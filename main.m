%% main.m - Enhanced Source Separation
% Luu y: Can co 2 file phu tro la 'my_istft.m' va 'fastica_robust.m'
clear all; close all; clc;
format compact;

disp('=== KHOI TAO HE THONG ===');

%% 1. CAI DAT THAM SO
rng(1,'twister'); 

fs_target = 16000; 
NFFT = 4096;
noverlap_ratio = 0.75;
win_type = 'hamming';
num_sources = 2; 

result_folder = 'result_enhanced';
if ~exist(result_folder, 'dir'), mkdir(result_folder); end

%% 2. TAI DU LIEU DAU VAO (MIC1 VA MIC2)
disp('--- Dang tai du lieu dau vao ---');
files = {'mic1.wav', 'mic2.wav'};
X_data = cell(length(files), 1);

try
    % Buoc 1: Tai tat ca cac file va luu vao cell array
    lengths = [];
    for i = 1:length(files)
        [x, fs_in] = audioread(files{i});
        if fs_in ~= fs_target
            x = resample(x, fs_target, fs_in);
        end
        % Chuyen thanh vector cot (column vector)
        if size(x, 2) > 1
            x = x(:, 1); % Lay kenh dau tien neu stereo
        end
        X_data{i} = x;
        lengths = [lengths, length(x)];
        disp(['  + ' files{i} ': ' num2str(fs_in) ' Hz, ' num2str(length(x)/fs_target) ' giay']);
    end
    
    % Buoc 2: Tim do dai ngan nhat
    min_len = min(lengths);
    
    % Buoc 3: Cat tat ca ve cung do dai va tao ma tran X
    X = zeros(length(files), min_len);
    for i = 1:length(files)
        X(i, :) = X_data{i}(1:min_len)';
    end
    
    disp(['Da load ' num2str(length(files)) ' file microphone tu dia.']);
    disp(['Do dai tin hieu: ' num2str(min_len/fs_target) ' giay']);
catch ME
    error(['Loi khi tai file: ' ME.message]);
end

% Luu tin hieu dau vao de kiem tra
audiowrite(fullfile(result_folder, 'input_mic1.wav'), X(1,:)', fs_target);
audiowrite(fullfile(result_folder, 'input_mic2.wav'), X(2,:)', fs_target);

%% 3. TIEN XU LY
disp('--- Tien xu ly (Wiener) ---');
[b_hp, a_hp] = butter(4, 50/(fs_target/2), 'high');
X_filt = filter(b_hp, a_hp, X, [], 2);

if strcmp(win_type, 'hamming'), WINDOW = hamming(NFFT); else, WINDOW = hanning(NFFT); end
NOVERLAP = floor(length(WINDOW) * noverlap_ratio);

X_denoised = zeros(size(X_filt));
for ch = 1:size(X_filt,1)
    [S_noisy, F, T] = spectrogram(X_filt(ch,:), WINDOW, NOVERLAP, NFFT, fs_target);
    
    % --- FIX LOI "min" TRUOC DO ---
    num_frames = size(S_noisy, 2);
    frames_to_use = min(10, num_frames); 
    noise_profile = mean(abs(S_noisy(:, 1:frames_to_use)).^2, 2);
    
    % Wiener Filter
    P_y = abs(S_noisy).^2;
    P_n = repmat(noise_profile, 1, size(P_y, 2));
    Gain = max(0, (P_y - P_n) ./ (P_y + eps));
    Gain = medfilt2(Gain, [3 1]); 
    S_clean = S_noisy .* Gain;
    
    % --- FIX LOI "Dimension Mismatch" TAI DAY ---
    % 1. Lay tin hieu ra mot bien tam
    temp_sig = my_istft(S_clean, WINDOW, NOVERLAP, NFFT, size(X_filt, 2));
    
    % 2. Tinh do dai hop le de gan (lay min giua tin hieu goc va tai tao)
    len_assign = min(length(temp_sig), size(X_denoised, 2));
    
    % 3. Chi gan phan du lieu hop le vao ma tran
    X_denoised(ch, 1:len_assign) = temp_sig(1:len_assign);
end
audiowrite(fullfile(result_folder, 'mix_denoised.wav'), X_denoised', fs_target);

%% 4. CLUSTERING
disp('--- Phuong phap 1: Clustering ---');
[S1, ~, ~] = spectrogram(X_denoised(1,:), WINDOW, NOVERLAP, NFFT, fs_target);
[S2, ~, ~] = spectrogram(X_denoised(2,:), WINDOW, NOVERLAP, NFFT, fs_target);

IPD = angle(S2 ./ (S1 + eps));
ILD = 20*log10((abs(S2) + eps) ./ (abs(S1) + eps));
mag_sum = abs(S1) + abs(S2);
threshold = prctile(mag_sum(:), 85); 
mask_active = mag_sum > threshold;
features = [IPD(mask_active), ILD(mask_active)/10]; 

try
    [idx_cluster, C] = kmeans(features, num_sources, 'Replicates', 3);
    figure; scatter(features(:,1), features(:,2)*10, 10, idx_cluster, 'filled'); title('Clustering');
    
    for i = 1:num_sources
        dist_map = sqrt((IPD - C(i,1)).^2 + ((ILD/10) - C(i,2)).^2);
        mask_i = exp(-dist_map.^2 / 0.5);
        
        % Goi ham my_istft
        x_sep = my_istft(S1 .* mask_i, WINDOW, NOVERLAP, NFFT, length(X_denoised));
        audiowrite(fullfile(result_folder, sprintf('output_cluster_src%d.wav', i)), x_sep/max(abs(x_sep)), fs_target);
    end
    disp('Da luu file cluster.');
catch
    disp('Loi Clustering. Bo qua.');
end

%% 5. ICA
disp('--- Phuong phap 2: ICA ---');
try
    % Goi ham fastica_robust
    [S_ica, A_est, W] = fastica_robust(X_denoised); 
    
    num_sources_ica = size(S_ica, 1);
    S_ica_normalized = zeros(size(S_ica));
    for i = 1:num_sources_ica
        S_ica_normalized(i, :) = S_ica(i, :) / (max(abs(S_ica(i, :))) + eps);
    end
    
    % Thoi gian
    t = (0:size(X_denoised, 2)-1) / fs_target;
    
    %% Figure: Waveform tin hieu dau vao
    figure('Name', 'Tin hieu dau vao (ICA)', 'Position', [100, 100, 1200, 600]);
    subplot(2,1,1);
    plot(t, X_denoised(1,:));
    xlabel('Thoi gian (s)');
    ylabel('Bien do');
    title('Microphone 1 - Tin hieu da tien xu ly');
    grid on;
    xlim([0, max(t)]);
    
    subplot(2,1,2);
    plot(t, X_denoised(2,:));
    xlabel('Thoi gian (s)');
    ylabel('Bien do');
    title('Microphone 2 - Tin hieu da tien xu ly');
    grid on;
    xlim([0, max(t)]);
    
    %% Figure: Waveform cac nguon da tach bang ICA
    figure('Name', 'Cac nguon da tach bang ICA', 'Position', [150, 150, 1200, 600]);
    for i = 1:num_sources_ica
        subplot(num_sources_ica, 1, i);
        plot(t, S_ica_normalized(i, :));
        xlabel('Thoi gian (s)');
        ylabel('Bien do');
        title(['Nguon ' num2str(i) ' - Tach bang BSS-ICA']);
        grid on;
        xlim([0, max(t)]);
    end
    
    %% Figure: Spectrogram tin hieu dau vao
    figure('Name', 'Spectrogram tin hieu dau vao (ICA)', 'Position', [200, 200, 1200, 600]);
    subplot(2,1,1);
    [S1_ica, F1_ica, T1_ica] = spectrogram(X_denoised(1,:), WINDOW, NOVERLAP, NFFT, fs_target);
    imagesc(T1_ica, F1_ica, 20*log10(abs(S1_ica) + eps));
    axis xy; colorbar;
    xlabel('Thoi gian (s)');
    ylabel('Tan so (Hz)');
    title('Microphone 1 - Spectrogram');
    ylim([0, min(8000, fs_target/2)]);
    
    subplot(2,1,2);
    [S2_ica, F2_ica, T2_ica] = spectrogram(X_denoised(2,:), WINDOW, NOVERLAP, NFFT, fs_target);
    imagesc(T2_ica, F2_ica, 20*log10(abs(S2_ica) + eps));
    axis xy; colorbar;
    xlabel('Thoi gian (s)');
    ylabel('Tan so (Hz)');
    title('Microphone 2 - Spectrogram');
    ylim([0, min(8000, fs_target/2)]);
    
    %% Figure: Spectrogram cac nguon da tach
    figure('Name', 'Spectrogram cac nguon da tach (ICA)', 'Position', [250, 250, 1200, 600]);
    for i = 1:num_sources_ica
        subplot(num_sources_ica, 1, i);
        [S_sep_ica, F_sep_ica, T_sep_ica] = spectrogram(S_ica_normalized(i, :), WINDOW, NOVERLAP, NFFT, fs_target);
        imagesc(T_sep_ica, F_sep_ica, 20*log10(abs(S_sep_ica) + eps));
        axis xy; colorbar;
        xlabel('Thoi gian (s)');
        ylabel('Tan so (Hz)');
        title(['Nguon ' num2str(i) ' - Spectrogram']);
        ylim([0, min(8000, fs_target/2)]);
    end
    
    %% Figure: So sanh tong quan
    figure('Name', 'So sanh tong quan (ICA)', 'Position', [300, 300, 1400, 800]);
    
    subplot(2,2,1);
    plot(t, X_denoised(1,:));
    xlabel('Thoi gian (s)');
    ylabel('Bien do');
    title('Dau vao: Microphone 1');
    grid on;
    xlim([0, max(t)]);
    
    subplot(2,2,2);
    plot(t, X_denoised(2,:));
    xlabel('Thoi gian (s)');
    ylabel('Bien do');
    title('Dau vao: Microphone 2');
    grid on;
    xlim([0, max(t)]);
    
    subplot(2,2,3);
    plot(t, S_ica_normalized(1, :));
    xlabel('Thoi gian (s)');
    ylabel('Bien do');
    title('Ket qua: Nguon 1 (ICA)');
    grid on;
    xlim([0, max(t)]);
    
    subplot(2,2,4);
    if num_sources_ica >= 2
        plot(t, S_ica_normalized(2, :));
        xlabel('Thoi gian (s)');
        ylabel('Bien do');
        title('Ket qua: Nguon 2 (ICA)');
        grid on;
        xlim([0, max(t)]);
    end
    
    % Luu file
    for i = 1:num_sources_ica
        sig = S_ica_normalized(i, :);
        audiowrite(fullfile(result_folder, sprintf('output_ica_src%d.wav', i)), sig', fs_target);
    end
    disp('Da luu file ICA.');
catch ME
    disp(['Loi ICA: ' ME.message]);
end

%% 6. GSC BEAMFORMER
disp('--- Phuong phap 3: GSC Beamformer ---');
try
    % Thử các hướng khác nhau để tách nguồn
    look_directions = [-45, -30, -15, 0, 15, 30, 45];  % Các hướng để thử (độ)
    filter_length = 128;
    mu = 0.01;
    
    gsc_outputs = [];
    for dir_idx = 1:length(look_directions)
        theta = look_directions(dir_idx);
        [y_enhanced, ~] = gsc_beamformer(X_denoised, fs_target, theta, filter_length, mu);
        gsc_outputs = [gsc_outputs; y_enhanced];
        disp(['  + Hướng ' num2str(theta) ' độ: hoàn thành']);
    end
    
    num_gsc_sources = size(gsc_outputs, 1);
    
    % Chuẩn hóa
    gsc_normalized = zeros(size(gsc_outputs));
    for i = 1:num_gsc_sources
        gsc_normalized(i, :) = gsc_outputs(i, :) / (max(abs(gsc_outputs(i, :))) + eps);
    end
    
    % Thoi gian
    if ~exist('t', 'var')
        t = (0:size(X_denoised, 2)-1) / fs_target;
    end
    
    %% Figure: Waveform cac nguon da tach bang GSC
    figure('Name', 'Cac nguon da tach bang GSC', 'Position', [400, 400, 1200, 600]);
    for i = 1:num_gsc_sources
        subplot(num_gsc_sources, 1, i);
        plot(t(1:length(gsc_normalized(i, :))), gsc_normalized(i, :));
        xlabel('Thoi gian (s)');
        ylabel('Bien do');
        title(['Nguon ' num2str(i) ' - GSC Beamformer (huong ' num2str(look_directions(i)) ' do)']);
        grid on;
        xlim([0, max(t)]);
    end
    
    %% Figure: Spectrogram cac nguon da tach bang GSC
    figure('Name', 'Spectrogram cac nguon da tach (GSC)', 'Position', [450, 450, 1200, 600]);
    for i = 1:num_gsc_sources
        subplot(num_gsc_sources, 1, i);
        [S_gsc, F_gsc, T_gsc] = spectrogram(gsc_normalized(i, :), WINDOW, NOVERLAP, NFFT, fs_target);
        imagesc(T_gsc, F_gsc, 20*log10(abs(S_gsc) + eps));
        axis xy; colorbar;
        xlabel('Thoi gian (s)');
        ylabel('Tan so (Hz)');
        title(['Nguon ' num2str(i) ' - GSC Beamformer (huong ' num2str(look_directions(i)) ' do)']);
        ylim([0, min(8000, fs_target/2)]);
    end
    
    % Luu file
    for i = 1:num_gsc_sources
        output_file = fullfile(result_folder, sprintf('output_gsc_src%d.wav', i));
        audiowrite(output_file, gsc_normalized(i, :)', fs_target);
    end
    disp('Da luu file GSC.');
catch ME
    disp(['Loi GSC: ' ME.message]);
end

disp('=== HOAN THANH ===');