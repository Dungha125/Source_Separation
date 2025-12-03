%% bss_ica_separate.m
% Script don gian de tach nguon am thanh tu mic1.wav va mic2.wav bang BSS-ICA
% Su dung: Chay script nay trong MATLAB
clear all; close all; clc;

disp('=== TACH NGUON AM THANH BANG BSS-ICA ===');
disp('');

%% 1. THAM SO
fs_target = 16000;  % Tan so lay mau mong muon
result_folder = 'result_bss_ica';
if ~exist(result_folder, 'dir'), mkdir(result_folder); end

%% 2. TAI FILE MICROPHONE
disp('--- Dang tai file microphone ---');
mic_files = {'mic1.wav', 'mic2.wav'};
X = [];

for i = 1:length(mic_files)
    if ~exist(mic_files{i}, 'file')
        error(['Khong tim thay file: ' mic_files{i}]);
    end
    
    [x, fs_in] = audioread(mic_files{i});
    disp(['  + ' mic_files{i} ': ' num2str(fs_in) ' Hz, ' num2str(length(x)/fs_in) ' giay']);
    
    % Resample neu can
    if fs_in ~= fs_target
        x = resample(x, fs_target, fs_in);
        disp(['    -> Resample ve ' num2str(fs_target) ' Hz']);
    end
    
    % Chuyen thanh mono neu can
    if size(x, 2) > 1
        x = mean(x, 2); % Lay trung binh 2 kenh
    end
    
    X = [X; x']; % Them vao ma tran (moi hang la mot microphone)
end

% Dua ve cung do dai
min_len = size(X, 2);
X = X(:, 1:min_len);
disp(['  -> Do dai chung: ' num2str(min_len/fs_target) ' giay']);
disp('');

%% 3. TIEN XU LY TOI UU CHO GIONG NOI
disp('--- Tien xu ly toi uu cho giong noi ---');

% 3.1. Pre-emphasis filter (tang cuong tan so cao cua giong noi)
alpha_pre = 0.97;
X_pre = zeros(size(X));
for i = 1:size(X, 1)
    X_pre(i, :) = filter([1, -alpha_pre], 1, X(i, :));
end

% 3.2. Bandpass filter cho tan so giong noi (300-3400 Hz)
speech_low = 300;
speech_high = 3400;
[b_bp, a_bp] = butter(4, [speech_low, speech_high]/(fs_target/2), 'bandpass');
X_filt = zeros(size(X_pre));
for i = 1:size(X_pre, 1)
    X_filt(i, :) = filter(b_bp, a_bp, X_pre(i, :));
end

% 3.3. High-pass filter de loai bo nhieu DC
[b_hp, a_hp] = butter(4, 80/(fs_target/2), 'high');
for i = 1:size(X_filt, 1)
    X_filt(i, :) = filter(b_hp, a_hp, X_filt(i, :));
end

% Chuan hoa
for i = 1:size(X_filt, 1)
    X_filt(i, :) = X_filt(i, :) / (max(abs(X_filt(i, :))) + eps);
end

disp('  -> Da loc va chuan hoa tin hieu (toi uu cho giong noi)');
disp('');

%% 4. TACH NGUON BANG BSS-ICA
disp('--- Tach nguon bang BSS-ICA ---');
try
    % Goi ham fastica_robust
    [S_ica, A_est, W] = fastica_robust(X_filt);
    
    num_sources = size(S_ica, 1);
    disp(['  -> Tim thay ' num2str(num_sources) ' nguon am thanh']);
    disp('');
    
    %% 5. HIEN THI VA LUU KET QUA
    disp('--- Hien thi va luu ket qua ---');
    
    % Chuan hoa cac nguon
    S_normalized = zeros(size(S_ica));
    for i = 1:num_sources
        S_normalized(i, :) = S_ica(i, :) / (max(abs(S_ica(i, :))) + eps);
    end
    
    % Thoi gian
    t = (0:min_len-1) / fs_target;
    
    %% Figure 1: Waveform tin hieu dau vao
    figure('Name', 'Tin hieu dau vao', 'Position', [100, 100, 1200, 600]);
    subplot(2,1,1);
    plot(t, X_filt(1,:));
    xlabel('Thoi gian (s)');
    ylabel('Bien do');
    title('Microphone 1 - Tin hieu dau vao');
    grid on;
    xlim([0, max(t)]);
    
    subplot(2,1,2);
    plot(t, X_filt(2,:));
    xlabel('Thoi gian (s)');
    ylabel('Bien do');
    title('Microphone 2 - Tin hieu dau vao');
    grid on;
    xlim([0, max(t)]);
    
    %% Figure 2: Waveform cac nguon da tach
    figure('Name', 'Cac nguon da tach', 'Position', [150, 150, 1200, 600]);
    for i = 1:num_sources
        subplot(num_sources, 1, i);
        plot(t, S_normalized(i, :));
        xlabel('Thoi gian (s)');
        ylabel('Bien do');
        title(['Nguon ' num2str(i) ' - Tach bang BSS-ICA']);
        grid on;
        xlim([0, max(t)]);
    end
    
    %% Figure 3: Spectrogram tin hieu dau vao
    figure('Name', 'Spectrogram tin hieu dau vao', 'Position', [200, 200, 1200, 600]);
    NFFT_plot = 2048;
    NOVERLAP_plot = floor(NFFT_plot * 0.75);
    WINDOW_plot = hamming(NFFT_plot);
    
    subplot(2,1,1);
    [S1, F1, T1] = spectrogram(X_filt(1,:), WINDOW_plot, NOVERLAP_plot, NFFT_plot, fs_target);
    imagesc(T1, F1, 20*log10(abs(S1) + eps));
    axis xy; colorbar;
    xlabel('Thoi gian (s)');
    ylabel('Tan so (Hz)');
    title('Microphone 1 - Spectrogram');
    ylim([0, min(8000, fs_target/2)]);
    
    subplot(2,1,2);
    [S2, F2, T2] = spectrogram(X_filt(2,:), WINDOW_plot, NOVERLAP_plot, NFFT_plot, fs_target);
    imagesc(T2, F2, 20*log10(abs(S2) + eps));
    axis xy; colorbar;
    xlabel('Thoi gian (s)');
    ylabel('Tan so (Hz)');
    title('Microphone 2 - Spectrogram');
    ylim([0, min(8000, fs_target/2)]);
    
    %% Figure 4: Spectrogram cac nguon da tach
    figure('Name', 'Spectrogram cac nguon da tach', 'Position', [250, 250, 1200, 600]);
    for i = 1:num_sources
        subplot(num_sources, 1, i);
        [S_sep, F_sep, T_sep] = spectrogram(S_normalized(i, :), WINDOW_plot, NOVERLAP_plot, NFFT_plot, fs_target);
        imagesc(T_sep, F_sep, 20*log10(abs(S_sep) + eps));
        axis xy; colorbar;
        xlabel('Thoi gian (s)');
        ylabel('Tan so (Hz)');
        title(['Nguon ' num2str(i) ' - Spectrogram']);
        ylim([0, min(8000, fs_target/2)]);
    end
    
    %% Figure 5: So sanh tong quan (subplot 2x2)
    figure('Name', 'So sanh tong quan', 'Position', [300, 300, 1400, 800]);
    
    % Waveform mic1
    subplot(2,2,1);
    plot(t, X_filt(1,:));
    xlabel('Thoi gian (s)');
    ylabel('Bien do');
    title('Dau vao: Microphone 1');
    grid on;
    xlim([0, max(t)]);
    
    % Waveform mic2
    subplot(2,2,2);
    plot(t, X_filt(2,:));
    xlabel('Thoi gian (s)');
    ylabel('Bien do');
    title('Dau vao: Microphone 2');
    grid on;
    xlim([0, max(t)]);
    
    % Waveform nguon 1
    subplot(2,2,3);
    plot(t, S_normalized(1, :));
    xlabel('Thoi gian (s)');
    ylabel('Bien do');
    title('Ket qua: Nguon 1');
    grid on;
    xlim([0, max(t)]);
    
    % Waveform nguon 2
    subplot(2,2,4);
    plot(t, S_normalized(2, :));
    xlabel('Thoi gian (s)');
    ylabel('Bien do');
    title('Ket qua: Nguon 2');
    grid on;
    xlim([0, max(t)]);
    
    % De-emphasis (dao nguoc pre-emphasis)
    alpha_pre = 0.97;
    S_deemph = zeros(size(S_normalized));
    for i = 1:num_sources
        S_deemph(i, :) = filter(1, [1, -alpha_pre], S_normalized(i, :));
        S_deemph(i, :) = S_deemph(i, :) / (max(abs(S_deemph(i, :))) + eps);
    end
    
    %% Luu file
    for i = 1:num_sources
        output_file = fullfile(result_folder, sprintf('speaker_%d_ica.wav', i));
        audiowrite(output_file, S_deemph(i, :)', fs_target);
        disp(['  + Luu: ' output_file ' (nguoi ' num2str(i) ')']);
    end
    
    % Luu tin hieu dau vao de so sanh
    audiowrite(fullfile(result_folder, 'input_mic1.wav'), X_filt(1,:)', fs_target);
    audiowrite(fullfile(result_folder, 'input_mic2.wav'), X_filt(2,:)', fs_target);
    
catch ME
    disp(['LOI ICA: ' ME.message]);
end

%% 6. TACH NGUON BANG GSC BEAMFORMER
disp('--- Tach nguon bang GSC Beamformer ---');
try
    % Thử nhiều hướng hơn để tìm giọng nói tốt nhất
    look_directions = [-60, -45, -30, -15, 0, 15, 30, 45, 60];  % Nhiều hướng hơn
    filter_length = 256;  % Tăng độ dài filter cho giọng nói
    mu = 0.005;  % Giảm step size để ổn định hơn cho giọng nói
    
    gsc_outputs = [];
    for dir_idx = 1:length(look_directions)
        theta = look_directions(dir_idx);
        [y_enhanced, ~] = gsc_beamformer(X_filt, fs_target, theta, filter_length, mu);
        gsc_outputs = [gsc_outputs; y_enhanced];
        disp(['  + Hướng ' num2str(theta) ' độ: hoàn thành']);
    end
    
    num_gsc_sources = size(gsc_outputs, 1);
    
    % Chọn 2 hướng có năng lượng cao nhất (giả định 2 người nói)
    gsc_energies = sum(gsc_outputs.^2, 2);
    [~, sorted_idx] = sort(gsc_energies, 'descend');
    best_directions = sorted_idx(1:min(2, length(sorted_idx)));
    
    gsc_selected = gsc_outputs(best_directions, :);
    num_gsc_sources = size(gsc_selected, 1);
    
    % Chuẩn hóa
    gsc_normalized = zeros(size(gsc_selected));
    for i = 1:num_gsc_sources
        gsc_normalized(i, :) = gsc_selected(i, :) / (max(abs(gsc_selected(i, :))) + eps);
    end
    
    % De-emphasis
    alpha_pre = 0.97;
    gsc_deemph = zeros(size(gsc_normalized));
    for i = 1:num_gsc_sources
        gsc_deemph(i, :) = filter(1, [1, -alpha_pre], gsc_normalized(i, :));
        gsc_deemph(i, :) = gsc_deemph(i, :) / (max(abs(gsc_deemph(i, :))) + eps);
    end
    
    % Thoi gian
    if ~exist('t', 'var')
        t = (0:size(X_filt, 2)-1) / fs_target;
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
    NFFT_gsc = 2048;
    NOVERLAP_gsc = floor(NFFT_gsc * 0.75);
    WINDOW_gsc = hamming(NFFT_gsc);
    
    for i = 1:num_gsc_sources
        subplot(num_gsc_sources, 1, i);
        [S_gsc, F_gsc, T_gsc] = spectrogram(gsc_normalized(i, :), WINDOW_gsc, NOVERLAP_gsc, NFFT_gsc, fs_target);
        imagesc(T_gsc, F_gsc, 20*log10(abs(S_gsc) + eps));
        axis xy; colorbar;
        xlabel('Thoi gian (s)');
        ylabel('Tan so (Hz)');
        title(['Nguon ' num2str(i) ' - GSC Beamformer (huong ' num2str(look_directions(i)) ' do)']);
        ylim([0, min(8000, fs_target/2)]);
    end
    
    % Luu file
    for i = 1:num_gsc_sources
        output_file = fullfile(result_folder, sprintf('speaker_%d_gsc_%ddeg.wav', i, look_directions(best_directions(i))));
        audiowrite(output_file, gsc_deemph(i, :)', fs_target);
        disp(['  + Luu: ' output_file ' (nguoi ' num2str(i) ', huong ' num2str(look_directions(best_directions(i))) ' do)']);
    end
    
    disp('Da luu file GSC.');
catch ME
    disp(['Loi GSC: ' ME.message]);
end

disp('');
disp('=== HOAN THANH ===');
disp(['Ket qua da duoc luu trong thu muc: ' result_folder]);

