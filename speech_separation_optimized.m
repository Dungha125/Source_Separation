%% speech_separation_optimized.m
% Script toi uu hoa de tach nguon noi cua moi nguoi tu mic1.wav va mic2.wav
% Muc tieu: Tach duoc giong noi cua tung nguoi ro rang nhat
clear all; close all; clc;

disp('=== TACH NGUON NOI TOI UU ===');
disp('Muc tieu: Tach duoc giong noi cua moi nguoi');
disp('');

%% 1. THAM SO TOI UU CHO GIONG NOI
fs_target = 16000;  % Tan so lay mau phu hop cho giong noi
result_folder = 'result_speech_separation';
if ~exist(result_folder, 'dir'), mkdir(result_folder); end

% Tham so cho giong noi
speech_low_freq = 300;    % Tan so thap nhat cua giong noi (Hz)
speech_high_freq = 3400;  % Tan so cao nhat cua giong noi (Hz)

%% 2. TAI FILE MICROPHONE
disp('--- Dang tai file microphone ---');
mic_files = {'mic1.wav', 'mic2.wav'};
X_data = cell(length(mic_files), 1);

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
    
    X_data{i} = x;
end

% Dua ve cung do dai
lengths = cellfun(@length, X_data);
min_len = min(lengths);
X = zeros(length(mic_files), min_len);
for i = 1:length(mic_files)
    X(i, :) = X_data{i}(1:min_len)';
end
disp(['  -> Do dai chung: ' num2str(min_len/fs_target) ' giay']);
disp('');

%% 3. TIEN XU LY TOI UU CHO GIONG NOI
disp('--- Tien xu ly toi uu cho giong noi ---');

% 3.1. Pre-emphasis filter (tang cuong tan so cao cua giong noi)
alpha_pre = 0.97;  % Tham so pre-emphasis
X_pre = zeros(size(X));
for i = 1:size(X, 1)
    X_pre(i, :) = filter([1, -alpha_pre], 1, X(i, :));
end
disp('  -> Da ap dung pre-emphasis filter');

% 3.2. Bandpass filter cho tan so giong noi
[b_bp, a_bp] = butter(4, [speech_low_freq, speech_high_freq]/(fs_target/2), 'bandpass');
X_filt = zeros(size(X_pre));
for i = 1:size(X_pre, 1)
    X_filt(i, :) = filter(b_bp, a_bp, X_pre(i, :));
end
disp(['  -> Da loc bandpass (' num2str(speech_low_freq) '-' num2str(speech_high_freq) ' Hz)']);

% 3.3. High-pass filter de loai bo nhieu DC va tan so rat thap
[b_hp, a_hp] = butter(4, 80/(fs_target/2), 'high');
for i = 1:size(X_filt, 1)
    X_filt(i, :) = filter(b_hp, a_hp, X_filt(i, :));
end

% 3.4. Chuan hoa
for i = 1:size(X_filt, 1)
    X_filt(i, :) = X_filt(i, :) / (max(abs(X_filt(i, :))) + eps);
end
disp('  -> Da chuan hoa tin hieu');
disp('');

%% 4. TACH NGUON BANG BSS-ICA (TOI UU CHO NOI)
disp('--- Tach nguon bang BSS-ICA (toi uu cho noi) ---');
try
    % Goi ham fastica_robust
    [S_ica, A_est, W] = fastica_robust(X_filt);
    
    num_sources = size(S_ica, 1);
    disp(['  -> Tim thay ' num2str(num_sources) ' nguon am thanh']);
    
    % Chuan hoa cac nguon
    S_normalized = zeros(size(S_ica));
    for i = 1:num_sources
        S_normalized(i, :) = S_ica(i, :) / (max(abs(S_ica(i, :))) + eps);
    end
    
    % Hau xu ly: De-emphasis (dao nguoc pre-emphasis)
    S_deemph = zeros(size(S_normalized));
    for i = 1:num_sources
        S_deemph(i, :) = filter(1, [1, -alpha_pre], S_normalized(i, :));
    end
    
    % Luu file
    for i = 1:num_sources
        output_file = fullfile(result_folder, sprintf('speaker_%d_ica.wav', i));
        audiowrite(output_file, S_deemph(i, :)', fs_target);
        disp(['  + Luu: ' output_file]);
    end
    
    disp('  -> Da luu ket qua BSS-ICA');
catch ME
    disp(['  Loi ICA: ' ME.message]);
    S_deemph = [];
end
disp('');

%% 5. TACH NGUON BANG GSC BEAMFORMER (TOI UU CHO NOI)
disp('--- Tach nguon bang GSC Beamformer (toi uu cho noi) ---');
try
    % Thử nhiều hướng hơn để tìm giọng nói tốt nhất
    look_directions = [-60, -45, -30, -15, 0, 15, 30, 45, 60];  % Nhiều hướng hơn
    filter_length = 256;  % Tăng độ dài filter cho giọng nói
    mu = 0.005;  % Giảm step size để ổn định hơn cho giọng nói
    
    gsc_outputs = [];
    gsc_energies = [];
    
    for dir_idx = 1:length(look_directions)
        theta = look_directions(dir_idx);
        [y_enhanced, ~] = gsc_beamformer(X_filt, fs_target, theta, filter_length, mu);
        
        % Tính năng lượng để chọn hướng tốt nhất
        energy = sum(y_enhanced.^2);
        gsc_outputs = [gsc_outputs; y_enhanced];
        gsc_energies = [gsc_energies, energy];
        
        disp(['  + Hướng ' num2str(theta) ' độ: nang luong = ' num2str(energy, '%.2e')]);
    end
    
    % Chọn 2 hướng có năng lượng cao nhất (giả định 2 người nói)
    [~, sorted_idx] = sort(gsc_energies, 'descend');
    best_directions = sorted_idx(1:min(2, length(sorted_idx)));
    
    num_gsc_sources = length(best_directions);
    gsc_selected = gsc_outputs(best_directions, :);
    
    % Chuẩn hóa
    gsc_normalized = zeros(size(gsc_selected));
    for i = 1:num_gsc_sources
        gsc_normalized(i, :) = gsc_selected(i, :) / (max(abs(gsc_selected(i, :))) + eps);
    end
    
    % De-emphasis
    gsc_deemph = zeros(size(gsc_normalized));
    for i = 1:num_gsc_sources
        gsc_deemph(i, :) = filter(1, [1, -alpha_pre], gsc_normalized(i, :));
    end
    
    % Luu file
    for i = 1:num_gsc_sources
        output_file = fullfile(result_folder, sprintf('speaker_%d_gsc_%ddeg.wav', i, look_directions(best_directions(i))));
        audiowrite(output_file, gsc_deemph(i, :)', fs_target);
        disp(['  + Luu: ' output_file ' (huong tot nhat: ' num2str(look_directions(best_directions(i))) ' do)']);
    end
    
    disp('  -> Da luu ket qua GSC');
catch ME
    disp(['  Loi GSC: ' ME.message]);
end
disp('');

%% 6. HAU XU LY TIN HIEU GIONG NOI
disp('--- Hau xu ly tin hieu giong noi ---');
try
    % Noise gate: Loại bỏ phần im lặng/nhiễu
    noise_gate_threshold = 0.01;  % Ngưỡng cho noise gate
    
    % Xử lý từng nguồn đã tách
    if exist('S_deemph', 'var') && ~isempty(S_deemph)
        num_sources_final = size(S_deemph, 1);
        S_enhanced = zeros(size(S_deemph));
        
        for i = 1:num_sources_final
            sig = S_deemph(i, :);
            
            % Noise gate
            sig_abs = abs(sig);
            gate_mask = sig_abs > noise_gate_threshold;
            sig_gated = sig .* gate_mask;
            
            % Soft thresholding để làm mượt transition
            sig_smooth = filter([0.2, 0.6, 0.2], 1, gate_mask);
            sig_gated = sig .* sig_smooth;
            
            % Normalize lại
            sig_gated = sig_gated / (max(abs(sig_gated)) + eps);
            
            S_enhanced(i, :) = sig_gated;
            
            % Luu file da hau xu ly
            output_file = fullfile(result_folder, sprintf('speaker_%d_ica_enhanced.wav', i));
            audiowrite(output_file, S_enhanced(i, :)', fs_target);
            disp(['  + Luu: ' output_file ' (da hau xu ly)']);
        end
    end
    
    disp('  -> Da hau xu ly tin hieu');
catch ME
    disp(['  Loi hau xu ly: ' ME.message]);
end
disp('');

%% 7. HIEN THI KET QUA
disp('--- Hien thi ket qua ---');
t = (0:min_len-1) / fs_target;

% Figure: So sanh cac phuong phap
figure('Name', 'So sanh ket qua tach giong noi', 'Position', [100, 100, 1400, 900]);

if exist('S_enhanced', 'var') && ~isempty(S_enhanced)
    num_show = size(S_enhanced, 1);
    
    for i = 1:num_show
        % Waveform
        subplot(num_show, 3, (i-1)*3 + 1);
        plot(t, S_enhanced(i, :));
        xlabel('Thoi gian (s)');
        ylabel('Bien do');
        title(['Nguoi ' num2str(i) ' - BSS-ICA (da toi uu)']);
        grid on;
        xlim([0, max(t)]);
        
        % Spectrogram
        subplot(num_show, 3, (i-1)*3 + 2);
        NFFT = 2048;
        NOVERLAP = floor(NFFT * 0.75);
        WINDOW = hamming(NFFT);
        [S_spec, F_spec, T_spec] = spectrogram(S_enhanced(i, :), WINDOW, NOVERLAP, NFFT, fs_target);
        imagesc(T_spec, F_spec, 20*log10(abs(S_spec) + eps));
        axis xy; colorbar;
        xlabel('Thoi gian (s)');
        ylabel('Tan so (Hz)');
        title(['Spectrogram - Nguoi ' num2str(i)]);
        ylim([0, min(4000, fs_target/2)]);
        
        % Phổ tần số trung bình
        subplot(num_show, 3, (i-1)*3 + 3);
        Pxx = mean(abs(S_spec).^2, 2);
        plot(F_spec, 10*log10(Pxx + eps));
        xlabel('Tan so (Hz)');
        ylabel('Cong suat (dB)');
        title(['Phổ tan so - Nguoi ' num2str(i)]);
        grid on;
        xlim([0, min(4000, fs_target/2)]);
    end
end

disp('  -> Da hien thi ket qua');
disp('');

%% 8. KET LUAN
disp('=== HOAN THANH ===');
disp(['Ket qua da duoc luu trong thu muc: ' result_folder]);
disp('');
disp('File ket qua:');
disp('  - speaker_1_ica.wav: Nguoi 1 (BSS-ICA)');
disp('  - speaker_2_ica.wav: Nguoi 2 (BSS-ICA)');
disp('  - speaker_1_ica_enhanced.wav: Nguoi 1 (BSS-ICA + hau xu ly)');
disp('  - speaker_2_ica_enhanced.wav: Nguoi 2 (BSS-ICA + hau xu ly)');
disp('  - speaker_X_gsc_Xdeg.wav: Nguoi X (GSC Beamformer)');
disp('');
disp('Luu y: Hay nghe thu cac file de chon phuong phap tot nhat!');

