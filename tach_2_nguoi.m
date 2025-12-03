%% tach_2_nguoi.m
% Script tach tieng 2 nguoi tu 2 file microphone - KET HOP NHIEU THUAT TOAN
% Dau vao: mic1.wav va mic2.wav (2 microphone thu cung 2 nguoi noi)
% Dau ra: nguoi_1.wav va nguoi_2.wav (tieng rieng cua tung nguoi)
clear all; close all; clc;

disp('=============================================================');
disp('    TACH TIENG 2 NGUOI - KET HOP NHIEU THUAT TOAN           ');
disp('=============================================================');
disp('');

%% 1. THAM SO
fs_target = 16000;
result_folder = 'output_tach_2_nguoi';
if ~exist(result_folder, 'dir'), mkdir(result_folder); end

%% 2. TAI FILE MICROPHONE
disp('BUOC 1: TAI FILE MICROPHONE');
disp('-----------------------------------------------------------');
mic_files = {'mic1.wav', 'mic2.wav'};
X_data = cell(2, 1);

for i = 1:2
    if ~exist(mic_files{i}, 'file')
        error(['Khong tim thay file: ' mic_files{i}]);
    end
    
    [x, fs_in] = audioread(mic_files{i});
    disp(['  + ' mic_files{i} ': ' num2str(fs_in) ' Hz, ' num2str(length(x)/fs_in) ' giay']);
    
    if fs_in ~= fs_target
        x = resample(x, fs_target, fs_in);
    end
    
    if size(x, 2) > 1
        x = mean(x, 2);
    end
    
    X_data{i} = x(:)';
end

len1 = length(X_data{1});
len2 = length(X_data{2});
min_len = min(len1, len2);

X = zeros(2, min_len);
X(1, :) = X_data{1}(1:min_len);
X(2, :) = X_data{2}(1:min_len);

disp(['  -> Da tai 2 file, do dai: ' num2str(min_len/fs_target, '%.2f') ' giay']);
disp('');

%% 3. TIEN XU LY
disp('BUOC 2: TIEN XU LY TIN HIEU');
disp('-----------------------------------------------------------');

alpha = 0.97;
X_pre = zeros(size(X));
for i = 1:2
    X_pre(i, :) = filter([1, -alpha], 1, X(i, :));
end

[b_bp, a_bp] = butter(6, [300, 3400]/(fs_target/2), 'bandpass');
X_bp = zeros(size(X_pre));
for i = 1:2
    X_bp(i, :) = filtfilt(b_bp, a_bp, X_pre(i, :));
end

for i = 1:2
    X_bp(i, :) = X_bp(i, :) / (max(abs(X_bp(i, :))) + eps);
end
disp('  + Pre-emphasis, Bandpass (300-3400 Hz), Chuan hoa');
disp('');

%% 4. PHUONG PHAP 1: BSS-ICA
disp('BUOC 3: TACH NGUON BANG NHIEU PHUONG PHAP');
disp('-----------------------------------------------------------');
disp('Phuong phap 1: BSS-ICA');

all_sources = {};
source_idx = 1;

try
    [S_ica, ~, ~] = fastica_robust(X_bp);
    
    if size(S_ica, 1) >= 2
        for i = 1:2
            sig = S_ica(i, :) / (max(abs(S_ica(i, :))) + eps);
            sig = filter(1, [1, -alpha], sig);
            sig = sig / (max(abs(sig)) + eps);
            
            all_sources{source_idx} = struct('signal', sig, 'method', 'ICA', 'index', i);
            source_idx = source_idx + 1;
        end
        disp('  + ICA: Tim thay 2 nguon');
    end
catch ME
    disp(['  ! Loi ICA: ' ME.message]);
end

%% 5. PHUONG PHAP 2: CLUSTERING (IPD/ILD)
disp('Phuong phap 2: Clustering (IPD/ILD)');

try
    NFFT = 2048;
    NOVERLAP = floor(NFFT * 0.75);
    WINDOW = hamming(NFFT);
    
    [S1, F, T] = spectrogram(X_bp(1,:), WINDOW, NOVERLAP, NFFT, fs_target);
    [S2, ~, ~] = spectrogram(X_bp(2,:), WINDOW, NOVERLAP, NFFT, fs_target);
    
    % Tinh IPD va ILD
    IPD = angle(S2 ./ (S1 + eps));
    ILD = 20*log10((abs(S2) + eps) ./ (abs(S1) + eps));
    
    % Chi lay cac diem co nang luong cao
    mag_sum = abs(S1) + abs(S2);
    threshold = prctile(mag_sum(:), 80);
    mask_active = mag_sum > threshold;
    
    features = [IPD(mask_active), ILD(mask_active)/10];
    
    % K-means clustering
    [idx_cluster, C] = kmeans(features, 2, 'Replicates', 5, 'MaxIter', 500);
    
    % Tao mask cho moi cluster
    for i = 1:2
        mask_full = zeros(size(IPD));
        
        % Gan lai cluster cho tat ca cac diem
        cluster_idx = 1;
        for f = 1:size(IPD, 1)
            for t = 1:size(IPD, 2)
                if mask_active(f, t)
                    if idx_cluster(cluster_idx) == i
                        mask_full(f, t) = 1;
                    end
                    cluster_idx = cluster_idx + 1;
                end
            end
        end
        
        % Tao soft mask tu khoang cach den center
        dist_map = sqrt((IPD - C(i,1)).^2 + ((ILD/10) - C(i,2)).^2);
        soft_mask = exp(-dist_map.^2 / 0.3);  % Mask mem hon
        
        % Ap dung mask
        S_masked = S1 .* soft_mask;
        
        % Chuyen ve mien thoi gian
        x_sep = my_istft(S_masked, WINDOW, NOVERLAP, NFFT, min_len);
        sig = x_sep(1:min(length(x_sep), min_len));
        
        if length(sig) < min_len
            sig = [sig, zeros(1, min_len - length(sig))];
        end
        
        sig = sig / (max(abs(sig)) + eps);
        
        all_sources{source_idx} = struct('signal', sig, 'method', 'Clustering', 'index', i);
        source_idx = source_idx + 1;
    end
    disp('  + Clustering: Tim thay 2 nguon');
catch ME
    disp(['  ! Loi Clustering: ' ME.message]);
end

%% 6. PHUONG PHAP 3: BEAMFORMING (THU NHIEU HUONG)
disp('Phuong phap 3: Beamforming');

try
    look_directions = [-60, -45, -30, 0, 30, 45, 60];
    beamform_results = [];
    beamform_energies = [];
    
    for theta = look_directions
        y = delay_and_sum_beamformer(X_bp, fs_target, theta);
        beamform_results = [beamform_results; y];
        beamform_energies = [beamform_energies, sum(y.^2)];
    end
    
    % Chon 2 huong co nang luong cao nhat
    [~, sorted_idx] = sort(beamform_energies, 'descend');
    
    for i = 1:2
        idx = sorted_idx(i);
        sig = beamform_results(idx, :);
        sig = filter(1, [1, -alpha], sig);
        sig = sig / (max(abs(sig)) + eps);
        
        all_sources{source_idx} = struct('signal', sig, 'method', 'Beamforming', 'index', look_directions(idx));
        source_idx = source_idx + 1;
    end
    disp('  + Beamforming: Tim thay 2 nguon');
catch ME
    disp(['  ! Loi Beamforming: ' ME.message]);
end

disp('');

%% 7. CHON 2 NGUON TOT NHAT
disp('BUOC 4: CHON 2 NGUON TOT NHAT');
disp('-----------------------------------------------------------');

num_sources = length(all_sources);
disp(['  + Co ' num2str(num_sources) ' nguon tu cac phuong phap']);

if num_sources < 2
    error('Khong du nguon de tach 2 nguoi!');
end

% Tinh correlation giua tat ca cac cap
best_corr = inf;
best_i = 1;
best_j = 2;

for i = 1:num_sources-1
    for j = i+1:num_sources
        sig1 = all_sources{i}.signal;
        sig2 = all_sources{j}.signal;
        
        % Dam bao cung do dai
        len = min(length(sig1), length(sig2));
        sig1 = sig1(1:len);
        sig2 = sig2(1:len);
        
        corr_val = corrcoef(sig1, sig2);
        if ~isnan(corr_val(1, 2))
            corr_abs = abs(corr_val(1, 2));
            
            % Muon correlation thap nhat (khac nhau nhat)
            if corr_abs < best_corr
                best_corr = corr_abs;
                best_i = i;
                best_j = j;
            end
        end
    end
end

disp(['  + Cap tot nhat: ' all_sources{best_i}.method ' va ' all_sources{best_j}.method]);
disp(['  + Correlation: ' num2str(best_corr, '%.3f')]);

nguoi_1 = all_sources{best_i}.signal;
nguoi_2 = all_sources{best_j}.signal;

% Dam bao cung do dai
if length(nguoi_1) < min_len
    nguoi_1 = [nguoi_1, zeros(1, min_len - length(nguoi_1))];
else
    nguoi_1 = nguoi_1(1:min_len);
end

if length(nguoi_2) < min_len
    nguoi_2 = [nguoi_2, zeros(1, min_len - length(nguoi_2))];
else
    nguoi_2 = nguoi_2(1:min_len);
end

disp('');

%% 8. LOAI BO PHAN DINH BANG WIENER MASKING
disp('BUOC 5: LOAI BO PHAN DINH BANG WIENER MASKING');
disp('-----------------------------------------------------------');

if best_corr > 0.3
    disp('  + Phat hien phan dinh, ap dung Wiener masking...');
    
    % Chuyen sang mien tan so
    NFFT = 2048;
    NOVERLAP = floor(NFFT * 0.75);
    WINDOW = hamming(NFFT);
    
    [S1, F, T] = spectrogram(nguoi_1, WINDOW, NOVERLAP, NFFT, fs_target);
    [S2, ~, ~] = spectrogram(nguoi_2, WINDOW, NOVERLAP, NFFT, fs_target);
    
    % Tinh cong suat
    P1 = abs(S1).^2;
    P2 = abs(S2).^2;
    P_total = P1 + P2 + eps;
    
    % Tao Wiener mask
    mask1 = P1 ./ P_total;
    mask2 = P2 ./ P_total;
    
    % Lam mem mask (soft masking)
    mask1 = mask1.^2;  % Binary masking: mask^2
    mask2 = mask2.^2;
    
    % Ap dung mask
    S1_clean = S1 .* mask1;
    S2_clean = S2 .* mask2;
    
    % Chuyen ve mien thoi gian
    nguoi_1_clean = my_istft(S1_clean, WINDOW, NOVERLAP, NFFT, min_len);
    nguoi_2_clean = my_istft(S2_clean, WINDOW, NOVERLAP, NFFT, min_len);
    
    nguoi_1_clean = nguoi_1_clean(1:min_len);
    nguoi_2_clean = nguoi_2_clean(1:min_len);
    
    % Kiem tra correlation
    corr_after = abs(corrcoef(nguoi_1_clean, nguoi_2_clean));
    if ~isnan(corr_after(1, 2))
        disp(['  + Correlation sau Wiener masking: ' num2str(corr_after(1, 2), '%.3f')]);
        
        % Chi ap dung neu cai thien
        if corr_after(1, 2) < best_corr
            nguoi_1 = nguoi_1_clean;
            nguoi_2 = nguoi_2_clean;
            disp('  + Da ap dung Wiener masking');
        else
            disp('  + Giu nguyen ket qua goc');
        end
    end
else
    disp('  + Correlation da thap, khong can Wiener masking');
end

disp('');

%% 9. HAU XU LY CUOI CUNG
disp('BUOC 6: HAU XU LY CUOI CUNG');
disp('-----------------------------------------------------------');

% Voice Activity Detection (VAD)
for i = 1:2
    if i == 1
        sig = nguoi_1;
    else
        sig = nguoi_2;
    end
    
    % Tinh energy theo frame
    frame_size = round(0.025 * fs_target);
    hop_size = round(0.010 * fs_target);
    num_frames = floor((length(sig) - frame_size) / hop_size) + 1;
    
    frame_energy = zeros(1, num_frames);
    for f = 1:num_frames
        start_idx = (f - 1) * hop_size + 1;
        end_idx = min(start_idx + frame_size - 1, length(sig));
        frame_energy(f) = sum(sig(start_idx:end_idx).^2);
    end
    
    % Nguong VAD adaptive
    noise_level = prctile(frame_energy, 20);
    vad_threshold = noise_level * 2.5;
    
    vad_mask = frame_energy > vad_threshold;
    
    % Mo rong mask
    vad_mask_expanded = vad_mask;
    for f = 2:num_frames-1
        if vad_mask(f-1) || vad_mask(f+1)
            vad_mask_expanded(f) = 1;
        end
    end
    
    % Ap dung VAD mask len tin hieu
    sig_vad = zeros(size(sig));
    for f = 1:num_frames
        start_idx = (f - 1) * hop_size + 1;
        end_idx = min(start_idx + hop_size - 1, length(sig));
        if vad_mask_expanded(f)
            sig_vad(start_idx:end_idx) = sig(start_idx:end_idx);
        end
    end
    
    % Smooth transition
    sig_vad = filter(ones(1, 30)/30, 1, sig_vad);
    
    if i == 1
        nguoi_1 = sig_vad;
    else
        nguoi_2 = sig_vad;
    end
end
disp('  + Voice Activity Detection (VAD)');

% Chuan hoa cuoi cung
nguoi_1 = nguoi_1 / (max(abs(nguoi_1)) + eps);
nguoi_2 = nguoi_2 / (max(abs(nguoi_2)) + eps);
disp('  + Chuan hoa cuoi cung');
disp('');

%% 10. LOAI BO PHAN DINH BANG TRU TRUC TIEP
disp('BUOC 7: LOAI BO PHAN DINH BANG TRU TRUC TIEP');
disp('-----------------------------------------------------------');

corr_check = abs(corrcoef(nguoi_1, nguoi_2));
if ~isnan(corr_check(1, 2)) && corr_check(1, 2) > 0.3
    disp(['  + Correlation truoc khi tru: ' num2str(corr_check(1, 2), '%.3f')]);
    
    % Phuong phap orthogonalization: lam cho 2 nguon vuong goc
    % nguoi_1_new = nguoi_1
    % nguoi_2_new = nguoi_2 - projection(nguoi_2 onto nguoi_1)
    
    projection_coef = (nguoi_2 * nguoi_1') / (nguoi_1 * nguoi_1' + eps);
    nguoi_2_ortho = nguoi_2 - projection_coef * nguoi_1;
    
    projection_coef2 = (nguoi_1 * nguoi_2_ortho') / (nguoi_2_ortho * nguoi_2_ortho' + eps);
    nguoi_1_ortho = nguoi_1 - projection_coef2 * nguoi_2_ortho;
    
    % Chuan hoa
    nguoi_1_ortho = nguoi_1_ortho / (max(abs(nguoi_1_ortho)) + eps);
    nguoi_2_ortho = nguoi_2_ortho / (max(abs(nguoi_2_ortho)) + eps);
    
    % Kiem tra correlation
    corr_ortho = abs(corrcoef(nguoi_1_ortho, nguoi_2_ortho));
    if ~isnan(corr_ortho(1, 2))
        disp(['  + Correlation sau khi orthogonalization: ' num2str(corr_ortho(1, 2), '%.3f')]);
        
        % Chi ap dung neu cai thien
        if corr_ortho(1, 2) < corr_check(1, 2) * 0.8
            nguoi_1 = nguoi_1_ortho;
            nguoi_2 = nguoi_2_ortho;
            disp('  + Da ap dung orthogonalization');
        else
            disp('  + Giu nguyen ket qua');
        end
    end
else
    disp('  + Correlation da thap, khong can xu ly them');
end

disp('');

%% 11. LUU KET QUA
disp('BUOC 8: LUU KET QUA');
disp('-----------------------------------------------------------');

output_file1 = fullfile(result_folder, 'nguoi_1.wav');
output_file2 = fullfile(result_folder, 'nguoi_2.wav');

audiowrite(output_file1, nguoi_1', fs_target);
audiowrite(output_file2, nguoi_2', fs_target);

disp(['  + Luu: ' output_file1]);
disp(['  + Luu: ' output_file2]);
disp('');

%% 12. HIEN THI KET QUA
disp('BUOC 9: HIEN THI KET QUA');
disp('-----------------------------------------------------------');

t = (0:min_len-1) / fs_target;

figure('Name', 'Ket qua tach 2 nguoi', 'Position', [100, 100, 1600, 900]);

% Row 1: Input
subplot(3, 3, 1);
plot(t, X(1, :));
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('DAU VAO: Microphone 1');
grid on; xlim([0, max(t)]);

subplot(3, 3, 2);
plot(t, X(2, :));
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('DAU VAO: Microphone 2');
grid on; xlim([0, max(t)]);

subplot(3, 3, 3);
[S_in, F_in, T_in] = spectrogram(X(1, :), WINDOW, NOVERLAP, NFFT, fs_target);
imagesc(T_in, F_in, 20*log10(abs(S_in) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Input');
ylim([0, 4000]);

% Row 2: Nguoi 1
subplot(3, 3, 4);
plot(t, nguoi_1);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('KET QUA: Nguoi 1');
grid on; xlim([0, max(t)]);
set(gca, 'Color', [1 0.95 0.95]);  % Background mau nhe

subplot(3, 3, 5);
[S1_final, F1_final, T1_final] = spectrogram(nguoi_1, WINDOW, NOVERLAP, NFFT, fs_target);
imagesc(T1_final, F1_final, 20*log10(abs(S1_final) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Nguoi 1');
ylim([0, 4000]);

subplot(3, 3, 6);
Pxx1 = mean(abs(S1_final).^2, 2);
plot(F1_final, 10*log10(Pxx1 + eps), 'LineWidth', 2);
xlabel('Tan so (Hz)'); ylabel('Cong suat (dB)');
title('Pho tan so - Nguoi 1');
grid on; xlim([0, 4000]);

% Row 3: Nguoi 2
subplot(3, 3, 7);
plot(t, nguoi_2);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('KET QUA: Nguoi 2');
grid on; xlim([0, max(t)]);
set(gca, 'Color', [0.95 0.95 1]);  % Background mau nhe

subplot(3, 3, 8);
[S2_final, F2_final, T2_final] = spectrogram(nguoi_2, WINDOW, NOVERLAP, NFFT, fs_target);
imagesc(T2_final, F2_final, 20*log10(abs(S2_final) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Nguoi 2');
ylim([0, 4000]);

subplot(3, 3, 9);
Pxx2 = mean(abs(S2_final).^2, 2);
plot(F2_final, 10*log10(Pxx2 + eps), 'LineWidth', 2);
xlabel('Tan so (Hz)'); ylabel('Cong suat (dB)');
title('Pho tan so - Nguoi 2');
grid on; xlim([0, 4000]);

savefig(fullfile(result_folder, 'ket_qua.fig'));
disp('  + Da hien thi ket qua');
disp('');

%% 13. DANH GIA CUOI CUNG
disp('BUOC 10: DANH GIA CUOI CUNG');
disp('-----------------------------------------------------------');

energy1 = sum(nguoi_1.^2);
energy2 = sum(nguoi_2.^2);
corr_final = abs(corrcoef(nguoi_1, nguoi_2));

disp(['  + Nang luong nguoi 1: ' num2str(energy1, '%.2e')]);
disp(['  + Nang luong nguoi 2: ' num2str(energy2, '%.2e')]);

if ~isnan(corr_final(1, 2))
    disp(['  + Correlation cuoi cung: ' num2str(corr_final(1, 2), '%.3f')]);
    
    if corr_final(1, 2) < 0.2
        disp('  ');
        disp('  *** KET QUA XUAT SAC: 2 nguoi da duoc tach ro rang! ***');
    elseif corr_final(1, 2) < 0.4
        disp('  ');
        disp('  *** KET QUA TOT: 2 nguoi da duoc tach, co the con it dinh. ***');
    elseif corr_final(1, 2) < 0.6
        disp('  ');
        disp('  *** KET QUA KHAM: 2 nguoi da duoc tach nhung con dinh nhieu. ***');
    else
        disp('  ');
        disp('  *** CANH BAO: Ket qua chua tot, 2 nguoi con dinh nhieu! ***');
        disp('  *** Goi y: Thu dieu chinh tham so hoac dung phuong phap khac. ***');
    end
end

disp('');
disp('=============================================================');
disp('                      HOAN THANH                             ');
disp('=============================================================');
disp('');
disp(['Ket qua: ' result_folder]);
disp(['  - nguoi_1.wav']);
disp(['  - nguoi_2.wav']);
disp(['  - ket_qua.fig']);
disp('');
disp('HAY NGHE THU 2 FILE VA XEM FIGURE DE KIEM TRA!');
disp('');
