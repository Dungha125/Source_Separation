%% integrated_separation.m
% THUAT TOAN TICH HOP: Tat ca cac phuong phap LONG VAO NHAU
% Khong phai chon, ma la KET HOP trong cung 1 quy trinh
clear all; close all; clc;

disp('=============================================================');
disp('  THUAT TOAN TICH HOP: ICA + BEAMFORMING + CLUSTERING       ');
disp('=============================================================');
disp('');

%% 1. TAI VA TIEN XU LY
disp('BUOC 1: TAI VA TIEN XU LY');
disp('-----------------------------------------------------------');

fs_target = 16000;
result_folder = 'output_integrated';
if ~exist(result_folder, 'dir'), mkdir(result_folder); end

mic_files = {'mic1.wav', 'mic2.wav'};
X_data = cell(2, 1);

for i = 1:2
    if ~exist(mic_files{i}, 'file'), error(['Khong tim thay: ' mic_files{i}]); end
    [x, fs_in] = audioread(mic_files{i});
    if fs_in ~= fs_target, x = resample(x, fs_target, fs_in); end
    if size(x, 2) > 1, x = mean(x, 2); end
    X_data{i} = x(:)';
end

min_len = min(cellfun(@length, X_data));
X = zeros(2, min_len);
X(1, :) = X_data{1}(1:min_len);
X(2, :) = X_data{2}(1:min_len);

alpha = 0.97;
X_pre = zeros(size(X));
for i = 1:2
    X_pre(i, :) = filter([1, -alpha], 1, X(i, :));
end

[b_bp, a_bp] = butter(6, [250, 4000]/(fs_target/2), 'bandpass');
X_bp = zeros(size(X_pre));
for i = 1:2
    X_bp(i, :) = filtfilt(b_bp, a_bp, X_pre(i, :));
    X_bp(i, :) = X_bp(i, :) / (max(abs(X_bp(i, :))) + eps);
end

disp(['  + Tien xu ly, do dai: ' num2str(min_len/fs_target, '%.2f') 's']);
disp('');

%% 2. PHAN TICH MIEN TAN SO
disp('BUOC 2: PHAN TICH TIME-FREQUENCY');
disp('-----------------------------------------------------------');

NFFT = 2048;
NOVERLAP = floor(NFFT * 0.75);
WINDOW = hamming(NFFT);

[S1_tf, F, T] = spectrogram(X_bp(1,:), WINDOW, NOVERLAP, NFFT, fs_target);
[S2_tf, ~, ~] = spectrogram(X_bp(2,:), WINDOW, NOVERLAP, NFFT, fs_target);

[num_freq, num_frame] = size(S1_tf);
disp(['  + STFT: ' num2str(num_freq) ' x ' num2str(num_frame)]);
disp('');

%% 3. BUOC 1: ICA TRONG MIEN THOI GIAN -> UOC TINH MIXING MATRIX
disp('BUOC 3: ICA - UOC TINH MIXING MATRIX');
disp('-----------------------------------------------------------');

[S_ica, ~, W_ica] = fastica_robust(X_bp);

if size(S_ica, 1) < 2, error('ICA chi tim thay 1 nguon!'); end

% Uoc tinh mixing matrix: A = W^(-1)
A_ica = inv(W_ica);

disp('  + Ma tran tron A (tu ICA):');
disp(['    A = [' num2str(A_ica(1,1), '%.3f') ', ' num2str(A_ica(1,2), '%.3f') ']']);
disp(['        [' num2str(A_ica(2,1), '%.3f') ', ' num2str(A_ica(2,2), '%.3f') ']']);
disp('  + Y nghia: X = A * S (Input = Ma tran tron * Nguon)');
disp('');

%% 4. BUOC 2: CLUSTERING IPD/ILD -> UOC TINH HUONG
disp('BUOC 4: CLUSTERING IPD/ILD - UOC TINH HUONG');
disp('-----------------------------------------------------------');

IPD = angle(S2_tf ./ (S1_tf + eps));
ILD = 20*log10((abs(S2_tf) + eps) ./ (abs(S1_tf) + eps));

mag_sum = abs(S1_tf) + abs(S2_tf);
threshold = prctile(mag_sum(:), 70);
mask_voice = mag_sum > threshold;

features = [IPD(mask_voice), ILD(mask_voice)/10];
[idx_cluster, C] = kmeans(features, 2, 'Replicates', 10, 'MaxIter', 1000);

% Uoc tinh goc tu IPD
d = 0.05; c = 343;
f_mid = 1000;

theta1_rad = asin(max(min(C(1,1) * c / (2*pi*f_mid*d), 1), -1));
theta2_rad = asin(max(min(C(2,1) * c / (2*pi*f_mid*d), 1), -1));

theta1_deg = rad2deg(theta1_rad);
theta2_deg = rad2deg(theta2_rad);

disp(['  + Clustering: 2 nhom tu IPD/ILD']);
disp(['  + Uoc tinh huong 1: ' num2str(theta1_deg, '%.1f') ' do']);
disp(['  + Uoc tinh huong 2: ' num2str(theta2_deg, '%.1f') ' do']);
disp('');

%% 5. BUOC 3: TAO SPATIAL WEIGHT KET HOP BEAMFORMER
disp('BUOC 5: TAO SPATIAL WEIGHT (KET HOP BEAMFORMER)');
disp('-----------------------------------------------------------');
disp('  => Long ghep: DAS + GSC + Differential');
disp('');

% 5.1. Delay-and-Sum Beamformer weight
disp('  5.1. Delay-and-Sum Beamforming...');
d_mic = 0.05;

try
    y1_das = delay_and_sum_beamformer(X_bp, fs_target, theta1_deg);
    y2_das = delay_and_sum_beamformer(X_bp, fs_target, theta2_deg);
    
    [S1_das_tf, ~, ~] = spectrogram(y1_das, WINDOW, NOVERLAP, NFFT, fs_target);
    [S2_das_tf, ~, ~] = spectrogram(y2_das, WINDOW, NOVERLAP, NFFT, fs_target);
    
    P1_das = abs(S1_das_tf).^2;
    P2_das = abs(S2_das_tf).^2;
    P_total_das = P1_das + P2_das + eps;
    
    weight_das_1 = P1_das ./ P_total_das;
    weight_das_2 = P2_das ./ P_total_das;
    disp('      + DAS: OK');
catch ME
    weight_das_1 = 0.5 * ones(num_freq, num_frame);
    weight_das_2 = 0.5 * ones(num_freq, num_frame);
    disp(['      ! DAS loi: ' ME.message]);
end

% 5.2. GSC Beamformer weight
disp('  5.2. GSC Beamforming...');
try
    [y1_gsc, ~] = gsc_beamformer(X_bp, fs_target, theta1_deg, 256, 0.003);
    [y2_gsc, ~] = gsc_beamformer(X_bp, fs_target, theta2_deg, 256, 0.003);
    
    [S1_gsc_tf, ~, ~] = spectrogram(y1_gsc, WINDOW, NOVERLAP, NFFT, fs_target);
    [S2_gsc_tf, ~, ~] = spectrogram(y2_gsc, WINDOW, NOVERLAP, NFFT, fs_target);
    
    P1_gsc = abs(S1_gsc_tf).^2;
    P2_gsc = abs(S2_gsc_tf).^2;
    P_total_gsc = P1_gsc + P2_gsc + eps;
    
    weight_gsc_1 = P1_gsc ./ P_total_gsc;
    weight_gsc_2 = P2_gsc ./ P_total_gsc;
    disp('      + GSC: OK');
catch ME
    weight_gsc_1 = 0.5 * ones(num_freq, num_frame);
    weight_gsc_2 = 0.5 * ones(num_freq, num_frame);
    disp(['      ! GSC loi: ' ME.message]);
end

% 5.3. Differential Microphone Array weight
disp('  5.3. Differential Array...');
try
    y_diff = differential_microphone_array(X_bp, 1);
    y_diff_inv = -y_diff;
    
    [S_diff_tf, ~, ~] = spectrogram(y_diff, WINDOW, NOVERLAP, NFFT, fs_target);
    [S_diff_inv_tf, ~, ~] = spectrogram(y_diff_inv, WINDOW, NOVERLAP, NFFT, fs_target);
    
    P_diff = abs(S_diff_tf).^2;
    P_diff_inv = abs(S_diff_inv_tf).^2;
    P_total_diff = P_diff + P_diff_inv + eps;
    
    weight_diff_1 = P_diff ./ P_total_diff;
    weight_diff_2 = P_diff_inv ./ P_total_diff;
    disp('      + Differential: OK');
catch ME
    weight_diff_1 = 0.5 * ones(num_freq, num_frame);
    weight_diff_2 = 0.5 * ones(num_freq, num_frame);
    disp(['      ! Differential loi: ' ME.message]);
end

% KET HOP CAC BEAMFORMER WEIGHT (LONG GHEP)
disp('');
disp('  => Ket hop beamformer weights:');

% He so cho tung beamformer (tong = 1.0)
w_das = 0.40;  % Delay-and-Sum
w_gsc = 0.40;  % GSC
w_diff = 0.20; % Differential

spatial_weight_1 = w_das * weight_das_1 + w_gsc * weight_gsc_1 + w_diff * weight_diff_1;
                   
spatial_weight_2 = w_das * weight_das_2 + w_gsc * weight_gsc_2 + w_diff * weight_diff_2;

% Normalize
spatial_sum = spatial_weight_1 + spatial_weight_2 + eps;
spatial_weight_1 = spatial_weight_1 ./ spatial_sum;
spatial_weight_2 = spatial_weight_2 ./ spatial_sum;

disp(['      DAS=' num2str(w_das) ' + GSC=' num2str(w_gsc) ' + Diff=' num2str(w_diff)]);
disp('  + Da ket hop beamformer weights');
disp('');

%% 6. BUOC 4: TAO ICA-BASED MASK
disp('BUOC 6: TAO ICA-BASED MASK');
disp('-----------------------------------------------------------');

% Chuyen ICA sources sang TF domain
S1_ica = S_ica(1, :) / (max(abs(S_ica(1, :))) + eps);
S2_ica = S_ica(2, :) / (max(abs(S_ica(2, :))) + eps);

[S1_ica_tf, ~, ~] = spectrogram(S1_ica, WINDOW, NOVERLAP, NFFT, fs_target);
[S2_ica_tf, ~, ~] = spectrogram(S2_ica, WINDOW, NOVERLAP, NFFT, fs_target);

% Wiener mask tu ICA
P1_ica = abs(S1_ica_tf).^2;
P2_ica = abs(S2_ica_tf).^2;
P_total_ica = P1_ica + P2_ica + eps;

mask_ica_1 = P1_ica ./ P_total_ica;
mask_ica_2 = P2_ica ./ P_total_ica;

disp('  + Da tao ICA-based Wiener mask');
disp('');

%% 7. BUOC 5: TAO CLUSTERING-BASED MASK
disp('BUOC 7: TAO CLUSTERING-BASED MASK');
disp('-----------------------------------------------------------');

mask_cluster_1 = zeros(num_freq, num_frame);
mask_cluster_2 = zeros(num_freq, num_frame);

for f = 1:num_freq
    for t = 1:num_frame
        dist1 = sqrt((IPD(f,t) - C(1,1))^2 + ((ILD(f,t)/10) - C(1,2))^2);
        dist2 = sqrt((IPD(f,t) - C(2,1))^2 + ((ILD(f,t)/10) - C(2,2))^2);
        
        w1 = exp(-dist1^2 / 0.5);
        w2 = exp(-dist2^2 / 0.5);
        
        total = w1 + w2 + eps;
        mask_cluster_1(f, t) = w1 / total;
        mask_cluster_2(f, t) = w2 / total;
    end
end

disp('  + Da tao Clustering-based mask tu IPD/ILD');
disp('');

%% 8. BUOC 6: KET HOP TAT CA CAC MASK (CORE - LONG VAO NHAU)
disp('BUOC 8: KET HOP TAT CA CAC MASK');
disp('-----------------------------------------------------------');
disp('  => KET HOP: ICA mask + Spatial weight + Clustering mask');
disp('');

% He so ket hop (co the dieu chinh)
w_ica = 0.4;        % Trong so ICA
w_spatial = 0.3;    % Trong so Spatial (tu beamforming)
w_cluster = 0.3;    % Trong so Clustering

% KET HOP: Average가重 cua 3 loai mask
mask_integrated_1 = w_ica * mask_ica_1 + w_spatial * spatial_weight_1 + w_cluster * mask_cluster_1;
mask_integrated_2 = w_ica * mask_ica_2 + w_spatial * spatial_weight_2 + w_cluster * mask_cluster_2;

% Normalize
mask_sum = mask_integrated_1 + mask_integrated_2 + eps;
mask_integrated_1 = mask_integrated_1 ./ mask_sum;
mask_integrated_2 = mask_integrated_2 ./ mask_sum;

disp(['  + He so ket hop: ICA=' num2str(w_ica) ', Spatial=' num2str(w_spatial) ', Cluster=' num2str(w_cluster)]);
disp('  + Da tao integrated mask');
disp('');

%% 9. ITERATIVE REFINEMENT voi INTEGRATED MASK
disp('BUOC 9: ITERATIVE REFINEMENT');
disp('-----------------------------------------------------------');

S1_current = S1_tf .* mask_integrated_1;
S2_current = S2_tf .* mask_integrated_2;

for iter = 1:5
    % Cap nhat Wiener mask tu ket qua hien tai
    P1_current = abs(S1_current).^2;
    P2_current = abs(S2_current).^2;
    P_total_current = P1_current + P2_current + eps;
    
    mask_wiener_new_1 = P1_current ./ P_total_current;
    mask_wiener_new_2 = P2_current ./ P_total_current;
    
    % KET HOP LAI voi spatial va clustering (LONG VAO)
    mask_iter_1 = 0.5 * mask_wiener_new_1 + 0.3 * spatial_weight_1 + 0.2 * mask_cluster_1;
    mask_iter_2 = 0.5 * mask_wiener_new_2 + 0.3 * spatial_weight_2 + 0.2 * mask_cluster_2;
    
    % Binary masking
    mask_alpha = 2.0 + iter * 0.2;  % Tang dan (binary hon theo iteration)
    mask_iter_1 = mask_iter_1.^mask_alpha;
    mask_iter_2 = mask_iter_2.^mask_alpha;
    
    % Normalize
    mask_sum_iter = mask_iter_1 + mask_iter_2 + eps;
    mask_iter_1 = mask_iter_1 ./ mask_sum_iter;
    mask_iter_2 = mask_iter_2 ./ mask_sum_iter;
    
    % Ap dung len SPECTROGRAM GOC (quan trong!)
    S1_current = S1_tf .* mask_iter_1;
    S2_current = S2_tf .* mask_iter_2;
    
    % Kiem tra
    y1_check = my_istft(S1_current, WINDOW, NOVERLAP, NFFT, min_len);
    y2_check = my_istft(S2_current, WINDOW, NOVERLAP, NFFT, min_len);
    len_check = min([length(y1_check), length(y2_check), min_len]);
    
    if len_check > 100
        corr_iter = corrcoef(y1_check(1:len_check), y2_check(1:len_check));
        if ~isnan(corr_iter(1, 2))
            disp(['  Iteration ' num2str(iter) ': Correlation = ' num2str(abs(corr_iter(1, 2)), '%.3f') ', mask_alpha = ' num2str(mask_alpha, '%.1f')]);
        end
    end
end

disp('  + Hoan thanh iterative refinement');
disp('');

%% 10. CHUYEN VE MIEN THOI GIAN
disp('BUOC 10: ISTFT - CHUYEN VE MIEN THOI GIAN');
disp('-----------------------------------------------------------');

nguoi_1_temp = my_istft(S1_current, WINDOW, NOVERLAP, NFFT, min_len);
nguoi_2_temp = my_istft(S2_current, WINDOW, NOVERLAP, NFFT, min_len);

len1 = length(nguoi_1_temp);
len2 = length(nguoi_2_temp);

if len1 < min_len
    nguoi_1 = [nguoi_1_temp, zeros(1, min_len - len1)];
else
    nguoi_1 = nguoi_1_temp(1:min_len);
end

if len2 < min_len
    nguoi_2 = [nguoi_2_temp, zeros(1, min_len - len2)];
else
    nguoi_2 = nguoi_2_temp(1:min_len);
end

disp('  + Da chuyen ve mien thoi gian');
disp('');

%% 11. HAU XU LY
disp('BUOC 11: HAU XU LY');
disp('-----------------------------------------------------------');

% De-emphasis
nguoi_1 = filter(1, [1, -alpha], nguoi_1);
nguoi_2 = filter(1, [1, -alpha], nguoi_2);

nguoi_1 = nguoi_1 / (max(abs(nguoi_1)) + eps);
nguoi_2 = nguoi_2 / (max(abs(nguoi_2)) + eps);

% Orthogonalization manh (lap 5 lan)
for orth_iter = 1:5
    proj = (nguoi_2 * nguoi_1') / (nguoi_1 * nguoi_1' + eps);
    nguoi_2 = nguoi_2 - 0.85 * proj * nguoi_1;
    
    proj2 = (nguoi_1 * nguoi_2') / (nguoi_2 * nguoi_2' + eps);
    nguoi_1 = nguoi_1 - 0.85 * proj2 * nguoi_2;
    
    nguoi_1 = nguoi_1 / (max(abs(nguoi_1)) + eps);
    nguoi_2 = nguoi_2 / (max(abs(nguoi_2)) + eps);
end

% VAD
for i = 1:2
    if i == 1, sig = nguoi_1; else, sig = nguoi_2; end
    
    frame_size = round(0.025 * fs_target);
    hop_size = round(0.010 * fs_target);
    num_frames = floor((length(sig) - frame_size) / hop_size) + 1;
    
    frame_energy = zeros(1, num_frames);
    for f = 1:num_frames
        start_idx = (f - 1) * hop_size + 1;
        end_idx = min(start_idx + frame_size - 1, length(sig));
        frame_energy(f) = sum(sig(start_idx:end_idx).^2);
    end
    
    vad_threshold = prctile(frame_energy, 15) * 3;
    vad_mask = frame_energy > vad_threshold;
    vad_mask_ext = conv(double(vad_mask), ones(1, 5), 'same') > 0;
    
    sig_vad = zeros(size(sig));
    for f = 1:num_frames
        start_idx = (f - 1) * hop_size + 1;
        end_idx = min(start_idx + hop_size - 1, length(sig));
        if vad_mask_ext(f)
            sig_vad(start_idx:end_idx) = sig(start_idx:end_idx);
        end
    end
    
    sig_vad = filter(ones(1, 50)/50, 1, sig_vad);
    
    if i == 1, nguoi_1 = sig_vad; else, nguoi_2 = sig_vad; end
end

nguoi_1 = nguoi_1 / (max(abs(nguoi_1)) + eps);
nguoi_2 = nguoi_2 / (max(abs(nguoi_2)) + eps);

disp('  + De-emphasis, Orthogonalization x5, VAD');
disp('');

%% 12. LUU KET QUA
disp('BUOC 12: LUU KET QUA');
disp('-----------------------------------------------------------');

audiowrite(fullfile(result_folder, 'nguoi_1.wav'), nguoi_1', fs_target);
audiowrite(fullfile(result_folder, 'nguoi_2.wav'), nguoi_2', fs_target);

disp('  + Luu: nguoi_1.wav');
disp('  + Luu: nguoi_2.wav');
disp('');

%% 13. HIEN THI
disp('BUOC 13: HIEN THI');
disp('-----------------------------------------------------------');

t = (0:min_len-1) / fs_target;

fig = figure('Name', 'Integrated Algorithm', 'Position', [50, 50, 1800, 1000]);

% Row 1: Input va Clustering
subplot(4, 4, 1);
plot(t, X_bp(1, :), 'k');
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Input: Mic 1');
grid on; xlim([0, max(t)]);

subplot(4, 4, 2);
plot(t, X_bp(2, :), 'k');
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Input: Mic 2');
grid on; xlim([0, max(t)]);

subplot(4, 4, 3);
scatter(features(:,1), features(:,2)*10, 3, idx_cluster, 'filled');
xlabel('IPD'); ylabel('ILD');
title('Clustering IPD/ILD');
hold on;
plot(C(1,1), C(1,2)*10, 'rx', 'MarkerSize', 15, 'LineWidth', 3);
plot(C(2,1), C(2,2)*10, 'bx', 'MarkerSize', 15, 'LineWidth', 3);
legend('Data', 'Center 1', 'Center 2');
grid on;

subplot(4, 4, 4);
imagesc(T, F, 20*log10(abs(S1_tf) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram Input');
ylim([0, 4000]);

% Row 2: Masks
subplot(4, 4, 5);
imagesc(mask_ica_1);
axis xy; colorbar; caxis([0 1]);
title('ICA Mask - Nguon 1');
xlabel('Frame'); ylabel('Tan so');

subplot(4, 4, 6);
imagesc(spatial_weight_1);
axis xy; colorbar; caxis([0 1]);
title('Spatial Weight - Nguon 1');
xlabel('Frame'); ylabel('Tan so');

subplot(4, 4, 7);
imagesc(mask_cluster_1);
axis xy; colorbar; caxis([0 1]);
title('Clustering Mask - Nguon 1');
xlabel('Frame'); ylabel('Tan so');

subplot(4, 4, 8);
imagesc(mask_iter_1);
axis xy; colorbar; caxis([0 1]);
title('INTEGRATED Mask - Nguon 1');
xlabel('Frame'); ylabel('Tan so');

% Row 3: Nguoi 1
subplot(4, 4, 9);
plot(t, nguoi_1, 'r', 'LineWidth', 1.5);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('OUTPUT: Nguoi 1');
grid on; xlim([0, max(t)]);
set(gca, 'Color', [1 0.95 0.95]);

subplot(4, 4, 10);
[S1_out, F1_out, T1_out] = spectrogram(nguoi_1, WINDOW, NOVERLAP, NFFT, fs_target);
imagesc(T1_out, F1_out, 20*log10(abs(S1_out) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Nguoi 1');
ylim([0, 4000]);

subplot(4, 4, 11);
Pxx1 = mean(abs(S1_out).^2, 2);
plot(F1_out, 10*log10(Pxx1 + eps), 'r', 'LineWidth', 2);
xlabel('Tan so (Hz)'); ylabel('Cong suat (dB)');
title('Pho tan so - Nguoi 1');
grid on; xlim([0, 4000]);

subplot(4, 4, 12);
plot(t, S1_ica, 'Color', [0.7 0.7 0.7]);
hold on;
plot(t, nguoi_1, 'r', 'LineWidth', 1.5);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('So sanh: ICA goc (xam) vs Ket qua (do)');
legend('ICA goc', 'Ket qua');
grid on; xlim([0, max(t)]);

% Row 4: Nguoi 2
subplot(4, 4, 13);
plot(t, nguoi_2, 'b', 'LineWidth', 1.5);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('OUTPUT: Nguoi 2');
grid on; xlim([0, max(t)]);
set(gca, 'Color', [0.95 0.95 1]);

subplot(4, 4, 14);
[S2_out, F2_out, T2_out] = spectrogram(nguoi_2, WINDOW, NOVERLAP, NFFT, fs_target);
imagesc(T2_out, F2_out, 20*log10(abs(S2_out) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Nguoi 2');
ylim([0, 4000]);

subplot(4, 4, 15);
Pxx2 = mean(abs(S2_out).^2, 2);
plot(F2_out, 10*log10(Pxx2 + eps), 'b', 'LineWidth', 2);
xlabel('Thoi gian (Hz)'); ylabel('Cong suat (dB)');
title('Pho tan so - Nguoi 2');
grid on; xlim([0, 4000]);

subplot(4, 4, 16);
plot(t, S2_ica, 'Color', [0.7 0.7 0.7]);
hold on;
plot(t, nguoi_2, 'b', 'LineWidth', 1.5);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('So sanh: ICA goc (xam) vs Ket qua (xanh)');
legend('ICA goc', 'Ket qua');
grid on; xlim([0, max(t)]);

savefig(fullfile(result_folder, 'ket_qua_integrated.fig'));
disp('  + Da hien thi');
disp('');

%% 14. DANH GIA
disp('BUOC 14: DANH GIA');
disp('-----------------------------------------------------------');

corr_final = corrcoef(nguoi_1, nguoi_2);
if ~isnan(corr_final(1, 2))
    disp(['  + Correlation giua 2 nguoi: ' num2str(abs(corr_final(1, 2)), '%.3f')]);
end

corr_in1 = corrcoef(nguoi_1, X_bp(1, :));
corr_in2 = corrcoef(nguoi_2, X_bp(2, :));

if ~isnan(corr_in1(1, 2))
    disp(['  + Corr nguoi 1 vs mic 1: ' num2str(abs(corr_in1(1, 2)), '%.3f')]);
end
if ~isnan(corr_in2(1, 2))
    disp(['  + Corr nguoi 2 vs mic 2: ' num2str(abs(corr_in2(1, 2)), '%.3f')]);
end



