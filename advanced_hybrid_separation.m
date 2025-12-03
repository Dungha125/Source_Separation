%% advanced_hybrid_separation.m
% THUAT TOAN LAI NANG CAO: Ket hop ICA + Tat ca cac phuong phap Beamforming
% Muc tieu: Tach duoc tieng 2 nguoi rieng biet, khong phai ban sao cua input
clear all; close all; clc;

disp('=============================================================');
disp('  THUAT TOAN LAI NANG CAO: ICA + MULTI-BEAMFORMING          ');
disp('=============================================================');
disp('');

%% 1. TAI VA TIEN XU LY
disp('BUOC 1: TAI VA TIEN XU LY');
disp('-----------------------------------------------------------');

fs_target = 16000;
result_folder = 'output_advanced_hybrid';
if ~exist(result_folder, 'dir'), mkdir(result_folder); end

mic_files = {'mic1.wav', 'mic2.wav'};
X_data = cell(2, 1);

for i = 1:2
    if ~exist(mic_files{i}, 'file')
        error(['Khong tim thay file: ' mic_files{i}]);
    end
    
    [x, fs_in] = audioread(mic_files{i});
    if fs_in ~= fs_target, x = resample(x, fs_target, fs_in); end
    if size(x, 2) > 1, x = mean(x, 2); end
    X_data{i} = x(:)';
end

min_len = min(cellfun(@length, X_data));
X = zeros(2, min_len);
X(1, :) = X_data{1}(1:min_len);
X(2, :) = X_data{2}(1:min_len);

% Tien xu ly
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

disp(['  + Da tien xu ly, do dai: ' num2str(min_len/fs_target, '%.2f') ' giay']);
disp('');

%% 2. CHUYEN SANG MIEN TAN SO
disp('BUOC 2: PHAN TICH MIEN TAN SO');
disp('-----------------------------------------------------------');

NFFT = 2048;
NOVERLAP = floor(NFFT * 0.75);
WINDOW = hamming(NFFT);

[S1_tf, F, T] = spectrogram(X_bp(1,:), WINDOW, NOVERLAP, NFFT, fs_target);
[S2_tf, ~, ~] = spectrogram(X_bp(2,:), WINDOW, NOVERLAP, NFFT, fs_target);

disp(['  + STFT: ' num2str(size(S1_tf, 1)) ' tan so x ' num2str(size(S1_tf, 2)) ' frame']);

% Tinh IPD/ILD
IPD = angle(S2_tf ./ (S1_tf + eps));
ILD = 20*log10((abs(S2_tf) + eps) ./ (abs(S1_tf) + eps));

disp('  + Da tinh IPD/ILD');
disp('');

%% 3. UOC TINH HUONG CUA 2 NGUON (Tu IPD/ILD clustering)
disp('BUOC 3: UOC TINH HUONG CUA 2 NGUON');
disp('-----------------------------------------------------------');

mag_sum = abs(S1_tf) + abs(S2_tf);
threshold = prctile(mag_sum(:), 75);
mask_active = mag_sum > threshold;

features = [IPD(mask_active), ILD(mask_active)/10];
[idx_cluster, C] = kmeans(features, 2, 'Replicates', 10, 'MaxIter', 1000);

% Tu center IPD, uoc tinh goc
% IPD ≈ (2*pi*f*d*sin(theta))/c
% => sin(theta) ≈ IPD * c / (2*pi*f*d)

d = 0.05;  % khoang cach mic (5cm)
c = 343;   % toc do am thanh

% Dung tan so trung binh (500-2000 Hz cho giong noi)
f_mid = 1000;  % Hz
theta1_est = asin(C(1,1) * c / (2*pi*f_mid*d));
theta2_est = asin(C(2,1) * c / (2*pi*f_mid*d));

% Chuyen sang do
theta1_deg = rad2deg(theta1_est);
theta2_deg = rad2deg(theta2_est);

disp(['  + Uoc tinh huong nguon 1: ' num2str(theta1_deg, '%.1f') ' do']);
disp(['  + Uoc tinh huong nguon 2: ' num2str(theta2_deg, '%.1f') ' do']);
disp('');

%% 4. AP DUNG TAT CA CAC BEAMFORMER CHO 2 HUONG
disp('BUOC 4: AP DUNG TAT CA CAC BEAMFORMER');
disp('-----------------------------------------------------------');

% Luu tat ca ket qua beamforming
beamform_results = {};
idx = 1;

% 4.1. Delay-and-Sum
try
    y1_das = delay_and_sum_beamformer(X_bp, fs_target, theta1_deg);
    y2_das = delay_and_sum_beamformer(X_bp, fs_target, theta2_deg);
    beamform_results{idx} = struct('signal', y1_das, 'method', 'DAS', 'theta', theta1_deg); idx = idx + 1;
    beamform_results{idx} = struct('signal', y2_das, 'method', 'DAS', 'theta', theta2_deg); idx = idx + 1;
    disp('  + Delay-and-Sum: 2 nguon');
catch ME
    disp(['  ! Loi DAS: ' ME.message]);
end

% 4.2. GSC
try
    [y1_gsc, ~] = gsc_beamformer(X_bp, fs_target, theta1_deg, 256, 0.003);
    [y2_gsc, ~] = gsc_beamformer(X_bp, fs_target, theta2_deg, 256, 0.003);
    beamform_results{idx} = struct('signal', y1_gsc, 'method', 'GSC', 'theta', theta1_deg); idx = idx + 1;
    beamform_results{idx} = struct('signal', y2_gsc, 'method', 'GSC', 'theta', theta2_deg); idx = idx + 1;
    disp('  + GSC Beamformer: 2 nguon');
catch ME
    disp(['  ! Loi GSC: ' ME.message]);
end

% 4.3. MVDR
try
    y1_mvdr = mvdr_beamformer(X_bp, fs_target, theta1_deg);
    y2_mvdr = mvdr_beamformer(X_bp, fs_target, theta2_deg);
    beamform_results{idx} = struct('signal', y1_mvdr, 'method', 'MVDR', 'theta', theta1_deg); idx = idx + 1;
    beamform_results{idx} = struct('signal', y2_mvdr, 'method', 'MVDR', 'theta', theta2_deg); idx = idx + 1;
    disp('  + MVDR Beamformer: 2 nguon');
catch ME
    disp(['  ! Loi MVDR: ' ME.message]);
end

% 4.4. LCMV
try
    % LCMV voi null constraint
    y1_lcmv = lcmv_beamformer(X_bp, fs_target, theta1_deg, theta2_deg);
    y2_lcmv = lcmv_beamformer(X_bp, fs_target, theta2_deg, theta1_deg);
    beamform_results{idx} = struct('signal', y1_lcmv, 'method', 'LCMV', 'theta', theta1_deg); idx = idx + 1;
    beamform_results{idx} = struct('signal', y2_lcmv, 'method', 'LCMV', 'theta', theta2_deg); idx = idx + 1;
    disp('  + LCMV Beamformer: 2 nguon');
catch ME
    disp(['  ! Loi LCMV: ' ME.message]);
end

% 4.5. Differential
try
    y_diff = differential_microphone_array(X_bp, 1);
    y_diff_inv = -y_diff;  % Dao nguoc
    beamform_results{idx} = struct('signal', y_diff, 'method', 'Differential', 'theta', 0); idx = idx + 1;
    beamform_results{idx} = struct('signal', y_diff_inv, 'method', 'Differential-Inv', 'theta', 0); idx = idx + 1;
    disp('  + Differential Array: 2 nguon');
catch ME
    disp(['  ! Loi Differential: ' ME.message]);
end

disp(['  => Tong cong: ' num2str(length(beamform_results)) ' ket qua beamforming']);
disp('');

%% 5. ICA CO BAN
disp('BUOC 5: ICA CO BAN');
disp('-----------------------------------------------------------');

try
    [S_ica, ~, ~] = fastica_robust(X_bp);
    
    if size(S_ica, 1) >= 2
        S1_ica = S_ica(1, :) / (max(abs(S_ica(1, :))) + eps);
        S2_ica = S_ica(2, :) / (max(abs(S_ica(2, :))) + eps);
        
        beamform_results{idx} = struct('signal', S1_ica, 'method', 'ICA', 'theta', 0); idx = idx + 1;
        beamform_results{idx} = struct('signal', S2_ica, 'method', 'ICA', 'theta', 0); idx = idx + 1;
        
        disp('  + ICA: 2 nguon');
    end
catch ME
    disp(['  ! Loi ICA: ' ME.message]);
end

disp(['  => Tong cong: ' num2str(length(beamform_results)) ' ket qua']);
disp('');

%% 6. CHON CAP KET QUA TOT NHAT (KHAC NHAU NHAT)
disp('BUOC 6: CHON CAP KET QUA TOT NHAT');
disp('-----------------------------------------------------------');

num_results = length(beamform_results);
if num_results < 2
    error('Khong du ket qua!');
end

% Dam bao tat ca ket qua co cung do dai
for i = 1:num_results
    len_i = length(beamform_results{i}.signal);
    if len_i < min_len
        beamform_results{i}.signal = [beamform_results{i}.signal, zeros(1, min_len - len_i)];
    elseif len_i > min_len
        beamform_results{i}.signal = beamform_results{i}.signal(1:min_len);
    end
end

% Tinh correlation matrix
corr_matrix = ones(num_results, num_results);
for i = 1:num_results
    for j = i+1:num_results
        sig1 = beamform_results{i}.signal;
        sig2 = beamform_results{j}.signal;
        
        corr_val = corrcoef(sig1, sig2);
        if ~isnan(corr_val(1, 2))
            corr_matrix(i, j) = abs(corr_val(1, 2));
            corr_matrix(j, i) = corr_matrix(i, j);
        end
    end
end

% Tim cap co correlation thap nhat VA khong giong input
min_corr = inf;
best_i = 1;
best_j = 2;

for i = 1:num_results
    for j = i+1:num_results
        % Kiem tra correlation voi input
        corr_input_i1 = abs(corrcoef(beamform_results{i}.signal, X_bp(1, :)));
        corr_input_i2 = abs(corrcoef(beamform_results{i}.signal, X_bp(2, :)));
        corr_input_j1 = abs(corrcoef(beamform_results{j}.signal, X_bp(1, :)));
        corr_input_j2 = abs(corrcoef(beamform_results{j}.signal, X_bp(2, :)));
        
        max_corr_input_i = max(corr_input_i1(1,2), corr_input_i2(1,2));
        max_corr_input_j = max(corr_input_j1(1,2), corr_input_j2(1,2));
        
        % Chi chon neu CA HAI deu khac input (corr < 0.9)
        if max_corr_input_i < 0.9 && max_corr_input_j < 0.9
            if corr_matrix(i, j) < min_corr
                min_corr = corr_matrix(i, j);
                best_i = i;
                best_j = j;
            end
        end
    end
end

disp(['  + Chon: ' beamform_results{best_i}.method ' va ' beamform_results{best_j}.method]);
disp(['  + Correlation: ' num2str(min_corr, '%.3f')]);

S1_selected = beamform_results{best_i}.signal;
S2_selected = beamform_results{best_j}.signal;

disp('');

%% 7. CAI THIEN BANG MASK TRONG MIEN TAN SO
disp('BUOC 7: CAI THIEN BANG TIME-FREQUENCY MASKING');
disp('-----------------------------------------------------------');

% Chuyen 2 nguon da chon sang mien tan so
[S1_sel_tf, ~, ~] = spectrogram(S1_selected, WINDOW, NOVERLAP, NFFT, fs_target);
[S2_sel_tf, ~, ~] = spectrogram(S2_selected, WINDOW, NOVERLAP, NFFT, fs_target);

% Tao Ideal Binary Mask (IBM) cung
P1 = abs(S1_sel_tf).^2;
P2 = abs(S2_sel_tf).^2;

% Binary mask: gan ve nguon co cong suat lon hon
IBM_1 = double(P1 > P2);
IBM_2 = double(P2 >= P1);

% Lam mem mask (soft binary mask)
IBM_1_soft = imgaussfilt(IBM_1, 1.5);  % Gaussian blur
IBM_2_soft = imgaussfilt(IBM_2, 1.5);

% Normalize
IBM_sum = IBM_1_soft + IBM_2_soft + eps;
IBM_1_soft = IBM_1_soft ./ IBM_sum;
IBM_2_soft = IBM_2_soft ./ IBM_sum;

% Ap dung mask len SPECTROGRAM GOC (quan trong!)
S1_masked = S1_tf .* IBM_1_soft;
S2_masked = S2_tf .* IBM_2_soft;

disp('  + Da tao va ap dung Ideal Binary Mask');
disp('');

%% 8. ITERATIVE REFINEMENT
disp('BUOC 8: ITERATIVE REFINEMENT (LAP 5 LAN)');
disp('-----------------------------------------------------------');

S1_refined = S1_masked;
S2_refined = S2_masked;

for iter = 1:5
    % Tinh lai cong suat
    P1_new = abs(S1_refined).^2;
    P2_new = abs(S2_refined).^2;
    P_total = P1_new + P2_new + eps;
    
    % Wiener mask
    mask_w1 = P1_new ./ P_total;
    mask_w2 = P2_new ./ P_total;
    
    % Ket hop voi spatial mask
    % Spatial mask tu IPD/ILD
    mask_spatial_1 = zeros(size(IPD));
    mask_spatial_2 = zeros(size(IPD));
    
    for f = 1:size(IPD, 1)
        for t = 1:size(IPD, 2)
            dist1 = sqrt((IPD(f,t) - C(1,1))^2 + ((ILD(f,t)/10) - C(1,2))^2);
            dist2 = sqrt((IPD(f,t) - C(2,1))^2 + ((ILD(f,t)/10) - C(2,2))^2);
            
            weight1 = exp(-dist1^2 / 0.5);
            weight2 = exp(-dist2^2 / 0.5);
            
            total = weight1 + weight2 + eps;
            mask_spatial_1(f, t) = weight1 / total;
            mask_spatial_2(f, t) = weight2 / total;
        end
    end
    
    % Ket hop: 70% Wiener + 30% Spatial
    mask_final_1 = 0.7 * mask_w1 + 0.3 * mask_spatial_1;
    mask_final_2 = 0.7 * mask_w2 + 0.3 * mask_spatial_2;
    
    % Binary masking manh
    mask_beta = 2.5;
    mask_final_1 = mask_final_1.^mask_beta;
    mask_final_2 = mask_final_2.^mask_beta;
    
    % Normalize
    mask_sum = mask_final_1 + mask_final_2 + eps;
    mask_final_1 = mask_final_1 ./ mask_sum;
    mask_final_2 = mask_final_2 ./ mask_sum;
    
    % Ap dung len spectrogram GOC
    S1_refined = S1_tf .* mask_final_1;
    S2_refined = S2_tf .* mask_final_2;
    
    % Kiem tra correlation
    y1_temp = my_istft(S1_refined, WINDOW, NOVERLAP, NFFT, min_len);
    y2_temp = my_istft(S2_refined, WINDOW, NOVERLAP, NFFT, min_len);
    
    len_check = min([length(y1_temp), length(y2_temp), min_len]);
    if len_check > 0
        corr_check = corrcoef(y1_temp(1:len_check), y2_temp(1:len_check));
        if ~isnan(corr_check(1, 2))
            disp(['    Iter ' num2str(iter) ': Correlation = ' num2str(abs(corr_check(1, 2)), '%.3f')]);
        end
    end
end

disp('  + Da hoan thanh iterative refinement');
disp('');

%% 9. CHUYEN VE MIEN THOI GIAN
disp('BUOC 9: CHUYEN VE MIEN THOI GIAN');
disp('-----------------------------------------------------------');

nguoi_1_temp = my_istft(S1_refined, WINDOW, NOVERLAP, NFFT, min_len);
nguoi_2_temp = my_istft(S2_refined, WINDOW, NOVERLAP, NFFT, min_len);

% Dam bao do dai
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

%% 10. HAU XU LY CUOI CUNG
disp('BUOC 10: HAU XU LY CUOI CUNG');
disp('-----------------------------------------------------------');

% De-emphasis
nguoi_1 = filter(1, [1, -alpha], nguoi_1);
nguoi_2 = filter(1, [1, -alpha], nguoi_2);

nguoi_1 = nguoi_1 / (max(abs(nguoi_1)) + eps);
nguoi_2 = nguoi_2 / (max(abs(nguoi_2)) + eps);

% Orthogonalization MANH (lap 5 lan)
for orth_iter = 1:5
    proj_coef = (nguoi_2 * nguoi_1') / (nguoi_1 * nguoi_1' + eps);
    nguoi_2 = nguoi_2 - 0.8 * proj_coef * nguoi_1;
    
    proj_coef2 = (nguoi_1 * nguoi_2') / (nguoi_2 * nguoi_2' + eps);
    nguoi_1 = nguoi_1 - 0.8 * proj_coef2 * nguoi_2;
    
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

disp('  + De-emphasis, Orthogonalization, VAD');
disp('');

%% 11. LUU KET QUA
disp('BUOC 11: LUU KET QUA');
disp('-----------------------------------------------------------');

output_file1 = fullfile(result_folder, 'nguoi_1.wav');
output_file2 = fullfile(result_folder, 'nguoi_2.wav');

audiowrite(output_file1, nguoi_1', fs_target);
audiowrite(output_file2, nguoi_2', fs_target);

disp(['  + Luu: nguoi_1.wav']);
disp(['  + Luu: nguoi_2.wav']);
disp('');

%% 12. HIEN THI KET QUA
disp('BUOC 12: HIEN THI KET QUA');
disp('-----------------------------------------------------------');

t = (0:min_len-1) / fs_target;

figure('Name', 'Advanced Hybrid Separation', 'Position', [50, 50, 1700, 950]);

% Row 1: Input
subplot(3, 4, 1);
plot(t, X_bp(1, :));
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Input: Mic 1');
grid on; xlim([0, max(t)]);

subplot(3, 4, 2);
plot(t, X_bp(2, :));
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Input: Mic 2');
grid on; xlim([0, max(t)]);

subplot(3, 4, 3);
imagesc(T, F, 20*log10(abs(S1_tf) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Input Mic1');
ylim([0, 4000]);

subplot(3, 4, 4);
scatter(features(:,1), features(:,2)*10, 5, idx_cluster, 'filled');
xlabel('IPD'); ylabel('ILD');
title(['Clustering: 2 nhom (theta1=' num2str(theta1_deg, '%.0f') ', theta2=' num2str(theta2_deg, '%.0f') ')']);
grid on;

% Row 2: Nguoi 1
subplot(3, 4, 5);
plot(t, nguoi_1);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('OUTPUT: Nguoi 1');
grid on; xlim([0, max(t)]);
set(gca, 'Color', [1 0.95 0.95]);

subplot(3, 4, 6);
[S1_final, F1, T1] = spectrogram(nguoi_1, WINDOW, NOVERLAP, NFFT, fs_target);
imagesc(T1, F1, 20*log10(abs(S1_final) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Nguoi 1');
ylim([0, 4000]);

subplot(3, 4, 7);
Pxx1 = mean(abs(S1_final).^2, 2);
plot(F1, 10*log10(Pxx1 + eps), 'r', 'LineWidth', 2);
xlabel('Tan so (Hz)'); ylabel('Cong suat (dB)');
title('Pho tan so - Nguoi 1');
grid on; xlim([0, 4000]);

subplot(3, 4, 8);
imagesc(IBM_1_soft);
axis xy; colorbar; caxis([0 1]);
xlabel('Frame'); ylabel('Tan so bin');
title('Binary Mask - Nguoi 1');

% Row 3: Nguoi 2
subplot(3, 4, 9);
plot(t, nguoi_2);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('OUTPUT: Nguoi 2');
grid on; xlim([0, max(t)]);
set(gca, 'Color', [0.95 0.95 1]);

subplot(3, 4, 10);
[S2_final, F2, T2] = spectrogram(nguoi_2, WINDOW, NOVERLAP, NFFT, fs_target);
imagesc(T2, F2, 20*log10(abs(S2_final) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Nguoi 2');
ylim([0, 4000]);

subplot(3, 4, 11);
Pxx2 = mean(abs(S2_final).^2, 2);
plot(F2, 10*log10(Pxx2 + eps), 'b', 'LineWidth', 2);
xlabel('Tan so (Hz)'); ylabel('Cong suat (dB)');
title('Pho tan so - Nguoi 2');
grid on; xlim([0, 4000]);

subplot(3, 4, 12);
imagesc(IBM_2_soft);
axis xy; colorbar; caxis([0 1]);
xlabel('Frame'); ylabel('Tan so bin');
title('Binary Mask - Nguoi 2');

savefig(fullfile(result_folder, 'ket_qua_advanced.fig'));
disp('  + Da hien thi ket qua');
disp('');

%% 13. DANH GIA CHI TIET
disp('BUOC 13: DANH GIA CHI TIET');
disp('-----------------------------------------------------------');

energy1 = sum(nguoi_1.^2);
energy2 = sum(nguoi_2.^2);
disp(['  + Nang luong nguoi 1: ' num2str(energy1, '%.2e')]);
disp(['  + Nang luong nguoi 2: ' num2str(energy2, '%.2e')]);

corr_final = corrcoef(nguoi_1, nguoi_2);
if ~isnan(corr_final(1, 2))
    corr_val = abs(corr_final(1, 2));
    disp(['  + Correlation giua 2 nguoi: ' num2str(corr_val, '%.3f')]);
end

% Correlation voi input
corr_in1 = corrcoef(nguoi_1, X_bp(1, :));
corr_in2 = corrcoef(nguoi_2, X_bp(2, :));

if ~isnan(corr_in1(1, 2))
    disp(['  + Correlation nguoi 1 vs input mic 1: ' num2str(abs(corr_in1(1, 2)), '%.3f')]);
end
if ~isnan(corr_in2(1, 2))
    disp(['  + Correlation nguoi 2 vs input mic 2: ' num2str(abs(corr_in2(1, 2)), '%.3f')]);
end

disp('  ');
if ~isnan(corr_final(1, 2))
    if corr_val < 0.15
        disp('  ✓✓✓ KET QUA XUAT SAC! ✓✓✓');
    elseif corr_val < 0.30
        disp('  ✓✓ KET QUA TOT! ✓✓');
    elseif corr_val < 0.50
        disp('  ✓ KET QUA KHAM ✓');
    else
        disp('  ✗ CANH BAO: Con dinh nhieu ✗');
    end
end

disp('');
disp('=============================================================');
disp('=============================================================');
disp('');
disp('THUAT TOAN LAI NANG CAO');
disp('');
disp('Cac buoc:');
disp('  1. Tien xu ly cho giong noi (pre-emphasis, bandpass)');
disp('  2. STFT - chuyen sang mien tan so');
disp('  3. Clustering IPD/ILD de uoc tinh huong 2 nguoi');
disp('  4. Ap dung TAT CA cac beamformer (DAS, GSC, MVDR, LCMV, Diff)');
disp('  5. Chon cap ket qua tot nhat (corr thap, khac input)');
disp('  6. Tao Ideal Binary Mask (IBM) tu cap da chon');
disp('  7. Iterative refinement (lap 5 lan):');
disp('     - Ket hop Wiener mask + Spatial mask');
disp('     - Binary masking manh (mask^2.5)');
disp('     - Ap dung len spectrogram GOC');
disp('  8. ISTFT - chuyen ve mien thoi gian');
disp('  9. Hau xu ly manh: Orthogonalization x5, VAD');
disp('');
disp('Diem manh:');
disp('  - Ket hop TAT CA cac phuong phap beamforming');
disp('  - Chon tu dong cap tot nhat (khac input, khac nhau)');
disp('  - Binary masking manh loai bo phan dinh');
disp('  - Iterative refinement cai thien dan');
disp('  - Orthogonalization lap lai loai bo dinh hoan toan');
disp('=============================================================');
disp('');

