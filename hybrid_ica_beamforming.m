%% hybrid_ica_beamforming.m
% Thuat toan lai: KET HOP BSS-ICA voi Beamforming va Masking
% Phat trien tu BSS-ICA, khong phai chay rieng roi chon
% Dau vao: mic1.wav va mic2.wav
% Dau ra: nguoi_1.wav va nguoi_2.wav

clear all; close all; clc;

disp('=============================================================');
disp('   THUAT TOAN LAI: BSS-ICA + BEAMFORMING + MASKING          ');
disp('=============================================================');
disp('');

%% 1. THAM SO
fs_target = 16000;
result_folder = 'output_hybrid_ica';
if ~exist(result_folder, 'dir'), mkdir(result_folder); end

%% 2. TAI VA TIEN XU LY
disp('BUOC 1: TAI VA TIEN XU LY');
disp('-----------------------------------------------------------');

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

[b_bp, a_bp] = butter(6, [300, 3400]/(fs_target/2), 'bandpass');
X_bp = zeros(size(X_pre));
for i = 1:2
    X_bp(i, :) = filtfilt(b_bp, a_bp, X_pre(i, :));
    X_bp(i, :) = X_bp(i, :) / (max(abs(X_bp(i, :))) + eps);
end

disp('  + Da tien xu ly');
disp('');

%% 3. BUOC 1: ICA CO BAN (Khoi tao)
disp('BUOC 2: ICA CO BAN (Khoi tao)');
disp('-----------------------------------------------------------');

[S_ica, A_ica, W_ica] = fastica_robust(X_bp);

if size(S_ica, 1) < 2
    error('ICA chi tim thay 1 nguon!');
end

% Lay 2 nguon dau
S1_init = S_ica(1, :) / (max(abs(S_ica(1, :))) + eps);
S2_init = S_ica(2, :) / (max(abs(S_ica(2, :))) + eps);

disp(['  + ICA khoi tao: ' num2str(size(S_ica, 1)) ' nguon']);

% Tinh mixing matrix A (gia dinh linear mixing)
% X = A * S => A = X * S' * (S * S')^(-1)
A_est = X_bp * [S1_init; S2_init]' * inv([S1_init; S2_init] * [S1_init; S2_init]' + eye(2)*1e-6);

disp(['  + Uoc tinh mixing matrix A']);
disp(['    A = [' num2str(A_est(1,1), '%.3f') ', ' num2str(A_est(1,2), '%.3f') ']']);
disp(['        [' num2str(A_est(2,1), '%.3f') ', ' num2str(A_est(2,2), '%.3f') ']']);
disp('');

%% 4. BUOC 2: SPATIAL INFORMATION (Tu Beamforming)
disp('BUOC 3: TRÍCH XUAT SPATIAL INFORMATION');
disp('-----------------------------------------------------------');

% Tinh IPD va ILD de uoc tinh huong cua moi nguon
NFFT = 2048;
NOVERLAP = floor(NFFT * 0.75);
WINDOW = hamming(NFFT);

[S1_tf, F, T] = spectrogram(X_bp(1,:), WINDOW, NOVERLAP, NFFT, fs_target);
[S2_tf, ~, ~] = spectrogram(X_bp(2,:), WINDOW, NOVERLAP, NFFT, fs_target);

IPD = angle(S2_tf ./ (S1_tf + eps));
ILD = 20*log10((abs(S2_tf) + eps) ./ (abs(S1_tf) + eps));

% Uoc tinh huong cua 2 nguon tu ket qua ICA
% Nguon co A(1,i) > A(2,i) => gan mic 1 (ben trai)
% Nguon co A(2,i) > A(1,i) => gan mic 2 (ben phai)

if abs(A_est(1,1)) > abs(A_est(2,1))
    % Nguon 1 gan mic 1 hon => ben trai
    theta1_est = -30;  % do (trai)
    theta2_est = 30;   % do (phai)
else
    theta1_est = 30;
    theta2_est = -30;
end

disp(['  + Uoc tinh huong nguon 1: ' num2str(theta1_est) ' do']);
disp(['  + Uoc tinh huong nguon 2: ' num2str(theta2_est) ' do']);
disp('');

%% 5. BUOC 3: SPATIAL MASK (Tu Beamforming + ILD/IPD)
disp('BUOC 4: TAO SPATIAL MASK');
disp('-----------------------------------------------------------');

% Tao spatial mask cho moi nguon dua tren IPD/ILD
% Nguon 1 (ben trai): IPD < 0, ILD > 0
% Nguon 2 (ben phai): IPD > 0, ILD < 0

% Clustering de xac dinh 2 vung spatial
mag_sum = abs(S1_tf) + abs(S2_tf);
threshold = prctile(mag_sum(:), 80);
mask_active = mag_sum > threshold;

features = [IPD(mask_active), ILD(mask_active)/10];
[idx_cluster, C] = kmeans(features, 2, 'Replicates', 5, 'MaxIter', 500);

% Tao soft mask
mask_spatial_1 = zeros(size(IPD));
mask_spatial_2 = zeros(size(IPD));

for f = 1:size(IPD, 1)
    for t = 1:size(IPD, 2)
        % Tinh khoang cach den 2 center
        dist1 = sqrt((IPD(f,t) - C(1,1))^2 + ((ILD(f,t)/10) - C(1,2))^2);
        dist2 = sqrt((IPD(f,t) - C(2,1))^2 + ((ILD(f,t)/10) - C(2,2))^2);
        
        % Soft mask dua tren khoang cach (closer = higher weight)
        weight1 = exp(-dist1^2 / 0.5);
        weight2 = exp(-dist2^2 / 0.5);
        
        % Normalize
        total_weight = weight1 + weight2 + eps;
        mask_spatial_1(f, t) = weight1 / total_weight;
        mask_spatial_2(f, t) = weight2 / total_weight;
    end
end

disp('  + Da tao spatial mask tu IPD/ILD');
disp('');

%% 6. BUOC 4: KET HOP ICA voi SPATIAL MASK
disp('BUOC 5: KET HOP ICA voi SPATIAL MASK');
disp('-----------------------------------------------------------');

% Chuyen ICA sources sang mien tan so
[S1_ica_tf, ~, ~] = spectrogram(S1_init, WINDOW, NOVERLAP, NFFT, fs_target);
[S2_ica_tf, ~, ~] = spectrogram(S2_init, WINDOW, NOVERLAP, NFFT, fs_target);

% Tinh Wiener mask cho moi nguon ICA
P1_ica = abs(S1_ica_tf).^2;
P2_ica = abs(S2_ica_tf).^2;
P_total = P1_ica + P2_ica + eps;

mask_wiener_1 = P1_ica ./ P_total;
mask_wiener_2 = P2_ica ./ P_total;

% KET HOP: Wiener mask (tu ICA) * Spatial mask (tu IPD/ILD)
mask_combined_1 = mask_wiener_1 .* mask_spatial_1;
mask_combined_2 = mask_wiener_2 .* mask_spatial_2;

% Lam mem mask (mask^alpha, alpha < 1 => mem hon)
mask_alpha = 1.5;  % > 1 => cung hon (binary), < 1 => mem hon
mask_combined_1 = mask_combined_1.^mask_alpha;
mask_combined_2 = mask_combined_2.^mask_alpha;

% Normalize lai mask
mask_sum = mask_combined_1 + mask_combined_2 + eps;
mask_combined_1 = mask_combined_1 ./ mask_sum;
mask_combined_2 = mask_combined_2 ./ mask_sum;

disp('  + Ket hop Wiener mask (ICA) + Spatial mask (IPD/ILD)');

% Ap dung mask len spectrogram goc
S1_separated = S1_tf .* mask_combined_1;
S2_separated = S2_tf .* mask_combined_2;

% Chuyen ve mien thoi gian
nguoi_1_temp = my_istft(S1_separated, WINDOW, NOVERLAP, NFFT, min_len);
nguoi_2_temp = my_istft(S2_separated, WINDOW, NOVERLAP, NFFT, min_len);

% Dam bao do dai dung
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

disp('  + Da ap dung combined mask');
disp('');

%% 7. HAU XU LY
disp('BUOC 6: HAU XU LY');
disp('-----------------------------------------------------------');

% De-emphasis
nguoi_1 = filter(1, [1, -alpha], nguoi_1);
nguoi_2 = filter(1, [1, -alpha], nguoi_2);

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
    
    vad_threshold = prctile(frame_energy, 25) * 2;
    vad_mask = frame_energy > vad_threshold;
    
    % Mo rong mask
    for f = 2:num_frames-1
        if vad_mask(f-1) || vad_mask(f+1)
            vad_mask(f) = 1;
        end
    end
    
    % Ap dung mask
    sig_vad = zeros(size(sig));
    for f = 1:num_frames
        start_idx = (f - 1) * hop_size + 1;
        end_idx = min(start_idx + hop_size - 1, length(sig));
        if vad_mask(f)
            sig_vad(start_idx:end_idx) = sig(start_idx:end_idx);
        end
    end
    
    if i == 1, nguoi_1 = sig_vad; else, nguoi_2 = sig_vad; end
end

disp('  + De-emphasis, VAD');

% Chuan hoa
nguoi_1 = nguoi_1 / (max(abs(nguoi_1)) + eps);
nguoi_2 = nguoi_2 / (max(abs(nguoi_2)) + eps);

% Orthogonalization de loai bo phan dinh con lai
corr_check = abs(corrcoef(nguoi_1, nguoi_2));
if ~isnan(corr_check(1, 2)) && corr_check(1, 2) > 0.2
    disp(['  + Correlation truoc orthogonalization: ' num2str(corr_check(1, 2), '%.3f')]);
    
    projection = (nguoi_2 * nguoi_1') / (nguoi_1 * nguoi_1' + eps);
    nguoi_2_ortho = nguoi_2 - projection * nguoi_1;
    
    projection2 = (nguoi_1 * nguoi_2_ortho') / (nguoi_2_ortho * nguoi_2_ortho' + eps);
    nguoi_1_ortho = nguoi_1 - projection2 * nguoi_2_ortho;
    
    nguoi_1 = nguoi_1_ortho / (max(abs(nguoi_1_ortho)) + eps);
    nguoi_2 = nguoi_2_ortho / (max(abs(nguoi_2_ortho)) + eps);
    
    corr_after = abs(corrcoef(nguoi_1, nguoi_2));
    if ~isnan(corr_after(1, 2))
        disp(['  + Correlation sau orthogonalization: ' num2str(corr_after(1, 2), '%.3f')]);
    end
end

disp('');

%% 8. LUU KET QUA
disp('BUOC 7: LUU KET QUA');
disp('-----------------------------------------------------------');

output_file1 = fullfile(result_folder, 'nguoi_1.wav');
output_file2 = fullfile(result_folder, 'nguoi_2.wav');

audiowrite(output_file1, nguoi_1', fs_target);
audiowrite(output_file2, nguoi_2', fs_target);

disp(['  + Luu: nguoi_1.wav']);
disp(['  + Luu: nguoi_2.wav']);
disp('');

%% 9. HIEN THI KET QUA
disp('BUOC 8: HIEN THI KET QUA');
disp('-----------------------------------------------------------');

t = (0:min_len-1) / fs_target;

figure('Name', 'Hybrid ICA + Beamforming + Masking', 'Position', [100, 100, 1600, 1000]);

% Row 1: Input + ICA initial
subplot(4, 3, 1);
plot(t, X_bp(1, :));
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Input: Mic 1');
grid on; xlim([0, max(t)]);

subplot(4, 3, 2);
plot(t, X_bp(2, :));
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Input: Mic 2');
grid on; xlim([0, max(t)]);

subplot(4, 3, 3);
imagesc(T, F, 20*log10(abs(S1_tf) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Input');
ylim([0, 4000]);

% Row 2: ICA khoi tao
subplot(4, 3, 4);
plot(t, S1_init);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('ICA khoi tao: Nguon 1');
grid on; xlim([0, max(t)]);

subplot(4, 3, 5);
plot(t, S2_init);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('ICA khoi tao: Nguon 2');
grid on; xlim([0, max(t)]);

subplot(4, 3, 6);
imagesc(mask_combined_1);
axis xy; colorbar;
xlabel('Frame'); ylabel('Tan so bin');
title('Combined Mask - Nguon 1');

% Row 3: Ket qua cuoi cung - Nguoi 1
subplot(4, 3, 7);
plot(t, nguoi_1);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('KET QUA: Nguoi 1');
grid on; xlim([0, max(t)]);
set(gca, 'Color', [1 0.95 0.95]);

subplot(4, 3, 8);
[S1_final, F1, T1] = spectrogram(nguoi_1, WINDOW, NOVERLAP, NFFT, fs_target);
imagesc(T1, F1, 20*log10(abs(S1_final) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Nguoi 1');
ylim([0, 4000]);

subplot(4, 3, 9);
Pxx1 = mean(abs(S1_final).^2, 2);
plot(F1, 10*log10(Pxx1 + eps), 'r', 'LineWidth', 2);
xlabel('Tan so (Hz)'); ylabel('Cong suat (dB)');
title('Pho tan so - Nguoi 1');
grid on; xlim([0, 4000]);

% Row 4: Ket qua cuoi cung - Nguoi 2
subplot(4, 3, 10);
plot(t, nguoi_2);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('KET QUA: Nguoi 2');
grid on; xlim([0, max(t)]);
set(gca, 'Color', [0.95 0.95 1]);

subplot(4, 3, 11);
[S2_final, F2, T2] = spectrogram(nguoi_2, WINDOW, NOVERLAP, NFFT, fs_target);
imagesc(T2, F2, 20*log10(abs(S2_final) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Nguoi 2');
ylim([0, 4000]);

subplot(4, 3, 12);
Pxx2 = mean(abs(S2_final).^2, 2);
plot(F2, 10*log10(Pxx2 + eps), 'b', 'LineWidth', 2);
xlabel('Tan so (Hz)'); ylabel('Cong suat (dB)');
title('Pho tan so - Nguoi 2');
grid on; xlim([0, 4000]);

savefig(fullfile(result_folder, 'ket_qua_hybrid.fig'));
disp('  + Da hien thi ket qua');
disp('');

%% 10. DANH GIA
disp('BUOC 9: DANH GIA KET QUA');
disp('-----------------------------------------------------------');

energy1 = sum(nguoi_1.^2);
energy2 = sum(nguoi_2.^2);
corr_final = corrcoef(nguoi_1, nguoi_2);

disp(['  + Nang luong nguoi 1: ' num2str(energy1, '%.2e')]);
disp(['  + Nang luong nguoi 2: ' num2str(energy2, '%.2e')]);

if ~isnan(corr_final(1, 2))
    corr_val = abs(corr_final(1, 2));
    disp(['  + Correlation giua 2 nguoi: ' num2str(corr_val, '%.3f')]);
    
    disp('  ');
    if corr_val < 0.15
        disp('  *** KET QUA XUAT SAC! ***');
    elseif corr_val < 0.3
        disp('  *** KET QUA TOT! ***');
    elseif corr_val < 0.5
        disp('  *** KET QUA KHAM ***');
    else
        disp('  *** CANH BAO: Con dinh nhieu ***');
    end
end

disp('');
disp('=============================================================');
disp('THUẬT TOÁN LAI: ICA + SPATIAL MASKING');
disp('');
disp('Cac buoc da thuc hien:');
disp('  1. ICA de tach nguon khoi tao');
disp('  2. Uoc tinh mixing matrix A');
disp('  3. Tinh spatial mask tu IPD/ILD (Clustering)');
disp('  4. Tinh Wiener mask tu ICA sources');
disp('  5. Ket hop: Combined mask = Wiener mask * Spatial mask');
disp('  6. Ap dung combined mask len spectrogram goc');
disp('  7. Hau xu ly: De-emphasis, VAD, Orthogonalization');
disp('');
disp('Uu diem:');
disp('  - Ket hop diem manh cua ICA (tach thong ke) va Beamforming (spatial)');
disp('  - Mask mem (soft masking) giu nguyen chi tiet tin hieu');
disp('  - Tu dong loai bo phan dinh');
disp('=============================================================');
disp('');

