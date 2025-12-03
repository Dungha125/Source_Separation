%% speech_separation_combined.m
% Script ket hop cac thuat toan de tach duoc tieng cua 2 nguoi rieng biet
% Muc tieu: Output la 2 file rieng biet cho 2 nguoi, su dung ket hop cac thuat toan
% Cac thuat toan: ICA, GSC, Clustering, Delay-and-Sum, Differential, MVDR, LCMV
clear all; close all; clc;

disp('=== TACH TIENG 2 NGUOI BANG KET HOP CAC THUAT TOAN ===');
disp('Muc tieu: Output 2 file rieng biet cho 2 nguoi');
disp('');

%% 1. THAM SO
fs_target = 16000;
result_folder = 'result_speech_combined';
if ~exist(result_folder, 'dir'), mkdir(result_folder); end

speech_low_freq = 300;
speech_high_freq = 3400;
alpha_pre = 0.97;
t = [];  % Sẽ được khởi tạo sau

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
    
    if fs_in ~= fs_target
        x = resample(x, fs_target, fs_in);
    end
    
    if size(x, 2) > 1
        x = mean(x, 2);
    end
    
    X_data{i} = x;
end

lengths = cellfun(@length, X_data);
min_len = min(lengths);
X = zeros(length(mic_files), min_len);
for i = 1:length(mic_files)
    X(i, :) = X_data{i}(1:min_len)';
end
t = (0:min_len-1) / fs_target;
disp(['  -> Do dai chung: ' num2str(min_len/fs_target) ' giay']);
disp('');

%% 3. TIEN XU LY
disp('--- Tien xu ly tin hieu ---');
% Pre-emphasis
X_pre = zeros(size(X));
for i = 1:size(X, 1)
    X_pre(i, :) = filter([1, -alpha_pre], 1, X(i, :));
end

% Bandpass filter
[b_bp, a_bp] = butter(4, [speech_low_freq, speech_high_freq]/(fs_target/2), 'bandpass');
X_filt = zeros(size(X_pre));
for i = 1:size(X_pre, 1)
    X_filt(i, :) = filter(b_bp, a_bp, X_pre(i, :));
end

% High-pass
[b_hp, a_hp] = butter(4, 80/(fs_target/2), 'high');
for i = 1:size(X_filt, 1)
    X_filt(i, :) = filter(b_hp, a_hp, X_filt(i, :));
end

% Chuan hoa
for i = 1:size(X_filt, 1)
    X_filt(i, :) = X_filt(i, :) / (max(abs(X_filt(i, :))) + eps);
end
disp('  -> Da tien xu ly');
disp('');

%% 4. TACH NGUON BANG NHIEU PHUONG PHAP
disp('=== TACH NGUON BANG NHIEU PHUONG PHAP ===');
all_results = struct();
result_idx = 1;

%% 4.1. BSS-ICA
disp('--- Phuong phap 1: BSS-ICA ---');
try
    [S_ica, ~, ~] = fastica_robust(X_filt);
    num_ica = size(S_ica, 1);
    
    S_ica_processed = zeros(num_ica, min_len);
    for i = 1:num_ica
        sig = S_ica(i, :) / (max(abs(S_ica(i, :))) + eps);
        sig = filter(1, [1, -alpha_pre], sig); % De-emphasis
        sig = sig / (max(abs(sig)) + eps);
        S_ica_processed(i, 1:min(length(sig), min_len)) = sig(1:min(length(sig), min_len));
        
        all_results(result_idx).method = 'ICA';
        all_results(result_idx).signal = sig(1:min(length(sig), min_len));
        all_results(result_idx).index = i;
        result_idx = result_idx + 1;
    end
    
    save_method_results('BSS_ICA', S_ica_processed, result_folder, fs_target, t, X_filt);
    disp(['  -> Tim thay ' num2str(num_ica) ' nguon bang ICA']);
catch ME
    disp(['  Loi ICA: ' ME.message]);
end
disp('');

%% 4.2. GSC BEAMFORMER
disp('--- Phuong phap 2: GSC Beamformer ---');
try
    look_directions = [-60, -45, -30, -15, 0, 15, 30, 45, 60];
    filter_length = 256;
    mu = 0.005;
    
    gsc_results = [];
    gsc_energies = [];
    
    for dir_idx = 1:length(look_directions)
        theta = look_directions(dir_idx);
        [y_enhanced, ~] = gsc_beamformer(X_filt, fs_target, theta, filter_length, mu);
        energy = sum(y_enhanced.^2);
        gsc_results = [gsc_results; y_enhanced];
        gsc_energies = [gsc_energies, energy];
    end
    
    % Chon 2 huong tot nhat
    [~, sorted_idx] = sort(gsc_energies, 'descend');
    best_idx = sorted_idx(1:min(2, length(sorted_idx)));
    
    gsc_selected = zeros(length(best_idx), min_len);
    for i = 1:length(best_idx)
        sig = gsc_results(best_idx(i), :);
        sig = sig / (max(abs(sig)) + eps);
        sig = filter(1, [1, -alpha_pre], sig);
        sig = sig / (max(abs(sig)) + eps);
        gsc_selected(i, 1:min(length(sig), min_len)) = sig(1:min(length(sig), min_len));
        
        all_results(result_idx).method = 'GSC';
        all_results(result_idx).signal = sig(1:min(length(sig), min_len));
        all_results(result_idx).index = look_directions(best_idx(i));
        result_idx = result_idx + 1;
    end
    
    save_method_results('GSC_Beamformer', gsc_selected, result_folder, fs_target, t, X_filt);
    disp(['  -> Tim thay ' num2str(length(best_idx)) ' nguon bang GSC']);
catch ME
    disp(['  Loi GSC: ' ME.message]);
end
disp('');

%% 4.3. CLUSTERING (IPD/ILD)
disp('--- Phuong phap 3: Clustering (IPD/ILD) ---');
try
    NFFT = 2048;
    NOVERLAP = floor(NFFT * 0.75);
    WINDOW = hamming(NFFT);
    
    [S1, ~, ~] = spectrogram(X_filt(1,:), WINDOW, NOVERLAP, NFFT, fs_target);
    [S2, ~, ~] = spectrogram(X_filt(2,:), WINDOW, NOVERLAP, NFFT, fs_target);
    
    IPD = angle(S2 ./ (S1 + eps));
    ILD = 20*log10((abs(S2) + eps) ./ (abs(S1) + eps));
    mag_sum = abs(S1) + abs(S2);
    threshold = prctile(mag_sum(:), 85);
    mask_active = mag_sum > threshold;
    features = [IPD(mask_active), ILD(mask_active)/10];
    
    [idx_cluster, C] = kmeans(features, 2, 'Replicates', 3);
    
    cluster_results = zeros(2, min_len);
    for i = 1:2
        dist_map = sqrt((IPD - C(i,1)).^2 + ((ILD/10) - C(i,2)).^2);
        mask_i = exp(-dist_map.^2 / 0.5);
        
        x_sep = my_istft(S1 .* mask_i, WINDOW, NOVERLAP, NFFT, min_len);
        sig = x_sep / (max(abs(x_sep)) + eps);
        cluster_results(i, 1:min(length(sig), min_len)) = sig(1:min(length(sig), min_len));
        
        all_results(result_idx).method = 'Clustering';
        all_results(result_idx).signal = sig(1:min(length(sig), min_len));
        all_results(result_idx).index = i;
        result_idx = result_idx + 1;
    end
    
    save_method_results('Clustering', cluster_results, result_folder, fs_target, t, X_filt);
    disp('  -> Tim thay 2 nguon bang Clustering');
catch ME
    disp(['  Loi Clustering: ' ME.message]);
end
disp('');

%% 4.4. DELAY-AND-SUM BEAMFORMER
disp('--- Phuong phap 4: Delay-and-Sum Beamformer ---');
try
    look_directions = [-60, -45, -30, -15, 0, 15, 30, 45, 60];
    
    das_results = [];
    das_energies = [];
    
    for dir_idx = 1:length(look_directions)
        theta = look_directions(dir_idx);
        y = delay_and_sum_beamformer(X_filt, fs_target, theta);
        energy = sum(y.^2);
        das_results = [das_results; y];
        das_energies = [das_energies, energy];
    end
    
    % Chon 2 huong tot nhat
    [~, sorted_idx] = sort(das_energies, 'descend');
    best_idx = sorted_idx(1:min(2, length(sorted_idx)));
    
    das_selected = zeros(length(best_idx), min_len);
    for i = 1:length(best_idx)
        sig = das_results(best_idx(i), :);
        sig = filter(1, [1, -alpha_pre], sig);
        sig = sig / (max(abs(sig)) + eps);
        das_selected(i, 1:min(length(sig), min_len)) = sig(1:min(length(sig), min_len));
        
        all_results(result_idx).method = 'DelayAndSum';
        all_results(result_idx).signal = sig(1:min(length(sig), min_len));
        all_results(result_idx).index = look_directions(best_idx(i));
        result_idx = result_idx + 1;
    end
    
    save_method_results('DelayAndSum', das_selected, result_folder, fs_target, t, X_filt);
    disp(['  -> Tim thay ' num2str(length(best_idx)) ' nguon bang Delay-and-Sum']);
catch ME
    disp(['  Loi Delay-and-Sum: ' ME.message]);
end
disp('');

%% 4.5. DIFFERENTIAL MICROPHONE ARRAY
disp('--- Phuong phap 5: Differential Microphone Array ---');
try
    diff_results = zeros(2, min_len);
    
    % First-order
    y1 = differential_microphone_array(X_filt, 1);
    y1 = filter(1, [1, -alpha_pre], y1);
    y1 = y1 / (max(abs(y1)) + eps);
    diff_results(1, 1:min(length(y1), min_len)) = y1(1:min(length(y1), min_len));
    
    all_results(result_idx).method = 'Differential';
    all_results(result_idx).signal = y1(1:min(length(y1), min_len));
    all_results(result_idx).index = 1;
    result_idx = result_idx + 1;
    
    % Second-order (nếu có thể)
    if size(X_filt, 1) >= 2
        y2 = -y1;  % Đảo ngược để có nguồn thứ 2
        y2 = y2 / (max(abs(y2)) + eps);
        diff_results(2, 1:min(length(y2), min_len)) = y2(1:min(length(y2), min_len));
        
        all_results(result_idx).method = 'Differential';
        all_results(result_idx).signal = y2(1:min(length(y2), min_len));
        all_results(result_idx).index = 2;
        result_idx = result_idx + 1;
    end
    
    save_method_results('Differential', diff_results, result_folder, fs_target, t, X_filt);
    disp('  -> Tim thay 2 nguon bang Differential');
catch ME
    disp(['  Loi Differential: ' ME.message]);
end
disp('');

%% 4.6. MVDR BEAMFORMER
disp('--- Phuong phap 6: MVDR Beamformer ---');
try
    look_directions = [-60, -45, -30, -15, 0, 15, 30, 45, 60];
    
    mvdr_results = [];
    mvdr_energies = [];
    
    for dir_idx = 1:length(look_directions)
        theta = look_directions(dir_idx);
        try
            y = mvdr_beamformer(X_filt, fs_target, theta);
            energy = sum(y.^2);
            mvdr_results = [mvdr_results; y];
            mvdr_energies = [mvdr_energies, energy];
        catch
            % Bỏ qua nếu lỗi
        end
    end
    
    if ~isempty(mvdr_results)
        % Chon 2 huong tot nhat
        [~, sorted_idx] = sort(mvdr_energies, 'descend');
        best_idx = sorted_idx(1:min(2, length(sorted_idx)));
        
        mvdr_selected = zeros(length(best_idx), min_len);
        for i = 1:length(best_idx)
            sig = mvdr_results(best_idx(i), :);
            sig = filter(1, [1, -alpha_pre], sig);
            sig = sig / (max(abs(sig)) + eps);
            mvdr_selected(i, 1:min(length(sig), min_len)) = sig(1:min(length(sig), min_len));
            
            all_results(result_idx).method = 'MVDR';
            all_results(result_idx).signal = sig(1:min(length(sig), min_len));
            all_results(result_idx).index = look_directions(best_idx(i));
            result_idx = result_idx + 1;
        end
        
        save_method_results('MVDR', mvdr_selected, result_folder, fs_target, t, X_filt);
        disp(['  -> Tim thay ' num2str(length(best_idx)) ' nguon bang MVDR']);
    end
catch ME
    disp(['  Loi MVDR: ' ME.message]);
end
disp('');

%% 4.7. LCMV BEAMFORMER
disp('--- Phuong phap 7: LCMV Beamformer ---');
try
    look_directions = [-60, -45, -30, -15, 0, 15, 30, 45, 60];
    
    lcmv_results = [];
    lcmv_energies = [];
    
    for dir_idx = 1:length(look_directions)
        theta = look_directions(dir_idx);
        % Tạo null constraints từ các hướng khác
        null_dirs = look_directions;
        null_dirs(null_dirs == theta) = [];
        null_dirs = null_dirs(1:min(2, length(null_dirs)));  % Chọn 2 hướng null
        
        try
            y = lcmv_beamformer(X_filt, fs_target, theta, null_dirs);
            energy = sum(y.^2);
            lcmv_results = [lcmv_results; y];
            lcmv_energies = [lcmv_energies, energy];
        catch
            % Bỏ qua nếu lỗi
        end
    end
    
    if ~isempty(lcmv_results)
        % Chon 2 huong tot nhat
        [~, sorted_idx] = sort(lcmv_energies, 'descend');
        best_idx = sorted_idx(1:min(2, length(sorted_idx)));
        
        lcmv_selected = zeros(length(best_idx), min_len);
        for i = 1:length(best_idx)
            sig = lcmv_results(best_idx(i), :);
            sig = filter(1, [1, -alpha_pre], sig);
            sig = sig / (max(abs(sig)) + eps);
            lcmv_selected(i, 1:min(length(sig), min_len)) = sig(1:min(length(sig), min_len));
            
            all_results(result_idx).method = 'LCMV';
            all_results(result_idx).signal = sig(1:min(length(sig), min_len));
            all_results(result_idx).index = look_directions(best_idx(i));
            result_idx = result_idx + 1;
        end
        
        save_method_results('LCMV', lcmv_selected, result_folder, fs_target, t, X_filt);
        disp(['  -> Tim thay ' num2str(length(best_idx)) ' nguon bang LCMV']);
    end
catch ME
    disp(['  Loi LCMV: ' ME.message]);
end
disp('');

%% 6. KET HOP VA CHON KET QUA TOT NHAT
disp('=== KET HOP VA CHON KET QUA TOT NHAT ===');
num_results = result_idx - 1;
disp(['Tong cong co ' num2str(num_results) ' ket qua tu cac phuong phap']);

if num_results < 2
    error('Khong du ket qua de tach 2 nguoi!');
end

% Chuan hoa do dai tat ca cac ket qua ve cung do dai
for i = 1:num_results
    len = length(all_results(i).signal);
    if len < min_len
        all_results(i).signal = [all_results(i).signal, zeros(1, min_len - len)];
    elseif len > min_len
        all_results(i).signal = all_results(i).signal(1:min_len);
    end
end

% Tinh cac tieu chi danh gia cho moi ket qua
disp('  -> Dang danh gia chat luong cac ket qua...');
for i = 1:num_results
    sig = all_results(i).signal;
    
    % 1. Nang luong (energy)
    all_results(i).energy = sum(sig.^2);
    
    % 2. Entropy (do phan tan - cao hon = nhieu thong tin hon)
    sig_power = sig.^2;
    sig_power = sig_power / (sum(sig_power) + eps);
    all_results(i).entropy = -sum(sig_power .* log2(sig_power + eps));
    
    % 3. Zero crossing rate (ZCR) - dac trung cho giong noi
    sig_diff = diff(sign(sig));
    all_results(i).zcr = sum(sig_diff ~= 0) / length(sig);
    
    % 4. Spectral centroid (tan so trung tam)
    NFFT = 2048;
    NOVERLAP = floor(NFFT * 0.75);
    WINDOW = hamming(NFFT);
    [S, F, ~] = spectrogram(sig, WINDOW, NOVERLAP, NFFT, fs_target);
    S_mag = abs(S).^2;
    S_sum = sum(S_mag, 1);
    S_sum(S_sum == 0) = eps;
    centroid = sum(F .* sum(S_mag, 2)) / sum(S_sum);
    all_results(i).spectral_centroid = centroid;
end

% Tinh do tuong quan va do khac biet giua cac ket qua
correlation_matrix = zeros(num_results, num_results);
difference_score = zeros(num_results, num_results);

for i = 1:num_results
    for j = 1:num_results
        if i ~= j
            sig1 = all_results(i).signal;
            sig2 = all_results(j).signal;
            
            % Correlation
            corr_val = corrcoef(sig1, sig2);
            if ~isnan(corr_val(1, 2))
                correlation_matrix(i, j) = abs(corr_val(1, 2));
            else
                correlation_matrix(i, j) = 0;
            end
            
            % Difference score (cao hon = khac nhau hon)
            % Kết hợp nhiều tiêu chí
            energy_diff = abs(all_results(i).energy - all_results(j).energy) / (max(all_results(i).energy, all_results(j).energy) + eps);
            entropy_diff = abs(all_results(i).entropy - all_results(j).entropy);
            zcr_diff = abs(all_results(i).zcr - all_results(j).zcr);
            centroid_diff = abs(all_results(i).spectral_centroid - all_results(j).spectral_centroid) / (max(all_results(i).spectral_centroid, all_results(j).spectral_centroid) + eps);
            
            % Tổng hợp difference score (càng cao càng tốt = càng khác nhau)
            difference_score(i, j) = energy_diff + entropy_diff + zcr_diff + centroid_diff;
        end
    end
end

% Tim cap ket qua tot nhat
% Tiêu chí: correlation thấp VÀ difference score cao VÀ đảm bảo 2 nguồn khác nhau
best_score = -inf;
best_pair = [1, 2];
candidate_pairs = [];

% Thu thập tất cả các cặp ứng viên
for i = 1:num_results-1
    for j = i+1:num_results
        corr_val = correlation_matrix(i, j);
        diff_val = difference_score(i, j);
        combined_score = (1 - corr_val) * diff_val;
        
        % Điều kiện nghiêm ngặt hơn: correlation < 0.5 và difference > 0.5
        if corr_val < 0.5 && diff_val > 0.5
            candidate_pairs = [candidate_pairs; i, j, corr_val, diff_val, combined_score];
        end
    end
end

% Nếu có ứng viên, chọn cặp tốt nhất
if ~isempty(candidate_pairs)
    [~, best_idx] = max(candidate_pairs(:, 5));  % Chọn combined_score cao nhất
    best_pair = [candidate_pairs(best_idx, 1), candidate_pairs(best_idx, 2)];
    disp(['  -> Tim thay ' num2str(size(candidate_pairs, 1)) ' cap ung vien']);
else
    % Nếu không có ứng viên tốt, thử điều kiện lỏng hơn
    disp('  -> Canh bao: Khong tim thay cap ket qua tot (corr < 0.5)');
    disp('  -> Thu dieu kien long hon (corr < 0.6)...');
    
    for i = 1:num_results-1
        for j = i+1:num_results
            corr_val = correlation_matrix(i, j);
            diff_val = difference_score(i, j);
            combined_score = (1 - corr_val) * diff_val;
            
            if corr_val < 0.6 && diff_val > 0.3 && combined_score > best_score
                best_score = combined_score;
                best_pair = [i, j];
            end
        end
    end
    
    % Nếu vẫn không có, chọn cặp có correlation thấp nhất
    if best_score == -inf
        disp('  -> Chon cap co correlation thap nhat...');
        min_corr = inf;
        for i = 1:num_results-1
            for j = i+1:num_results
                if correlation_matrix(i, j) < min_corr
                    min_corr = correlation_matrix(i, j);
                    best_pair = [i, j];
                end
            end
        end
    end
end

disp(['  -> Chon cap ket qua tot nhat: ' all_results(best_pair(1)).method ' va ' all_results(best_pair(2)).method]);
disp(['  -> Do tuong quan: ' num2str(correlation_matrix(best_pair(1), best_pair(2)), '%.3f')]);
disp(['  -> Difference score: ' num2str(difference_score(best_pair(1), best_pair(2)), '%.3f')]);
disp(['  -> Nang luong nguoi 1: ' num2str(all_results(best_pair(1)).energy, '%.2e')]);
disp(['  -> Nang luong nguoi 2: ' num2str(all_results(best_pair(2)).energy, '%.2e')]);

% Lay ket qua da duoc chuan hoa
speaker1 = all_results(best_pair(1)).signal;
speaker2 = all_results(best_pair(2)).signal;

% Validation: Kiểm tra xem kết quả có giống input không
corr_with_input1 = abs(corrcoef(speaker1, X_filt(1, :)));
corr_with_input2 = abs(corrcoef(speaker2, X_filt(2, :)));

if ~isnan(corr_with_input1(1, 2)) && corr_with_input1(1, 2) > 0.8
    disp('  -> Canh bao: Nguoi 1 co ve giong voi input microphone 1!');
    disp('  -> Dang tim ket qua khac...');
    
    % Tìm kết quả khác có correlation thấp với input
    best_alt = best_pair(1);
    min_corr_input = corr_with_input1(1, 2);
    
    for i = 1:num_results
        if i ~= best_pair(2)
            corr_temp = abs(corrcoef(all_results(i).signal, X_filt(1, :)));
            if ~isnan(corr_temp(1, 2)) && corr_temp(1, 2) < min_corr_input
                min_corr_input = corr_temp(1, 2);
                best_alt = i;
            end
        end
    end
    
    if best_alt ~= best_pair(1)
        speaker1 = all_results(best_alt).signal;
        disp(['  -> Da thay doi nguoi 1 sang: ' all_results(best_alt).method]);
    end
end

if ~isnan(corr_with_input2(1, 2)) && corr_with_input2(1, 2) > 0.8
    disp('  -> Canh bao: Nguoi 2 co ve giong voi input microphone 2!');
    disp('  -> Dang tim ket qua khac...');
    
    best_alt = best_pair(2);
    min_corr_input = corr_with_input2(1, 2);
    
    for i = 1:num_results
        if i ~= best_pair(1)
            corr_temp = abs(corrcoef(all_results(i).signal, X_filt(2, :)));
            if ~isnan(corr_temp(1, 2)) && corr_temp(1, 2) < min_corr_input
                min_corr_input = corr_temp(1, 2);
                best_alt = i;
            end
        end
    end
    
    if best_alt ~= best_pair(2)
        speaker2 = all_results(best_alt).signal;
        disp(['  -> Da thay doi nguoi 2 sang: ' all_results(best_alt).method]);
    end
end

% Chuan hoa cuoi cung
speaker1 = speaker1 / (max(abs(speaker1)) + eps);
speaker2 = speaker2 / (max(abs(speaker2)) + eps);

% Kiểm tra lại correlation giữa 2 speaker
final_corr = abs(corrcoef(speaker1, speaker2));
if ~isnan(final_corr(1, 2))
    disp(['  -> Correlation giua 2 nguoi truoc hau xu ly: ' num2str(final_corr(1, 2), '%.3f')]);
end

disp('');

%% 7. HAU XU LY VA LOAI BO PHAN DINH
disp('--- Hau xu ly tin hieu va loai bo phan dinh ---');

% Tính correlation giữa 2 speaker để phát hiện phần dính
corr_speakers = corrcoef(speaker1, speaker2);
if ~isnan(corr_speakers(1, 2))
    corr_val = abs(corr_speakers(1, 2));
    disp(['  -> Correlation giua 2 speaker: ' num2str(corr_val, '%.3f')]);
    
    if corr_val > 0.3
        disp('  -> Phat hien phan dinh, dang loai bo...');
        
        % Phương pháp 1: Loại bỏ phần tương quan cao trong miền thời gian
        window_size = round(0.1 * fs_target);  % 100ms window
        overlap = round(window_size / 2);
        
        % Tính correlation trong từng cửa sổ
        num_windows = floor((min_len - window_size) / (window_size - overlap)) + 1;
        corr_windows = zeros(1, num_windows);
        
        for w = 1:num_windows
            start_idx = (w - 1) * (window_size - overlap) + 1;
            end_idx = min(start_idx + window_size - 1, min_len);
            
            if end_idx > start_idx
                seg1 = speaker1(start_idx:end_idx);
                seg2 = speaker2(start_idx:end_idx);
                corr_temp = corrcoef(seg1, seg2);
                if ~isnan(corr_temp(1, 2))
                    corr_windows(w) = abs(corr_temp(1, 2));
                end
            end
        end
        
        % Tìm các cửa sổ có correlation cao (phần dính)
        high_corr_threshold = 0.5;
        high_corr_mask = corr_windows > high_corr_threshold;
        
        % Loại bỏ phần dính: giảm biên độ ở các cửa sổ có correlation cao
        for w = 1:num_windows
            if high_corr_mask(w)
                start_idx = (w - 1) * (window_size - overlap) + 1;
                end_idx = min(start_idx + window_size - 1, min_len);
                
                if end_idx > start_idx
                    % Giảm biên độ của speaker có năng lượng thấp hơn trong cửa sổ này
                    energy1 = sum(speaker1(start_idx:end_idx).^2);
                    energy2 = sum(speaker2(start_idx:end_idx).^2);
                    
                    % Tính reduction factor mạnh hơn
                    reduction_factor = 1 - corr_windows(w) * 0.8;  % Giảm 80% nếu correlation = 1
                    reduction_factor = max(reduction_factor, 0.05);  % Giữ lại ít nhất 5%
                    
                    if energy1 > energy2
                        % Giảm speaker2 (người có năng lượng thấp hơn)
                        speaker2(start_idx:end_idx) = speaker2(start_idx:end_idx) * reduction_factor;
                    else
                        % Giảm speaker1
                        speaker1(start_idx:end_idx) = speaker1(start_idx:end_idx) * reduction_factor;
                    end
                end
            end
        end
    end
end

% Noise gate và smoothing
noise_threshold = 0.01;
for i = 1:2
    if i == 1
        sig = speaker1;
    else
        sig = speaker2;
    end
    
    % Noise gate
    sig_abs = abs(sig);
    gate_mask = sig_abs > noise_threshold;
    sig_smooth = filter([0.2, 0.6, 0.2], 1, double(gate_mask));
    sig = sig .* sig_smooth;
    
    % Normalize
    sig = sig / (max(abs(sig)) + eps);
    
    if i == 1
        speaker1 = sig;
    else
        speaker2 = sig;
    end
end

% Phương pháp 2: Loại bỏ phần dính trong miền tần số
% Nếu vẫn còn correlation cao, thử loại bỏ trong miền tần số
corr_after = corrcoef(speaker1, speaker2);
if ~isnan(corr_after(1, 2))
    corr_after_val = abs(corr_after(1, 2));
    disp(['  -> Correlation sau khi xu ly (buoc 1): ' num2str(corr_after_val, '%.3f')]);
    
    if corr_after_val > 0.4
        disp('  -> Tiep tuc loai bo phan dinh trong mien tan so...');
        
        % Chuyển sang miền tần số
        NFFT = 2048;
        NOVERLAP = floor(NFFT * 0.75);
        WINDOW = hamming(NFFT);
        
        [S1, F, T] = spectrogram(speaker1, WINDOW, NOVERLAP, NFFT, fs_target);
        [S2, ~, ~] = spectrogram(speaker2, WINDOW, NOVERLAP, NFFT, fs_target);
        
        % Tính correlation trong miền tần số
        S1_mag = abs(S1);
        S2_mag = abs(S2);
        
        % Tìm các tần số có correlation cao
        freq_corr = zeros(size(S1, 1), 1);
        for f_idx = 1:size(S1, 1)
            freq_corr(f_idx) = abs(corrcoef(S1_mag(f_idx, :), S2_mag(f_idx, :)));
            if isnan(freq_corr(f_idx))
                freq_corr(f_idx) = 0;
            end
        end
        
        % Loại bỏ các tần số có correlation cao ở speaker có năng lượng thấp hơn
        high_corr_freqs = find(freq_corr > 0.6);
        if ~isempty(high_corr_freqs)
            % Tính năng lượng tổng của mỗi speaker
            total_energy1 = sum(S1_mag(:).^2);
            total_energy2 = sum(S2_mag(:).^2);
            
            if total_energy1 > total_energy2
                % Giảm S2 ở các tần số có correlation cao
                for f_idx = high_corr_freqs
                    reduction = 1 - freq_corr(f_idx) * 0.7;
                    S2(f_idx, :) = S2(f_idx, :) * reduction;
                end
            else
                % Giảm S1
                for f_idx = high_corr_freqs
                    reduction = 1 - freq_corr(f_idx) * 0.7;
                    S1(f_idx, :) = S1(f_idx, :) * reduction;
                end
            end
            
            % Chuyển lại miền thời gian
            speaker1 = my_istft(S1, WINDOW, NOVERLAP, NFFT, min_len);
            speaker2 = my_istft(S2, WINDOW, NOVERLAP, NFFT, min_len);
            
            speaker1 = speaker1(1:min_len);
            speaker2 = speaker2(1:min_len);
            
            disp(['  -> Da xu ly ' num2str(length(high_corr_freqs)) ' tan so co correlation cao']);
        end
    end
end

% Kiểm tra lại correlation sau khi xử lý
corr_final = corrcoef(speaker1, speaker2);
if ~isnan(corr_final(1, 2))
    disp(['  -> Correlation cuoi cung: ' num2str(abs(corr_final(1, 2)), '%.3f')]);
end

disp('  -> Da hau xu ly');
disp('');

%% 8. LUU KET QUA CUOI CUNG
disp('=== LUU KET QUA CUOI CUNG ===');
output_file1 = fullfile(result_folder, 'nguoi_1.wav');
output_file2 = fullfile(result_folder, 'nguoi_2.wav');

audiowrite(output_file1, speaker1', fs_target);
audiowrite(output_file2, speaker2', fs_target);

disp(['  + Luu: ' output_file1]);
disp(['  + Luu: ' output_file2]);
disp('');

%% 9. HIEN THI KET QUA CUOI CUNG
disp('--- Hien thi ket qua cuoi cung ---');
figure('Name', 'Ket qua cuoi cung: Tach tieng 2 nguoi', 'Position', [100, 100, 1400, 900]);

subplot(3, 2, 1);
plot(t, X_filt(1, :));
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Dau vao: Microphone 1');
grid on; xlim([0, max(t)]);

subplot(3, 2, 2);
plot(t, X_filt(2, :));
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Dau vao: Microphone 2');
grid on; xlim([0, max(t)]);

subplot(3, 2, 3);
plot(t, speaker1);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Ket qua: Nguoi 1');
grid on; xlim([0, max(t)]);

subplot(3, 2, 4);
NFFT = 2048;
NOVERLAP = floor(NFFT * 0.75);
WINDOW = hamming(NFFT);
[S1_spec, F1_spec, T1_spec] = spectrogram(speaker1, WINDOW, NOVERLAP, NFFT, fs_target);
imagesc(T1_spec, F1_spec, 20*log10(abs(S1_spec) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Nguoi 1');
ylim([0, min(4000, fs_target/2)]);

subplot(3, 2, 5);
plot(t, speaker2);
xlabel('Thoi gian (s)'); ylabel('Bien do');
title('Ket qua: Nguoi 2');
grid on; xlim([0, max(t)]);

subplot(3, 2, 6);
[S2_spec, F2_spec, T2_spec] = spectrogram(speaker2, WINDOW, NOVERLAP, NFFT, fs_target);
imagesc(T2_spec, F2_spec, 20*log10(abs(S2_spec) + eps));
axis xy; colorbar;
xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
title('Spectrogram - Nguoi 2');
ylim([0, min(4000, fs_target/2)]);

% Luu figure cuoi cung
savefig(fullfile(result_folder, 'final_results.fig'));
disp('  -> Da hien thi ket qua');
disp('');

%% 10. KET LUAN
disp('=== HOAN THANH ===');
disp(['Ket qua da duoc luu trong thu muc: ' result_folder]);
disp('');
disp('OUTPUT CUOI CUNG:');
disp('  - nguoi_1.wav: Tieng nguoi thu nhat');
disp('  - nguoi_2.wav: Tieng nguoi thu hai');
disp('');
disp('Cac phuong phap da su dung:');
for i = 1:num_results
    disp(['  + ' all_results(i).method ' (index: ' num2str(all_results(i).index) ')']);
end
disp('');
disp('Cap ket qua duoc chon:');
disp(['  - Nguoi 1: ' all_results(best_pair(1)).method]);
disp(['  - Nguoi 2: ' all_results(best_pair(2)).method]);
disp('');
disp('Ket qua chi tiet cua tung phuong phap duoc luu trong cac thu muc con!');
disp('Hay nghe thu 2 file de kiem tra chat luong!');
