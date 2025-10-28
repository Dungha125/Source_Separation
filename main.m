%% main_enhanced.m
% Underdetermined source separation with enhanced frequency domain processing
clear all; close all; clc;
format compact;

dis = 1;
if dis, disp('=== KHOI TAO HE THONG ==='); end

rng('default');
rng(1,'twister');

%% === Tham so he thong ===
M = 1;
u = 0.5;
N = 2;
th = 2.5;
stopthresholdini = 8000;
TC1 = 0.2;
TC2 = 0.08;
numlags = 2;
thepow = 15;
minpow = 25;
evalu = 1;

%% === CAI TIEN: Tham so STFT nang cao ===
winnumber = 3;
NFFT = 4096;  % Tang do phan giai tan so
k = 4;
switch winnumber
    case 1, WINDOW = hanning(NFFT/k);
    case 2, WINDOW = hann(NFFT/k);
    case 3, WINDOW = hamming(NFFT/k);
    case 4, WINDOW = bartlett(NFFT/k);
    case 5, WINDOW = triang(NFFT/k);
    case 6, WINDOW = blackman(NFFT/k);
    case 7, WINDOW = rectwin(NFFT/k);
    otherwise, WINDOW = hamming(NFFT/k);
end
noverlapfactor = 0.75;
NOVERLAP = length(WINDOW)*noverlapfactor;

%% === Tao thu muc ket qua ===
result_folder = 'result_enhanced';
if ~exist(result_folder, 'dir')
    mkdir(result_folder);
end

%% === Load audio sources ===
if evalu
    file_list = {'sp1.wav', 'sp2.wav'};
    
    max_len = 0;
    for i = 1:length(file_list)
        info = audioinfo(file_list{i});
        if info.TotalSamples > max_len
            max_len = info.TotalSamples;
        end
    end
    
    s = zeros(max_len, length(file_list));
    fs = 0;
    for i = 1:length(file_list)
        if i == 1
            [audio_data, fs_temp] = audioread(file_list{i});
            fs = fs_temp;
        else
            [audio_data, ~] = audioread(file_list{i});
        end
        s(1:length(audio_data), i) = audio_data;
    end
    
    labelvec = {'A','B'};
    
    NS = 1:size(s,2);
    Ns = zeros(1,N);
    for i=1:N
        rnd = ceil(rand(1)*(size(s,2)-i+1));
        Ns(i) = NS(rnd);
        NS(rnd) = [];
    end
    N = length(Ns);
    S = s(:,Ns)';
    
    % Hien thi nguon goc
    figure('Name', 'Tin hieu nguon goc', 'NumberTitle', 'off');
    for i = 1:N
        subplot(N, 1, i);
        plot((1:size(S,2))/fs, S(i,:));
        title(sprintf('Nguon %d - %s', i, labelvec{Ns(i)}));
        xlabel('Thoi gian (s)'); ylabel('Bien do'); grid on;
    end
    saveas(gcf, fullfile(result_folder, 'source_waveforms.png'));
    
    % Tao mixing matrix
    Npos = max(N,7);
    all_theta = linspace(0,pi,Npos);
    theta = zeros(1,N);
    for i=1:N
        rnd = ceil(rand(1)*(Npos-i+1));
        theta(i) = all_theta(rnd);
        all_theta(rnd) = [];
    end
    
    A = calcA(theta,u);
    X = A * S;
    
    % Them nhieu
    snr_dB = 25;
    X_noisy = zeros(size(X));
    for ch = 1:size(X,1)
        X_noisy(ch,:) = awgn(X(ch,:), snr_dB, 'measured');
    end
    
    audiowrite(fullfile(result_folder,'stereomix.wav'), X_noisy', fs);
    audiowrite(fullfile(result_folder,'stereomix_clean.wav'), X', fs);
    
    % Hien thi mixture
    figure('Name', 'Tin hieu hon hop', 'NumberTitle', 'off');
    subplot(2,2,1); plot((1:size(X,2))/fs, X(1,:));
    title('Kenh trai (Clean)'); xlabel('Thoi gian (s)'); ylabel('Bien do'); grid on;
    subplot(2,2,2); plot((1:size(X,2))/fs, X(2,:));
    title('Kenh phai (Clean)'); xlabel('Thoi gian (s)'); ylabel('Bien do'); grid on;
    subplot(2,2,3); plot((1:size(X_noisy,2))/fs, X_noisy(1,:));
    title(sprintf('Kenh trai (Noisy - %d dB)', snr_dB));
    xlabel('Thoi gian (s)'); ylabel('Bien do'); grid on;
    subplot(2,2,4); plot((1:size(X_noisy,2))/fs, X_noisy(2,:));
    title(sprintf('Kenh phai (Noisy - %d dB)', snr_dB));
    xlabel('Thoi gian (s)'); ylabel('Bien do'); grid on;
    saveas(gcf, fullfile(result_folder, 'mixture_waveforms.png'));
    
    % === CAI TIEN: Phan tich STFT chi tiet ===
    if dis, disp('Phan tich STFT cua tin hieu hon hop...'); end
    [S_left, F, T] = spectrogram(X(1,:), WINDOW, NOVERLAP, NFFT, fs);
    [S_right, ~, ~] = spectrogram(X(2,:), WINDOW, NOVERLAP, NFFT, fs);
    
    % Hien thi spectrogram
    figure('Name', 'Spectrogram - Kenh trai', 'NumberTitle', 'off');
    imagesc(T, F, 20*log10(abs(S_left)+eps));
    axis xy; colorbar; colormap jet;
    title('Spectrogram - Kenh trai');
    xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
    caxis([max(max(20*log10(abs(S_left)+eps)))-60, max(max(20*log10(abs(S_left)+eps)))]);
    saveas(gcf, fullfile(result_folder, 'spectrogram_left.png'));
    
    figure('Name', 'Spectrogram - Kenh phai', 'NumberTitle', 'off');
    imagesc(T, F, 20*log10(abs(S_right)+eps));
    axis xy; colorbar; colormap jet;
    title('Spectrogram - Kenh phai');
    xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
    caxis([max(max(20*log10(abs(S_right)+eps)))-60, max(max(20*log10(abs(S_right)+eps)))]);
    saveas(gcf, fullfile(result_folder, 'spectrogram_right.png'));
    
    % === CAI TIEN: Phan tich goc pha va IPD ===
    if dis, disp('Tinh toan Inter-channel Phase Difference (IPD)...'); end
    IPD = angle(S_left .* conj(S_right));
    
    figure('Name', 'Inter-channel Phase Difference', 'NumberTitle', 'off');
    imagesc(T, F, IPD);
    axis xy; colorbar; colormap hsv;
    title('Inter-channel Phase Difference (IPD)');
    xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
    caxis([-pi, pi]);
    saveas(gcf, fullfile(result_folder, 'ipd_analysis.png'));
    
    % === CAI TIEN: Phan tich ILD (Inter-channel Level Difference) ===
    if dis, disp('Tinh toan Inter-channel Level Difference (ILD)...'); end
    ILD = 20*log10(abs(S_left)+eps) - 20*log10(abs(S_right)+eps);
    
    figure('Name', 'Inter-channel Level Difference', 'NumberTitle', 'off');
    imagesc(T, F, ILD);
    axis xy; colorbar; colormap jet;
    title('Inter-channel Level Difference (ILD)');
    xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
    caxis([-20, 20]);
    saveas(gcf, fullfile(result_folder, 'ild_analysis.png'));
    
    % === CAI TIEN: Clustering trong khong gian IPD-ILD ===
    if dis, disp('Thuc hien clustering trong khong gian IPD-ILD...'); end
    
    % Chon cac diem co nang luong cao
    magnitude = abs(S_left) + abs(S_right);
    threshold = prctile(magnitude(:), 75); % Chi lay 25% diem manh nhat
    
    valid_idx = magnitude > threshold;
    ipd_vec = IPD(valid_idx);
    ild_vec = ILD(valid_idx);
    
    % K-means clustering
    if length(ipd_vec) > 100
        features = [ipd_vec(:), ild_vec(:)];
        [idx_cluster, centroids] = kmeans(features, N, 'Replicates', 10);
        
        figure('Name', 'IPD-ILD Clustering', 'NumberTitle', 'off');
        scatter(ipd_vec, ild_vec, 10, idx_cluster, 'filled', 'MarkerFaceAlpha', 0.3);
        hold on;
        plot(centroids(:,1), centroids(:,2), 'kx', 'MarkerSize', 15, 'LineWidth', 3);
        hold off;
        xlabel('IPD (rad)'); ylabel('ILD (dB)');
        title('K-means Clustering trong khong gian IPD-ILD');
        colorbar; grid on;
        saveas(gcf, fullfile(result_folder, 'ipd_ild_clustering.png'));
    end
    
    % Mat na ly tuong
    for i=1:N
        vd = zeros(1,N); vd(i)=1; ivd=(vd-1)*(-1);
        [imaskL{i},imaskR{i},SNRiL(i),SNRiR(i)] = idealmask(A*diag(vd)*S, A*diag(ivd)*S, fs, NFFT, WINDOW, NOVERLAP);
    end
    [cmL,cmR] = colorimask(imaskL,imaskR,fs);
    
    figure('Name', 'Mat na ly tuong', 'NumberTitle', 'off');
    subplot(1,2,1); imshow(cmL); title('Kenh trai');
    subplot(1,2,2); imshow(cmR); title('Kenh phai');
    saveas(gcf, fullfile(result_folder, 'ideal_mask.png'));
else
    [X,fs] = audioread(fullfile(result_folder,'stereomix.wav'));
    X = X';
end

powpow = 10*log10((sum(X(1,:).^2)+sum(X(2,:).^2))/(2*size(X,2)));
thE = powpow - thepow;
minpower = powpow - minpow;

%% === CAI TIEN: Tien xu ly trong mien tan so ===
if dis, disp('=== TIEN XU LY MIEN TAN SO ==='); end

% Loc thong cao de loai DC offset
[b_hp, a_hp] = butter(4, 50/(fs/2), 'high');
X_filtered = zeros(size(X));
for ch = 1:size(X,1)
    X_filtered(ch,:) = filtfilt(b_hp, a_hp, X(ch,:));
end

% Giam nhieu bang Wiener filtering
if dis, disp('Ap dung Wiener filtering...'); end
X_denoised = zeros(size(X_filtered));
OriginalLength = size(X,2); % Kich thuoc muc tieu
for ch = 1:size(X_filtered,1)
    [S_noisy, F_vec, T_vec] = spectrogram(X_filtered(ch,:), WINDOW, NOVERLAP, NFFT, fs);
    
    % Uoc luong pho nhieu tu cac frame dau
    noise_frames = 1:min(10, size(S_noisy,2));
    noise_psd = mean(abs(S_noisy(:,noise_frames)).^2, 2);
    
    % Wiener filter
    signal_psd = abs(S_noisy).^2;
    % Broadcast noise_psd thanh ma tran cung kich thuoc voi signal_psd
    noise_psd_matrix = repmat(noise_psd, 1, size(signal_psd, 2));
    wiener_gain = max(0, 1 - noise_psd_matrix ./ (signal_psd + eps));
    S_clean = S_noisy .* wiener_gain;
    
    % Tong hop lai (Line 252 - Thay the code bi loi)
    temp_denoised = real(overlapadd(S_clean, WINDOW, NOVERLAP));
    
    % Dam bao do dai: Chi lay cac mau co san va dien them 0 neu can (Line 254 - Sua loi)
    ActualLength = length(temp_denoised);
    
    % Chi so index toi da co the lay tu temp_denoised
    MaxIndex = min(OriginalLength, ActualLength);
    
    % Gan phan co san
    X_denoised(ch, 1:MaxIndex) = temp_denoised(1:MaxIndex);
    
    % Neu temp_denoised ngan hon, phan con lai cua X_denoised(ch,:) van la 0 (do da khoi tao)
end

% Cap nhat X
X = X_denoised;

%% === Main separation loop ===
x = {X};
mask = {[]};
fmask = {[]};
enermask = {[]};
delete_me_again = sg(X(1,:),NFFT,fs,WINDOW,NOVERLAP);
lastremmask = zeros(size(delete_me_again));
clear delete_me_again

countmax = 20;
finalcnt = 1;
enercnt = 1;
exitcnt = 1;
cnt = 0;
lastfmasklength = 0;

if dis, disp('=== BAT DAU THUAT TOAN TACH ==='); end

while cnt < countmax
    sx = size(x,2);
    cnt = cnt + 1;
    if dis, fprintf('Lap %d/%d...\n', cnt, countmax); end
    
    stopthreshold = stopthresholdini;
    xbuffer = {};
    maskbuffer = {};
    
    for n = 1:sx
        Xn = x{n};
        Xm = Xn - repmat(mean(Xn, 2), 1, size(Xn, 2));
        C = cov(Xm.');
        [E,D] = eig(C);
        d = diag(D);
        d(d<=0) = eps;
        Wwhite = inv(sqrt(D)) * E';
        Xwhite = Wwhite * Xm;
        
        try
            [y{n}, Aest] = icaML(Xwhite);
        catch ME
            warning('ICA failed: %s', ME.message);
            [y{n}, Aest] = icaML(Xn);
        end
        
        for m=1:2
            den = 10*sqrt(var(y{n}(m,:))) + eps;
            y{n}(m,:) = y{n}(m,:) / den;
        end
        
        if evalu
            [newX{1}(1,:),newX{1}(2,:),newX{2}(1,:),newX{2}(2,:), msk{1}, msk{2}] = ...
                applymasks(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', mask{n}, fs, th, NFFT, WINDOW, NOVERLAP, cmR, cmL);
        else
            [newX{1}(1,:),newX{1}(2,:),newX{2}(1,:),newX{2}(2,:), msk{1}, msk{2}] = ...
                applymasks(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', mask{n}, fs, th, NFFT, WINDOW, NOVERLAP);
        end
        
        for m=1:2
            condi = oneortwo_cond(newX{m}(1,:), newX{m}(2,:), fs);
            est = enerstop(newX{m}(1,:), newX{m}(2,:), thE, minpower);
            
            if est == 2
                if dis, disp('  -> Khong phai tin hieu: nang luong qua thap'); end
            elseif est == 1
                [L,R,enermask{enercnt}] = getfinalmask(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', mask{n}, fs, th, m, NFFT, WINDOW, NOVERLAP, 1);
                audiowrite(fullfile(result_folder, sprintf('enerstereo%d.wav', enercnt)), [L, R], fs);
                enercnt = enercnt + 1;
            elseif condi > stopthreshold
                [L,R,fmask{finalcnt}] = getfinalmask(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', mask{n}, fs, 1, m, NFFT, WINDOW, NOVERLAP, 0);
                audiowrite(fullfile(result_folder, sprintf('finalstereo%d.wav', finalcnt)), [L, R], fs);
                
                % === CAI TIEN: Phan tich chi tiet tin hieu tach ===
                if ~isempty(fmask{finalcnt}) && isnumeric(fmask{finalcnt})
                    [S_sep, F_sep, T_sep] = spectrogram(L, WINDOW, NOVERLAP, NFFT, fs);
                    
                    figure('Name', sprintf('Phan tich tin hieu %d', finalcnt), 'NumberTitle', 'off');
                    
                    % Spectrogram
                    subplot(2,2,1);
                    imagesc(T_sep, F_sep, 20*log10(abs(S_sep)+eps));
                    axis xy; colorbar; colormap jet;
                    title(sprintf('Spectrogram - Tin hieu %d', finalcnt));
                    xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
                    
                    % Pho nang luong
                    subplot(2,2,2);
                    energy_profile = mean(abs(S_sep).^2, 2);
                    plot(F_sep, 10*log10(energy_profile+eps));
                    title('Pho nang luong');
                    xlabel('Tan so (Hz)'); ylabel('Cong suat (dB)');
                    grid on;
                    
                    % Mat na
                    subplot(2,2,3);
                    imagesc(fmask{finalcnt});
                    colorbar; title('Mat na thoi gian-tan so');
                    xlabel('Khung thoi gian'); ylabel('Tan so bin');
                    
                    % Waveform
                    subplot(2,2,4);
                    plot((1:length(L))/fs, L);
                    title('Dang song');
                    xlabel('Thoi gian (s)'); ylabel('Bien do');
                    grid on;
                    
                    saveas(gcf, fullfile(result_folder, sprintf('analysis_signal_%d.png', finalcnt)));
                end
                
                finalcnt = finalcnt + 1;
            else
                xbuffer = {xbuffer{:}, newX{m}};
                maskbuffer = {maskbuffer{:}, msk{m}};
            end
        end
    end
    
    x = xbuffer;
    mask = maskbuffer;
    
    if isempty(xbuffer)
        if dis, disp('Dung thuat toan tach'); end
        
        if length(fmask) ~= lastfmasklength
            fmask = multisigcheck(fmask, X(1,:)', X(2,:)', TC1, fs, NFFT, WINDOW, NOVERLAP, numlags);
            lastfmasklength = length(fmask);
            finalcnt = lastfmasklength + 1;
        end
        
        if ~isempty(enermask{1})
            fmask = nosigcorr(fmask, enermask, X(1,:)', X(2,:)', TC2, fs, NFFT, WINDOW, NOVERLAP, numlags, result_folder);
        end
        
        enercnt = 1;
        enermask = {[]};
        
        if isempty(fmask{1})
            if dis, disp('Khong tach duoc tin hieu nao.'); end
            exitcnt = 4; break;
        else
            if evalu
                [L,R,remainingmask] = getremainingmask(X(1,:)', X(2,:)', fmask, fs, NFFT, WINDOW, NOVERLAP, cmL, cmR);
            else
                [L,R,remainingmask] = getremainingmask(X(1,:)', X(2,:)', fmask, fs, NFFT, WINDOW, NOVERLAP);
            end
        end
        
        if isequal(lastremmask, remainingmask)
            audiowrite(fullfile(result_folder, 'remaining.wav'), [L,R], fs);
            break;
        else
            lastremmask = remainingmask;
        end
        
        audiowrite(fullfile(result_folder, 'remaining.wav'), [L,R], fs);
        
        if exitcnt >= 4, break; end
        exitcnt = exitcnt + 1;
    end
end

%% === Evaluation ===
if evalu
    [valL,valR,e1L,e1R,e2L,e2R,q,lbl,cflag] = comparemasks(fmask, imaskL, imaskR, labelvec(Ns), fs, length(s));
    
    for i=1:N
        Xalone(:,:,i) = A(:,i) * S(i,:);
    end
    
    [PLEL,PREL,PLNR,PRNR,SNRL,SNRR,SNRiLi,SNRiRi,SNRxL,SNRxR] = ...
        calcELNR(e1L,e1R,e2L,e2R,imaskL,imaskR,q,NFFT,WINDOW,NOVERLAP,Xalone,lbl);
    
    % Luu ket qua
    mEL = 100*(mean(PLEL)+mean(PREL))/2;
    mNR = 100*(mean(PLNR)+mean(PRNR))/2;
    mSNRi = (mean(SNRiLi)+mean(SNRiRi))/2;
    mSNRo = (mean(SNRL)+mean(SNRR))/2;
    mSNR = mSNRo - mSNRi;
    mSNRx = (mean(SNRxL)+mean(SNRxR))/2;
    
    flag = (length(fmask)~=N || cflag==1);
    
    save(fullfile(result_folder,'results.mat'), 'mEL','mNR','mSNRi','mSNRo','mSNR','mSNRx','flag');
    
    % Bao cao ket qua
    fprintf('\n=== KET QUA DANH GIA ===\n');
    fprintf('Energy Loss: %.2f%%\n', mEL);
    fprintf('Noise Reduction: %.2f%%\n', mNR);
    fprintf('SNR dau vao: %.2f dB\n', mSNRi);
    fprintf('SNR dau ra: %.2f dB\n', mSNRo);
    fprintf('Cai thien SNR: %.2f dB\n', mSNR);
    fprintf('So tin hieu: %d (Mong doi: %d)\n', length(fmask), N);
end

if dis, disp('=== HOAN THANH ==='); end