%% main.m
% Chuong trinh chinh: Underdetermined source separation voi hien thi figure
clear all; close all; clc;
format compact;

dis = 1;
if dis, disp('Initialisation...'); end

rng('default');
rng(1,'twister');

M = 1;
u = 0.5;
N = 2;
th = 2.5;            % Tang manh nguong nhan dang
stopthresholdini = 8000; % Tang nguong dung thuat toan
TC1 = 0.2;           % Tang nguong kiem tra da tin hieu
TC2 = 0.08;          % Tang nguong loai bo nhieu
numlags = 2;         % Tang so lag kiem tra tuong quan
thepow = 15;         % Giam nguong nang luong
minpow = 25;         % Giam nguong cong suat toi thieu
evalu = 1;

winnumber = 3;
NFFT = 2048;
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

%% === Tao thu muc luu ket qua ===
result_folder = 'result';
if ~exist(result_folder, 'dir')
    mkdir(result_folder);
end

%% === Load / create sources and stereo mix ===
if evalu
    file_list = {
        'sp1.wav', 
        'sp2.wav',
   
    };

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

    % Hien thi tin hieu nguon goc
    figure('Name', 'Tin hieu nguon goc (Waveform)', 'NumberTitle', 'off');
    for i = 1:N
        subplot(N, 1, i);
        plot((1:size(S,2))/fs, S(i,:));
        title(sprintf('Nguon %d - %s', i, labelvec{Ns(i)}));
        xlabel('Thoi gian (s)');
        ylabel('Bien do');
        grid on;
    end
    % === [THEM MOI] Luu hinh ===
    saveas(gcf, fullfile(result_folder, 'source_waveforms.png'));
    % ==========================

    Npos = max(N,7);
    all_theta = linspace(0,pi,Npos);
    theta = zeros(1,N);
    for i=1:N
        rnd = ceil(rand(1)*(Npos-i+1));
        theta(i) = all_theta(rnd);
        all_theta(rnd) = [];
    end

    A = calcA(theta,u);
    X = A * S; % Hon hop S?CH (Clean)
    snr_dB = 30;
    X_noisy = zeros(size(X));
    for ch = 1:size(X,1)
        % awgn t? ??ng ?i?u ch?nh công su?t ?? ??t SNR yêu c?u
        X_noisy(ch,:) = awgn(X(ch,:), snr_dB, 'measured'); % Hon hop CO NHIEU (Noisy)
    end

    % ---- Ghi mixture có nhi?u ----
    audiowrite('stereomix.wav', X_noisy', fs);

    % N?u mu?n gi? c? b?n s?ch ?? so sánh:
    audiowrite('stereomix_clean.wav', X', fs);
        % audiowrite('stereomix.wav', X', fs);

    audiowrite('stereomix.wav', X_noisy', fs);
    audiowrite('stereomix_clean.wav', X', fs);
    audiowrite(fullfile(result_folder,'stereomix.wav'), X', fs); % Luu ban SACH vao result
    
    % === [CHINH SUA] Hien thi tin hieu hon hop stereo (Clean vs Noisy) ===
    figure('Name', 'Tin hieu hon hop stereo (Clean vs Noisy)', 'NumberTitle', 'off');
    
    % Tin hieu S?CH (Clean)
    subplot(2,2,1);
    plot((1:size(X,2))/fs, X(1,:));
    title('Kenh trai (Clean)');
    xlabel('Thoi gian (s)'); ylabel('Bien do'); grid on;
    
    subplot(2,2,2);
    plot((1:size(X,2))/fs, X(2,:));
    title('Kenh phai (Clean)');
    xlabel('Thoi gian (s)'); ylabel('Bien do'); grid on;

    % Tin hieu CO NHIEU (Noisy)
    subplot(2,2,3);
    plot((1:size(X_noisy,2))/fs, X_noisy(1,:));
    title(sprintf('Kenh trai (Noisy - %d dB SNR)', snr_dB));
    xlabel('Thoi gian (s)'); ylabel('Bien do'); grid on;

    subplot(2,2,4);
    plot((1:size(X_noisy,2))/fs, X_noisy(2,:));
    title(sprintf('Kenh phai (Noisy - %d dB SNR)', snr_dB));
    xlabel('Thoi gian (s)'); ylabel('Bien do'); grid on;

    saveas(gcf, fullfile(result_folder, 'mixture_waveforms_noisy_vs_clean.png'));
    % ====================================================================
else
    if ~exist(fullfile(result_folder,'stereomix.wav'),'file')
        error('Khong tim thay stereo mix: result/stereomix.wav');
    end
    [X,fs] = audioread(fullfile(result_folder,'stereomix.wav'));
    X = X';
end

powpow = 10*log10((sum(X(1,:).^2)+sum(X(2,:).^2))/(2*size(X,2)));
thE = powpow - thepow;
minpower = powpow - minpow;

if evalu
    for i=1:N
        vd = zeros(1,N); vd(i)=1; ivd=(vd-1)*(-1);
        [imaskL{i},imaskR{i},SNRiL(i),SNRiR(i)] = idealmask(A*diag(vd)*S, A*diag(ivd)*S, fs, NFFT, WINDOW, NOVERLAP);
    end
    [cmL,cmR] = colorimask(imaskL,imaskR,fs);

    % Hien thi mat na mau ly tuong
    figure('Name', 'Mat na mau ly tuong', 'NumberTitle', 'off');
    subplot(1,2,1);
    imshow(cmL);
    title('Mat na mau - Kenh trai');
    
    subplot(1,2,2);
    imshow(cmR);
    title('Mat na mau - Kenh phai');
    
    % Luu hinh
    saveas(gcf, fullfile(result_folder, 'color_mask_ideal.png'));
end

x = {X};
mask = {[]};
fmask = {[]};
enermask = {[]};
delete_me_again = sg(X(1,:),NFFT,fs,WINDOW,NOVERLAP);
lastremmask = zeros(size(delete_me_again));
clear delete_me_again

countmax = 20;   % Giam so lan lap de tranh tach qua muc
finalcnt = 1;
enercnt = 1;
exitcnt = 1;
cnt = 0;
lastfmasklength = 0;

if dis, disp('Starting algorithm...'); end

%% === Main separation loop ===
while cnt < countmax
    sx = size(x,2);
    cnt = cnt + 1;
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
            warning('ICA failed, using unwhitened data: %s', ME.message);
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
                if dis, disp('Not a speech signal - too low energy'); end
            elseif est == 1
                [L,R,enermask{enercnt}] = getfinalmask(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', mask{n}, fs, th, m, NFFT, WINDOW, NOVERLAP, 1);
                stestr = fullfile(result_folder, sprintf('enerstereo%d.wav', enercnt));
                audiowrite(stestr, [L, R], fs);
                
                % Hien thi mat na nang luong
                if ~isempty(enermask{enercnt}) && isnumeric(enermask{enercnt})
                    figure('Name', sprintf('Mat na nang luong %d', enercnt), 'NumberTitle', 'off');
                    imagesc(enermask{enercnt});
                    colorbar;
                    title(sprintf('Mat na nang luong - Tin hieu %d', enercnt));
                    xlabel('Khung thoi gian');
                    ylabel('Tan so bin');
                    saveas(gcf, fullfile(result_folder, sprintf('enermask_%d.png', enercnt)));
                end
                
                enercnt = enercnt + 1;

            elseif condi > stopthreshold
                [L,R,fmask{finalcnt}] = getfinalmask(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', mask{n}, fs, 1, m, NFFT, WINDOW, NOVERLAP, 0);
                stestr = fullfile(result_folder, sprintf('finalstereo%d.wav', finalcnt));
                audiowrite(stestr, [L, R], fs);
                
                % Hien thi mat na cuoi cung
                if ~isempty(fmask{finalcnt}) && isnumeric(fmask{finalcnt})
                    figure('Name', sprintf('Mat na cuoi cung %d', finalcnt), 'NumberTitle', 'off');
                    imagesc(fmask{finalcnt});
                    colorbar;
                    title(sprintf('Mat na cuoi cung - Tin hieu %d', finalcnt));
                    xlabel('Khung thoi gian');
                    ylabel('Tan so bin');
                    saveas(gcf, fullfile(result_folder, sprintf('finalmask_%d.png', finalcnt)));
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
        if dis, disp('Stopping separation algorithm'); end

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
            if dis, disp('No signals segregated.'); end
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
            
            % Hien thi mat na con lai
            if ~isempty(remainingmask) && isnumeric(remainingmask)
                figure('Name', 'Mat na con lai', 'NumberTitle', 'off');
                imagesc(remainingmask);
                colorbar;
                title('Mat na tin hieu con lai');
                xlabel('Khung thoi gian');
                ylabel('Tan so bin');
                saveas(gcf, fullfile(result_folder, 'remaining_mask.png'));
            end
            
            break;
        else
            lastremmask = remainingmask;
        end

        audiowrite(fullfile(result_folder, 'remaining.wav'), [L,R], fs);

        if exitcnt >= 4, break; end
        exitcnt = exitcnt + 1;
    end
end

if isempty(fmask)
    flag = 1; mEL=0; mNR=0; mSNRi=0; mSNRo=0; mSNRx=0;
else
    if length(fmask) ~= lastfmasklength
        fmask = multisigcheck(fmask, X(1,:)', X(2,:)', TC1, fs, NFFT, WINDOW, NOVERLAP, numlags);
    end
    if ~isempty(enermask{1})
        fmask = nosigcorr(fmask, enermask, X(1,:)', X(2,:)', TC2, fs, NFFT, WINDOW, NOVERLAP, numlags, result_folder);
    end
end

if dis, disp('Separation done.'); end

%% === Evaluation ===
if evalu
    [valL,valR,e1L,e1R,e2L,e2R,q,lbl,cflag] = comparemasks(fmask, imaskL, imaskR, labelvec(Ns), fs, length(s));

    % Hien thi so sanh mat na
    if ~iscell(valL) && isnumeric(valL)
        figure('Name', 'So sanh mat na - Kenh trai', 'NumberTitle', 'off');
        imagesc(valL); 
        colorbar;
        title('So sanh mat na - Kenh trai');
        xlabel('Mat na uoc luong');
        ylabel('Mat na ly tuong');
        saveas(gcf, fullfile(result_folder, 'compare_mask_left.png'));
    end
    
    if ~iscell(valR) && isnumeric(valR)
        figure('Name', 'So sanh mat na - Kenh phai', 'NumberTitle', 'off');
        imagesc(valR); 
        colorbar;
        title('So sanh mat na - Kenh phai');
        xlabel('Mat na uoc luong');
        ylabel('Mat na ly tuong');
        saveas(gcf, fullfile(result_folder, 'compare_mask_right.png'));
    end

    for i=1:N
        Xalone(:,:,i) = A(:,i) * S(i,:);
    end

    [PLEL,PREL,PLNR,PRNR,SNRL,SNRR,SNRiLi,SNRiRi,SNRxL,SNRxR] = ...
        calcELNR(e1L,e1R,e2L,e2R,imaskL,imaskR,q,NFFT,WINDOW,NOVERLAP,Xalone,lbl);

    % Hien thi ket qua SNR
    figure('Name', 'Ket qua SNR', 'NumberTitle', 'off');
    subplot(2,1,1);
    bar([SNRiLi; SNRL]');
    legend('SNR dau vao', 'SNR dau ra');
    title('So sanh SNR - Kenh trai');
    xlabel('Tin hieu');
    ylabel('SNR (dB)');
    grid on;
    
    subplot(2,1,2);
    bar([SNRiRi; SNRR]');
    legend('SNR dau vao', 'SNR dau ra');
    title('So sanh SNR - Kenh phai');
    xlabel('Tin hieu');
    ylabel('SNR (dB)');
    grid on;
    saveas(gcf, fullfile(result_folder, 'snr_comparison.png'));

    datafile = fullfile(result_folder, 'data.mat');
    save(datafile, 'PLEL','PREL','PLNR','PRNR','SNRL','SNRR','SNRiL','SNRiR','lbl','th','stopthresholdini','TC1','TC2','thepow','minpow','NFFT','winnumber','k','NOVERLAP','Ns','theta');

    mEL = 100*(mean(PLEL)+mean(PREL))/2;
    mNR = 100*(mean(PLNR)+mean(PRNR))/2;
    mSNRi = (mean(SNRiLi)+mean(SNRiRi))/2;
    mSNRo = (mean(SNRL)+mean(SNRR))/2;
    mSNRx = (mean(SNRxL)+mean(SNRxR))/2;

    if length(fmask)~=N || cflag==1
        flag = 1;
    else
        flag = 0;
    end

    mSNR = mSNRo - mSNRi;
    save(fullfile(result_folder,'mydata.mat'), 'mEL','mNR','mSNRi','mSNRo','mSNR','mSNRx','flag');
    
    % Hien thi ket qua tong hop
    figure('Name', 'Ket qua tong hop', 'NumberTitle', 'off');
    metrics = [mEL, mNR, mSNRi, mSNRo, mSNR];
    bar(metrics);
    set(gca, 'XTickLabel', {'EL (%)', 'NR (%)', 'SNRi (dB)', 'SNRo (dB)', 'Delta SNR (dB)'});
    title('Cac chi so danh gia');
    ylabel('Gia tri');
    grid on;
    saveas(gcf, fullfile(result_folder, 'metrics_summary.png'));
    
    fprintf('\n=== KET QUA DANH GIA ===\n');
    fprintf('Energy Loss (EL): %.2f%%\n', mEL);
    fprintf('Noise Reduction (NR): %.2f%%\n', mNR);
    fprintf('SNR dau vao: %.2f dB\n', mSNRi);
    fprintf('SNR dau ra: %.2f dB\n', mSNRo);
    fprintf('Cai thien SNR: %.2f dB\n', mSNR);
    fprintf('So tin hieu tach duoc: %d (Mong doi: %d)\n', length(fmask), N);
    
    if length(fmask) > N
        fprintf('\n*** CANH BAO: Tach duoc nhieu tin hieu hon mong doi! ***\n');
        fprintf('Nguyen nhan co the:\n');
        fprintf('  - Mot nguon bi tach thanh nhieu phan\n');
        fprintf('  - Phat hien nhieu hoac thanh phan khong mong muon\n');
        fprintf('  - Can tang cac tham so: th, stopthresholdini, TC1, TC2\n');
    elseif length(fmask) < N
        fprintf('\n*** CANH BAO: Tach duoc it tin hieu hon mong doi! ***\n');
        fprintf('Nguyen nhan co the:\n');
        fprintf('  - Hai nguon bi gop chung thanh mot\n');
        fprintf('  - Can giam cac tham so: th, stopthresholdini\n');
    else
        fprintf('\n==> Ket qua TOT: So tin hieu tach dung bang so nguon!\n');
    end
    
    % Hien thi spectrogram cua cac tin hieu tach duoc
    figure('Name', 'Spectrogram cac tin hieu tach duoc', 'NumberTitle', 'off');
    num_sigs = min(length(fmask), 4);
    for i = 1:num_sigs
        subplot(2, 2, i);
        if ~isempty(fmask{i}) && isnumeric(fmask{i})
            imagesc(fmask{i});
            colorbar;
            title(sprintf('Tin hieu tach %d', i));
            xlabel('Khung thoi gian');
            ylabel('Tan so bin');
        end
    end
    saveas(gcf, fullfile(result_folder, 'separated_signals_spec.png'));
    
    % Kiem tra tuong quan giua cac tin hieu tach duoc
    if length(fmask) >= 2
        figure('Name', 'Phan tich tuong quan', 'NumberTitle', 'off');
        
        % Doc cac file audio da tach
        corr_matrix = zeros(length(fmask));
        for i = 1:length(fmask)
            for j = 1:length(fmask)
                try
                    file1 = fullfile(result_folder, sprintf('finalstereo%d.wav', i));
                    file2 = fullfile(result_folder, sprintf('finalstereo%d.wav', j));
                    
                    if exist(file1, 'file') && exist(file2, 'file')
                        [sig1, ~] = audioread(file1);
                        [sig2, ~] = audioread(file2);
                        
                        % Tinh tuong quan
                        min_len = min(length(sig1), length(sig2));
                        c = corrcoef(sig1(1:min_len,1), sig2(1:min_len,1));
                        corr_matrix(i,j) = abs(c(1,2));
                    end
                catch
                    corr_matrix(i,j) = 0;
                end
            end
        end
        
        % Hien thi ma tran tuong quan
        imagesc(corr_matrix);
        colorbar;
        caxis([0 1]);
        title('Ma tran tuong quan giua cac tin hieu');
        xlabel('Tin hieu');
        ylabel('Tin hieu');
        
        % Them gia tri len hinh
        for i = 1:size(corr_matrix,1)
            for j = 1:size(corr_matrix,2)
                if corr_matrix(i,j) > 0
                    text(j, i, sprintf('%.2f', corr_matrix(i,j)), ...
                        'HorizontalAlignment', 'center', ...
                        'Color', 'white', 'FontWeight', 'bold');
                end
            end
        end
        
        saveas(gcf, fullfile(result_folder, 'correlation_matrix.png'));
        
        % Canh bao neu tuong quan qua cao
        fprintf('\n=== PHAN TICH TUONG QUAN ===\n');
        high_corr_found = false;
        for i = 1:size(corr_matrix,1)
            for j = i+1:size(corr_matrix,2)
                if corr_matrix(i,j) > 0.7
                    fprintf('CANH BAO: Tin hieu %d va %d co tuong quan cao (%.2f)\n', ...
                        i, j, corr_matrix(i,j));
                    fprintf('  => Co the la cung mot nguon bi tach thanh 2 phan!\n');
                    high_corr_found = true;
                elseif corr_matrix(i,j) > 0.3
                    fprintf('Luu y: Tin hieu %d va %d co tuong quan trung binh (%.2f)\n', ...
                        i, j, corr_matrix(i,j));
                end
            end
        end
        
        if ~high_corr_found && length(fmask) >= 2
            fprintf('TOT: Cac tin hieu tach duoc co tuong quan thap (< 0.3)\n');
            fprintf('  => Tach thanh cong %d nguon doc lap!\n', length(fmask));
        end
        
        % Phan tich pho nang luong
        figure('Name', 'Phan tich pho nang luong', 'NumberTitle', 'off');
        for i = 1:min(length(fmask), 2)
            subplot(2, 1, i);
            if ~isempty(fmask{i}) && isnumeric(fmask{i})
                % Tinh nang luong trung binh theo tan so
                energy_profile = mean(fmask{i}, 2);
                plot(energy_profile);
                title(sprintf('Pho nang luong - Tin hieu %d', i));
                xlabel('Tan so bin');
                ylabel('Nang luong trung binh');
                grid on;
                
                % Tim dinh pho
                [peaks, locs] = findpeaks(energy_profile, 'MinPeakHeight', max(energy_profile)*0.3);
                hold on;
                plot(locs, peaks, 'rv', 'MarkerFaceColor', 'r');
                legend('Pho nang luong', 'Dinh pho');
                hold off;
            end
        end
        saveas(gcf, fullfile(result_folder, 'energy_spectrum.png'));
        
        % So sanh dac trung pho
        fprintf('\n=== SO SANH DAC TRUNG PHO ===\n');
        for i = 1:min(length(fmask), 2)
            if ~isempty(fmask{i}) && isnumeric(fmask{i})
                energy_profile = mean(fmask{i}, 2);
                [peaks, locs] = findpeaks(energy_profile, 'MinPeakHeight', max(energy_profile)*0.3);
                
                fprintf('Tin hieu %d:\n', i);
                fprintf('  - So dinh pho: %d\n', length(peaks));
                if ~isempty(locs)
                    num_peaks = min(3, length(locs));
                    fprintf('  - Vi tri dinh chinh (bin): %s\n', mat2str(locs(1:num_peaks)));
                    num_peaks = min(3, length(locs));
                    fprintf('  - Tan so tuong ung (Hz): %s\n', ...
                        mat2str(round(locs(1:num_peaks)*fs/NFFT)));
                end
                fprintf('  - Nang luong trung binh: %.4f\n', mean(energy_profile));
                fprintf('  - Nang luong toi da: %.4f\n', max(energy_profile));
            end
        end
    end
end

% === [THEM MOI] Hien thi waveform cac tin hieu da tach ===
if dis, disp('Generating separated waveform plots...'); end
if finalcnt > 1
    num_separated = finalcnt - 1;
    figure('Name', 'Tin hieu da tach (Waveform)', 'NumberTitle', 'off');
    
    % Tinh toan so hang, so cot cho subplot
    num_rows = ceil(sqrt(num_separated));
    num_cols = ceil(num_separated / num_rows);
    
    max_len = 0;
    signals = {};
    fs_sig = fs; % Gia su cung fs

    % Doc tat ca cac file truoc de tim do dai lon nhat
    for i = 1:num_separated
        fname = fullfile(result_folder, sprintf('finalstereo%d.wav', i));
        if exist(fname, 'file')
            [sig, fs_sig] = audioread(fname);
            signals{i} = sig;
            if length(sig) > max_len
                max_len = length(sig);
            end
        else
            signals{i} = [];
        end
    end
    
    if max_len > 0
        t = (0:max_len-1) / fs_sig;

        % Ve subplot
        for i = 1:num_separated
            if ~isempty(signals{i})
                sig = signals{i};
                % Pad zero (them so 0) vao cuoi neu can
                if length(sig) < max_len
                    sig(end+1:max_len, :) = 0;
                end
                
                subplot(num_rows, num_cols, i);
                plot(t, sig(:,1)); % Chi ve kenh Trai (Left)
                title(sprintf('Tin hieu tach %d', i));
                xlabel('Thoi gian (s)');
                ylabel('Bien do');
                grid on;
                xlim([0 t(end)]); % Dam bao truc x giong nhau
            end
        end
        
        saveas(gcf, fullfile(result_folder, 'separated_waveforms.png'));
    else
        if dis, disp('No valid separated audio files found to plot.'); end
    end
else
    if dis, disp('No final signals were saved to plot (finalcnt <= 1).'); end
end
% ========================================================


if dis, disp('Done'); end