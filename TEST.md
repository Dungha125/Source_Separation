%% main_enhanced.m
% Enhanced Underdetermined Source Separation
% Improvements: Better preprocessing, adaptive thresholds, spectral processing

clear all; close all; clc;
format compact;

dis = 1;
if dis, disp('=== Enhanced USS Initialization ==='); end

% Reproducible random seed
rng('default');
rng(1,'twister');

%% ========== ENHANCED PARAMETERS ==========
M = 2;            % number of microphones
u = 0.5;          % cardioid control
N = 3;            % number of sources
th = 1;

% Adaptive stopping criteria
stopthresholdini = 3000;
adaptive_threshold = true;  % NEW: Enable adaptive thresholding

% Enhanced correlation thresholds
TC1 = 0.12;       % Increased from 0.1 for better duplicate detection
TC2 = 0.035;      % Slightly increased from 0.03
numlags = 2;      % Increased from 1 for better correlation estimation

% Energy thresholds
thepow = 18;      % Reduced from 20 for more sensitive detection
minpow = 28;      % Reduced from 30

evalu = 1;

%% ========== ENHANCED STFT PARAMETERS ==========
winnumber = 3;
NFFT = 4096;      % Increased from 2048 for better frequency resolution
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

%% ========== LOAD/CREATE SOURCES ==========
if evalu
    s = [];
    [s(:,1),fs]  = audioread('ukma.wav');
    [s(:,2),~]   = audioread('frma.wav');
    [s(:,3),~]   = audioread('itfe.wav');
    [s(:,4),~]   = audioread('cnfe.wav');
    [s(:,5),~]   = audioread('rufe.wav');
    [s(:,6),~]   = audioread('gema.wav');
    [s(:,7),~]   = audioread('nlma.wav');
    [s(:,8),~]   = audioread('jpfe.wav');
    [s(:,9),~]   = audioread('brfe.wav');
    [s(:,10),~]  = audioread('esma.wav');
    [s(:,11),~]  = audioread('dkma.wav');
    [s(:,12),~]  = audioread('ukfe.wav');
    labelvec = {'A','B','C','D','E','F','G','H','I','J','K','L'};

    % Random source selection
    NS = 1:size(s,2);
    Ns = zeros(1,N);
    for i=1:N
        rnd = ceil(rand(1)*(size(s,2)-i+1));
        Ns(i) = NS(rnd);
        NS(rnd) = [];
    end
    N = length(Ns);
    S = s(:,Ns)';

    % Random directions
    Npos = max(N,7);
    all_theta = linspace(0,pi,Npos);
    theta = zeros(1,N);
    for i=1:N
        rnd = ceil(rand(1)*(Npos-i+1));
        theta(i) = all_theta(rnd);
        all_theta(rnd) = [];
    end

    % Mixing
    A = calcA(theta,u);
    X = A * S;
    audiowrite('stereomix.wav', X', fs);
else
    if ~exist('stereomix.wav','file')
        error('Không tìm thấy stereo mix: stereomix.wav');
    end
    [X,fs] = audioread('stereomix.wav');
    X = X';
end

%% ========== NEW: ENHANCED PREPROCESSING ==========
% Apply high-pass filter to remove DC and low-frequency noise
[b_hp, a_hp] = butter(4, 80/(fs/2), 'high');
X_filtered = zeros(size(X));
for ch = 1:size(X,1)
    X_filtered(ch,:) = filtfilt(b_hp, a_hp, X(ch,:));
end
X = X_filtered;

% Adaptive energy estimation
powpow = 10*log10((sum(X(1,:).^2)+sum(X(2,:).^2))/(2*size(X,2)));
thE = powpow - thepow;
minpower = powpow - minpow;

%% ========== IDEAL MASKS (if evaluation) ==========
if evalu
    for i=1:N
        vd = zeros(1,N); vd(i)=1; ivd=(vd-1)*(-1);
        [imaskL{i},imaskR{i},SNRiL(i),SNRiR(i)] = ...
            idealmask(A*diag(vd)*S, A*diag(ivd)*S, fs, NFFT, WINDOW, NOVERLAP);
    end
    [cmL,cmR] = colorimask(imaskL,imaskR,fs);
end

%% ========== INITIALIZATION ==========
x = {X};
mask = {[]};
fmask = {[]};
enermask = {[]};

delete_me = sg(X(1,:),NFFT,fs,WINDOW,NOVERLAP);
lastremmask = zeros(size(delete_me));
clear delete_me

countmax = 30;
finalcnt = 1;
enercnt = 1;
exitcnt = 1;
cnt = 0;
lastfmasklength = 0;

% NEW: Track separation quality history
separation_quality_history = [];

if dis, disp('=== Starting Enhanced Separation ==='); end

%% ========== MAIN ENHANCED SEPARATION LOOP ==========
while cnt < countmax
    sx = size(x,2);
    cnt = cnt + 1;
    
    % NEW: Adaptive threshold adjustment
    if adaptive_threshold && cnt > 1
        if ~isempty(separation_quality_history)
            recent_quality = mean(separation_quality_history(max(1,end-2):end));
            if recent_quality < 0.5
                stopthreshold = stopthresholdini * 0.8;  % Lower threshold
            else
                stopthreshold = stopthresholdini * 1.1;  % Raise threshold
            end
        else
            stopthreshold = stopthresholdini;
        end
    else
        stopthreshold = stopthresholdini;
    end
    
    if dis, fprintf('Iteration %d/%d (threshold: %.1f)\n', cnt, countmax, stopthreshold); end

    xbuffer = {};
    maskbuffer = {};

    for n = 1:sx
        %% ========== ENHANCED PCA WHITENING ==========
        Xn = x{n};
        
        % Remove mean
        Xmean = mean(Xn, 2);
        Xm = Xn - repmat(Xmean, 1, size(Xn,2));
        
        % Robust covariance estimation
        C = cov(Xm.');
        
        % Eigenvalue decomposition
        [E, D] = eig(C);
        
        % NEW: More robust eigenvalue handling
        d = diag(D);
        d_original = d;
        
        % Set very small eigenvalues to small positive value
        threshold_eig = max(d) * 1e-6;
        d(d < threshold_eig) = threshold_eig;
        
        % Whitening matrix with regularization
        Wwhite = diag(1./sqrt(d)) * E';
        Xwhite = Wwhite * Xm;
        
        % Additional normalization
        for ch = 1:size(Xwhite,1)
            Xwhite(ch,:) = Xwhite(ch,:) / (std(Xwhite(ch,:)) + eps);
        end

        %% ========== ICA WITH ERROR HANDLING ==========
        try
            [y{n}, Aest] = icaML(Xwhite);
            ica_success = true;
        catch ME
            warning('ICA failed: %s. Using fallback method.', ME.message);
            try
                % Fallback to simpler whitening
                [y{n}, Aest] = icaML(Xm);
                ica_success = true;
            catch
                warning('Both ICA attempts failed. Skipping this component.');
                ica_success = false;
                continue;
            end
        end
        
        if ~ica_success
            continue;
        end

        %% ========== ENHANCED NORMALIZATION ==========
        for m = 1:2
            % Robust variance estimation (using MAD)
            signal_mad = mad(y{n}(m,:), 1);
            robust_std = signal_mad * 1.4826;  % Convert MAD to STD estimate
            den = 10 * robust_std + eps;
            y{n}(m,:) = y{n}(m,:) / den;
        end

        %% ========== MASK ESTIMATION ==========
        if dis, disp('  -> Estimating masks...'); end
        
        if evalu
            [newX{1}(1,:),newX{1}(2,:),newX{2}(1,:),newX{2}(2,:), msk{1}, msk{2}] = ...
                applymasks(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', ...
                          mask{n}, fs, th, NFFT, WINDOW, NOVERLAP, cmR, cmL);
        else
            [newX{1}(1,:),newX{1}(2,:),newX{2}(1,:),newX{2}(2,:), msk{1}, msk{2}] = ...
                applymasks(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', ...
                          mask{n}, fs, th, NFFT, WINDOW, NOVERLAP);
        end

        %% ========== ENHANCED STOPPING CRITERIA ==========
        for m = 1:2
            condi = oneortwo_cond(newX{m}(1,:), newX{m}(2,:), fs);
            est = enerstop(newX{m}(1,:), newX{m}(2,:), thE, minpower);
            
            % NEW: Calculate separation quality metric
            sep_quality = calculate_separation_quality(newX{m}, X, NFFT, WINDOW, NOVERLAP);
            separation_quality_history = [separation_quality_history, sep_quality];
            
            if est == 2
                if dis, disp('  -> Rejected: Energy too low'); end
            elseif est == 1
                if dis, disp('  -> Low quality signal detected'); end
                [L,R,enermask{enercnt}] = getfinalmask(X(1,:)', X(2,:)', ...
                    y{n}(1,:)', y{n}(2,:)', mask{n}, fs, th, m, NFFT, WINDOW, NOVERLAP, 1);
                stestr = sprintf('enerstereo%d.wav', enercnt);
                audiowrite(stestr, [L, R], fs);
                if dis, fprintf('  -> Saved: %s\n', stestr); end
                enercnt = enercnt + 1;
            elseif condi > stopthreshold && sep_quality > 0.3  % NEW: Added quality check
                if dis, disp('  -> HIGH QUALITY - Saving as final signal'); end
                [L,R,fmask{finalcnt}] = getfinalmask(X(1,:)', X(2,:)', ...
                    y{n}(1,:)', y{n}(2,:)', mask{n}, fs, 1, m, NFFT, WINDOW, NOVERLAP, 0);
                stestr = sprintf('finalstereo%d.wav', finalcnt);
                audiowrite(stestr, [L, R], fs);
                if dis, fprintf('  -> Saved: %s (Quality: %.2f)\n', stestr, sep_quality); end
                finalcnt = finalcnt + 1;
            else
                xbuffer = {xbuffer{:}, newX{m}};
                maskbuffer = {maskbuffer{:}, msk{m}};
            end
        end
    end % for n

    x = xbuffer;
    mask = maskbuffer;

    %% ========== CHECK FOR COMPLETION ==========
    if isempty(xbuffer)
        if dis, disp('=== No more components to process ==='); end

        % Merge duplicates
        if length(fmask) ~= lastfmasklength
            fmask = multisigcheck(fmask, X(1,:)', X(2,:)', TC1, fs, NFFT, WINDOW, NOVERLAP, numlags);
            lastfmasklength = length(fmask);
            finalcnt = lastfmasklength + 1;
        end

        % Check correlations with low-energy signals
        if ~isempty(enermask{1})
            fmask = nosigcorr(fmask, enermask, X(1,:)', X(2,:)', TC2, fs, NFFT, WINDOW, NOVERLAP, numlags);
        end

        enercnt = 1;
        enermask = {[]};

        if isempty(fmask{1})
            if dis, disp('No signals segregated. Stopping.'); end
            exitcnt = 4; 
            break;
        else
            if dis, disp('Computing remaining mask...'); end
            if evalu
                [L,R,remainingmask] = getremainingmask(X(1,:)', X(2,:)', fmask, fs, NFFT, WINDOW, NOVERLAP, cmL, cmR);
            else
                [L,R,remainingmask] = getremainingmask(X(1,:)', X(2,:)', fmask, fs, NFFT, WINDOW, NOVERLAP);
            end
        end

        % Check if mask changed
        if isequal(lastremmask, remainingmask)
            if dis, disp('Mask unchanged. Stopping.'); end
            exitcnt = 4;
            audiowrite('remaining.wav', [L,R], fs);
            break;
        else
            lastremmask = remainingmask;
        end

        audiowrite('remaining.wav', [L,R], fs);
        
        % Check remaining signal quality
        if enerstop(L', R', thE, minpower, fs) > 0
            if dis, disp('Remaining signal has insufficient energy.'); end
            break;
        elseif oneortwo_cond(L', R', fs) > stopthreshold
            stestr = sprintf('finalstereo%d.wav', finalcnt);
            audiowrite(stestr, [L,R], fs);
            if dis, fprintf('Saved remaining as %s\n', stestr); end
            fmask{finalcnt} = remainingmask;
            break;
        else
            mask = {remainingmask};
            x = {[L'; R']};
        end

        if exitcnt >= 4, break; end
        exitcnt = exitcnt + 1;
    end
end % while

%% ========== FINAL CLEANUP ==========
if isempty(fmask)
    flag = 1; 
    mEL=0; mNR=0; mSNRi=0; mSNRo=0; mSNRx=0;
else
    if length(fmask) ~= lastfmasklength
        fmask = multisigcheck(fmask, X(1,:)', X(2,:)', TC1, fs, NFFT, WINDOW, NOVERLAP, numlags);
    end
    if ~isempty(enermask{1})
        fmask = nosigcorr(fmask, enermask, X(1,:)', X(2,:)', TC2, fs, NFFT, WINDOW, NOVERLAP, numlags);
    end
end

if dis, disp('=== Separation Complete ==='); end

%% ========== EVALUATION ==========
if evalu
    if dis, disp('=== Evaluating Results ==='); end

    clear x mask
    [valL,valR,e1L,e1R,e2L,e2R,q,lbl,cflag] = ...
        comparemasks(fmask, imaskL, imaskR, labelvec(Ns), fs, length(s));

    % Ground truth signals
    for i=1:N
        Xalone(:,:,i) = A(:,i) * S(i,:);
    end

    [PLEL,PREL,PLNR,PRNR,SNRL,SNRR,SNRiLi,SNRiRi,SNRxL,SNRxR] = ...
        calcELNR(e1L,e1R,e2L,e2R,imaskL,imaskR,q,NFFT,WINDOW,NOVERLAP,Xalone,lbl);

    % Save results
    datafile = 'data_enhanced.mat';
    save(datafile, 'PLEL','PREL','PLNR','PRNR','SNRL','SNRR','SNRiL','SNRiR',...
         'lbl','th','stopthresholdini','TC1','TC2','thepow','minpow',...
         'NFFT','winnumber','k','NOVERLAP','Ns','theta','separation_quality_history');

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
    save mydata_enhanced mEL mNR mSNRi mSNRo mSNR mSNRx flag
    
    % Display results
    fprintf('\n=== PERFORMANCE METRICS ===\n');
    fprintf('Mean Error (EL): %.2f%%\n', mEL);
    fprintf('Mean Noise Residual (NR): %.2f%%\n', mNR);
    fprintf('Input SNR: %.2f dB\n', mSNRi);
    fprintf('Output SNR: %.2f dB\n', mSNRo);
    fprintf('SNR Improvement: %.2f dB\n', mSNR);
    fprintf('Number of sources detected: %d (Expected: %d)\n', length(fmask), N);
end

if dis, disp('=== All Done ==='); end

%% ========== NEW HELPER FUNCTION ==========
function quality = calculate_separation_quality(separated, original, NFFT, WINDOW, NOVERLAP)
    % Calculate separation quality based on spectral distance
    % Returns value between 0 (poor) and 1 (excellent)
    
    % Compute spectrograms
    [~,~,~,P_sep] = spectrogram(separated(1,:), WINDOW, NOVERLAP, NFFT);
    [~,~,~,P_orig] = spectrogram(original(1,:), WINDOW, NOVERLAP, NFFT);
    
    % Compute spectral distance (simplified)
    if size(P_sep,2) ~= size(P_orig,2)
        % Resize if needed
        min_len = min(size(P_sep,2), size(P_orig,2));
        P_sep = P_sep(:,1:min_len);
        P_orig = P_orig(:,1:min_len);
    end
    
    % Normalized cross-correlation in frequency domain
    P_sep_norm = abs(P_sep) ./ (sum(abs(P_sep(:))) + eps);
    P_orig_norm = abs(P_orig) ./ (sum(abs(P_orig(:))) + eps);
    
    % Calculate quality metric
    correlation = sum(P_sep_norm(:) .* P_orig_norm(:));
    sparsity = sum(abs(P_sep(:)) > 0.1*max(abs(P_sep(:)))) / numel(P_sep);
    
    quality = 0.6 * (1 - correlation) + 0.4 * sparsity;
    quality = max(0, min(1, quality));
end
