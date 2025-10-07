%% main.m
% Ch??ng trình chính: Underdetermined source separation
% T??ng thích MATLAB 2016
clear all; close all; clc;
format compact;

dis = 1;
if dis, disp('Initialisation...'); end

% C?u hình ng?u nhiên ?? l?p l?i ???c k?t qu?
rng('default');
rng(1,'twister');

% Tham s? chung
M = 2;            % s? microphone
u = 0.5;          % cardioid control
N = 3;            % s? ngu?n ?? t?o mixture khi evalu=1
th = 1;
stopthresholdini = 3000;
TC1 = 0.1;
TC2 = 0.03;
numlags = 1;
thepow = 20;
minpow = 30;

evalu = 1; % ch?y ch? ?? tách + ?ánh giá

% STFT params
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

%% === Load / create sources and stereo mix ===
if evalu
    s = [];
    [s(:,1),fs]  = audioread('ukma.wav'); % A
    [s(:,2),~]   = audioread('frma.wav'); % B
    [s(:,3),~]   = audioread('itfe.wav'); % C
    [s(:,4),~]   = audioread('cnfe.wav'); % D
    [s(:,5),~]   = audioread('rufe.wav'); % E
    [s(:,6),~]   = audioread('gema.wav'); % F
    [s(:,7),~]   = audioread('nlma.wav'); % G
    [s(:,8),~]   = audioread('jpfe.wav'); % H
    [s(:,9),~]   = audioread('brfe.wav'); % I
    [s(:,10),~]  = audioread('esma.wav'); % J
    [s(:,11),~]  = audioread('dkma.wav'); % K
    [s(:,12),~]  = audioread('ukfe.wav'); % L
    labelvec = {'A','B','C','D','E','F','G','H','I','J','K','L'};

    % Random pick N sources from pool
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

    % Mixing with two cardioid microphones
    A = calcA(theta,u);   % b?n ph?i có file calcA.m
    X = A * S;            % 2 x T
    audiowrite('stereomix.wav', X', fs);
else
    if ~exist('stereomix.wav','file')
        error('Không tìm th?y stereo mix: stereomix.wav');
    end
    [X,fs] = audioread('stereomix.wav');
    X = X';
end

% estimate energy thresholds
powpow = 10*log10((sum(X(1,:).^2)+sum(X(2,:).^2))/(2*size(X,2)));
thE = powpow - thepow;
minpower = powpow - minpow;

% Create ideal masks if evalu
if evalu
    for i=1:N
        vd = zeros(1,N); vd(i)=1; ivd=(vd-1)*(-1);
        [imaskL{i},imaskR{i},SNRiL(i),SNRiR(i)] = idealmask(A*diag(vd)*S, A*diag(ivd)*S, fs, NFFT, WINDOW, NOVERLAP);
    end
    [cmL,cmR] = colorimask(imaskL,imaskR,fs);
end

% buffers
x = {X};
mask = {[]};
fmask = {[]};
enermask = {[]};
delete_me_again = sg(X(1,:),NFFT,fs,WINDOW,NOVERLAP);
lastremmask = zeros(size(delete_me_again));
clear delete_me_again

% counters
countmax = 30;
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
        % --- PCA whitening (prewhiten) to improve ICA stability ---
        % x{n} expected size 2 x T
        Xn = x{n};
        % remove mean
        Xm = Xn - repmat(mean(Xn,2), 1, size(Xn,2));
        % covariance
        C = cov(Xm.');
        [E,D] = eig(C);
        % avoid tiny eigenvalues
        d = diag(D);
        d(d<=0) = eps;
        % whitening matrix
        Wwhite = inv(sqrt(D)) * E';
        Xwhite = Wwhite * Xm; % 2 x T whitened

        % Pass whitened data to ICA
        try
            [y{n}, Aest] = icaML(Xwhite); % b?n c?n file icaML.m
            % y{n} returned is 2 x T
            % Transform back to original space if needed (not necessary for masks)
        catch ME
            warning('ICA failed, using unwhitened data: %s', ME.message);
            [y{n}, Aest] = icaML(Xn);
        end

        % normalization
        for m=1:2
            den = 10*sqrt(var(y{n}(m,:))) + eps;
            y{n}(m,:) = y{n}(m,:) / den;
        end

        % estimate mask
        if dis, disp('Finding mask...'); end
        if evalu
            [newX{1}(1,:),newX{1}(2,:),newX{2}(1,:),newX{2}(2,:), msk{1}, msk{2}] = ...
                applymasks(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', mask{n}, fs, th, NFFT, WINDOW, NOVERLAP, cmR, cmL);
        else
            [newX{1}(1,:),newX{1}(2,:),newX{2}(1,:),newX{2}(2,:), msk{1}, msk{2}] = ...
                applymasks(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', mask{n}, fs, th, NFFT, WINDOW, NOVERLAP);
        end

        % binaural stopping criteria
        for m=1:2
            condi = oneortwo_cond(newX{m}(1,:), newX{m}(2,:), fs);
            est = enerstop(newX{m}(1,:), newX{m}(2,:), thE, minpower);
            if est == 2
                if dis, disp('Not a speech signal - too low energy'); end
            elseif est == 1
                if dis, disp('Not a good quality speech signal'); end
                [L,R,enermask{enercnt}] = getfinalmask(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', mask{n}, fs, th, m, NFFT, WINDOW, NOVERLAP, 1);
                stestr = sprintf('enerstereo%d.wav', enercnt);
                if dis, disp(stestr); end
                audiowrite([stestr], [L, R], fs);
                enercnt = enercnt + 1;
            elseif condi > stopthreshold
                if dis, disp('Save as final signal...'); end
                [L,R,fmask{finalcnt}] = getfinalmask(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', mask{n}, fs, 1, m, NFFT, WINDOW, NOVERLAP, 0);
                stestr = sprintf('finalstereo%d.wav', finalcnt);
                if dis, disp(stestr); end
                audiowrite([stestr], [L, R], fs);
                finalcnt = finalcnt + 1;
            else
                xbuffer = {xbuffer{:}, newX{m}};
                maskbuffer = {maskbuffer{:}, msk{m}};
            end
        end
    end % for n

    x = xbuffer;
    mask = maskbuffer;

    if isempty(xbuffer)
        if dis, disp('Stopping separation algorithm'); end

        % merge duplicates and check correlations
        if length(fmask) ~= lastfmasklength
            fmask = multisigcheck(fmask, X(1,:)', X(2,:)', TC1, fs, NFFT, WINDOW, NOVERLAP, numlags);
            lastfmasklength = length(fmask);
            finalcnt = lastfmasklength + 1;
        end

        if ~isempty(enermask{1})
            fmask = nosigcorr(fmask, enermask, X(1,:)', X(2,:)', TC2, fs, NFFT, WINDOW, NOVERLAP, numlags);
        end

        enercnt = 1;
        enermask = {[]};

        if isempty(fmask{1})
            if dis, disp('No signals segregated. Stopping separation procedure.'); end
            exitcnt = 4; break;
        else
            if dis, disp('Finding remaining mask...'); end
            if evalu
                [L,R,remainingmask] = getremainingmask(X(1,:)', X(2,:)', fmask, fs, NFFT, WINDOW, NOVERLAP, cmL, cmR);
            else
                [L,R,remainingmask] = getremainingmask(X(1,:)', X(2,:)', fmask, fs, NFFT, WINDOW, NOVERLAP);
            end
        end

        if isequal(lastremmask, remainingmask)
            if dis, disp('No changes. Stopping separation procedure.'); end
            exitcnt = 4;
            audiowrite('remaining.wav', [L,R], fs);
            break;
        else
            lastremmask = remainingmask;
        end

        audiowrite('remaining.wav', [L,R], fs);
        condi = oneortwo_cond(newX{m}(1,:), newX{m}(2,:), fs);
        if enerstop(L', R', thE, minpower, fs) > 0
            if dis, disp('Not a speech signal'); end
            break;
        elseif condi > stopthreshold
            stestr = sprintf('finalstereo%d.wav', finalcnt);
            audiowrite(stestr, [L,R], fs);
            if dis, disp(['save remaining as ' stestr]); end
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

if isempty(fmask)
    flag = 1; mEL=0; mNR=0; mSNRi=0; mSNRo=0; mSNRx=0;
else
    if length(fmask) ~= lastfmasklength
        fmask = multisigcheck(fmask, X(1,:)', X(2,:)', TC1, fs, NFFT, WINDOW, NOVERLAP, numlags);
    end
    if ~isempty(enermask{1})
        fmask = nosigcorr(fmask, enermask, X(1,:)', X(2,:)', TC2, fs, NFFT, WINDOW, NOVERLAP, numlags);
    end
end

if dis, disp('Separation done.'); end

%% === Evaluation & save results (n?u evalu) ===
if evalu
    if dis, disp('Evaluation of outputs'); end

    clear x mask
    [valL,valR,e1L,e1R,e2L,e2R,q,lbl,cflag] = comparemasks(fmask, imaskL, imaskR, labelvec(Ns), fs, length(s));

    % standalone signals (ground truth)
    for i=1:N
        Xalone(:,:,i) = A(:,i) * S(i,:);
    end

    [PLEL,PREL,PLNR,PRNR,SNRL,SNRR,SNRiLi,SNRiRi,SNRxL,SNRxR] = ...
        calcELNR(e1L,e1R,e2L,e2R,imaskL,imaskR,q,NFFT,WINDOW,NOVERLAP,Xalone,lbl);

    datafile = 'data.mat';
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
    save mydata mEL mNR mSNRi mSNRo mSNR mSNRx flag
end

if dis, disp('Done'); end
