%% main.m
% Ch??ng tr�nh ch�nh: Underdetermined source separation
clear all; close all; clc;
format compact;

dis = 1;
if dis, disp('Initialisation...'); end

rng('default');
rng(1,'twister');

M = 1;
u = 0.5;
N = 2;
th = 1;
stopthresholdini = 3000;
TC1 = 0.1;
TC2 = 0.03;
numlags = 1;
thepow = 20;
minpow = 30;
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

%% === T?o th? m?c l?u k?t qu? ===
result_folder = 'result';
if ~exist(result_folder, 'dir')
    mkdir(result_folder);
end

%% === Load / create sources and stereo mix ===
if evalu
    file_list = {
        'sp1.wav', 
        'sp2.wav',
        'sp3.wav',
        'sp4.wav'
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

    labelvec = {'A','B','C','D'};

    NS = 1:size(s,2);
    Ns = zeros(1,N);
    for i=1:N
        rnd = ceil(rand(1)*(size(s,2)-i+1));
        Ns(i) = NS(rnd);
        NS(rnd) = [];
    end
    N = length(Ns);
    S = s(:,Ns)';

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
    audiowrite(fullfile(result_folder,'stereomix.wav'), X', fs);
else
    if ~exist(fullfile(result_folder,'stereomix.wav'),'file')
        error('Kh�ng t�m th?y stereo mix: result/stereomix.wav');
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

    %% === L?u h�nh th? c�ng (b? % n?u mu?n l?u) ===
    % figure; imshow(cmL); title('Color Mask Left');
    % exportgraphics(gcf, fullfile(result_folder, 'color_mask_left.png'), 'Resolution', 300);
    % figure; imshow(cmR); title('Color Mask Right');
    % exportgraphics(gcf, fullfile(result_folder, 'color_mask_right.png'), 'Resolution', 300);
end

x = {X};
mask = {[]};
fmask = {[]};
enermask = {[]};
delete_me_again = sg(X(1,:),NFFT,fs,WINDOW,NOVERLAP);
lastremmask = zeros(size(delete_me_again));
clear delete_me_again

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
                enercnt = enercnt + 1;

                %% === L?u h�nh th? c�ng n?u h�m c� hi?n th? h�nh ===
                % exportgraphics(gcf, fullfile(result_folder, sprintf('enermask_%d.png', enercnt)), 'Resolution', 300);

            elseif condi > stopthreshold
                [L,R,fmask{finalcnt}] = getfinalmask(X(1,:)', X(2,:)', y{n}(1,:)', y{n}(2,:)', mask{n}, fs, 1, m, NFFT, WINDOW, NOVERLAP, 0);
                stestr = fullfile(result_folder, sprintf('finalstereo%d.wav', finalcnt));
                audiowrite(stestr, [L, R], fs);
                finalcnt = finalcnt + 1;

                % exportgraphics(gcf, fullfile(result_folder, sprintf('finalmask_%d.png', finalcnt)), 'Resolution', 300);
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

    %% === L?u h�nh th? c�ng khi so s�nh m?t n? ===
    % figure; imagesc(valL); title('Left Mask Comparison'); colorbar;
    % exportgraphics(gcf, fullfile(result_folder, 'compare_mask_left.png'), 'Resolution', 300);
    % figure; imagesc(valR); title('Right Mask Comparison'); colorbar;
    % exportgraphics(gcf, fullfile(result_folder, 'compare_mask_right.png'), 'Resolution', 300);

    for i=1:N
        Xalone(:,:,i) = A(:,i) * S(i,:);
    end

    [PLEL,PREL,PLNR,PRNR,SNRL,SNRR,SNRiLi,SNRiRi,SNRxL,SNRxR] = ...
        calcELNR(e1L,e1R,e2L,e2R,imaskL,imaskR,q,NFFT,WINDOW,NOVERLAP,Xalone,lbl);

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
end

if dis, disp('Done'); end
