function [PLEL,PREL,PLNR,PRNR,SNRL,SNRR,SNRiL,SNRiR,SNRoL,SNRoR] = ...
    calcELNR(e1L,e1R,e2L,e2R,imaskL,imaskR,q,NFFT,WINDOW,NOVERLAP,Xalone, labelvec)
% calcELNR - tính EL/NR và SNR, l?u file WAV và hi?n th? waveform/spectrogram
% Compatible MATLAB 2016
% input:
%   e1L,e1R,e2L,e2R : cell masks
%   imaskL,imaskR   : ideal masks
%   q               : mapping vector
%   NFFT, WINDOW, NOVERLAP : STFT params
%   Xalone          : ground truth standalone signals (2 x T x N)
%   labelvec        : optional nhãn ngu?n (cell)
%
% outputs: các ch? s? (m?ng)

% ??c microphone mix
[wholesig, Fs] = audioread('stereomix.wav');
if size(wholesig,2) == 1
    wholesig = [wholesig, zeros(size(wholesig))];
end

% n?u có remaining & ideal files, ??c (n?u không, b? qua)
hasRemaining = exist('remaining.wav','file') == 2;
if hasRemaining
    [remsig, Fs_rem] = audioread('remaining.wav');
    if Fs_rem ~= Fs, remsig = resample(remsig, Fs, Fs_rem); end
end

% compute STFT-like spectrogram matrices (sg expects column vector)
ywholeL = sg(wholesig(:,1), NFFT, Fs, WINDOW, NOVERLAP);
ywholeR = sg(wholesig(:,2), NFFT, Fs, WINDOW, NOVERLAP);

numOutputs = length(e1L);
% initialize outputs
PLEL = zeros(1,numOutputs); PREL = zeros(1,numOutputs);
PLNR = zeros(1,numOutputs); PRNR = zeros(1,numOutputs);
SNRL = zeros(1,numOutputs); SNRR = zeros(1,numOutputs);
SNRiL = zeros(1,numOutputs); SNRiR = zeros(1,numOutputs);
SNRoL = zeros(1,numOutputs); SNRoR = zeros(1,numOutputs);

for i=1:numOutputs
    % l?y output file n?u có
    outname = sprintf('finalstereo%d.wav', i);
    if exist(outname,'file')==2
        [O_n, Fs_o] = audioread(outname);
        if Fs_o ~= Fs, O_n = resample(O_n, Fs, Fs_o); end
        if size(O_n,2)==1, O_n = [O_n, zeros(size(O_n))]; end
    else
        % N?u không có file, tái t?o t? mask (n?u có) ho?c dùng zeros
        % T?o tín hi?u b?ng invspecgram t? mask (n?u e1L/e1R t?n t?i)
        try
            speL = ywholeL .* e1L{i};
            speR = ywholeR .* e1R{i};
            left = invspecgram(speL, NFFT, Fs, WINDOW, NOVERLAP);
            right = invspecgram(speR, NFFT, Fs, WINDOW, NOVERLAP);
            O_n = [left(:), right(:)];
        catch
            O_n = zeros(size(wholesig));
        end
    end

    % apply masks to whole spectrogram and invert (estimated components)
    spe1L = ywholeL .* e1L{i};
    spe1R = ywholeR .* e1R{i};
    spe2L = ywholeL .* e2L{i};
    spe2R = ywholeR .* e2R{i};
    spIL_n = ywholeL .* imaskL{q(i)};
    spIR_n = ywholeR .* imaskR{q(i)};

    e1L_n = invspecgram(spe1L, NFFT, Fs, WINDOW, NOVERLAP);
    e1R_n = invspecgram(spe1R, NFFT, Fs, WINDOW, NOVERLAP);
    e2L_n = invspecgram(spe2L, NFFT, Fs, WINDOW, NOVERLAP);
    e2R_n = invspecgram(spe2R, NFFT, Fs, WINDOW, NOVERLAP);
    IL_n   = invspecgram(spIL_n, NFFT, Fs, WINDOW, NOVERLAP);
    IR_n   = invspecgram(spIR_n, NFFT, Fs, WINDOW, NOVERLAP);
    whL    = invspecgram(ywholeL, NFFT, Fs, WINDOW, NOVERLAP);
    whR    = invspecgram(ywholeR, NFFT, Fs, WINDOW, NOVERLAP);

    % tránh chia cho 0
    if sum(IL_n.^2)==0, IL_n = IL_n + eps; end
    if sum(IR_n.^2)==0, IR_n = IR_n + eps; end
    if sum(O_n(:,1).^2)==0, O_n(:,1)=O_n(:,1)+eps; end
    if sum(O_n(:,2).^2)==0, O_n(:,2)=O_n(:,2)+eps; end

    PLEL(i) = sum(e1L_n.^2)/sum(IL_n.^2);
    PREL(i) = sum(e1R_n.^2)/sum(IR_n.^2);
    PLNR(i) = sum(e2L_n.^2)/sum(O_n(:,1).^2);
    PRNR(i) = sum(e2R_n.^2)/sum(O_n(:,2).^2);

    SNRL(i) = 10*log10(sum(IL_n.^2)/sum((IL_n - O_n(:,1)).^2 + eps));
    SNRR(i) = 10*log10(sum(IR_n.^2)/sum((IR_n - O_n(:,2)).^2 + eps));
    SNRiL(i) = 10*log10(sum(IL_n.^2)/sum((IL_n - whL).^2 + eps));
    SNRiR(i) = 10*log10(sum(IR_n.^2)/sum((IR_n - whR).^2 + eps));

    % output vs ground-truth Xalone if provided
    try
        XA1 = Xalone(1,:, q(i));
        XA2 = Xalone(2,:, q(i));
        XA1 = XA1(:)'; XA2 = XA2(:)';
        minlength = min([length(O_n), length(XA1)]);
        XA1 = XA1(1:minlength); XA2 = XA2(1:minlength);
        On1 = O_n(1:minlength,1); On2 = O_n(1:minlength,2);
        XA1n = XA1 / (sqrt(var(XA1))+eps);
        XA2n = XA2 / (sqrt(var(XA2))+eps);
        On1n = On1 / (sqrt(var(On1))+eps);
        On2n = On2 / (sqrt(var(On2))+eps);
        SNRoL(i) = 10*log10(sum(On1n.^2)/sum((XA1n' - On1n).^2 + eps));
        SNRoR(i) = 10*log10(sum(On2n.^2)/sum((XA2n' - On2n).^2 + eps));
    catch
        SNRoL(i)=NaN; SNRoR(i)=NaN;
    end

    % Keep last e1/e2/O_n for plotting/saving outside loop
    if i==1
        last_e1L = e1L_n; last_e1R = e1R_n;
        last_e2L = e2L_n; last_e2R = e2R_n;
        last_O = O_n;
    end
end

% ===============================
% LUU FILE WAV: mic, speaker1, speaker2, separated
% ===============================
disp('Saving WAV files: mic.wav, speaker1.wav, speaker2.wav, separated.wav ...');
try
    audiowrite('mic.wav', wholesig, Fs);
    try
        s1 = [last_e1L(:), last_e1R(:)];
    catch
        s1 = zeros(size(wholesig));
    end
    try
        s2 = [last_e2L(:), last_e2R(:)];
    catch
        s2 = zeros(size(wholesig));
    end
    L = size(wholesig,1);
    if size(s1,1) < L, s1(end+1:L,1)=0; s1(end+1:L,2)=0; end
    if size(s2,1) < L, s2(end+1:L,1)=0; s2(end+1:L,2)=0; end
    if size(last_O,1) < L, last_O(end+1:L,1)=0; last_O(end+1:L,2)=0; end

    audiowrite('speaker1.wav', s1(1:L,:), Fs);
    audiowrite('speaker2.wav', s2(1:L,:), Fs);
    audiowrite('separated.wav', last_O(1:L,:), Fs);
    disp('Saved successfully.');
catch ME
    warning('Khong the luu WAV: %s', ME.message);
end

% ===============================
% VE WAVEFORM CHO CAC FILE
% ===============================
files = {'mic.wav','speaker1.wav','speaker2.wav','separated.wav'};
titles = {'Tin hieu ghi duoc tai micro', 'Tin hieu nguoi noi thu nhat', ...
          'Tin hieu nguoi noi thu hai', 'Tin hieu sau xu ly / tach nguon'};

for i=1:length(files)
    if exist(files{i},'file')==2
        [y, fs_y] = audioread(files{i});
        T = (0:size(y,1)-1)/fs_y;
        figure('Name',titles{i},'NumberTitle','off');
        subplot(2,1,1);
        plot(T, y(:,1));
        title([titles{i} ' - Kenh trai']);
        xlabel('Thoi gian (s)'); ylabel('Bien do');
        grid on;
        subplot(2,1,2);
        plot(T, y(:,2));
        title([titles{i} ' - Kenh phai']);
        xlabel('Thoi gian (s)'); ylabel('Bien do');
        grid on;
        annotation('textbox',[0 0.95 1 0.05],'String',titles{i},'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');
    else
        warning('File %s khong ton tai => bo qua hinh', files{i});
    end
end

% ===============================
% (Tuy chon) ve spectrogram cho separated.wav
% ===============================
if exist('separated.wav','file')==2
    [ysep, fssep] = audioread('separated.wav');
    if size(ysep,2) >= 1
        figure('Name','Spectrogram separated (kenh trai)','NumberTitle','off');
        specgram(ysep(:,1), 1024, fssep);
        title('Spectrogram - separated (kenh trai)');
    end
end

end
