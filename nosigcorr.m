function fmask = nosigcorr(fmask, enermask, wholesigL, wholesigR, tc, Fs, NFFT, WINDOW, NOVERLAP, numlags)
%==========================================================================
%  H�M: nosigcorr
%--------------------------------------------------------------------------
%  M?c ?�ch:
%   - Ki?m tra t??ng quan gi?a c�c t�n hi?u n?ng l??ng th?p (enermask)
%     v� c�c t�n hi?u ch�nh (fmask).
%   - N?u t�n hi?u n?ng l??ng th?p t??ng quan cao v?i t�n hi?u ch�nh,
%     th� h?p nh?t (merge) m?t n? t?n s? t??ng ?ng.
%   - Sau ?� t�i t?o l?i t�n hi?u �m thanh t? m?t n? ?� h?p nh?t.
%
%  Tham s? v�o:
%   + fmask:     Cell array ch?a m?t n? t?n s? c?a c�c t�n hi?u cu?i.
%   + enermask:  Cell array ch?a m?t n? n?ng l??ng th?p.
%   + wholesigL: T�n hi?u g?c k�nh tr�i (Left).
%   + wholesigR: T�n hi?u g?c k�nh ph?i (Right).
%   + tc:        Ng??ng t??ng quan (threshold coefficient).
%   + Fs:        T?n s? l?y m?u (sample rate).
%   + NFFT:      K�ch th??c FFT (d�ng trong ph�n t�ch ph?).
%   + WINDOW:    C?a s? s? d?ng trong FFT.
%   + NOVERLAP:  S? m?u ch?ng l?n gi?a c�c khung.
%   + numlags:   S? l??ng ?? tr? khi t�nh h�m t??ng quan ch�o.
%
%  Gi� tr? tr? v?:
%   + fmask:     B? m?t n? sau khi h?p nh?t.
%
%  Ghi ch�:
%   MATLAB R2015 tr? ?i ?� lo?i b? 'wavread' v� 'wavwrite',
%   thay v�o ?� d�ng 'audioread' v� 'audiowrite'.
%==========================================================================

% --- T�nh ph? t�n hi?u g?c k�nh tr�i & ph?i ---
ywholeL = sg(wholesigL, NFFT, Fs, WINDOW, NOVERLAP);
ywholeR = sg(wholesigR, NFFT, Fs, WINDOW, NOVERLAP);

Nsig = length(fmask);     % S? t�n hi?u ch�nh
Nesig = length(enermask); % S? t�n hi?u n?ng l??ng th?p

% --- ??c c�c file finalstereo{i}.wav (t�n hi?u ch�nh) ---
for i = 1:Nsig
    str = strcat('finalstereo', num2str(i), '.wav');
    [s, fs_wav] = audioread(str);

    % ??m b?o t?n s? m?u ??ng nh?t
    if fs_wav ~= Fs
        s = resample(s, Fs, fs_wav);
    end

    % L?u k�nh tr�i/ph?i
    x(:, i) = s(:, 1);  % Ch? l?y k�nh tr�i
    y(:, i) = s(:, 2);  % K�nh ph?i
end

% --- Chu?n b? bi?n ---
todelete = [];   % Danh s�ch enermask b? g?p
changed  = [];   % Danh s�ch fmask ?� thay ??i
xco = zeros(Nsig, Nesig); % Ma tr?n l?u h? s? t??ng quan

% --- T�nh t??ng quan gi?a t?ng fmask v� enermask ---
for i = 1:Nsig
    for j = 1:Nesig
        str = strcat('enerstereo', num2str(j), '.wav');
        [s, fs_wav] = audioread(str);

        % ??m b?o c�ng t?n s?
        if fs_wav ~= Fs
            s = resample(s, Fs, fs_wav);
        end

        xe = s(:, 1); % L?y k�nh tr�i
        % T�nh h? s? t??ng quan ch�o t?i ?a
        xco(i, j) = max(xcorr(x(:, i), xe, numlags, 'coeff'));
    end
end

% --- Ki?m tra ?i?u ki?n t??ng quan v� h?p nh?t m?t n? ---
for i = 1:Nsig
    for j = 1:Nesig
        % N?u h? s? t??ng quan cao v� l� gi� tr? l?n nh?t ? c?t ?�
        if xco(i, j) > tc && xco(i, j) == max(xco(:, j))
            if issparse(enermask{j})
                enermask{j} = full(enermask{j});
            end

            % H?p nh?t m?t n? fmask{i} v� enermask{j}
            fmask{i} = fmask{i} + enermask{j} - fmask{i} .* enermask{j};

            % T�i t?o l?i t�n hi?u
            yL = ywholeL .* fmask{i};
            yR = ywholeR .* fmask{i};

            % Ghi nh?n c�c thay ??i
            todelete = [todelete, j];
            changed  = [changed, i];

            % D�ng h�m ngh?ch ??o ph? ?? t�i t?o t�n hi?u th?i gian
            x(:, i) = invspecgram(yL, NFFT, Fs, WINDOW, NOVERLAP);
            y(:, i) = invspecgram(yR, NFFT, Fs, WINDOW, NOVERLAP);

            % Chuy?n l?i v? sparse ?? ti?t ki?m b? nh?
            enermask{j} = sparse(enermask{j});
        end
    end
end

% --- Hi?n th? th�ng tin g? l?i ---
disp('Ma tr?n t??ng quan xco:');
disp(xco);
disp('C�c ch? s? enermask b? x�a:');
disp(todelete);
disp('C�c fmask ?� thay ??i:');
disp(changed);

% --- Ghi l?i c�c file finalstereo{i}.wav ?� c?p nh?t ---
% --- Ghi l?i c�c file finalstereo{i}.wav ?� c?p nh?t ---
for i = 1:length(fmask)
    stestr = strcat('finalstereo', int2str(i), '.wav');
    audiowrite(stestr, [x(:, i), y(:, i)], Fs);
end

% --- X�a c�c file n?ng l??ng th?p ?� h?p nh?t ---
for j = 1:Nesig
    stestr = strcat('enerstereo', int2str(j), '.wav');
    if exist(stestr, 'file') == 2   % <== ki?m tra file t?n t?i
        delete(stestr);
    end
end
end

