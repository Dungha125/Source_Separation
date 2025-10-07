function fmask = nosigcorr(fmask, enermask, wholesigL, wholesigR, tc, Fs, NFFT, WINDOW, NOVERLAP, numlags)
%==========================================================================
%  HÀM: nosigcorr
%--------------------------------------------------------------------------
%  M?c ?ích:
%   - Ki?m tra t??ng quan gi?a các tín hi?u n?ng l??ng th?p (enermask)
%     và các tín hi?u chính (fmask).
%   - N?u tín hi?u n?ng l??ng th?p t??ng quan cao v?i tín hi?u chính,
%     thì h?p nh?t (merge) m?t n? t?n s? t??ng ?ng.
%   - Sau ?ó tái t?o l?i tín hi?u âm thanh t? m?t n? ?ã h?p nh?t.
%
%  Tham s? vào:
%   + fmask:     Cell array ch?a m?t n? t?n s? c?a các tín hi?u cu?i.
%   + enermask:  Cell array ch?a m?t n? n?ng l??ng th?p.
%   + wholesigL: Tín hi?u g?c kênh trái (Left).
%   + wholesigR: Tín hi?u g?c kênh ph?i (Right).
%   + tc:        Ng??ng t??ng quan (threshold coefficient).
%   + Fs:        T?n s? l?y m?u (sample rate).
%   + NFFT:      Kích th??c FFT (dùng trong phân tích ph?).
%   + WINDOW:    C?a s? s? d?ng trong FFT.
%   + NOVERLAP:  S? m?u ch?ng l?n gi?a các khung.
%   + numlags:   S? l??ng ?? tr? khi tính hàm t??ng quan chéo.
%
%  Giá tr? tr? v?:
%   + fmask:     B? m?t n? sau khi h?p nh?t.
%
%  Ghi chú:
%   MATLAB R2015 tr? ?i ?ã lo?i b? 'wavread' và 'wavwrite',
%   thay vào ?ó dùng 'audioread' và 'audiowrite'.
%==========================================================================

% --- Tính ph? tín hi?u g?c kênh trái & ph?i ---
ywholeL = sg(wholesigL, NFFT, Fs, WINDOW, NOVERLAP);
ywholeR = sg(wholesigR, NFFT, Fs, WINDOW, NOVERLAP);

Nsig = length(fmask);     % S? tín hi?u chính
Nesig = length(enermask); % S? tín hi?u n?ng l??ng th?p

% --- ??c các file finalstereo{i}.wav (tín hi?u chính) ---
for i = 1:Nsig
    str = strcat('finalstereo', num2str(i), '.wav');
    [s, fs_wav] = audioread(str);

    % ??m b?o t?n s? m?u ??ng nh?t
    if fs_wav ~= Fs
        s = resample(s, Fs, fs_wav);
    end

    % L?u kênh trái/ph?i
    x(:, i) = s(:, 1);  % Ch? l?y kênh trái
    y(:, i) = s(:, 2);  % Kênh ph?i
end

% --- Chu?n b? bi?n ---
todelete = [];   % Danh sách enermask b? g?p
changed  = [];   % Danh sách fmask ?ã thay ??i
xco = zeros(Nsig, Nesig); % Ma tr?n l?u h? s? t??ng quan

% --- Tính t??ng quan gi?a t?ng fmask và enermask ---
for i = 1:Nsig
    for j = 1:Nesig
        str = strcat('enerstereo', num2str(j), '.wav');
        [s, fs_wav] = audioread(str);

        % ??m b?o cùng t?n s?
        if fs_wav ~= Fs
            s = resample(s, Fs, fs_wav);
        end

        xe = s(:, 1); % L?y kênh trái
        % Tính h? s? t??ng quan chéo t?i ?a
        xco(i, j) = max(xcorr(x(:, i), xe, numlags, 'coeff'));
    end
end

% --- Ki?m tra ?i?u ki?n t??ng quan và h?p nh?t m?t n? ---
for i = 1:Nsig
    for j = 1:Nesig
        % N?u h? s? t??ng quan cao và là giá tr? l?n nh?t ? c?t ?ó
        if xco(i, j) > tc && xco(i, j) == max(xco(:, j))
            if issparse(enermask{j})
                enermask{j} = full(enermask{j});
            end

            % H?p nh?t m?t n? fmask{i} và enermask{j}
            fmask{i} = fmask{i} + enermask{j} - fmask{i} .* enermask{j};

            % Tái t?o l?i tín hi?u
            yL = ywholeL .* fmask{i};
            yR = ywholeR .* fmask{i};

            % Ghi nh?n các thay ??i
            todelete = [todelete, j];
            changed  = [changed, i];

            % Dùng hàm ngh?ch ??o ph? ?? tái t?o tín hi?u th?i gian
            x(:, i) = invspecgram(yL, NFFT, Fs, WINDOW, NOVERLAP);
            y(:, i) = invspecgram(yR, NFFT, Fs, WINDOW, NOVERLAP);

            % Chuy?n l?i v? sparse ?? ti?t ki?m b? nh?
            enermask{j} = sparse(enermask{j});
        end
    end
end

% --- Hi?n th? thông tin g? l?i ---
disp('Ma tr?n t??ng quan xco:');
disp(xco);
disp('Các ch? s? enermask b? xóa:');
disp(todelete);
disp('Các fmask ?ã thay ??i:');
disp(changed);

% --- Ghi l?i các file finalstereo{i}.wav ?ã c?p nh?t ---
% --- Ghi l?i các file finalstereo{i}.wav ?ã c?p nh?t ---
for i = 1:length(fmask)
    stestr = strcat('finalstereo', int2str(i), '.wav');
    audiowrite(stestr, [x(:, i), y(:, i)], Fs);
end

% --- Xóa các file n?ng l??ng th?p ?ã h?p nh?t ---
for j = 1:Nesig
    stestr = strcat('enerstereo', int2str(j), '.wav');
    if exist(stestr, 'file') == 2   % <== ki?m tra file t?n t?i
        delete(stestr);
    end
end
end

