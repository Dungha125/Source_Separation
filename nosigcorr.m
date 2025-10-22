function fmask = nosigcorr(fmask, enermask, wholesigL, wholesigR, tc, Fs, NFFT, WINDOW, NOVERLAP, numlags, result_folder)
% nosigcorr: ki?m tra t??ng quan gi?a fmask và enermask, merge n?u c?n.
% Tham s? b? sung:
%   result_folder (tùy ch?n) - th? m?c ch?a finalstereo*.wav và enerstereo*.wav
%
% G?i trong main:
%   fmask = nosigcorr(fmask, enermask, X(1,:)', X(2,:)', TC2, fs, NFFT, WINDOW, NOVERLAP, numlags, result_folder);

if nargin < 11 || isempty(result_folder)
    result_folder = '.'; % fallback n?u ch?a truy?n
end

disp('--- Running nosigcorr ---');

% Tính spectrogram g?c (dùng ?? nhân mask và tái t?o)
ywholeL = sg(wholesigL, NFFT, Fs, WINDOW, NOVERLAP);
ywholeR = sg(wholesigR, NFFT, Fs, WINDOW, NOVERLAP);

Nsig = length(fmask);
Nesig = length(enermask);

% Kh?i t?o x,y r?ng (s? m? r?ng khi ??c file)
x = [];
y = [];

% --- ??c các file finalstereo{i}.wav (tín hi?u chính) ---
for i = 1:Nsig
    fname = sprintf('finalstereo%d.wav', i);
    fullpath = fullfile(result_folder, fname);
    if exist(fullpath, 'file') ~= 2
        warning('nosigcorr: file %s không t?n t?i — b? qua (i=%d).', fullpath, i);
        continue;
    end

    [s_sig, fs_wav] = audioread(fullpath);

    % ??m b?o stereo
    if isempty(s_sig)
        warning('nosigcorr: file %s r?ng — b? qua.', fullpath);
        continue;
    end
    if size(s_sig,2) == 1
        s_sig = [s_sig, s_sig];
    end
    if fs_wav ~= Fs
        s_sig = resample(s_sig, Fs, fs_wav);
    end

    len = size(s_sig,1);

    % m? r?ng ma tr?n x và y thích h?p
    if size(x,1) < len
        x(len, max(size(x,2),1)) = 0; %#ok<AGROW>
    end
    if size(y,1) < len
        y(len, max(size(y,2),1)) = 0; %#ok<AGROW>
    end
    if size(x,2) < i
        x(:, i) = 0; %#ok<AGROW>
    end
    if size(y,2) < i
        y(:, i) = 0; %#ok<AGROW>
    end

    x(1:len, i) = s_sig(:,1); % left
    y(1:len, i) = s_sig(:,2); % right
end

% N?u không có d? li?u -> thoát s?m
if isempty(x)
    warning('nosigcorr: Không có tín hi?u finalstereo h?p l? — không x? lý ti?p.');
    return;
end

% --- Chu?n b? ---
todelete = [];
changed  = [];
xco = zeros(Nsig, Nesig);

% --- Tính t??ng quan gi?a t?ng fmask và enermask ---
for i = 1:Nsig
    for j = 1:Nesig
        fname = sprintf('enerstereo%d.wav', j);
        fullpath = fullfile(result_folder, fname);
        if exist(fullpath, 'file') ~= 2
            % không có file enermask j => gi? 0
            xco(i,j) = 0;
            continue;
        end

        [s_e, fs_wav] = audioread(fullpath);
        if isempty(s_e)
            xco(i,j) = 0;
            continue;
        end
        if size(s_e,2) == 1
            s_e = [s_e, s_e];
        end
        if fs_wav ~= Fs
            s_e = resample(s_e, Fs, fs_wav);
        end

        xe = s_e(:,1);
        % ch?n ?? dài chung nh? nh?t
        lenx = min(length(x(:,i)), length(xe));
        if lenx <= 1
            xco(i,j) = 0;
            continue;
        end
        xi = x(1:lenx, i);
        xe = xe(1:lenx);

        xi = xi - mean(xi);
        xe = xe - mean(xe);

        try
            c = xcorr(xi, xe, numlags, 'coeff');
            xco(i, j) = max(c);
        catch
            xco(i, j) = 0;
        end
    end
end

% --- Ki?m tra ?i?u ki?n t??ng quan và h?p nh?t m?t n? ---
for i = 1:Nsig
    for j = 1:Nesig
        if xco(i, j) > tc && xco(i, j) == max(xco(:, j))
            if issparse(enermask{j})
                enermask{j} = full(enermask{j});
            end

            % H?p nh?t mask
            fmask{i} = fmask{i} + enermask{j} - fmask{i} .* enermask{j};

            % Tái t?o tín hi?u t? mask
            yL = ywholeL .* fmask{i};
            yR = ywholeR .* fmask{i};

            todelete = [todelete, j];
            changed  = [changed, i];

            xrec = invspecgram(yL, NFFT, Fs, WINDOW, NOVERLAP);
            yrec = invspecgram(yR, NFFT, Fs, WINDOW, NOVERLAP);

            lenrec = length(xrec);
            if size(x,1) < lenrec
                x(lenrec, max(size(x,2),1)) = 0; %#ok<AGROW>
                y(lenrec, max(size(y,2),1)) = 0; %#ok<AGROW>
            end
            if size(x,2) < i
                x(:, i) = 0; %#ok<AGROW>
            end
            if size(y,2) < i
                y(:, i) = 0; %#ok<AGROW>
            end

            x(1:lenrec, i) = xrec;
            y(1:lenrec, i) = yrec;

            enermask{j} = sparse(enermask{j});
        end
    end
end

% --- Hi?n th? log ---
disp('Ma tran tuong quan xco:'); disp(xco);
disp('Cac chi so enermask bi xoa:'); disp(todelete);
disp('Cac fmask ?a thay ?oi:'); disp(changed);

% --- Ghi l?i các file finalstereo{i}.wav ?ã c?p nh?t ---
for i = 1:size(x,2)
    fname = sprintf('finalstereo%d.wav', i);
    fullpath = fullfile(result_folder, fname);
    xi = x(:, i);
    yi = y(:, i);
    L = max(length(xi), length(yi));
    if length(xi) < L, xi(L,1) = 0; end
    if length(yi) < L, yi(L,1) = 0; end
    audiowrite(fullpath, [xi, yi], Fs);
end

% --- Xóa file enerstereo ?ã h?p nh?t ---
for j = 1:Nesig
    fname = sprintf('enerstereo%d.wav', j);
    fullpath = fullfile(result_folder, fname);
    if exist(fullpath, 'file') == 2
        delete(fullpath);
    end
end

disp('? nosigcorr hoàn t?t.');

end
