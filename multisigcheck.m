function fmask = multisigcheck(fmask, L, R, TC1, fs, NFFT, WINDOW, NOVERLAP, numlags)
% multisigcheck - Ki?m tra xem m?t m?t n? th?i-t?n có ch?a nhi?u ngu?n hay không
%
% INPUT:
%   fmask     - m?t n? th?i-t?n (ma tr?n kích th??c freq x time) ho?c vector (tr??ng h?p ??c bi?t)
%   L, R      - tín hi?u kênh trái và ph?i (??u là vector hàng ho?c c?t)
%   TC1       - ng??ng (tham s? nh?, ví d? 0.1). Hàm s? so sánh v?i ?? l?ch ILD (dB)
%   fs        - sampling rate
%   NFFT      - kích th??c FFT dùng cho spectrogram
%   WINDOW    - c?a s? (window) dùng cho spectrogram (d?ng vector ho?c chi?u dài)
%   NOVERLAP  - s? m?u ch?ng l?p gi?a các c?a s?
%   numlags   - (không dùng ? b?n ??n gi?n này, gi? ?? t??ng thích)
%
% OUTPUT:
%   fmask     - m?t n? ?ã ???c c?p nh?t (n?u phát hi?n ?a ngu?n thì gi?m tr?ng s? vùng t??ng ?ng)
%
% Ghi chú: Hàm không m? h?p tho?i, ch? x? lý d? li?u truy?n vào.

% ----- chu?n hóa d?ng vector tín hi?u -----
if size(L,1) > 1 && size(L,2) == 1
    L = L(:)'; % hàng
elseif size(L,2) > size(L,1)
    L = L(1,:); % ?ã là hàng
else
    L = L(:)'; % fallback
end
if size(R,1) > 1 && size(R,2) == 1
    R = R(:)'; % hàng
elseif size(R,2) > size(R,1)
    R = R(1,:);
else
    R = R(:)';
end

% ----- Tính spectrogram cho 2 kênh -----
% N?u WINDOW truy?n vào là chi?u (s?) thì t?o window t??ng ?ng
if isscalar(WINDOW)
    winVec = hann(WINDOW);
else
    winVec = WINDOW;
end

[S_L, F, T] = spectrogram(L, winVec, NOVERLAP, NFFT, fs);
[S_R, ~, ~] = spectrogram(R, winVec, NOVERLAP, NFFT, fs);

% Kích th??c TF
[fn, tn] = size(S_L);

% ----- N?u fmask có kích th??c gi?ng spectrogram thì x? lý theo ô TF -----
if isequal(size(fmask), [fn, tn])
    maskTF = fmask ~= 0; % ch?n ô th?i-t?n ?ang n?m trong mask
    if ~any(maskTF(:))
        % không có ô nào ???c mask: tr? v? không thay ??i
        return;
    end

    % L?y magnitude trên ô ???c mask
    magL = abs(S_L(maskTF));
    magR = abs(S_R(maskTF));

    % Tính ILD (dB) cho các ô mask
    ILD = 20*log10( (magL + eps) ./ (magR + eps) );

    % Th?ng kê: ?? l?ch chu?n ILD (dB) và d?i (range)
    stdILD = std(ILD);
    rangeILD = max(ILD) - min(ILD);

    % Ng??ng quy v? dB: n?u ?? l?ch chu?n l?n h?n 3 dB (tùy nghi?m) -> ?a ngu?n
    % Dùng TC1 ?? ?i?u ch?nh: TC1 nh? -> ng??ng cao h?n; TC1 m?c ??nh 0.1
    % Chúng ta thi?t l?p ng??ng_dB = 3 + 20*TC1 (linh ho?t)
    thresh_dB = 3 + 20 * TC1;

    if stdILD > thresh_dB || rangeILD > 2*thresh_dB
        % Phát hi?n ?a ngu?n trong vùng mask -> gi?m tr?ng s? mask ? các ô ?ó
        if nargout > 0
            fmask(maskTF) = fmask(maskTF) * 0.5;
        else
            fmask(maskTF) = fmask(maskTF) * 0.5;
        end
        if isscalar(TC1)
            disp(['[multisigcheck] Phát hi?n kh? n?ng ?a ngu?n trong vùng mask — gi?m tr?ng s? (stdILD=', num2str(stdILD,'%.2f'), ' dB).']);
        else
            disp(['[multisigcheck] Phát hi?n kh? n?ng ?a ngu?n — gi?m tr?ng s? vùng mask.']);
        end
    else
        % Không th?y b?ng ch?ng ?a ngu?n -> không thay ??i
        if nargout == 0
            % nothing
        end
    end

    return;
end

% ----- N?u fmask không ph?i ma tr?n th?i-t?n (fallback): dùng t??ng quan chéo -----
% ?ây là tình hu?ng ph?: th?c hi?n ki?m tra toàn c?c trên tín hi?u
xc = xcorr(L, R, numlags, 'coeff');
[maxval, idx] = max(abs(xc));
lag = idx - (numlags + 1);
% N?u t??ng quan c?c ??i th?p (ví d? < 0.6) -> kh? n?ng nhi?u ngu?n (kém ??ng pha)
if maxval < (1 - TC1)
    % n?u fmask là vector ho?c cell ch?a mask d?ng r?i, gi?m tr?ng s?
    try
        maskIdx = fmask ~= 0;
        fmask(maskIdx) = fmask(maskIdx) * 0.5;
    catch
        % n?u fmask không th? c?p nh?t thì b? qua
    end
    disp(['[multisigcheck] Fallback: t??ng quan chéo th?p (', num2str(maxval,'%.2f'), ') -> gi?m tr?ng s? mask.']);
else
    % t??ng quan cao -> gi? nguyên
end

end
