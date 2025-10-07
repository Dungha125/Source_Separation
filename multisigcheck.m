function fmask = multisigcheck(fmask, L, R, TC1, fs, NFFT, WINDOW, NOVERLAP, numlags)
% multisigcheck - Ki?m tra xem m?t m?t n? th?i-t?n c� ch?a nhi?u ngu?n hay kh�ng
%
% INPUT:
%   fmask     - m?t n? th?i-t?n (ma tr?n k�ch th??c freq x time) ho?c vector (tr??ng h?p ??c bi?t)
%   L, R      - t�n hi?u k�nh tr�i v� ph?i (??u l� vector h�ng ho?c c?t)
%   TC1       - ng??ng (tham s? nh?, v� d? 0.1). H�m s? so s�nh v?i ?? l?ch ILD (dB)
%   fs        - sampling rate
%   NFFT      - k�ch th??c FFT d�ng cho spectrogram
%   WINDOW    - c?a s? (window) d�ng cho spectrogram (d?ng vector ho?c chi?u d�i)
%   NOVERLAP  - s? m?u ch?ng l?p gi?a c�c c?a s?
%   numlags   - (kh�ng d�ng ? b?n ??n gi?n n�y, gi? ?? t??ng th�ch)
%
% OUTPUT:
%   fmask     - m?t n? ?� ???c c?p nh?t (n?u ph�t hi?n ?a ngu?n th� gi?m tr?ng s? v�ng t??ng ?ng)
%
% Ghi ch�: H�m kh�ng m? h?p tho?i, ch? x? l� d? li?u truy?n v�o.

% ----- chu?n h�a d?ng vector t�n hi?u -----
if size(L,1) > 1 && size(L,2) == 1
    L = L(:)'; % h�ng
elseif size(L,2) > size(L,1)
    L = L(1,:); % ?� l� h�ng
else
    L = L(:)'; % fallback
end
if size(R,1) > 1 && size(R,2) == 1
    R = R(:)'; % h�ng
elseif size(R,2) > size(R,1)
    R = R(1,:);
else
    R = R(:)';
end

% ----- T�nh spectrogram cho 2 k�nh -----
% N?u WINDOW truy?n v�o l� chi?u (s?) th� t?o window t??ng ?ng
if isscalar(WINDOW)
    winVec = hann(WINDOW);
else
    winVec = WINDOW;
end

[S_L, F, T] = spectrogram(L, winVec, NOVERLAP, NFFT, fs);
[S_R, ~, ~] = spectrogram(R, winVec, NOVERLAP, NFFT, fs);

% K�ch th??c TF
[fn, tn] = size(S_L);

% ----- N?u fmask c� k�ch th??c gi?ng spectrogram th� x? l� theo � TF -----
if isequal(size(fmask), [fn, tn])
    maskTF = fmask ~= 0; % ch?n � th?i-t?n ?ang n?m trong mask
    if ~any(maskTF(:))
        % kh�ng c� � n�o ???c mask: tr? v? kh�ng thay ??i
        return;
    end

    % L?y magnitude tr�n � ???c mask
    magL = abs(S_L(maskTF));
    magR = abs(S_R(maskTF));

    % T�nh ILD (dB) cho c�c � mask
    ILD = 20*log10( (magL + eps) ./ (magR + eps) );

    % Th?ng k�: ?? l?ch chu?n ILD (dB) v� d?i (range)
    stdILD = std(ILD);
    rangeILD = max(ILD) - min(ILD);

    % Ng??ng quy v? dB: n?u ?? l?ch chu?n l?n h?n 3 dB (t�y nghi?m) -> ?a ngu?n
    % D�ng TC1 ?? ?i?u ch?nh: TC1 nh? -> ng??ng cao h?n; TC1 m?c ??nh 0.1
    % Ch�ng ta thi?t l?p ng??ng_dB = 3 + 20*TC1 (linh ho?t)
    thresh_dB = 3 + 20 * TC1;

    if stdILD > thresh_dB || rangeILD > 2*thresh_dB
        % Ph�t hi?n ?a ngu?n trong v�ng mask -> gi?m tr?ng s? mask ? c�c � ?�
        if nargout > 0
            fmask(maskTF) = fmask(maskTF) * 0.5;
        else
            fmask(maskTF) = fmask(maskTF) * 0.5;
        end
        if isscalar(TC1)
            disp(['[multisigcheck] Ph�t hi?n kh? n?ng ?a ngu?n trong v�ng mask � gi?m tr?ng s? (stdILD=', num2str(stdILD,'%.2f'), ' dB).']);
        else
            disp(['[multisigcheck] Ph�t hi?n kh? n?ng ?a ngu?n � gi?m tr?ng s? v�ng mask.']);
        end
    else
        % Kh�ng th?y b?ng ch?ng ?a ngu?n -> kh�ng thay ??i
        if nargout == 0
            % nothing
        end
    end

    return;
end

% ----- N?u fmask kh�ng ph?i ma tr?n th?i-t?n (fallback): d�ng t??ng quan ch�o -----
% ?�y l� t�nh hu?ng ph?: th?c hi?n ki?m tra to�n c?c tr�n t�n hi?u
xc = xcorr(L, R, numlags, 'coeff');
[maxval, idx] = max(abs(xc));
lag = idx - (numlags + 1);
% N?u t??ng quan c?c ??i th?p (v� d? < 0.6) -> kh? n?ng nhi?u ngu?n (k�m ??ng pha)
if maxval < (1 - TC1)
    % n?u fmask l� vector ho?c cell ch?a mask d?ng r?i, gi?m tr?ng s?
    try
        maskIdx = fmask ~= 0;
        fmask(maskIdx) = fmask(maskIdx) * 0.5;
    catch
        % n?u fmask kh�ng th? c?p nh?t th� b? qua
    end
    disp(['[multisigcheck] Fallback: t??ng quan ch�o th?p (', num2str(maxval,'%.2f'), ') -> gi?m tr?ng s? mask.']);
else
    % t??ng quan cao -> gi? nguy�n
end

end
