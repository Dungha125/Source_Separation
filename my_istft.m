function x_rec = my_istft(S, win, noverlap, nfft, target_len)
    % Ham Inverse STFT rieng biet
    [n_freq, n_frames] = size(S);
    hop_size = length(win) - noverlap;
    
    if mod(nfft, 2) == 0
        S_full = [S; conj(flipud(S(2:end-1, :)))];
    else
        S_full = [S; conj(flipud(S(2:end, :)))];
    end
    
    x_frames = real(ifft(S_full, nfft));
    x_len_est = (n_frames-1) * hop_size + nfft;
    x_rec = zeros(1, x_len_est);
    win_norm = zeros(1, x_len_est);
    
    for i = 1:n_frames
        idx = (i-1)*hop_size + 1;
        ii = idx : idx + nfft - 1;
        x_rec(ii) = x_rec(ii) + x_frames(:, i)';
        win_norm(ii) = win_norm(ii) + win';
    end
    
    win_norm(win_norm < 1e-6) = 1; 
    x_rec = x_rec ./ win_norm;
    
    if nargin >= 5
        x_rec = x_rec(1:min(length(x_rec), target_len));
    end
end