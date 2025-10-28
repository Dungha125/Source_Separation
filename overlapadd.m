function x = overlapadd(S, window, noverlap)
% OVERLAPADD Tong hop lai tin hieu tu STFT bang phuong phap overlap-add
%
% Inputs:
%   S        - Ma tran STFT (nfreq x nframes)
%   window   - Cua so phan tich
%   noverlap - So mau chong lap
%
% Output:
%   x        - Tin hieu thoi gian da tong hop

    [nfreq, nframes] = size(S);
    nfft = 2*(nfreq-1);
    hop = length(window) - noverlap;
    
    % Khoi tao tin hieu dau ra
    x = zeros(1, (nframes-1)*hop + length(window));
    
    % Overlap-add cho tung frame
    for i = 1:nframes
        % IFFT de chuyen ve mien thoi gian
        % Tao lai phan tan so am bang cach lay conjugate va dao nguoc
        full_spectrum = [S(:,i); conj(flipud(S(2:end-1,i)))];
        frame = real(ifft(full_spectrum, nfft));
        
        % Chi lay phan dau bang do dai cua so
        frame = frame(1:length(window)) .* window;
        
        % Cong don vao tin hieu dau ra
        start_idx = (i-1)*hop + 1;
        end_idx = start_idx + length(window) - 1;
        x(start_idx:end_idx) = x(start_idx:end_idx) + frame';
    end
end