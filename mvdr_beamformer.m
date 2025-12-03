function y = mvdr_beamformer(X, fs, look_direction, noise_cov)
% MVDR_BEAMFORMER - Minimum Variance Distortionless Response Beamformer
% 
% Inputs:
%   X - Tín hiệu đầu vào (M x N): M = số microphone, N = số mẫu
%   fs - Tần số lấy mẫu (Hz)
%   look_direction - Hướng nhìn (góc độ)
%   noise_cov - Ma trận covariance của nhiễu (optional)
%
% Output:
%   y - Tín hiệu đã được beamforming

[M, N] = size(X);

% Khoảng cách giữa các microphone
d = 0.05;  % mét
c = 343;   % Tốc độ âm thanh (m/s)

% Chuyển sang miền tần số
NFFT = 2048;
NOVERLAP = floor(NFFT * 0.75);
WINDOW = hamming(NFFT);

% Tính STFT cho mỗi microphone
[S1, F, T] = spectrogram(X(1, :), WINDOW, NOVERLAP, NFFT, fs);
num_freqs = size(S1, 1);
num_frames = size(S1, 2);
S = zeros(M, num_freqs, num_frames);
S(1, :, :) = S1;
for m = 2:M
    [S_temp, ~, ~] = spectrogram(X(m, :), WINDOW, NOVERLAP, NFFT, fs);
    S(m, :, :) = S_temp;
end

% Steering vector
theta_look = deg2rad(look_direction);

% Tính steering vector cho mỗi tần số
a = zeros(M, num_freqs);
for f_idx = 1:num_freqs
    omega = 2 * pi * freqs(f_idx);
    for m = 1:M
        mic_pos = (m - 1) * d;
        phase = omega * mic_pos * sin(theta_look) / c;
        a(m, f_idx) = exp(1j * phase);
    end
end

% Tính covariance matrix (nếu không có noise_cov)
if nargin < 4 || isempty(noise_cov)
    % Ước tính từ dữ liệu
    R_xx = zeros(M, M, num_freqs);
    for f_idx = 1:num_freqs
        S_f = squeeze(S(:, f_idx, :));
        R_xx(:, :, f_idx) = (S_f * S_f') / num_frames;
    end
    R_nn = R_xx;  % Giả định R_xx = R_nn (đơn giản hóa)
else
    R_nn = noise_cov;
end

% Tính MVDR weights
Y = zeros(num_freqs, num_frames);
for f_idx = 1:num_freqs
    a_f = a(:, f_idx);
    R_nn_f = squeeze(R_nn(:, :, f_idx));
    
    % MVDR: w = R_nn^(-1) * a / (a' * R_nn^(-1) * a)
    R_inv = inv(R_nn_f + eye(M) * 1e-6);  % Regularization
    w = (R_inv * a_f) / (a_f' * R_inv * a_f);
    
    % Áp dụng weights
    S_f = squeeze(S(:, f_idx, :));
    Y(f_idx, :) = w' * S_f;
end

% Inverse STFT
y = my_istft(Y, WINDOW, NOVERLAP, NFFT, N);
y = y(1:N);

% Chuẩn hóa
y = y / (max(abs(y)) + eps);

end

