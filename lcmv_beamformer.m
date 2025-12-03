function y = lcmv_beamformer(X, fs, look_direction, constraint_directions)
% LCMV_BEAMFORMER - Linearly Constrained Minimum Variance Beamformer
% 
% Inputs:
%   X - Tín hiệu đầu vào (M x N): M = số microphone, N = số mẫu
%   fs - Tần số lấy mẫu (Hz)
%   look_direction - Hướng nhìn chính (góc độ)
%   constraint_directions - Các hướng cần null (góc độ, optional)
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

% Tính STFT
[S1, F, T] = spectrogram(X(1, :), WINDOW, NOVERLAP, NFFT, fs);
num_freqs = size(S1, 1);
num_frames = size(S1, 2);
S = zeros(M, num_freqs, num_frames);
S(1, :, :) = S1;
for m = 2:M
    [S_temp, ~, ~] = spectrogram(X(m, :), WINDOW, NOVERLAP, NFFT, fs);
    S(m, :, :) = S_temp;
end

% Tính covariance matrix
R_xx = zeros(M, M, num_freqs);
for f_idx = 1:num_freqs
    S_f = squeeze(S(:, f_idx, :));
    R_xx(:, :, f_idx) = (S_f * S_f') / num_frames;
end

% Constraint matrix
theta_look = deg2rad(look_direction);
if nargin < 4 || isempty(constraint_directions)
    constraint_directions = [];  % Chỉ có constraint cho look direction
end

% Tính steering vectors
a_look = zeros(M, num_freqs);
for f_idx = 1:num_freqs
    omega = 2 * pi * freqs(f_idx);
    for m = 1:M
        mic_pos = (m - 1) * d;
        phase = omega * mic_pos * sin(theta_look) / c;
        a_look(m, f_idx) = exp(1j * phase);
    end
end

% Constraint: C' * w = f
% C = [a_look, a_null1, a_null2, ...]
% f = [1, 0, 0, ...]
num_constraints = 1 + length(constraint_directions);
C = zeros(M, num_constraints, num_freqs);
f = [1; zeros(num_constraints-1, 1)];

for f_idx = 1:num_freqs
    C(:, 1, f_idx) = a_look(:, f_idx);
    
    % Null constraints
    for null_idx = 1:length(constraint_directions)
        theta_null = deg2rad(constraint_directions(null_idx));
        omega = 2 * pi * freqs(f_idx);
        a_null = zeros(M, 1);
        for m = 1:M
            mic_pos = (m - 1) * d;
            phase = omega * mic_pos * sin(theta_null) / c;
            a_null(m) = exp(1j * phase);
        end
        C(:, 1 + null_idx, f_idx) = a_null;
    end
end

% Tính LCMV weights
Y = zeros(num_freqs, num_frames);
for f_idx = 1:num_freqs
    R_f = squeeze(R_xx(:, :, f_idx));
    C_f = squeeze(C(:, :, f_idx));
    
    % LCMV: w = R^(-1) * C * (C' * R^(-1) * C)^(-1) * f
    R_inv = inv(R_f + eye(M) * 1e-6);  % Regularization
    w = R_inv * C_f * inv(C_f' * R_inv * C_f) * f;
    
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

