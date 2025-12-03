function [y_enhanced, y_blocked] = gsc_beamformer(X, fs, look_direction, filter_length, mu)
% GSC_BEAMFORMER - Generalized Sidelobe Canceller Beamformer
% 
% Inputs:
%   X - Tín hiệu đầu vào (M x N): M = số microphone, N = số mẫu
%   fs - Tần số lấy mẫu (Hz)
%   look_direction - Hướng nhìn (góc độ, 0 = thẳng, 90 = bên phải)
%   filter_length - Độ dài bộ lọc thích nghi (mặc định: 64)
%   mu - Hệ số học (step size) cho NLMS (mặc định: 0.01)
%
% Outputs:
%   y_enhanced - Tín hiệu đã được tăng cường từ hướng look_direction
%   y_blocked - Tín hiệu đã bị chặn (nhiễu từ các hướng khác)
%
% GSC bao gồm 3 phần:
%   1. Fixed beamformer (Delay-and-Sum)
%   2. Blocking matrix (chặn tín hiệu từ hướng chính)
%   3. Adaptive filter (loại bỏ nhiễu)

[M, N] = size(X);

% Tham số mặc định
if nargin < 3, look_direction = 0; end
if nargin < 4, filter_length = 64; end
if nargin < 5, mu = 0.01; end

% Khoảng cách giữa các microphone (giả định 5cm cho 2 mic)
d = 0.05;  % mét
c = 343;   % Tốc độ âm thanh (m/s)

% Chuyển góc sang radian
theta_look = deg2rad(look_direction);

%% 1. FIXED BEAMFORMER (Delay-and-Sum)
% Tính toán delay cho mỗi microphone
delays_samples = zeros(M, 1);
for m = 1:M
    % Giả định microphone 1 ở gốc, microphone 2 ở vị trí d
    mic_pos = (m - 1) * d;  % Vị trí của microphone thứ m
    delay_time = mic_pos * sin(theta_look) / c;
    delays_samples(m) = round(delay_time * fs);
end

% Điều chỉnh delay (thêm padding để tránh index âm)
max_delay = max(abs(delays_samples));
X_padded = [zeros(M, max_delay), X, zeros(M, max_delay)];

% Áp dụng delay và tổng hợp
y_fixed = zeros(1, N);
for n = 1:N
    idx = n + max_delay;
    sum_val = 0;
    for m = 1:M
        sample_idx = idx - delays_samples(m);
        if sample_idx > 0 && sample_idx <= size(X_padded, 2)
            sum_val = sum_val + X_padded(m, sample_idx);
        end
    end
    y_fixed(n) = sum_val / M;
end

%% 2. BLOCKING MATRIX
% Tạo ma trận chặn để loại bỏ tín hiệu từ hướng look_direction
% Với 2 microphone, blocking matrix đơn giản là: B = [1, -1]
if M == 2
    % Blocking matrix: tạo tín hiệu null ở hướng look_direction
    % Cân bằng theo biên độ và delay
    delay_diff = delays_samples(2) - delays_samples(1);
    
    % Tạo tín hiệu chênh lệch (blocking matrix đơn giản)
    % Với 2 mic, blocking matrix tạo null ở hướng look_direction
    y_blocked = X(2, :) - X(1, :);
else
    % Trường hợp nhiều microphone hơn
    y_blocked = X(2, :) - X(1, :);
end

%% 3. ADAPTIVE FILTER (NLMS)
% Sử dụng NLMS để loại bỏ nhiễu từ y_blocked trong y_fixed

% Khởi tạo bộ lọc
w = zeros(filter_length, 1);
y_enhanced = zeros(1, N);
alpha = 1e-6;  % Regularization parameter cho NLMS

% Tạo ma trận đầu vào cho adaptive filter (sử dụng y_blocked làm reference)
u_buffer = zeros(filter_length, 1);

for n = 1:N
    % Cập nhật buffer
    u_buffer = [y_blocked(n); u_buffer(1:end-1)];
    
    % Tín hiệu lỗi
    e = y_fixed(n) - w' * u_buffer;
    
    % Cập nhật trọng số NLMS
    u_norm = u_buffer' * u_buffer + alpha;
    w = w + (mu / u_norm) * e * u_buffer;
    
    % Đầu ra
    y_enhanced(n) = e;
end

% Chuẩn hóa
y_enhanced = y_enhanced / (max(abs(y_enhanced)) + eps);
y_blocked = y_blocked / (max(abs(y_blocked)) + eps);

end

