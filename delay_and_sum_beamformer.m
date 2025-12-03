function y = delay_and_sum_beamformer(X, fs, look_direction)
% DELAY_AND_SUM_BEAMFORMER - Delay and Sum Beamformer
% 
% Inputs:
%   X - Tín hiệu đầu vào (M x N): M = số microphone, N = số mẫu
%   fs - Tần số lấy mẫu (Hz)
%   look_direction - Hướng nhìn (góc độ, 0 = thẳng, 90 = bên phải)
%
% Output:
%   y - Tín hiệu đã được beamforming

[M, N] = size(X);

% Khoảng cách giữa các microphone (giả định 5cm cho 2 mic)
d = 0.05;  % mét
c = 343;   % Tốc độ âm thanh (m/s)

% Chuyển góc sang radian
theta_look = deg2rad(look_direction);

% Tính toán delay cho mỗi microphone
delays_samples = zeros(M, 1);
for m = 1:M
    mic_pos = (m - 1) * d;
    delay_time = mic_pos * sin(theta_look) / c;
    delays_samples(m) = round(delay_time * fs);
end

% Điều chỉnh delay (thêm padding)
max_delay = max(abs(delays_samples));
X_padded = [zeros(M, max_delay), X, zeros(M, max_delay)];

% Áp dụng delay và tổng hợp
y = zeros(1, N);
for n = 1:N
    idx = n + max_delay;
    sum_val = 0;
    for m = 1:M
        sample_idx = idx - delays_samples(m);
        if sample_idx > 0 && sample_idx <= size(X_padded, 2)
            sum_val = sum_val + X_padded(m, sample_idx);
        end
    end
    y(n) = sum_val / M;
end

% Chuẩn hóa
y = y / (max(abs(y)) + eps);

end

