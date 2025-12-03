function y = differential_microphone_array(X, order)
% DIFFERENTIAL_MICROPHONE_ARRAY - Differential Microphone Array
% 
% Inputs:
%   X - Tín hiệu đầu vào (M x N): M = số microphone, N = số mẫu
%   order - Bậc của differential array (1 = first-order, 2 = second-order)
%
% Output:
%   y - Tín hiệu đã được xử lý

[M, N] = size(X);

if M < 2
    error('Can it nhat 2 microphone cho differential array');
end

if order == 1
    % First-order differential: y = x1 - x2
    y = X(1, :) - X(2, :);
elseif order == 2
    % Second-order differential: y = x1 - 2*x2 + x3 (neu co 3 mic)
    if M >= 3
        y = X(1, :) - 2*X(2, :) + X(3, :);
    else
        % Với 2 mic, dùng first-order
        y = X(1, :) - X(2, :);
    end
else
    % Default: first-order
    y = X(1, :) - X(2, :);
end

% Chuẩn hóa
y = y / (max(abs(y)) + eps);

end

