function [S, A_est, W] = fastica_robust(X)
    % FastICA Implementation (Fixed-point algorithm)
    % X: Mixed signals (Rows = sensors, Cols = samples)
    
    [M, N] = size(X);
    
    % --- SUA LOI TAI DAY (Dung bsxfun de truong thich moi phien ban MATLAB) ---
    % 1. Centering: Tru di gia tri trung binh
    mu = mean(X, 2);
    X = bsxfun(@minus, X, mu);
    % ------------------------------------------------------------------------
    
    % 2. Whitening
    C = (X * X') / (N - 1); % Tinh Covariance thu cong de dam bao chinh xac
    [E, D] = eig(C);
    d = diag(D);
    
    % Sap xep eigenvalues de dam bao on dinh
    [d, idx] = sort(d, 'descend');
    E = E(:, idx);
    
    % Tranh chia cho 0 hoac so phuc
    inv_sqrt_D = diag(1 ./ (sqrt(d) + eps)); 
    V = inv_sqrt_D * E';
    Z = V * X;
    
    % 3. FastICA Algorithm
    max_iter = 1000;
    tol = 1e-4;
    W = randn(M, M);
    W = orth(W')'; % Orthogonalize initial weights
    
    for i = 1:M
        w = W(i,:)';
        for k = 1:max_iter
            w_old = w;
            
            % Update rule: w = E[z*g(w'z)] - E[g'(w'z)]*w
            wx = w' * Z;
            g_wx = tanh(wx);
            dg_wx = 1 - g_wx.^2;
            
            % Tinh trung binh (Mean)
            w = (Z * g_wx')/N - mean(dg_wx) * w;
            
            % Decorrelation (Gram-Schmidt)
            if i > 1
                w = w - W(1:i-1,:)' * W(1:i-1,:) * w;
            end
            w = w / norm(w);
            
            % Check convergence
            if abs(abs(w' * w_old) - 1) < tol
                break;
            end
        end
        W(i,:) = w';
    end
    
    % 4. Source Estimation
    S = W * Z;
    A_est = inv(W * V);
end