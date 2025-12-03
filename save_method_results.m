function save_method_results(method_name, signals, result_folder, fs, t, X_filt)
% SAVE_METHOD_RESULTS - Luu ket qua va figure cho mot phuong phap
% 
% Inputs:
%   method_name - Ten phuong phap
%   signals - Ma tran ket qua (num_sources x signal_length)
%   result_folder - Thu muc ket qua chinh
%   fs - Tan so lay mau
%   t - Vector thoi gian
%   X_filt - Tin hieu dau vao da loc

    method_folder = fullfile(result_folder, method_name);
    if ~exist(method_folder, 'dir'), mkdir(method_folder); end
    
    num_sources = size(signals, 1);
    
    % Luu file audio
    for i = 1:num_sources
        sig = signals(i, :);
        output_file = fullfile(method_folder, sprintf('source_%d.wav', i));
        audiowrite(output_file, sig', fs);
    end
    
    % Tao figure
    fig = figure('Name', ['Ket qua: ' method_name], 'Position', [100, 100, 1400, 900], 'Visible', 'off');
    
    % Input signals
    subplot(3, 2, 1);
    plot(t, X_filt(1, :));
    xlabel('Thoi gian (s)'); ylabel('Bien do');
    title('Dau vao: Microphone 1');
    grid on; xlim([0, max(t)]);
    
    subplot(3, 2, 2);
    plot(t, X_filt(2, :));
    xlabel('Thoi gian (s)'); ylabel('Bien do');
    title('Dau vao: Microphone 2');
    grid on; xlim([0, max(t)]);
    
    % Output signals
    for i = 1:min(2, num_sources)
        subplot(3, 2, 2*i + 1);
        plot(t, signals(i, :));
        xlabel('Thoi gian (s)'); ylabel('Bien do');
        title(['Ket qua: Nguon ' num2str(i)]);
        grid on; xlim([0, max(t)]);
        
        % Spectrogram
        subplot(3, 2, 2*i + 2);
        NFFT = 2048;
        NOVERLAP = floor(NFFT * 0.75);
        WINDOW = hamming(NFFT);
        [S_spec, F_spec, T_spec] = spectrogram(signals(i, :), WINDOW, NOVERLAP, NFFT, fs);
        imagesc(T_spec, F_spec, 20*log10(abs(S_spec) + eps));
        axis xy; colorbar;
        xlabel('Thoi gian (s)'); ylabel('Tan so (Hz)');
        title(['Spectrogram - Nguon ' num2str(i)]);
        ylim([0, min(4000, fs/2)]);
    end
    
    % Luu figure
    fig_file = fullfile(method_folder, 'results_figure.fig');
    savefig(fig, fig_file);
    close(fig);
end

