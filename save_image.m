function save_image(data, filename, title_str)
    if ~isnumeric(data)
        data = double(data);
    end
    h = figure;
    imagesc(data); title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
    colorbar; axis image; colormap jet;
    try
        saveas(h, filename);
    catch
        imwrite(mat2gray(data), filename);
    end
    close(h);
end
