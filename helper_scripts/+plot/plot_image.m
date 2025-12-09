function h = plot_image(img, colored_by, units, save_path)
    % No export unless save path is specified
    if nargin < 4
        save_path = '';
    end
    
    colored_by = string(colored_by);
    plot_title = "Vessel skeleton color-mapped to " + colored_by;
    cmap = [0, 0, 0; turbo(255)];
    h = imagesc(img);
    axis image;
    axis off;
    colormap(cmap);
    cb = colorbar;
    cb.Label.String = colored_by + " (" + string(units) + ")";
    cb.Label.FontSize = 12;
    title(plot_title);

    if ~isempty(save_path)
        set(gcf, 'Visible', 'off');
        exportgraphics(gcf, save_path, 'Resolution', 600);
        close(gcf);
        h = [];
    else
        set(gcf, 'Visible', 'on');
        drawnow;
    end
end