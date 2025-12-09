function h = plot_panel(img, plot_title)
    h = imagesc(img);
    axis image;
    axis off;
    title(plot_title,'Interpreter','none');
end