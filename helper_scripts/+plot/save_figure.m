function save_figure(figH, parameters)
    % save_figure  Prompt user for a filename and save the figure.
    %
    %   save_figure(figH) pops up a “Save As” dialog for common image
    %   formats and only returns once the user picks a file (or cancels).
    %   If the file already exists, it asks to confirm overwrite.
    %
    %   parameters:
    %     figH: handle to the figure you want to save (default gcf)
    %     parameters: to launch GUI in results dir
    %
    if nargin<1 || isempty(figH) || ~ishghandle(figH,'figure')
        figH = gcf;
    end

    % 1) Prompt for filename & path
    [fname, pth, ~] = uiputfile( ...
        {'*.svg', 'SVG vector file (*.svg)'; ...
        '*.tif;*.tiff', 'TIF image (*.tif,*.tiff)'; ...
         '*.png','PNG image (*.png)'; ...
         '*.jpg;*.jpeg','JPEG image (*.jpg,*.jpeg)'; ...
         '*.pdf','PDF vector file (*.pdf)'}, ...
        'Save Figure As', parameters.output_dir);

    % 2) Handle cancel
    if isequal(fname,0) || isequal(pth,0)
        disp('Save cancelled by user.');
        return;
    end

    fullfileName = fullfile(pth, fname);

    % 3) Save using your preferred method
    [~,~,ext] = fileparts(fname);
    try
        switch lower(ext)
            case '.svg'
                exportgraphics(figH, fullfileName, 'ContentType','vector', 'Resolution',600);
            case {'.tif','.tiff'}
                exportgraphics(figH, fullfileName, 'Resolution',600);
            case '.png'
                exportgraphics(figH, fullfileName, 'Resolution',600);
            case {'.jpg','.jpeg'}
                exportgraphics(figH, fullfileName, 'Resolution',600);
            case '.pdf'
                exportgraphics(figH, fullfileName,'ContentType','vector','Resolution',600);
            otherwise
                warning('Unknown extension "%s", using saveas.', ext);
                saveas(figH, fullfileName);
        end
        fprintf('Saved figure to:\n  %s\n', fullfileName);
    catch ME
        warning('Failed to save figure');
    end
end
