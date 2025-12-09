function hAxes = prepare_canvas(N)
% prepareCanvas  Create a multi‐panel canvas for N plots (max 3 columns)
%
%   hAxes = prepareCanvas(N) creates a new figure with a tiled layout
%   that has up to 3 columns and as many rows as needed (ceil(N/3)).
%   It returns a 1×N array of axes handles hAxes, one for each plot slot.
%
%   Parameters:
%     N: positive integer, number of plots
% 
%   Outputs:
%     hAxes: 1×N array of axes handles

    narginchk(1,1);
    validateattributes(N, {'numeric'}, ...
        {'scalar','integer','positive'}, mfilename, 'N', 1);

    % Maximum number of columns
    nCols = min(3, N);

    % Compute required number of rows
    nRows = ceil(N / nCols);

    if nRows > 10
        error("Cannot generate figure with more than 10 rows")
    end

    % Create figure and tiled layout
    figure;
    t = tiledlayout(nRows, nCols, ...
        'Padding',    'compact', ...
        'TileSpacing','compact');

    % Preallocate axes handle array
    hAxes = gobjects(1, N);

    % Create each tile
    for k = 1:N
        hAxes(k) = nexttile(t);
    end

    % Turn off any extra tiles
    totalTiles = nRows * nCols;
    if totalTiles > N
        for k = (N+1):totalTiles
            ax = nexttile(t, k);
            axis(ax, 'off');
        end
    end

    linkprop(hAxes, 'Clim');
end
