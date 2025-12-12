function summaryCell = summarise_morph_measurements(parameters, ff, measurements, varargin)
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'log_to_console', true, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});
    log_to_console = p.Results.log_to_console;

    meanDiameter = measurements.meanDiameter;
    meanLength = measurements.meanLength;
    meanDensity = measurements.meanDensity;
    fracDimension = measurements.fracDimension;

    summaryCell = {meanDiameter, meanLength, meanDensity, fracDimension};

    if log_to_console
        % Print the quantification results
        fprintf("    %d. %s:  Diameter: %.3f,  Length: %.3f,  Density: %.3f,  Fractal dimension: %.3f\n", ...
            ff, parameters.filelist(ff).name, meanDiameter, meanLength, meanDensity, fracDimension);
    end
end
