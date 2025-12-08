function log_morph_measurements(parameters, ff, measurements)
    meanDiameter = measurements.meanDiameter;
    meanLength = measurements.meanLength;
    meanDensity = measurements.meanDensity;
    fracDimension = measurements.fracDimension;

    % Print the quantification results
    fprintf("    %d. %s:  Diameter: %.3f,  Length: %.3f,  Density: %.3f,  Fractal dimension: %.3f\n", ...
        ff, parameters.filelist(ff).name, meanDiameter, meanLength, meanDensity, fracDimension);
end
