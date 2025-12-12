%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 3: Quantification of Blood Vessel Morphology
%%%%
%%%% Version:     1.3
%%%% Date:        25/02/2024
%%%%
%%%% Authors:     Tricia Loo (DBI-Infra IACF, tricia.loo@sund.ku.dk)
%%%%              Julia Mertesdorf (DBI-Infra IACF, jume@di.ku.dk)
%%%%              Peidi Xu (DBI-Infra IACF, peidi.xu@sund.ku.dk)
%%%%
%%%% Description: This third section of the pipeline performs quantitative
%%%%              analysis on a series of cropped OCT images stored in the
%%%%              provided input directory. It calculates four
%%%%              morphological metrics for each image: mean vessel
%%%%              diameter, vessel length, vessel density and fractal
%%%%              dimension, from a mean-intensity-projected image of
%%%%              'quantify_range_um's around a given 'quantify_frame'.
%%%%              All results are saved in a CSV file.
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ================= Initialization =======================================
close all
clear all
addpath('helper_scripts');

%%% Specify result_dir (path where the desired 'Parameters.mat' is stored)
result_dir = 'results\yymmdd_hhmmss\';

%%% Load exisiting user configurations
load(fullfile(result_dir, 'Parameters.mat'));

%%% Modify any existing configuration below by setting:
% ParameterName = newValue;
% e.g. thresholding.method = "test_all";
% e.g. thresholding.sensitivity = [0.1, 0.2, 0.3, 0.4, 0.5];
% e.g. quantify_frame = 156;
% e.g. quantify_range_um = 5;


%% ================ Parse Inputs and Allocate Storage =====================
% Update modified variables
existingVars = whos("-file",fullfile(result_dir, 'Parameters.mat'));
save(fullfile(result_dir, 'Parameters.mat'), existingVars.name)

try
    % Load image information and find any pre-processed image
    imgInfo = readtable(fullfile(result_dir, 'ImageSummary.csv'), 'ReadRowNames', true, 'Delimiter', ',');
    if exist(fullfile(result_dir, 'ProcessedImages'), 'dir')
        image_dir = fullfile(result_dir, 'ProcessedImages');
    end

    % Calculate image sizes
    quantify_range  = round(quantify_range_um/pixel_size(3));
    img_reso        = pixel_size(1:2);

    % Preallocate storage
    num_images    = length(filelist);
    meanDiameter  = zeros(num_images, 1);
    meanLength    = zeros(num_images, 1);
    meanDensity   = zeros(num_images, 1);
    fracDimension = zeros(num_images, 1);
    readErrorIdx  = cell(num_images, 2);
    [readErrorIdx{:, 1}] = deal(false);

    detailed_result_path = fullfile(result_dir, 'DetailedResults');

catch ME
    fprintf("Check configuration file: ")
    rethrow(ME);
end

% Create folder for saving detailed results
if save_detailed_results && ~isfolder(detailed_result_path)
    mkdir(detailed_result_path);
end

% ================= Quantify all images in folder =========================
fprintf("\nProcessing %d images:\n", num_images);
for ff = 1:length(filelist)

    % 1) Parse image information
    img_path   = fullfile(image_dir, filelist(ff).name);

    if isscalar(quantify_frame)
        central_frame = quantify_frame;
    else
        central_frame  = imgInfo{filelist(ff).name, quantify_frame};
        if isnan(central_frame)
            meanDiameter(ff)    = NaN;
            meanLength(ff)      = NaN;
            meanDensity(ff)     = NaN;
            fracDimension(ff)   = NaN;

            readErrorIdx{ff, 1} = true;
            readErrorIdx{ff, 2} = NaN;
            continue
        end
    end

    quantify_slices = [central_frame-quantify_range, central_frame+quantify_range];

    % Read image around chosen depth
    try
        ImageStack = tiffreadVolume(img_path, 'PixelRegion', {[1 inf], [1 inf], quantify_slices});
    catch
        ImageStack = tiffreadVolume(img_path, 'PixelRegion', {[1 inf], [1 inf]});

        quantify_slices(2) = size(ImageStack, 3);
        readErrorIdx{ff, 1} = true;
        readErrorIdx{ff, 2} = [quantify_slices, diff(quantify_slices)*pixel_size(3)];
    end

    if numel(size(ImageStack)) == 4
        ImageStack = ImageStack(:,:,:,1);
    end

    % Create mean-intensity-projection and skeletonize image
    StackMIP = mean(ImageStack, 3);
    if max(StackMIP(:)) <= 1
        StackMIP = StackMIP * 255;
    end
    StackMIP = double(int16(squeeze(StackMIP)));
    [skeleton, binary_img] = skeletonization(StackMIP, median_filter_size, frangi_opts, thresholding, false, "");

    % 2) Quantify the blood vessel network morphology
    % Split skeleton into branches
    [labeled_skeleton, ~] = splitbranches(skeleton);

    % Calculate vessel diameter at each point along skeleton
    skeleton_diameters = double(bwdist(~binary_img).*skeleton*2);

    % Make Vessel Measurements per Branch
    vessel_morph = regionprops(labeled_skeleton, skeleton_diameters, 'MeanIntensity', 'Area');
    vessel_label = regionprops(labeled_skeleton, labeled_skeleton, 'MeanIntensity');
    vessel_measurements = table([vessel_label.MeanIntensity]', ...
                                [vessel_morph.Area]', [vessel_morph.Area]'*img_reso(1), ...
                                [vessel_morph.MeanIntensity]', [vessel_morph.MeanIntensity]' *img_reso(1), ...
                                'VariableNames', {'Label', 'Length_px', 'Length_um', 'MeanDiameter_px', 'MeanDiameter_um'});
    vessel_measurements = vessel_measurements([vessel_measurements.Label]>1,:);
    num_vessels = size(vessel_measurements, 1);

    % Average Vessel Branch Length in um
    meanLength(ff) = mean([vessel_measurements.Length_um]);

    % Average Vessel Diameter in um (including at branch points)
    avg_diameter = mean(skeleton_diameters(skeleton));
    meanDiameter(ff) = avg_diameter*img_reso(1);

    % Average Vessel Density per mm2
    img_reso_mm = img_reso*(10^-3);
    img_area = prod(size(binary_img).*img_reso_mm);
    meanDensity(ff) = num_vessels/img_area;

    % Fractal Dimension
    [n, r] = boxcount2D(skeleton);
    fracDimension(ff) = -1*fit(log(r)', log(n)', 'poly1').p1;

    % Print the quantification results
    fprintf("    %d. %s:  Diameter: %.3f,  Length: %.3f,  Density: %.3f,  Fractal dimension: %.3f\n", ...
        ff, filelist(ff).name, meanDiameter(ff), meanLength(ff), meanDensity(ff), fracDimension(ff));

    % Save detailed results
    if save_detailed_results
        [~, filename, ~] = fileparts(filelist(ff).name);
        save_path = fullfile(detailed_result_path, filename);
        writetable(vessel_measurements, [save_path '_perBranchMeasurements.csv'])

        % Write Integer Output
        imwrite(uint16(binary_img), [save_path '_vessels.tif'])
        imwrite(uint16(skeleton), [save_path '_skeleton.tif'])
        imwrite(uint16(labeled_skeleton), [save_path '_labeledSkeleton.tif'])

        % Plot mapped measurements
        turbo_on_black = [0, 0, 0; turbo(255)];

        imagesc(skeleton_diameters), axis image, axis off,
        colormap(turbo_on_black), colorbar
        title("Vessel skeleton color-mapped to vessel diameter")
        exportgraphics(gcf, [save_path '_vesselDiameter.png'], 'Resolution', 600);

        skeleton_branchdiameter = labelmapper(labeled_skeleton, [vessel_measurements.Label], [vessel_measurements.MeanDiameter_um]);
        imagesc(skeleton_branchdiameter), axis image, axis off,
        colormap(turbo_on_black), colorbar
        title("Vessel skeleton color-mapped to mean branch diameter")
        exportgraphics(gcf, [save_path '_branchMeanDiameter.png'], 'Resolution', 600);

        skeleton_branchlength = labelmapper(labeled_skeleton, [vessel_measurements.Label], [vessel_measurements.Length_um]);
        imagesc(skeleton_branchlength), axis image, axis off,
        colormap(turbo_on_black), colorbar
        title("Vessel skeleton color-mapped to branch length")
        exportgraphics(gcf, [save_path '_branchLength.png'], 'Resolution', 600);

        close
    end

end

% ================= Write results to file =================================
% Printout warnings for files where error occured
warningIdx  = [readErrorIdx{:, 1}];
if any(warningIdx)
    warningHead = {"The following files had errors in their SPD detection", "Frames read as SPD (NaN or start_frame, stop_frame, depth(um))"};
    warningBody = [{filelist(warningIdx).name}', readErrorIdx(warningIdx, 2)];
    warningLog  = [warningHead; warningBody];
    writecell(warningLog, fullfile(result_dir, 'WarningLog.txt'))
end

% Write final parameters to table
if exist("prep_opts")
    prep_params = [fieldnames(prep_opts), struct2cell(prep_opts)];
else
    prep_params = [];
end

param_table = [{"image_dir"}, {image_dir};
               {"output_dir"}, {output_dir};
               {"pixel_size"}, {pixel_size};
               prep_params;
               {"median_filter_size"}, {median_filter_size};
               fieldnames(frangi_opts), struct2cell(frangi_opts);
               fieldnames(thresholding), struct2cell(thresholding);
               {"quantify_frame"}, {quantify_frame};
               {"quantify_range_um"}, {quantify_range_um};
               {"save_detailed_results"}, {save_detailed_results}];
writecell(param_table, fullfile(result_dir, "Parameters.csv"))

% Morphology calculations
morphTable = table(meanDiameter, meanLength, meanDensity, fracDimension);
morphTable.Properties.RowNames = {filelist.name};
morphTable.Properties.VariableNames ...
    = {'Mean_Diameter (um)', 'Mean_Branch_Length (um)', ...
       'Vessel_Density (vessel/mm2)', 'Fractal_Dimension'};
writetable(morphTable, fullfile(result_dir, 'MorphologyResults.csv'), WriteRowNames=true);
fprintf("\n"); disp(morphTable);
fprintf("\nOCTA Script 3: Quantification of Blood Vessel Morphology is DONE\n");