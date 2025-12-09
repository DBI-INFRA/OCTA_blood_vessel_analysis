%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 3: Quantification of Blood Vessel Morphology
%%%%
%%%% Version:     1.3
%%%% Date:        25/02/2024
%%%%
%%%% Authors:     Tricia Loo (DBI-Infra IACF, tricia.loo@sund.ku.dk)
%%%%              Julia Mertesdorf (DBI-Infra IACF, jume@di.ku.dk)
%%%%              Peidi Xu (DBI-Infra IACF, peidi.xu@sund.ku.dk)
%%%%              Jesko Wagner (DBI-Infra IACF, jesko.wagner@sund.ku.dk)
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
parameters = load(fullfile(result_dir, 'Parameters.mat'));

%%% Modify any existing configuration below by setting:
% ParameterName = newValue;
% e.g. thresholding.method = "test_all";
% e.g. thresholding.sensitivity = [0.1, 0.2, 0.3, 0.4, 0.5];
% e.g. quantify_frame = 156;
% e.g. quantify_range_um = 5;


%% ================ Parse Inputs and Allocate Storage =====================
% Update modified variables
existingVars = fieldnames(parameters);
for i = 1:numel(existingVars)
    eval([existingVars{i} ' = parameters.' existingVars{i} ';']);
end
save(fullfile(result_dir, 'Parameters.mat'), existingVars{:});

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
    morphTable = cell(num_images, 4);

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
    [morph, readErrorIdx] = plot.quant_morph(ff, parameters, imgInfo, readErrorIdx);
    morphTable(ff, :) = plot.summarise_morph_measurements(parameters, ff, morph.avg, log_to_console=true);

    % Save detailed results
    if save_detailed_results
        [~, filename, ~] = fileparts(filelist(ff).name);
        save_path = fullfile(detailed_result_path, filename);
        writetable(morph.object, [save_path '_perBranchMeasurements.csv'])

        % Write Integer Output
        imwrite(~morph.image.binary, [save_path '_vessels.tif'], compression="none")
        imwrite(~morph.image.skeleton, [save_path '_skeleton.tif'], compression="none")
        imwrite(morph.image.labeled_skeleton, [save_path '_labeledSkeleton.tif'], compression="none")

        % Plot mapped measurements
        plot.plot_image(morph.image.skeleton_diameters, ...
            "vessel diameter", "µm", ...
            [save_path '_vesselDiameter.png']);

        skeleton_branchdiameter = labelmapper(morph.image.labeled_skeleton, [morph.object.Label], [morph.object.MeanDiameter_um]);
        plot.plot_image(skeleton_branchdiameter, ...
            "mean branch diameter", "µm", ...
            [save_path '_branchMeanDiameter.png']);

        skeleton_branchlength = labelmapper(morph.image.labeled_skeleton, [morph.object.Label], [morph.object.Length_um]);
        plot.plot_image(skeleton_branchlength, ...
            "branch length", "µm", ...
            [save_path '_branchLength.png']);

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
morphTable = cell2table(morphTable);
morphTable.Properties.RowNames = {filelist.name};
morphTable.Properties.VariableNames ...
    = {'Mean_Diameter (um)', 'Mean_Branch_Length (um)', ...
       'Vessel_Density (vessel/mm2)', 'Fractal_Dimension'};
writetable(morphTable, fullfile(result_dir, 'MorphologyResults.csv'), WriteRowNames=true);
fprintf("\n"); disp(morphTable);
fprintf("\nOCTA Script 3: Quantification of Blood Vessel Morphology is DONE\n");