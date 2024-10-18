%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 3: Quantification of Blood Vessel Morphology
%%%%
%%%% Version:     1.2
%%%% Date:        26/07/2024
%%%%
%%%% Authors:     Tricia Loo (DBI-Infra IACF, tricia.loo@sund.ku.dk)
%%%%              Julia Mertesdorf (DBI-Infra IACF, jume@di.ku.dk)
%%%%
%%%% Description: This third section of the pipeline performs quantitative  
%%%%              analysis on a series of cropped OCT images stored in the  
%%%%              provided input directory. At the depth of the superficial 
%%%%              plexus (SPD) of the blood vessel network, it calculates 
%%%%              various morphological metrics for each image, including 
%%%%              the mean vessel diameter, vessel length, vessel density  
%%%%              and fractal dimension. All results are saved in a CSV file.
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ================= Initialization ========================================
close all
clear all
addpath('helper_scripts');

% ================= User Parameters =======================================
% Input image directory containing the cropped OCT images

result_path = 'results/240730_170157';

% Manually set to the SPD frame number if you already have the number, otherwise
% leave it to -1 so that it will be inferred from the previous step. Note
% that it is in the number of frame, not number of micrometers

SPD_frame = -1;

% Set to true if you want to test & view different thresholding methods
test_thresholding_methods = false;
% Choose different values for the sensitivity parameter in local adaptive
% thresholding. The <test_threshold>-method will then apply local adaptive
% thresholding with each different sensitivity parameter, and create a new
% subplot for each version that can be found under the following path:
% <3_MIP_skeletonization_results/Thresholding_Comparision>
% Note: If test_thresholding is enabled, this script will apply local
% adaptive thresholding with sensitivity=0.3 as the default
local_adaptive_thresholding_sensitivities = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];

% Parameters for skeletonization (median-filter & Frangi-filter)
save_skeleton = true;

% Parameter for the depth range around the estimated SPD in microns
SPD_range_um = 30;

% The method to apply for thresholding, i.e. dividing the image into foreground (vessels) & background
% The selectable options are: "local_adaptive_thresholding", "fuzzy_thresholding", "otsu_thresholding"
thresholding_method = "local_adaptive_thresholding";

% The sensitivity of the local adaptive thresholding algorithm (higher values will
% pick up more & smaller vessels, but potentially also more noise).
thresholding_sensitivity = 0.3;


% ================= Parse input and preallocate storage ===================
fileList    = dir(fullfile(result_path, "1_AlignedImages", "*.tif*"));
params_path = fullfile(result_path, 'InputParameters.mat');
load(params_path);

% ================= The below parameters are set different from OCTA_2 and OCTA_3 ===================
thresholding.sensitivity = thresholding_sensitivity;
thresholding.method = thresholding_method;

% If you want to disable morphological opening, uncomment the following line
% thresholding.opening_size = 0; 

if test_thresholding_methods
    thresholding.method = "test_all";
    thresholding.test_sensitivities = local_adaptive_thresholding_sensitivities;
end

SPD_range   = round(SPD_range_um/pixel_size(3));
img_reso    = pixel_size(1:2);
imgInfo     = readtable(fullfile(result_path, 'ImageSummary.csv'), 'ReadRowNames',true);

num_images    = length(fileList);
meanDiameter  = zeros(num_images, 1);
meanLength    = zeros(num_images, 1);
meanDensity   = zeros(num_images, 1);
fracDimension = zeros(num_images, 1);
readErrorIdx  = cell(num_images, 2);
[readErrorIdx{:, 1}] = deal(false);

skeleton_result_path = fullfile(result_path, '3_MIP_skeletonization_results');

if save_skeleton && ~isfolder(skeleton_result_path)
    mkdir(skeleton_result_path);
end

% Add new user parameters to parameters table and save to csv file
save(params_path, "SPD_range_um", "-append");
table_path = fullfile(result_path, 'InputParameters.csv');
params_value = [{pixel_size}; {SPD_range_um}; struct2cell(preprocess_opts); ...
                {median_filter_size}; struct2cell(frangi_opts); struct2cell(thresholding)];
params_names = [{"pixel_size"}; {"SPD_range_um"}; fieldnames(preprocess_opts);...
                {"median_filter_size"}; fieldnames(frangi_opts); fieldnames(thresholding)];
params_table = [params_names, params_value];
fprintf("\n"); disp(params_table);
writecell(params_table, table_path);

% ================= Quantify all images in folder =========================
fprintf("\nProcessing %d images:\n", num_images);
for ff = 1:length(fileList)
    
    % 1) Parse image information
    img_path   = fullfile(fileList(ff).folder, fileList(ff).name);

    if (SPD_frame == -1)
        SPD_frame  = imgInfo{fileList(ff).name, 'SPD_Frame'};
        if isnan(SPD_frame)
            meanDiameter(ff)    = NaN;
            meanLength(ff)      = NaN;
            meanDensity(ff)     = NaN;
            fracDimension(ff)   = NaN;
    
            readErrorIdx{ff, 1} = true;
            readErrorIdx{ff, 2} = NaN;
            continue
        end
    end

    SPD_slices = [SPD_frame-SPD_range, SPD_frame+SPD_range];

    % Read image around SPD depth
    try
        ImageStack = tiffreadVolume(img_path, 'PixelRegion', {[1 inf], [1 inf], SPD_slices});
    catch
        ImageStack = tiffreadVolume(img_path, 'PixelRegion', {[1 inf], [1 inf]});
        
        SPD_slices(2) = size(ImageStack, 3);
        readErrorIdx{ff, 1} = true; 
        readErrorIdx{ff, 2} = [SPD_slices, diff(SPD_slices)*pixel_size(3)];
    end
    
    % Create mean-intensity-projection and skeletonize image
    StackMIP = mean(ImageStack, 3);
    if max(StackMIP(:)) <= 1
        StackMIP = StackMIP * 255;
    end
    StackMIP = double(int16(squeeze(StackMIP)));
    save_path = fullfile(skeleton_result_path, fileList(ff).name);
    [skeleton, binary_img] = Skeletonization(StackMIP, median_filter_size, frangi_opts, thresholding, save_skeleton, save_path);
    
    %imwrite(binary_img, strcat(save_path, '_binary.png'))

    % 2) Quantify the blood vessel network morphology
    % Average Vessel Diameter in um
    avg_diameter = vessel_diameter(binary_img, skeleton);
    meanDiameter(ff) = avg_diameter*img_reso(1);

    % Average Vessel Branch Length in um
    [avg_length, num_vessels] = vessel_length(skeleton);
    meanLength(ff) = avg_length*img_reso(1);

    % Average Vessel Density per mm2
    img_reso_mm = img_reso*(10^-3);
    img_area = prod(size(binary_img).*img_reso_mm);
    meanDensity(ff) = num_vessels/img_area;

    % Fractal Dimension
    [n, r] = BoxCount2D(skeleton);
    fracDimension(ff) = -1*fit(log(r)', log(n)', 'poly1').p1;

    % Print the quantification results
    fprintf("    %d. %s:  Diameter: %.3f,  Length: %.3f,  Density: %.3f,  Fractal dimension: %.3f\n", ...
        ff, fileList(ff).name, meanDiameter(ff), meanLength(ff), meanDensity(ff), fracDimension(ff));
end

% ================= Write results to file =================================
% Printout warnings for files where error occured
warningHead = {"The following files had errors in their SPD detection", "Frames read as SPD (NaN or start_frame, stop_frame, depth(um))"};
warningIdx  = [readErrorIdx{:, 1}];
warningBody = [{fileList(warningIdx).name}', readErrorIdx(warningIdx, 2)];
warningLog  = [warningHead; warningBody];
writecell(warningLog, fullfile(result_path, 'WarningLog.txt'))

% Morphology calculations
morphTable = table(meanDiameter, meanLength, meanDensity, fracDimension);
morphTable.Properties.RowNames = {fileList.name};
morphTable.Properties.VariableNames ...
    = {'Mean_Diameter (um)', 'Mean_Branch_Length (um)', ...
       'Vessel_Density (vessel/mm2)', 'Fractal_Dimension'};
writetable(morphTable, fullfile(result_path, 'MorphologyResults.csv'), WriteRowNames=true);
fprintf("\n"); disp(morphTable);
fprintf("\nOCTA Script 3: Quantification of Blood Vessel Morphology is DONE\n");
