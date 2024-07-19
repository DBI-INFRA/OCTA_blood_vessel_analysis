%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 3: Quantification of Blood Vessel Morphology
%%%%
%%%% Version:     1.0
%%%% Date:        16/07/2024
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
result_path = 'results/240719_141653';

% Parameters for skeletonization (median-filter & Frangi-filter)
median_filter_size = 5;
frangi_opts.sigmarange = [1 6];
frangi_opts.sigmastepsize = 2;
frangi_opts.correctionconst1 = 0.8;
frangi_opts.correctionconst2 = 15;
thresholding_method = "fuzzy_thresholding";  % OPTIONS: "fuzzy_thresholding", "local_adaptive_thresholding"
save_skeleton = true;

% Parameter for the depth range around the estimated SPD (for mean intensity projection)
SPD_range = 30;

% ================= Parse input and preallocate storage ===================
fileList  = dir(fullfile(result_path, "1_AlignedImages", "*.tif*"));
imgInfo  = readtable(fullfile(result_path, 'ImageSummary.csv'), 'ReadRowNames',true);

num_images = length(fileList);
meanDiameter  = zeros(num_images, 1);
meanLength    = zeros(num_images, 1);
meanDensity   = zeros(num_images, 1);
fracDimension = zeros(num_images, 1);

skeleton_result_path = fullfile(result_path, '3_MIP_skeletonization_results');

if save_skeleton && ~isfolder(skeleton_result_path)
    mkdir(skeleton_result_path);
end

% ================= Quantify all images in folder =========================
fprintf("\nProcessing %d images:\n", num_images);
for ff = 1:length(fileList)
    
    % 1) Parse image information
    img_path   = fullfile(fileList(ff).folder, fileList(ff).name);

    % Resolution stored as pixels per centimeters in Tiff metadata 
    % Convert into um/pixel (1cm = 10e4 um)
    img_reso   = (10^4)./[imgInfo{fileList(ff).name,'XRes'}, imgInfo{fileList(ff).name,'YRes'}];
    SPD_depth  = imgInfo{fileList(ff).name, 'SPD_Depth'};
    SPD_slices = [SPD_depth-SPD_range, SPD_depth+SPD_range];

    % Read image around SPD depth
    ImageStack = tiffreadVolume(img_path, 'PixelRegion', {[1 inf], [1 inf], SPD_slices});
    
    % Create mean-intensity-projection and skeletonize image
    StackMIP = mean(ImageStack, 3);
    if max(StackMIP(:)) <= 1
        StackMIP = StackMIP * 255;
    end
    StackMIP = double(int16(squeeze(StackMIP)));
    save_path = fullfile(skeleton_result_path, fileList(ff).name);
    [skeleton, binary_img] = Skeletonization(StackMIP, median_filter_size, frangi_opts, thresholding_method, save_skeleton, save_path);

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
morphTable = table(meanDiameter, meanLength, meanDensity, fracDimension);
morphTable.Properties.RowNames = {fileList.name};
morphTable.Properties.VariableNames ...
    = {'Mean_Diameter (um)', 'Mean_Branch_Length (um)', ...
       'Vessel_Density (vessel/mm2)', 'Fractal_Dimension'};
writetable(morphTable, fullfile(result_path, 'MorphologyResults.csv'), WriteRowNames=true);
fprintf("\n"); disp(morphTable);
fprintf("\nOCTA Script 3: Quantification of Blood Vessel Morphology is DONE\n");
