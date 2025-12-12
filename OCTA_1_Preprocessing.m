%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 1: Preprocessing & Skin-Border Alignment
%%%%
%%%% Version:     1.3
%%%% Date:        25/02/2025
%%%%
%%%% Authors:     Tricia Loo (DBI-Infra IACF, tricia.loo@sund.ku.dk)
%%%%              Julia Mertesdorf (DBI-Infra IACF, jume@di.ku.dk)
%%%%              Peidi Xu (DBI-Infra IACF, peidi.xu@sund.ku.dk)
%%%%
%%%% Description: The first part of the OCTA analysis pipeline performs
%%%%              image pre-processing by automatically cropping, filtering
%%%%              and aligning images to the skin-border.
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ================= Initialisation =======================================
close all
clear all
addpath('helper_scripts');

%%% Specify result_dir (path where the desired 'Parameters.mat' is stored)
result_dir = 'results\yymmdd_hhmmss\';

%%% Load exisiting user configurations
load(fullfile(result_dir, 'Parameters.mat'));

%%% Modify any existing configuration below by setting 
% ParameterName = newValue;
% e.g. WT_orient = 'horizontal';

%% ================ Parse Inputs and Allocate Storage =====================
% Update modified variables
existingVars = whos("-file",fullfile(result_dir, 'Parameters.mat'));
save(fullfile(result_dir, 'Parameters.mat'), existingVars.name)

try
    % Load user configurations and image information
    load(fullfile(result_dir, 'Parameters.mat'))
    imgInfo = readtable(fullfile(result_dir, 'ImageSummary.csv'), 'Delimiter', ',');
    
    % Default (optimised) parameters for wavelet filtering (not recommended to change)
    % See waveletfft3D.m for details
    gamma = 10; order = 4; wname = 'db20';             
    
    % Preallocate storage
    fileslice = ones(length(filelist), 2);
    processed_image_dir = fullfile(result_dir, 'ProcessedImages');

    % Print user parameters & image path
    fprintf("\nIMAGE FOLDER PATH: \n  <%s>\n", image_dir);
    fprintf("\nOUTPUT FOLDER PATH: \n  <%s>\n", processed_image_dir);
    fprintf(['\nSELECTED PARAMETERS:\n' ...
             '  Cropping: %d, Sensitivity = %.2f \n' ...
             '  Wavelet Filter Orientation, Dimensions, Reps: %s, %s, %d \n' ...
             '  Alignment: %d \n' ...
             '  Write Z Displacement: %d \n'], ...
             prep_opts.AutoCrop, prep_opts.CropSensitivity, ...
             prep_opts.WT_orient, prep_opts.WT_dim, prep_opts.WT_reps, ...
             prep_opts.AutoAlign, prep_opts.WriteZDisplacement);

catch ME
    fprintf("Check configuration file, preprocess_image must be set to true: ")
    rethrow(ME);
end

% Create folder for processed image
mkdir(processed_image_dir);

%% ================= Preprocess all images in input folder ================
fprintf('\nPROCESSING FILES: \n');
for ff = 1:length(filelist)
    img_path = fullfile(image_dir, filelist(ff).name);
    fprintf("  Image: <%s>\n", filelist(ff).name);

    % 1) Preprocessing steps
    % Read and crop stacks to relevant (non-noise) data
    OCTStack = tiffreadVolume(img_path);
    OCTStack = OCTStack(:,:,:,1);
    if prep_opts.AutoCrop
        [CroppedOCTStack, crop_range] = autocropOCTstack(OCTStack, prep_opts.CropSensitivity);
        fileslice(ff,:) = crop_range;
    else 
        CroppedOCTStack = OCTStack;
        fileslice(ff,:) = [1, size(OCTStack,3)];
    end

    % Perform wavelet transform

    FilteredStack = waveletfft3D(CroppedOCTStack, gamma, order, wname,...
                        prep_opts.WT_orient, prep_opts.WT_reps, prep_opts.WT_dim);

    % Z alignment of image
    if prep_opts.AutoAlign
        [ResultStack, z_shift] = zalignstack(FilteredStack, 'median');
    else
        ResultStack = double(FilteredStack);
    end

    % 2) Write Results to file
    % Write aligned image
    imwrite(ResultStack(:,:,1), fullfile(processed_image_dir, filelist(ff).name));
    for ii = 2:size(ResultStack, 3)
         imwrite(ResultStack(:,:,ii), fullfile(processed_image_dir, filelist(ff).name), "Writemode", "append");
    end
    
    % Write zdisplacement
    if prep_opts.WriteZDisplacement & prep_opts.AutoAlign
        if ~exist(fullfile(processed_image_dir, 'Displacement'), "file")
            mkdir(fullfile(processed_image_dir, 'Displacement'));
        end
        imwrite(double(z_shift)/255, fullfile(processed_image_dir, 'Displacement', filelist(ff).name));
    end

end

%% ================= Write image processing summary =======================
imgInfo.("First Slice") = fileslice(:, 1);
imgInfo.("Last Slice") = fileslice(:, 2);
writetable(imgInfo, fullfile(result_dir, 'ImageSummary.csv'));

fprintf("\nOCTA Script 1: Preprocessing & Skin-Border Alignment DONE\n\n");