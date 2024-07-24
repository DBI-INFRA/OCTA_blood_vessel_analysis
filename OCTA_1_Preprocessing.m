%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 1: Preprocessing & Skin-Border Alignment
%%%%
%%%% Version:     1.0
%%%% Date:        16/07/2024
%%%%
%%%% Authors:     Tricia Loo (DBI-Infra IACF, tricia.loo@sund.ku.dk)
%%%%              Julia Mertesdorf (DBI-Infra IACF, jume@di.ku.dk)
%%%%
%%%% Description: The first part of the OCTA analysis pipeline performs
%%%%              image preprocessing by automatically cropping, filtering
%%%%              and aligning images to the skin-border.
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ================= Initialisation =======================================
close all
clear all
addpath('helper_scripts');

%% ================= User Parameters ======================================
% Input image directory containing the original OCT images.
% Note: The image data should be in local drive for fast read/write speed
input_dir = 'C:\Users\Tricia\Desktop\NinaMK\ADSample1';

% Result directory & optional result
result_dir = 'results';

% Pixel size in micrometers
pixel_size = [6.49, 6.49, 2.45];

% Optional results to write
writeZDisplacement = true;

% Automated image cropping
% The sample data has already been cropped to size and should have AutoCrop
% set to false. If using uncropped data, please set AutoCrop to true.
AutoCrop = true;
% Value between [0, 1] to control amount of end-cropping for the OCT stack. 
% 0 means more data will be kept and 1 means more data will be cropped. 
% Can be set to 'None' to prevent cropping from end of dataset
CropSensitivity = 0.5;

% Wavelet Transform Parameters
gamma = 10; order = 4; wname = 'db20'; reps = 2; orient = 'both'; dim = 3;

%% ================= Parse input and preallocate storage ==================

filelist  = dir(fullfile(input_dir, '*.tif*'));

fileslice = ones(length(filelist), 2);
fileinfo  = cell(length(filelist), 3);

timenow = char(datetime("now"), "yyMMdd_HHmmss");
output_dir  = fullfile(result_dir, timenow);
imwrite_dir = fullfile(output_dir, '1_AlignedImages');
mkdir(imwrite_dir);

% Save all user parameters in a table
params_path = fullfile(output_dir, 'InputParameters.mat');
preprocess_opts = struct("AutoCrop", AutoCrop, "CropSensitivity", CropSensitivity, ...
                    "WT_gamma", gamma, "WT_order", order, "WT_wname", wname, ...
                    "WT_reps", reps, "WT_orient", orient, "WT_dim", dim);
save(params_path, "pixel_size", "preprocess_opts");

%% ================= Preprocess all images in input folder ================
for ff = 1:length(filelist)
    img_path = fullfile(input_dir, filelist(ff).name);
    img_info = imfinfo(img_path);
    fileinfo(ff,:) = {img_info(1).XResolution, img_info(1).YResolution, ... 
                      img_info(1).ResolutionUnit};

    % 1) Preprocessing steps
    % Read and crop stacks to relevant (non-noise) data
    OCTStack = tiffreadVolume(img_path);
    OCTStack = OCTStack(:,:,:,1);
    if AutoCrop
        [CroppedOCTStack, crop_range] = AutoCropOCTStack(OCTStack, CropSensitivity);
        fileslice(ff,:) = crop_range;
    else 
        CroppedOCTStack = OCTStack;
        fileslice(ff,:) = [1, size(OCTStack,3)];
    end

    % Perform wavelet transform
    FilteredStack = WaveletFFT3D(CroppedOCTStack, ...
                                 gamma, order, wname, orient, reps, dim);

    % Z alignment of image
    [AlignedStack, z_shift] = ZAlignStack(FilteredStack, 'median');

    % 2) Write Results to file
    % Write aligned image
    imwrite(AlignedStack(:,:,1), fullfile(imwrite_dir, filelist(ff).name));
    for ii = 2:size(AlignedStack, 3)
         imwrite(AlignedStack(:,:,ii), fullfile(imwrite_dir, filelist(ff).name), "Writemode", "append");
    end
    
    % Write zdisplacement
    if writeZDisplacement
        if ~exist(fullfile(imwrite_dir, 'Displacement'), "file")
            mkdir(fullfile(imwrite_dir, 'Displacement'));
        end
        imwrite(double(z_shift)/255, fullfile(imwrite_dir, 'Displacement', filelist(ff).name));
    end

end

%% ================= Write summary of image alignment & resolution ========
T = array2table(fileslice, "RowNames", {filelist.name}, "VariableNames", {'First Slice', 'Last Slice'});
writetable(T, fullfile(output_dir, 'ImageSummary.csv'), 'WriteRowNames', true);
fprintf("\nOCTA Script 1: Preprocessing & Skin-Border Alignment DONE\n");
