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
input_dir = 'sample_data';

% Result directory & optional result
% result_dir = 'N:\SUN-BMI-DBI-DATA-NOBACKUP\NinaMK\Results';
result_dir = 'results';

% Optional results to write
writeZDisplacement = true;

% Automated image cropping
AutoCrop = false; %true;

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
        [CroppedOCTStack, crop_range] = AutoCropOCTStack(OCTStack);
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
T = horzcat(cell2table(fileinfo), array2table(fileslice));
T.Properties.RowNames = {filelist.name};
T.Properties.VariableNames = {'XRes', 'YRes', 'ResUnit', ...
                               'First Slice', 'Last Slice'};
writetable(T, fullfile(output_dir, 'ImageSummary.csv'), 'WriteRowNames', true);
fprintf("\nOCTA Script 1: Preprocessing & Skin-Border Alignment DONE\n");
