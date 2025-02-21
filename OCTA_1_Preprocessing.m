%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 1: Preprocessing & Skin-Border Alignment
%%%%
%%%% Version:     1.2
%%%% Date:        26/07/2024
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

% Input/Output Directory
input_dir = 'sample_data';  % Input image directory containing the original 
                            % OCT images. (Note: The image data should be
                            % in a local drive for fast read/write speed)

result_dir = 'results';     % Result directory, will be made if one does 
                            % not exist.

                            
% Pixel size in micrometers
pixel_size = [6.49, 6.49, 2.45];


% Automated image cropping in Z-direction
AutoCrop = false;           %% AutoCrop is used to remove low signal data 
                            %  from the front and enf of the stacj
                            %  The sample data has already been cropped to 
                            %  size and should have AutoCrop set to false. 
                            %  To enable, set AutoCrop to true.

CropSensitivity = 0.5;      %% Select a value between [0, 1] to control the
                            %  amount of cropping from the end of the stack. 
                            %  0 means more data will be kept and 1 means 
                            %  more data will be cropped. 
                            %  Can be set to 'None' to prevent cropping from
                            %  end of stack.


% Wavelet Filtering Parameters
orientation = 'both';  	    %% Direction of stripe artefact in image, choose
                            %  between 'horizontal', 'vertical' or 'both' 
dimension = 'z';            %% Direction to apply filter if image is 3D 
                            %  stack, choose between 'x', 'y' or 'z'.
                            %  'z' means that wavelet filter
% Additional parameters for wavelet filtering            
gamma = 10; order = 4;      %  Not recommended to change
wname = 'db20'; reps = 2;   %  See SuppWaveletFFT.m for details
                           

% Image Alignment to air-skin boundary
AutoAlign = true;           % Option to align image stack to air-skin boundary
writeZDisplacement = true;  % Option to write the z-displacement for image into a tiff file

%% ================= Parse input and preallocate storage ==================

filelist  = dir(fullfile(input_dir, '*.tif*'));

fileslice = ones(length(filelist), 2);
timenow = char(datetime("now"), "yyMMdd_HHmmss");
output_dir  = fullfile(result_dir, timenow);
imwrite_dir = fullfile(output_dir, '1_AlignedImages');
mkdir(imwrite_dir);

% Save all user parameters in a table
params_path = fullfile(output_dir, 'InputParameters.mat');
preprocess_opts = struct("AutoCrop", AutoCrop, "CropSensitivity", CropSensitivity, ...
                    "WT_gamma", gamma, "WT_order", order, "WT_wname", wname, ...
                    "WT_reps", reps, "WT_orient", orientation, "WT_dim", dimension);
save(params_path, "pixel_size", "preprocess_opts");

%% ================= Preprocess all images in input folder ================
for ff = 1:length(filelist)
    img_path = fullfile(input_dir, filelist(ff).name);
    img_info = imfinfo(img_path);

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
                                 gamma, order, wname, orientation, reps, dimension);

    % Z alignment of image
    if AutoAlign
        [ResultStack, z_shift] = ZAlignStack(FilteredStack, 'median');
    else
        ResultStack = double(FilteredStack);
    end

    % 2) Write Results to file
    % Write aligned image
    imwrite(ResultStack(:,:,1), fullfile(imwrite_dir, filelist(ff).name));
    for ii = 2:size(ResultStack, 3)
         imwrite(ResultStack(:,:,ii), fullfile(imwrite_dir, filelist(ff).name), "Writemode", "append");
    end
    
    % Write zdisplacement
    if writeZDisplacement & AutoAlign
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
