%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow configuration file
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========== Image Information and Directories ===========================
image_dir = 'sample_data';  %% Image directory containing .tif files of 
                            % OCT images to be processed.

output_dir = 'results';     %% Result directory will be made if it does not exist

pixel_size  = [6.49, 6.49, 2.45];

%% ========== Image Pre-processing Parameters (OCTA_1 only) ==================
preprocess_image = true;    %% Set to false to ignore all image preprocessing
                            %  parameters in this section

% Automated image cropping in Z-direction
AutoCrop = false;           %% AutoCrop is used to remove low signal data 
                            %  from the front and end of the stack
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
WT_orient = 'both';  	%% Direction of stripe artefact in image, choose
                        %  between 'horizontal', 'vertical' or 'both' 
WT_dim = 'z';           %% Direction to apply filter if image is 3D stack, 
                        %  choose between 'x', 'y' or 'z'.
                        %  'z' means that wavelet filter is applied
                        %  slice-by-slice in the 3rd dimension
WT_reps = 2;            %% Number of filtering repetitions, higher values
                        % will increase amounts of filters.

% Image Alignment to air-skin boundary
AutoAlign = true;           % Option to align image stack to air-skin boundary
WriteZDisplacement = false; % Option to write the z-displacement for image into a tiff file

%% ========== Vessel Enhancement Parameters (OCTA_2 AND 3) ================
% Parameter to adjust the median filter size (higher values = more smoothing)
median_filter_size = 5;

% Parameters for the Frangi filter for vessel enhancement
frangi_opts.sigmarange = [1 6];
frangi_opts.sigmastepsize = 1;
frangi_opts.correctionconst1 = 0.8;
frangi_opts.correctionconst2 = 15; 

%% ========== Vessel Quantification Parameters (OCTA_3 only) ==============
% Thresholding parameters, i.e. dividing the image into foreground (vessels) & background
thresholding.method = "local_adaptive_thresholding";
                                %%The selectable methods are: 
                                % "local_adaptive_thresholding", 
                                % "fuzzy_thresholding", 
                                % "otsu_thresholding", or "test_all"

thresholding.sensitivity = 0.3; %%The sensitivity of the local adaptive 
                                % thresholding algorithm (higher values 
                                % will pick up more & smaller vessels, 
                                % but potentially also more noise). If 
                                % method = 'test_all', multiple
                                % sensitivity levels can be given as a
                                % 1D array (e.g. [0.1, 0.2, 0.3, 0.4]).
                                % If method = 'local_adaptive_thresholding', 
                                % sensitivity will be first value given.

thresholding.opening_size = 75; %% Morphological opening removes some 
                                % small islands after segmentation, 
                                % based on area size (in pixels).

quantify_frame = 'SPD_Frame';   %% Set central frame for quantification. 
                                % Value can be either SINGLE FRAME NUMBER, or a 
                                % VARIABLE NAME in 'ImageSummary.csv'.
                                %
                                % E.g. To automatically measure morphology 
                                % around the SPD frame calculated in OCTA_2, 
                                % set quantify_frame = 'SPD_Frame'; to 
                                % automatically measure a set of user- 
                                % specified frames, set quantify_frame = 
                                % 'My_Frame' and create a column
                                % titled 'My_Frame' in 'ImageSummary.csv' 
                                % with frames to use for each image.

quantify_range_um = 30;         %% Depth around the selected frame to 
                                % include in quantification (microns) 

save_detailed_results = true; 	%% Saves detailed morphological results 
                                % including image of vessel skeleton,
                                % binary image of vessel, image of labeled
                                % skeleton, and table of per-branch
                                % measurements.

%% ================= Save configuration parameters ========================
% Create results folder
timenow = char(datetime("now"), "yyMMdd_HHmmss");
output_dir = fullfile(output_dir, timenow);
mkdir(output_dir)

% Create image summary file
filelist  = dir(fullfile(image_dir, '*.tif*'));
T = array2table({filelist.name}', "VariableNames", {'Filename'});
writetable(T, fullfile(output_dir, 'ImageSummary.csv'), 'WriteRowNames', true);

% Create parameter file as .m
params_path = fullfile(output_dir, 'Parameters.mat');
save(params_path, "image_dir", "output_dir", "pixel_size", "filelist")
if preprocess_image
    prep_opts = struct("AutoCrop", AutoCrop, "CropSensitivity", CropSensitivity, ...
                "WT_orient", WT_orient, "WT_dim", WT_dim, "WT_reps", WT_reps, ...
                "AutoAlign", AutoAlign, "WriteZDisplacement", WriteZDisplacement);
    save(params_path, "prep_opts", "-append");
end
save(params_path, "median_filter_size", "frangi_opts", "thresholding", "-append");
save(params_path, "quantify_frame", "quantify_range_um", "save_detailed_results", "-append");