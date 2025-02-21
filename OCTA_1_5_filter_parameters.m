
result_path = 'results/250221_134531';

% Parameter to adjust the median filter size (higher values = more smoothing)
median_filter_size = 5;

% Parameters for the Frangi filter for vessel enhancement
frangi_opts.sigmarange = [1 6];
frangi_opts.sigmastepsize = 1;
frangi_opts.correctionconst1 = 0.8;
frangi_opts.correctionconst2 = 15; 

% Morphological opening will remove some small islands after segmentation, based on 
% the area size (in pixels). This helps to remove noise. Adjust if necessary
thresholding.opening_size = 75;

% ================= Print parameters & parse input ========================
% Print the selected parameters & image path
fprintf("\nIMAGE FOLDER PATH: \n  <%s>\n", result_path);
fprintf('\nSELECTED PARAMETERS:\n  Median filter size: %d \n  Frangi filter: sigma range: [%d, %d], sigma stepsize: %d, correctionconst1: %.2f, correctionconst 2: %d\n\n', median_filter_size, frangi_opts.sigmarange(1), frangi_opts.sigmarange(2), frangi_opts.sigmastepsize, frangi_opts.correctionconst1, frangi_opts.correctionconst2);

% IMPORTANT NOTE: The thresholding method + sensitivity is now fixed for OCTA_2, since this
%                 variant worked best for accurate CLD & SPD detection.
% The method to apply for thresholding, i.e. dividing the image into foreground (vessels) & background
% The selectable options are: "local_adaptive_thresholding", "fuzzy_thresholding", "otsu_thresholding"
thresholding.method = "local_adaptive_thresholding";

% The sensitivity of the local adaptive thresholding algorithm (higher values will
% pick up more & smaller vessels, but potentially also more noise).
% IMPORTANT: The thresholding sensitivity should be set to 0 for OCTA_2 to improve  
% the accuracy of the CLD-detection. Note that it is however not necessarily recommended  
% to use sensitivity=0 for OCTA_3 - the ideal value should be selected based on visual 
% inspection of the thresholding results for different sensitivity values.
thresholding.sensitivity = 0;

% Add new user parameters to parameters table
params_path = fullfile(result_path, 'InputParameters.mat');
load(params_path, "pixel_size");
save(params_path, "median_filter_size", "frangi_opts", "thresholding", "-append");
