%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 2: Estimation of Capillary Loop Depth (CLD) 
%%%%                           and Superficial Plexus Depth (SPD)
%%%%
%%%% Version:     1.3
%%%% Date:        26/02/2025
%%%%
%%%% Authors:     Tricia Loo (DBI-Infra IACF, tricia.loo@sund.ku.dk)
%%%%              Julia Mertesdorf (DBI-Infra IACF, jume@di.ku.dk)
%%%%
%%%% Description: This second part of the pipeline estimates the Capillary  
%%%%              Loop Depth (CLD) and Superficial Plexus Depth (SPD) from  
%%%%              z-aligned Optical Coherence Tomography (OCT) images.
%%%%              The SPD & CLD depth are calculated by first skeletonizing  
%%%%              and then counting the number of independent blood vessel 
%%%%              networks per depth slice. The CLD-depth is defined as the  
%%%%              maximum number of unconnected capillary loops. The 
%%%%              SPD-depth is the depth position where the vessel network   
%%%%              is fully connected, i.e. the depth at which the gradient 
%%%%              of the curve of independent vessel networks reaches 0 for 
%%%%              the first time after the CLD-depth position. The script 
%%%%              adds the CLD- & SPD- depth as new columns to the existing 
%%%%              immage summary table (Note: Both depth numbers are given
%%%%              w.r.t the cropped & skin-aligned OCT-images, not the
%%%%              original ones). Furthermore, it saves the full results
%%%%              of the number of independent blood vessel networks per
%%%%              depth slice both in a csv-table and as a graph (the
%%%%              estimated CLD- & SPD-depths are marked in red).
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ====================== Initialization ==================================
close all
clear all
addpath('helper_scripts');

%%% Specify result_dir (path where the desired 'Parameters.mat' is stored)
result_dir = 'results\yymmdd_hhmmss\';

%%% Load exisiting user configurations
load(fullfile(result_dir, 'Parameters.mat'));

%%% Modify any existing configuration below by setting:
% ParameterName = newValue;
% e.g. median_filter_size = 2;

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
    
    % Default (optimised) thresholding parameters for CLD- & SPD- detection
    thresholding.opening_size = 75;
    thresholding.method = "local_adaptive_thresholding";
    thresholding.sensitivity = 0;
    
    % Preallocate storage
    CLD_depths = zeros(height(imgInfo), 1);
    SPD_depths = zeros(height(imgInfo), 1);
    CLD_SPD_result_path = fullfile(result_dir, 'CLD_SPD_estimation');

    % Print user parameters & image path
    fprintf("\nIMAGE FOLDER PATH: \n  <%s>\n", image_dir);
    fprintf("\nOUTPUT FOLDER PATH: \n  <%s>\n", CLD_SPD_result_path);
    fprintf(['\nSELECTED PARAMETERS:\n' ...
             '  Median filter size: %d \n' ...
             '  Frangi filter: sigma range: [%d, %d], ' ...
             'sigma stepsize: %d, correctionconst1: %.2f, correctionconst 2: %d\n\n'], ...
             median_filter_size, frangi_opts.sigmarange(1), frangi_opts.sigmarange(2), ...
             frangi_opts.sigmastepsize, frangi_opts.correctionconst1, frangi_opts.correctionconst2);
catch ME
    fprintf("Check configuration file: ")
    rethrow(ME);
end

% Create new folder for the results of CLD- & SPD-depth computation
mkdir(CLD_SPD_result_path);

%% ================ Compute the CLD- & SPD-depth for all images ===========
% Iterate over the image folder & compute the CLD- & SPD-depth for each image
fprintf("ESTIMATE THE CLD-DEPTH & SPD-DEPTH:\n")
for ff = 1:length(filelist)
    image_path = fullfile(filelist(ff).folder, filelist(ff).name);
    [CLD_depth, SPD_depth] = CLD_SPD_estimation(result_dir, image_path, median_filter_size, frangi_opts, thresholding, CLD_SPD_result_path, filelist(ff).name, pixel_size(3));
    CLD_depths(ff) = CLD_depth;
    SPD_depths(ff) = SPD_depth;
    CLD_depth_um = (CLD_depth-1)*pixel_size(3);
    if ~isnan(SPD_depth)
        SPD_depth_um = (SPD_depth-1)*pixel_size(3);
    else
        SPD_depth_um = NaN;
    end
    fprintf("  Image: <%s>:  CLD_depth = %.1f,  SPD_depth = %.1f\n", filelist(ff).name, CLD_depth_um, SPD_depth_um);
end

% Add the CLD- & SPD-depth info as a new column to the table
imgInfo.("CLD_Frame") = CLD_depths;
imgInfo.("SPD_Frame") = SPD_depths;
imgInfo.("CLD_Depth_um") = (CLD_depths-1)*pixel_size(3);
imgInfo.("SPD_Depth_um") = (SPD_depths-1)*pixel_size(3);
fprintf("\n"); disp(imgInfo);
writetable(imgInfo, fullfile(result_dir, 'ImageSummary.csv'), 'WriteRowNames', true);

fprintf("\nOCTA Script 2: Estimation of CLD-depth & SPD-depth DONE\n\n");