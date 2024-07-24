%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 2: Estimation of Capillary Loop Depth (CLD) 
%%%%                           and Superficial Plexus Depth (SPD)
%%%%
%%%% Version:     1.1
%%%% Date:        24/07/2024
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

% ================= Initialization ========================================
close all
clear all
addpath('helper_scripts');

% ================= User Parameters =======================================
% Input image directory containing the cropped OCT images & result directory

result_path = 'results/240723_201053';

% Parameter to adjust the median filter size (higher values = more smoothing)
median_filter_size = 5;

% Parameters for the Frangi filter for vessel enhancement
frangi_opts.sigmarange = [1 6];
frangi_opts.sigmastepsize = 1;
frangi_opts.correctionconst1 = 0.8;
frangi_opts.correctionconst2 = 15;

% The method to apply for thresholding, i.e. dividing the image into foreground (vessels) & background
% The selectable options are: "local_adaptive_thresholding", "fuzzy_thresholding", "otsu_thresholding"
thresholding.method = "local_adaptive_thresholding";

% The sensitivity of the local adaptive thresholding algorithm (higher values will
% pick up more of the vessels & smaller vessels, but potentially also more noise)
% Note: This parameter is only relevant if you selected "local_adaptive_thresholding"
thresholding.sensitivity = 0.3;

% ================= Print parameters & parse input ========================
% Print the selected parameters & image path
fprintf("\nIMAGE FOLDER PATH: \n  <%s>\n", result_path);
fprintf('\nSELECTED PARAMETERS:\n  Median filter size: %d \n  Frangi filter: sigma range: [%d, %d], sigma stepsize: %d, correctionconst1: %.2f, correctionconst 2: %d\n\n', median_filter_size, frangi_opts.sigmarange(1), frangi_opts.sigmarange(2), frangi_opts.sigmastepsize, frangi_opts.correctionconst1, frangi_opts.correctionconst2);

% Add new user parameters to parameters table
params_path = fullfile(result_path, 'InputParameters.mat');
load(params_path, "pixel_size");
save(params_path, "median_filter_size", "frangi_opts", "thresholding", "-append");

% Create new folder for the results of CLD- & SPD-depth computation
CLD_SPD_result_path = fullfile(result_path, '2_CLD_SPD_estimation');
mkdir(CLD_SPD_result_path);

% Parse input directory
fileList = dir(fullfile(result_path, '1_AlignedImages', '*.tif*'));
imgInfo = readtable(fullfile(result_path, 'ImageSummary.csv'), 'ReadRowNames', true, 'Delimiter', ',');

% ================= Compute the CLD- & SPD-depth for all images ===========
% Inizialize empty arrays to store the CLD-depth and SPD-depth per image
CLD_depths = zeros(height(imgInfo), 1);
SPD_depths = zeros(height(imgInfo), 1);

% Iterate over the image folder & compute the CLD- & SPD-depth for each image
fprintf("ESTIMATE THE CLD-DEPTH & SPD-DEPTH:\n")
for ff = 1:length(fileList)
    image_path = fullfile(fileList(ff).folder, fileList(ff).name);
    [CLD_depth, SPD_depth] = CLD_SPD_Estimation(result_path, image_path, median_filter_size, frangi_opts, thresholding, CLD_SPD_result_path, fileList(ff).name, pixel_size(3));
    CLD_depths(ff) = CLD_depth;
    SPD_depths(ff) = SPD_depth;
    CLD_depth_um = (CLD_depth-1)*pixel_size(3);
    if ~isnan(SPD_depth)
        SPD_depth_um = (SPD_depth-1)*pixel_size(3);
    else
        SPD_depth_um = NaN;
    end
    fprintf("  Image: <%s>:  CLD_depth = %.1f,  SPD_depth = %.1f\n", fileList(ff).name, CLD_depth_um, SPD_depth_um);
end

% Add the CLD- & SPD-depth info as a new column to the table
imgInfo.("CLD_Frame") = CLD_depths;
imgInfo.("SPD_Frame") = SPD_depths;
imgInfo.("CLD_Depth_um") = (CLD_depths-1)*pixel_size(3);
imgInfo.("SPD_Depth_um") = (SPD_depths-1)*pixel_size(3);
fprintf("\n"); disp(imgInfo);
table_path = fullfile(result_path, 'ImageSummary.csv');
writetable(imgInfo, table_path, 'WriteRowNames',true);

fprintf("\nOCTA Script 2: Estimation of CLD-depth & SPD-depth DONE\n");


% ================= Function for CLD- & SPD-depth estimation ==============
%% CLD-depth & SPD-depth Estimation
function [CLD_depth, SPD_depth] = CLD_SPD_Estimation(result_path, image_path, median_filter_size, frangi_opts, thresholding, CLD_SPD_result_path, img_name, pixel_size)

% Load the selected image
image = tiffreadVolume(image_path);

% Convert pixels to range [0, 255]
if max(image(:)) <= 1
    image = image * 255;
end
image = double(int16(squeeze(image)));

% Iterate over the depth of the image & compute the skeleton for each x-y-plane
%fprintf("COMPUTE THE NUMBER OF INDEPENDENT BLOOD VESSEL NETWORKS\n")
depth = size(image, 3);
num_segments_array = zeros(1, depth);
VISUALIZE = false;
for z = 1:depth
    depth_slice = image(:, :, z);
    [skeletonized_slice, ~] = Skeletonization(depth_slice, median_filter_size, frangi_opts, thresholding, VISUALIZE, "none");

    % Find connected segments & count the number of independent segments
    segments = bwconncomp(skeletonized_slice);
    num_segments = segments.NumObjects;
    num_segments_array(z) = num_segments;
    %fprintf("  Slice %d: #independent segments = %d\n", z, num_segments);
end

% Smooth the number of segments curve with a moving average filter (in the paper, the window is 20um)
window_size = 20;
smoothed_array = movmean(num_segments_array, window_size);

% Calculate the position of the CLD (maximum number of unconnected
% capillary loops) in first half of image
[~,locs,~,p] = findpeaks(smoothed_array(1:round(depth/2)));
[~, p_idx] = max(p);
CLD_depth = locs(p_idx);

% Calculate the position of the SPD (depth position where vessel network is fully
% connected, i.e. the gradient of the curve of independent segments vs depth reaches 0)
gradient = diff(smoothed_array);
smoothed_gradient = movmean(gradient, 10);
[~, min_grad_idx] = min(smoothed_gradient(CLD_depth:end));
min_grad_idx = CLD_depth + min_grad_idx;
SPD_depth = find(smoothed_gradient(min_grad_idx:end) >= 0, 1, 'first') + min_grad_idx;

% Check if any errors occured for the estimation of the CLD & SPD depth
if isempty(SPD_depth)
    [~, img_name, ~] = fileparts(image_path);
    fprintf('  ERROR: SPD depth is empty for image: <%s>\n', img_name);
    % Check if the error log file already exists - if not, create a new file
    SPD_error_filepath = fullfile(result_path, 'SPD_not_found.txt');
    if exist(SPD_error_filepath) == 2
        fid = fopen(SPD_error_filepath, 'a');
    else
        fid = fopen(SPD_error_filepath, 'w');
        fprintf(fid, ['ERROR: SPD_depth could not be calculated, please ' ...
            'manually select the SPD depth and save it into ImageSummary.csv. ' ...
            '\nThe SPD could not be calculated for the following images:\n\n']);
    end
    fprintf(fid, '%s\n', img_name);
    fclose(fid);

    % Set SPD-depth to NaN for manual correction
    SPD_depth = NaN;
end

% Adapt metrics to micron-scale
CLD_depth_um = round((CLD_depth-1) .* pixel_size, 1);
SPD_depth_um = round((SPD_depth-1) .* pixel_size, 1);
depth_in_microns = (0:depth-1) * pixel_size;

% Create a graph that plots the number of independent segments per depth slice
figure;
subplot(2, 1, 1);  % Subplot for original data
plot(depth_in_microns, num_segments_array, '-', 'LineWidth', 1.5);
xlabel('Depth in microns');
ylabel('#Independent segments');
title('Original Data');
grid on;
hold on; % Adding annotations
plot(CLD_depth_um, num_segments_array(CLD_depth), 'ro', 'MarkerSize', 10);
text(CLD_depth_um, num_segments_array(CLD_depth), sprintf('  CLD depth: %.1f', CLD_depth_um), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
if ~isnan(SPD_depth)
    plot(SPD_depth_um, num_segments_array(SPD_depth), 'ro', 'MarkerSize', 10);
    text(SPD_depth_um, num_segments_array(SPD_depth), sprintf('  SPD depth: %.1f', SPD_depth_um), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end
hold off;

subplot(2, 1, 2);  % Subplot for smoothed data
plot(depth_in_microns, smoothed_array, '-', 'LineWidth', 1.5);
xlabel('Depth in microns');
ylabel('#Independent segments');
title(['Smoothed Data (Window Size: ', num2str(window_size), ')']);
grid on;
hold on; % Adding annotations
plot(CLD_depth_um, smoothed_array(CLD_depth), 'ro', 'MarkerSize', 10);
text(CLD_depth_um, smoothed_array(CLD_depth), sprintf('  CLD depth: %.1f', CLD_depth_um), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
if ~isnan(SPD_depth)
    plot(SPD_depth_um, smoothed_array(SPD_depth), 'ro', 'MarkerSize', 10);
    text(SPD_depth_um, smoothed_array(SPD_depth), sprintf('  SPD depth: %.1f', SPD_depth_um), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end
hold off;

[~, fileName, ~] = fileparts(img_name);
saveas(gcf, fullfile(CLD_SPD_result_path, strcat(img_name, '_Independent_segments_per_depth_graph.png')));
close

% Interpolate data to micron-depth scale
new_micron_depth = round(linspace(0, (depth-1) * pixel_size, (depth-1) * pixel_size + 1)');
num_segments_interp = interp1(depth_in_microns, num_segments_array, new_micron_depth, 'linear');
smoothed_segments_interp = interp1(depth_in_microns, smoothed_array, new_micron_depth, 'linear');

% Save independent segments per depth array (original & smoothed) to a csv-file.
% Additionally, save the estimated CLD- & SPD-depth in separate columns
CLD_info_col = zeros(length(num_segments_interp), 1); CLD_info_col(1) = CLD_depth_um;
SPD_info_col = zeros(length(num_segments_interp), 1); SPD_info_col(1) = SPD_depth_um;
data_table = table(new_micron_depth, num_segments_interp, smoothed_segments_interp, CLD_info_col, SPD_info_col, 'VariableNames', {'Depth_um', 'num_segments', 'smoothed_num_segments', 'CLD_depth', 'SPD_depth'});
file_path = fullfile(CLD_SPD_result_path, strcat(fileName, '_Independent_segments_per_depth_data.csv'));
writetable(data_table, file_path);
%fprintf("\nNumber of segments per depth data saved to <%s>\n", file_path);
end
