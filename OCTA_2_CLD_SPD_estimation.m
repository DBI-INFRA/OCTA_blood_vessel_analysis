%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 2: Estimation of Capillary Loop Depth (CLD) 
%%%%                           and Superficial Plexus Depth (SPD)
%%%%
%%%% Version:     1.0
%%%% Date:        16/07/2024
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
result_path = 'results/240719_141653';

% Parameters for skeletonization (median-filter & Frangi-filter)
median_filter_size = 5;
frangi_opts.sigmarange = [1 6];
frangi_opts.sigmastepsize = 2;
frangi_opts.correctionconst1 = 0.8;
frangi_opts.correctionconst2 = 15;
thresholding_method = "fuzzy_thresholding";  % OPTIONS: "fuzzy_thresholding", "local_adaptive_thresholding"
VISUALIZE = false;

% ================= Print parameters & parse input ========================
% Print the selected parameters & image path
fprintf("\nIMAGE FOLDER PATH: \n  <%s>\n", result_path);
fprintf('\nSELECTED PARAMETERS:\n  Median filter size: %d \n  Frangi filter: sigma range: [%d, %d], sigma stepsize: %d, correctionconst1: %.2f, correctionconst 2: %d\n\n', median_filter_size, frangi_opts.sigmarange(1), frangi_opts.sigmarange(2), frangi_opts.sigmastepsize, frangi_opts.correctionconst1, frangi_opts.correctionconst2);

% Create new folder for the results of CLD- & SPD-depth computation
CLD_SPD_result_path = fullfile(result_path, '2_CLD_SPD_estimation');
mkdir(CLD_SPD_result_path);

% Parse input directory
fileList = dir(fullfile(result_path, '1_AlignedImages', '*.tif*'));
imgInfo = readtable(fullfile(result_path, 'ImageSummary.csv'), 'ReadRowNames', true);

% ================= Compute the CLD- & SPD-depth for all images ===========
% Inizialize empty arrays to store the CLD-depth and SPD-depth per image
CLD_depths = zeros(height(imgInfo), 1);
SPD_depths = zeros(height(imgInfo), 1);

% Iterate over the image folder & compute the CLD- & SPD-depth for each image
fprintf("ESTIMATE THE CLD-DEPTH & SPD-DEPTH:\n")
for ff = 1:length(fileList)
    image_path = fullfile(fileList(ff).folder, fileList(ff).name);
    [CLD_depth, SPD_depth] = CLD_SPD_Estimation(image_path, median_filter_size, frangi_opts, thresholding_method, VISUALIZE, CLD_SPD_result_path, fileList(ff).name);
    CLD_depths(ff) = CLD_depth;
    SPD_depths(ff) = SPD_depth;
    fprintf("  Image: <%s>:  CLD_depth = %d,  SPD_depth = %d\n", fileList(ff).name, CLD_depth, SPD_depth);
end

% Add the CLD- & SPD-depth info as a new column to the table
imgInfo.("CLD_Depth") = CLD_depths;
imgInfo.("SPD_Depth") = SPD_depths;
fprintf("\n"); disp(imgInfo);
table_path = fullfile(result_path, 'ImageSummary.csv');
writetable(imgInfo, table_path, 'WriteRowNames',true);

fprintf("\nOCTA Script 2: Estimation of CLD-depth & SPD-depth DONE\n");


% ================= Function for CLD- & SPD-depth estimation ==============
%% CLD-depth & SPD-depth Estimation
function [CLD_depth, SPD_depth] = CLD_SPD_Estimation(image_path, median_filter_size, frangi_opts, thresholding_method, VISUALIZE, CLD_SPD_result_path, img_name)

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
for z = 1:depth
    depth_slice = image(:, :, z);
    [skeletonized_slice, ~] = Skeletonization(depth_slice, median_filter_size, frangi_opts, thresholding_method, VISUALIZE, "none");

    % Find connected segments & count the number of independent segments
    segments = bwconncomp(skeletonized_slice);
    num_segments = segments.NumObjects;
    num_segments_array(z) = num_segments;
    %fprintf("  Slice %d: #independent segments = %d\n", z, num_segments);
end

% Smooth the number of segments curve with a moving average filter (in the paper, the window is 20um)
window_size = 20;
smoothed_array = movmean(num_segments_array, window_size);

% Calculate the position of the CLD (maximum number of unconnected capillary loops)
[~, max_idx] = max(smoothed_array);
CLD_depth = max_idx;

% Calculate the position of the SPD (depth position where vessel network is fully
% connected, i.e. the gradient of the curve of independent segments vs depth reaches 0)
gradient = diff(smoothed_array);
smoothed_gradient = movmean(gradient, 10);
[~, min_grad_idx] = min(smoothed_gradient(CLD_depth:end));
min_grad_idx = CLD_depth + min_grad_idx;
SPD_depth = find(smoothed_gradient(min_grad_idx:end) >= 0, 1, 'first') + min_grad_idx;

% Check if any errors occured for the estimation of the CLD & SPD depth
if isempty(SPD_depth) || isempty(CLD_depth) 
    error('CLD depth or SPD_depth is empty. Please check if the input is correct');
end

% Create a graph that plots the number of independent segments per depth slice
figure;
subplot(2, 1, 1);  % Subplot for original data
plot(1:depth, num_segments_array, '-', 'LineWidth', 1.5);
xlabel('Depth');
ylabel('#Independent segments');
title('Original Data');
grid on;
hold on; % Adding annotations
plot(CLD_depth, num_segments_array(CLD_depth), 'ro', 'MarkerSize', 10);
text(CLD_depth, num_segments_array(CLD_depth), sprintf('  CLD depth: %d', CLD_depth), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(SPD_depth, num_segments_array(SPD_depth), 'ro', 'MarkerSize', 10);
text(SPD_depth, num_segments_array(SPD_depth), sprintf('  SPD depth: %d', SPD_depth), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
hold off;

subplot(2, 1, 2);  % Subplot for smoothed data
plot(1:depth, smoothed_array, '-', 'LineWidth', 1.5);
xlabel('Depth');
ylabel('#Independent segments');
title(['Smoothed Data (Window Size: ', num2str(window_size), ')']);
grid on;
hold on; % Adding annotations
plot(CLD_depth, smoothed_array(CLD_depth), 'ro', 'MarkerSize', 10);
text(CLD_depth, smoothed_array(CLD_depth), sprintf('  CLD depth: %d', CLD_depth), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(SPD_depth, smoothed_array(SPD_depth), 'ro', 'MarkerSize', 10);
text(SPD_depth, smoothed_array(SPD_depth), sprintf('  SPD depth: %d', SPD_depth), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
hold off;

[~, fileName, ~] = fileparts(img_name);
saveas(gcf, fullfile(CLD_SPD_result_path, strcat(img_name, '_Independent_segments_per_depth_graph.png')));

close

% Save independent segments per depth array (original & smoothed) to a csv-file.
% Additionally, save the estimated CLD- & SPD-depth in separate columns
CLD_info_col = zeros(length(num_segments_array), 1); CLD_info_col(1) = CLD_depth;
SPD_info_col = zeros(length(num_segments_array), 1); SPD_info_col(1) = SPD_depth;
data_table = table((1:depth)', num_segments_array', smoothed_array', CLD_info_col, SPD_info_col, 'VariableNames', {'Depth', 'num_segments', 'smoothed_num_segments', 'CLD_depth', 'SPD_depth'});
file_path = fullfile(CLD_SPD_result_path, strcat(fileName, '_Independent_segments_per_depth_data.csv'));
writetable(data_table, file_path);
%fprintf("\nNumber of segments per depth data saved to <%s>\n", file_path);
end
