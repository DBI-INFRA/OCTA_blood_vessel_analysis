function [CLD_depth, SPD_depth] = CLD_SPD_estimation(result_path, image_path, median_filter_size, frangi_opts, thresholding, CLD_SPD_result_path, img_name, pixel_size)

% -------------------------------------------------------------------------
% Description:
%   This function computes:
%     - CLD (Capillary Loop Density) depth: the depth index at which the
%       number of independent vessel segments (unconnected capillary loops)
%       reaches a local maximum in the first half of the volume.
%     - SPD (Superficial Plexus Depth): the depth index where the vessel
%       network becomes fully connected again (determined from the gradient
%       of the smoothed number-of-segments curve).

%   The function reads the input volume using `tiffreadVolume`, normalizes
%   intensities to [0,255], computes a skeleton per depth slice using
%   `skeletonization`, counts connected components per slice, smooths the
%   resulting curve, locates the CLD and SPD indices, produces summary
%   plots, and writes a CSV file with per-micron interpolated data.

% Parameters:
%   - result_path:          path where high-level results / logs are stored
%   - image_path:           full path to the input multi-page TIFF or volume
%   - median_filter_size:   size (in pixels) of the median filter used in preprocessing
%   - frangi_opts:          options passed to the Frangi vesselness filter (see `frangi_2Dfilter`)
%   - thresholding:         thresholding settings used by `skeletonization`
%   - CLD_SPD_result_path:  directory where the per-image CLD/SPD CSV and plots are saved
%   - img_name:             filename (or identifier) used for saving outputs
%   - pixel_size:           size of one pixel in microns (um/pixel)
%
% Returns:
%   - CLD_depth:            estimated CLD depth expressed as a 1-based image slice index
%   - SPD_depth:            estimated SPD depth expressed as a 1-based image slice index (NaN if not found)
%
% Notes:
%   - The returned `CLD_depth` and `SPD_depth` values are returned as image
%     slice indices (1-based). If you need micron values, the function also
%     saves micron-scaled columns (`Depth_um`) and writes the CLD/SPD in
%     microns to the CSV file saved in `CLD_SPD_result_path`.
%   - If SPD cannot be determined the function logs the image name to
%     `SPD_not_found.txt` within `result_path` and returns `SPD_depth = NaN`.
%
% Authors: Tricia Loo (DBI-Infra IACF, tricia.loo@sund.ku.dk)
%          Julia Mertesdorf (DBI-Infra IACF, jume@di.ku.dk)
% Date:    26/02/2025


% Load the selected image
image = tiffreadVolume(image_path);
[~, fileName, ~] = fileparts(img_name);

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

SPD_range_um = 10;
SPD_range   = round(SPD_range_um/pixel_size);

for z = 1+SPD_range: depth-SPD_range
    
    % depth_slice = image(:, :, z);

    depth_slice = mean(image(:, :, z-SPD_range: z+SPD_range), 3);

    [skeletonized_slice, ~] = skeletonization(depth_slice, median_filter_size, frangi_opts, thresholding, VISUALIZE, "none");
    % Find connected segments & count the number of independent segments
    segments = bwconncomp(skeletonized_slice);
    num_segments = segments.NumObjects;
    num_segments_array(z) = num_segments;
    %fprintf("  Slice %d: #independent segments = %d\n", z, num_segments);
end

% imwrite(skeleton_results_image, skeleton_results_image_path); 

num_segments_array = num_segments_array(SPD_range+1: end-SPD_range);

% Smooth the number of segments curve with a moving average filter
window_size = round(20 / pixel_size);
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

CLD_depth = CLD_depth + SPD_range;
if ~isnan(SPD_depth)
    SPD_depth = SPD_depth + SPD_range;
end

% Adapt metrics to micron-scale
CLD_depth_um = round((CLD_depth-1) .* pixel_size, 1);
SPD_depth_um = round((SPD_depth-1) .* pixel_size, 1);
depth_in_microns = (0:depth-1) * pixel_size;


num_segments_array = [ones(1, SPD_range) * num_segments_array(1), num_segments_array, ones(1, SPD_range) * num_segments_array(end)];
smoothed_array = [ones(1, SPD_range) * smoothed_array(1), smoothed_array, ones(1, SPD_range) * smoothed_array(end)];

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