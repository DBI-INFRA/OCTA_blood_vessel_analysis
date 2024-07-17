function [skeletonized_img, binarized_img] = Skeletonization(image, median_filter_size, frangi_opts, VISUALIZE)
% ---------------------------------------------------------------------------
% Description:
%    This function segments & skeletonizes the blood vessel network of the
%    input image. Before thresholding, it performs some filtering with a
%    median and frangi filter to reduce noise and enhance vessel structures
%
% Parameters:
%    - image: Input 2D image
%    - median_filter_size: Size of the median filter to be applied
%    - frangi_opts: Structure containing the following filter options:
%        - sigmarange & sigmastepsize
%        - correctionconst1 & correctionconst2
%    - VISUALIZE: Boolean to indicate whether to visualize intermediate steps
%
% Returns:
%    - skeletonized_img: Skeletonized binary image of the blood vessels
%    - binarized_img: Binarized image after fuzzy thresholding
% ---------------------------------------------------------------------------

% Import required functions & define path to scripts
% addpath('helper_scripts');
import fuzzylogic.*
set(0,'DefaultFigureWindowStyle' , 'normal');

% 1) PREPROCESSING
% Median-filter the image
image_median = medfilt2(image,[median_filter_size median_filter_size]);

% Apply Frangi filter to image
image_frangi = frangi_2Dfilter(image_median, frangi_opts);

% 2) SEGMENTATION & SKELETONIZATION
% Apply fuzzy thresholding to binarize & segment the image
binarized_img = fuzzy_thresholding(image_frangi, 2, 3) - 1; % nth = 2 clusters

% Skeletonize the binary image
skeletonized_img = bwmorph(binarized_img, 'skel', Inf);

if VISUALIZE
    % Display the processed image at the different steps
    F1 = figure(1);
    subplot(2,3,1); imshow(image); title('Original image');
    subplot(2,3,2); imshow(image_median./255); title('2. Image after median filtering');
    subplot(2,3,3); imshow(image_frangi); title('3. Image after frangi filtering');
    subplot(2,3,4); imshow(binarized_img); title('4. Image after fuzzy thresholding');
    subplot(2,3,5); imshow(skeletonized_img); title('5. Image after skeletonization');
    F1.WindowState = 'maximized';
end
end
