function [skeletonized_img, binarized_img] = Skeletonization(image, median_filter_size, frangi_opts, VISUALIZE, save_path)
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

% 
% options = struct('FrangiScaleRange', [1 6], 'FrangiScaleRatio', 2, ...
%                  'FrangiBetaOne', 0.5, 'FrangiBetaTwo', 15, ...
%                  'BlackWhite', true); % Adjust these parameters as needed
% 
% [J, Scale, Direction] = FrangiFilter2D(image);
% 
% figure;
% subplot(1, 2, 1);
% imshow(image_median./255);
% title('Original Image');
% subplot(1, 2, 2);
% imshow(J, []);
% title('Frangi Filtered Image');
% level = graythresh(image_frangi); % Compute Otsu's threshold level
% binarized_img_otsu = imbinarize(image_frangi, level); % Binarize the image using the computed threshold


% Skeletonize the binary image
skeletonized_img = bwmorph(binarized_img, 'skel', Inf);

% image_frangi = (image_frangi - min(image_frangi(:))) / (max(image_frangi(:)) - min(image_frangi(:)));

% imwrite(uint8(image_median), 'image_median.png');
% imwrite(uint8(image_frangi* 255), 'image_frangi.png');
% imwrite(uint8(binarized_img * 255), 'binarized_img.png');

if VISUALIZE
    % Display the processed image at the different steps
    F1 = figure(1);

    subplot(1,4,1); imshow(image_median./255); title('1. Median filtering');
    subplot(1,4,2); imshow(image_frangi); title('2. Frangi filtering');
    subplot(1,4,3); imshow(binarized_img); title('3. Fuzzy thresholding');
    subplot(1,4,4); imshow(skeletonized_img); title('4. Skeletonization');

    % subplot(2,4,1); imshow(image./255); title('1. Original image');
    % subplot(2,4,2); imshow(image_median./255); title('2. Median filtering');
    % subplot(2,4,3); imshow(image_frangi); title('3. Frangi filtering');
    % subplot(2,4,4); imshow(binarized_img); title('4. Fuzzy thresholding');
    % subplot(2,4,5); imshow(skeletonized_img); title('5. Skeletonization');

    % subplot(2,4,6); imshow(binarized_img_otsu); title('6. binarized_img_otsu');

    F1.WindowState = 'maximized';

    save_path2 = save_path(1:end-4) + "_skeletonized.png";
    set(F1, 'PaperPositionMode', 'auto');
    print(F1, save_path2, '-dpng', '-r0', '-painters');
    %saveas(gcf, save_path2);
    close(F1);
end
end
