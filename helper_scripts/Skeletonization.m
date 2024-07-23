function [skeletonized_img, binarized_img] = Skeletonization(image, median_filter_size, frangi_opts, thresholding, VISUALIZE, save_path)
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

% pad the image to avoid boundary artefacts, currently hard-coded, may
% move to function parameter
pad = true;
if pad
    image_orig = image;
    % Define padding size
    padSize = 10;  % Adjust the padding size as needed
    % Pad the image symmetrically
    image = padarray(image_orig, [padSize padSize], 'symmetric');
end

% Median-filter the image
image_median = medfilt2(image,[median_filter_size median_filter_size]);

% Apply Frangi filter to image
image_frangi = frangi_2Dfilter(image_median, frangi_opts);

% Crop back if pad == true
if pad
    image_frangi = image_frangi(padSize+1:end-padSize, padSize+1:end-padSize);
end

% 2) SEGMENTATION & SKELETONIZATION
switch thresholding.method
    case 'fuzzy_thresholding'
        % Apply fuzzy thresholding to binarize & segment the image
        binarized_img = fuzzy_thresholding(image_frangi, 2, 3) - 1; % nth = 2 clusters
    case 'local_adaptive_thresholding'
        adaptive_threshold = adaptthresh(image_frangi, thresholding.sensitivity); % 2nd parameter is the sensitivity
        binarized_img = imbinarize(image_frangi, adaptive_threshold);
    case 'otsu_thresholding'
        threshold = graythresh(image_frangi);
        binarized_img = imbinarize(image_frangi, threshold);
    case 'test_all'
        disp("You selected <test_all>. This will generate a plot for " + ...
            "each available thresholding method, including local adaptive " + ...
            "thresholding with 3 different options for the sensitivity (0.2, 0.4, 0.6)")
        binarized_img_fuzzy = fuzzy_thresholding(image_frangi, 2, 3) - 1;
        threshold = graythresh(image_frangi);
        binarized_img_otsu = imbinarize(image_frangi, threshold);
        adaptive_tr1 = adaptthresh(image_frangi, 0.2);
        binarized_img_ad_tr1 = imbinarize(image_frangi, adaptive_tr1);
        adaptive_tr2 = adaptthresh(image_frangi, 0.4);
        binarized_img_ad_tr2 = imbinarize(image_frangi, adaptive_tr2);
        adaptive_tr3 = adaptthresh(image_frangi, 0.6);
        binarized_img_ad_tr3 = imbinarize(image_frangi, adaptive_tr3);
        binarized_img = binarized_img_fuzzy; % Use fuzzy tr. as the default

        if VISUALIZE
        % Create figure for different thresholding results
            Ft = figure(1);
            subplot(2,3,1); imshow(image_median./255); title('Median filtered image');
            subplot(2,3,2); imshow(binarized_img_fuzzy); title('Fuzzy thresholding');
            subplot(2,3,3); imshow(binarized_img_otsu); title('Otsu thresholding');
            subplot(2,3,4); imshow(binarized_img_ad_tr1); title('Adaptive local tr., sensitivity = 0.2');
            subplot(2,3,5); imshow(binarized_img_ad_tr2); title('Adaptive local tr., sensitivity = 0.4');
            subplot(2,3,6); imshow(binarized_img_ad_tr3); title('Adaptive local tr., sensitivity = 0.6');
            % Ft.WindowState = 'maximized';
            save_path2 = strcat(save_path(1:end-4), "_all_thresholding_methods.png");
            set(Ft, 'PaperPositionMode', 'auto');
            print(Ft, save_path2, '-dpng', '-r0', '-painters');
            saveas(gcf, save_path2);
            close(Ft);
        end

    otherwise
        disp(['The selected thresholding method does not exist. ' ...
            'Please select one from the following: {fuzzy_thresholding, ' ...
            'local adaptive thresholding, otsu thresholding}']);
        disp("Use default method <Fuzzy thresholding> instead");
        binarized_img = fuzzy_thresholding(image_frangi, 2, 3) - 1;
end

% Skeletonize the binary image
skeletonized_img = bwmorph(binarized_img, 'skel', Inf);


if VISUALIZE
    % Display the processed image at the different steps
    F1 = figure(1);
    subplot(1,4,1); imshow(image_median./255); title('1. Median filtering');
    subplot(1,4,2); imshow(image_frangi); title('2. Frangi filtering');
    subplot(1,4,3); imshow(binarized_img); title('3. Fuzzy thresholding');
    subplot(1,4,4); imshow(skeletonized_img); title('4. Skeletonization');
    % F1.WindowState = 'maximized';
    save_path2 = save_path(1:end-4) + "_skeletonized_" + join(struct2array(thresholding), "_") + ".png";
    set(F1, 'PaperPositionMode', 'auto');
    print(F1, save_path2, '-dpng', '-r0', '-painters');
    saveas(gcf, save_path2);
    close(F1);
end
end
