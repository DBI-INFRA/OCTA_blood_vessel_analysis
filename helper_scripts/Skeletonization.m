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
    padSize = 20;  % Adjust the padding size as needed
    % Pad the image symmetrically
    image = padarray(image_orig, [padSize padSize], 'symmetric');
end

% Median-filter the image
image_median = medfilt2(image,[median_filter_size median_filter_size]);

% Apply Frangi filter to image
image_frangi = frangi_2Dfilter(image_median, frangi_opts);

% Crop back if pad == true
if pad
    image = image(padSize+1:end-padSize, padSize+1:end-padSize);
    image_median = image_median(padSize+1:end-padSize, padSize+1:end-padSize);
    image_frangi = image_frangi(padSize+1:end-padSize, padSize+1:end-padSize);
end

[path, img_name, ~] = fileparts(save_path);

% 2) SEGMENTATION & SKELETONIZATION
switch thresholding.method
    case 'fuzzy_thresholding'
        % Apply fuzzy thresholding to binarize & segment the image
        binarized_img = fuzzy_thresholding(image_frangi, 2, 3) - 1; % nth = 2 clusters
    case 'local_adaptive_thresholding'
        adaptive_threshold = adaptthresh(image_frangi, thresholding.sensitivity(1)); % 2nd parameter is the sensitivity
        binarized_img = imbinarize(image_frangi, adaptive_threshold);
    case 'otsu_thresholding'
        threshold = graythresh(image_frangi);
        binarized_img = imbinarize(image_frangi, threshold);
    case 'test_all'
        sensitivities_str = sprintf('%.1f, ', thresholding.sensitivity);
        fprintf("  You selected <test_all>. This will generate a plot for " + ...
            "each available thresholding method, including local adaptive " + ...
            "thresholding with the following values for the sensitivity: %s\n", ...
            sensitivities_str);
        % Apply fuzzy & otsu thresholding
        binarized_img_fuzzy = fuzzy_thresholding(image_frangi, 2, 3) - 1;
        threshold = graythresh(image_frangi);
        binarized_img_otsu = imbinarize(image_frangi, threshold);

        % Apply adaptive thresholding for each chosen sensitivity value
        num_sensitivities = numel(thresholding.sensitivity);
        binarized_images_adaptive = cell(1, num_sensitivities);
        for i = 1:num_sensitivities
            adaptive_tr = adaptthresh(image_frangi, thresholding.sensitivity(i));
            binarized_images_adaptive{i} = imbinarize(image_frangi, adaptive_tr);
        end
        
        % Create figure for different thresholding results
        if VISUALIZE
            threshold_res_path = fullfile(path, "Thresholding_Comparision");
            if ~exist(threshold_res_path, 'dir')
                mkdir(threshold_res_path);
            end
            Ft = figure(1);

            num_rows = ceil((num_sensitivities + 3) / 3);
            subplot(num_rows, 3, 1); imshow(image_median./255); title('Median filtered image');
            subplot(num_rows, 3, 2); imshow(binarized_img_fuzzy); title('Fuzzy thresholding');
            subplot(num_rows, 3, 3); imshow(binarized_img_otsu); title('Otsu thresholding');
            for i = 1:num_sensitivities
                subplot(num_rows, 3, i + 3);
                imshow(binarized_images_adaptive{i});
                title(sprintf('Adaptive thresholding, s=%.1f', thresholding.sensitivity(i)));
            end
            
            % Save results as a matplot-figure
            save_path2 = fullfile(threshold_res_path, img_name);
            exportgraphics(gcf, strcat(save_path2, ".pdf"), 'ContentType', 'vector');
            saveas(gcf, strcat(save_path2, '.fig'));
            close(Ft);
        end

        % Set local adaptive thresholding with sensitivity = 0.3 as the 
        % default, which is used for all of the remaining computations
        adaptive_tr_def = adaptthresh(image_frangi, 0.3);
        binarized_img = imbinarize(image_frangi, adaptive_tr_def);

    otherwise
        fprintf(['The selected thresholding method <%s> does not exist. ' ...
            'Please select one from the following: {"fuzzy_thresholding", ' ...
            '"local adaptive thresholding", "otsu thresholding"}\n'], thresholding.method);
        disp("Use default method <Local adaptive thresholding> " + ...
            "with sensitivity = 0.3 instead");
        adaptive_tr_def = adaptthresh(image_frangi, 0.3);
        binarized_img = imbinarize(image_frangi, adaptive_tr_def);
end

if thresholding.opening_size > 0
    binarized_img = bwareaopen(binarized_img, thresholding.opening_size);
end

skeletonized_img = bwmorph(binarized_img, 'skel', Inf);

% Save skeleton to result folder
skeletonization_path = fullfile(path, "Skeletonization");
if ~exist(skeletonization_path, 'dir')
    mkdir(skeletonization_path);
end
skeletonized_img_path = fullfile(skeletonization_path, strcat(img_name, "_skeleton.tif"));
imwrite(skeletonized_img, skeletonized_img_path);

if VISUALIZE
    % Display the processed image at the different steps
    F1 = figure(1);
    subplot(1,4,1); imshow(image_median./255); title('1. Median filtering');
    subplot(1,4,2); imshow(image_frangi); title('2. Frangi filtering');
    if strcmp(thresholding.method, 'otsu_thresholding') || strcmp(thresholding.method, 'fuzzy_thresholding')
        threshold_title = char(strcat('3. ', char(strrep(thresholding.method, "_", " "))));
    elseif strcmp(thresholding.method, 'local_adaptive_thresholding')
        threshold_title = char(strcat("3. Adaptive thresholding, s=", num2str(thresholding.sensitivity(1))));
    else  % if "test_all" or a non-existing method was selected, use the default
        threshold_title = char("3. Adaptive thresholding, s=0.3");
    end
    subplot(1,4,3); imshow(binarized_img); title(threshold_title);
    subplot(1,4,4); imshow(skeletonized_img); title('4. Skeletonization');

    threshold_title_save = char(strrep(threshold_title(3:end), " ", "_"));
    save_path2 = strcat(save_path(1:end-4), "_", threshold_title_save);
    exportgraphics(gcf, strcat(save_path2, ".pdf"), 'ContentType', 'vector');
    saveas(gcf, strcat(save_path2, '.fig'));
    close(F1);
end
end
