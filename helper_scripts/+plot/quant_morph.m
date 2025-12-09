function [out, readErrorIdx] = quant_morph(ff, parameters, imgInfo, readErrorIdx)

% 1) Parse image information
img_path = fullfile(parameters.image_dir, parameters.filelist(ff).name);

avg_measurements = struct('meanDiameter', NaN, 'meanLength', NaN, 'meanDensity', NaN, 'fracDimension', NaN);

% Check validity of frame at which to quantify
if isscalar(parameters.quantify_frame)
    central_frame = parameters.quantify_frame;
else
    central_frame  = imgInfo{parameters.filelist(ff).name, parameters.quantify_frame};
end

if isnan(central_frame)
    readErrorIdx{ff, 1} = true;
    readErrorIdx{ff, 2} = NaN;
    return;
end

img_reso        = parameters.pixel_size(1:2);
quantify_range  = round(parameters.quantify_range_um/parameters.pixel_size(3));
quantify_slices = [central_frame-quantify_range, central_frame+quantify_range];

% Read image around chosen depth
try
    ImageStack = tiffreadVolume(img_path, 'PixelRegion', {[1 inf], [1 inf], quantify_slices});
catch
    ImageStack = tiffreadVolume(img_path, 'PixelRegion', {[1 inf], [1 inf]});
    
    quantify_slices(2) = size(ImageStack, 3);
    readErrorIdx{ff, 1} = true; 
    readErrorIdx{ff, 2} = [quantify_slices, diff(quantify_slices)*parameters.pixel_size(3)];
end

if numel(size(ImageStack)) == 4
    ImageStack = ImageStack(:,:,:,1);
end

% Create mean-intensity-projection and skeletonize image
StackMIP = mean(ImageStack, 3);
if max(StackMIP(:)) <= 1
    StackMIP = StackMIP * 255;
end
StackMIP = double(int16(squeeze(StackMIP)));
[skeleton, binary_img] = skeletonization(StackMIP, parameters.median_filter_size, parameters.frangi_opts, parameters.thresholding, false, "");


% 2) Quantify the blood vessel network morphology
% Split skeleton into branches
[labeled_skeleton, ~] = splitbranches(skeleton);

% Calculate vessel diameter at each point along skeleton
skeleton_diameters = double(bwdist(~binary_img).*skeleton*2);

% Make Vessel Measurements per Branch
vessel_morph = regionprops(labeled_skeleton, skeleton_diameters, 'MeanIntensity', 'Area');
vessel_label = regionprops(labeled_skeleton, labeled_skeleton, 'MeanIntensity');
per_vessel_measurements = table([vessel_label.MeanIntensity]', ...
                            [vessel_morph.Area]', [vessel_morph.Area]'*img_reso(1), ...
                            [vessel_morph.MeanIntensity]', [vessel_morph.MeanIntensity]' *img_reso(1), ...
                            'VariableNames', {'Label', 'Length_px', 'Length_um', 'MeanDiameter_px', 'MeanDiameter_um'}); 
per_vessel_measurements = per_vessel_measurements([per_vessel_measurements.Label]>1,:);
num_vessels = size(per_vessel_measurements, 1);

% Average Vessel Branch Length in um
avg_measurements.meanLength = mean([per_vessel_measurements.Length_um]);

% Average Vessel Diameter in um (including at branch points)
avg_diameter = mean(skeleton_diameters(skeleton));
avg_measurements.meanDiameter = avg_diameter*img_reso(1);

% Average Vessel Density per mm2
img_reso_mm = img_reso*(10^-3);
img_area = prod(size(binary_img).*img_reso_mm);
avg_measurements.meanDensity = num_vessels/img_area;

% Fractal Dimension
[n, r] = boxcount2D(skeleton);
avg_measurements.fracDimension = -1*fit(log(r)', log(n)', 'poly1').p1;

out = struct(); % Initialize the output struct
out.avg = avg_measurements;
out.object = per_vessel_measurements;
out.image = struct();
out.image.binary = binary_img;
out.image.skeleton = skeleton;
out.image.labeled_skeleton = uint16(labeled_skeleton);
out.image.skeleton_diameters = skeleton_diameters;
