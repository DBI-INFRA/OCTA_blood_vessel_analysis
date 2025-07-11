function [AlignedStack, ZDisplacement] = zalignstack(stack, method)
% -------------------------------------------------------------------------
% Description:
%    This function aligns a 3D stack along the z-dimension using a 
%    displacement field determined by finding the peak intensity slice and 
%    applying a smoothing method. 
%
% Parameters:
%    - stack: Input 3D stack to be aligned
%    - method: Method for z-alignment ('median' or 'plane'). Options:
%       - Median: Applies median filtering to smooth z-displacement
%       - Plane: Fits a plane to the z-positions and calculates the 
%                z-displacement based on the fitted plane equation
%
% Returns:
%    - AlignedStack: Aligned 3D stack
%    - ZDisplacement: Displacement field used for alignment
% -------------------------------------------------------------------------

stack_size = size(stack);

% Find peak in z-value
[~, zpeak_idx] = max(stack, [], 3, 'linear');
[~, ~, k] = ind2sub(stack_size, zpeak_idx);
k = filloutliers(k, "center", "percentiles", [0, 25], 2);

% Calculate smooth alignement
switch method
    case 'median'
    % Applies median filtering to create smooth z_displacement 
    ZDisplacement = int16(medfilt2(k, [25 25], 'symmetric'));

    case 'plane'
    % Calculated equation of plane based on image center
    [xx, yy] = meshgrid(1:stack_size(2), 1:stack_size(1));
    zv = ones(size(k));
    X = [xx(:), yy(:), zv(:)];
    C = X\k(:);
    fittedPlane = round(reshape(X*C, stack_size(1), stack_size(2)));
    ZDisplacement = fittedPlane;

end

disp_max = max(ZDisplacement, [], 'all');
d_field = zeros([stack_size, 3]);
d_field(:,:,:,3) = repmat(ZDisplacement, 1, 1, stack_size(3));
AlignedStack = double(imwarp(stack, d_field));
AlignedStack = AlignedStack(:,:,1:end-disp_max);



