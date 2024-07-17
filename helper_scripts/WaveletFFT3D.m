function [image] = WaveletFFT3D(image, gamma, order, wname, orient, reps, dim)
% ---------------------------------------------------------------------------
% Description:
%    This function extends SuppWaveletFFT.m for vertical/horizontal wavelet 
%    filtering to 3D. It applies the wavelet filtering to the selected
%    dimension "dim" of the 3D-volume.
%
% Parameters:
%    - image: Input 3D image
%    - gamma: Parameter controlling the filter strength
%    - order: The order of the wavelet decomposition
%    - wname: Wavelet name used for decomposition
%    - orient: Orientation for filtering ('h' = horizontal, 'v' = vertical, 'b' = both)
%    - reps: Number of iterations for wavelet filtering
%    - dim: Dimension to process directionally (1, 2, or 3)
%
% Returns:
%    - image: Filtered 3D image after wavelet reconstruction
% ---------------------------------------------------------------------------

% 1) Parse input
% Convert image into single
image = squeeze(single(image)./255);
n = ndims(image);
image_size = size(image);

% Parse image dimension for directional processing
if n > 3
    error('Image must have at most 3 dimensions');
elseif n == 2
    image_size = [image_size, 1];
elseif n ~= dim
    dimorder = 1:3; 
    dimorder(dim) = []; dimorder = [dimorder, dim];
    image = permute(image, dimorder);
    image_size = size(image);
end

% Parse orientation
switch orient(1)
    case 'h'
        O = 1;
    case 'v'
        O = 2;
    case 'b'
        O = 3;
    otherwise
        error('Unknown orientation input for WaveletFFT3D');
end

% 2) Perform wavelet transform iteratively on image
for ss = 1:image_size(end)
    slice = image(:,:,ss);
    for rr = 1:reps
        imageFiltered = SuppWaveletFFT(slice,gamma,order,wname,O);
        slice = min(imageFiltered,slice); % Take the minimum between the filtered and old image
    end
    image(:,:,ss) = slice;
end

if exist("dimorder", "var")
    image = ipermute(image, dimorder);
end