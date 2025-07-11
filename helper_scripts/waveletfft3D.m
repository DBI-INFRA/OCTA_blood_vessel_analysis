function image = waveletfft3D(image, gamma, order, wname, orient, reps, dim)
% ---------------------------------------------------------------------------
% Description:
%    This function applies the subfunction SuppWaveletFFT for vertical/
%    horizontal wavelet filtering in 3D. It applies the wavelet filtering 
%    to the selected dimension "dim" of the 3D-volume.
%
% Parameters:
%    - image: Input 3D image
%    - gamma: Parameter controlling the filter strength
%    - order: The order of the wavelet decomposition
%    - wname: Wavelet name used for decomposition
%    - orient: Orientation for filtering ('h' = horizontal, 'v' = vertical, 'b' = both)
%    - reps: Number of iterations for wavelet filtering
%    - dim: Dimension to process directionally (x, y, or z)
%
% Returns:
%    - image: Filtered 3D image after wavelet reconstruction
% ---------------------------------------------------------------------------

% 1) Parse input
% Convert image into single
image = squeeze(single(image)./255);
n = ndims(image);
image_size = size(image);

% Parse direction to apply filters
switch dim
    case 'y'
        D = 1;
    case 'x'
        D = 2;
    case 'z'
        D = 3;
    otherwise
        error('Unknown orientation input for WaveletFFT3D');
end

% Parse image dimension for directional processing
if n > 3
    error('Image must have at most 3 dimensions');
elseif n == 2
    image_size = [image_size, 1];
elseif n ~= D
    dimorder = 1:3; 
    dimorder(D) = []; dimorder = [dimorder, D];
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
        imageFiltered = SuppWaveletFFT(slice, gamma, order, wname, O);
        slice = min(imageFiltered,slice); % Take the minimum between the filtered and old image
    end
    image(:,:,ss) = slice;
end

if exist("dimorder", "var")
    image = ipermute(image, dimorder);
end
end 


function image = SuppWaveletFFT(image,gamma,order,wname,O)
% ---------------------------------------------------------------------------
% Description:
%    This function performs vertical / horizontal wavelet filtering of the 
%    input image. It applies wavelet decomposition followed by FFT-based 
%    suppression in the specified orientation(s).
%
% Parameters:
%    - image: Input 2D image
%    - gamma: Parameter controlling the filter strength
%    - order: The order of the wavelet decomposition
%    - wname: Wavelet name used for decomposition
%    - O: Orientation for filtering (1 = horizontal, 2 = vertical, 3 = both)
%
%    Returns:
%    - image: Filtered image after wavelet reconstruction.
% ---------------------------------------------------------------------------

detail = single(image);
x1 = size(detail,2); y1 = size(detail,1);

% Perform wavelet decomposition
for n = 1:order
    [detail,horiz{n},verti{n},diag{n}] = dwt2(detail,wname);
end

% Perform FFT on the vertical and/or horizontal frequency bands
for n = 1:order
    
    if O == 1 || O == 3
    fhoriz = fftshift(fft(horiz{n},[],2),2);
    [y,x] = size(fhoriz);
    hmask = repmat((1-exp(-[-floor(x/2):-floor(x/2)+x-1].^2/gamma)),y,1);
    fhoriz = fhoriz .* hmask;
    horiz{n} = ifft(ifftshift(fhoriz,2),[],2);
    end
    
    if O == 2 || O == 3
    fverti = fftshift(fft(verti{n},[],1),1);
    [y,x] = size(fverti);
    vmask = repmat((1-exp(-[-floor(y/2):-floor(y/2)+y-1].^2/gamma))',1,x);
    fverti = fverti .* vmask;
    verti{n} = ifft(ifftshift(fverti,1),[],1);
    end

end

% Reconstruct wavelets.
for n = order:-1:1 % Count backwards
    detail = detail(1:size(horiz{n},1),1:size(horiz{n},2));
    detail = idwt2(detail,horiz{n},verti{n},diag{n},wname);
end

image = detail(1:y1,1:x1);

end