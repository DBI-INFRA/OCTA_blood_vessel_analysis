function [n, r] = BoxCount2D(inputSkel)
% -------------------------------------------------------------------------
% Description:
%    This function calculates fractal dimensions for a 2D input skeleton 
%    using the box-counting method. It returns two arrays which can be used
%    for calculating the fractal dimension D, from the relation n = r^D. 
%    The fractal dimension D can be determined from the slope of the 
%    straight line fitted to the log-log plot of log n = D * log r.
%
% Parameters:
%    - inputSkel: Skeletonized binary image of blood vessels
%
% Returns:
%    - n: The number of boxes for each width r
%    - r: The widths of the boxes
% -------------------------------------------------------------------------

stack_size = size(inputSkel);
box_power  = nextpow2(stack_size);
pad_target = 2.^nextpow2(stack_size);

% Pad skeleton image to size
start_idx  = floor((pad_target-stack_size)./2)+1;
stop_idx   = start_idx+stack_size-1;
PaddedSkel = zeros(pad_target);
PaddedSkel(start_idx(1):stop_idx(1), start_idx(2):stop_idx(2)) = inputSkel;

% Initialise box counting and iterate over box sizes
max_box   = min(box_power);
r = 2.^(0:max_box);
n = zeros(1, max_box+1);

pyr_img   = dlarray(PaddedSkel, 'SS');
n(1) = sum(pyr_img>0, 'all');
for bx = 2:max_box+1
    pyr_img = avgpool(pyr_img, [2, 2], "Stride", 2)*4; % sum-pool a 2x2 box accross the image
    n(bx) = sum(pyr_img>0, 'all');
end
