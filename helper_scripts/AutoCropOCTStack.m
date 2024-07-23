function [stack_out, crop_range] = AutoCropOCTStack(stack_in, sensitivity)
% -------------------------------------------------------------------------
% Description:
%    This function automatically crops an OCT stack to the region where the 
%    signal is strongest. It identifies start and stop positions based on 
%    signal prominence and local minima.
%
% Parameters:
%    - stack_in: Input 3D OCT stack
%    - sensitivity: Value between [0, 1] to control amount of cropping in
%    deeper parts of OCT stack. 0 means less data will be cropped from end 
%    of image and 1 means more data will be cropped from the end. Default
%    = 0.5. Can be set to 'None' to prevent cropping from end of dataset
%
% Returns:
%    - stack_out: Cropped 3D OCT stack
%    - crop_range: Range of slices retained [start_pos, stop_pos]
% -------------------------------------------------------------------------

stack_mean = squeeze(mean(stack_in, [1 2]));
start_pos  = find(islocalmin(stack_mean,'MinProminence',10), 1);

if ~exist("sensitivity", "var")
    sensitivity = 0.5;    
end

if sensitivity == "None"
    stop_pos = size(stack_in, 3);
elseif (0 <= sensitivity) && (sensitivity <= 1)
    stop_range = sort([stack_mean(start_pos), median(stack_mean(end-9:end))]);
    stop_val   = sensitivity*diff(stop_range)+stop_range(1);
    stop_pos   = find(stack_mean(start_pos:end) < stop_val, 1) + start_pos;
else
    error("Sensitivity must either be a value between 0 and 1, or 'None'.")
end

crop_range = [start_pos, stop_pos];
stack_out = stack_in(:,:,start_pos:stop_pos);
end