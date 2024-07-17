function [stack_out, crop_range] = AutoCropOCTStack(stack_in)
% -------------------------------------------------------------------------
% Description:
%    This function automatically crops an OCT stack to the region where the 
%    signal is strongest. It identifies start and stop positions based on 
%    signal prominence and local minima.
%
% Parameters:
%    - stack_in: Input 3D OCT stack
%
% Returns:
%    - stack_out: Cropped 3D OCT stack
%    - crop_range: Range of slices retained [start_pos, stop_pos]
% -------------------------------------------------------------------------

stack_mean = squeeze(mean(stack_in, [1 2]));

start_pos = find(islocalmin(stack_mean,'MinProminence',10), 1);
stop_val = mean([stack_mean(start_pos), median(stack_mean(end-9:end))]);
stop_pos = find(stack_mean(start_pos:end) < stop_val, 1) + start_pos;

crop_range = [start_pos, stop_pos];
stack_out = stack_in(:,:,start_pos:stop_pos);
end