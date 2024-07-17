function [avg_length, num_segments] = vessel_length(skeletonized_img)
% -------------------------------------------------------------------------
% Description:
%    This function computes the average length of blood vessel segments
%    between all branching points, as well as the number of independent 
%    segments in a skeletonized image.
%
% Parameters:
%    - skeletonized_img: Skeletonized binary image of blood vessels
%
% Returns:
%    - avg_length: Average length of blood vessel segments
%    - num_segments: Number of blood vessel segments
% -------------------------------------------------------------------------

% Compute branching points and endpoints
branchpoints = bwmorph(skeletonized_img, 'branchpoints');
%fprintf('Branchpoints: [%s]\n', join(string(branchpoints), ','));

% Define structuring element in order to set the neighboring pixels of
% branch-points to zero.
cross_se = strel([0 1 0; 1 1 1; 0 1 0]);

% Dilate the branch points to include their 4-connected neighbors
dilated_branch_points = imdilate(branchpoints, cross_se);

% Subtract the dilated branch points from the skeletonized image (set to 0)
skeleton_segments = skeletonized_img & ~dilated_branch_points;

% Give each connected component a separate label
[segments_labeled, num_segments] = bwlabel(skeleton_segments);
%disp(num_segments);
%fprintf('Labeled Segments: [%s]\n', join(string(segments_labeled), ','));

% Compute & print the length of all blood vessel segments
lengths = regionprops(segments_labeled, 'Area');
%disp([lengths(:).Area]);

% Compute & return the average blood vessel length
avg_length = mean([lengths.Area]);

end

