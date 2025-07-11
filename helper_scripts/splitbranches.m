function [labeled_skeleton, num_segments] = splitbranches(skeletonized_img)
% -------------------------------------------------------------------------
% Description:
%    This function splits the blood vessels into labelled segments between 
%    branching points, and returns the number of independent 
%    segments in a skeletonized image.
%
% Parameters:
%    - skeletonized_img: Skeletonized binary image of blood vessels
%
% Returns:
%    - labeled_skeleton: Labeled skeleton image, with 0 as background, 1 as
%                        branch points and labels >=2 as vessel segments
%    - num_segments: Number of blood vessel segments
% -------------------------------------------------------------------------

% Compute branching points and endpoints
branch_points = bwmorph(skeletonized_img, 'branchpoints');

% Define structuring element in order to set the neighboring pixels of
% branch-points to zero.
cross_se = strel([0 1 0; 1 1 1; 0 1 0]);

% Dilate the branch points to include their 4-connected neighbors
dilated_branch_points = imdilate(branch_points, cross_se);

% Subtract the dilated branch points from the skeletonized image (set to 0)
skeleton_segments = skeletonized_img & ~dilated_branch_points;

% Give each connected component a separate label
[labeled_segments, num_segments] = bwlabel(skeleton_segments);

% Reconnect branch ends that were removed for CC labeling
branch_ends = dilated_branch_points & ~branch_points & skeletonized_img;
full_branch = labeled_segments + branch_ends;
labeled_branch = ordfilt2(full_branch, 9, ones(3,3)) .* logical(full_branch);
labeled_skeleton = imclearborder(labeled_branch + skeletonized_img, 4);

end

