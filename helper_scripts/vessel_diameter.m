function avg_diameter = vessel_diameter(binarized_img, skeletonized_img)
% --------------------------------------------------------------------------
% Description:
%    This function computes the average vessel diameter from a binary vessel 
%    image and its skeletonized version. It utilizes the Euclidean distance  
%    transform to measure distances from the skeleton to vessel boundaries.
%
% Parameters:
%    - binarized_img: Segmented binary image of blood vessels
%    - skeletonized_img: Skeletonized version of the segmented vessels
%
% Returns:
%    - avg_diameter: The average vessel diameter
% --------------------------------------------------------------------------

% Compute the Euclidean distance transform of the binary vessel image
distance_transform = bwdist(~binarized_img);

% Sample the distances along the skeletonized vessels
skeleton_distances = distance_transform(skeletonized_img);

% Calculate the average of these distances (i.e., the average radius) &
% multiply by 2 to get the average diameter.
average_radius = mean(skeleton_distances);
avg_diameter = 2 * average_radius;

% Toggle on visualization to inspect the results
visualize = false;

if visualize
    figure;
    
    % Plot the original binary image
    subplot(2, 2, 1);
    imshow(binarized_img);
    title('Binary Vessel Image');
    
    % Plot the skeletonized image
    subplot(2, 2, 2);
    imshow(skeletonized_img);
    title('Skeletonized Image');
    
    % Plot the distance transform
    subplot(2, 2, 3);
    imagesc(distance_transform);
    axis image;
    colorbar;
    hold on;
    [row, col] = find(skeletonized_img);
    plot(col, row, 'r.');
    title('Skeleton Overlaid on Distance Transform');
    hold off;
    
    % Plot the skeleton distances as a histogram
    subplot(2, 2, 4);
    histogram(skeleton_distances);
    title('Histogram of Skeleton Distances');
    xlabel('Distance (pixels)');
    ylabel('Frequency');
end
end
