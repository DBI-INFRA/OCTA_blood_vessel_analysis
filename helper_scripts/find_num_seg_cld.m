
%%%% Description: An optional script to print number of vessel segments at
%%%% CLD for all the images (so that you do not need to go to the csv file 
% for each image)

% ================= User Parameters =======================================
% Input image directory containing the cropped OCT images

result_dir = 'results/240730_170157';


% ================= Parse input and preallocate storage ===================
root_path = fullfile(result_dir, "2_CLD_SPD_estimation");

files = dir(fullfile(root_path, '*_data.csv'));
for k = 1:length(files)
    if files(k).name(1) ~= '.'
        cur_path = fullfile(root_path, files(k).name);
        a = readtable(cur_path);
        max_num_elements = round(max(a.smoothed_num_segments));
        
        file_name_parts = strsplit(files(k).name, '_');
        file_name = strcat(file_name_parts{1}, '_', file_name_parts{2});
        
        fprintf('max_num_elements for %s = %d\n', file_name, max_num_elements);
    end
end
