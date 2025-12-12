%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 4: Comparing Morphology Across Images
%%%%
%%%% Version:     1.0
%%%% Date:        09/12/2025
%%%%
%%%% Authors:     Jesko Wagner (DBI-Infra IACF, jesko.wagner@sund.ku.dk)
%%%%
%%%% Description: This fourth section of the pipeline performs viusual  
%%%%              comparison of measurements performed in the previous step.
%%%%              Specifically, it allows plotting metrics measured on
%%%%              separate images with a unified color bar, so they can be
%%%%              compared qualitatively. The resulting figures may also be
%%%%              exported for further use in publications.
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ================= Initialization =======================================
close all
clear all
addpath('helper_scripts');

%%% Specify result_dir (path where the desired 'Parameters.mat' is stored)
result_dir = 'results\yymmdd_hhmmss\';

%%% Load exisiting user configurations
parameters = load(fullfile(result_dir, 'Parameters.mat'));

%% ================ Parse Inputs and Allocate Storage =====================
% Update modified variables
existingVars = fieldnames(parameters);
for i = 1:numel(existingVars)
    eval([existingVars{i} ' = parameters.' existingVars{i} ';']);
end
save(fullfile(result_dir, 'Parameters.mat'), existingVars{:});


try
    % Load image information and find any pre-processed image
    imgInfo = readtable(fullfile(result_dir, 'ImageSummary.csv'), 'ReadRowNames', true, 'Delimiter', ',');
    if exist(fullfile(result_dir, 'ProcessedImages'), 'dir')
        image_dir = fullfile(result_dir, 'ProcessedImages');
    end

    % Contains fullMorphResults
    load(fullfile(result_dir, 'MorphologyResults_full.mat'));

catch ME
    fprintf("Check configuration file: ")
    rethrow(ME);
end

% File order depends on the OS
[files, path] = uigetfile({'*.tif;*.tiff', 'TIFF Files (*.tif, *.tiff)'}, ...
    'Select TIFF files', result_dir, 'MultiSelect', 'on');


if isequal(files, 0)
    error('No file was selected. Exiting.');
elseif isa(files, 'char') 
    files = cell({files}); % single image has been passed
end

% Subset morphology results by selected files
data = plot.get_struct_by_name(fullMorphResults, files);

% data.image options:
% binary = binary vessel segmentation
% skeleton = binary vessel skeleton segmentation
% skeleton_diameters = vessel skeleton colored by diameter
% labeled_skeleton = vessel skeleton labeled by which vessel branch it belongs to

options = ["Vessel diameter" ...
           "Mean branch diameter" "Branch length"];

           % decide which metric to plot
chosen_metric = menu('Select what to plot', options);

% prepare canvas to right dimensions
canvas = plot.prepare_canvas(length(files));

% plot each image with correct metric overlay
for k = 1:length(files)
    % activate axis
    axes(canvas(k));
    [~, plot_title] = fileparts(files(k));
    plot_data = data(k);
    if chosen_metric == 1
        plot_data = data(k).image.skeleton_diameters;
    elseif chosen_metric == 2
        plot_data = labelmapper(data(k).image.labeled_skeleton, [data(k).object.Label], [data(k).object.MeanDiameter_um]);
    elseif chosen_metric == 3
        plot_data = labelmapper(data(k).image.labeled_skeleton, [data(k).object.Label], [data(k).object.Length_um]);
    end
    plot.plot_panel(plot_data, plot_title);
end

cmap = [0, 0, 0; turbo(255)];
colormap(cmap);

cb = colorbar;
cb.Label.String = options(chosen_metric) + " (µm)";
cb.Label.FontSize = 12;
% export image to a directory
plot.save_figure(gcf, parameters)