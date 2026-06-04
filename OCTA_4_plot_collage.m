%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OCTA Workflow Pipeline 4: Comparing Morphology Across Images
%%%%
%%%% Version:     1.1
%%%% Date:        04/06/2026
%%%%
%%%% Authors:     Jesko Wagner (DBI-Infra IACF, jesko.wagner@sund.ku.dk)
%%%%
%%%% Description: This fourth section of the pipeline performs visual
%%%%              comparison of measurements performed in the previous step.
%%%%              Specifically, it allows plotting metrics measured on
%%%%              separate images with a unified color bar, so they can be
%%%%              compared qualitatively. The resulting figures may also be
%%%%              exported for further use in publications.
%%%%
%%%%              Comparison works across independent runs: the user selects
%%%%              one or more 'MorphologyResults_full.mat' files (possibly
%%%%              from different result folders), their measurements are
%%%%              combined, and the user then picks which images to compare.
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ================= Initialization =======================================
close all
clear all
addpath('helper_scripts');

%% ================= Select MorphologyResults_full.mat files ==============
% Let the user pick .mat files, one folder at a time, so that results from
% independent runs (which live in different result folders) can be combined.
% uigetfile returns files from a single folder per call, so we loop and keep
% accumulating until the user cancels.
matlist   = {};             % full paths to all selected .mat files
startPath = 'results';      % folder the next dialog opens in
while true
    [f, p] = uigetfile({'*.mat', 'Morphology results (*.mat)'}, ...
        'Select MorphologyResults_full.mat file(s) - Cancel when done', ...
        startPath, 'MultiSelect', 'on');
    if isequal(f, 0)
        break;              % user is done selecting
    end
    if ischar(f)
        f = {f};            % single file selected
    end
    matlist   = [matlist, fullfile(p, f)];
    startPath = p;          % reopen next dialog in the same folder
end

if isempty(matlist)
    error('No .mat files were selected. Exiting.');
end

%% ================= Combine results across selected files ================
% Read each .mat file and concatenate its fullMorphResults struct array. A
% 'source' field (the parent result folder name) is added to each entry so
% that images with identical filenames from different runs stay
% distinguishable in the selection dialog below.
allResults = [];
for m = 1:numel(matlist)
    loaded = load(matlist{m});
    if ~isfield(loaded, 'fullMorphResults') || isempty(loaded.fullMorphResults)
        warning('Skipping "%s": no fullMorphResults found.', matlist{m});
        continue;
    end
    these = loaded.fullMorphResults;
    [~, runName] = fileparts(fileparts(matlist{m}));   % result folder name
    [these.source] = deal(runName);

    if isempty(allResults)
        allResults = these;
    else
        these = orderfields(these, allResults);        % match field order
        allResults = [allResults, these];
    end
end

if isempty(allResults)
    error('None of the selected .mat files contained usable results. Exiting.');
end

%% ================= Select which images to compare =======================
% Present every image found across the selected files and let the user pick
% which ones to compare. Selection is by index, so duplicate filenames from
% different runs are handled unambiguously.
labels = arrayfun(@(s) sprintf('%s  |  %s', s.source, s.name), ...
    allResults, 'UniformOutput', false);
[sel, ok] = listdlg( ...
    'PromptString',   {'Select images to compare', '(Ctrl/Shift-click for multiple)'}, ...
    'ListString',     labels, ...
    'SelectionMode',  'multiple', ...
    'ListSize',       [450 300]);
if ~ok || isempty(sel)
    error('No images were selected. Exiting.');
end
data = allResults(sel);

% data.image options:
% binary = binary vessel segmentation
% skeleton = binary vessel skeleton segmentation
% skeleton_diameters = vessel skeleton colored by diameter
% labeled_skeleton = vessel skeleton labeled by which vessel branch it belongs to

options = ["Vessel diameter" ...
           "Mean branch diameter" "Branch length"];

% decide which metric to plot (listdlg centres on screen, unlike menu)
[chosen_metric, ok] = listdlg( ...
    'PromptString',  'Select what to plot', ...
    'ListString',    cellstr(options), ...
    'SelectionMode', 'single', ...
    'ListSize',      [220 90]);
if ~ok || isempty(chosen_metric)
    error('No metric was selected. Exiting.');
end

% prepare canvas to right dimensions
canvas = plot.prepare_canvas(numel(data));

% plot each image with correct metric overlay
for k = 1:numel(data)
    % activate axis
    axes(canvas(k));
    [~, plot_title] = fileparts(data(k).name);
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
cb.Layout.Tile = 'east';   % span the right of the whole layout, not just the last panel
cb.Label.String = options(chosen_metric) + " (µm)";
cb.Label.FontSize = 12;

% export image to a directory (defaults to the first selected file's folder)
parameters.output_dir = fileparts(matlist{1});
plot.save_figure(gcf, parameters);
