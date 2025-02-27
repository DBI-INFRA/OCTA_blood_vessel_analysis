# OCTA_blood_vessel_analysis
Multi-component pipeline for automated analysis of blood vessel morphology imaged using angiographic optical coherence tomography

## MATLAB Toolboxes required
- Curve Fitting Toolbox
- Deep Learning Toolbox
- Fuzzy Logic Toolbox
- Signal Processing Toolbox
- Wavelet Toolbox

## Usage and code structures

This package consists of a "config.m" file, and three "OCTA_... .m" scripts in a main folder, and a subfolder "helper_scripts" which includes MATLAB functions that are called by the three main scripts. 
The "config.m" file should be run before any of the main scripts are run. The three main scripts are expected to be run in order for analysis of blood vessel morphology at the Superficial Plexus Depth (SPD), but can be operated independently.

- **Expected Input**: An entire folder of OCTA images. We provide an additional folder of pre-cropped sample images in the folder `sample_data` as an example of expected input. Please keep `AutoCrop` set to `false` for these test images, since they have already been cropped for size.
- **Expected Output**: Each time "config.m" is run, a subfolder with the current date-time in the form of `yyMMdd_HHmmss` (e.g., "240722_154054") will be created in a chosen result folder "result_dir" along with three files "ImageSummary.csv", "Parameters.m", "Parameters.csv". The pipeline will save all intermediate and final outputs into the subfolder and will also expect to receive intermediate data from the subfolder.

Detailed usage of the main scripts are given below.

**The working directory should be set to the main folder where the three "OCTA_".m files are located.**

___
### config.m
This script allows the user to define parameters to be used in the three main scripts. 

**You MUST run this script, before running any of the other scripts.**

#### User Input Parameters
<u>Image Information and Directories:</u>
- **`image_dir`** : Path to folder where images to be processed are stored. (Default = "sample_data")
- **`output_dir`**: Path to folder where all results should be stored. A time-stamped subfolder will be automatically created here each time you run this script.
- **`pixel_size`**: \[x, y, z] size of each pixel in micrometers (Sample data has size \[6.49, 6.49, 2.45])

<u>Image Pre-processing Parameters (used in OCTA_1 only):</u>
- **`preprocess_image`**: true/false toggle. Set to false to ignore all image preprocessing parameters in this section. OCTA_1 will not work if this is set to false.
- **`AutoCrop`**: true/false toggle. Set to true to automatically crop noisy z-slices away from the dataset (see point 1 above). Set to false for testing the default ```sample_data``` folder.
- **`CropSensitivity`**: Value between \[0, 1] or 'None', default = 0.5. Controls where the automatic cropping of the image stops, set to 1 to crop away data more agressively and reduce file size; set to 0 for less agressive removal of slices, which may result in larger files and slower processing times later on in the pipeline. Can also be set to 'None' so that AutoCrop will keep data all the way to the end of the image.
- **`Wavelet Transform Parameters`**: Parameters used to remove motion artefacts in the wavelet-FFT filter. Choose direction of stripe artefaces to remove and how many times to apply filter.
- **`AutoAlign`**: true/false toggle. Set to true to automatically align image along air-skin border in z-direction (see point 3 above).
- **`writeZDisplacement`**: true/false toggle. When set to true, the pipeline will write the z-displacement field as a tiff file, which can be used for debugging issues with skin alignment.

<u>Vessel Enhancement Parameters (used in OCTA_2 AND 3):</u>
- **`median_filter_size`** (default = 5): The size of the median filter used to smooth the image before further processing. Each pixel in the output image will have the median value of its neighborhood (defined by `size x size`) in the input image. If your image has a lot of noise, you might want to increase this value to get smoother results. However, too much smoothing can cause loss of detail, e.g. small blood vessels may be smoothed out and not visible anymore. See the [Matlab documentation for the median filter](https://de.mathworks.com/help/images/ref/medfilt2.html) for more details.
- **`frangi_opts`**: These parameters control the Frangi filter, which enhances vessel-like structures in the image. Adjusting these values can fine-tune the vessel detection process. If the default parameters don't work well, try adjusting the sigmarange first, this parameter has the biggest influence on the resulting filtered image. To understand how the filter and parameters work in detail, see [Frangi et al. (1998)](https://doi.org/10.1007/BFb0056195).
  - **`.sigmarange`** (default = [1 6]): Controls the range of scales at which vessels will be detected, the first number defining the lower bound (smallest vessels) and the second the upper bound (biggest vessels). If the vessels in your image are very thin, reduce the lower bound (e.g., `sigmarange = [0.5 6]`). If your vessels are thick, increase the upper bound (e.g., `sigmarange = [1 8]`). Generally, if the vessels in the image vary a lot in size, the range should be large (e.g. `sigmarange = [0.5 8]`), while if they are all e.g. very small, both the lower and upper bound should be rather small (e.g. `sigmarange = [0.5 3]`).
  - **`.sigmastepsize`** (default = 1): The step size between the scales used in the Frangi filter. A smaller step size can provide more precise detection, but it increases computation time. A step size of `1` for a range of `[1 6]` will detect vessels at scale 1, 2, 3, 4, 5 and 6, a step size of `0.5` will compute vessels at scale 1, 1.5, 2, 2.5 etc. and a step size of `2` will skip over every second vessel scale, reducing the computation time.
  - **`.correctionconst1`** (default = 0.8): A constant used for correcting vesselness calculation. If the vessels in your images are under-detected, try reducing this value (e.g. `correctionconst1 = 0.6`).
  - **`.correctionconst2`** (default = 15): Another constant used for correcting vesselness calculation. If larger vessels aren't well detected, try increasing this value (e.g. `correctionconst2 = 20`).

<u>Vessel Quantification Parameters (used in OCTA_3 only):</u>
Note: Thresholding parameters are set to optimised defaults for OCTA_2 and cannot be changed from config file.
- **`thresholding`**: The following two parameters control the image segmentation step after filtering. Both the chosen method and sensitivity have a significant influence on how the image is separated into vessels and background. See for instance [Peter Bankhead's e-book chapter on thresholding](https://bioimagebook.github.io/chapters/2-processing/3-thresholding/thresholding.html) for an introduction to thresholding and different thresholding methods. 
  - **`.method`**: The thresholding method applied to segment the image into foreground (vessels) and background. The options available are `local_adaptive_thresholding` (the default method, uses[MATLAB's adaptthresh function](https://se.mathworks.com/help/images/ref/adaptthresh.html)), `otsu_thresholding` (uses [MATLAB's graythresh function](https://se.mathworks.com/help/images/ref/graythresh.html)) and `fuzzy_thresholding` (*helper_scripts/fuzzy_thresholding.m* from the [OCTAVA toolbox](https://github.com/GUntracht/OCTAVA/tree/main))
  - **`.sensitivity`**: The sensitivity parameter for the local adaptive thresholding algorithm. Higher values will pick up more details like smaller vessels, but potentially also more noise. The sensitivity is only relevant if local_adaptive_thresholding is selected under `.method`.
  - **`.opening_size`**: Morphological opening removes some small islands after segmentation, based on area size (in pixels)
- **`quantify_frame`**: Central frame for quantifying morphology. The value can be either SINGLE FRAME NUMBER, or a VARIABLE NAME (Column Name) in 'ImageSummary.csv'.
- **`quantify_range_um`**: Depth above and below the selected frame to include in quantification (microns)                             
- **`save_skeleton`**: Save skeletonisation image to file.

#### Expected Outputs (stored in "`results_dir`/`yyMMdd_HHmmss`/")
- **"ImageSummary.csv"**: List of all images from input_dir folder.
- **"Parameters.mat"**: .mat file storing above user input variables.

---
### OCTA_1_Preprocessing.m
The first part of the OCTA analysis pipeline performs image preprocessing. For each image in the `input_dir` folder, this script will:
1. (If `AutoCrop == true`) Crop away slices from the start and the end of the z-stack where the image signal is faint. This is done by taking the mean intensity of each z-plane and calculating an expected background noise from the mean intensity of the last 10 z-slices. Cropping begins when a local minima is first detected in the mean z-intensity signal, and stops when the mean z-intensity returns to a level close to the background noise (determined by the `CropSensitivity` parameter).
2. Remove motion artefacts slice-by-slice using a wavelet-FFT filter, as described in [Byers et al. (2017)](https://doi.org/10.1364%2FBOE.8.004551).
3. Detect the air-skin border by first finding the peak in the image intensity for each pixel along its z-dimension, then removing outliers greater than the 25th percentile, and finally by applying a 25x25 median filter. The resulting 3D surface was used to axially shift the OCT data in the z-dimension and align the entire image at the air-skin boundary.

#### User Input Parameters
<u>Required:</u>
- **`result_dir`**: You will need to modify the value of `result_dir` based on the subfolder generated by "config.m". e.g., Change from `result_path = 'results/240722_154054';` to `result_path = 'results/240722_172132';`.

<u>Optional:</u>
- Use the Image Pre-processing Parameters section in config.m to set pre-processing behaviour 
- You can overwrite existing config.m parameters in this script using "\<ParameterName> = \<newValue>;"

#### Expected Outputs (stored in "`output_dir`/`yyMMdd_HHmmss`/")
- **"ProcessedImages/"**: Folder of all processed images which have been cropped, filtered and aligned. Images written as tiff.
- **"ProcessedImages/Displacement/"**: Folder of z-displacement fields (if writeZDisplacement == true), written as tiff.
- **"ImageSummary.csv"**: Updated with FirstSlice and LastSlice specifying the z-frames used from the original image data after cropping.

___
### OCTA_2_CLD_SPD_estimation.m
The second part of the pipeline estimates the Capillary Loop Depth (CLD) and Superficial Plexus Depth (SPD) by performing blood vessel segmentation and skeletonisation for each z-slice in the image, and counting the number of independent skeletons (i.e., the number of unconnected blood vessel networks) through the depth of the 3D-image. The CLD is identified at the depth slice position with the maximum number of unconnected blood vessel networks, and the SPD at the depth slice position where the gradient of the curve of independent vessel networks reaches 0 for the first time after the detected CLD-depth.
The script will preferentially use preprocessed images from "`results_dir`/`yyMMdd_HHmmss`/ProcessedImages/", if it exists, or else will default to images in `input_dir`

#### User Input Parameters
<u>Required:</u>
- **`result_dir`**: You will need to modify the value of `result_dir` based on the subfolder generated by "config.m". e.g., Change from `result_path = 'results/240722_154054';` to `result_path = 'results/240722_172132';`.

<u>Optional:</u>
- Use the Vessel Enhancement Parameters section in config.m to set vessel processing behaviour 
- You can overwrite any existing parameters in 'config.m' using "\<ParameterName> = \<newValue>;"
- Note: Segmentation parameters for OCTA_2 are set to optimised defaults and cannot be changed from config file.

#### Expected Outputs (stored in "`results_dir`/`yyMMdd_HHmmss`/")
- **"CLD_SPD_estimation/"**: Folder of graphs and corresponding raw data excel sheets listing the number of independent blood vessel networks found at each image depth, with the identified CLD and SPD locations marked on each graph with a red circle.
- **"ImageSummary.csv"**: Updated with the automatically calculated `CLD_Frame`, `SPD_Frame`, `CLD_Depth_um` and `SPD_Depth_um` appended to the table for each processed image. To manually change the SPD depth that will be used in the third script "OCTA_3_MorphQuant.m", please edit the value listed under the `SPD_Frame` column in this file. The corresponding `SPD_Depth_um` will also have to be manually corrected.
- **"SPD_not_found.txt"**: Text file that will log all images where the SPD depth could not be automatically calculated. In these cases, you can manually select the SPD depth by examining the graphs found in "CLD_SPD_estimation/" and editing "ImageSummary.csv" as specified above.
___
### OCTA_3_MorphQuant.m
The third section of the pipeline performs quantitative analysis of a mean-intensity-projected blood vessel network around a user specified depth-range. It calculates various morphological metrics for each image, including the mean vessel diameter, vessel length, vessel density and fractal dimension (using the box counting method). All results are saved in a CSV file. 
The depth can be specified as a specific frame or as a Variable_Name in 'ImageSummary.csv'. 

#### User Input Parameters
<u>Required:</u>
- **`result_dir`**: You will need to modify the value of `result_dir` based on the subfolder generated by "config.m". e.g., Change from `result_path = 'results/240722_154054';` to `result_path = 'results/240722_172132';`.

<u>Optional:</u>
- Use the Vessel Quantification Parameters section in config.m to set vessel quantification behaviour 
- You can overwrite any existing parameters in 'config.m' using "\<ParameterName> = \<newValue>;"

#### Expected Outputs (stored in "`results_dir`/`yyMMdd_HHmmss/`")
- **"MIP_skeletonization_results/"**: Folder of the visualization of skeletonisation and thresholding (segmentation) results written as ".fig" and ".pdf" files.
- **"MorphologyResults.csv"**: A table showing the morphology results of the skeletonized blood vessel network at the SPD-depth of each image, including mean diameter, length, vessel density, fractal dimension.
- **"Parameters.csv"**: All user input variables saved in "Parameters.mat" are written out into a .csv file for ease of reading.
- **"WarningLog.txt"**: Text file to log all images where the full SPD range (`[SPD_depth-round(SPD_range_um/pixel_size(3)), SPD_depth+round(SPD_range_um/pixel_size(3))]`) cannot be read in, either because the SPD was not found (SPD_frame = "NaN") or because the specified SPD range goes beyond the slices available in the corresponding "1_AlignedImage/" file. In these cases, the user may choose to manually specify an SPD (see expected outputs for "OCTA_2_CLD_SPD_estimation.m") or re-crop the images with "OCTA_1_Preprocessing" with a lower `CropSensitivity` setting, for these specific images with errors.
