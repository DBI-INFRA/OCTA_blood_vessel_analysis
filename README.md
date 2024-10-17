# OCTA_blood_vessel_analysis
Multi-component pipeline for automated analysis of blood vessel morphology imaged using angiographic optical coherence tomography

## MATLAB Toolboxes required
- Fuzzy Logic Toolbox
- Curve Fitting Toolbox
- Wavelet Toolbox
- Deep Learning Toolbox
- Signal Processing Toolbox

## Usage and code structures

This package consists of three "OCTA_... .m" scripts in a main folder that needs to be run sequentially, and a subfolder "helper_scripts" which includes MATLAB functions that are called by these three main scripts. 

- **Expected Input**: An entire folder of OCTA images. We provide an additional folder of pre-cropped sample images in the folder `sample_data` as an example of expected input. Please keep `AutoCrop` set to `false` for these test images, since they have already been cropped for size.
- **Expected Output**: Each time "OCTA_1_Preprocessing.m" is run, a subfolder with the current date-time in the form of `yyMMdd_HHmmss` (e.g., "240722_154054") will be created in a chosen result folder `results_dir`. The pipeline will save all intermediate and final outputs into the subfolder and will also expect to receive intermediate data from the subfolder.

Detailed usage of the main scripts are given below.

**The working directory should be set to the main folder where the three "OCTA_".m files are located.**

---
### OCTA_1_Preprocessing.m

The first part of the OCTA analysis pipeline performs image preprocessing. For each image in the `input_dir` folder, this script will:
1. (If `AutoCrop == true`) Crop away slices from the start and the end of the z-stack where the image signal is faint. This is done by taking the mean intensity of each z-plane and calculating an expected background noise from the mean intensity of the last 10 z-slices. Cropping begins when a local minima is first detected in the mean z-intensity signal, and stops when the mean z-intensity returns to a level close to the background noise (determined by the `CropSensitivity` parameter).
2. Remove motion artefacts slice-by-slice using a wavelet-FFT filter, as described in [Byers et al. (2017)](https://doi.org/10.1364%2FBOE.8.004551).
3. Detect the air-skin border by first finding the peak in the image intensity for each pixel along its z-dimension, then removing outliers greater than the 25th percentile, and finally by applying a 25x25 median filter. The resulting 3D surface was used to axially shift the OCT data in the z-dimension and align the entire image at the air-skin boundary.

#### User Input Parameters
- **`input_dir`**: Path to folder where images to be processed are stored. (Default = "sample_data")
- **`result_dir`**: Path to folder where all results should be stored. A result subfolder will be automatically created here each time you run this script.
- **`pixel_size`**: \[x, y, z] size of each pixel in micrometers
- **`writeZDisplacement`**: true/false toggle. When set to true, the pipeline will write the z-displacement field (see point 3 above) as a tiff file, which can be used for debugging issues with skin alignment.
- **`AutoCrop`**: true/false toggle. Set to true to automatically crop noisy z-slices away from the dataset (see point 1 above). Set to false for testing the default ```sample_data``` folder.
- **`CropSensitivity`**: Value between \[0, 1] or 'None', default = 0.5. Controls where the automatic cropping of the image stops, set to 1 to crop away data more agressively and reduce file size; set to 0 for less agressive removal of slices, which may result in larger files and slower processing times later on in the pipeline. Can also be set to 'None' so that AutoCrop will keep data all the way to the end of the image.
- **`Wavelet Transform Parameters`**: Parameters used to remove motion artefacts in the wavelet-FFT filter. (See helper_scripts/SuppWaveletFFT.m for details.)

#### Expected Outputs (stored in "`results_dir`/`yyMMdd_HHmmss`/")
- **"1_AlignedImages/"**: Folder of all processed images which have been cropped, filtered and aligned. Images written as tiff.
- **"1_AlignedImages/Displacement/"**: Folder of z-displacement fields (if writeZDisplacement == true), written as tiff.
- **"ImageSummary.csv"**: List of all images from input_dir folder, with FirstSlice and LastSlice specifying the z-frames used from the original image data after cropping.
- **"InputParameters.mat"**: .mat file storing user input variables, including "pixel_size" and all preprocessing options ("AutoCrop", "CropSensitivity" and all Wavelet Transform Parameters).


___
### OCTA_2_CLD_SPD_estimation.m
The second part of the pipeline estimates the Capillary Loop Depth (CLD) and Superficial Plexus Depth (SPD) from the z-aligned Optical Coherence Tomography (OCT) images produced in "OCTA_1_Preprocessing.m". 

For each aligned image in `results_dir`/`yyMMdd_HHmmss/1_AlignedImages/`, the script will perform blood vessel segmentation and skeletonisation for each z-slice in the image, and count the number of independent skeletons (i.e., the number of unconnected blood vessel networks) through the depth of the 3D-image. The CLD is identified at the depth slice position with the maximum number of unconnected blood vessel networks, and the SPD at the depth slice position where the gradient of the curve of independent vessel networks reaches 0 for the first time after the detected CLD-depth.

#### User Input Parameters
<u>Parameters that **must** be adjusted:</u>
- **`result_path`**: This parameter needs to be modified to match the subfolder generated by the first script. Each time you run the first script, a new subfolder will be created, and the name of this subfolder needs to be inserted here. E.g. if the first script generated a subfolder named `240722_172132`, update the line as follows: `result_path = 'results/240722_172132';`.

<u>Parameters that can be adjusted if needed (**optional**):</u>
- **`median_filter_size`** (default = 5): The size of the median filter used to smooth the image before further processing. Each pixel in the output image will have the median value of its neighborhood (defined by `size x size`) in the input image. If your image has a lot of noise, you might want to increase this value to get smoother results. However, too much smoothing can cause loss of detail, e.g. small blood vessels may be smoothed out and not visible anymore. See the [Matlab documentation for the median filter](https://de.mathworks.com/help/images/ref/medfilt2.html) for more details.
- **`frangi_opts`**: These parameters control the Frangi filter, which enhances vessel-like structures in the image. Adjusting these values can fine-tune the vessel detection process. If the default parameters don't work well, try adjusting the sigmarange first, this parameter has the biggest influence on the resulting filtered image. To understand how the filter and parameters work in detail, see [Frangi et al. (1998)](https://doi.org/10.1007/BFb0056195).
  - **`.sigmarange`** (default = [1 6]): Controls the range of scales at which vessels will be detected, the first number defining the lower bound (smallest vessels) and the second the upper bound (biggest vessels). If the vessels in your image are very thin, reduce the lower bound (e.g., `sigmarange = [0.5 6]`). If your vessels are thick, increase the upper bound (e.g., `sigmarange = [1 8]`). Generally, if the vessels in the image vary a lot in size, the range should be large (e.g. `sigmarange = [0.5 8]`), while if they are all e.g. very small, both the lower and upper bound should be rather small (e.g. `sigmarange = [0.5 3]`).
  - **`.sigmastepsize`** (default = 1): The step size between the scales used in the Frangi filter. A smaller step size can provide more precise detection, but it increases computation time. A step size of `1` for a range of `[1 6]` will detect vessels at scale 1, 2, 3, 4, 5 and 6, a step size of `0.5` will compute vessels at scale 1, 1.5, 2, 2.5 etc. and a step size of `2` will skip over every second vessel scale, reducing the computation time.
  - **`.correctionconst1`** (default = 0.8): A constant used for correcting vesselness calculation. If the vessels in your images are under-detected, try reducing this value (e.g. `correctionconst1 = 0.6`).
  - **`.correctionconst2`** (default = 15): Another constant used for correcting vesselness calculation. If larger vessels aren't well detected, try increasing this value (e.g. `correctionconst2 = 20`).

<u>Parameters that are not recommended to adjust:</u>

The following parameters were adjustable in previous versions of this script, though we have optimized them for the detection of CLD & SPD and recommend not to change them, unless the resulting graphs and estimated CLD + SPD depths look really off.

- **`thresholding`**: The following two parameters control the image segmentation step after filtering. Both the chosen method and sensitivity have a significant influence on how the image is separated into vessels and background. See for instance [Peter Bankhead's e-book chapter on thresholding](https://bioimagebook.github.io/chapters/2-processing/3-thresholding/thresholding.html) for an introduction to thresholding and different thresholding methods. 
  - **`.method`**: The thresholding method applied to segment the image into foreground (vessels) and background. The options available are `local_adaptive_thresholding` (the default method, uses[MATLAB's adaptthresh function](https://se.mathworks.com/help/images/ref/adaptthresh.html)), `otsu_thresholding` (uses [MATLAB's graythresh function](https://se.mathworks.com/help/images/ref/graythresh.html)) and `fuzzy_thresholding` (*helper_scripts/fuzzy_thresholding.m* from the [OCTAVA toolbox](https://github.com/GUntracht/OCTAVA/tree/main))
  - **`.sensitivity`**: The sensitivity parameter for the local adaptive thresholding algorithm. Higher values will pick up more details like smaller vessels, but potentially also more noise. The sensitivity is only relevant if local_adaptive_thresholding is selected under `.method`.
  
#### Expected Outputs (stored in "`results_dir`/`yyMMdd_HHmmss`/")
- **"2_CLD_SPD_estimation/"**: Folder of graphs and corresponding raw data excel sheets listing the number of independent blood vessel networks found at each image depth, with the identified CLD and SPD locations marked on each graph with a red circle.
- **"ImageSummary.csv"**: Updated ImageSummary file from "OCTA_1", now with the automatically calculated `CLD_Frame`, `SPD_Frame`, `CLD_Depth_um` and `SPD_Depth_um` appended to the table for each processed image. To manually change the SPD depth that will be used in the third script "OCTA_3_MorphQuant.m", please edit the value listed under the `SPD_Frame` column in this file. The corresponding `SPD_Depth_um` will also have to be manually corrected.
- **"InputParameters.mat"**: .mat file with all new user input parameters appended to the existing user parameters file, including the `median_filter_size`, `frangi_opts`, `thresholding.method` and `thresholding.sensitivity`.
- **"SPD_not_found.txt"**: Text file that will log all images where the SPD depth could not be automatically calculated. In these cases, you can manually select the SPD depth by examining the graphs found in "2_CLD_SPD_estimation/" and editing "ImageSummary.csv" as specified above.

___
### OCTA_3_MorphQuant.m
The third section of the pipeline performs quantitative analysis of the blood vessel network at the depth of the SPD (the script will include signal from the frames below and above the SPD to get more robust results. This range of frames around the SPD can be specified by the user). It calculates various morphological metrics for each image, including the mean vessel diameter, vessel length, vessel density and fractal dimension (using the box counting method). All results are saved in a CSV file.

#### User Input Parameters
- **`result_path`**: You will need to go to modify the value of `result_path` based on the subfolder generated by the first script. e.g., Change from `result_path = 'results/240722_154054';` to `result_path = 'results/240722_172132';`.
- **`SPD_frame`**: Manually specified SPD frame number. Default value is `-1` indicating that the SPD_frame will be inferred from output of previous step. Set it to a positive integer value to manually specify the SPD frame number.
- **`thresholding_sensitivity`**: Sensitivity parameter for local adaptive thresholding algorithm (higher values will pick up more and smaller vessels, but potentially also more noise). Only relevant if local_adaptive_thresholding is selected. This value is set different between OCTA_2 and OCTA_3.
- **`test_thresholding_methods`**: true/false toggle. Set to true to test and view the following different thresholding methods --- "fuzzy_thresholding", "otsu_thresholding", and "local_adaptive_thresholding" with the user's specified sensitivities (See below). If enabled, the script will continue to calculate morphological metrics with the default thresholding method ("local_adaptive_thresholding", sensitivity = 0.3). If not enabled, the script will use the thresholding parameters specified in "OCTA_2_CLD_SPD_estimation.m" and saved in "InputParameters.mat".
  - **`local_adaptive_thresholding_sensitivities`**: User-specified local_adaptive_thresholding sensitivities to test, given as an array. The script will generate a subplot of the local adaptive thresholding result for each different given sensitivity parameter in this array.
- **`save_skeleton`**: true/false toggle to write calculated image skeleton to file.
- **`SPD_range_um`**: Depth around the calculated/specified SPD_Depth listed in "ImageSummary.csv" to consider when quantifying vessel morphology, in micrometers. Default value is 30 micrometers, following [Byers et al. (2017)](https://doi.org/10.1364%2FBOE.8.004551). Morphology quantification will then be performed on the mean intensity projection of slices `[SPD_depth-round(SPD_range_um/pixel_size(3)), SPD_depth+round(SPD_range_um/pixel_size(3))]` (SPD_depth $\pm$ SPD_range_um/pixel_size(3)), where `pixel_size(3)` is the voxel size in the z-dimension in microns.

#### Expected Outputs (stored in "`results_dir`/`yyMMdd_HHmmss/`")
- **"3_MIP_skeletonization_results/"**: Folder of the visualization of skeletonisation and thresholding (segmentation) results written as ".fig" and ".pdf" files.
- **"MorphologyResults.csv"**: A table showing the morphology results of the skeletonized blood vessel network at the SPD-depth of each image, including mean diameter, length, vessel density, fractal dimension.
- **"InputParameters.csv"**: All user input variables saved in "InputParameters.mat" are written out into a .csv file for ease of reading.

- **"WarningLog.txt"**: Text file to log all images where the full SPD range (`[SPD_depth-round(SPD_range_um/pixel_size(3)), SPD_depth+round(SPD_range_um/pixel_size(3))]`) cannot be read in, either because the SPD was not found (SPD_frame = "NaN") or because the specified SPD range goes beyond the slices available in the corresponding "1_AlignedImage/" file. In these cases, the user may choose to manually specify an SPD (see expected outputs for "OCTA_2_CLD_SPD_estimation.m") or re-crop the images with "OCTA_1_Preprocessing" with a lower `CropSensitivity` setting, for these specific images with errors.
