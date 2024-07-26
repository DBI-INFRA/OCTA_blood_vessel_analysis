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

For each aligned image in `results_dir`/`yyMMdd_HHmmss/1_AlignedImages/`, the script will perform blood vessel segmentation and skeletonisation for each z-slice in the image, and count the number of independent skeletons (i.e., the number of unconnected blood vessel networks) through the depth of the 3D-image. The CLD is identified at the slice position with the maximum number of unconnected blood vessel networks, and the SPD at the slice position where the gradient of the curve of independent vessel networks reaches 0 for the first time after the CLD-depth position.

#### User Input Parameters
- **`result_path`**: You will need to go to modify the value of `result_path` based on the subfolder generated by the first script. e.g., Change from `result_path = 'results/240722_154054';` to `result_path = 'results/240722_172132';`.
- **`median_filter_size`**: Size of median filter to apply to image (higher values = more smoothing);
- **`frangi_opts`**: Parameters for the Frangi filter for vessel enhancement (see [Frangi et al. (1998)](https://doi.org/10.1007/BFb0056195))
- **`thresholding`**: Method for image segmentation after filtering step. See usage guide for "OCTA_3_MorphQuant.m" below to view more options for testing the available thresholding methods.
  - **`.method`**: Choose between "local_adaptive_thresholding" (default, MATLAB adaptthresh), "otsu_thresholding" (MATLAB graythresh) and "fuzzy_thresholding" ("helper_scripts/fuzzy_thresholding.m" from the [OCTAVA toolbox](https://github.com/GUntracht/OCTAVA/tree/main))
  - **`.sensitivity`**: Sensitivity parameter for local adaptive thresholding algorithm (higher values will pick up more and smaller vessels, but potentially also more noise). Only relevant if local_adaptive_thresholding is selected.
  
#### Expected Outputs (stored in "`results_dir`/`yyMMdd_HHmmss`/")
- **"2_CLD_SPD_estimation/"**: Folder of graphs and corresponding raw data listing the number of independent blood vessel networks found at each image depth, with the identified CLD and SPD locations marked on each graph with a red circle.
- **"ImageSummary.csv"**: Updated ImageSummary file from "OCTA_1", now with the automatically calculated CLD_Frame, SPD_Frame, CLD_Depth_um and SPD_Depth_um appended for each image processed. To manually change the the SPD used for "OCTA_3_MorphQuant.m", please edit the the value listed under the "SPD_Frame" column in this file. The corresponding SPD_Depth_um will also have to be manually corrected.
- **"InputParameters.mat"**: .mat file with all new user input variables appended to the existing file, including "median_filter_size", "frangi_opts", and "thresholding").
- **"SPD_not_found.txt"**: Text file to log all images where the SPD_depth could not be calculated. User may then manually select the SPD_depth by examining the graphs in "2_CLD_SPD_estimation/" and editing "ImageSummary.csv" as specified above.

___
### OCTA_3_MorphQuant.m
The third section of the pipeline performs quantitative analysis of the blood vessel network at the depth of the SPD (the script will include signal from the frames below and above the SPD to get more robust results. This range of frames around the SPD can be specified by the user). It calculates various morphological metrics for each image, including the mean vessel diameter, vessel length, vessel density and fractal dimension (using the box counting method). All results are saved in a CSV file.

#### User Input Parameters
- **`result_path`**: You will need to go to modify the value of `result_path` based on the subfolder generated by the first script. e.g., Change from `result_path = 'results/240722_154054';` to `result_path = 'results/240722_172132';`.
- **`thresholding_sensitivity`**: Sensitivity parameter for local adaptive thresholding algorithm (higher values will pick up more and smaller vessels, but potentially also more noise). Only relevant if local_adaptive_thresholding is selected. This value is set different between OCTA_2 and OCTA_3.
- **`test_thresholding_methods`**: true/false toggle. Set to true to test and view the following different thresholding methods --- "fuzzy_thresholding", "otsu_thresholding", and "local_adaptive_thresholding" with the user's specified sensitivities (See below). If enabled, the script will continue to calculate morphological metrics with the default thresholding method ("local_adaptive_thresholding", sensitivity = 0.3). If not enabled, the script will use the thresholding parameters specified in "OCTA_2_CLD_SPD_estimation.m" and saved in "InputParameters.mat".
  - **`local_adaptive_thresholding_sensitivities`**: User-specified local_adaptive_thresholding sensitivities to test, given as an array. The script will generate a subplot of the local adaptive thresholding result for each different given sensitivity parameter in this array.
- **`save_skeleton`**: true/false toggle to write calculated image skeleton to file.
- **`SPD_range_um`**: Depth around the calculated/specified SPD_Depth listed in "ImageSummary.csv" to consider when quantifying vessel morphology, in micrometers. Morphology quantification will then be performed on the mean intensity projection of slices `[SPD_depth-round(SPD_range_um/pixel_size(3)), SPD_depth+round(SPD_range_um/pixel_size(3))]` (SPD_depth $\pm$ SPD_range_um/pixel_size(3)), where `pixel_size(3)` is the voxel size in the z-dimension in microns.

#### Expected Outputs (stored in "`results_dir`/`yyMMdd_HHmmss/`")
- **"3_MIP_skeletonization_results/"**: Folder of the visualization of skeletonisation and thresholding (segmentation) results written as ".fig" and ".pdf" files.
- **"MorphologyResults.csv"**: A table showing the morphology results of the skeletonized blood vessel network at the SPD-depth of each image, including mean diameter, length, vessel density, fractal dimension.
- **"InputParameters.csv"**: All user input variables saved in "InputParameters.mat" are written out into a .csv file for ease of reading.

- **"WarningLog.txt"**: Text file to log all images where the full SPD range (`[SPD_depth-round(SPD_range_um/pixel_size(3)), SPD_depth+round(SPD_range_um/pixel_size(3))]`) cannot be read in, either because the SPD was not found (SPD_frame = "NaN") or because the specified SPD range goes beyond the slices available in the corresponding "1_AlignedImage/" file. In these cases, the user may choose to manually specify an SPD (see expected outputs for "OCTA_2_CLD_SPD_estimation.m") or re-crop the images with "OCTA_1_Preprocessing" with a lower `CropSensitivity` setting, for these specific images with errors.
