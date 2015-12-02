mincrmnms - implementation of RMNMS

This is a C implementation of RMNMS (Rotation Invariant Multi-modal
Non-local MS lesion segmentation) inspired by the original implem-
entation of BEaST by Simon F. Eskildsen (Simon F. ESkildsen, Pierrick
Coupé, Vladimir Fonov, José V. Manjón, Kelvin K. Leung, Nicolas 
Guizard, Shafik N. Wassef, Lasse R. Østergaard, D. Louis Collins, 
and The Alzheimer's Disease Neuroimaging Initiative, BEaST: Brain 
extraction based on nonlocal segmentation technique, NeuroImage, 
vol. 59(3), pp. 2362-2373. ISSN 1053-8119, 
10.1016/j.neuroimage.2011.09.012.


The latest RMNMS source code can be found at:
https://github.com/BIC-MNI/RMNMS

mincrmnms works with MINC1 and MINC2 images. However, experimental
support for the NIfTI format has been added. This has not yet been
tested thoroughly.

mincrmnms needs a library of priors to work (see below).

Compiling
---------
mincrmnms requires either MINC or NIfTI libraries. mincrmnms has been
tested on Debian type Linux systems, such as Ubuntu.
To configure type:


mkdir build && cd build

ccmake ../CMakeLists.txt

and set the right paths. Then

or

cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt/minc/ -DLIBMINC_DIR:PATH=/opt/minc/lib/ -DBUILD_TESTING:BOOL=OFF ../


make
make install

Troubleshooting:
- NIFTI_ROOT should be set to /usr if you installed NIfTI libraries
using the package libnifti-dev
- If the compiler cannot find hdf5.h you probably need to install
  libhdf5-serial-dev
- If you get the message: "Could not find module FindLIBMINC.cmake or
  a configuration file for package LIBMINC.", you must point to the
  directory containing either FindLIBMINC.cmake or
  LIBMINCConfig.cmake. If you have installed MINC Tool Kit,
  http://www.bic.mni.mcgill.ca/ServicesSoftware/ServicesSoftwareMincToolKit,
  the directory is most likely /opt/minc/lib

Library
-------
The library folder MUST contain these files:
library.masks.1mm
library.t2w.1mm
library.flr.1mm


mincrmnms will try to access these six files.

The library.{t2w,flr}.* files contain filenames of the normalized images. 
It is important that the filenames are in the same order across the
library files. mincrmnms uses the line number in the files to link
images at different resolutions, as well as linking the images to the
segmentations.
Similarly, the library.masks.* files contain filenames of the "expert"
segmentations.

In addition, mincrmnms will by default try to access two binary masks
named margin_mask.mnc and intersection_mask.mnc. These specify respectively
which voxels to process and which voxels are automatically included in the
output segmentation. These can be manually set using -mask and -positive
and disabled using -no_mask and -no_positive.

mincrmnms assumes that all images are in the same space and have the
same origin. This is not checked at runtime and will lead to errors if
it is not fulfilled.

mincrmnms uses a simple intensity based comparison metric. Thus, it is
very important that the intensities of the library images have been
normalized.

See README.library for how to install existing libraries.

Usage
-----
mincrmnms [options] <library dir> <input t2w> <input flr> <output> 

<library>: path to the library
<input>: input image
<output>: output segmentation


Examples
--------
Suppose your library resides in ~/rmnms/ and your normalized input
image is named t2w.mnc and flr.mnc, this command will provide a 
segmentation:

mincrmnms ~/rmnms/lib/ t2w.mnc flr.mnc output.mnc -conf conf/default.1mm.conf -fill -mask mask.mnc

Explanation of options
----------------------
 -probability:      Output the probability map instead of crisp mask.
 -flip:             Flip images around the mid-sagittal plane to increase patch count.
 -load_moments:     Do not calculate moments instead use precalculated library moments. (for optimization purposes)
 -fill:             Fill holes in the binary output.
 		    Just in case we get errors inside the mask.
 -median:           Apply a median filter on the probability map.
 		    Makes the segmentation slightly more robust, but may limit the accuracy.
 -nlm_filter:       Apply an NLM filter on the probability map (experimental).
 -verbose:          Enable verbose output.
 -clobber:          Clobber output files
 -abspath:          File paths in the library are absolute (default is relative to library root).
 -voxel_size:       Specify voxel size for calculations (4, 2, or 1). Assumes no multiscale. Use configuration file for multiscale.
                Default value: 4
 -patch_size:       Specify patch size for single scale approach.
                Default value: 1
 -search_area:      Specify size of search area for single scale approach.
                Default value: 2
 -alpha:            Specify confidence level Alpha.
                Default value: 0.5
 -beta:             Specify smoothness factor Beta.
                Default value: 0.25
 -threshold:        Specify threshold for patch selection.
                Default value: 0.95
 -selection_num:    Specify number of selected images.
                Default value: 20
 -positive:         Specify mask of positive segmentation (inside mask) instead of the default mask.
 		    This will be added to the final segmentation
 -output_selection: Specify file to output selected files.
 -count:            Specify file to output the patch count.
 -configuration:    Specify configuration file.
 		    See the 'conf' folder for example configurations.
 -mask:             Specify a segmentation mask instead of the the default mask.
 -same_resolution:  Output final mask with the same resolution as input file.
 -no_mask:          Do not apply a segmentation mask. Perform the segmentation over the entire image.
 -no_positive:      Do not apply a positive mask.


Tuning the parameters
---------------------
The main parameters are those in the configuration files. You should
always set them in the configuration file instead of at the command
line. For example, the configuration file (default.1mm.conf) looks
like this:

# voxelsize patchsize searcharea alpha beta threshold num_selected
1 1 5 0.5 0.25 0.95 40

The header shows the parameter names. "voxelsize" is on of 1, 2, or 4
indicating the scale step. RMNMS only suport native scale. The
propagation is controlled by "alpha", which determine how much
information is propagated. Here alpha is not 0.2 for the 4mm scale, which
means that probabilities in the range 0.2 - 0.8 are propagated, while
<0.2 are considered background and >0.8 are considered foreground
(object). If you increase "alpha", you trust your lowres segmentation
more and propagate less to the next scale. For the final scale "alpha"
is usually 0.5, because this is the threshold for the final
segmentation. You can adjust this final threshold to control
consistent over/under-segmentation.

"patchsize" is the size of the patch when comparing structures across
the library. 1 means a patch size of 3x3x3, 2 means a patch size of
5x5x5, and so on. Increasing this may give better results, but also
seriously increases the computational time.

"searcharea" is the size of the spatial neighborhood in which to look
for similar patches. This is similar to the "patchsize" in that e.g.
2 ~ 5x5x5 and 4 ~ 9x9x9.

"beta" is a smoothness parameter. Usually in the range 0-1. Larger
beta means more smooth.

"threshold" is a parameter controlling the preselection of
patches. Preselection makes sure to only include patches that have
some similarity. The range is 0-1. The higher the more strict (fewer
patches selected).

"num_selected" is the number of images to select from the library when
looking for similar structures. Usually the higher the
better. However, the improvement is asymptotic and the memory usage
quickly rises with larger N.

For the command line options, the important ones are:

-same_res which simply makes sure that the output has the same
 resolution as the input no matter the configuration file.

-median applies a median filter on the probability map before
 propagation. This also smoothes the result and should be used when
 the library is not perfect (i.e. from another population than the
 image to segment).

-fill can always be used as this just morphologically fills any holes
 in the segmentation.

The remaining command line options are not really relevant for
tuning. However, you may want to output the probability map (using
-probability) when determining the best "alpha".

Populating the library
----------------------
The best way to improve your results is to populate the library with
images/masks from the same scanner as the the images you are trying to
segment. One way to do this is to run RMNMS with the default ICBM/ADNI
images and then select the best masks among the results, possibly
perform some manual corrections, and put them into the library. Then
run RMNMS again. This bootstrapping method can be performed
iteratively with increasing performance improvements.

For more on populating the library, please see README.library

Reference
---------
Please cite RMNMS as:

Nicolas Guizard, Pierrick Coupé, Vladimir Fonov, Jose V. Manjon, 
Douglas L. Arnold, D. Louis Collins, “Rotation-invariant multi-
contrast non-local means for MS lesion segmentation”, Neuroimage: 
Clinical, 2015, doi:10.1016/j.nicl.2015.05.001 

Contact
-------
For questions and feedback, please contact
Nicolas GUIZARD <n.guizard@gmail.com>
