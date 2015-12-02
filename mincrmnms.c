/*  rmnms.c
 *
 *  Copyright 2015  Nicolas Guizard Simon Fristed Eskildsen, Vladimir Fonov,
 *                Pierrick Coupé, Jose V. Manjon
 *
 *  This file is part of rmnms.
 *
 *  rmnms is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  rmnms is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with rmnms.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  For questions and feedback, please contact:
 *  Nicolas Guizard <n.guizard@gmail.com>
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif //HAVE_CONFIG_H

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#include "ParseArgv.h"
#include "array_alloc.h"
#include "nlmseg.h"
#include "rmnms.h"
#include "label.h"

const char LICENSE[]="Copyright (C) 2015\tNicolas Guizard, Simon Fristed Eskildsen, Vladimir Fonov, \n\
\t\t\tPierrick Coupe, Jose V. Manjon\n\n\
This program comes with ABSOLUTELY NO WARRANTY; for details type 'cat COPYING'. \n\
This is free software, and you are welcome to redistribute it under certain\n\
conditions; type 'cat COPYING' for details.\n\
";


int main(int argc, char  *argv[] )
{
  char *input_file1,*input_file2,*output_file,*libdir, history_label[1024];
  char imagelist1[FILENAMELENGTH],imagelist2[FILENAMELENGTH], masklist[FILENAMELENGTH],meanlist1[FILENAMELENGTH],meanlist2[FILENAMELENGTH],varlist1[FILENAMELENGTH],varlist2[FILENAMELENGTH];
  char ***images1, ***images2, ***masks, ***means1, ***means2, ***vars1, ***vars2;
  int num_images,i,sizes[3][5],tmpsizes[5],volumesize,*selection,steps=3,filled=0;
  float *imagedata1,*imagedata2,*maskdata,*meandata1,*meandata2,*vardata1, *vardata2, **subject1, **subject2,**mask,**positivemask=NULL,**segsubject,**patchcount,**filtered;
  float max,min;
  float **segmented;
  float *tempdata;
  int scale,scaledvolumesize,scales[3] = {1,2,4};
  int masksize=0,initialscale,targetscale,scalesteps;
  rmnms_conf input_conf[3],configuration[3];
  image_metadata **meta;
  image_metadata *mask_meta;
  int targetvoxelsize=1;

  VIO_BOOL outputprob = FALSE;
  VIO_BOOL flipimages = FALSE;
  VIO_BOOL load_moments = FALSE;
  VIO_BOOL fill_output = FALSE;
  VIO_BOOL verbose = FALSE;
  VIO_BOOL medianfilter = FALSE;
  VIO_BOOL patchfilter = FALSE;
  VIO_BOOL abspath = FALSE;
  VIO_BOOL same_res = FALSE;
  VIO_BOOL clobber  = FALSE;
  VIO_BOOL nomask = FALSE;
  VIO_BOOL nopositive  = FALSE;

  int voxelsize=4;
  int sizepatch = 1;
  int searcharea = 2;
  double alpha = 0.5;
  double beta = 0.25;
  double threshold = 0.95;
  int selectionsize = 20;
  time_t timer;

  char *positive_file=NULL;
  char *selection_file=NULL;
  char *count_file=NULL;
  char *conf_file=NULL;
  char *mask_file=NULL;

  char *default_rmnms_library=RMNMS_LIBRARY_PREFIX;
  char *default_rmnms_mask=RMNMS_LIBRARY_PREFIX"/margin_mask.mnc";
  char *default_rmnms_positive_file=RMNMS_LIBRARY_PREFIX"/intersection_mask.mnc";
  char *default_rmnms_config=RMNMS_LIBRARY_PREFIX"/default.2mm.conf";

  /* Argument table */
  ArgvInfo argTable[] = {
    {
      "-probability", ARGV_CONSTANT, (char *) TRUE, (char *) &outputprob,
      "Output the probability map instead of crisp mask."
    },
    {
      "-flip", ARGV_CONSTANT, (char *) TRUE, (char *) &flipimages,
      "Flip images around the mid-sagittal plane to increase patch count."
    },
    {
      "-load_moments", ARGV_CONSTANT, (char *) TRUE, (char *) &load_moments,
      "Do not calculate moments instead use precalculated library moments. (for optimization purposes)"
    },
    {
      "-fill", ARGV_CONSTANT, (char *) TRUE, (char *) &fill_output,
      "Fill holes in the binary output."
    },
    {
      "-median", ARGV_CONSTANT, (char *) TRUE, (char *) &medianfilter,
      "Apply a median filter on the probability map."
    },
    {
      "-nlm_filter", ARGV_CONSTANT, (char *) TRUE, (char *) &patchfilter,
      "Apply an NLM filter on the probability map (experimental)."
    },
    {
      "-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
      "Enable verbose output."
    },
    {
      "-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
      "Clobber output files"
    },
    {
      "-abspath", ARGV_CONSTANT, (char *) TRUE, (char *) &abspath,
      "File paths in the library are absolute (default is relative to library root)."
    },

    {
      "-voxel_size", ARGV_INT, (char *) 1, (char *) &voxelsize,
      "Specify voxel size for calculations (4, 2, or 1). Assumes no multiscale. Use configuration file for multiscale."
    },
    {
      "-patch_size", ARGV_INT, (char *) 1, (char *) &sizepatch,
      "Specify patch size for single scale approach."
    },
    {
      "-search_area", ARGV_INT, (char *) 1, (char *) &searcharea,
      "Specify size of search area for single scale approach."
    },
    {
      "-alpha", ARGV_FLOAT, (char *) 1, (char *) &alpha,
      "Specify confidence level Alpha."
    },
    {
      "-beta", ARGV_FLOAT, (char *) 1, (char *) &beta,
      "Specify smoothness factor Beta."
    },
    {
      "-threshold", ARGV_FLOAT, (char *) 1, (char *) &threshold,
      "Specify threshold for patch selection."
    },
    {
      "-selection_num", ARGV_INT, (char *) 1, (char *) &selectionsize,
      "Specify number of selected images."
    },

    {
      "-positive", ARGV_STRING, (char *) 1, (char *) &positive_file,
      "Specify mask of positive segmentation (inside mask) instead of the default mask."
    },
    {
      "-output_selection", ARGV_STRING, (char *) 1, (char *) &selection_file,
      "Specify file to output selected files."
    },
    {
      "-count", ARGV_STRING, (char *) 1, (char *) &count_file,
      "Specify file to output the patch count."
    },
    {
      "-configuration", ARGV_STRING, (char *) 1, (char *) &conf_file,
      "Specify configuration file."
    },
    {
      "-mask", ARGV_STRING, (char *) 1, (char *) &mask_file,
      "Specify a segmentation mask instead of the the default mask."
    },
    {
      "-same_resolution", ARGV_CONSTANT, (char *) TRUE, (char *) &same_res,
      "Output final mask with the same resolution as input file."
    },
    {
      "-no_mask", ARGV_CONSTANT, (char *) TRUE, (char *) &nomask,
      "Do not apply a segmentation mask. Perform the segmentation over the entire image."
    },
    {
      "-no_positive", ARGV_CONSTANT, (char *) TRUE, (char *) &nopositive,
      "Do not apply a positive mask."
    },

    {NULL, ARGV_END, NULL, NULL, NULL}
  };
  
  fprintf(stderr,"\nrmnms --\t\tan implementation of RMNMS (Rotation-invariant Multi-modal non-local MS) \n\t\t\tusing non-local Segmentation Technique) version %s\n\n",PACKAGE_VERSION);

  /* Get the time, overwriting newline */
  timer = time(NULL);

  /* make minc-type stamp*/
  sprintf(history_label,"%s>>>%s",ctime(&timer),argv[0]);
  for (i=1; i<argc; i++)
    sprintf(history_label,"%s %s",history_label,argv[i]);

  /* Get arguments */
  if (ParseArgv(&argc, argv, argTable, 0) || (argc < 5)) {
    fprintf(stderr,LICENSE);
    fprintf(stderr,
            "\nUsage: %s [options] <library dir> <input T1W> <input FLAIR> <output>\n",
            argv[0]);
    fprintf(stderr,"       %s -help\n\n", argv[0]);

    exit(STATUS_ERR);
  }

  /* if(argc>3) */
  /* { */
  libdir = argv[argc-4];
  input_file1 = argv[argc-3];
  input_file2 = argv[argc-2];
  output_file = argv[argc-1];
  /* } else { */
  /*   libdir = default_rmnms_library; */
  /*   input_file = argv[argc-2]; */
  /*   output_file = argv[argc-1]; */

  /*   positive_file=default_rmnms_positive_file; */
  /*   conf_file=default_rmnms_config; */
  /*   medianfilter=TRUE; */
  /*   fill_output=TRUE; */
  /*   same_res=TRUE; */
  /*   fprintf(stderr,"WARNING: Running rmnms with default parameters:\n -median -configuration %s -positive %s -fill -same_resolution\n",conf_file,positive_file); */
  /* } */

  if (mask_file==NULL) {
    mask_file=malloc((strlen(libdir)+20)*sizeof(*mask_file));
    sprintf(mask_file,"%s/margin_mask.mnc",libdir);
  }
  if ((!nopositive) && (positive_file==NULL)) {
    positive_file=malloc((strlen(libdir)+30)*sizeof(*positive_file));
    sprintf(positive_file,"%s/intersection_mask.mnc",libdir);
  }

  if(!clobber) {
    if(!access(output_file,F_OK)) {
      fprintf(stderr,"ERROR! File exists: %s , run with -clobber\n",output_file);
      return STATUS_ERR;
    }
    if(count_file && !access(count_file,F_OK)) {
      fprintf(stderr,"ERROR! File exists: %s , run with -clobber\n",count_file);
      return STATUS_ERR;
    }
  }

  if ((voxelsize>4) || (voxelsize<1) || (voxelsize==3)) {
    fprintf(stderr,"ERROR! Initial voxel size must be either 4, 2, or 1\n");
    return STATUS_ERR;
  }

  meta = (image_metadata **)malloc(3*sizeof(image_metadata*));

  meta[0] = read_volume(input_file1, &tempdata, sizes[0]);
  if (meta[0] == NULL) {
    fprintf(stderr,"ERROR! Image not read: %s\n",input_file1);
    return STATUS_ERR;
  }
  volumesize=sizes[0][0]*sizes[0][1]*sizes[0][2];

  subject1 = alloc_2d_float(3,volumesize*sizeof(**subject1));
  cp_volume(tempdata, subject1[0], sizes[0]);
  free(tempdata);

  meta[0] = read_volume(input_file2, &tempdata, sizes[0]);
  if (meta[0] == NULL) {
    fprintf(stderr,"ERROR! Image not read: %s\n",input_file2);
    return STATUS_ERR;
  }
  volumesize=sizes[0][0]*sizes[0][1]*sizes[0][2];

  subject2 = alloc_2d_float(3,volumesize*sizeof(**subject2));
  cp_volume(tempdata, subject2[0], sizes[0]);
  free(tempdata);

  if (read_volume(mask_file, &tempdata, tmpsizes) == NULL) {
    fprintf(stderr,"ERROR! Image not read: %s\n",mask_file);
    return STATUS_ERR;
  }

  if ((tmpsizes[0]!=sizes[0][0]) || (tmpsizes[1]!=sizes[0][1]) || (tmpsizes[2]!=sizes[0][2])) {
    fprintf(stderr,"ERROR! Mask dimension does not match image dimension!\n");
    return STATUS_ERR;
  }
  mask = alloc_2d_float(3,volumesize*sizeof(**mask));
  cp_volume(tempdata, mask[0], sizes[0]);
  free(tempdata);

  if (nomask) {
    /* option for no segmentation mask - set the mask to all ones */
    wipe_data(mask[0],sizes[0],1.0);
  }

  if (positive_file!=NULL) {
    image_metadata *positive_meta;
    if ((positive_meta=read_volume(positive_file, &tempdata, tmpsizes)) == NULL) {
      fprintf(stderr,"ERROR! Image not read: %s\n",positive_file);
      return STATUS_ERR;
    }
    if ((tmpsizes[0]!=sizes[0][0]) || (tmpsizes[1]!=sizes[0][1]) || (tmpsizes[2]!=sizes[0][2])) {
      fprintf(stderr,"ERROR! Positive mask dimension does not match image dimension!\n");
      return STATUS_ERR;
    }
    positivemask = alloc_2d_float(3,volumesize*sizeof(**mask));
    cp_volume(tempdata, positivemask[0], sizes[0]);
    free(tempdata);
    free_meta(positive_meta);

    down_sample(positivemask[0], positivemask[1], 2, sizes[0]);
    down_sample(positivemask[0], positivemask[2], 4, sizes[0]);
  }

  segmented = alloc_2d_float(3,volumesize*sizeof(**segmented));

  /* downsample the subject and mask */
  down_sample(subject1[0], subject1[1], 2, sizes[0]);
  down_sample(subject1[0], subject1[2], 4, sizes[0]);
  down_sample(subject2[0], subject2[1], 2, sizes[0]);
  down_sample(subject2[0], subject2[2], 4, sizes[0]);
  down_sample(mask[0], mask[1], 2, sizes[0]);
  down_sample(mask[0], mask[2], 4, sizes[0]);

  /* populate the entire configuration table for compatibility reasons */
  for (i=0; i<3; i++) {
    configuration[i].voxelsize = voxelsize;
    configuration[i].patchsize = sizepatch;
    configuration[i].searcharea = searcharea;
    configuration[i].alpha = alpha;
    configuration[i].beta = beta;
    configuration[i].threshold = threshold;
    configuration[i].selectionsize = selectionsize;
  }


  if (conf_file != NULL) {
    steps=read_configuration(conf_file, input_conf);
    if (steps==STATUS_ERR) {
      fprintf(stderr,"Error in configuration file. Values outside limits.\n");
      return STATUS_ERR;
    }
    initialscale=-1;
    targetscale=4;
    for (i=0; i<steps; i++) {
      scale=(int)(input_conf[i].voxelsize/2);
      configuration[scale].voxelsize=input_conf[i].voxelsize;
      configuration[scale].patchsize=input_conf[i].patchsize;
      configuration[scale].searcharea=input_conf[i].searcharea;
      configuration[scale].alpha=input_conf[i].alpha;
      configuration[scale].beta=input_conf[i].beta;
      configuration[scale].threshold=input_conf[i].threshold;
      configuration[scale].selectionsize=input_conf[i].selectionsize;
      if (scale>initialscale)
        initialscale=scale;
      if (scale<targetscale)
        targetscale=scale;
    }
  } else {
    /* if no configuration file, apply user settings for single scale */
    targetscale=initialscale=(int)(voxelsize/2);
  }

  scalesteps=initialscale-targetscale+1;

  fprintf(stderr,"%d scale steps:\n",scalesteps);

  for (i=initialscale; i>=targetscale; i--) {
    fprintf(stderr,"%d %d %d %4.2lf %4.2lf %4.2lf %d\n",configuration[i].voxelsize,configuration[i].patchsize,configuration[i].searcharea,configuration[i].alpha,configuration[i].beta,configuration[i].threshold,configuration[i].selectionsize);
  }

  images1 = alloc_3d_char(3,MAXLIBSIZE, FILENAMELENGTH);
  images2 = alloc_3d_char(3,MAXLIBSIZE, FILENAMELENGTH);
  masks = alloc_3d_char(3,MAXLIBSIZE, FILENAMELENGTH);
  means1 = alloc_3d_char(3,MAXLIBSIZE, FILENAMELENGTH);
  means2 = alloc_3d_char(3,MAXLIBSIZE, FILENAMELENGTH);
  vars1 = alloc_3d_char(3,MAXLIBSIZE, FILENAMELENGTH);
  vars2 = alloc_3d_char(3,MAXLIBSIZE, FILENAMELENGTH);

  /*for (scale=initialscale;scale>=0;scale--){*/
  for (scale=0; scale>=0; scale--) {

    sprintf(imagelist1,"%s/library.t1w.%dmm",libdir,scales[scale]);
    sprintf(imagelist2,"%s/library.flr.%dmm",libdir,scales[scale]);
    sprintf(masklist,"%s/library.masks.%dmm",libdir,scales[scale]);
    if (load_moments) {
      sprintf(meanlist1,"%s/library.meanst1w.%dmm",libdir,scales[scale]);
      sprintf(meanlist2,"%s/library.meansflr.%dmm",libdir,scales[scale]);
      sprintf(varlist1,"%s/library.varst1w.%dmm",libdir,scales[scale]);
      sprintf(varlist2,"%s/library.varsflair.%dmm",libdir,scales[scale]);
    }

    num_images=read_list(imagelist1,images1[scale],abspath?"":libdir);
    
    if (read_list(masklist,masks[scale],abspath?"":libdir)!=num_images) {
      fprintf(stderr,"ERROR! Number of t1w images and masks does not match!\n");
      return STATUS_ERR;
    }

    num_images=read_list(imagelist2,images2[scale],abspath?"":libdir);
    if (read_list(masklist,masks[scale],abspath?"":libdir)!=num_images) {
      fprintf(stderr,"ERROR! Number of flair images and masks does not match!\n");
      return STATUS_ERR;
    }

    if (num_images<configuration[scale].selectionsize) {
      fprintf(stderr,"ERROR! Cannot select more images than in the library!\n\tlibrary images: %d\n\tselection: %d\n",num_images,configuration[scale].selectionsize);
      return STATUS_ERR;
    }

    if (load_moments) {
      if (read_list(meanlist1,means1[scale],abspath?"":libdir)!=num_images) {
        fprintf(stderr,"ERROR! Number of images and means does not match!\n");
        return STATUS_ERR;
      }
      if (read_list(meanlist2,means2[scale],abspath?"":libdir)!=num_images) {
        fprintf(stderr,"ERROR! Number of images and means does not match!\n");
        return STATUS_ERR;
      }
      if (read_list(varlist1,vars1[scale],abspath?"":libdir)!=num_images) {
        fprintf(stderr,"ERROR! Number of images and vars does not match!\n");
        return STATUS_ERR;
      }
      if (read_list(varlist2,vars2[scale],abspath?"":libdir)!=num_images) {
        fprintf(stderr,"ERROR! Number of images and vars does not match!\n");
        return STATUS_ERR;
      }
    }
  }

  if ((mask_meta=read_volume(mask_file, &tempdata, tmpsizes)) == NULL) {
    fprintf(stderr,"ERROR! Image not read: %s\n",mask_file);
    return STATUS_ERR;
  }

  if ((tmpsizes[0]!=sizes[0][0]) || (tmpsizes[1]!=sizes[0][1]) || (tmpsizes[2]!=sizes[0][2])) {
    fprintf(stderr,"ERROR! Image dimension does not match library image dimension!\n");
    return STATUS_ERR;
  }
  free(tempdata);
  free_meta(mask_meta);

  fprintf(stderr,"\nrmnms -- debugging: %s\n", images1[0][0]);

  meta[1] = read_volume(images1[0][0], &tempdata, sizes[1]);
  if (meta[1] == NULL) {
    fprintf(stderr,"ERROR! Image not read: %s\n",images1[0][0]);
    return STATUS_ERR;
  }
  free(tempdata);

  meta[2] = read_volume(images2[0][0], &tempdata, sizes[1]);
  if (meta[2] == NULL) {
    fprintf(stderr,"ERROR! Image not read: %s\n",images2[0][0]);
    return STATUS_ERR;
  }
  free(tempdata);
  
    /* make the downsampled masks crisp */
  threshold_data(mask[1],sizes[1],0.5);
  
  segsubject = alloc_2d_float(3,volumesize*sizeof(**segsubject));
  patchcount = alloc_2d_float(3,volumesize*sizeof(**patchcount));
  filtered   = alloc_2d_float(3,volumesize*sizeof(**filtered));

  if (verbose) fprintf(stderr,"Initial voxel size: %d\nTarget voxel size: %d\n",scales[initialscale],scales[targetscale]);

  for (scale=initialscale; scale>=targetscale; scale--) {

    selection = (int *)malloc(configuration[scale].selectionsize*sizeof(*selection));
    pre_selection4d(subject1[scale], subject2[scale], mask[scale], images1[scale], images2[scale], sizes[scale], num_images, configuration[scale].selectionsize, selection, selection_file,verbose);

    if (verbose) fprintf(stderr,"Performing segmentation at %dmm resolution\nReading files ",scales[scale]);

    scaledvolumesize = sizes[scale][0]*sizes[scale][1]*sizes[scale][2];

    imagedata1 = (float *)malloc(configuration[scale].selectionsize*scaledvolumesize*sizeof(*imagedata1));
    imagedata2 = (float *)malloc(configuration[scale].selectionsize*scaledvolumesize*sizeof(*imagedata2));
    maskdata =  (float *)malloc(configuration[scale].selectionsize*scaledvolumesize*sizeof(*maskdata));
    meandata1 =  (float *)malloc(configuration[scale].selectionsize*scaledvolumesize*sizeof(*meandata1));
    meandata2 =  (float *)malloc(configuration[scale].selectionsize*scaledvolumesize*sizeof(*meandata2));
    vardata1 =   (float *)malloc(configuration[scale].selectionsize*scaledvolumesize*sizeof(*vardata1));
    vardata2 =   (float *)malloc(configuration[scale].selectionsize*scaledvolumesize*sizeof(*vardata2));

    /* read the library images, masks, and moments */
    for (i=0; i<configuration[scale].selectionsize; i++) {
      image_metadata *_meta;
      if (verbose) fprintf(stderr,".");
      if ((_meta=read_volume(images1[scale][selection[i]], &tempdata, tmpsizes)) == NULL) {
        fprintf(stderr,"ERROR! Image not read: %s\n",images1[scale][selection[i]]);
        return STATUS_ERR;
      }
      cp_volume(tempdata, imagedata1+i*scaledvolumesize, tmpsizes);
      free(tempdata);
      free_meta(_meta);
            
      if (verbose) fprintf(stderr,".");
      if ((_meta=read_volume(images2[scale][selection[i]], &tempdata, tmpsizes)) == NULL) {
        fprintf(stderr,"ERROR! Image not read: %s\n",images2[scale][selection[i]]);
        return STATUS_ERR;
      }
      cp_volume(tempdata, imagedata2+i*scaledvolumesize, tmpsizes);
      free(tempdata);
      free_meta(_meta);
     
    }
    if (verbose) fprintf(stderr,"*");
    for (i=0; i<configuration[scale].selectionsize; i++) {
      image_metadata *_meta;
      if (verbose) fprintf(stderr,".");
      if ((_meta=read_volume(masks[scale][selection[i]], &tempdata, tmpsizes)) == NULL) {
        fprintf(stderr,"ERROR! Image not read: %s\n",masks[scale][selection[i]]);
        return STATUS_ERR;
      }
      cp_volume(tempdata, maskdata+i*scaledvolumesize, tmpsizes);
      free(tempdata);
      free_meta(_meta);
    }
    if (verbose) fprintf(stderr,"*");

    if (!load_moments) {
      /* calculate the mean and variance for the library images */
      /* this must be done if the selected patch size is different from the one used in the precalculation */
      for (i=0; i<configuration[scale].selectionsize; i++) {
        if (verbose) fprintf(stderr,"c");
        ComputeFirstMoment(imagedata1+i*scaledvolumesize, meandata1+i*scaledvolumesize, sizes[scale], configuration[scale].patchsize, &min, &max);
        ComputeSecondMoment(imagedata1+i*scaledvolumesize, meandata1+i*scaledvolumesize, vardata1+i*scaledvolumesize, sizes[scale], configuration[scale].patchsize, &min, &max);
        ComputeFirstMoment(imagedata2+i*scaledvolumesize, meandata2+i*scaledvolumesize, sizes[scale], configuration[scale].patchsize, &min, &max);
        ComputeSecondMoment(imagedata2+i*scaledvolumesize, meandata2+i*scaledvolumesize, vardata2+i*scaledvolumesize, sizes[scale], configuration[scale].patchsize, &min, &max);
      }
    } else {
      for (i=0; i<configuration[scale].selectionsize; i++) {
        image_metadata *_meta;
        if (verbose) fprintf(stderr,".");
        if ((_meta=read_volume(means1[scale][selection[i]], &tempdata, tmpsizes)) == NULL) {
          fprintf(stderr,"ERROR! Image not read: %s\n",means1[scale][selection[i]]);
          return STATUS_ERR;
        }
        cp_volume(tempdata, meandata1+i*scaledvolumesize, tmpsizes);
        free(tempdata);
        free_meta(_meta);
      }
      if (verbose) fprintf(stderr,"*");
      for (i=0; i<configuration[scale].selectionsize; i++) {
        image_metadata *_meta;
        fprintf(stderr,".");
        if ((_meta=read_volume(vars1[scale][selection[i]], &tempdata, tmpsizes)) == NULL) {
          fprintf(stderr,"ERROR! Image not read: %s\n",vars1[scale][selection[i]]);
          return STATUS_ERR;
        }
        cp_volume(tempdata, vardata1+i*scaledvolumesize, tmpsizes);
        free(tempdata);
        free_meta(_meta);
                
        if (verbose) fprintf(stderr,".");
        if ((_meta=read_volume(means2[scale][selection[i]], &tempdata, tmpsizes)) == NULL) {
          fprintf(stderr,"ERROR! Image not read: %s\n",means2[scale][selection[i]]);
          return STATUS_ERR;
        }
        cp_volume(tempdata, meandata2+i*scaledvolumesize, tmpsizes);
        free(tempdata);
        free_meta(_meta);
      }
      if (verbose) fprintf(stderr,"*");
      for (i=0; i<configuration[scale].selectionsize; i++) {
        image_metadata *_meta;
        fprintf(stderr,".");
        if ((_meta=read_volume(vars2[scale][selection[i]], &tempdata, tmpsizes)) == NULL) {
          fprintf(stderr,"ERROR! Image not read: %s\n",masks[scale][selection[i]]);
          return STATUS_ERR;
        }
        cp_volume(tempdata, vardata2+i*scaledvolumesize, tmpsizes);
        free(tempdata);
        free_meta(_meta);
      }
    }
    if (verbose) fprintf(stderr,"\n");
    /* end of reading files */

    /* remove any disconnected parts */
    masksize = getLargestObject_float(mask[scale], sizes[scale], 1, 0);

    if (verbose) fprintf(stderr,"Mask size: %d\nAlpha: %f\n",masksize,configuration[scale].alpha);

    /* make sure we starting from a clean slate */
    wipe_data(segsubject[scale],sizes[scale],0.0);

    if (flipimages) {
      /* doubling the library selection by flipping images along the mid-sagittal plane */
      imagedata1 = (float *)realloc(imagedata1,configuration[scale].selectionsize*2*scaledvolumesize*sizeof(*imagedata1));
      imagedata2 = (float *)realloc(imagedata2,configuration[scale].selectionsize*2*scaledvolumesize*sizeof(*imagedata2));
      maskdata = (float *)realloc(maskdata,configuration[scale].selectionsize*2*scaledvolumesize*sizeof(*maskdata));
      meandata1 = (float *)realloc(meandata1,configuration[scale].selectionsize*2*scaledvolumesize*sizeof(*meandata1));
      meandata2 = (float *)realloc(meandata2,configuration[scale].selectionsize*2*scaledvolumesize*sizeof(*meandata2));
      vardata1 = (float *)realloc(vardata1,configuration[scale].selectionsize*2*scaledvolumesize*sizeof(*vardata1));
      vardata2 = (float *)realloc(vardata2,configuration[scale].selectionsize*2*scaledvolumesize*sizeof(*vardata2));

      for (i=0; i<configuration[scale].selectionsize; i++) {
        flip_data(imagedata1+i*scaledvolumesize, imagedata1+(configuration[scale].selectionsize+i)*scaledvolumesize, sizes[scale]);
        flip_data(imagedata2+i*scaledvolumesize, imagedata2+(configuration[scale].selectionsize+i)*scaledvolumesize, sizes[scale]);
        flip_data(maskdata+i*scaledvolumesize, maskdata+(configuration[scale].selectionsize+i)*scaledvolumesize, sizes[scale]);
        flip_data(meandata1+i*scaledvolumesize, meandata1+(configuration[scale].selectionsize+i)*scaledvolumesize, sizes[scale]);
        flip_data(meandata2+i*scaledvolumesize, meandata2+(configuration[scale].selectionsize+i)*scaledvolumesize, sizes[scale]);
        flip_data(vardata1+i*scaledvolumesize, vardata1+(configuration[scale].selectionsize+i)*scaledvolumesize, sizes[scale]);
        flip_data(vardata2+i*scaledvolumesize, vardata2+(configuration[scale].selectionsize+i)*scaledvolumesize, sizes[scale]);
      }

      max = nlmsegFuzzy4Drmnmns(subject1[scale],subject2[scale], imagedata1, imagedata2, maskdata, meandata1, meandata2, vardata1, vardata2, mask[scale], configuration[scale].patchsize, configuration[scale].searcharea, configuration[scale].beta, configuration[scale].threshold, sizes[scale], configuration[scale].selectionsize*2, segsubject[scale], patchcount[scale]);
    } else {
      max = nlmsegFuzzy4Drmnmns(subject1[scale],subject2[scale], imagedata1, imagedata2, maskdata, meandata1, meandata2, vardata1, vardata2, mask[scale], configuration[scale].patchsize, configuration[scale].searcharea, configuration[scale].beta, configuration[scale].threshold, sizes[scale], configuration[scale].selectionsize, segsubject[scale], patchcount[scale]);
    }

    free(imagedata1);
    free(imagedata2);
    free(maskdata);
    free(meandata1);
    free(meandata2);
    free(vardata1);
    free(vardata2);


/*    if (positive_file!=NULL) {
      /* add the certain positive segmentation (inside mask) */
/*      add_mask_data(segsubject[scale],positivemask[scale],sizes[scale]);      
/*    }

    /* add the certain segmentation from the previous scale */
    /*add_mask_data(segsubject[scale],segmented[scale],sizes[scale]);*/

    if (medianfilter) {
      median_filter(segsubject[scale], sizes[scale], 3);
    }

    /* the patch filter is experimental */
/*    if (patchfilter) {
      wipe_data(filtered[scale],sizes[scale],0.0);
      wipe_data(patchcount[scale],sizes[scale],0.0);
      max = nlmfilter(subject[scale], mask[scale], segsubject[scale], 2*configuration[scale].patchsize, 2*configuration[scale].searcharea, configuration[scale].beta, configuration[scale].threshold, sizes[scale], filtered[scale], patchcount[scale]);
      combine_maps(segsubject[scale], filtered[scale], mask[scale], sizes[scale]);
    }*/

    if (scale > targetscale) {
      /* if performing a higher resolution step, upsample the result and create new mask */
      resize_trilinear(segsubject[scale], sizes[scale], sizes[scale-1], segsubject[scale-1]);
      masksize=update_mask(segsubject[scale-1], mask[scale-1], segmented[scale-1], sizes[scale-1], configuration[scale].alpha, 1.0-configuration[scale].alpha);
    }

    free(selection);
  } // for each scale


  if (count_file!=NULL) {
    write_volume_generic(count_file, patchcount[targetscale], meta[targetscale],FALSE);
  }

  if(targetscale!=0 && same_res) { /* need to upsample final output */
    if (verbose) fprintf(stderr,"Upsampling to input resolution, %dx%dx%d\n",sizes[0][0],sizes[0][1],sizes[0][2]);
    resize_trilinear(segsubject[targetscale], sizes[targetscale], sizes[0], segsubject[0]);
    masksize=update_mask(segsubject[0], mask[0], segmented[0], sizes[0], configuration[targetscale].alpha, 1.0-configuration[targetscale].alpha);
    targetscale=0;
    configuration[targetscale].alpha = alpha;
  }

  if (!outputprob) {
    if (verbose) fprintf(stderr,"Thresholding estimator at %f\n",configuration[targetscale].alpha);
    threshold_data(segsubject[targetscale], sizes[targetscale], configuration[targetscale].alpha);
    /*getLargestObject_float(segsubject[targetscale], sizes[targetscale], 1, 0);*/

    if (fill_output) {
      wipe_data(mask[targetscale], sizes[targetscale], 1.0);
      filled = flood_fill_float(segsubject[targetscale], mask[targetscale], sizes[targetscale], 0, 0, 0, 0, 6);
      //segsubject[targetscale]=mask[targetscale];
      cp_volume(mask[targetscale],segsubject[targetscale],sizes[targetscale]);
    }
  }

  meta[targetscale]->history=strdup(history_label);
  if(write_volume_generic(output_file, segsubject[targetscale], meta[targetscale],!outputprob)) {
    fprintf(stderr,"Can't save output to %s\n",output_file);
    return STATUS_ERR;
  }

  free_2d_float(mask);
  free_2d_float(subject1);
  free_2d_float(subject2);
  if (positive_file!=NULL)
    free_2d_float(positivemask);
  
  free_2d_float(filtered);
  free_2d_float(segmented);
  free_2d_float(segsubject);
  free_2d_float(patchcount);

  free_3d_char(images1);
  free_3d_char(images2);
  free_3d_char(masks);
  free_3d_char(means1);
  free_3d_char(means2);
  free_3d_char(vars1);
  free_3d_char(vars2);
  
  free_meta(meta[2]);
  free_meta(meta[1]);
  
  free(meta);

  return STATUS_OK;
}
