/*  rmnms.h
 *
 *  Copyright 2011  Simon Fristed Eskildsen, Vladimir Fonov,
 *   	      	    Pierrick Coupé, Jose V. Manjon
 *
 *  This file is part of mincrmnms.
 *
 *  mincrmnms is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  mincrmnms is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with mincrmnms.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  For questions and feedback, please contact:
 *  Simon Fristed Eskildsen <eskild@gmail.com> 
 */


#ifndef RMNMS_H
#define RMNMS_H

#include <stdio.h>

#ifdef HAVE_MINC
#ifndef DEF_VOLUME_IO
#include "mincio.h"
#endif
#endif //HAVE_MINC

#include "basic.h"

#define MAXLIBSIZE 1000
#define FILENAMELENGTH 255

#define VOXELSIZEMAX 4
#define VOXELSIZEMIN 1
#define PATCHSIZEMAX 10
#define PATCHSIZEMIN 0
#define SEARCHAREAMAX 10
#define SEARCHAREAMIN 0
#define ALPHAMAX 1.0
#define ALPHAMIN 0.0
#define BETAMAX 5.0
#define BETAMIN 0.0
#define THRESHOLDMAX 1.0
#define THRESHOLDMIN 0.0

typedef struct {
  int index;
  float ssd;
} ssd_t;

typedef struct {
  int voxelsize;
  int patchsize;
  int searcharea;
  double alpha;
  double beta;
  double threshold;
  int selectionsize;
} rmnms_conf;

int fgetline(FILE *fp, char line[], int max);

int median_filter(float *volume, int *sizes, int filtersize);

int trilinear_interpolant(float *volume, int *sizes, point3D coord, float *result);
int resize_volume(float *input, int *sizes, int *sizes2, float *result);
int resize_trilinear(float *input, int *sizes, int *sizes2, float *result);

void cp_volume(float *data, float *copy, int *sizes);

int flip_data(float *data, float *result, int *sizes);
int combine_maps(float *data, float *map, float *mask, int *sizes);

int down_sample(float *subject, float *result, int factor, int *sizes);
int up_sample(float *subject, float *result, int factor, int *sizes, int *targetsizes);

int threshold_data(float *data, int *sizes, float threshold);

int add_mask_data(float *data1, float *mask, int *sizes);
int wipe_data(float *data1, int *sizes, float value);
int update_mask(float *subject, float *mask, float *segmented, int *sizes, float min, float max);

int flood_fill_float(float *data, float *output, int *sizes, int sx, int sy, int sz, float fill_value, int connectivity);

int pre_selection(float *subject, float *mask, char **images, int *sizes, int librarysize, int num_selected, int *selection, char *outfile, VIO_BOOL verbose);
int pre_selection4d(float *subject1, float *subject2, float *mask, char **images1, char **images2, int *sizes, int librarysize, int num_selected, int *selection, char *outfile, VIO_BOOL verbose);

int read_configuration(char *filename, rmnms_conf *conf);
int read_list(char *filename, char **list,char *basedir);

image_metadata * read_volume(char *filename, float **data, int *sizes);
int write_volume_generic(char *filename, float *data, image_metadata *meta,VIO_BOOL binary_mask );

#endif
