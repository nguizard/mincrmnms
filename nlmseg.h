/*  nlmseg.h
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


#ifndef NLMSEG_H
#define NLMSEG_H

void ComputeFirstMoment(float* ima, float* means, const int* dims, int f, float *min, float *max);
void ComputeSecondMoment(float* ima, float* means, float* variance, const int* dims,int f, float *min, float *max);
void ComputeFirstMoment4D(float* ima,float* atlas, float* means, float* Ameans,const int* dims, int f);
void ComputeSecondMoment4D(float* ima,float* atlas, float* means, float* Ameans, float* variance, float* Avariance,const int* dims,int f);

float SSDPatch(float* PatchImg, float* PatchTemplate, int f);
void ExtractPatch(float* ima, float* Patch, int x, int y, int z, int size, int sx, int sy, int sz);
void ExtractPatch4D(float* ima, float* Patch, int x,int y, int z, int t,int size,int sx,int sy,int sz);

float nlmsegFuzzy4D(float *subject, float *imagedata, float *maskdata, float *meandata, float *vardata, float *mask, int sizepatch, int searcharea, float alpha, float threshold, int sizes[3], int librarysize, float *SegSubject, float *PatchCount);
float nlmsegFuzzy4Drmnmns(float *subject1, float *subject2, float *imagedata1, float *imagedata2, float *maskdata, float *meandata1, float *meandata2, float *vardata1, float *vardata2, float *mask, int sizepatch, int searcharea, float alpha, float threshold, int sizes[3], int librarysize, float *SegSubject, float *PatchCount);
float nlmfilter(float *subject, float *mask, float *maskdata, int sizepatch, int searcharea, float beta, float threshold, int dims[3], float *SegSubject, float *PatchCount);

#endif
