// checkcel - check the validity of an Affymetrix CEL file.
// Copyright (C) 2016 Alastair Droop, The Leeds MRC Medical Bioinformatics Centre <a.p.droop@leeds.ac.uk>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __checkcel_celdata_h
#define __checkcel_celdata_h

//Define the maximum CEL mean value:
#define MAX_INTENSITY_VALUE 65535

// Define the struct that holds data from a CELfile:
typedef struct {
  char valid;
  char type;
  char intensity_stats_calculated;
  char *array;
  char *algorithm;
  int32_t rows;
  int32_t cols;
  int32_t cell_margin;
  u_int32_t outliers;
  u_int32_t masked;
  float intensity_min;
  float intensity_max;
  int intensity_n_unique;
  int intensity_n_invalid;
} CELdata;


void init_CELdata(CELdata *d);
void free_CELdata(CELdata *d);
void print_CELdata(CELdata *d);

//Extract the chip name from a given string:
void extract_chipname(char *str, CELdata *cel_data);

// Calculate statistics from an array of intensity values:
void calculate_intensity_stats(float *data, size_t n, float *max_value, float *min_value, int *unique, int *invalid);

#endif
