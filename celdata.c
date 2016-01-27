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

#include <stdio.h>
#include "cel.h"

void init_CELdata(CELdata *d){
  d->valid = 0;
  d->type = CEL_TYPE_UNKNOWN;
  d->intensity_stats_calculated = 0;
  d->array = NULL;
  d->algorithm = NULL;
  d->rows = 0;
  d->cols = 0;
  d->cell_margin = 0;
  d->outliers = 0;
  d->masked = 0;
  d->intensity_min = 0;
  d->intensity_max = 0;
  d->intensity_n_unique = 0;
  d->intensity_n_invalid = 0;
}

void free_CELdata(CELdata *d){
  d->valid = 0;
  d->type = CEL_TYPE_UNKNOWN;
  d->intensity_stats_calculated = 0;
  if(d->array != NULL){
    free(d->array);
    d->array = NULL;
  }
  if(d->algorithm != NULL){
    free(d->algorithm);
    d->algorithm = NULL;
  }
  d->rows = 0;
  d->cols = 0;
  d->cell_margin = 0;
  d->outliers = 0;
  d->masked = 0;
  d->intensity_min = 0;
  d->intensity_max = 0;
  d->intensity_n_unique = 0;
  d->intensity_n_invalid = 0;
}
#define CEL_TYPE_UNKNOWN 100
#define CEL_TYPE_BINARY 101
#define CEL_TYPE_CALVIN 102
#define CEL_TYPE_TEXT 103

void print_CELdata(CELdata *d){
  char *type_str = "unknown";
  if((d->valid != 1) || (d->type == CEL_TYPE_UNKNOWN)){
    printf("(invalid)\n");
    return;
  }
  if(d->type == CEL_TYPE_BINARY) type_str = "binary";
  else if(d->type == CEL_TYPE_CALVIN) type_str = "calvin";
  else if(d->type == CEL_TYPE_TEXT) type_str = "text";
  if(d->intensity_stats_calculated != 1) printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n", type_str, d->array, d->algorithm, d->rows, d->cols, d->cell_margin, d->outliers, d->masked);
  else printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%0.0f\t%0.0f\t%d\t%d\n", type_str, d->array, d->algorithm, d->rows, d->cols, d->cell_margin, d->outliers, d->masked, d->intensity_min, d->intensity_max, d->intensity_n_unique, d->intensity_n_invalid);
}

void extract_chipname(char *str, CELdata *d){
  size_t array_length;
  char *pointer, *str_end;
  if(d->array != NULL) free(d->array);
  str_end = strstr((char*)str, ".1sq");
  if(str_end != NULL){
    pointer = str_end;
    while(pointer[0] != ' ') pointer--;
    pointer++;
    array_length = str_end - pointer + 1;
  } else {
    pointer = "unknown";
    array_length = strlen(pointer) + 1;
  }
  d->array = malloc(array_length * sizeof(char));
  memcpy(d->array, pointer, array_length);
  d->array[array_length - 1] = '\0';
}

void calculate_intensity_stats(float *data, size_t n, float *max_value, float *min_value, int *unique, int *invalid){
  int i;
  int unique_values[MAX_INTENSITY_VALUE + 1];
  float curr_value;
  if(data == NULL) return;
  *min_value = MAX_INTENSITY_VALUE + 1;
  *max_value = -1;
  *unique = 0;
  *invalid = 0;
  memset(&unique_values, 0, (MAX_INTENSITY_VALUE + 1) * sizeof(int));
  for(i=0; i < n; i++){
    curr_value = data[i];
    if(curr_value < 0){
      (*invalid)++;
      continue;
    }
    if(curr_value > MAX_INTENSITY_VALUE){
      (*invalid)++;
      continue;
    }
    if(curr_value < *min_value) *min_value = curr_value;
    if(curr_value > *max_value) *max_value = curr_value;
    unique_values[(int)round(curr_value)] ++;    
  }
  for(i=0; i<MAX_INTENSITY_VALUE + 1; i++) if(unique_values[i] != 0) (*unique)++;
  if(*invalid == n){
    *unique = 0;
    *min_value = 0.0;
    *max_value = 0.0;
  }
}
