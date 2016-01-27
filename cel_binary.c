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

char is_CELbinary(CELfile f){
  char bitflip = 0;
  if(check_endian() == MACHINE_BIG_ENDIAN) bitflip = 1;
  char result;
  int32_t magic_number, version;
  reset_CELfile(f);
  // Check the magic number:
  result = readCEL_int32(&magic_number, 1, f, bitflip);
  if((result != CEL_READ_VALUE_OK) || (magic_number != 64)) return 0;
  // Check the magic number:
  result = readCEL_int32(&version, 1, f, bitflip);
  if((result != CEL_READ_VALUE_OK) || (version != 4)) return 0;
  // All initial checks look good, to reset the file and return success:
  reset_CELfile(f);
  return 1;
}

char readCELbinary(CELfile f, CELdata *d, char read_intensity, char verbose){
  char bitflip = 0;
  char result;
  int32_t i, magic_number, version, cells, subgrids;
  size_t n;
  float *intensities;
  CELbinary_spotdata *spotdata, *p;
  char *header = NULL;
  char *parameters = NULL;
  //Sort out the endianness of the machine we're on:
  if(check_endian() == MACHINE_BIG_ENDIAN) bitflip = 1;
  reset_CELfile(f);
  // Initially, set the data to invalid:
  d->valid = 0;
  // Check the magic number:
  result = readCEL_int32(&magic_number, 1, f, bitflip);
  if((result != CEL_READ_VALUE_OK) || (magic_number != 64)) return 1;
  // Check the magic number:
  result = readCEL_int32(&version, 1, f, bitflip);
  if((result != CEL_READ_VALUE_OK) || (version != 4)) return 1;
  // Get the chip dimensions:
  if(readCEL_int32(&d->cols, 1, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  if(readCEL_int32(&d->rows, 1, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  if(readCEL_int32(&cells, 1, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  if(cells != d->rows * d->cols) return 1;
  // Read in the file header:
  if(readCEL_str(&header, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  extract_chipname(header, d);
  if(verbose == 1) printf("header: \"%s\"\n", header);
  free(header); //APD 10/08
  header = NULL; //APD 10/08
  //Read in the algorithm:
  if(readCEL_str(&d->algorithm, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  for(i=0; i<strlen(d->algorithm); i++) d->algorithm[i] = tolower(d->algorithm[i]);
  if(verbose == 1) printf("algorithm: \"%s\"\n", d->algorithm);
  //Read in the parameters:
  if(readCEL_str(&parameters, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  if(verbose == 1) printf("parameters: \"%s\"\n", parameters);
  free(parameters); //APD 10/08
  parameters = NULL; //APD 10/08
  // Read in the cell margin:
  if(readCEL_int32(&d->cell_margin, 1, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  // Read in the outlier number:
  if(readCEL_uint32(&d->outliers, 1, f, bitflip) != CEL_READ_VALUE_OK) return 1;  
  // Read in the masked cell count:
  if(readCEL_uint32(&d->masked, 1, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  // Read in the subgrid number:
  if(readCEL_int32(&subgrids, 1, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  // Read in the intensity data if needed:
  if(read_intensity ==1){
    intensities = (float*)malloc(cells * sizeof(float));
    if(intensities == NULL) return 1;
    spotdata = (CELbinary_spotdata*)malloc(cells * sizeof(CELbinary_spotdata));
    if(spotdata == NULL){
      free(intensities);
      intensities = NULL;
      return 1;
    }
    n = fread(spotdata, sizeof(CELbinary_spotdata), cells, f.handle);
    if(n != cells){
      free(intensities);
      intensities = NULL;
      free(spotdata);
      spotdata = NULL;
      return 1;
    }
    p = spotdata;
    for(i=0; i<cells; i++){
      intensities[i] = p->intensity;
      p++;
    }
    free(spotdata);
    spotdata = NULL;
    calculate_intensity_stats(intensities, cells, &d->intensity_max, &d->intensity_min, &d->intensity_n_unique, &d->intensity_n_invalid);
    d->intensity_stats_calculated = 1;
    free(intensities);
    intensities = NULL;    
  }
  // Mark the CEL data as valid:
  d->type = CEL_TYPE_BINARY;
  d->valid = 1;
  return 1;
}
