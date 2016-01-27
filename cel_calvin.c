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

char readCELcalvin_parameter(CELcalvin_parameter *p, CELfile f, char bitflip){
  char result;
  char* type = NULL;
  result = readCEL_wstr(&(p->name), f, bitflip);
  if(result != CEL_READ_VALUE_OK){
    free(p->name);
    p->name = NULL;
    return CEL_READ_VALUE_FAILED;    
  }
  result = readCEL_int32(&(p->value_length), 1, f, bitflip);
  if(result != CEL_READ_VALUE_OK){
    free(p->value);
    free(p->name);
    p->name = NULL;
    p->value = NULL;
    return CEL_READ_VALUE_FAILED;    
  }
  result = readCEL_char(&(p->value), p->value_length, f);
  if(result != CEL_READ_VALUE_OK){
    free(p->value);
    free(p->name);
    p->name = NULL;
    p->value = NULL;
    return CEL_READ_VALUE_FAILED;    
  }
  result = readCEL_wstr(&type, f, bitflip);
  if(result != CEL_READ_VALUE_OK){
    free(p->value);
    free(p->name);
    p->name = NULL;
    p->value = NULL;
    free(type);
    type = NULL;
    return CEL_READ_VALUE_FAILED;    
  }
  p->type = CEL_CALVIN_MIMETYPE_UNKNOWN;
  if(strcmp(type, "text/x-calvin-integer-8") == 0) p->type = CEL_CALVIN_MIMETYPE_INT8;
  else if(strcmp(type, "text/x-calvin-unsigned-integer-8") == 0) p->type = CEL_CALVIN_MIMETYPE_UINT8;
  else if(strcmp(type, "text/x-calvin-integer-16") == 0) p->type = CEL_CALVIN_MIMETYPE_INT16;
  else if(strcmp(type, "text/x-calvin-unsigned-integer-16") == 0) p->type = CEL_CALVIN_MIMETYPE_UINT16;
  else if(strcmp(type, "text/x-calvin-integer-32") == 0) p->type = CEL_CALVIN_MIMETYPE_INT32;
  else if(strcmp(type, "text/x-calvin-unsigned-integer-32") == 0) p->type = CEL_CALVIN_MIMETYPE_UINT32;
  else if(strcmp(type, "text/x-calvin-float") == 0) p->type = CEL_CALVIN_MIMETYPE_FLOAT;
  else if(strcmp(type, "text/plain") == 0) p->type = CEL_CALVIN_MIMETYPE_PLAINTEXT;
  else if(strcmp(type, "text/ascii") == 0) p->type = CEL_CALVIN_MIMETYPE_ASCII;
  free(type);
  type = NULL;
  return CEL_READ_VALUE_OK;
}

void freeCELcalvin_parameter(CELcalvin_parameter *p){
  if(p->name != NULL){
    free(p->name);
    p->name = NULL;
  }
  if(p->value != NULL){
    free(p->value);
    p->value = NULL;
  }
  p->type = CEL_CALVIN_MIMETYPE_UNKNOWN;
  p->value_length = 0;
}

int8_t decode_CELcalvin_parameter_int8(CELcalvin_parameter *p){
  int8_t decoded_value;
  memcpy(&decoded_value, p->value, sizeof(int8_t));
  return decoded_value;
}

u_int8_t decode_CELcalvin_parameter_uint8(CELcalvin_parameter *p){
  u_int8_t decoded_value;
  memcpy(&decoded_value, p->value, sizeof(u_int8_t));
  return decoded_value;  
}

int16_t decode_CELcalvin_parameter_int16(CELcalvin_parameter *p, char bitflip){
  int16_t decoded_value;
  memcpy(&decoded_value, p->value, sizeof(int16_t));
  if(bitflip == 1) decoded_value = (((decoded_value>>8)&0xff) | ((decoded_value&0xff)<<8));
  return decoded_value;
}

u_int16_t decode_CELcalvin_parameter_uint16(CELcalvin_parameter *p, char bitflip){
  int16_t decoded_value;
  memcpy(&decoded_value, p->value, sizeof(int16_t));
  if(bitflip == 1) decoded_value = (((decoded_value>>8)&0xff) | ((decoded_value&0xff)<<8));
  return decoded_value;
}

int32_t decode_CELcalvin_parameter_int32(CELcalvin_parameter *p, char bitflip){
  int32_t decoded_value;
  memcpy(&decoded_value, p->value, sizeof(int32_t));
  if(bitflip == 1) decoded_value = (((decoded_value>>24)&0xff) | ((decoded_value&0xff)<<24) | ((decoded_value>>8)&0xff00) | ((decoded_value&0xff00)<<8));
  return decoded_value;
}

u_int32_t decode_CELcalvin_parameter_uint32(CELcalvin_parameter *p, char bitflip){
  u_int32_t decoded_value;
  memcpy(&decoded_value, p->value, sizeof(u_int32_t));
  if(bitflip == 1) decoded_value = (((decoded_value>>24)&0xff) | ((decoded_value&0xff)<<24) | ((decoded_value>>8)&0xff00) | ((decoded_value&0xff00)<<8));
  return decoded_value;  
}

float decode_CELcalvin_parameter_float(CELcalvin_parameter *p, char bitflip){
  u_int32_t data;
  float decoded_value;
  memcpy(&data, p->value, sizeof(u_int32_t));
  if(bitflip == 1) data = (((data>>24)&0xff) | ((data&0xff)<<24) | ((data>>8)&0xff00) | ((data&0xff00)<<8));
  memcpy(&decoded_value, &data, sizeof(u_int32_t));
  return decoded_value;
}

void decode_CELcalvin_parameter_plaintext(CELcalvin_parameter *p, char **s){
  int i;
  *s = (char*)malloc(((p->value_length / 2) + 1) * sizeof(char));
  memset(*s, 0, ((p->value_length / 2) + 1) * sizeof(char));
  for(i=0; i<(p->value_length / 2); i++) (*s)[i] = p->value[1 + (i * 2)];
}

void printCELcalvin_parameter(CELcalvin_parameter *p){
  char *buffer = NULL;
  if(p->type == CEL_CALVIN_MIMETYPE_INT8) printf("parameter (%d) \"%s\" = %d\n", p->type, p->name, decode_CELcalvin_parameter_int8(p));
  else if(p->type == CEL_CALVIN_MIMETYPE_UINT8) printf("parameter (%d) \"%s\" = %d\n", p->type, p->name, decode_CELcalvin_parameter_uint8(p));
  else if(p->type == CEL_CALVIN_MIMETYPE_INT16) printf("parameter (%d) \"%s\" = %d\n", p->type, p->name, decode_CELcalvin_parameter_int16(p, 1));
  else if(p->type == CEL_CALVIN_MIMETYPE_UINT16) printf("parameter (%d) \"%s\" = %d\n", p->type, p->name, decode_CELcalvin_parameter_uint16(p, 1));
  else if(p->type == CEL_CALVIN_MIMETYPE_INT32) printf("parameter (%d) \"%s\" = %d\n", p->type, p->name, decode_CELcalvin_parameter_int32(p, 1));
  else if(p->type == CEL_CALVIN_MIMETYPE_UINT32) printf("parameter (%d) \"%s\" = %d\n", p->type, p->name, decode_CELcalvin_parameter_uint32(p, 1));
  else if(p->type == CEL_CALVIN_MIMETYPE_FLOAT) printf("parameter (%d) \"%s\" = %f\n", p->type, p->name, decode_CELcalvin_parameter_float(p, 1));
  else if(p->type == CEL_CALVIN_MIMETYPE_PLAINTEXT){
    decode_CELcalvin_parameter_plaintext(p, &buffer);
    printf("parameter (%d) \"%s\" = \"%s\"\n", p->type, p->name, buffer);
    free(buffer);
    buffer = NULL;
  }
  else printf("parameter (%d) \"%s\" = \"%s\"\n", p->type, p->name, p->value);
}

char readCELcalvin_datagroup(CELcalvin_datagroup *g, CELfile f, char bitflip){
  if(readCEL_uint32(&g->next_pos, 1, f, bitflip) != CEL_READ_VALUE_OK) return CEL_READ_VALUE_FAILED;
  if(readCEL_uint32(&g->first_dataset_pos, 1, f, bitflip) != CEL_READ_VALUE_OK) return CEL_READ_VALUE_FAILED;
  if(readCEL_int32(&g->dataset_number, 1, f, bitflip) != CEL_READ_VALUE_OK) return CEL_READ_VALUE_FAILED;
  if(readCEL_wstr(&g->name, f, bitflip) != CEL_READ_VALUE_OK){
    free(g->name);
    g->name = NULL;
    return CEL_READ_VALUE_FAILED;
  }
  return CEL_READ_VALUE_OK;
}

void free_CELcalvin_datagroup(CELcalvin_datagroup *g){
  g->next_pos = 0;
  g->first_dataset_pos = 0;
  g->dataset_number = 0;
  if(g->name != NULL){
    free(g->name);
    g->name = NULL;
  }
}

char readCELcalvin_dataset_column(CELcalvin_dataset_column *c, CELfile f, char bitflip){
  if(readCEL_wstr(&c->name, f, bitflip) != CEL_READ_VALUE_OK){
    free(c->name);
    c->name = NULL;
    return CEL_READ_VALUE_FAILED;
  }
  if(readCEL_int8(&c->type, 1, f) != CEL_READ_VALUE_OK) return CEL_READ_VALUE_FAILED;
  if(readCEL_int32(&c->size, 1, f, bitflip) != CEL_READ_VALUE_OK) return CEL_READ_VALUE_FAILED;
  return CEL_READ_VALUE_OK;
}

void free_CELcalvin_dataset_column(CELcalvin_dataset_column *c){
  c->type = 0;
  c->size = 0;
  if(c->name != NULL){
    free(c->name);
    c->name = NULL;
  }
}

char readCELcalvin_dataset(CELcalvin_dataset *g, CELfile f, char bitflip){
  int i;
  if(readCEL_uint32(&g->first_element_pos, 1, f, bitflip) == CEL_READ_VALUE_FAILED) return CEL_READ_VALUE_FAILED;
  if(readCEL_uint32(&g->next_dataset_pos, 1, f, bitflip) == CEL_READ_VALUE_FAILED) return CEL_READ_VALUE_FAILED;
  if(readCEL_wstr(&g->name, f, bitflip) == CEL_READ_VALUE_FAILED){
    free(g->name);
    g->name = NULL;
    return CEL_READ_VALUE_FAILED;
  }
  if(readCEL_int32(&g->parameter_number, 1, f, bitflip) == CEL_READ_VALUE_FAILED) return CEL_READ_VALUE_FAILED;
  g->parameters = (CELcalvin_parameter*)malloc(g->parameter_number * sizeof(CELcalvin_parameter));
  if(g->parameters == NULL){
    free(g->name);
    g->name = NULL;
    return CEL_READ_VALUE_FAILED;    
  }
  for(i=0; i<g->parameter_number; i++) readCELcalvin_parameter(&(g->parameters[i]), f, bitflip);
  if(readCEL_uint32(&(g->column_number), 1, f, bitflip) == CEL_READ_VALUE_FAILED){
    for(i=0; i<g->parameter_number; i++) freeCELcalvin_parameter(&(g->parameters[i]));
    free(g->parameters);
    g->parameters = NULL;
    free(g->name);
    g->name = NULL;
    return CEL_READ_VALUE_FAILED;
  }
  g->columns = (CELcalvin_dataset_column*)malloc(g->column_number * sizeof(CELcalvin_dataset_column));
  if(g->columns == NULL){
    for(i=0; i<g->parameter_number; i++) freeCELcalvin_parameter(&(g->parameters[i]));
    free(g->parameters);
    g->parameters = NULL;
    free(g->name);
    g->name = NULL;
    return CEL_READ_VALUE_FAILED;
  }
  for(i=0; i<g->column_number; i++) readCELcalvin_dataset_column(&(g->columns[i]), f, bitflip);
  if(readCEL_uint32(&(g->row_number), 1, f, bitflip) == CEL_READ_VALUE_FAILED) return CEL_READ_VALUE_FAILED;
  return CEL_READ_VALUE_OK;
}

void free_CELcalvin_dataset(CELcalvin_dataset *g){
  int i;
  if(g->name != NULL){
    free(g->name);
    g->name = NULL;
  }
  if(g->parameters != NULL){
    for(i=0; i<g->parameter_number; i++) freeCELcalvin_parameter(&(g->parameters[i]));
    free(g->parameters);
    g->parameters = NULL;
  }
  if(g->columns != NULL){
    for(i=0; i<g->column_number; i++) free_CELcalvin_dataset_column(&(g->columns[i]));
    free(g->columns);
    g->columns = NULL;
  }
}

char is_CELcalvin(CELfile f){
  char bitflip = 0;
  if(check_endian() == MACHINE_LITTLE_ENDIAN) bitflip = 1;
  char result;
  u_int8_t magic_number, version;
  int32_t group_number;
  u_int32_t first_group_offset;
  char *data_type = NULL;
  reset_CELfile(f);
  // Check the magic number:
  result = readCEL_uint8(&magic_number, 1, f);
  if((result != CEL_READ_VALUE_OK) || (magic_number != 59)){
    reset_CELfile(f);
    return 0;
  }
  // Check the version number:
  result = readCEL_uint8(&version, 1, f);
  if((result != CEL_READ_VALUE_OK) || (version != 1)){
    reset_CELfile(f);
    return 0;
  }
  // Read in the group count:
  if(readCEL_int32(&group_number, 1, f, bitflip) != CEL_READ_VALUE_OK){
    reset_CELfile(f);
    return 0;
  }
  // Read in the offset of the first group:
  if(readCEL_uint32(&first_group_offset, 1, f, bitflip) != CEL_READ_VALUE_OK){
    reset_CELfile(f);
    return 0;
  }
  // Check the data type identifier:
  result = readCEL_str(&data_type, f, bitflip);
  if(result != CEL_READ_VALUE_OK){
    free(data_type);
    data_type = NULL;
    reset_CELfile(f);
    return 0;
  }
  if(strcmp(data_type, "affymetrix-calvin-intensity") != 0){
    free(data_type);
    data_type = NULL;
    reset_CELfile(f);
    return 0;
  }
  free(data_type);
  data_type = NULL;
  // All initial checks look good, to reset the file and return success:
  reset_CELfile(f);
  return 1;
}

char readCELcalvin(CELfile f, CELdata *d, char read_intensity, char verbose){
  char bitflip = 0;
  if(check_endian() == MACHINE_LITTLE_ENDIAN) bitflip = 1;
  char result;
  char *data_type = NULL;
  char *file_id = NULL;
  char *date = NULL;
  char *locale = NULL;
  float *intensities;
  u_int8_t magic_number, version;
  int32_t group_number, parameter_number;
  u_int32_t first_group_offset;
  CELcalvin_parameter parameter;
  int i, j;
  reset_CELfile(f);
  // Initially, set the data to invalid:
  d->valid = 0;
  // Check the magic number:
  result = readCEL_uint8(&magic_number, 1, f);
  if((result != CEL_READ_VALUE_OK) || (magic_number != 59)) return 1;
  // Check the version number:
  result = readCEL_uint8(&version, 1, f);
  if((result != CEL_READ_VALUE_OK) || (version != 1)) return 1;
  // Read in the group count:
  if(readCEL_int32(&group_number, 1, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  // Read in the offset of the first group:
  if(readCEL_uint32(&first_group_offset, 1, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  // Check the data type identifier:
  result = readCEL_str(&data_type, f, bitflip);
  if(result != CEL_READ_VALUE_OK){
    free(data_type);
    data_type = NULL;
    return 1;
  }
  if(strcmp(data_type, "affymetrix-calvin-intensity") != 0){
    free(data_type);
    data_type = NULL;
    return 1;
  }
  free(data_type);
  data_type = NULL;
  // Read in the unique file ID:
  if(readCEL_str(&file_id, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  if(verbose == 1) printf("file ID: %s\n", file_id);
  free(file_id);
  file_id = NULL;
  // Read in the date:
  result = readCEL_wstr(&date, f, bitflip);
  free(date);
  date = NULL;
  if(result != CEL_READ_VALUE_OK) return 1;
  // Read in the locale:
  result = readCEL_wstr(&locale, f, 1);
  if(verbose == 1) printf("locale: \"%s\"\n", locale);
  free(locale);
  locale = NULL;
  if(result != CEL_READ_VALUE_OK) return 1;
  // Read in the number of parameters:
  if(readCEL_int32(&parameter_number, 1, f, bitflip) != CEL_READ_VALUE_OK) return 1;
  if(verbose == 1) printf("parameter count: %d\n", parameter_number);
  // Read in each parameter in turn, and extract the data we need:
  for(i=0; i<parameter_number; i++){
    readCELcalvin_parameter(&parameter, f, bitflip);
    if(verbose == 1) printCELcalvin_parameter(&parameter);
    if(strcmp(parameter.name, "affymetrix-array-type") == 0) decode_CELcalvin_parameter_plaintext(&parameter, &(d->array));
    else if(strcmp(parameter.name, "affymetrix-algorithm-param-CellIntensityCalculationType") == 0){
      decode_CELcalvin_parameter_plaintext(&parameter, &(d->algorithm));
      for(j=0; j<strlen(d->algorithm); j++) d->algorithm[j] = tolower(d->algorithm[j]);
    }
    if(strcmp(parameter.name, "affymetrix-cel-rows") == 0) d->rows = decode_CELcalvin_parameter_int32(&parameter, bitflip);
    if(strcmp(parameter.name, "affymetrix-cel-cols") == 0) d->cols = decode_CELcalvin_parameter_int32(&parameter, bitflip);
    if(strcmp(parameter.name, "affymetrix-algorithm-param-CellMargin") == 0) d->cell_margin = decode_CELcalvin_parameter_int32(&parameter, bitflip);
    freeCELcalvin_parameter(&parameter);
  }
  //  Read in the single data group:
  fseek(f.handle, first_group_offset, SEEK_SET);
  CELcalvin_datagroup data_group;
  readCELcalvin_datagroup(&data_group, f, bitflip);
  if(verbose == 1) printf("first data group \"%s\" contains %d datasets:\n", data_group.name, data_group.dataset_number);
  // Go to the start of the first data set:
  fseek(f.handle, data_group.first_dataset_pos, SEEK_SET);
  CELcalvin_dataset data_set;
  for(i=0; i<data_group.dataset_number; i++){
    readCELcalvin_dataset(&data_set, f, bitflip);
    if(verbose == 1) printf(" dataset [%d] \"%s\" contains %d parameter(s), %d column(s) and %d row(s)\n", i, data_set.name, data_set.parameter_number, data_set.column_number, data_set.row_number);
    if(strcmp(data_set.name, "Outlier") == 0) d->outliers = data_set.row_number;
    if(strcmp(data_set.name, "Mask") == 0) d->masked = data_set.row_number;
    if((strcmp(data_set.name, "Intensity") == 0) && (read_intensity == 1)){
      fseek(f.handle, data_set.first_element_pos, SEEK_SET);
      intensities = (float*)malloc(data_set.row_number * sizeof(float));
      if(intensities == NULL){
        free_CELcalvin_dataset(&data_set);
        free_CELcalvin_datagroup(&data_group);
        return 1;
      }
      result = readCEL_float(intensities, data_set.row_number, f, bitflip);
      if(result != CEL_READ_VALUE_OK){
        free(intensities);
        intensities = NULL;
        free_CELcalvin_dataset(&data_set);
        free_CELcalvin_datagroup(&data_group);
        return 1;
      }
      calculate_intensity_stats(intensities, data_set.row_number, &d->intensity_max, &d->intensity_min, &d->intensity_n_unique, &d->intensity_n_invalid);
      d->intensity_stats_calculated = 1;
      free(intensities);
      intensities = NULL;
    }
    fseek(f.handle, data_set.next_dataset_pos, SEEK_SET);
    free_CELcalvin_dataset(&data_set);
  }
  free_CELcalvin_datagroup(&data_group);
  // No issues, so this must be valid:
  d->type = CEL_TYPE_CALVIN;
  d->valid = 1;
  return 0;
}
