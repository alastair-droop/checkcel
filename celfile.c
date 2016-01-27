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

char check_endian(){
  int x=1;
  if(*(char*)&x == 1) return MACHINE_LITTLE_ENDIAN;
  return MACHINE_BIG_ENDIAN;
}

CELfile open_CELfile(char* path){
  CELfile f;
  // Set default values for the structure:
  f.open = 0;
  f.path = NULL;
  f.name = NULL;
  // Allocate memory for the full path:
  f.path = malloc((strlen(path) + 1) * sizeof(char));
  if(f.path == NULL){
    f.open = 0;
    f.path = NULL;
    return f;
  }
  // Copy the full path to the object:
  strcpy(f.path, path);
  // Create a pointer to the file name part of the full file path:
  f.name = strrchr(f.path, '/');
  if(f.name == NULL) f.name = f.path;
  else f.name ++;
  //  Attempt to open the file:
  f.handle = fopen(f.path, "r");
  if(f.handle == NULL){
    fclose(f.handle);
    f.open = 0;
    return f;
  }
  // Set the file status to open:
  f.open = 1;
  return f;
}

void close_CELfile(CELfile f){
  if(f.path != NULL){
    free(f.path);
    fclose(f.handle);
  }
  f.open = 0;
  f.path = NULL;
  f.name = NULL;
}

void reset_CELfile(CELfile f){
  if(f.open == 1) rewind(f.handle);
}

char readCEL_int8(int8_t *value, size_t n, CELfile f){
  if(fread(value, sizeof(int8_t), n, f.handle) == n) return CEL_READ_VALUE_OK;
  return CEL_READ_VALUE_FAILED;
}

char readCEL_int16(int16_t *value, size_t n, CELfile f, char bitflip){
  unsigned char *buffer;
  int i;
  unsigned int a, b;
  buffer = (unsigned char*)malloc(n * 2 * sizeof(unsigned char));
  if(buffer == NULL) return CEL_READ_VALUE_FAILED;
  if(fread(buffer, sizeof(char), n * 2, f.handle) != n * 2) return CEL_READ_VALUE_FAILED;
  for(i=0; i<n; i++){
    a = (unsigned int)buffer[0 + (i * 2)];
    b = (unsigned int)buffer[1 + (i * 2)];
    if(bitflip == 0) value[i] = a | (b << 8);
    else value[i] = b | (a << 8);
  }
  free(buffer);
  return CEL_READ_VALUE_OK;
}

char readCEL_int32(int32_t *value, size_t n, CELfile f, char bitflip){
  unsigned char *buffer;
  int i;
  unsigned int a, b, c, d;
  buffer = (unsigned char*)malloc(n * 4 * sizeof(unsigned char));
  
  if(buffer == NULL) return CEL_READ_VALUE_FAILED;
  if(fread(buffer, sizeof(char), n * 4, f.handle) != n * 4) return CEL_READ_VALUE_FAILED;
  for(i=0; i<n; i++){
    a = (unsigned int)buffer[0 + (i * 4)];
    b = (unsigned int)buffer[1 + (i * 4)];
    c = (unsigned int)buffer[2 + (i * 4)];
    d = (unsigned int)buffer[3 + (i * 4)];
    if(bitflip == 0) value[i] = a | (b << 8) | (c << 16) | (d << 24);
    else value[i] = d | (c << 8) | (b << 16) | (a << 24);
  }
  free(buffer);
  return CEL_READ_VALUE_OK;
}

char readCEL_uint8(u_int8_t *value, size_t n, CELfile f){
  if(fread(value, sizeof(u_int8_t), n, f.handle) == n) return CEL_READ_VALUE_OK;
  return CEL_READ_VALUE_FAILED;
}

char readCEL_uint16(u_int16_t *value, size_t n, CELfile f, char bitflip){
  unsigned char *buffer;
  int i;
  unsigned int a, b;
  buffer = (unsigned char*)malloc(n * 2 * sizeof(unsigned char));
  if(buffer == NULL) return CEL_READ_VALUE_FAILED;
  if(fread(buffer, sizeof(char), n * 2, f.handle) != n * 2) return CEL_READ_VALUE_FAILED;
  for(i=0; i<n; i++){
    a = (unsigned int)buffer[0 + (i * 2)];
    b = (unsigned int)buffer[1 + (i * 2)];
    if(bitflip == 0) value[i] = a | (b << 8);
    else value[i] = b | (a << 8);
  }
  free(buffer);
  return CEL_READ_VALUE_OK;  
}

char readCEL_uint32(u_int32_t *value, size_t n, CELfile f, char bitflip){
  unsigned char *buffer;
  int i;
  unsigned int a, b, c, d;
  buffer = (unsigned char*)malloc(n * 4 * sizeof(unsigned char));
  if(buffer == NULL) return CEL_READ_VALUE_FAILED;
  if(fread(buffer, sizeof(char), n * 4, f.handle) != n * 4) return CEL_READ_VALUE_FAILED;
  for(i=0; i<n; i++){
    a = (unsigned int)buffer[0 + (i * 4)];
    b = (unsigned int)buffer[1 + (i * 4)];
    c = (unsigned int)buffer[2 + (i * 4)];
    d = (unsigned int)buffer[3 + (i * 4)];
    if(bitflip == 0) value[i] = a | (b << 8) | (c << 16) | (d << 24);
    else value[i] = d | (c << 8) | (b << 16) | (a << 24);
  }
  free(buffer);
  return CEL_READ_VALUE_OK;
}

char readCEL_float(float *value, size_t n, CELfile f, char bitflip){
  return readCEL_int32((int32_t *) value, n, f, bitflip);
}

char readCEL_char(char **c, size_t n, CELfile f){
  *c = (char*)malloc((n + 1) * sizeof(char));
  if(*c == NULL) return CEL_READ_VALUE_FAILED;
  memset(*c, 0, (n + 1) * sizeof(char));
  fread(*c, sizeof(char), n, f.handle);
  (*c)[n] = '\0';
  return CEL_READ_VALUE_OK;
}

char readCEL_str(char **s, CELfile f, char bitflip){
  int32_t string_length;
  if(readCEL_int32(&string_length, 1, f, bitflip) != CEL_READ_VALUE_OK) return CEL_READ_VALUE_FAILED;
  if(readCEL_char(s, string_length, f) == CEL_READ_VALUE_FAILED){
    free(*s);
    *s = NULL;
    return CEL_READ_VALUE_FAILED;
  }
  return CEL_READ_VALUE_OK;
}

char readCEL_wstr(char **s, CELfile f, char bitflip){
  int32_t i, string_length;
  char *buffer;
  if(readCEL_int32(&string_length, 1, f, bitflip) != CEL_READ_VALUE_OK) return CEL_READ_VALUE_FAILED;
  if(readCEL_char(&buffer, string_length * 2, f) == CEL_READ_VALUE_FAILED){
    free(buffer);
    buffer = NULL;
    return CEL_READ_VALUE_FAILED;
  }
  *s = (char*)malloc((string_length + 1) * sizeof(char));
  if(*s == NULL){
    free(buffer);
    buffer = NULL;
    return CEL_READ_VALUE_FAILED;    
  }
  memset(*s, 0, (string_length + 1) * sizeof(char));
  for(i=0; i < string_length; i++)(*s)[i] = buffer[1 + (i * 2)];
  free(buffer);
  buffer = NULL;
  return CEL_READ_VALUE_OK;
}

char check_CELtype(CELfile f){
  if(is_CELcalvin(f) == 1) return CEL_TYPE_CALVIN;
  if(is_CELbinary(f) == 1) return CEL_TYPE_BINARY;
  if(is_CELtext(f) == 1) return CEL_TYPE_TEXT;
  return CEL_TYPE_UNKNOWN;
}

char readCEL(CELfile f, CELdata *d, char read_intensity, char verbose){
  char type;
  type = check_CELtype(f);
  init_CELdata(d);
  if(type == CEL_TYPE_CALVIN) return readCELcalvin(f, d, read_intensity, verbose);
  if(type == CEL_TYPE_BINARY) return readCELbinary(f, d, read_intensity, verbose);
  if(type == CEL_TYPE_TEXT) return readCELtext(f, d, read_intensity, verbose);
  return 1;
}

