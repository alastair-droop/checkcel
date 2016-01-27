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

#ifndef __checkcel_cel_calvin_h
#define __checkcel_cel_calvin_h

//Define the calvin parameter MIME types:
#define CEL_CALVIN_MIMETYPE_UNKNOWN 0
#define CEL_CALVIN_MIMETYPE_INT8 1
#define CEL_CALVIN_MIMETYPE_UINT8 2
#define CEL_CALVIN_MIMETYPE_INT16 3
#define CEL_CALVIN_MIMETYPE_UINT16 4
#define CEL_CALVIN_MIMETYPE_INT32 5
#define CEL_CALVIN_MIMETYPE_UINT32 6
#define CEL_CALVIN_MIMETYPE_FLOAT 7
#define CEL_CALVIN_MIMETYPE_PLAINTEXT 8
#define CEL_CALVIN_MIMETYPE_ASCII 9

// Structure to hold Calvin parameter object:
typedef struct {
  char* name;
  char* value;
  int value_length;
  char type;
} CELcalvin_parameter;

char readCELcalvin_parameter(CELcalvin_parameter *p, CELfile f, char bitflip);
void freeCELcalvin_parameter(CELcalvin_parameter *p);
void printCELcalvin_parameter(CELcalvin_parameter *p);

// Structure to hold a Calvin data group object:
typedef struct {
  u_int32_t next_pos;
  u_int32_t first_dataset_pos;
  int32_t dataset_number;
  char* name;
} CELcalvin_datagroup;

void free_CELcalvin_datagroup(CELcalvin_datagroup *g);
char readCELcalvin_datagroup(CELcalvin_datagroup *g, CELfile f, char bitflip);

// Structure to hold a Calvin dataset column:
typedef struct {
  char* name;
  int8_t type;
  int32_t size;
} CELcalvin_dataset_column;

char readCELcalvin_dataset_column(CELcalvin_dataset_column *c, CELfile f, char bitflip);
void free_CELcalvin_dataset_column(CELcalvin_dataset_column *c);

// Structure to hold Calvin dataset objects:  
typedef struct {
  u_int32_t first_element_pos;
  u_int32_t next_dataset_pos;
  char* name;
  int32_t parameter_number;
  CELcalvin_parameter *parameters;
  u_int32_t column_number;
  CELcalvin_dataset_column *columns;
  u_int32_t row_number;
} CELcalvin_dataset;

char readCELcalvin_dataset(CELcalvin_dataset *g, CELfile f, char bitflip);
void free_CELcalvin_dataset(CELcalvin_dataset *g);


// Functions to decode the data stored in a calvin parameter object:
int8_t decode_CELcalvin_parameter_int8(CELcalvin_parameter *p);
u_int8_t decode_CELcalvin_parameter_uint8(CELcalvin_parameter *p);
int16_t decode_CELcalvin_parameter_int16(CELcalvin_parameter *p, char bitflip);
u_int16_t decode_CELcalvin_parameter_uint16(CELcalvin_parameter *p, char bitflip);
int32_t decode_CELcalvin_parameter_int32(CELcalvin_parameter *p, char bitflip);
u_int32_t decode_CELcalvin_parameter_uint32(CELcalvin_parameter *p, char bitflip);
float decode_CELcalvin_parameter_float(CELcalvin_parameter *p, char bitflip);
void decode_CELcalvin_parameter_plaintext(CELcalvin_parameter *p, char **s);

char is_CELcalvin(CELfile f);
char readCELcalvin(CELfile f, CELdata *d, char read_intensity, char verbose);

#endif
