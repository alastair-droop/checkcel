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

#ifndef __checkcel_celfile_h
#define __checkcel_celfile_h

//Define the value read status:
#define CEL_READ_VALUE_OK 0
#define CEL_READ_VALUE_FAILED 1

//Define the CEL file types:
#define CEL_TYPE_UNKNOWN 100
#define CEL_TYPE_BINARY 101
#define CEL_TYPE_CALVIN 102
#define CEL_TYPE_TEXT 103

//A function to check the endianness of the current machine:
#define MACHINE_LITTLE_ENDIAN 0
#define MACHINE_BIG_ENDIAN 1
char check_endian();

// Define the struct to hold a CELfile connection:
typedef struct {
  char open;
  char type;
  char *path;
  char *name;
  FILE *handle;
} CELfile;

// Functions to manipulate the CELfile connection:
CELfile open_CELfile(char* path);
void close_CELfile(CELfile f);
void reset_CELfile(CELfile f);

// Functions to read signed integers from a CELfile:
char readCEL_int8(int8_t *value, size_t n, CELfile f);
char readCEL_int16(int16_t *value, size_t n, CELfile f, char bitflip);
char readCEL_int32(int32_t *value, size_t n, CELfile f, char bitflip);

// Functions to read unsigned integers from a CELfile:
char readCEL_uint8(u_int8_t *value, size_t n, CELfile f);
char readCEL_uint16(u_int16_t *value, size_t n, CELfile f, char bitflip);
char readCEL_uint32(u_int32_t *value, size_t n, CELfile f, char bitflip);

// Functions to read floating point values from a CELfile:
char readCEL_float(float *value, size_t n, CELfile f, char bitflip);

// Functions to read strings from a CELfile:
char readCEL_char(char **c, size_t n, CELfile f);
char readCEL_str(char **s, CELfile f, char bitflip);
char readCEL_wstr(char **s, CELfile f, char bitflip);

// Function to test which type a file is:
char check_CELtype(CELfile f);

// Function to open an arbitrary CEL file:
char readCEL(CELfile f, CELdata *d, char read_intensity, char verbose);

#endif
