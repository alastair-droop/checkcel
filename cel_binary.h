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

#ifndef __checkcel_cel_binary_h
#define __checkcel_cel_binary_h

//Define the structure holding the spot-level data:
#pragma pack(1)
typedef struct {
  float intensity;
  float sd;
  int16_t pixels;
} CELbinary_spotdata;
#pragma pack()

char is_CELbinary(CELfile f);
char readCELbinary(CELfile f, CELdata *d, char read_intensity, char verbose);

#endif
