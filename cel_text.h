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

#ifndef __checkcel_cel_text_h
#define __checkcel_cel_text_h

//Define the maximum line size:
#define CEL_TEXT_MAX_LINE 10000
#define CEL_TEXT_HEADER_MAX 250

// Define the different line types:
#define CEL_TEXT_UNKNOWN_LINE 0
#define CEL_TEXT_HEADER_LINE 1
#define CEL_TEXT_TAG_LINE 2
#define CEL_TEXT_FAILED 3

typedef struct {
  char line[CEL_TEXT_MAX_LINE + 1];
  char line_type;
  char section[CEL_TEXT_HEADER_MAX + 1];
  char tag[CEL_TEXT_HEADER_MAX + 1];
  char data[CEL_TEXT_MAX_LINE + 1];
} CELtext_current_state;

void readCELtext_line(CELfile f, CELtext_current_state *state);

char is_CELtext(CELfile f);
char readCELtext(CELfile f, CELdata *d, char read_intensity, char verbose);

#endif
