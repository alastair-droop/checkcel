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

void readCELtext_line(CELfile f, CELtext_current_state *state){
  char *result;
  int n;
 
  result = fgets(state->line, CEL_TEXT_MAX_LINE, f.handle);
  if(result == NULL){
    state->line_type = CEL_TEXT_FAILED;
    return;
  }
    
  //Read in the header if possible:
  n = sscanf(state->line, "[%[ABCDEFGHIJKLMNOPQRSTUVWXYZ]]", state->section);
  if(n == 1){
    state->line_type = CEL_TEXT_HEADER_LINE;
    return;
  }
  
  //Read in the tag-key data if possible:
  n = sscanf(state->line, "%[^=]=%s", state->tag, state->data);
  if(n == 2){
    state->line_type = CEL_TEXT_TAG_LINE;
    return;
  }
  
  //Non standard line:
  state->line_type = CEL_TEXT_UNKNOWN_LINE;
}

char is_CELtext(CELfile f){
  CELtext_current_state state;
  reset_CELfile(f);
  readCELtext_line(f, &state);
  if((state.line_type != CEL_TEXT_HEADER_LINE) || (strcmp(state.section, "CEL") != 0)) return 0;
  readCELtext_line(f, &state);
  if((state.line_type != CEL_TEXT_TAG_LINE) || (strcmp(state.tag, "Version") != 0) || (strcmp(state.data, "3") != 0)) return 0;
  // All initial checks look good, to reset the file and return success:
  reset_CELfile(f);
  return 1;
}

char readCELtext(CELfile f, CELdata *d, char read_intensity, char verbose){
  unsigned int i, x, y, intensity_number;
  float *intensities;
  char data_line[CEL_TEXT_MAX_LINE + 1];
  CELtext_current_state state;
  char *p, *result;
  d->valid = 0;
  while(1){
    readCELtext_line(f, &state);
    if(state.line_type == CEL_TEXT_FAILED) break;
    
    if((state.line_type == CEL_TEXT_TAG_LINE) && (strcmp(state.section, "HEADER") == 0)){
      // A tag line in the header section:s
      if(strcmp(state.tag, "Rows") == 0) sscanf(state.data, "%d", &d->rows);
      if(strcmp(state.tag, "Cols") == 0) sscanf(state.data, "%d", &d->cols);
      if(strcmp(state.tag, "Algorithm") == 0){
        d->algorithm = (char*)malloc((strlen(state.data) + 1) * sizeof(char));
        if(d->algorithm == NULL){
          free(d->algorithm);
          d->algorithm = NULL;
          return 1;
        }
        sscanf(state.data, "%s", d->algorithm);
        for(i=0; i<strlen(d->algorithm); i++) d->algorithm[i] = tolower(d->algorithm[i]);
      }
      if(strcmp(state.tag, "DatHeader") == 0) extract_chipname(state.line, d);
      if(strcmp(state.tag, "AlgorithmParameters") == 0){
        p = strstr(state.data, "CellMargin:");
        if(p != NULL){
          p += 11;
          sscanf(p, "%d", &d->cell_margin);
        }
      }
    }
    
    if((state.line_type == CEL_TEXT_HEADER_LINE) && (strcmp(state.section, "INTENSITY") == 0)){
      intensity_number = 0;
      readCELtext_line(f, &state);
      if(strcmp(state.tag, "NumberCells") != 0) return 1;
      sscanf(state.data, "%d", &intensity_number);
      if (read_intensity == 1) {
        intensities = (float*)malloc(intensity_number * sizeof(float));
        if(intensities == NULL) return 1;
        result = fgets(data_line, CEL_TEXT_MAX_LINE, f.handle);
        if(result == NULL){
          free(intensities);
          intensities = NULL;
          return 1;
        }
        for(i=0; i<intensity_number; i++){
          result = fgets(data_line, CEL_TEXT_MAX_LINE, f.handle);
          if(result == NULL){
            free(intensities);
            intensities = NULL;
            return 1;
          }
          sscanf(data_line, "%d%d%f", &x, &y, &intensities[i]);
        }
        calculate_intensity_stats(intensities, intensity_number, &d->intensity_max, &d->intensity_min, &d->intensity_n_unique, &d->intensity_n_invalid);
        d->intensity_stats_calculated = 1;
        free(intensities);
        intensities = NULL;
      } else {
        for(i=0; i<intensity_number; i++) if(fgets(data_line, CEL_TEXT_MAX_LINE, f.handle) == NULL) return 1;
      }
      continue;
    }

    if((state.line_type == CEL_TEXT_HEADER_LINE) && (strcmp(state.section, "MASKS") == 0)){
      readCELtext_line(f, &state);
      if(strcmp(state.tag, "NumberCells") != 0) return 1;
      sscanf(state.data, "%d", &d->masked);
      for(i=0; i<d->masked; i++) if(fgets(data_line, CEL_TEXT_MAX_LINE, f.handle) == NULL) return 1;
      continue;
    }
    
    if((state.line_type == CEL_TEXT_HEADER_LINE) && (strcmp(state.section, "OUTLIERS") == 0)){
      readCELtext_line(f, &state);
      if(strcmp(state.tag, "NumberCells") != 0) return 1;
      sscanf(state.data, "%d", &d->outliers);
      for(i=0; i<d->outliers; i++) if(fgets(data_line, CEL_TEXT_MAX_LINE, f.handle) == NULL) return 1;
      continue;
    }
  }  
  // No issues, so this must be valid:
  d->type = CEL_TYPE_TEXT;
  d->valid = 1;
  return 0;
}
