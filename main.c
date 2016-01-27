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
#include <getopt.h>
#include <glob.h>
#include <string.h>
#include "cel.h"

void print_usage(){
  printf("usage: checkcel [-cfvh] file [...]\n");
}

void print_version(){
  printf("checkcel 1.5.0 (2016-01-27)\n");
}

int main (int argc, const char * argv[])
{
  int i, j, option;
  char check_spot_level_data, filter_bad_files;
  glob_t glob_data;
  CELfile f;
  CELdata cel_data;

  // Sort out the command line options:
  check_spot_level_data = 0;
  filter_bad_files = 0;
  while ((option = getopt(argc, (char* const*)argv, "cfvh")) != -1){
    switch (option){
      case 'c':
        check_spot_level_data = 1;
        break;
      case 'f':
        filter_bad_files = 1;
        break;
      case 'v':
        print_version();
        return 0;
      case 'h':
        print_usage();
        printf("Options:\n");
        printf("-c: calculate and display intensity statistics\n");
        printf("-f: filter out invalid CEL files\n");
        printf("-h: display this help information\n");
        printf("-h: display version\n");
        printf("\nOutput columns:\n");
        printf("  : CEL file name\n");
        printf("  : file format\n");
        printf("  : chip ID\n");
        printf("  : creation algorithm\n");
        printf("  : row count\n");
        printf("  : column count\n");
        printf("  : margin\n");
        printf("  : outier cell count\n");
        printf("  : masked cell count\n");
        printf("\nIntensity statistics:\n");
        printf("  : minimum intensity value\n");
        printf("  : maximum intensity value\n");
        printf("  : unique value count\n");
        printf("  : invalid value count\n");
        return 0;
      default:
        print_usage();
        return 1;
    }
  }

  // Loop over the remaining command line arguments:
  for(i=optind; i<argc; i++){
    //  Expand the wildcard file listing to get a list of valid files to process:
    glob(argv[i], 0, NULL, &glob_data);
    if(glob_data.gl_matchc < 1){
      printf("no matching file\n");
      return 1;
    }
    // Initialise the CEL file dataset:
    init_CELdata(&cel_data);
    //  Run through each file in turn, processing it:
    for(j=0; j<glob_data.gl_matchc; j++){
      f = open_CELfile(glob_data.gl_pathv[j]);
      readCEL(f, &cel_data, check_spot_level_data, 0);
      if(cel_data.valid == 1){
        printf("%s\t", f.name);
        print_CELdata(&cel_data);
      } else if(filter_bad_files != 1) printf("%s\tunknown\n", f.name);
      free_CELdata(&cel_data);
      close_CELfile(f);
    }
  }
  return 0;
}
