/*
This source code is part of specs, an application for
protein residue conservation scoring.
Written by Ivana Mihalek.
Copyright (C) 2007-2013 Ivana Mihalek.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see<http://www.gnu.org/licenses/>.

Contact: ivana.mihalek@gmail.com.
*/
# include "specs.h"

int read_epitope (Options * options, Protein * protein, int *epitope){

    char buf [LONGSTRING];
    char aux [PDB_ATOM_RES_NO_LEN+2];
    char res_type;
    int  found, total, res_ctr;
    FILE * fptr;
   
    fptr   = efopen ( options->epi_file_name, "r" );
    if (! fptr ) return 1;
    
    memset ( epitope, 0, protein->length*sizeof(int));
    total = 0;
    while ( ! feof ( fptr) ) {
 	if ( fgets (buf, LONGSTRING, fptr ) && strlen(buf) > 1 ) {
	    if ( buf[0] == '#' ) continue;
	    if (! (sscanf (buf, "%s %c", aux, &res_type)  ) ) {
		fprintf (stderr, "Error reading %s.\n", options->epi_file_name);
		fclose (fptr);
		return 1;
	    }
	    found = 0;
	    for ( res_ctr=0; res_ctr < protein->length; res_ctr++) {
		if ( ! strcmp ( protein->sequence[res_ctr].pdb_id, aux) ) {
		    if ( res_type !=  protein->sequence[res_ctr].res_type_short ) {
			fprintf (stderr, "While reading %s: Type mismatch for res id %s (pdb:%c, here:%c) .\n",
				 options->epi_file_name, aux, protein->sequence[res_ctr].res_type_short, res_type);
			fclose (fptr);
			return 1;
		    }
		    found = 1;
		    epitope[res_ctr]= 1;
		    total++;
		}
	    }
	    if ( ! found ) {
		fprintf (stderr, "While reading %s: Residue with the id %s not found.\n",
			 options->epi_file_name, aux);
		fclose (fptr);
		return 1;
	    }
	}
    }
    
    if ( !total ) {
	fprintf (stderr, "No positions matching pdb found in %s.\n", options->epi_file_name);
	return 1;
    }

    fclose (fptr);
    return 0;
}
