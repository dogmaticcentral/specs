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

int afa_out (Options * options, Alignment * alignment) {

    int i, pos;
    char * seq;
    FILE * fptr;
    char filename[BUFFLEN];
    
    sprintf (filename, "%s.patched.afa", options->outname);

    fptr = efopen (filename, "w");
    if (!fptr) return 1;
    
    for (i=0; i<alignment->number_of_seqs; i++) {
	
	fprintf ( fptr,  ">%s\n", alignment->name[i]);
	seq = alignment->sequence[i];
	
	for (pos=0; pos<alignment->length; pos++) {
	    if ( seq [pos] == '.' ) {
		/* seaview doesn't like dots */
		fprintf ( fptr,  "-");
	    } else {
		fprintf ( fptr,  "%c",  seq [pos]);
	    }
	    if ( pos%50 == 49 ) fprintf (fptr, "\n");
	}
	if ( pos%50 )  fprintf (fptr, "\n");
  
    }
    fclose (fptr);
    
    return 0;
}
