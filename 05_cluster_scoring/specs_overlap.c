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

int epitope_detection_specificity ( Protein *protein,  int * res_rank, int * int_cvg,
				    int  surface, int *epitope, 
				    double * sensitivity,  double * specificity) {

    int pos, ctr, epitope_size = 0, epitope_complement_size = 0;
    int L = protein->length;
    
    if ( ! epitope ) return 0;
    
    for (pos = 0; pos < L; pos ++ ) {
	if ( surface &&  ! protein->sequence[pos].solvent_accessible ) continue;
	epitope_size += epitope[pos];
	epitope_complement_size += 1-epitope[pos];
    }
    if ( ! epitope_size || !epitope_complement_size ) {
	fprintf (stderr, "Epitope defintion error.\n");
	exit (1);	    
    }
    
    for (ctr=0; ctr < L  && int_cvg[ctr]; ctr++ ) { 
	sensitivity[ctr] = 0; 
	for (pos = 0; pos < L; pos ++ ) {
	    if ( surface &&  ! protein->sequence[pos].solvent_accessible ) continue;
	    if  (  epitope[pos]  ) {
		if ( res_rank[pos] <= int_cvg[ctr] ) 	sensitivity[ctr]++;
	    } else {
		if ( res_rank[pos] >  int_cvg[ctr] ) 	specificity[ctr]++;
	    }
	}
	sensitivity[ctr] /= epitope_size;
	specificity[ctr] /= epitope_complement_size;
    }

    return 0;
}
