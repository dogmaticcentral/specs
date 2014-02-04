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


/* res_cvg array maps each residue into cvg it belongs to */
/* int_cvg is an array of all coverages (expressed as integers)  implied by the score array  */

int coverage ( Protein * protein, Alignment * alignment, int * almt2prot,
	       double * score, int almt_length,
	       int * res_rank, int * number_of_positions_at_cvg ,
	       int restrict2structure, int surface) {

    double * protein_score, prev_score; /* "score" refers to alignment positions */
    int * sorted_res;
    int pos, ctr, ctr2;
    int first, cvg_ctr, surface_size;
    int new_length =0;
    
		
    restrict2structure = almt2prot && restrict2structure; /* paranoid programing ?*/
    
    if (restrict2structure) {
	new_length = protein->length;
    } else {
	new_length = alignment->number_of_protected_positions;
    }

    /*allocate */
    protein_score = (double *) emalloc (new_length*sizeof(double));
    if (!protein_score ) return 1;
    sorted_res    =    (int *) emalloc (new_length*sizeof(int));
    if (!sorted_res )    return 1;

    
    /* remove gapped positions from the score array */
    pos = 0;
    for (ctr=0; ctr < almt_length; ctr ++ ) {
	if ( restrict2structure ) {
	    pos =  almt2prot[ctr];
	    if ( pos >= 0 ) {
		protein_score[pos] = score[ctr];
	    }

	    /* the reference and pdb  sequences are protected */   
	} else if ( alignment->protected_position[ctr] ){
	    
	    if ( pos >= new_length ) {
		fprintf (stderr, "Error in %s:%d: pos ctr >= allocated length (%d) \n",
			 __FILE__, __LINE__, new_length);
		exit (1);
	    }
	    protein_score[pos] = score[ctr];
	    pos ++;
	}
    }
    
    
    /* sort protein residues according to the new array */
    for (pos=0; pos < new_length; pos++) sorted_res[pos] = pos;
    array_qsort ( sorted_res, protein_score, new_length);

     /* find the lowest score in the game */
    prev_score = 0;
    for (ctr=0; ctr < new_length; ctr++) {
	prev_score = protein_score[ sorted_res[ctr] ];
	break; 
    }
    
    
    /* turn the sorted array to coverage info */
    first   = 0;
    cvg_ctr = 0;
    number_of_positions_at_cvg[cvg_ctr] = 0;
    surface_size     = 0;
    for (ctr=0; ctr < new_length; ctr++) {
	surface_size ++;
	if ( protein_score[ sorted_res[ctr] ] <= prev_score ) {
	    number_of_positions_at_cvg[cvg_ctr] ++;
	} else {
	    prev_score  = protein_score[ sorted_res[ctr] ];
	    for (ctr2=first; ctr2 <ctr; ctr2++ ) {
		res_rank[ sorted_res[ctr2] ] = number_of_positions_at_cvg[cvg_ctr];
	    }
	    first = ctr;
	    cvg_ctr ++;
	    number_of_positions_at_cvg[cvg_ctr] =  number_of_positions_at_cvg[cvg_ctr-1] + 1;
	}
	
    }
    for (ctr2=first; ctr2 <ctr; ctr2++ ) {
	res_rank[ sorted_res[ctr2] ] =  number_of_positions_at_cvg[cvg_ctr];
    }
    
   
    /* sanity : */
    if (  number_of_positions_at_cvg[cvg_ctr] != surface_size ) {
	fprintf (stderr, "Error: surf size %d, counted %d \n", surface_size, number_of_positions_at_cvg[cvg_ctr] );
	exit (1);
    }

    /* free */
    free (protein_score);
    free (sorted_res);

    
    return 0;
}
