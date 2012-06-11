# include "specs.h"


/* res_cvg array maps each residue into cvg it belongs to */
/* int_cvg is an array of all coverages (expressed as integers)  implied by the score array  */

int coverage ( Protein * protein, Alignment * alignment, int * almt2prot,
	       double * score, int almt_length,
	       int * res_rank, int * int_cvg ,
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
	new_length = alignment->length - alignment->refseq_gaps;
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
	    
	} else if ( alignment->refseq &&  alignment->refseq[0][ctr] != '.' ){
	    
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

    /* turn the sorted array to coverage info */
    
    /* find the lowest score in the game */
    prev_score = 0;
    for (ctr=0; ctr < new_length; ctr++) {
	if ( surface &&  ! protein->sequence[sorted_res[ctr]].solvent_accessible ) continue; 
	prev_score = protein_score[ sorted_res[ctr] ];
	break; 
    }
    
    
    first   = 0;
    cvg_ctr = 0;
    int_cvg[cvg_ctr] = 0;
    surface_size     = 0;
    for (ctr=0; ctr < new_length; ctr++) {
	if ( surface &&  ! protein->sequence[sorted_res[ctr]].solvent_accessible ) continue;
	surface_size ++;
	if ( protein_score[ sorted_res[ctr] ] <= prev_score ) {
	    int_cvg[cvg_ctr] ++;
	} else {
	    prev_score  = protein_score[ sorted_res[ctr] ];
	    for (ctr2=first; ctr2 <ctr; ctr2++ ) {
		if ( surface &&  ! protein->sequence[sorted_res[ctr2]].solvent_accessible ) continue;
		res_rank[ sorted_res[ctr2] ] = int_cvg[cvg_ctr];
	    }
	    first = ctr;
	    cvg_ctr ++;
	    int_cvg[cvg_ctr] =  int_cvg[cvg_ctr-1] + 1;
	}
    }
    for (ctr2=first; ctr2 <ctr; ctr2++ ) {
	if ( surface &&  ! protein->sequence[sorted_res[ctr2]].solvent_accessible ) continue;
	res_rank[ sorted_res[ctr2] ] =  int_cvg[cvg_ctr];
    }
    
    
    /* sanity : */
    if (  int_cvg[cvg_ctr] != surface_size ) {
	fprintf (stderr, "Error: surf size %d, counted %d \n", surface_size, int_cvg[cvg_ctr] );
	exit (1);
    }

    /* free */
    free (protein_score);
    free (sorted_res);

    
    return 0;
}
