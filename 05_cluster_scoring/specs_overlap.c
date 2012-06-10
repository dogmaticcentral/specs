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
