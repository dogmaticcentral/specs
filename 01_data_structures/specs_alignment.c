# include "specs.h"

/************************************************************/

int process_almt (Options *options, Alignment *alignment) {
    
    int retval;
    int count_gaps (Options *options, Alignment * alignment);
    int seq_pw_dist(Alignment * alignment);
    /* gaps */
    count_gaps (options, alignment);

    alignment->seq_dist = NULL;

    /*allocate space for various indicators of sequence similarity*/
    alignment->seq_dist =
	dmatrix ( alignment->number_of_seqs, alignment->number_of_seqs);
    
    alignment->aligned_sites =
	intmatrix ( alignment->number_of_seqs, alignment->number_of_seqs);
    if ( ! alignment->aligned_sites ) return 1;
    
    alignment->identical_sites =
	intmatrix ( alignment->number_of_seqs, alignment->number_of_seqs);
    if ( ! alignment->identical_sites ) return 1;
    
    alignment->similar_sites =
	intmatrix ( alignment->number_of_seqs, alignment->number_of_seqs);
    if ( ! alignment->similar_sites ) return 1;
    
 
 

    
    retval   = seq_pw_dist (alignment);
    if ( retval) return retval;

  
    return 0;
}

/*****************************************************************/
int count_gaps (Options *options, Alignment * alignment) {

    int s, c;
    alignment->seq_gaps    = (int *) emalloc (alignment->number_of_seqs*sizeof(int));
    if (!alignment->seq_gaps) return 1;
    alignment->column_gaps = (int *) emalloc (alignment->length*sizeof(int));
    if (!alignment->column_gaps) return 1;
    
    for ( s=0; s<alignment->number_of_seqs; s++ ) {
	for ( c=0; c<alignment->length; c++) {
	    if ( alignment->sequence[s][c] == '.' ) {
		alignment->column_gaps[c] ++;
		alignment->seq_gaps[s] ++;
	    }
	}
    }
    alignment->refseq_gaps = 0;
    if (options->no_refseqs) {
	int pos = 0;
	for ( c=0; c<alignment->length; c++) {
	    if ( alignment->refseq[0][c] == '.' ) {
		alignment->refseq_gaps ++;
	    } else {
		pos ++;
	    
	    }
	}
    }
    
    return 0;
}
