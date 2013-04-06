# include "specs.h"

/* 'Normalization' here ferers to using only non-gapped postions to say something
   about the conservation. [Rather than tresting gaps as 21st amino acid type.
   Apparently Valdar is 'nomralized' by definition. */

int  scoring ( Options *options, Alignment * alignment, Tree * tree,
	       int * similar_to, double * score, int score_ctr) {
    
    int  carbone (Alignment * alignment, Tree * tree, int * similar_to,
		  int normalize, double *score);
    int  entropy ( Alignment * alignment, int * similar_to,  int normalize, double *score);
    int  hybrid (Alignment * alignment, Tree * tree,
		 int * similar_to, int normalize,  double *score);
    int  int_trace (Alignment * alignment, Tree * tree, int * similar_to,
		    int normalize,  double *score);
    int  majority_fraction ( Alignment * alignment, int * similar_to,
			     int normalize, double *score);
    int  pheno_init (char * table_file);
    int  pheno (Alignment * alignment, double *score, int blosum_ctr);
    int  hybrid_pheno (Alignment * alignment, double *score, Tree * tree);
    int  pheno_over_nodes (Alignment * alignment, double *score, Tree * tree);
    int  valdar ( Alignment * alignment, double *score);
    
    static int table_init =0;
    int normalize = 1;

    
    if (!table_init  && 
	    ( options->method[score_ctr] == PHENO ||
	      options->method[score_ctr] == PHENO_HYBRID ||
	      options->method[score_ctr] == PHENO_NODES)){
	if ( pheno_init(PROBS) ) exit (1);
	table_init = 1;
    }
    switch  ( options->method[score_ctr]) {
    case VALDAR:
	valdar ( alignment, score);
	break;
    case ENTROPY:
	entropy ( alignment, NULL,  normalize, score);
	break;
    case ENTR_W_SIM:
	entropy ( alignment, similar_to, normalize, score);
	break;
    case IVET:
	int_trace ( alignment, tree, NULL, normalize, score);
	break;
    case CARBONE:
	carbone ( alignment, tree, NULL, normalize, score);
	break;
    case IV_W_SIM:
	int_trace ( alignment, tree, similar_to, normalize, score);
	break;
    case MAJORITY:
	majority_fraction ( alignment, NULL, normalize, score);
	break;
    case MAJORITY_W_SIM:
	majority_fraction ( alignment, similar_to, normalize, score);
	break;
    case RVET:
	hybrid ( alignment, tree, NULL, 0,  score);
	break;
    case RVET_NORM:
	hybrid ( alignment, tree, NULL, 1,  score);
	break;
    case RV_W_SIM:
	hybrid ( alignment, tree, similar_to, normalize,  score);
	break;
    case PHENO:
	pheno (alignment, score, 1);
	break;
    case PHENO_HYBRID:
	hybrid_pheno ( alignment, score, tree);
	break;
    case PHENO_NODES:
	pheno_over_nodes ( alignment, score, tree);
	break;
	
    }
    /* sink the gapped positions to the bottom, if requested */
    if ( options->max_gaps ) {
	int ctr;
	double max_score = -1;
	for ( ctr=0; ctr < alignment->length; ctr++ ) {
	    if ( max_score < score[ctr] ) max_score = score[ctr];
	}
	for ( ctr=0; ctr < alignment->length; ctr++ ) {
	    if ((double)alignment->column_gaps[ctr]/alignment->number_of_seqs> options->max_gaps ) {
		score[ctr] = max_score;
	    }
	}
    }

    return 0;
}
