# include "specs.h"

int determinants (double sig_cutoff, Alignment * alignment,  Tree * tree,
		  int queried_node_id, double *score)
{
    /* root and leaf should correspond to the tree built on the alignment */ 
    int pos, node_ctr, no_seqs = alignment->number_of_seqs;
    int kstr;
    int int_max_gaps = MAX_GAPS*no_seqs;
    int retval;
    double df, chsq, prob;
    Node  *first_inner_node, *root, *leaf,  *node;
    int chisquare ( int * pool_population, int * node_population, int no_bins,
		    int kstr, double *df, double *chsq, double *prob);

    memset (score, 0, alignment->length*sizeof(double) );

    root = tree->root;
    leaf = tree->leaf;    
    first_inner_node = leaf + no_seqs;

    for ( pos=0; pos < alignment->length; pos++) {
	
	//printf ("pos %4d \n", pos+1); fflush ( stdout);
	score [pos] = 1.0;
	/* skip the column if gapped  */
	if ( (double)	alignment->column_gaps[pos] > int_max_gaps ) continue;
	
	/* set the population arrays */
	entropy_recursive (root,  pos, NULL, 0);
	
	/* look for the queried node: */ 
	for ( node_ctr= 0; node_ctr < no_seqs-1; node_ctr++ ) {
	    node = first_inner_node + node_ctr;
	    
	    if ( node->id != queried_node_id ) continue;
	    if ( node->number_of_leaves < MIN_NO_LEAVES ) continue;
	    if ( node->number_of_leaves == no_seqs ) continue;

	    /* difference btw the distribution of aa types in the initial pool,
	       and in the subtree */
	    retval =  chisquare (root->population, node->population,
				 ASCII, kstr=0, &df, &chsq, &prob);
	    if ( retval ) { /* error handling */
		return 1;
	    }
	    score [pos] = (prob < sig_cutoff) ? prob : 1.0;
  	}
	printf ( "** %4d  %8.3le \n", pos, prob);
    }

    return 0;
    
}

/******************************************************************/
/* find dterminants at each bifurcation of the tree               */
int path_determinants (double sig_cutoff, Alignment * alignment,  Tree * tree,
		       Node ** path_to_refseq,
		       double  **per_node_prob, int **description)
{
    Node * refseq_leaf;
    Node * node, *left, *right;
    int retval;
    int path_ctr, pos, ctr;
    int kstr;
    double df, chsq, prob;
    int chisquare ( int * pool_population, int * node_population, int no_bins,
		    int kstr, double *df, double *chsq, double *prob);
     
    /* find path to the refseq */
    /* find refseq, to begin with */
    refseq_leaf = tree->leaf;
    while ( refseq_leaf < tree->leaf + alignment->number_of_seqs &&
	    strcmp (alignment-> refseq_name[0], refseq_leaf ->name) ) refseq_leaf++;

    if ( refseq_leaf == tree->leaf + alignment->number_of_seqs) {
	fprintf ( stderr, "Refseq %s not found in the tree\n",
		  alignment-> refseq_name[0]);
	// return 1;
    } else {
	printf ( "Refseq %s found at node %d.\n", refseq_leaf->name, refseq_leaf->id);
    }

    for (path_ctr=0; path_ctr <alignment->number_of_seqs; path_ctr++) {
	path_to_refseq[path_ctr] = NULL;
    }
    
    node = refseq_leaf;
    path_ctr  = 0;
    do {
	path_to_refseq[path_ctr] = node;
	path_ctr ++;
    } while ( node->type != ROOT  && (node = node->parent));

    
    /* for each position, for each node on the path */
    for (pos = 0; pos<alignment->length; pos++) {
	/* set the population arrays */
	entropy_recursive (tree->root,  pos, NULL, 0);
	
	path_ctr  = 1;
	while (path_to_refseq[path_ctr]) {
	    /* compare the distributions on the on-path and off-path nodes */
	    
	    left = path_to_refseq[path_ctr]->left;
	    right =  path_to_refseq[path_ctr]->right;
	    retval =  chisquare (left->population, right->population,
				 ASCII, kstr=0, &df, &chsq, &prob);
	    if ( retval ) { /* error handling */
		return 1;
	    }
	    /* if there is a statistical difference btw the distribution,
	       mark the position as gain-of-function, loss-of-function,
	       or as a discriminant */

	    per_node_prob[pos][path_ctr] = prob;
	    if ( prob > sig_cutoff ) {
		description[pos][path_ctr] = INS;
	    } else {
		/* if entropy is small left and right,
		   call it a discriminant */
		/* in the lack of better, call the  entropy
		   small if it is less than 30% of its max value */
		int no_types_left = 0;
		int no_types_right = 0;
		for (ctr=0; ctr < ASCII; ctr++ ) {
		    no_types_left  += (left->population[ctr]);
		    no_types_right += (right->population[ctr]);
		}
		if ( left->entropy < 0.3*log (no_types_left)  &&
		     right->entropy < 0.3*log (no_types_right)  ) {
		    description[pos][path_ctr] = DISCR;
		} else {
		    /* else, depending on which one is on the path
		       to the refseq, call it LOF or GOF */
		    if ( left == path_to_refseq [path_ctr-1] ) {
			if ( left->entropy < 0.3*log (no_types_left) ) {
			    description[pos][path_ctr] = GOF;
			} else {
			    description[pos][path_ctr] = LOF;
			}
		    } else {
			if ( left->entropy < 0.3*log (no_types_left) ) {
			    description[pos][path_ctr] = LOF;
			} else {
			    description[pos][path_ctr] = GOF;
			}
		    }
		}
	    }
	    path_ctr ++;
	}
    }
    

    return 0;
}
