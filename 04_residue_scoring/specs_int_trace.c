# include "specs.h"

int int_trace (Alignment * alignment, Tree * tree, int * similar_to,  int normalize, double *score){
    
    int pos, rank, group, done;
    int int_trace_recursive_part ( Node *node,  int * similar_to,  int normalize, int pos );
  
    /* calculate the score */
    for ( pos=0; pos < alignment->length; pos++) {
 	int_trace_recursive_part  ( tree->root, similar_to, normalize, pos);
 	score[pos] =  tree->no_of_leaves;
	done = 0;
	for ( rank = tree->no_of_leaves - 1; rank >=1  && !done; rank-- ) {
	    for ( group=0; group<rank && ! done; group++) {
		if ( ! tree->group_root[rank][group]->consensus ) {
		    score[pos] = tree->group_root[rank][group]-> id + 1;
		    done = 1;
		}
	    }
	}
	if ( ! done ) score[pos] = 1;
    }
    return 0;
}

/***********************************************************************/
int int_trace_recursive_part ( Node *node,  int * similar_to, int normalize, int pos ) {

    if ( node->type == LEAF ) {
	int aa;
	aa = node->seq[pos];
	if ( similar_to ) aa = similar_to[aa];
	node->consensus = aa;
	
    } else {
	int_trace_recursive_part (node->left, similar_to, normalize, pos);
	int_trace_recursive_part (node->right, similar_to, normalize, pos);

	if ( node->left->consensus ==  node->right->consensus ) {
	    node->consensus = node->left->consensus;
	    
	} else if (normalize) {
	    if  (node->left->consensus == '.') {
		node->consensus = node->right->consensus;
	    } else if (node->right->consensus == '.') {
		node->consensus = node->right->consensus;
	    } else {
		node->consensus = '\0';
	    }
	    
	} else {
	    node->consensus = '\0';
	}
    }

    return 0;
}

/***********************************************************************/
int carbone (Alignment * alignment, Tree * tree, int * similar_to,  int normalize, double *score){
    
    int pos, rank, group, done;
    int conserved_groups;
    int int_trace__recursive_part ( Node *node,  int * similar_to, int pos );
    /* calculate the score */
    for ( pos=0; pos < alignment->length; pos++) {
	
 	int_trace_recursive_part  ( tree->root, similar_to, normalize, pos);
	
	if ( tree->group_root[1][0]->consensus ) {
	    score[pos] = 1;
	    
	} else {
	    
	    score[pos] =  tree->no_of_leaves;
	    done = 0;
	    conserved_groups = 0;
	    
	    for ( rank = 2; rank < tree->no_of_leaves  && !done; rank++ ) {
		for ( group=0; group<rank && ! done; group++) {
		    if (  tree->group_root[rank][group]->consensus ) {
			conserved_groups ++;
		    }
		}
		if ( conserved_groups >= 2 ){
		    score[pos] = rank;
		    done = 1;
		}
	    }
	}
    }
    return 0;
}
