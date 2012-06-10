# include "specs.h"
int upgma_and_nj (Options * options, Alignment * alignment, Tree * tree) {

    double  distance1, distance2;
    int retval;
    int node_ctr, no_seqs = alignment->number_of_seqs, tree_size;
    int closest1, closest2,* current_list;
    int upper, last_ctr;
    Node * node;
    int ( *closest_in_curr_list ) (Alignment * alignment, Node * node, int tree_size, int new_node, 
				   int *  current_list, int * closest1, int *closest2,
				   double * dist_1_ptr, double * dist_2_ptr );
    int closest_in_curr_list_nj (Alignment *alignment, Node * node, int tree_size, int new_node, 
			     int *  current_list, int * closest1, int *closest2,
				 double * dist_1_ptr, double * dist_2_ptr );
    int closest_in_curr_list_consensus_upgma (Alignment *alignment, Node * node,
					  int tree_size, int dummy, 
					  int *  current_list, int * closest1,
					  int *closest2, double * dist_1_ptr,
					  double * dist_2_ptr );
    /* I need dummy here to make upgma and nj defs look formally the same */ 
    int closest_in_curr_list_upgma (Alignment *alignment, Node * node, int tree_size, int dummy, 
			     int *  current_list, int * closest1, int *closest2,
				    double * dist_1_ptr, double * dist_2_ptr );
    int fill_consensus_array (int seq_length, char * seq1, char * seq2, char *consensus);
    switch (options->tree_method) {
    case UPGMA:
	closest_in_curr_list = closest_in_curr_list_upgma;
	break;
    case CONSENSUS_UPGMA:
	closest_in_curr_list = closest_in_curr_list_consensus_upgma;
	break;
    case NEIGHBOR_JOINING:
	closest_in_curr_list = closest_in_curr_list_nj;
	break;
    default:
	fprintf ( stderr,"Unrecognized tree building method.\n");
	return 1;
    }

    /*****************************************/
    /*****************************************/
    /*          build tree                   */
    /*****************************************/
    /*****************************************/
    
    /* allocate space */
    tree_size = 2*no_seqs -1;
    node = (Node *) emalloc ( tree_size*sizeof(Node) );
    /* initialize leaves to point to sequences from the alignment*/
    for ( node_ctr=0; node_ctr < no_seqs; node_ctr++ ) {
	node[node_ctr].id   = node_ctr;
	node[node_ctr].type = LEAF;
	node[node_ctr].seq  = alignment->sequence[node_ctr];
        node[node_ctr].name = alignment->name[node_ctr];
	node[node_ctr].number_of_leaves = 1;
    }
    /* initialize current list of nodes whose distance needs to be compared */
    current_list = (int*) emalloc ( tree_size*sizeof(int));
    for ( node_ctr=0; node_ctr < no_seqs; node_ctr++ ) {
	current_list[node_ctr] = 1; 
    }

    /* find children for each of the remaining nodes */
    upper = (options->tree_method == NEIGHBOR_JOINING) ? tree_size - 1 : tree_size;
    for ( node_ctr=no_seqs; node_ctr < upper; node_ctr++ ) {
	retval = closest_in_curr_list (alignment, node, tree_size, node_ctr,
				       current_list, &closest1, &closest2,
				       &distance1, &distance2);
	if (retval) return retval;
	/* hack, so I can order the nodes in the tree: */
	if ( distance1 <= 0 ) distance1 = 0.001;
	if ( distance2 <= 0 ) distance2 = 0.001;
	
	/* fill in the new  node fields */ 
	node[node_ctr].left             = &node[closest1];
	node[node_ctr].dist_to_left     = distance1;
	node[node_ctr].right            = &node[closest2];
	node[node_ctr].dist_to_right    = distance2;
	node[node_ctr].type             = INNER;
	node[node_ctr].id               = tree_size - node_ctr; /* this will serve as the "rank" */
	node[node_ctr].number_of_leaves = node[closest1].number_of_leaves +
	    node[closest2].number_of_leaves;
	if ( ! 	node[node_ctr].seq ) {
	    node[node_ctr].seq   = emalloc ( alignment->length*sizeof(char) );
	    if ( ! node[node_ctr].seq ) exit (1);
	}
	node[node_ctr].consensus_length =
	    fill_consensus_array( alignment->length, node[closest1].seq,
				 node[closest2].seq, node[node_ctr].seq);
	    
	node[closest1].parent = & node[node_ctr];
	node[closest2].parent = & node[node_ctr];
	node[closest1].dist_to_parent = distance1;
	node[closest2].dist_to_parent = distance2;
	/* remove the two from the current list, and replace them with the parent */ 
	current_list[closest1] = 0;
	current_list[closest2] = 0;
	current_list[node_ctr] = 1;
    }

    last_ctr = node_ctr - 1;
    tree->leaf = node;
    (node + tree_size -1)->type = ROOT;
    tree->root = node + tree_size - 1;

    
    if( options->tree_method == NEIGHBOR_JOINING) {
	
	int place_root ( Node *leaf, int no_of_nodes, Node * root) ;
	int rank_inner_nodes (Node *root, int  no_of_nodes);
	/* make the last two nodes  parents of each other */
	/* place the root */
	for ( node_ctr=0; node_ctr < upper; node_ctr++ ) {
	    if ( current_list[node_ctr] ) printf ( " %d ", node_ctr);
	}
	printf ("\n");
	for ( node_ctr=0; node_ctr < upper; node_ctr++ ) {
	    if ( current_list[node_ctr] ) {
		node[node_ctr].parent =  & node[last_ctr];
		node[last_ctr].parent =  & node[node_ctr];
		break;
	    }
	}
	place_root ( tree->leaf,  tree_size-1, tree->root);
	rank_inner_nodes ( tree->root,  tree_size);
    }

    tree->size = tree_size;
    tree->no_of_leaves = no_seqs; 
# if 0
    print_tree (stdout, tree->root);
    exit(0);
# endif
    
    free (current_list);
   
    return 0;

}
/********************************************************************************/
/********************************************************************************/
int closest_in_curr_list_upgma (Alignment * alignment, Node * node, int tree_size, int dummy, 
				int *  current_list, int * closest1, int *closest2,
				double * dist_1_ptr, double * dist_2_ptr ){
    int ctr1, ctr2;
    double distance, min_distance;
    double **seq_dist = alignment->seq_dist;
    double  node_distance ( double **seq_dist, Node* node1, Node* node2 );

    min_distance = 1000;
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
        if ( ! current_list[ctr1]) continue;
	for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    distance = node_distance ( seq_dist,  node+ctr1, node+ctr2 );
	    //exit (0);
	    /* we need the average */ 
	    distance /= (node+ctr1)->number_of_leaves*(node+ctr2)->number_of_leaves;
	    if ( distance > 1 ) {
		printf ( "closest_in_curr_list_upgma(): %d   %d  %d   %d  %8.2lf\n",
			 ctr1, ctr2,  (node+ctr1)->number_of_leaves,
			 (node+ctr2)->number_of_leaves, distance);
		exit (0);
	    }
	    if ( distance < min_distance) {
		min_distance = distance;
		*closest1 = ctr1;
		*closest2 = ctr2;
	    }
	}
    }
    *dist_1_ptr = min_distance/2;
    *dist_2_ptr = min_distance/2;
    return 0;
}
/********************************************************************************/
/********************************************************/
int fill_consensus_array (int seq_length, char * seq1, char * seq2, char *consensus) {
    
    int i, ctr = 0;
    
    for (i=0; i<seq_length; i++ ) {
	if( seq1[i] && seq1[i] != '.' &&  seq1[i] ==  seq2[i]) {
	    consensus[i] = seq1[i];
	    ctr ++;
	} else {
	   consensus[i]   = '\0';
	}
    }

    return ctr;
}
