# include "specs.h"

/********************************************************************************/
int build_tree (Options * options, Alignment * alignment, Tree * tree) {
    
# define NUMBER_OF_BINS 20

    int rank, group_tag, retval;
    int upgma_and_nj (Options * options, Alignment * alignment, Tree * tree);
   
    switch (options->tree_method) {
    case UPGMA:
    case CONSENSUS_UPGMA:
    case NEIGHBOR_JOINING:
	retval =  upgma_and_nj (options, alignment, tree);
	 /* should fill the tree->size value*/
	if ( retval ) return retval;
	break;
    default:
	fprintf ( stderr,"Unrecognized tree building method.\n");
	return 1;
    }
    
    /* find groups (by similarity) */
    tree->group_root_sim = node_matrix (NUMBER_OF_BINS, tree->no_of_leaves);
    for ( rank = 1; rank < tree->no_of_leaves; rank++ ) {
	group_tag = -1;
	find_group_roots_sim (tree->root, tree->no_of_leaves-1, NUMBER_OF_BINS,  tree->group_root_sim);
    }
    
    /* find groups */
    tree->group_root = node_matrix (tree->no_of_leaves, tree->no_of_leaves);
    for ( rank = 1; rank < tree->no_of_leaves; rank++ ) {
	group_tag = -1;
	find_group_roots (tree->root,  tree->group_root,  rank, & group_tag );
    }

    return 0;
}


/********************************************************************************/
double  node_distance ( double **seq_dist, Node* node1, Node* node2 ){

    
    if (  (node1->type == LEAF) &&  (node2->type == LEAF) ) {
	return seq_dist[node1->id][node2->id];
    } else {
	double distance = 0.0;
	if ( node1->type != LEAF ) { /* left side depth first */
	    distance += node_distance( seq_dist, node1->left,  node2);
	    distance += node_distance( seq_dist, node1->right, node2);
	} else 	if ( node2->type != LEAF ) {
	    distance += node_distance( seq_dist, node1,  node2->left);
	    distance += node_distance( seq_dist, node1,  node2->right);
	}
	return distance;
    }
}
/********************************************************************************/
/********************************************************************************/
int print_tree (FILE * fptr,Node * node ) {
    
    if ( node->type != LEAF ) {
	fprintf ( fptr, "(");
	/* down the left branch*/ 
	print_tree ( fptr, node->left);
	fprintf ( fptr,","); 
	/* down the right  branch*/ 
	print_tree ( fptr, node->right);
	/* handle unrooted tree: */
	if (node->type != ROOT && node->parent->type==LEAF ) {
	    fprintf ( fptr,","); 
	    print_tree ( fptr, node->parent);
	}
	
	fprintf ( fptr, ")");
	fprintf ( fptr,"%d", node->id);
	if ( node->dist_to_parent > 0.0001 ) {
	    fprintf ( fptr,":%f ", node->dist_to_parent );
	} else {
	    fprintf ( fptr,":0.0001");
	}
	
    } else {
        fprintf ( fptr,"%s", node->name);
	if ( node->dist_to_parent > 0.0001 ) {
	    fprintf ( fptr,":%f ", node->dist_to_parent);
	} else {
	    fprintf ( fptr,":0.0001");
	}
    }

    return 0;
}


/********************************************************/
void print_debug_tree (FILE * fptr, Node *node) {

  
    if (node->type != LEAF ) {

	fprintf ( fptr, "(");
	/* down the left branch*/ 
	print_debug_tree ( fptr, node->left);
	fprintf ( fptr,","); 
	/* down the right  branch*/ 
	print_debug_tree ( fptr, node->right);
	/* handle unrooted tree: */
	if (node->type != ROOT && node->parent->type==LEAF ) {
	    fprintf ( fptr,","); 
	    print_debug_tree ( fptr, node->parent);
	}
	fprintf ( fptr, ")");
	
	fprintf ( fptr,"%d", node->id);
	if ( node->dist_to_parent > 0.0001 ) {
	    fprintf ( fptr,":%f \n", node->dist_to_parent );
	} else {
	    fprintf ( fptr,":0.0001\n");
	}
	
    } else {
	
        fprintf ( fptr,"%s [%d] ", node->name, node->id);
	if ( node->dist_to_parent > 0.0001 ) {
	    fprintf ( fptr,":%f \n", node->dist_to_parent);
	} else {
	    fprintf ( fptr,":0.0001\n");
	}
	
    }
    
    return ;
}


/**********************************************************************************************/
int place_root ( Node *leaf, int no_of_nodes, Node * root) {
    
    int edge_ctr = 0, longest_ctr = 0, i;
    int find_path ( Node ** path, int * edge_ctr, Node ** longest_path,  int * longest_ctr,Node * node, Node* from) ;

    Node** path, ** longest_path;
    if ( ! (path = emalloc (no_of_nodes*sizeof(Node*) ) ) ) exit (1);
    if ( ! (longest_path = emalloc (no_of_nodes*sizeof(Node*) ) ) ) exit (1);
    /* print all possible paths: */
    /* find any leaf to start fom: */
    longest_ctr = 0;
    edge_ctr = 0;
    for (i = 0; i < no_of_nodes; i++ ) {
	if ( (leaf + i) -> type != LEAF  ) continue;
	find_path (path, &edge_ctr, longest_path, &longest_ctr, leaf + i, NULL);
   }
    /* printf (" %4d  %s --->  %s \n", longest_ctr, longest_path[0]->name,  longest_path[longest_ctr-1]->name); */
    
    /* place the root above the middle node on that path */
    {
	int insert_node ( Node * new_node, Node * old_node_1, Node * old_node_2);
	int reorient_nodes ( Node * node, Node * from );
	root->type = ROOT;
	insert_node (root, longest_path[longest_ctr/2-1], longest_path[longest_ctr/2]);
 	reorient_nodes ( root, NULL);
    }
    return 0;
}

/**********************************************************************************************/
int insert_node ( Node * new_node, Node * old_node_1, Node * old_node_2) {
    new_node->left = old_node_1;
    new_node->right= old_node_2;

    if  ( old_node_1 ==  old_node_2->left ) {
	old_node_2->left = new_node;
    } else if  ( old_node_1 ==  old_node_2->right ){
	old_node_2->right = new_node;
    } else if  ( old_node_1 ==  old_node_2->parent ){
	old_node_2->parent = new_node;
	old_node_2->dist_to_parent = old_node_1->dist_to_parent =  old_node_2->dist_to_parent/2;
	
    } else {
  	fprintf ( stderr, "Error inserting node.\n");
	exit (1);
    }

    if  ( old_node_2 ==  old_node_1->left ) {
	old_node_1->left = new_node;
    } else if  ( old_node_2 ==  old_node_1->right ){
	old_node_1->right = new_node;
    } else if  ( old_node_2 ==  old_node_1->parent ){
	old_node_1->parent = new_node;
 	old_node_2->dist_to_parent = old_node_1->dist_to_parent =  old_node_1->dist_to_parent/2;
   } else {
  	fprintf ( stderr, "Error inserting node.\n");
	exit (1);
    }

    new_node->dist_to_left  =  old_node_1 -> dist_to_parent;
    new_node->dist_to_right =  old_node_2 -> dist_to_parent;
    return 0;
}
/**********************************************************************************************/
int reorient_nodes ( Node * node, Node * from ) {

    Node *left, *right;
    double dist_to_left, dist_to_right, dist_to_parent;
    /* decide on "left" and "right" */
    if (  from ==  node->parent  ) {
	
	left = node->left;
	right = node->right;
	
	dist_to_left   = node->dist_to_left;
	dist_to_right  = node->dist_to_right;
	dist_to_parent =  node->dist_to_parent;
	
    } else if (from  ==  node->left  ) {
	
	right = node->right;
	left  = node->parent;
	
	dist_to_right  = node->dist_to_right;
	dist_to_left   = node->dist_to_parent;
	dist_to_parent = node->dist_to_left;
	
    } else if ( from == node->right  ) {
	
	left   = node->left;
	right  = node->parent;
	
	dist_to_left   = node->dist_to_left;
	dist_to_right  = node->dist_to_parent;
	dist_to_parent = node->dist_to_right;

	
	
    } else {
	fprintf ( stderr, "Error traversing tree in reorient_nodes.\n"); 
	exit (1);
    }
    node->parent = from; 
    node->left   = left; 
    node->right  = right; 
    node->dist_to_left   = dist_to_left;
    node->dist_to_right  = dist_to_right;
    node->dist_to_parent = dist_to_parent;
 
    if ( node->type != LEAF ) { 
	reorient_nodes (left,  node); 
	reorient_nodes (right, node);
    }
    return 0;
}
/**********************************************************************************************/
int find_path ( Node ** path, int * edge_ctr, Node ** longest_path,  int * longest_ctr,Node * node, Node* from) {

    if ( ! node )  {
	fprintf ( stderr, "Error traversing tree in find_path (incomming node = 0x0).\n");
	exit (1);
    }
    
    /* add yourself to the path */
    path[*edge_ctr] = node;
    /* increase the counter     */
    (*edge_ctr) ++;
    if (  node->type != LEAF ) {
	Node *left, *right;
	/* decide on "left" and "right" */
	if (  node->parent == from ) {
	    left = node->left;
	    right = node->right;
	} else if ( node->left == from ) {
	    right = node->right;
	    left  = node->parent;
	} else if ( node->right == from ) {
	    right = node->parent;
	    left  = node->left;
	} else {
	    fprintf ( stderr, "Error traversing tree in find_path.\n");
	    exit (1);
	}
	/* go down the left */
	find_path (path, edge_ctr, longest_path, longest_ctr,  left, node);
	/* go down the   right */
	find_path (path, edge_ctr, longest_path, longest_ctr,  right, node);
    } else {
	if ( *edge_ctr == 1 ) { /* we are at the beginning */
	    find_path (path, edge_ctr, longest_path, longest_ctr,node->parent, node);
	} else {  /* we are at the end */
	    if ( *edge_ctr > *longest_ctr ) {
		*longest_ctr = *edge_ctr;
		memcpy ( longest_path, path, (*longest_ctr)*sizeof(Node*));
	    }
	}
    }
    /* get yourself off the path */
    path[*edge_ctr] = NULL;
    /* decrease the counter */
    (*edge_ctr) --;
  
    return 0;
}


/**********************************************************************************************/
int rank_inner_nodes (Node *root, int  no_of_nodes) {

    int i, inner_node_ctr;
    int  * node_label_sorted;
    double *distance;
    Node ** label2node;
    int find_leaves_dist ( Node *node, int *inner_node_ctr, int * node_label_sorted, Node ** label2node, double * distance );
    
    /* allocate */
    if ( ! ( node_label_sorted = emalloc (no_of_nodes*sizeof(int)))) return 1;
    if ( ! ( distance  = emalloc (no_of_nodes*sizeof(double)))) return 1;
    if ( ! ( label2node  = emalloc (no_of_nodes*sizeof(Node*)))) return 1;
    
    /*find distance to the root from each inner node */
    inner_node_ctr = 0;
    find_leaves_dist ( root, &inner_node_ctr, node_label_sorted, label2node, distance );
   

    /* sort inner nodes by that distance */
    array_qsort ( node_label_sorted, distance, inner_node_ctr);
    /* assign them their sorted number */
    for (i = 0; i < inner_node_ctr; i++ ) {
	label2node[node_label_sorted[i]]->id = i+1;
    }
   
    free (node_label_sorted);
    free (label2node);
    free (distance);

    return 0;
}
/**********************************************************************************************/
int find_leaves_dist ( Node *node, int *inner_node_ctr, int * node_label_sorted, Node ** label2node, double * distance ){
    double distance_to_parent (Node * node);
 
    if ( node -> type == LEAF ) {
    } else {
	node_label_sorted[*inner_node_ctr] = *inner_node_ctr;
	label2node[*inner_node_ctr] = node;
	distance[*inner_node_ctr] = distance_to_parent (node);
 	(*inner_node_ctr) ++;
	find_leaves_dist ( node->left, inner_node_ctr, node_label_sorted, label2node, distance );
	find_leaves_dist ( node->right, inner_node_ctr, node_label_sorted, label2node, distance );
    }

    return 0;
}



/**********************************************************************************************/
double distance_to_parent (Node * node){

    Node * current = node ;
    double distance= 0.0;
    
    while ( current->parent ) {
	distance += current -> dist_to_parent;
	current = current->parent;
    };

    return distance;
    
}

/******************************************************************************/
int find_group_roots (Node *node, Node*** group_root,  int rank, int * group_tag ) {
    if ( node->type==LEAF ) {
    } else {
	
	/* determining the roots for all the groups at the rank "rank" */
	/* remember:  at rank N, N is the smallest group root possible */
	if ( node->id == rank ) {
	    ++(*group_tag);
	    group_root[rank][*group_tag] = node;
	} else if( node->id < rank ){
	    
	    Node *left, *right;
	    left  = node->left;
	    right = node->right;

	    if ( left->id >= rank  || left->type==LEAF) {
		++(*group_tag);
		group_root[rank][*group_tag] = left;
	    } else {
		find_group_roots ( left, group_root, rank, group_tag);
	    }
	    if ( right->id >= rank  || right->type==LEAF) {
		++(*group_tag);
		group_root[rank][*group_tag] = right;
	    } else {
		find_group_roots ( right, group_root, rank, group_tag);
	    }
	}
	
    }
    return 0;
}


/******************************************************************************/
int find_group_roots_sim (Node *root, int no_of_nodes, int no_of_bins, Node*** group_root) {

    int ctr,  *bin_ctr;
    double  bin_size = 1.0/no_of_bins;
    int find_group_roots_sim_rec_part (Node *node,  Node*** group_root, int * bin_ctr, double bin_size) ;
    
    if ( ! (bin_ctr = (int*) emalloc (no_of_bins*sizeof(int))) ) exit (1);
    for (ctr=0; ctr< no_of_bins; ctr++) {
	memset (group_root[ctr], 0, no_of_nodes);
    }
    find_group_roots_sim_rec_part (root, group_root, bin_ctr, bin_size);

# if 0
    int g;
    for (ctr=0; ctr< no_of_bins; ctr++) {
	printf (" %d \n", ctr);
	for (g=0; g < no_of_nodes && group_root[ctr][g]; g++) {
	    printf ("\t %2d  %3d \n", g, group_root[ctr][g]->id );
	}
    }
    exit(1);
# endif
    free(bin_ctr);
    return 0;
}


/******************************************************************************/
int find_group_roots_sim_rec_part (Node *node,  Node*** group_root, int * bin_ctr, double  bin_size) {

    int bin;
    double sim;
    /* if leaf, asign bin value to -1 and return*/
    if ( node->type == LEAF ) {
	node->bin = -1;
	return 0;
    }

    /* take care of itself first, than of the children */
    /* find bin */
    sim = node -> avg_sim;
    bin = floor ( (int)(sim/bin_size) );
    /* is bin different than parents? */
    if ( node->type ==ROOT || node->parent->bin != bin) {
	/* yes - assign */
	group_root[bin][bin_ctr[bin]] = node;
	node->bin = bin;
	bin_ctr[bin] ++;
    } else { /* still asign node->bin, so that children know */
	node->bin = bin;
    }
    /* children */
    find_group_roots_sim_rec_part ( node->left, group_root, bin_ctr, bin_size);
    find_group_roots_sim_rec_part ( node->right, group_root, bin_ctr, bin_size);
    
    return 0;
}
/********************************************************************************/
/********************************************************************************/
int  set_precedence_table ( Node* node, int number_of_nodes, int **is_precedent){

    int node_id;
    
    if (node->left->type != LEAF )  {
	is_precedent[node->id][node->left->id] = 1;
	set_precedence_table ( node->left, number_of_nodes, is_precedent);
	for ( node_id = node->left->id + 1; node_id <= number_of_nodes; node_id++ ) {
	    is_precedent[node->id][node_id] |= is_precedent[node->left->id][node_id];
	}
    }
    if (node->right->type != LEAF ) {
	is_precedent[node->id][node->right->id] = 1;
	set_precedence_table ( node->right, number_of_nodes, is_precedent);
	for ( node_id = node->right->id + 1; node_id <= number_of_nodes; node_id++ ) {
	    is_precedent[node->id][node_id] |= is_precedent[node->right->id][node_id];
	}
    }
     
    return 0;
}


/********************************************************************************/
int  print_leaves (FILE * fptr, Node * node) {

    if ( node->type == LEAF ) {
	fprintf (fptr,"\t\t%s\n", node->name);
    } else {
	print_leaves (fptr, node->left);
	print_leaves (fptr, node->right);
    }
    
    return 0;
}
