# include "specs.h"

int output_specs_score ( char * base_filename, Protein * protein, Alignment * alignment, Node * leaf,
			 double ** score, double ** complement_score, double cutoff, int ** is_precedent) {

    int node_id, pos, no_seqs = alignment->number_of_seqs;
    int ctr, neg_ctr, complement_ctr;
    int number_of_nodes_above_cutoff,  number_of_nodes_complement_above_cutoff, number_of_nodes_below_neg_cutoff;
    int ctr1, ctr2, node_id_1, node_id_2;
    int node_ctr;
    int node_pos_ctr, ** node_related_position;
    Node * node_ptr, *first_inner_node;
    int * above_cutoff = NULL;
    int * complement_above_cutoff = NULL;
    int * below_neg_cutoff = NULL;
    char filename[BUFFLEN];
    FILE * fptr;
    int  print_leaves (FILE * fptr, Node * node);

    if ( ! ( above_cutoff = emalloc (no_seqs*sizeof(int) ) ) )  exit (1);
    if ( ! ( complement_above_cutoff = emalloc (no_seqs*sizeof(int) ) ) )  exit (1);
    if ( ! ( below_neg_cutoff = emalloc (no_seqs*sizeof(int) ) ) )  exit (1);
    if ( ! (  node_related_position   = intmatrix ( no_seqs, alignment->length+1) ) )   exit (1);

# if 1
    pos = 44;
    for ( node_id= 1; node_id < no_seqs; node_id++ ) {
	//if ( complement_score[pos][node_id] > cutoff ||
	//score[pos][node_id] <  -2.5 ) 
	    fprintf (stdout, "%4d     %8.3f    %8.3f \n", node_id,
		     score[pos][node_id], complement_score[pos][node_id] );
    }
    exit (0);
# endif
    
    sprintf (filename, "%s.specs", base_filename);
    fptr = efopen (filename, "w");
    if ( !fptr) return 1;

    for ( pos=0; pos < alignment->length; pos++) {
	
	number_of_nodes_above_cutoff = 0;
	ctr = 0; neg_ctr = 0;complement_ctr = 0;
	for ( node_id= 1; node_id < no_seqs; node_id++ ) {

	    node_ctr = 2*no_seqs-1 - node_id;
	    if ( (leaf+node_ctr)-> number_of_leaves < 5) continue;
	    
	    if ( score[pos][node_id] > cutoff ){
		above_cutoff [ctr] = node_id;
		ctr ++;
	    } else if (score[pos][node_id] < -1*cutoff){
		below_neg_cutoff [ctr] = node_id;
		neg_ctr ++;
	    }
	    if ( complement_score[pos][node_id] > cutoff ){
		complement_above_cutoff [ctr] = node_id;
		complement_ctr ++;
	    } 	    
	}
	if ( !ctr && ! neg_ctr) continue;
	above_cutoff [ctr] = -1;
	below_neg_cutoff [neg_ctr] = -1;
	complement_above_cutoff [ctr] = -1;
	
	number_of_nodes_above_cutoff = ctr;
	number_of_nodes_below_neg_cutoff = neg_ctr;
	number_of_nodes_complement_above_cutoff = complement_ctr;

	/********************************  above cutoff ****************************************************/
	for (ctr1=0; above_cutoff [ctr1] >=  0; ctr1++ ) {
	    node_id_1 = above_cutoff [ctr1];
	    if ( ! node_id_1 ) continue;
	    node_ctr = 2*no_seqs-1 - node_id_1;
	    if ( (leaf+node_ctr)-> number_of_leaves < 5) continue;
	    
	    for (ctr2 = ctr1+1; above_cutoff [ctr2] >= 0 ; ctr2++ ) {
		node_id_2 = above_cutoff [ctr2];
		if ( ! node_id_2 ) continue;
		node_ctr = 2*no_seqs-1 - node_id_2;
		if ( (leaf+node_ctr)-> number_of_leaves < 5) continue;

		/* go for the highest z */
		/* if ( is_precedent[node_id_1][node_id_2] || is_precedent[node_id_2][node_id_1] ) {
		    number_of_nodes_above_cutoff --;
		    
		    if ( score[pos][node_id_2] < score[pos][node_id_1] ) {
			above_cutoff [ctr2] = 0;
		    } else {
			above_cutoff [ctr1] = 0;
			break;
		    }
		    }
		*/
		
		/* alternatively, go for the ancestor */
		if ( is_precedent[node_id_1][node_id_2] ) {
		    number_of_nodes_above_cutoff --;
		    above_cutoff [ctr2] = 0;
		} else if ( is_precedent[node_id_2][node_id_1] ) {
		    number_of_nodes_above_cutoff --;
		    above_cutoff [ctr1] = 0;
		    break;
		}
	    }
	}
	if (  number_of_nodes_above_cutoff < 0) {
	    printf (" %d foul\n", pos);
	    exit (1);
	}
	/********************************  complement above cutoff ********************************************/
	for (ctr1=0; complement_above_cutoff [ctr1] >=  0; ctr1++ ) {
	    node_id_1 = complement_above_cutoff [ctr1];
	    if ( ! node_id_1 ) continue;
	    node_ctr = 2*no_seqs-1 - node_id_1;
	    if ( (leaf+node_ctr)-> number_of_leaves < 5) continue;
	    
	    for (ctr2 = ctr1+1; complement_above_cutoff [ctr2] >= 0 ; ctr2++ ) {
		node_id_2 = complement_above_cutoff [ctr2];
		if ( ! node_id_2 ) continue;
		node_ctr = 2*no_seqs-1 - node_id_2;
		if ( (leaf+node_ctr)-> number_of_leaves < 5) continue;
		
		/* go for the ancestor */
		if ( is_precedent[node_id_1][node_id_2] ) {
		    number_of_nodes_complement_above_cutoff --;
		    complement_above_cutoff [ctr2] = 0;
		} else if ( is_precedent[node_id_2][node_id_1] ) {
		    number_of_nodes_complement_above_cutoff --;
		    complement_above_cutoff [ctr1] = 0;
		    break;
		}
	    }
	}
	if (  number_of_nodes_complement_above_cutoff < 0) {
	    printf (" %d foul\n", pos);
	    exit (1);
	}

	/*************************** far below the cutoff: *****************************************************/ 
	for (ctr1=0; below_neg_cutoff [ctr1] >=  0; ctr1++ ) {
	    node_id_1 = below_neg_cutoff [ctr1];
	    if ( ! node_id_1 ) continue;
	    node_ctr = 2*no_seqs-1 - node_id_1;
	    if ( (leaf+node_ctr)-> number_of_leaves < 5) continue;
	    
	    for (ctr2 = ctr1+1; below_neg_cutoff [ctr2] >= 0 ; ctr2++ ) {
		node_id_2 = below_neg_cutoff [ctr2];
		if ( ! node_id_2 ) continue;
		node_ctr = 2*no_seqs-1 - node_id_2;
		if ( (leaf+node_ctr)-> number_of_leaves < 5) continue;

		/* alternatively, go for the ancestor */
		if ( is_precedent[node_id_1][node_id_2] ) {
		    number_of_nodes_below_neg_cutoff --;
		    below_neg_cutoff [ctr2] = 0;
		} else if ( is_precedent[node_id_2][node_id_1] ) {
		    number_of_nodes_below_neg_cutoff --;
		    below_neg_cutoff [ctr1] = 0;
		    break;
		}
	    }
	}
	if (  number_of_nodes_below_neg_cutoff < 0) {
	    printf (" %d foul\n", pos);
	    exit (1);
	}

	/*********************** the actual output, finally ****************************************************/
	if ( number_of_nodes_above_cutoff > 0 ) {
	    if ( number_of_nodes_above_cutoff > 1 ) {
		fprintf (fptr, "%4d  %5s:", pos + 1, "discr" );
	    } else {
		fprintf (fptr, "%4d  %5s:", pos + 1, "det" );
	    }
	    for (ctr1=0; above_cutoff [ctr1] >=  0; ctr1++ ) {
		node_id_1 = above_cutoff [ctr1];
		if ( ! node_id_1 ) continue;
		node_ctr = 2*no_seqs-1 - node_id_1;
		if ( (leaf+node_ctr)-> number_of_leaves < 5) continue;
		fprintf (fptr, " %4d", node_id_1);
		node_ctr = 2*no_seqs-1 -  node_id_1;
		(leaf+node_ctr) -> marked = 1;
		node_related_position[node_id_1][0] ++;
		node_pos_ctr =  node_related_position[node_id_1][0]; /* use the zeroth pos as a counter */
		node_related_position[node_id_1][ node_pos_ctr ] = pos + 1;
	    }
	    fprintf (fptr, "\n");
	}
	if ( number_of_nodes_complement_above_cutoff > 0 ) {
	    fprintf (fptr, "\t %4d  %5s:", pos + 1, "compl" );
	    for (ctr1=0; complement_above_cutoff [ctr1] >=  0; ctr1++ ) {
		node_id_1 = complement_above_cutoff [ctr1];
		if ( ! node_id_1 ) continue;
		node_ctr = 2*no_seqs-1 - node_id_1;
		if ( (leaf+node_ctr)-> number_of_leaves < 5) continue;
		fprintf (fptr, " %4d", node_id_1);
		node_ctr = 2*no_seqs-1 -  node_id_1;
		(leaf+node_ctr) -> marked = 1;
		node_related_position[node_id_1][0] ++;
		node_pos_ctr =  node_related_position[node_id_1][0]; /* use the zeroth pos as a counter */
		node_related_position[node_id_1][ node_pos_ctr ] = pos + 1;
	    }
	    fprintf (fptr, "\n");
	}
	if ( number_of_nodes_below_neg_cutoff > 0 ) {
	    fprintf (fptr, "\t %4d  %5s:", pos + 1, "neg discr" );
	    for (ctr1=0; below_neg_cutoff [ctr1] >=  0; ctr1++ ) {
		// printf (" ** %d \n", below_neg_cutoff [ctr1] );
		node_id_1 = below_neg_cutoff [ctr1];
		if ( ! node_id_1 ) continue;
		node_ctr = 2*no_seqs-1 - node_id_1;
		if ( (leaf+node_ctr)-> number_of_leaves < 5) continue;
		fprintf (fptr, " %4d", node_id_1);
		node_ctr = 2*no_seqs-1 -  node_id_1;
		(leaf+node_ctr) -> marked = 1;
		node_related_position[node_id_1][0] ++;
		node_pos_ctr =  node_related_position[node_id_1][0]; /* use the zeroth pos as a counter */
		node_related_position[node_id_1][ node_pos_ctr ] = pos + 1;
	    }
	    fprintf (fptr, "\n");
	}


	

    }
    



    fprintf ( fptr, "\n\n");
    fprintf ( fptr, "nodes:\n=============\n\n");
    
    first_inner_node = leaf + no_seqs;

    for ( node_ctr = 0; node_ctr < no_seqs-1; node_ctr ++ ) {
	node_ptr = first_inner_node + node_ctr;
	if ( node_ptr->marked ) {
	    fprintf ( fptr, "\n\t node %4d\n", node_ptr->id );
	    fprintf ( fptr, "\t positions:");
	    for ( node_pos_ctr=1; node_pos_ctr<= node_related_position[node_ptr->id][0]; node_pos_ctr++ ) {
		fprintf ( fptr, " %4d",  node_related_position [node_ptr->id][ node_pos_ctr ]);
	    }
	    fprintf (fptr, "\n");
	    fprintf ( fptr, "\t sequences\n");
	    print_leaves ( fptr, node_ptr);
	}
    }
    fprintf (fptr, "\n");
    
    free (above_cutoff);
    free (complement_above_cutoff);
    free (below_neg_cutoff);
    free_matrix ( (void*) node_related_position);
    
//fclose (fptr);
    return 0;
}
