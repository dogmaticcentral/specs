# include "specs.h"

/***********************************************************************************************************/

# define MAX_OVERLAP 0.3
# define MAX_PROB 1.e-6

int output_specs_score ( Options *options, Protein * protein, Alignment * alignment,
			 int *almt2prot, Node * leaf,
			 double ** score, double ** complement_score, double ** p_value,
			 double **probability, double **overlap, 
			 double cutoff, int ** is_precedent) {
# if 0
    char *base_filename = options ->outname;
    int node_id, pos, no_seqs = alignment->number_of_seqs;
    int node_id_2;
    int node_ctr,* print, first;
    int ** related;
    int int_max_gaps = MAX_GAPS*no_seqs;
    Node *root;
    int  refseq_number,  refseq_2_number = -1;
    char filename[BUFFLEN];
    char in_qry2;
    FILE * fptr;
    int print_leaves (FILE * fptr, Node * node);
    int recursive_cleanup (Node *node, int * related, double * probability);
    
    if ( ! ( related   = intmatrix (alignment->length+1, no_seqs ) ) )   exit (1);
    if ( ! ( print = emalloc (no_seqs*sizeof(int)  ) ) )  exit (1);
    printf ("in output\n");

#if 1
    sprintf (filename, "%s.specs", base_filename);
    fptr = efopen (filename, "w");
    if ( !fptr) return 1;
# else
    fptr = stdout;
# endif

    if ( !  options->refseq ) {
	fprintf ( stderr, "Please define refseq.\n");
	exit (0);
    }
    /* locate the refseq */
    refseq_number = find_refseq (alignment, options->refseq);
    if ( options->qry2[0]) refseq_2_number = find_refseq (alignment, options->qry2);
    
    
    /* associate nodes and positions */
    for ( node_id = 2; node_id < no_seqs; node_id++ ) {
	node_ctr = 2*no_seqs -1 - node_id;
	if ( (leaf+node_ctr)-> number_of_leaves < MIN_NO_LEAVES) continue;
	for ( pos=0; pos < alignment->length; pos++) {
	    if ( alignment->column_gaps[pos] > int_max_gaps) continue;
	    if ( overlap[pos][node_id] < MAX_OVERLAP &&  probability[pos][node_id]  <= MAX_PROB) {
		related[pos][node_id] = 1;
	    }
	}
    }
    /* recursive cleanup (so parents and children do not appear for the same reason */
    for ( pos=0; pos < alignment->length; pos++) {
	if (	alignment->column_gaps[pos] > int_max_gaps ) continue;
	recursive_cleanup ( root = leaf+2*no_seqs-1 - 1, related[pos], probability[pos] ) ;
    }
    
    /* which nodes should be printed? */
    for ( node_id = 2; node_id < no_seqs; node_id++ ) {
 	for ( pos=0; pos < alignment->length; pos++) {
	    if ( alignment->column_gaps[pos] > int_max_gaps ) continue;
	    if ( related[pos][node_id] ) {
		print[node_id] = 1;
		break;
	    }
	}
    }
    for ( node_id = 2; node_id < no_seqs; node_id++ ) {
	node_ctr = 2*no_seqs-1 - node_id;

	if ( ! print[node_id] ) continue;
	fprintf ( fptr, "\n\nnode %d\n=============\n", node_id);
	
	fprintf ( fptr, " %4s  %4s  %4s %4s %5s  %5s  %5s   %5s   other nodes\n",
		    "pos", "pdbid", "qry1", "qry2", "entr", "c_entr", "prob", "ovlp");
	
	for ( pos=0; pos < alignment->length; pos++) {
	    if (  ! related[pos][node_id] ) continue;
	    in_qry2 =  (options->qry2[0]) ? alignment->sequence[refseq_2_number][pos]: '-';
	    fprintf ( fptr, " %4d   %4s   %1c   %1c  %5.1lf    %5.1lf   %5.1le   %5.2lf   ",
		      pos+1, protein->sequence[ almt2prot[pos]].pdb_id, 
		      alignment->sequence[refseq_number][pos],
		      in_qry2,
		      score[pos][node_id], complement_score[pos][node_id],
		      probability[pos][node_id], overlap[pos][node_id] );

	    /* other nodes which have this same pos as discriminant  */
	    first = 1;
	    for ( node_id_2 = 2; node_id_2  < no_seqs; node_id_2 ++ ) {
		if ( node_id_2 == node_id ) continue;
		if ( related[pos][node_id_2] ) {
		    if ( ! first ) {
			fprintf ( fptr, ",");
		    }
		    fprintf ( fptr, "%d",node_id_2); 
		    first = 0;
		}
	    }
	   
	    fprintf ( fptr, "\n");
	}
    }

    fclose (fptr);

    /* print out leaves */
    sprintf (filename, "%s.leaves", base_filename);
    fptr = efopen (filename, "w");
    if ( !fptr) return 1;
    

    for ( node_id = 2; node_id < no_seqs; node_id++ ) {
	if ( ! print[node_id] ) continue;
	node_ctr = 2*no_seqs -1 - node_id;
	fprintf ( fptr, "\n\nnode %d\n=============\n", node_id);
	print_leaves ( fptr, leaf+node_ctr);
	fprintf (fptr, "\n");
   }
    
    
    fclose (fptr);

# endif    
    return 0;
}


/*****************************************************************/
int recursive_cleanup (Node *node, int * related, double * probability) {

    if ( node->number_of_leaves < MIN_NO_LEAVES || node->type == LEAF ) {
	return 0;
    }
    recursive_cleanup (node->left,  related, probability);   
    recursive_cleanup (node->right, related, probability);

    if ( ! related [node->left->id] && ! related [node->right->id] ) {
	return 0;
    }
    /* if both related */
    if (  related [node->left->id] &&  related [node->right->id] ) {
	int distr_overlap ( int * pool_population, int * node_population, int no_bins, double * ovlp );
	double ovlp;
	distr_overlap ( node->left->population, node->right->population, ASCII, & ovlp );
	
	if ( ovlp > 0.8) { /* if the distributions overlap, keep the parent */
	    related [node->right->id] = 0;
	    related [node->left->id] = 0;
	    related [node->id] = 1;
	} else {  /* keep the children and drop the parent */
	     related [node->id] = 0;
	}
    } else {
	//return 0;
	int child_id;
	/* if one related, compare with parent, and keep the one with the smaller p*/
	child_id =  related [node->left->id] ? related [node->left->id] : related [node->right->id];
	related [node->id] = 0;
	related [child_id] = 1;
	    
/* 	if ( probability[node->id] > probability [child_id] ) { */
/* 	    related [node->id] = 0; */
/* 	    related [child_id] = 1; */
/* 	} else { */
/* 	    related [node->id] = 1; */
/* 	    related [child_id] = 0; */
/* 	} */
    }

    return 0;
}
