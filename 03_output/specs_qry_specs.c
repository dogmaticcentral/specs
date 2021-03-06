/*
This source code is part of specs, an application for
protein residue conservation scoring.
Written by Ivana Mihalek.
Copyright (C) 2007-2013 Ivana Mihalek.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see<http://www.gnu.org/licenses/>.

Contact: ivana.mihalek@gmail.com.
*/
# include "specs.h"
int output_specs_for_refseq ( Options * optin, Protein * protein, Alignment * alignment,
			      int * almt2protein,  Node * leaf,
			      double ** score, double cutoff , int ** is_precedent) {

    int node_id, pos, no_seqs = alignment->number_of_seqs;
    int ctr,number_of_nodes_above_cutoff ;
    int ctr1, ctr2, node_id_1, node_id_2;
    int node_ctr;
    int node_pos_ctr, ** node_related_position;
    Node * node_ptr, * first_inner_node;
    int * above_cutoff = NULL;
    char filename[BUFFLEN];
    FILE * fptr;
    int  print_leaves (FILE * fptr, Node * node);

    if ( ! ( above_cutoff = emalloc (no_seqs*sizeof(int) ) ) )  exit (1);
    if ( ! (  node_related_position   = intmatrix ( no_seqs, alignment->length+1) ) )   exit (1);
    sprintf (filename, "%s.qry_specs", options->outname);
    fptr = efopen (filename, "w");
    if ( !fptr) return 1;

    for ( pos=0; pos < alignment->length; pos++) {
	
	number_of_nodes_above_cutoff = 0;
	ctr = 0;
	for ( node_id= 1; node_id < no_seqs; node_id++ ) {
	    if ( score[pos][node_id] > cutoff ){
		above_cutoff [ctr] = node_id;
		ctr ++;
	    }
	}
	if ( !ctr ) continue;
	above_cutoff [ctr] = -1;
	
	number_of_nodes_above_cutoff = ctr;
	
	for (ctr1=0; above_cutoff [ctr1] >=  0; ctr1++ ) {
	    node_id_1 = above_cutoff [ctr1];
	    if ( ! node_id_1 ) continue;
	    
	    for (ctr2 = ctr1+1; above_cutoff [ctr2] >= 0 ; ctr2++ ) {
		node_id_2 = above_cutoff [ctr2];
		if ( ! node_id_2 ) continue;

		/* go for the highest z */
		/* if ( is_precedent[node_id_1][node_id_2] || is_precedent[node_id_2][node_id_1] ) {
		    number_of_nodes_above_cutoff --;
		    
		    if ( score[pos][node_id_2] < score[pos][node_id_1] ) {
			above_cutoff [ctr2] = 0;
		    } else {
			above_cutoff [ctr1] = 0;
			break;
		    }
		    }*/
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
	if (  number_of_nodes_above_cutoff <= 0) {
	    printf (" %d foul\n", pos);
	    for (ctr1=0; above_cutoff [ctr1] >=  0; ctr1++ ) {
		printf (" ** %d \n", above_cutoff [ctr1] );
	    }
	    exit (1);
	}
	if ( number_of_nodes_above_cutoff > 1 ) {
	    fprintf (fptr, "%4d  %5s:", pos + 1, "discr" );
	} else {
	    fprintf (fptr, "%4d  %5s:", pos + 1, "det" );
	}
	for (ctr1=0; above_cutoff [ctr1] >=  0; ctr1++ ) {
	    // printf (" ** %d \n", above_cutoff [ctr1] );
	    node_id_1 = above_cutoff [ctr1];
	    if ( ! node_id_1 ) continue;
	    fprintf (fptr, " %4d", node_id_1);
	    node_ctr = 2*no_seqs-1 -  node_id_1;
	    (leaf+node_ctr) -> marked = 1;
	    node_related_position[node_id_1][0] ++;
	    node_pos_ctr =  node_related_position[node_id_1][0]; /* use the zeroth pos as a counter */
	    node_related_position[node_id_1][ node_pos_ctr ] = pos + 1;
	}
	fprintf (fptr, "\n");

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
    free_matrix ( (void*) node_related_position);
    
//fclose (fptr);
    return 0;
}
