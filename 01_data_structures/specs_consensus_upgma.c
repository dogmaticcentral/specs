# include "specs.h"
/********************************************************************************/
/********************************************************************************/
int closest_in_curr_list_consensus_upgma (Alignment *alignment, Node * node,
					  int tree_size, int dummy, 
					  int *  current_list, int * closest1,
					  int *closest2, double * dist_1_ptr,
					  double * dist_2_ptr ){
    int ctr1, ctr2;
    int consensus_length, max_consensus_length;
    int seq_length = alignment->length;
    int node_consensus_length (int seq_length, Node * node1, Node* node2);
    double distance, distance_closest;
    
    max_consensus_length = -1;
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
        if ( ! current_list[ctr1]) continue;
	for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    consensus_length = node_consensus_length ( seq_length, node+ctr1, node+ctr2 );
	 
	    if ( consensus_length > max_consensus_length) {
		printf ( " ***  %4d  %4d   %4d \n",
			 seq_length, max_consensus_length, consensus_length);
		max_consensus_length = consensus_length;
		*closest1 = ctr1;
		*closest2 = ctr2;
		distance_closest = node_distance (alignment->seq_dist,  node+ctr1, node+ctr2 );
		distance_closest /= (node+ctr1)->number_of_leaves*(node+ctr2)->number_of_leaves;
		
	    } else if ( consensus_length == max_consensus_length ) {
		/* tie break: */
		/* look for smaller avg distance - like in the regular upgma */
		distance = node_distance ( alignment->seq_dist,  node+ctr1, node+ctr2 );
		distance /= (node+ctr1)->number_of_leaves*(node+ctr2)->number_of_leaves;
		if ( distance < distance_closest) {
		    distance_closest = distance;
		    *closest1 = ctr1;
		    *closest2 = ctr2;
		} else {
		    printf ( "tie ...\n");
		}
	    }
	}
    }
    *dist_1_ptr =  *dist_2_ptr = (1.- (double)max_consensus_length/seq_length)/2;

    return 0;
}
/********************************************************/
int node_consensus_length (int seq_length, Node * node1, Node* node2) {

    int i, ctr =0;

    for (i=0; i <seq_length; i++) {
	ctr += ( node1->seq[i] &&  node1->seq[i]!= '.' &&  node1->seq[i] == node2->seq[i] );
    }

    return ctr;
    
}
