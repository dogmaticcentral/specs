# include "specs.h"
int closest_in_curr_list_nj (Alignment *alignment, Node * node, int tree_size, int new_node, 
			     int *  current_list, int * closest1, int *closest2,
			     double * dist_1_ptr, double * dist_2_ptr ){
    /* actually, for nj, they are not the closest, but rather the pair which,
       if joined next, gives the minimum total sum of branch lengths */
    /* new_node  is where I'll put the parent for the closest pair I find */
    
    int ctr1, ctr2, curr_list_length;
    int no_leaves = (tree_size+1)/2;
    double sum, min_sum;
    double avg1, avg2, avg1_min, avg2_min;
    double **seq_dist = alignment->seq_dist;
    static double ** dist_table = NULL;
    double  nj_sum_of_branch_lengths (double ** dist_table,  int * current_list,
				   int  node_ctr1, int node_ctr2,
				   int tree_size, double *avg1_ptr, double *avg2_ptr);

    if ( ! dist_table ) { /* initialize distance table */
	/* allocate */
	dist_table = dmatrix (tree_size, tree_size);
	/* copy distances btw the leaves */
	for (ctr1=0; ctr1 < no_leaves; ctr1++ ) {
	    for (ctr2=ctr1+1; ctr2 < no_leaves; ctr2++ ) {
		dist_table[ctr1][ctr2] = dist_table[ctr2][ctr1] = seq_dist[ctr1][ctr2];
	    }
	}
    }
    
    curr_list_length = 0;
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	curr_list_length += current_list[ctr1];
    }

    avg1_min = 0;
    avg2_min = 0;
    min_sum = 10000;
    
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	if ( ! current_list[ctr1]) continue;
	for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    sum = nj_sum_of_branch_lengths ( dist_table, current_list,
					     ctr1, ctr2, tree_size, &avg1, &avg2);
	    if ( sum < min_sum) {
		min_sum   = sum;
		*closest1 = ctr1;
		*closest2 = ctr2;
		avg1_min  = avg1; 
		avg2_min  = avg2; 
	    }
	}
    }
    
    printf ( "closest:  %3d  %3d  \n\n", *closest1, *closest2);
    
    *dist_1_ptr = ( dist_table[*closest1][*closest2]+ avg1_min - avg2_min)/2;
    *dist_2_ptr = ( dist_table[*closest1][*closest2]+ avg2_min - avg1_min)/2;

    dist_table[new_node][*closest1] =  dist_table[*closest1][new_node] = *dist_1_ptr;
    dist_table[new_node][*closest2] =  dist_table[*closest2][new_node] = *dist_2_ptr;

    /* update the distance table */
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	if ( ! current_list[ctr1]) continue;
	if ( ctr1 == new_node )   continue;
	if ( ctr1 == *closest1 ||  ctr1 == *closest2) continue;
	dist_table[new_node][ctr1] = dist_table[ctr1][new_node] =
	    ( dist_table[*closest1][ctr1] + dist_table[*closest2][ctr1]
	      -  dist_table[*closest2][*closest1] )* 0.5;
    }
	
	
 
    if ( curr_list_length == 3 ) {
	for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	    for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
		printf ( " %3d  %3d  %8.4lf \n", ctr1, ctr2,  dist_table[ctr1][ctr2]);
	    }
	    printf ( "\n");
	}
    }
    /* printf (" %4d    %4d  %4d   %8.3lf   %8.3lf   %8.3lf   %8.3lf   %8.3lf  %8.3lf   \n", */
    /*  curr_list_length,  *closest1, *closest2, dist_table[*closest1][*closest2],
	min_distance, avg1_min, avg2_min, *dist_1_ptr, *dist_2_ptr); */
   
    return 0;
}


/***************************************************************************************/
double  nj_sum_of_branch_lengths (double ** dist_table,  int * current_list,
				  int  node_ctr1, int node_ctr2, int tree_size,
				  double *avg1_ptr, double *avg2_ptr) {
    double sum, term;
    double avg[2];
    int current_pair[2] = {node_ctr1,node_ctr2};
    int ctr1, ctr2, n;
    
    sum =  dist_table[node_ctr1][node_ctr2]/2;


    for (ctr1 = 0; ctr1 < 2; ctr1++ ) {
	n =0;
	avg[ctr1] = 0;
	for (ctr2= 0; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    if ( ctr2 == node_ctr1 || ctr2  == node_ctr2 )
		continue;
	    avg[ctr1] += dist_table[ current_pair[ctr1]][ ctr2];
	    n ++;
	}
	if ( !n ) {
	    printf ( "** %d %d \n", node_ctr1, node_ctr2);
	    printf ( "available:  \n");
	    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
		if ( ! current_list[ctr1]) continue;
		printf ( "%4d", ctr1);
	    }
	    printf ( "\n");
	    exit (1);
	}
    
 	    
	avg[ctr1] /= n;
    }

    sum += (avg[0] + avg[1])/2;
    
    *avg1_ptr = avg[0];
    *avg2_ptr = avg[1];
    
    term = 0;
    n = 0;
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	if ( ! current_list[ctr1]) continue;
	if ( ctr1 == node_ctr1 ||  ctr2  == node_ctr2 )
	    continue;
	for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    if ( ctr2 == node_ctr1 ||  ctr2 == node_ctr2 )
		continue;
	    term += dist_table[ctr1][ctr2];
	    n ++;
	}
    }
    sum += term/n;
		
	    
    return sum;
}
