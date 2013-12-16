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

int hybrid (Alignment * alignment, Tree * tree, int * similar_to,
	    int normalize,  double *score){
    
    int pos, no_seqs = alignment->number_of_seqs;
    int rank, g;
    double subsum, rho;

    /* calculate the score */
    for ( pos=0; pos < alignment->length; pos++) {
	entropy_recursive(tree->root, pos, similar_to, normalize);
	rho = 1.0;
	for ( rank = 1; rank < no_seqs; rank++ ) {
	    subsum = 0.0;
	    for ( g=0; g<rank; g++) {
		subsum  +=  tree->group_root[rank][g]->entropy;
	    }
	    rho += subsum/rank;
	}
 	score[pos] =  rho;
    }

    return 0;
    
}

/******************************************************************************/
int entropy_recursive (Node *node,  int pos, int * similar_to, int normalize) {

    if ( node->type == LEAF ) {
	int  aa;
	memset (node->population, 0, ASCII*sizeof(int));
	aa = (int) node->seq[pos];
	if ( similar_to) aa = similar_to[aa];
	node->population[ aa ] = 1;
	node->entropy = 0.0;
    } else {
	int ctr;
	int no_of_leaves = node->number_of_leaves;
	int no_of_types = 0;
	double fr;
	entropy_recursive ( node->left, pos, similar_to, normalize);
	entropy_recursive (node->right, pos, similar_to, normalize);
	/* for each aa type: population = population left + pop right */
	node->entropy = 0.0;
	if ( normalize ) {
	    ctr = '.';
	    node->population[ctr] = node->left->population[ctr]
		+ node->right->population[ctr];
	    no_of_leaves -= node->population[ctr];
	}
	
	for (ctr=0; ctr < ASCII; ctr++ ) {
	    if (normalize && ctr == '.') continue;
	    node->population[ctr] = node->left->population[ctr]
		+ node->right->population[ctr];
	    /* frequency is pop/number of seqs */
	    if ( node->population[ctr] ) {
		if ( normalize ) no_of_types++;
		fr =  (double)node->population[ctr]/no_of_leaves;
		node->entropy -= fr*log(fr);
	    }
	}
	if (normalize && no_of_types>1)
	    node->entropy /= log( no_of_types) ;
   }

    return 0;
}

/******************************************************************************/
int entropy_complement (Node *node, Node * root,  int pos) {

    if ( node->type == LEAF ) {
	node->entropy_complement = 0.0;
    } else {
	int ctr, complement_population;
	int complement_size = root->number_of_leaves - node->number_of_leaves;
	double fr;
	entropy_complement (node->left, root, pos);
	entropy_complement (node->right, root, pos);
	/* for each aa type: population = population left + pop right */
	node->entropy_complement = 0.0;
	
	if ( node->type != ROOT ) {
	    for (ctr=0; ctr < ASCII; ctr++ ) {
		complement_population = root->population[ctr] - node->population[ctr]; 
		/* frequency is pop/number of seqs */
		if ( complement_population ) {
		    fr =  (double)complement_population/complement_size;
		    node->entropy_complement -= fr*log(fr);
		}
	    }
	}
    }

    return 0;
}



