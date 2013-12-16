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

int alignment_tree_sanity_check (Alignment *alignment );
int associate_seqs (Node *node, Alignment *alignment );
int count_leaves ( Node * node);
int find_tree_size (char *filename );
int recursive_tree_scan (FILE * fptr , Node *leaf,  int *node_ctr, Node * parent, int * rooted);
int rank_rooted_tree ( Node * leaf, int no_of_nodes);

/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/

int  read_tree (Options *options,  Alignment *alignment, Tree * tree) {

    int no_of_nodes, tree_size, no_of_leaves;
    int node_ctr, rooted, retval;
    int rank, group_tag;
    FILE *fptr;

    
    /* find the tree size */
    no_of_nodes =  find_tree_size (options->tree_file_name);

    /* allocate */
    tree_size  = no_of_nodes;
    tree->leaf = (Node *) emalloc ( (tree_size+1)*sizeof(Node) ); /* use extra place for the root */
    for ( node_ctr = 0; node_ctr < tree_size; node_ctr++ ) {
	if ( ! ( (tree->leaf+node_ctr)->name = emalloc (NAME_LENGTH))) return 1;
    }

    /* read in the tree */
    if ( ! (fptr=efopen (options->tree_file_name, "r" )))  exit (1);
    node_ctr = rooted = 0;
    retval = recursive_tree_scan (fptr, tree->leaf, &node_ctr, NULL, & rooted);
    fclose (fptr);
    if ( retval ) return retval;
    
    /* assign ranks to the nodes */
    if (rooted) {
	rank_rooted_tree (tree->leaf, no_of_nodes);
	tree->root = tree->leaf;
	tree->size = no_of_nodes;
	
    } else {
	int place_root (Node *leaf, int no_of_nodes, Node *root);
	int rank_inner_nodes (Node *root, int  no_of_nodes);
	int count_leaves (Node * node);
	
	tree->root = tree->leaf + no_of_nodes;  /* the root will be placed in the last position*/
	place_root (tree->leaf, no_of_nodes, tree->root);
	rank_inner_nodes (tree->root, no_of_nodes + 1); /* the root at the last place */
	count_leaves (tree->root);
	tree->size = no_of_nodes + 1;
	
    }

    /* associate sequences with the leaves */
     tree->no_of_leaves = no_of_leaves = (tree->size+1)/2;
    if ( no_of_leaves != alignment->number_of_seqs ) {
	fprintf (stderr, "Alignment/tree size mismatch: no_of_leaves %d number_of_seqs %d\n",
		 no_of_leaves,  alignment->number_of_seqs);
	exit (1);
    }

    associate_seqs (tree->root, alignment);
    alignment_tree_sanity_check(alignment); /* check that all seqs in the alignment
					       have a node associated with it */
   
    /* find groups */
    tree->group_root = node_matrix (tree->no_of_leaves, tree->no_of_leaves);
    for ( rank = 1; rank < tree->no_of_leaves; rank++ ) {
	group_tag = -1;
	find_group_roots (tree->root, tree->group_root, rank, &group_tag);
    }
   
    return 0;
}
    

/**********************************************************************************************/
/**********************************************************************************************/
int alignment_tree_sanity_check (Alignment *alignment ) {
    int s;
    for (s=0; s<alignment->number_of_seqs; s++) {
	if (!alignment->node[s]) {
	    fprintf (stderr, "Sequence %s not found in the tree\n", alignment->name[s]);
	    exit (1);
	    
	}
    }
    return 0;
}

/**********************************************************************************************/
int associate_seqs (Node *node, Alignment *alignment ) {
    
    if ( node->type == LEAF ) {
	
	int s, found = 0;
	for (s=0; s<alignment->number_of_seqs; s++) {
	    if ( ! strcmp( node->name, alignment->name[s] ) ) {
		found = 1;
		node->seq = alignment->sequence[s];
		alignment->node[s] = node;
	    }
	}
	if ( ! found ) {
	    fprintf (stderr, "Leaf name %s not found in the alignment\n", node->name);
	    exit (1);
	}
	
 	
    } else {
	associate_seqs (node->left, alignment);
	associate_seqs (node->right, alignment);
    }
    return 0;
}



/**********************************************************************************************/
int  recursive_tree_scan (FILE * fptr , Node *leaf,  int *node_ctr, Node * parent, int * rooted) {
    
    char next_char;
    char aux_str[NAME_LENGTH];
    int ctr;
    double aux_double;
    Node * node  =  leaf + *node_ctr;
   
    while ( isspace (next_char = getc (fptr)));
    if ( next_char == '(' ) {  /* otherwise we have hit a leaf */
	
	node->type = INNER;

	/* left subtree */
	(*node_ctr)++;
	node -> left = leaf + *node_ctr;
	
	recursive_tree_scan (fptr, leaf, node_ctr, node, rooted);
	

	
	while ( isspace (next_char = getc (fptr)) );

	if ( next_char != ',' ) {
	    fprintf (stderr, "error reading nhx (1)\n");
	    exit (1);
	}
	
	/* right subtree */
	(*node_ctr)++;
	node->right = leaf + *node_ctr;
	
	recursive_tree_scan (fptr, leaf, node_ctr, node, rooted);
		
	
	/* special case  of an unrooted tree */
	while ( isspace (next_char = getc (fptr)) );

	if ( next_char == ',') {
	    (*node_ctr)++;
	    node->parent = leaf + *node_ctr;
	    
	    recursive_tree_scan (fptr, leaf, node_ctr, node, rooted);
	    
	    /* I am my own grandpa */
	    node -> dist_to_parent = node -> parent -> dist_to_parent;

	    *rooted = 0;
	} else if ( next_char == ')' ) {
	    *rooted = 1;
	} else {
	    fprintf (stderr, "error reading nhx (2)\n");
	    exit (2);
	}
    } else {
	/*this is a leaf */
	ungetc  ( next_char, fptr);
	node->type = LEAF;
    }

    /* see if there is a node description ( "name:dist_to_parent" ) */ 
    while ( isspace (next_char = getc (fptr))  );
    ungetc  ( next_char, fptr);
    
    if (next_char == ',' || next_char == ')'|| next_char == ';' ) {
	/* no node description */
	return 0;
    }
   
    /* here */
    memset (aux_str, 0, NAME_LENGTH);
    ctr = 0;
    while (! feof (fptr) &&  (next_char = getc(fptr)) != ':' &&  next_char  != ',' &&  next_char  != ')') {
	aux_str[ctr] = next_char;
	ctr ++;
    }
    if ( next_char == ')' || next_char  == ',' )
	ungetc  ( next_char, fptr);
    if ( ctr ) {
	/* the name is present */
	/* we take that it means a leaf */
	memset ( node -> name, 0, NAME_LENGTH);
	sprintf (node -> name, "%s", aux_str);
    }
    memset (aux_str, 0, NAME_LENGTH);
    ctr = 0;
    while ( isdigit (next_char = getc(fptr) ) ||  next_char=='.' ||  next_char=='-') {
	if ( ctr >= NAME_LENGTH ) {
	    fprintf ( stderr, "Error reading tree: name too long.");
	    return 1; 
	}
 	aux_str[ctr] = next_char;
	ctr ++;
    }
    ungetc ( next_char, fptr);
   
    /* clustalw gives me negative distances (??) --> neighbor joining method does */
    aux_double = atof (aux_str);
    
    node -> dist_to_parent = (aux_double > 0.0) ? aux_double : 0.0 ;    
    node -> parent         = parent;

   
    return 0;
}

/****************************************************************************************/
/****************************************************************************************/

/* assign ranks to nodes in the rooted tree*/
int rank_rooted_tree ( Node * leaf, int no_of_nodes) {

    int  node_ctr, next;
    int min_ctr;
    double dist, min_dist;
    Node * node;

    /*******************************/
    /* assign ranks to inner nodes */ 
    /*******************************/
    /* number of inner nodes =  (no_of_nodes )/2  (integer division!)*/
    for ( next = (no_of_nodes )/2; next >= 1; next -- ) {
	
	/* Run through all the unmarked nodes in the tree. */
	/* The unmarked node with  exactly 2 marked neighbors, */
	/* with the minimal neighbor-to-neighbor distance gets assigned  the `next` value. */
	min_dist = 1000.0;
	min_ctr  = 2;
	for (node_ctr =0; node_ctr <  no_of_nodes; node_ctr++) {

	    node = leaf + node_ctr;
	    if ( node  -> type == LEAF) continue;
	    if ( node -> id) continue;
	    
	    dist = node->left ->dist_to_parent + node->right ->dist_to_parent;
		
	    if ( dist < min_dist) {
		min_dist        = dist;
		min_ctr         = node_ctr;
	    }
	}
	
	(leaf + min_ctr) -> id     = next;
	(leaf + min_ctr) -> type   = (next > 1 ) ? INNER : ROOT;
    }
							   


    return 0;
}

/****************************************************************************************/
int count_leaves ( Node * node) {

    if ( node->type == LEAF ) {
	node->number_of_leaves = 1;
    } else {
	count_leaves ( node->left);
	count_leaves ( node->right);
 	node->number_of_leaves = node->left->number_of_leaves +  node->right->number_of_leaves;
    }

    return 0;
}

/****************************************************************************************/
# define BUFFLEN_TREE  50000

int  find_tree_size (char *filename ) {
    FILE *fptr;
    int node_ctr = 0;
    int bra = 0, ket = 0, comma =  0;
    char line[BUFFLEN_TREE] = {'\0'}, *lineptr;
    if ( ! (fptr=efopen (filename, "r" ) ) )  exit (1);
    /* every open bracket means one iner node*/
    while(fgets(line, BUFFLEN_TREE, fptr)!= NULL){
	lineptr = line;
	while ( *lineptr != '\n' ) {
	    switch ( *lineptr ) {
	    case '(':
		bra++;
		break;
	    case ')':
		ket++;
		break;
	    case ',':
		comma++;
	    }
	    
	    lineptr++;
	}
    }
    fclose (fptr);


    if ( bra != ket ) {
	fprintf ( stderr, "Error reading in the tree: unbalanced brackets.\n");
	fprintf ( stderr, "bra:  %4d     ket:  %4d     comma: %4d \n", bra, ket, comma);
	exit (1);
    }

    if  ( ket < comma ) {
	printf( "Unrooted tree.\n");
    } else {
	printf( "Rooted tree.\n");
    }
    node_ctr = ket + comma + 1;
    printf ("Number of nodes: %d.\n", node_ctr);
    
    return node_ctr; 
}
