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
/*********************************************************************************************/
/*********************************************************************************************/

double  pair_prob[16][128][128]; /* 16 blosum matrices: 30, 35, .., 95, and 100*/


int pheno (Alignment * alignment, double *score, int blosum_ctr) {
    int col, seq, ctr, ctr2;
    int freq[ASCII]; /* ASCII == 128,  ASCII size  */
    int aux;
    double p, entropy, norm;
    for (col = 0; col < alignment->length; col++){
	/* find frequencies */
	memset (freq, 0, ASCII*sizeof(int));
	for (seq=0; seq< alignment->number_of_seqs; seq++){
	    if ( alignment->sequence[seq][col] == '.' ) continue;
	    if ( strchr ( "BXZ" , alignment->sequence[seq][col] )) continue;
	    freq [ (int) alignment->sequence[seq][col] ] ++;
	}
	/* find entropy */
	aux = alignment->number_of_seqs - alignment->column_gaps[col];
	if ( aux < 2 ) continue;
	norm = 1.0 /(aux*(aux-1));
	entropy = 0.0;
	for ( ctr=0; ctr < 128; ctr++) {
	    if ( ! freq[ctr] ) continue;
	    for ( ctr2=ctr; ctr2 < 128; ctr2++) {
		if ( ! freq[ctr2] ) continue;
    
		p = (double)freq[ctr]*freq[ctr2]*norm;
		entropy -=  p*log(p/pair_prob[blosum_ctr][ctr][ctr2]);
	    }
	}
	score[col] = entropy;
    }

    return 0;
}
/*********************************************************************************************/

int pheno_recursive (Node *node, int pos) {
    
    if ( node->type == LEAF ) {
	memset (node->population, 0, ASCII*sizeof(int));
	node->population[(int) node->seq[pos] ] = 1;
	node->entropy = 0.0;
    } else {
	int ctr, ctr2;
	int no_of_leaves = node->number_of_leaves;
	int blosum_ctr = 0; /* <====== calculate */ 
	double p, aux, norm;
	pheno_recursive (node->left, pos);
	pheno_recursive (node->right, pos);
	/* for each aa type: population = population left + pop right */
	for (ctr=0; ctr < ASCII; ctr++ ) {
	    node->population[ctr] = node->left->population[ctr] + node->right->population[ctr];
	}
	blosum_ctr = 1;
	//blosum_ctr = (lround (node->avg_sim*10)*10-30)/5;
	//if ( blosum_ctr < 0)  blosum_ctr = 0;
	node->entropy = 0.0;	
	aux = no_of_leaves - node->population['.'];
	if ( aux < 2 ) return 0;
	norm = 1.0 /(aux*(aux-1));
	for (ctr=0; ctr < ASCII; ctr++ ) {
	    if ( ctr == (int)'.' ) continue;
	    if ( strchr ( "BXZ" , (char)ctr)) continue;
	    if ( ! node->population[ctr] ) continue;
	    for (ctr2=ctr; ctr2 < ASCII; ctr2++ ) {
		if ( ctr2 == (int)'.' ) continue;
		if ( strchr ("BXZ" , (char)ctr2) ) continue;
		if ( ! node->population[ctr2] ) continue;
		p = (double)node->population[ctr] * node->population[ctr2]*norm;
		if ( ! pair_prob[blosum_ctr][ctr][ctr2] ) {
		    printf ( "blah.\n");
		}
		node->entropy -=  p*log(p/pair_prob[blosum_ctr][ctr][ctr2] );
	    }
	}
   }

    return 0;
    
}
    
/********************************************************************************/
int hybrid_pheno (Alignment * alignment, double *score, Tree * tree){
    
    int pos,  no_seqs = alignment->number_of_seqs;
    int rank, g;
    double subsum, rho;
    int pheno_recursive (Node *node, int pos) ;
    
    /* calculate the score */
    for ( pos=0; pos < alignment->length; pos++) {
	pheno_recursive(tree->root, pos);
	rho = 0.0;
	for ( rank = 1; rank <no_seqs ; rank++ ) {
	    subsum = 0.0;
	    for ( g=0; g<rank; g++) {
		subsum  +=  tree->group_root[rank][g]->entropy;
	    }	    
	    rho += subsum/(rank);
	}
 	score[pos] =  rho;
    }
    return 0;
}

/********************************************************************************/
int pheno_over_nodes_1 (Alignment * alignment, double *score, Tree * tree){
    int pos,  no_seqs = alignment->number_of_seqs;
    int node_ctr;
    double  rho;
    int pheno_recursive (Node *node, int pos) ;
    
    /* calculate the score */
    for ( pos=0; pos < alignment->length; pos++) {
	pheno_recursive(tree->root, pos);
	rho = 0.0;
	for (node_ctr = no_seqs; node_ctr <  2*no_seqs-1; node_ctr++) {
	    rho += (tree->leaf+node_ctr)->entropy;
	}
 	score[pos] =  rho;
    }
    return 0;
}
/********************************************************************************/
int pheno_over_nodes_2 (Alignment * alignment, double *score, Tree * tree){
    int pos,  no_seqs = alignment->number_of_seqs;
    int node_ctr, id;
    double  rho;
    int pheno_recursive (Node *node, int pos) ;
    
    /* calculate the score */
    for ( pos=0; pos < alignment->length; pos++) {
	pheno_recursive(tree->root, pos);
	rho = 0.0;
	for (node_ctr = no_seqs; node_ctr <  2*no_seqs-1; node_ctr++) {
	    id = (tree->leaf+node_ctr)->id;
	    rho += (tree->leaf+node_ctr)->entropy/id;
	    //printf ( "%8d %4d   %6.2lf    %4d    %6.2lf  \n",pos, node_ctr, (tree->leaf+node_ctr)->entropy,
	    // (tree->leaf+node_ctr)->id, (tree->leaf+node_ctr)->entropy/(tree->leaf+node_ctr)->id);
	}
 	score[pos] =  rho;
    }

    return 0;
}
/***********************************************************************/
int pheno_over_nodes (Alignment * alignment, double *score, Tree * tree){
    return pheno_over_nodes_2 ( alignment, score, tree);
}

/***********************************************************************/
# define MAX_BLOSUM_INDEX 210

int pheno_init (char * table_file) {
    
    char *amino_acid_order = "ARNDCQEGHILKMFPSTWYV";
    char comment_char;
    char line[LONGSTRING];
    char token[MAX_TOK][MEDSTRING] = {{'\0'}};
    double blosum_probs[15][MAX_BLOSUM_INDEX];
    int  max_token, tok_ctr, retval;
    int line_ctr;
    int i, j, ctr, blosum_ctr;
    int char_i, char_j;
    int aao_strlen = strlen(amino_acid_order);
    FILE * fptr;
    int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING], char * fmt, char * warnstr) ;
    /* open table file */
    if ( !(fptr = efopen (table_file, "r" )) ) return 1;
    memset ( line, 0, LONGSTRING);
    ctr = blosum_ctr =  line_ctr = 0;
    while(fgets(line,LONGSTRING,fptr)!=NULL){
	line_ctr++;
	retval = tokenize ( token, &max_token, line, comment_char= '!' );
	switch ( retval ) {
	case  TOK_TOOMNY:
	    errmsg ( stderr, line_ctr, line, "\t\t %s\n", "Too many tokens.");
	    return 1;
	    break;
	case TOK_TOOLONG:
	    errmsg ( stderr, line_ctr, line, "\t\t %s\n", "Token too long.");
	    return 1;
	    break;
	}
	if ( max_token < 0 ) continue;
	if (  ! strncmp (token[0], "check", 5)  ) continue;
	if (  ! strncmp (token[0], "blosum", 6)  ) {
	    if ( max_token < 1 ) {
		fprintf ( stderr, "Error reading blosum tables.\n");
		return 1;
	    }  else {
		blosum_ctr = atoi ( token[1]);
		blosum_ctr = (blosum_ctr - 30)/5;
		ctr = 0;
	    }
	} else {
	    for (tok_ctr=0; tok_ctr <= max_token; tok_ctr ++ ) {
		if ( ctr >=  MAX_BLOSUM_INDEX) {
		    fprintf ( stderr, "Error reading blosum tables (line %d).\n", line_ctr);
		    return 1;
		}
		blosum_probs[blosum_ctr][ctr] = atof( token[tok_ctr] );
		ctr++;
	    }
	}
	
    }
    
    fclose (fptr);
    
    /*get the prob matrix to a usable form */
    for (blosum_ctr=0; blosum_ctr<16; blosum_ctr++) {
	ctr = 0;
	for(i=0;i<aao_strlen;i++){
	    char_i = (int) amino_acid_order [i];
	    for (j=0;j<=i;j++){
		char_j = (int) amino_acid_order [j];
		pair_prob[blosum_ctr][char_i][char_j] = pair_prob[blosum_ctr][char_j][char_i] = blosum_probs[blosum_ctr][ctr];
		ctr++;
	    }
	}
    }
    return 0;
}

/*********************************************************************************************/
