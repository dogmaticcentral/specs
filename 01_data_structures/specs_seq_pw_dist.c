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

/* pairwise distances between sequences */

/*********************************************************************************************/
/*********************************************************************************************/

/*********************************************************************************************/
/*********************************************************************************************/

int seq_pw_dist (Alignment * alignment) {
    char *amino_acid_order = "ABCDEFGHIKLMNPQRSTVWXYZ";

    int blosum62[]={
	4,
	-2,  4,
	0, -3,  9,
	-2,  4, -3,  6,
	-1,  1, -4,  2,  5,
	-2, -3, -2, -3, -3,  6,
	0, -1, -3, -1, -2, -3,  6,
	-2,  0, -3, -1,  0, -1, -2,  8,
	-1, -3, -1, -3, -3,  0, -4, -3,  4,
	-1,  0, -3, -1,  1, -3, -2, -1, -3,  5,
	-1, -4, -1, -4, -3,  0, -4, -3,  2, -2,  4,
	-1, -3, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5,
	-2,  3, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6,
	-1, -2, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7,
	-1,  0, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,
	-1, -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5,
	1,  0, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,
	0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,
	0, -3, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4,
	-3, -4, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,
	0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1,  0,  0, -1, -2, -1,
	-2, -3, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2, -1,  7,
	-1,  1, -3,  1,  4, -3, -2,  0, -3,  1, -3, -1,  0, -1,  3,  0,  0, -1, -2, -3, -1, -2,  4};
    
     int i,j,pos, ctr, off_diag;
    int aao_strlen;
    int number_of_similar_sites;
    int number_of_common_sites;
    int number_of_id_sites;
    int similarity[128][128];
    char * seq_i, *seq_j;
    int char_i, char_j;
    double avg;


    /*get the similarity matrix to a usable form */
    aao_strlen = strlen(amino_acid_order);
    
    ctr = 0;
    avg = 0;
    off_diag = 0;
    for(i=0;i<aao_strlen;i++){
	char_i = (int) amino_acid_order [i];
	for (j=0;j<=i;j++){
	    char_j = (int) amino_acid_order [j];
	    similarity[char_i][char_j] = similarity[char_j][char_i] = blosum62[ctr];
	    if ( i != j  && blosum62[ctr] > 0) {
		avg +=  blosum62[ctr];
		off_diag ++;
	    }
	    ctr++;
	}
    }
    
    avg /= off_diag;
    avg = ceil(avg);
    
    int length_i, length_j, shorter;

    /* sanity */
    for (i=0; i<alignment->number_of_seqs; i++) {
	seq_i = alignment->sequence[i];
	length_i = alignment->length - alignment->seq_gaps[i];
	if ( ! length_i ) {
	    fprintf (stderr, "Error in seq_pw_dist(): %s is all gaps (?!).\n",
		     alignment->name[i]);
	    exit (1);
	}
    }

    
    for (i=0; i<alignment->number_of_seqs; i++) {
	seq_i = alignment->sequence[i];
	length_i = alignment->length - alignment->seq_gaps[i];
	for (j=i+1; j<alignment->number_of_seqs; j++) {
	    length_j = alignment->length - alignment->seq_gaps[j];
	    seq_j = alignment->sequence[j];
	    
	    number_of_common_sites  = 0;
	    number_of_id_sites      = 0;
	    number_of_similar_sites = 0;
	    
	    for (pos=0; pos<alignment->length; pos++) {
		
		if ( seq_i[pos] == '.' ) continue;
		if ( seq_j[pos] == '.' ) continue;
		
		number_of_common_sites ++;
		if ( seq_i[pos] == seq_j[pos] ) number_of_id_sites++;
		if ( similarity [(int)seq_i[pos]] [(int) seq_j[pos] ] >= avg ) {
		    number_of_similar_sites ++;
		}
	    }
	    shorter = (length_i<length_j) ? length_i : length_j;
	    alignment->seq_dist[i][j] = 1 - (double) number_of_similar_sites/shorter;
	    /* alignment->seq_dist[i][j] = 1 - (double) number_of_similar_sites/number_of_common_sites; */
	    alignment->seq_dist[j][i] = alignment->seq_dist[i][j];

	    alignment->aligned_sites[i][j]   =
		alignment->aligned_sites[j][i]   =  number_of_common_sites;
	    alignment->identical_sites[i][j] =
		alignment->identical_sites[j][i] =  number_of_id_sites;
	    alignment->similar_sites[i][j]   =
		alignment->similar_sites[j][i]   =  number_of_similar_sites;
	    
	}
    }
# if 0
    for (i=0; i<alignment->number_of_seqs; i++) {
	for (j=i+1; j<alignment->number_of_seqs; j++) {
	    printf (" %4d  %4d  %10s %10s  %8.3lf\n",
		    i+1, j+1, alignment->name[i],  alignment->name[j], alignment->seq_dist[i][j]);
	}
    }
	exit (0);
# endif   
    return 0;
}
