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

int valdar ( Alignment * alignment, double *score){

    char *amino_acid_order = "ARNDCQEGHILKMFPSTWYV";
    int pet91[] = {
        2,
       -1, 5,
        0, 0, 3,
        0, -1, 2, 5,
       -1, -1, -1, -3, 11,
       -1, 2, 0, 1, -3, 5,
       -1, 0, 1, 4, -4, 2, 5,
        1, 0, 0, 1, -1, -1, 0, 5,
       -2, 2, 1, 0, 0, 2, 0, -2, 6,
        0, -3, -2, -3, -2, -3, -3, -3, -3, 4,
       -1, -3, -3, -4, -3, -2, -4, -4, -2, 2, 5,
       -1, 4, 1, 0, -3, 2, 1, -1, 1, -3, -3, 5,
       -1, -2, -2, -3, -2, -2, -3, -3, -2, 3, 3, -2, 6,
       -3, -4, -3, -5, 0, -4, -5, -5, 0, 0, 2, -5, 0, 8,
        1, -1, -1, -2, -2, 0, -2, -1, 0, -2, 0, -2, -2, -3, 6,
        1, -1, 1, 0, 1, -1, -1, 1, -1, -1, -2, -1, -1, -2, 1, 2,
        2, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 0, -2, 1, 1, 2,
       -4, 0, -5, -5, 1, -3, -5, -2, -3, -4, -2, -3, -3, -1, -4, -3, -4, 15,
       -3, -2, -1, -2, 2, -2, -4, -4, 4, -2, -1, -3, -2, 5, -3, -1, -3, 0, 9,
        1, -3, -2, -2, -2, -3, -2, -2, -3, 4, 2, -3, 2, 0, -1, -1, 0, -3, -3, 4};
    int i,j, ctr, diag;
    int char_i, char_j;
    int aao_strlen;
    double diag_avg, min, max;
    double similarity[128][128] = {{0.0}};
    
    /*get the similarity matrix to a usable form */
    aao_strlen = strlen(amino_acid_order);
    
    ctr = 0;
    diag_avg = 0;
    diag = 0;
    for(i=0;i<aao_strlen;i++){
	char_i = (int) amino_acid_order [i];
	for (j=0;j<=i;j++){
	    char_j = (int) amino_acid_order [j];
	    similarity[char_i][char_j] = similarity[char_j][char_i] = pet91[ctr];
	    if ( i == j ) {
		diag_avg +=  pet91[ctr];
		diag ++;
	    }
	    ctr++;
	}
    }
    diag_avg /= diag;
    for(i=0;i<128;i++) similarity[i][i] = diag_avg;
    
    max = -100; min = 100;
    for(i=0;i<128;i++){
	for(j=0;j<128;j++){
	    if ( max < similarity[i][j] ) max =  similarity[i][j] ;
	    if ( min > similarity[i][j] ) min =  similarity[i][j] ;
	}
    }
    for(i=0;i<128;i++){
	for(j=0;j<128;j++){
	    similarity[i][j] = (double)(similarity[i][j] - min )/(max-min);
	}
    }
    

    char *seq_i, *seq_j;
    int number_of_aligned_sites, pos;
    double sim_sum;
    double **seq_dist = dmatrix (alignment->number_of_seqs, alignment->number_of_seqs);
    if ( ! seq_dist) return 1;
    
    /* set sequence distance */
    for (i=0; i<alignment->number_of_seqs; i++) {
	seq_i = alignment->sequence[i];
	seq_dist[i][i] = 0.0;
	for (j=i+1; j<alignment->number_of_seqs; j++) {
	    seq_j = alignment->sequence[j];
	    number_of_aligned_sites = 0;
	    sim_sum = 0.0;
	    for (pos=0; pos<alignment->length; pos++) {
		if ( seq_i[pos] == '.' && seq_j[pos] == '.' ) continue;
		number_of_aligned_sites ++;
		if  (seq_i[pos] != '.' &&  seq_j[pos] != '.' ) {
		    sim_sum += similarity [(int)seq_i[pos]] [(int) seq_j[pos] ] ;
		}
	    }
	    seq_dist[i][j] = 1.0 - sim_sum/number_of_aligned_sites;
	    seq_dist[j][i] = seq_dist[i][j];
	    /*printf ( " %3d  %3d  %4.2lf   %4.2lf\n",
	      i, j, alignment->seq_dist[i][j], seq_dist[i][j]);*/
	    /* the resolution of this distance is rather poor  */
	}
    }

    /* sequence weight */
    double * seq_weight = emalloc (alignment->number_of_seqs *sizeof(double));
    double norm;
    if ( ! seq_weight) return 1;
    for (i=0; i<alignment->number_of_seqs; i++) {
	seq_weight[i] = 0.0;
	for (j=0; j<alignment->number_of_seqs; j++) {
	    seq_weight[i] += seq_dist[i][j];
	}
	seq_weight[i] /= (alignment->number_of_seqs-1);
    }

    free ( seq_dist[0] );
    free ( seq_dist);

    norm = 0.0;
    for (i=0; i<alignment->number_of_seqs; i++) {
	for (j=i+1; j<alignment->number_of_seqs; j++) {
	    norm += seq_weight[i]*seq_weight[j];
	}
    }

    norm = sqrt (norm);
    for (i=0; i<alignment->number_of_seqs; i++) {
	seq_weight[i] /= norm;
    }


    for (pos=0; pos<alignment->length; pos++) {
	score[pos] = 0.0;
    }

    double sim;
    for (i=0; i<alignment->number_of_seqs; i++) {
	seq_i = alignment->sequence[i];
	for (j=i+1; j<alignment->number_of_seqs; j++) {
	    seq_j = alignment->sequence[j];
	    for (pos=0; pos<alignment->length; pos++) {
		if ( seq_i[pos] != '.' && seq_j[pos] != '.' ) {
		    sim = similarity [(int)seq_i[pos]] [(int) seq_j[pos] ] ;
		    score[pos] += seq_weight[i]*seq_weight[j]*sim;
		}
	    }
	}

    }
    /* take the complement of Valdars's score, so
       it assigns 0 to the most conserved and 1 to the least */
    for (pos=0; pos<alignment->length; pos++) {
	score[pos]  = 1 - score[pos];
    }
    free ( seq_weight);
    
    return 0;
} 
