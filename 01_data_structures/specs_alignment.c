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

/************************************************************/

int process_almt (Options *options, Alignment *alignment) {
    
    int retval;
    int process_exons (Options *options, Alignment * alignment);
    int count_gaps    (Options *options, Alignment * alignment);
    int protected_positions (Options *options, Alignment * alignment);
    int seq_pw_dist(Alignment * alignment);

    /* store the position of exons, and replace them with gaps in the alignment */
    process_exons (options, alignment);

    /* gaps */
    count_gaps (options, alignment);

    /* protected positions */
    protected_positions (options, alignment);

    alignment->seq_dist = NULL;

    /*allocate space for various indicators of sequence similarity*/
    alignment->seq_dist =
	dmatrix ( alignment->number_of_seqs, alignment->number_of_seqs);
    
    alignment->aligned_sites =
	intmatrix ( alignment->number_of_seqs, alignment->number_of_seqs);
    if ( ! alignment->aligned_sites ) return 1;
    
    alignment->identical_sites =
	intmatrix ( alignment->number_of_seqs, alignment->number_of_seqs);
    if ( ! alignment->identical_sites ) return 1;
    
    alignment->similar_sites =
	intmatrix ( alignment->number_of_seqs, alignment->number_of_seqs);
    if ( ! alignment->similar_sites ) return 1;
    
 
 

    
    retval   = seq_pw_dist (alignment);
    if ( retval) return retval;

  
    return 0;
}

/*****************************************************************/
int process_exons (Options *options, Alignment * alignment) {

    int ctr;
    
    /* check if we have exons marked in the sequence */
    int exon_bdry_found = 0;
    for (ctr=0; ctr < alignment->number_of_seqs; ctr++ ) {
	if ( strchr (alignment->sequence[ctr], 'Z') ) {
	    exon_bdry_found = 1;
	    break;
	}
    }

    if (  !exon_bdry_found )
	return 0; // we're done here
    
    /* how many exons are there in each reference sequence? */
    int refseq_ctr;
    alignment->number_of_exon_bdries = emalloc (options->no_refseqs*sizeof(int) );
    if (!alignment->number_of_exon_bdries) return 1;
    
    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  {
	char * seq      = alignment->refseq[refseq_ctr];
	char * last_pos = alignment->refseq[refseq_ctr] + alignment->length;
	int z_count = 0;
	while (  seq<last_pos &&  (seq = strchr (seq,'Z')) )   {
	    seq ++;
	    z_count++;
	}
	alignment->number_of_exon_bdries[refseq_ctr] = z_count;
    }
       
    /* allocate space to store exon bdries */
    if (! (alignment->exon_bdry = chmatrix(options->no_refseqs, alignment->length) ) ) return 1;
    
    /* store the exon positions in a separate storage place */
    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  {
	char * seq        = alignment->refseq[refseq_ctr];
	int pos;
	for (pos=0; pos < alignment->length; pos++) {
	    alignment->exon_bdry[refseq_ctr][pos] = (seq[pos] == 'Z');
	}
    }
    /* in all sequences replace Z by '.', which I am unfortunately using throughout the code for the gap*/
    for (ctr=0; ctr < alignment->number_of_seqs; ctr++ ) {
	char * seq      = alignment->sequence[ctr];
	char * last_pos = alignment->sequence[ctr] + alignment->length;
	while (  seq<last_pos &&  (seq = strchr (seq,'Z')) )   {
	    (*seq) = '.';
	    seq ++;
	}
    }	
    return 0;
}

/*****************************************************************/
int protected_positions (Options *options, Alignment * alignment) {

    int s, c;
    alignment->protected_position = (int *) emalloc (alignment->length*sizeof(int));
    if (!alignment->protected_position) return 1;

    if (alignment->no_refseqs) {
	alignment->number_of_protected_positions = alignment->length;
	for ( c=0; c<alignment->length; c++) {
	    alignment->protected_position[c] = 1;
	}
	
    } else {

	alignment->number_of_protected_positions = 0;
	for ( c=0; c<alignment->length; c++)  {
	    for (s=0; s<alignment->no_refseqs; s++) {
		if ( alignment->refseq[s][c] != '-' ) {
		    alignment->protected_position[c] = 1;
		    alignment->number_of_protected_positions ++;
		    break;
		}
	    } 
	}
    }
    
    return 0;
}

/*****************************************************************/
int count_gaps (Options *options, Alignment * alignment) {

    int s, c;
    alignment->seq_gaps    = (int *) emalloc (alignment->number_of_seqs*sizeof(int));
    if (!alignment->seq_gaps) return 1;
    alignment->column_gaps = (int *) emalloc (alignment->length*sizeof(int));
    if (!alignment->column_gaps) return 1;
    
    for ( s=0; s<alignment->number_of_seqs; s++ ) {
	for ( c=0; c<alignment->length; c++) {
	    if ( alignment->sequence[s][c] == '.' ) {
		alignment->column_gaps[c] ++;
		alignment->seq_gaps[s] ++;
	    }
	}
    }
    alignment->refseq_gaps = 0;
    if (options->no_refseqs) {
	int pos = 0;
	for ( c=0; c<alignment->length; c++) {
	    if ( alignment->refseq[0][c] == '.' ) {
		alignment->refseq_gaps ++;
	    } else {
		pos ++;
	    
	    }
	}
    }
    
    return 0;
}
