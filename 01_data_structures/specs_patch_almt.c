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


typedef struct {
    int seqno;
    int aligned_length;
    double pct_id;
    double pct_sim;
} Dist_descr;


int patch_almt ( Options *options, Alignment * alignment) {

    int i,j, k, p, no_seqs = alignment->number_of_seqs;
    char logname [BUFFLEN];
    int pos, patched;
    int number_of_seqs_to_ignore = 0;
    int * ignore_seq;
    double pct_id, pct_sim;
    char *seq_i, *seq_j;
    char ** new_sequence;
    Dist_descr *nbr;
    FILE * log;
    
    
    sprintf (logname, "%s.patchlog", options->outname);
    if ( ! (log = efopen( logname, "w"))  ) return 1;

    if ( !(ignore_seq =  (int*) emalloc (alignment->number_of_seqs*sizeof(int)) ) ) return 1;
    
    
    /* for each sequence, except the reference sequence, */
    /* sort the remaining sequences in the order of preference  */
    /* for patching */
    if ( ! (nbr=emalloc (no_seqs*sizeof (Dist_descr) )) ) return 1;

    
    /* the default: the new seqs are the same old seqs */
    new_sequence = chmatrix (alignment->number_of_seqs, alignment->length);
    if ( !new_sequence ) return 1;
    memcpy ( new_sequence[0],  alignment->sequence[0],
	     alignment->number_of_seqs*alignment->length*sizeof(char) );

    /* now inspect sequences one by one - i runs over sequences (not over positions) */
    for (i=0; i<no_seqs; i++) {
	if ( !alignment->seq_gaps[i] ) continue;

	int skip = 0;
	if (options->no_refseqs) {
	    int refseq_ctr;
	    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++) {
		if (! strcmp (alignment->name[i], options->refseq_name[refseq_ctr]) ) {
		    skip = 1;
		    break;
		}
	    }
	}
	if ( skip) continue;


	/* careful with the pdbseq here: */
	if (  options->pdbseq_name && !options->skip_pdbseq ) {
	    for (p=0; p<no_seqs; p++) {
		if (! strcmp ( alignment->name[p], options->pdbseq_name) ) {
		    skip = 1;
		    break;
		}
	    }
	}
	if ( skip) continue;

	/* here we have a  sequence that has gaps, and is not protected from patching */
	/* find the closest neighbor to patch from                                    */
	nbr[0].seqno = -1; nbr[0].aligned_length = 0.0; nbr[0].pct_id = 0.0; nbr[0].pct_sim = 0.0;
	
	for (j=0; j<no_seqs; j++) {
	    if (i==j) continue;
	    int aligned_length = alignment->aligned_sites[i][j];
	    pct_id  = (double)alignment->identical_sites[j][i]/aligned_length;
	    pct_sim =   (double)alignment->similar_sites[j][i]/aligned_length;

	    /* find place in the nbr array */
	    k = 0;
	    while ( k<no_seqs && aligned_length <nbr[k].aligned_length) k++;
	    while ( k<no_seqs && pct_id <nbr[k].pct_id) k++;
	    while ( k<no_seqs && pct_sim<nbr[k].pct_sim) k++;

	    /* move no_seqs-1-k elements of the nbr array;
	       and write this sequence in position k */
	    if (  k < no_seqs) {
		if ( k < no_seqs-1 ) {
		    memmove ( nbr+k+1, nbr+k, (no_seqs-1-k)*sizeof(Dist_descr) ); 
		}
		nbr[k].seqno = j;
		nbr[k].pct_sim = pct_sim;
		nbr[k].pct_id  = pct_id;
		nbr[k].aligned_length = aligned_length;
	    }
	}

	seq_i = alignment->sequence[i];
	fprintf (log, "%s, patching with sim cutoff %6.2lf:\n",
		 alignment->name[i], options->patch_sim_cutoff);
		
	patched = 0;
	for (pos=0; pos < alignment->length; pos ++) {
	    /* we are not patching postns which are not gappped,
	        nor the positions which are gap for most of the alignment */
	    if (  seq_i[pos] != '.' ||
		  (double) alignment->column_gaps[pos]/alignment->number_of_seqs > options->max_gaps )
		continue;
	    /* patch from the closest relative possible */
	    for (k=0; k<no_seqs; k++) {
		if (nbr[k].pct_sim < options->patch_sim_cutoff) break;
		if (nbr[k].aligned_length < options->patch_min_length*alignment->length) continue;
		j = nbr[k].seqno;
		if ( j < 0) break; /* how and why was this supposed to happen? */
		seq_j = alignment->sequence[j];
		if ( seq_j[pos] != '.' ) {
		    new_sequence[i][pos] = seq_j[pos];
		    patched = 1;
		    fprintf (log, "\tpos %3d from %s (aligned length: %3d   sim: %6.3lf   id: %6.3lf)\n",
			     pos+1, alignment->name[j], nbr[k].aligned_length, nbr[k].pct_sim, nbr[k].pct_id);
		    break;
		}
	    }
	}
	if ( !patched )  {
	    fprintf (log, "\t not patched\n");
	    int nongaps_after_patching = 0;
	    double fract_length        = 1;
	    for (pos=0; pos < alignment->length; pos ++)  nongaps_after_patching += (new_sequence[i][pos]!='.');
	    fract_length = (double)nongaps_after_patching/alignment->length;
	    if ( fract_length <  options->min_fragment_length) {
		/* fraction of the alignment length that this seqeunce covers is small */
		/* is there a close neighbor that covers the full length though?       */
		/* (for example, if all fish seqeunces are shorter thatn 2/3 of the alignmen, we do not want to throw them out;
		   but, if one is a fragment while the other seems to be whole, we want to get rid of the fragment) */
		double fract_length_j = 0;
		for (k=0; k<no_seqs; k++) {
		    
		    if (nbr[k].pct_sim < 0.8) break;  /* <<<< have some hardcoded numbers here! */
		    j = nbr[k].seqno;
		    if ( j < 0) break; /* how and why was this supposed to happen? */
		    fract_length_j = (double)(alignment->length - alignment->seq_gaps[j])/alignment->length;
		    if (fract_length_j > 0.8) {
			
			fprintf (log, "\t fraction of the non-gapped length (%4.2f) is smaller than allowed (%4.2f)\n",
				 fract_length, options->min_fragment_length);
			fprintf (log, "\t however, its neighbor  %s  (fract sim: %4.2f) covers %4.2f of the alignment length \n",
				 alignment->name[j], nbr[k].pct_sim, fract_length_j);
			fprintf (log, "\t %s will be ignored in the further consideration\n", alignment->name[i]);
			ignore_seq[i] = 1;
			number_of_seqs_to_ignore ++;
			break;
		    }
		}
	    }
	}
    }

    if ( !number_of_seqs_to_ignore) {
	free_cmatrix (alignment->sequence);
	alignment->sequence = new_sequence;
	
    } else {
	
	/* what changes ?*/
	/* sequences themselves & their names*/
	int new_seq_ctr = 0;
	char **tmp_name  = (char **) emalloc(no_seqs*sizeof(char*));
	if (!tmp_name) return 1;
	
	memset (alignment->sequence[0], 0, no_seqs*alignment->length*sizeof(char));
	for (i=0; i < no_seqs; i++ ) {
	    tmp_name[i] = alignment->name[i];
	}
	for (i=0; i < no_seqs; i++ ) {
	    if ( ignore_seq[i] ) continue;
	    memcpy (alignment->sequence[new_seq_ctr], new_sequence[i], alignment->length*sizeof(char));
	    alignment->name[new_seq_ctr] = tmp_name[i];
	    new_seq_ctr ++;
	}
 	/* number of sequences in the alignment */
	alignment->number_of_seqs = new_seq_ctr;
	
	free_cmatrix (new_sequence);
	free (tmp_name);
    }
	
    /* make sure that alignment->refseq and alignment->pdbseq
       (which are aliases) point to the location in the new array */
    for (i=0; i < no_seqs; i++ ) {
	int refseq_i;
	for (refseq_i=0; refseq_i<options->no_refseqs; refseq_i++)  {
	    if (! strcmp ( alignment->name[i], alignment->refseq_name[refseq_i]) ) {
		alignment->refseq[refseq_i] = alignment->sequence[i];
		break;
	    }
	}
    }

    if (  options->pdbname[0] && !options->skip_pdbseq ) {
	for (i=0; i < no_seqs; i++ ) {
	    if (! strcmp ( alignment->name[i], options->pdbseq_name) ) {
		alignment->pdbseq = alignment->sequence[i];
		break;
	    }
	}
     }

     /* the number of gaps and the pairwise distances beteen seqs are now different */
     /* gaps */
     count_gaps (options, alignment);
     /* various indicators of sequence similarity */
     if (! seq_pw_dist (alignment)) return 1;

     
    
     free (nbr);
     fclose (log);
    
     return 0;
}

