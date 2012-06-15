# include "specs.h"


typedef struct {
    int seqno;
    int aligned;
    double pct_id;
    double pct_sim;
} Dist_descr;


int patch_almt ( Options *options, Alignment * alignment) {

    int i,j, k, p, no_seqs = alignment->number_of_seqs;
    char logname [BUFFLEN];
    int pos, patched;
    double pct_id, pct_sim;
    char *seq_i, *seq_j;
    char ** new_sequence;
    Dist_descr *nbr;
    FILE * log;
    
    
    sprintf (logname, "%s.patchlog", options->outname);
    if ( ! (log = efopen( logname, "w"))  ) return 1;
    
    /* for each sequence, except the reference sequence, */
    /* sort the remaining sequences in the order of preference  */
    /* for patching */
    if ( ! (nbr=emalloc (no_seqs*sizeof (Dist_descr) )) ) return 1;
    new_sequence = chmatrix (alignment->number_of_seqs, alignment->length);
    if ( !new_sequence ) return 1;
  
    memcpy ( new_sequence[0],  alignment->sequence[0],
	     alignment->number_of_seqs*alignment->length*sizeof(char) );
   
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

	
	
	
	nbr[0].seqno = -1; nbr[0].pct_id = 0.0; nbr[0].pct_sim = 0.0;
	
	for (j=0; j<no_seqs; j++) {
	    if (i==j) continue;
	    
	    pct_id  = (double)alignment->identical_sites[j][i]/alignment->aligned_sites[i][j];
	    pct_sim = (double)alignment->similar_sites[j][i]/alignment->aligned_sites[i][j];

	    /* find place in the nbr array */
	    k = 0;
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
		nbr[k].aligned = alignment->aligned_sites[i][j];
	    }
	}

	seq_i = alignment->sequence[i];
	fprintf (log, "%s, patched with sim cutoff %6.2lf\n",
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
		if (nbr[k].aligned < options->patch_min_length*alignment->length) continue;
		j = nbr[k].seqno;
		if ( j < 0) break; /* how and why was this supposed to happen? */
		seq_j = alignment->sequence[j];
		if ( seq_j[pos] != '.' ) {
		    new_sequence[i][pos] = seq_j[pos];
		    patched = 1;
		    fprintf (log, "\tpos %3d from %s (aligned length: %3d   sim: %6.3lf   id: %6.3lf)\n",
			     pos+1, alignment->name[j], nbr[k].aligned, nbr[k].pct_sim, nbr[k].pct_id);
		    break;
		}
	    }
	}
	if ( !patched )  fprintf (log, "\t not patched\n");

    }


   
    free_cmatrix (alignment->sequence);
    alignment->sequence = new_sequence;


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


    
    
    free (nbr);
    fclose (log);
    
    return 0;
}

