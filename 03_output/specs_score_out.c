# include "specs.h"


int  output_score ( Options * options, Protein * protein, Alignment * alignment,
		    int * almt2prot, double ** score,  int **res_rank, int surface){

    int almt_pos, score_ctr, length, seqctr, ctr, refseq_ctr;
    int pos[options->no_refseqs], restrict2structure = almt2prot && options->restrict2structure;
    double freq[ASCII];
    int printed[ASCII];
    char aux_str[22]; /*20 amino acid ypes + gap + newline */
    FILE * fptr;
    char filename[BUFFLEN], pdbid[PDB_ATOM_RES_NO_LEN+1] = {'\0'};
    double cvg;
    char aa = '\0';

    sprintf (filename, "%s.score", options->outname);

    fptr = efopen (filename, "w");
    if ( !fptr) return 1;
    
    fprintf (fptr, "%%%6s %4s", "almt", "pdb ");
    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  {
	fprintf (fptr, " %15s ",  options->refseq_name[refseq_ctr]);
    }
    fprintf (fptr, "%8s", "gaps");

    
    for ( score_ctr=0; score_ctr<options->number_of_methods; score_ctr++) {
	if ( options->raw_score ) {
	    fprintf (fptr, " %15s", options->method_name[score_ctr]);
	} else {
	    fprintf (fptr, " %6s", options->method_name[score_ctr]);
	}
    }
    if ( surface ) {
      fprintf (fptr, "%6s", "surf");
    }
    fprintf (fptr, "   substitutions\n");

    if ( surface ) {
	int resctr;
	length = 0;
	for ( resctr=0; resctr< protein->length; resctr++) {
	    length += protein->sequence[resctr].solvent_accessible;
	}
    } else {
	if ( restrict2structure) {
	    length = protein->length;
	} else {
	    length = alignment->length - alignment->refseq_gaps;
	}
    }

    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  pos[refseq_ctr] = -1;
    
    for (almt_pos = 0; almt_pos < alignment->length; almt_pos++) {

	for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)
	    if (alignment->refseq[refseq_ctr][almt_pos] != '.' )  pos[refseq_ctr]++;
	
	/*  sequential id, pdb id, aa type, frac of gaps */
	if ( almt2prot && almt2prot[almt_pos] >= 0 ) {
	    sprintf (pdbid, "%s", protein->sequence[ almt2prot[almt_pos] ].pdb_id );
	} else {
	    sprintf (pdbid, "%s", "-");
	}
	fprintf (fptr, "%6d %5s", almt_pos+1, pdbid);

	
	for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++) {
	    if ( almt2prot && almt2prot[almt_pos] >= 0 ) {
		aa = protein->sequence[ almt2prot[almt_pos] ].res_type_short;
	    } else {
		aa = alignment->refseq[refseq_ctr][almt_pos];
	    }
	    fprintf (fptr, "%6c",  aa);
	}
	
	fprintf (fptr, "%8.2lf", (double)alignment->column_gaps[almt_pos]/alignment->number_of_seqs);

	/* scores */
	for ( score_ctr=0; score_ctr<options->number_of_methods; score_ctr++) {
	    if ( restrict2structure ) {
		cvg =  ( almt2prot[almt_pos] >= 0 ) ? 
		(double)res_rank[score_ctr] [ almt2prot[almt_pos] ]/length : 1;
	    } else {
		/* the first refseq is especially dear to our heart: */
		if (options->no_refseqs && alignment->refseq[0][almt_pos] == '.' ) {
		    cvg = 1;
		} else {
		    cvg = (double)res_rank[score_ctr][pos[0]]/length;
 		}
	    }
	    if (options->raw_score) {
		fprintf (fptr, " %8.3lf %6.2lf", score[score_ctr][almt_pos], cvg);
	    } else {
		fprintf (fptr, " %6.2lf",  cvg);
	    }
	}
	
	/* surface accessibility */
	if ( surface && almt2prot ) {
	  int acc =  ( almt2prot[almt_pos] < 0 ) ?  
	    -1: protein->sequence[almt2prot[almt_pos]].solvent_accessible;
	    fprintf ( fptr," %6d",  acc);
	}

	/* substitutions seen in the alignment */
	memset ( freq, 0, ASCII*sizeof(double));
	for (seqctr=0; seqctr < alignment->number_of_seqs; seqctr++) {
	    if ( strchr ( "XZB", alignment->sequence[seqctr][almt_pos])) continue;
	    freq [ (int) alignment->sequence[seqctr][almt_pos] ] ++;
	}
	memset ( printed, 0, ASCII*sizeof(int));
	memset ( aux_str, 0, 22);
	ctr = 0;
	for (seqctr=0; seqctr < alignment->number_of_seqs; seqctr++) {
	    if ( strchr ( "XZB", alignment->sequence[seqctr][almt_pos])) continue;
	    if ( printed[ (int) alignment->sequence[seqctr][almt_pos] ] ) continue;
	    aux_str[ctr] =  alignment->sequence[seqctr][almt_pos];
	    printed[ (int) alignment->sequence[seqctr][almt_pos] ] = 1;
	    ctr ++;
	}
	fprintf (fptr, "  %s", aux_str);

	
	fprintf (fptr, "\n");
    }

    fclose (fptr);

    return 0;
}
