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

    /******************************************************************************/
    /* normalization for the fractional coverage (= fract coverage  of what)      */
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
	    length = alignment->number_of_protected_positions;
	}
    }
    /******************************************************************************/
    /* print out the header                                                       */
    fprintf (fptr, "%%%6s", "almt");
    
    fprintf (fptr, "%8s", "gaps(%)");
    
    for ( score_ctr=0; score_ctr<options->number_of_methods; score_ctr++) {
	if ( options->raw_score ) {
	    fprintf (fptr, " %15s", options->method_name[score_ctr]);
	} else {
	    fprintf (fptr, " %6s", options->method_name[score_ctr]);
	}
    }

    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  {
	fprintf (fptr, " %15s ",  options->refseq_name[refseq_ctr]);
	fprintf (fptr, " pos_in_%s ",   options->refseq_name[refseq_ctr]);
    }

    
    if ( almt2prot ) {
	fprintf (fptr, " %8s ", "pdb_aa");
	fprintf (fptr, " %8s ", "pdb_id");
    }

    if ( options->dssp_file_name[0] ) {
      fprintf (fptr, "%6s", "surf");
    }
    fprintf (fptr, "   substitutions\n");


    /******************************************************************************/
    /* print out the rest                                                         */
    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  pos[refseq_ctr] = -1;
    
    for (almt_pos = 0; almt_pos < alignment->length; almt_pos++) {

	for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)
	    if (alignment->refseq[refseq_ctr][almt_pos] != '.' )  pos[refseq_ctr]++;
	
	/*  sequential id*/
	fprintf (fptr, " %6d ", almt_pos+1);
	
	/* percent gaps in the alignment */
	fprintf (fptr, " %3d ", (int)(100*(double)alignment->column_gaps[almt_pos]/alignment->number_of_seqs));
	
	
	/* scores */
	for ( score_ctr=0; score_ctr<options->number_of_methods; score_ctr++) {
	    if ( restrict2structure ) {
		cvg =  ( almt2prot[almt_pos] >= 0 ) ? 
		(double)res_rank[score_ctr] [ almt2prot[almt_pos] ]/length : 1;
	    } else {
		/* the first refseq is especially dear to our heart: */
		if ( alignment->protected_position[almt_pos] && pos[0] >= 0 && pos[0] < alignment->length) {
		    cvg = (double)res_rank[score_ctr][pos[0]]/length;
 		} else {
		    cvg = 1;
		}
	    }
	    if (options->raw_score) {
		fprintf (fptr, " %8.3lf %6.2lf", score[score_ctr][almt_pos], cvg);
	    } else {
		fprintf (fptr, " %6.2lf",  cvg);
	    }
	}
	
	/* reference sequence(s) */
	for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++) {
	    aa = alignment->refseq[refseq_ctr][almt_pos];
	    if ( aa == '.' ) {
		fprintf (fptr, "%6c  %6c",  '.',  '.');
	    } else {
		fprintf (fptr, "%6c  %6d ",  aa, pos[refseq_ctr]+1);
	    }
	}
	

	/* pdb aa, and position (identifier ~= position) */
	if ( almt2prot ) {
	    if ( almt2prot[almt_pos] >= 0 ) {
		sprintf (pdbid, "%s", protein->sequence[ almt2prot[almt_pos] ].pdb_id );
		aa =  protein->sequence[ almt2prot[almt_pos] ].res_type_short;
	    } else {
	        sprintf (pdbid, "%s", ".");
		aa = '.';
	    }
	    fprintf (fptr, " %5c  %5s", aa, pdbid);
	}
	
	/* surface accessibility */
	if ( options->dssp_file_name[0] ) {
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
