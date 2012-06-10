# include "specs.h"

 int output_path_determinants (Options * options,  Protein * protein, Alignment * alignment,
			       int * almt2prot, Node ** path_to_refseq,
			       double  **per_node_prob, int **description){

     FILE * fptr;
     char filename[BUFFLEN], pdbid[PDB_ATOM_RES_NO_LEN];   
     int path_ctr, almt_pos;
     char aa = '\0';
     int output_raw = 1;
     int output_flag[alignment->number_of_seqs];
     
     sprintf (filename, "%s.path", options->outname);

     fptr = efopen (filename, "w");
     if ( !fptr) return 1;

     /* don't output the nodes for which there is nothing
	statistically significant happening  - first figure
        out which nodes are in question */
     for (path_ctr=1; path_ctr < alignment->number_of_seqs
	      && path_to_refseq[path_ctr]; path_ctr++ ) {
	 output_flag[path_ctr] = 0;
	 for (almt_pos = 0; almt_pos < alignment->length; almt_pos++) {
	     if ( description[almt_pos][path_ctr] != INS ) {
		 output_flag[path_ctr] = 1;
		 break;
	     }
	 }
     }
     fprintf (fptr, "%%%6s %5s%5s%8s", "almt", "pdb ", "aa ",  "gaps");
     for (path_ctr=1; path_ctr < alignment->number_of_seqs
	      && path_to_refseq[path_ctr]; path_ctr++ ) {
	 if ( ! output_flag[path_ctr] ) continue;
	 if (output_raw) {
	     fprintf (fptr," %12d", path_to_refseq[path_ctr]->id);
	     
	 } else {
	     fprintf (fptr," %4d", path_to_refseq[path_ctr]->id);
	 }
     }
     fprintf (fptr, "\n");
     
     for (almt_pos = 0; almt_pos < alignment->length; almt_pos++) {

	
	/*  sequential id, pdb id, aa type, frac of gaps */
	if ( almt2prot && almt2prot[almt_pos] >= 0 ) {
	    sprintf (pdbid, "%s", protein->sequence[ almt2prot[almt_pos] ].pdb_id );
	    aa = protein->sequence[ almt2prot[almt_pos] ].res_type_short;
	} else {
	    sprintf (pdbid, "%s", "-");
	    aa = alignment->refseq[almt_pos];
	}
	fprintf (fptr, "%6d%6s%6c%8.2lf", almt_pos+1, pdbid, aa,
		 (double)alignment->column_gaps[almt_pos]/alignment->number_of_seqs);
	for (path_ctr=1; path_ctr < alignment->number_of_seqs
		 && path_to_refseq[path_ctr]; path_ctr++ ) {
	    
	    if ( ! output_flag[path_ctr] ) continue;
	    
	    if (output_raw) {
		switch ( description[almt_pos][path_ctr] ) {
		case INS:
		    /* fprintf (fptr," %4s", "INS"); */
		    fprintf (fptr,"  %6s %4s ", "-",  "-");
		    break;
		case LOF:
		    fprintf (fptr,"  %3.1le %4s",
			     per_node_prob[almt_pos][path_ctr], "LOF");
		    break;
		case GOF:
		    fprintf (fptr,"  %3.1le %4s",
			     per_node_prob[almt_pos][path_ctr], "GOF");
		    break;
		case DISCR:
		    fprintf (fptr,"  %3.1le %4s",
			     per_node_prob[almt_pos][path_ctr], "DIS");
		    break;
		default:
		    fprintf (fptr,"  %3.1le %4s",
			     per_node_prob[almt_pos][path_ctr], " ? ");
		}
	    } else {
		switch ( description[almt_pos][path_ctr] ) {
		case INS:
		    /* fprintf (fptr," %4s", "INS"); */
		    fprintf (fptr," %4s", "-");
		    break;
		case LOF:
		    fprintf (fptr," %4s", "LOF");
		    break;
		case GOF:
		    fprintf (fptr," %4s", "GOF");
		    break;
		case DISCR:
		    fprintf (fptr," %4s", "DIS");
		    break;
		default:
		    fprintf (fptr," %4s", " ? ");
		}
	    }
	    
	}
	fprintf (fptr, "\n");

     }
     

     return 0;

 }
