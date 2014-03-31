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

int main ( int argc, char * argv[]) {

    Options options;
    Protein protein;
    Alignment alignment;
    Tree tree;
    int    *prot2almt = NULL, *almt2prot = NULL;
    int retval;
    int surface = 0;
    int * epitope = NULL;
    int patch_almt ( Options *options, Alignment * alignment);
    int scoring_module (Options *options, Alignment * alignment, Tree * tree,
			Protein * protein, int * prot2almt, int * almt2prot,
			int  surface, int * epitope);    
    /* command file is required */
    if ( argc < 2 ) {
	fprintf ( stderr, "Usage: %s <command file>.\n", argv[0]);
	exit (0);
    }
    retval = read_cmd_file ( argv[1], &options);
    if (retval) exit(retval);
    retval = logger (&options, INTRO, "");
    if (retval) exit(retval);
   

    /*******************************************/
    /*                                         */
    /*  alignment input, tree building         */
    /*                                         */
    /*******************************************/

    /* read in the alignment */
    retval = read_alignment (&options, &alignment);
    if (retval) exit(retval);
    /* find  number of gaps, seq distances*/
    retval =  process_almt(&options, &alignment);
    if (retval) exit(retval);
    /* "patch" the alignment */

    if (options.patch_sim_cutoff > -1) {
	if (patch_almt (&options, &alignment)) return 1;
	/* output the patched alignment (afa will do)*/
	/* output gap as a "-" to make seaview happy */
	afa_out (&options, &alignment);
    }

    /* build the seq similarity tree */
    memset (&tree, 0, sizeof(Tree));

    if ( options.tree_file_name[0] )  {
	retval  = read_tree(&options,  &alignment, &tree);
	if (retval) return retval;
	
    } else {

	retval  = build_tree(&options, &alignment, &tree);
	if (retval) return retval;
 
    }
  
    retval = output_tree ( options.outname, tree.root);
    if (retval) exit(retval);


    /*******************************************/
    /*                                         */
    /*  mapping on the structure, if provided  */
    /*                                         */
    /*******************************************/
    if ( options.pdbname[0]) {
	/* warn if no chain given */
	if ( !options.chain) {
	    retval = logger (&options, WARN,
			     "No chain specified. Using the first one.");
	    if ( retval) exit (1);
	}
	if (retval) exit(retval);
	/* read in the structure */
	retval = read_pdb (options.pdbname, &protein, options.chain);
	if (retval) exit(retval);

	/* find mapping between the structure and the alignment*/
	if ( ! (prot2almt = (int *) emalloc (protein.length*sizeof(int))) ) exit (1);
	if ( ! (almt2prot = (int *) emalloc (alignment.length*sizeof(int))) )exit (1);
	retval = struct_almt_mapping (&protein, &alignment, prot2almt, almt2prot);
	if (retval) exit(retval);
	
	/*  read in the epitope, if provided       */
	if ( options.epi_file_name[0]) {
	    if ( ! (epitope = (int *) emalloc (protein.length*sizeof(int))) )exit (1);
	    retval  =  read_epitope (&options, &protein, epitope);
	    if (retval) exit(retval);
	} else {
	    epitope = NULL;
	}
    }

    /*******************************************/
    /*                                         */
    /*  surface, if provided                   */
    /*                                         */
    /*******************************************/
    if ( options.dssp_file_name[0] ) {
	/* surface = 0; */
        /* let's not couple the two - rather, the surface
	   could be made protected, but that's yet to be implemented*/
	retval  = read_dssp (&options, &protein);
	if (retval) exit(retval);
    }
    
    
    /*******************************************/
    /*                                         */
    /*         conservation scores             */
    /*                                         */
    /*******************************************/
    if ( options.path ) {
	double ** per_node_prob;  
	int ** description;
	Node ** path_to_refseq;
	path_to_refseq = emalloc (alignment.number_of_seqs*sizeof(Node*) );
	if (! path_to_refseq)  exit (1);
	if (!(per_node_prob = dmatrix  (alignment.length,
					   alignment.number_of_seqs)))  exit (1);
	if (!(description   = intmatrix  (alignment.length,
					   alignment.number_of_seqs)))  exit (1);
	path_determinants (options.significance_cutoff, &alignment,  &tree,
			   path_to_refseq, per_node_prob, description);
	if (retval) exit(retval);
	output_path_determinants (&options, &protein, &alignment, almt2prot,
				  path_to_refseq,  per_node_prob, description);
	if (retval) exit(retval);
    }

    if ( options.do_scoring ) {
	retval = scoring_module (&options, &alignment, &tree, &protein,
				 prot2almt, almt2prot, surface, epitope);
	if (retval) exit(retval);
    }
    logger ( &options, NOTE, "");

    return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

/* scoring */

int scoring_module (Options *options, Alignment * alignment,
		    Tree * tree, Protein * protein, int * prot2almt,
		    int * almt2prot, int  surface, int * epitope ) {

    double ** score, ** sensitivity, ** specificity;
    double* clustering_score; /* will be used only if pdb given*/
    char outfilename[BUFFLEN] = {'\0'};
    int ** res_rank, **int_cvg;
    int retval;
    int similar_to[ASCII] = {'\0'};
    int score_ctr;
    int number_of_scores = options->number_of_methods;
    /* set up the similarity array */
    set_similarity ( options, similar_to );
    
    if (! (score       = dmatrix  (number_of_scores, alignment->length)))  exit (1);
    if (! (res_rank    = intmatrix(number_of_scores, alignment->length)))  exit (1);
    if (! (int_cvg     = intmatrix(number_of_scores, alignment->length)))  exit (1);
    if (! (sensitivity = dmatrix  (number_of_scores, alignment->length)))  exit (1);
    if (! (specificity = dmatrix  (number_of_scores, alignment->length)))  exit (1);
    
    if ( options->pdbname[0]) {
	clustering_score  =  (double*) emalloc (protein->length*sizeof(double));
	if (!clustering_score) exit(1);
    }
	
    for (score_ctr = 0; score_ctr< options->number_of_methods; score_ctr++) {

	if ( options->method[score_ctr] == EXTERN ) {
            /*  thrid party score, if provided         */
	     retval = read_score_almt ( options->scorename, alignment, score[score_ctr] ); 
	} else if ( options->method[score_ctr] == DET ) {
	     retval = determinants (options->significance_cutoff, alignment, tree,
				    options->queried_node, score[score_ctr]);
	} else {
	     retval = scoring  (options, alignment,  tree, similar_to,
				score[score_ctr], score_ctr);
	}
	if (retval) exit(retval);
	
	coverage (protein, alignment, almt2prot, score[score_ctr], alignment->length,
		  res_rank[score_ctr], int_cvg[score_ctr], options->restrict2structure, surface );
	
	if ( options->pdbname[0]) {
	    epitope_detection_specificity ( protein, res_rank[score_ctr],
					    int_cvg[score_ctr], surface, epitope,
					    sensitivity[score_ctr], specificity[score_ctr] );
	    /* find the clustering score */
	    if ( ! epitope ) { /* I am hacking if I'm looking for
				  epitope finding specificity - residues
				  scored only on surface => clustering meaningless) */
		if ( ! options->restrict2structure ) {
		    printf ("Note: in this version need restrict2structure to do clustering."
			    "(This should be fixed.)\n");
		} else {
		    clustering ( protein,  res_rank[score_ctr],
				 int_cvg[score_ctr], clustering_score);

		    /* clustering info */
		    memset  (outfilename, 0, BUFFLEN*sizeof(char));
		    sprintf (outfilename, "%s.%s",options->outname,
			     options->method_name[score_ctr]);
		    retval = output_clustering (outfilename, int_cvg[score_ctr],
						clustering_score, protein->length );
		    if (retval) exit(retval);
		}
	    }
	}
    }
    

    /***********************************/
    /* output                          */
    retval = output_score (options, protein, alignment, almt2prot,
			   score, res_rank, surface);
    if (retval) return (retval);
    
    if ( epitope ) {
	retval = output_epi_spec (options,  protein, int_cvg,
				  sensitivity, specificity, surface);
    }
    if (retval) return (retval);

    
    return 0;
}
