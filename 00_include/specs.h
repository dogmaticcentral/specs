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
/*******************************************/
/* by  Ivana Mihalek, 2006-2009            */
/*******************************************/

# ifndef _BC_H
# define _BC_H
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include <time.h>
# include "specs_pdb.h"
# include "specs_geometry.h"
# include "specs_node.h"
# include "specs_alignment.h"
# include "specs_utils.h"

# define VERSION "2012"
# define PROBS   "/home/ivanam/baubles/specs/08_data/probs"


/*****************************/
/*   log   modes             */
/*****************************/
# define INTRO 1
# define WARN  2
# define STATUS  3
# define NOTE  4

# define BUFFLEN  250
/******************************/
/*   tokenizer                */
/******************************/

# define TOK_TOOMNY  1 /* tokenazier error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */
# define LONGSTRING  250
# define MEDSTRING   100
# define SHORTSTRING  25

/******************************/
/*   geometry                 */
/******************************/
/* distance to be called a nighbor in the protein structure */
# define CUTOFF_DIST 4.0

/******************************/
/* tree building methods      */
/******************************/
# define UPGMA 1
# define NEIGHBOR_JOINING 2
# define CONSENSUS_UPGMA 3

/******************************/
/*   minimal interesting      */
/*     group of sequences     */
/******************************/
# define  MIN_NO_LEAVES 5
# define  MAX_GAPS  0.3

/******************************/
/*  path discriminants        */
/******************************/
# define SIGNIFICANCE_CUTOFF  1.e-4
typedef enum { INS, LOF, GOF, DISCR } Description; 

/******************************/
/*   user options:            */
/******************************/
/* please update echo options int the logger */
typedef struct {
    char pdbname [BUFFLEN];
    char almtname[BUFFLEN];
    char pdbseq_name [BUFFLEN];  /* the name of the pdb sequence in the alignment */
    char outname [BUFFLEN];
    char epi_file_name  [BUFFLEN];
    char dssp_file_name [BUFFLEN];
    char scorename [BUFFLEN]; /* file with the score already provided */
    char tree_file_name[BUFFLEN]; /* input the tree */
    int do_scoring;
    int path; /* determinants on the path to the query */
    int raw_score;
    int restrict2structure; /* output option: rerank the score so the "coverage"
			       refers to input structure only */
    int skip_pdbseq;  /* ship pdbseq when doing the scoring - so the scores can be
			 mapped on an unrelated structures */
    int number_of_methods;
    int *method;
    char ** method_name;
    int tree_method;
    int queried_node;
    int no_refseqs;
    char ** refseq_name; /* reference sequence(s) -in particular
			    - gaps in the refseq[0] do not count when calculating the coverage*/
    char chain;
    double max_gaps;
    double acc_cutoff;
    double significance_cutoff;
    double patch_sim_cutoff;
    double patch_min_length;
    
} Options; 

/* scoring methods */
typedef enum { ENTROPY,    VALDAR,    RVET,   RVET_NORM,  MAJORITY, PHENO,  IVET, CARBONE, 
	       ENTR_W_SIM, RV_W_SIM, IV_W_SIM, MAJORITY_W_SIM,  PHENO_HYBRID, PHENO_NODES,
	       EXTERN, DET} Method_type;
# define NR_IMPLEMENTED_METHODS  14 /* EXTERNal is not an implemented method */


/* number of rounds in the simulation*/
# define ROUNDS_PER_SITE 100

/******************************/
/* function declarations :    */
/******************************/
int  afa_out (Options * options, Alignment * alignment);
int  build_tree (Options * options, Alignment * alignment, Tree * tree);
int  clustering (Protein *protein,  int* res_rank, int * int_cvg,
		  double *clustering_score);
int  coverage (Protein * protein, Alignment * alignment,
		int * almt2prot, double * score,
		int almt_length, int * res_rank, int * int_cvg,
	       int restrict2structure, int surface);

int  determinants (double sig_cutoff, Alignment * alignment,  Tree * tree,
		  int queried_node_id, double *score);

int  determine_adj_matrix ( int ** adj_matrix, Residue * sequence,
			    int no_res, double cutoff_dist);
int  entropy_complement (Node *node, Node * root, int pos) ;
int  entropy_recursive (Node *node,  int pos, int * similar_to, int normalize);
int  epitope_detection_specificity ( Protein *protein,  int * res_rank,
				     int * int_cvg,
				    int  surface, int *epitope, 
				     double * sensitivity,  double * specificity);
int find_group_roots     (Node *node, Node*** group_root,  int rank, int * group_tag );
int find_group_roots_sim (Node *first_node, int no_of_nodes, int no_of_bins, Node*** group_root);
int  find_query (Alignment * alignment, char * query_name);
int  logger (Options * options, int mode, char *msg );
int  output_clustering (char * base_filename,  int * int_cvg,
			double *clustering_score, int length);
int  output_epi_spec ( Options * options, Protein * protein, int ** int_cvg,
		       double **sensitivity, double **specificity, int surface);

int  output_path_determinants (Options * options,  Protein * protein, Alignment * alignment,
			       int * almt2prot, Node ** path_to_query,
			       double  **per_node_prob, int **description);
int  output_score ( Options * options, Protein * protein, Alignment * alignment,
		    int * almt2prot, double ** score,  int **res_rank, int surface);

int  output_specs_score ( Options *options, Protein * protein, Alignment * alignment,
			  int *almt2prot, Node * leaf,
			 double ** score, double ** complement_score, double ** p_value,
			 double **probability, double **overlap, 
			 double cutoff, int ** is_precedent);
int  output_tree (char * filename ,Node * node );
int  path_determinants (double sig_cutoff, Alignment * alignment,
			Tree * tree, Node ** path_to_query,
		       double  **per_node_prob, int **description);
int  print_tree (FILE * fptr,Node * node );
int  process_almt (Options *options, Alignment *alignment) ;
int  read_alignment (Options * options, Alignment * alignment);
int  read_cmd_file (char *filename, Options * options);
int  read_dssp (Options * options, Protein *protein);
int  read_epitope (Options * options, Protein * protein, int *epitope);
int  read_pdb      ( char * pdbname, Protein   * protein,  char chain);
int  read_score ( char *filename, Protein * protein, int * score2prot, double *score );
int  read_score_almt ( char *filename, Alignment * alignment, double *score );
int  read_tree (Options *options,  Alignment *alignment, Tree * tree);
int  scoring ( Options *options, Alignment * alignment, Tree *tree, 
	       int * similar_to, double * score, int score_ctr);
int  set_similarity ( Options * options, int * similar_to );
char single_letter ( char code[]);
int  struct_almt_mapping (Protein * protein, Alignment * alignment,
			  int * prot2almt, int * almt2prot);

int tokenize ( char  token[MAX_TOK][MEDSTRING], int * max_token, char * line, char comment_char);
# endif
