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


int read_cmd_file (char *filename, Options * options) {
        FILE * fptr, *log = stdout;
        char line[LONGSTRING];
        char token[MAX_TOK][MEDSTRING] = {{'\0'}};
        char comment_char;
        int max_token, current_char = 0;
        int line_ctr, token_ctr, retval;
        int method_ctr;
        int would_use_extern = 0;
        int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
                     char * fmt, char * warnstr);
        int find_method ( char *token, int * method, char * method_name);


        fptr   = efopen ( filename, "r" );
        if (!fptr ) return 1;

        memset (options, 0, sizeof(Options));
        options->do_scoring = 1;
        options->tree_method = UPGMA;
        options->patch_sim_cutoff = -1;
        options->patch_min_length  = 0.9;
        options->min_fragment_length = 0.66;

        line_ctr = 0;
        memset ( line, 0, LONGSTRING);
        while(fgets(line,LONGSTRING,fptr)!=NULL) {
                line_ctr++;
                /* tokenize */
                retval = tokenize ( token, &max_token, line, comment_char= '!' );
                switch ( retval ) {
                case  TOK_TOOMNY:
                        errmsg ( log, line_ctr, line, "\t\t %s\n", "Too many tokens.");
                        fclose (log);
                        break;
                case TOK_TOOLONG:
                        errmsg ( log, line_ctr, line, "\t\t %s\n", "Token too long.");
                        fclose (log);
                        break;
                }
                if ( max_token < 0 ) continue;


                /* turn  first token to lowercase */
                current_char = 0;
                while ( token[0][current_char] &&  (current_char < SHORTSTRING) ) {
                        if ( token[0][current_char] >= 65 && token[0][current_char] <= 90)
                                token[0][current_char] += 32;
                        current_char++;
                }

                /* check first token for meaning (first token should be a keyword)*/
                if (  !strncmp (token[0], "align", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by a file name (full path).\n",
                                         token[0]);
                                return 1;
                        }
                        sprintf ( options->almtname, "%s", token[1]);
                } else if  (  !strncmp (token[0], "acc", 3)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by the accessibility cutoff.\n",
                                         token[0]);
                                return 1;
                        }
                        options->acc_cutoff = atoi(token[1]);

                } else if  (  !strncmp (token[0], "chai", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by a chain identifier.\n",
                                         token[0]);
                                return 1;
                        }
                        options->chain = token[1][0];

                } else if  (  !strncmp (token[0], "dssp", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by a file name (full path).\n",
                                         token[0]);
                                return 1;
                        }
                        sprintf ( options->dssp_file_name, "%s", token[1]);

                } else if  (  !strncmp (token[0], "epitope", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by a file name (full path).\n",
                                         token[0]);
                                return 1;
                        }
                        sprintf ( options->epi_file_name, "%s", token[1]);

                } else if  (  !strncmp (token[0], "r4s", 3)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by the r4s output file name (full path).\n",
                                         token[0]);
                                return 1;
                        }
                        sprintf ( options->scorename, "%s", token[1]);

                } else if  (  !strncmp (token[0], "intree", 6)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by a file name (full path).\n",
                                         token[0]);
                                return 1;
                        }
                        sprintf ( options->tree_file_name, "%s", token[1]);

                } else if (  !strncmp (token[0], "meth", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by a (at least one) name.\n",
                                         token[0]);
                                return 1;
                        }
                        options->number_of_methods = max_token;
                        options->method = emalloc (options->number_of_methods*sizeof(int));
                        if ( !options->method) return 1;
                        options->method_name = chmatrix (options->number_of_methods, SHORTSTRING);
                        if ( !options->method_name) return 1;
                        method_ctr = 0;
                        for (token_ctr = 1; token_ctr<= max_token; token_ctr++) {
                                /* turn  second token to lowercase */
                                current_char = 0;
                                while ( token[token_ctr][current_char] &&  (current_char < SHORTSTRING) ) {
                                        if ( token[token_ctr][current_char] < 90)
                                                token[token_ctr][current_char] += 32;
                                        current_char++;
                                }
                                retval = find_method ( token[token_ctr], options->method+method_ctr,
                                                       options->method_name[method_ctr]);
                                if (retval) { /* unrecognized method */
                                        options->number_of_methods--;
                                } else {
                                        if ( options->method[method_ctr] == EXTERN ) {
                                                /* this way the extern word may or may not be used,
                                                   as long as the external file is provided */
                                                would_use_extern = 1;
                                        }
                                        method_ctr++;
                                }
                        }

                } else if (  !strncmp (token[0], "min_frag_len", 12)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by min allowed fragment length.\n",
                                         token[0]);
                                return 1;
                        }
                        options->min_fragment_length = atof (token[1]);

                } else if (  !strncmp (token[0], "node", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by the node id.\n", token[0]);
                                return 1;
                        } else if ( !atoi (token[1]) ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Node id should be larger than 0.\n", "");
                        }
                        options->queried_node = atoi (token[1]);

                } else if (  !strncmp (token[0], "noscore", 7)  ) {
                        options->do_scoring = 0;

                } else if (  !strncmp (token[0], "outn", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by a name.\n", token[0]);
                                return 1;
                        }
                        sprintf ( options->outname, "%s", token[1]);

                } else if (  !strncmp (token[0], "patch_sim_cutoff", 16)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by min fraction seq. similarity to be used in patching.\n",
                                         token[0]);
                                return 1;
                        }
                        options->patch_sim_cutoff = atof (token[1]);
                } else if (  !strncmp (token[0], "patch_min_length", 16)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by min common seq. length fraction  to be used in patching.\n",
                                         token[0]);
                                return 1;
                        }
                        options->patch_min_length = atof (token[1]);


                } else if (  !strncmp (token[0], "path", 4)  ) {
                        options->path = 1;
                } else if  (  !strncmp (token[0], "pdbf", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by a file name (full path).\n",
                                         token[0]);
                                return 1;
                        }
                        sprintf ( options->pdbname, "%s", token[1]);

                } else if (  !strncmp (token[0], "pdbseq", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by a name.\n", token[0]);
                                return 1;
                        }
                        sprintf ( options->pdbseq_name, "%s", token[1]);

                } else if (  !strncmp (token[0], "raw", 3)  ) {
                        options->raw_score = 1;

                } else if (  !strncmp (token[0], "refseq", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by (at least one) name.\n", token[0]);
                                return 1;
                        }
                        options->no_refseqs = max_token;
                        options->refseq_name = chmatrix (options->no_refseqs, ALMT_NAME_LENGTH);
                        if ( !options->refseq_name) return 1;
                        for (token_ctr = 1; token_ctr<= max_token; token_ctr++) {
                                sprintf ( options->refseq_name[token_ctr-1], "%s", token[token_ctr]);
                        }

                } else if (  !strncmp (token[0], "restr", 5)  ) {
                        options->restrict2structure = 1;

                } else if (  !strncmp (token[0], "sign", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed  by a significance cutoff\n",
                                         token[0]);
                                return 1;
                        }
                        options->significance_cutoff = atof (token[1]);
                } else if (  !strncmp (token[0], "sink", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed  by a percentage of gaps which \"sinks\" a position.\n",
                                         token[0]);
                                return 1;
                        }
                        options->max_gaps = atof (token[1]);
                } else if (  !strncmp (token[0], "skip", 4)  ) {
                        options->skip_pdbseq = 1;
                } else if (  !strncmp (token[0], "tree", 4)  ) {
                        if ( max_token < 1 ) {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Keyord %s should be followed by a name.\n", token[0]);
                                return 1;
                        }
                        if ( !strncmp (  token[1], "nj", 2) ) {
                                options->tree_method = NEIGHBOR_JOINING;
                        } else if ( !strncmp (  token[1], "consensus", 5) )  {
                                options->tree_method = CONSENSUS_UPGMA;
                        } else if ( !strncmp (  token[1], "upgma", 5) )  {
                                options->tree_method = UPGMA;
                        } else {
                                errmsg ( log, line_ctr, line,
                                         "\t\t Unrecognozed tree building method: %s.\n",
                                         token[0]);
                                return 1;
                        }
                } else {
                        errmsg ( log, line_ctr, line, "\t\t Keyword %s not recognized.\n", token[0]);
                        return 1;

                }
                memset (line, 0, LONGSTRING);
        }
        fclose (fptr);

        if ( !options->pdbseq_name[0] && options->pdbname[0] ) {
                fprintf ( stderr, "Pdb sequence should be part of the alignment, ");
                fprintf ( stderr, "and the name of the pdb sequence provided to do the mapping.\n");
                return 1;
        }
        if ( !options->scorename[0] && !options->almtname[0] ) {
                fprintf ( stderr, "Neither score nor the alignment are provided.");
                fprintf ( stderr, "(What am I supposed to do?).\n");
                return 1;
        }
        if ( !options->outname[0]) sprintf (options->outname, "%s", "specs");
        if ( !options->acc_cutoff ) options->acc_cutoff = 10;
        if ( !options->significance_cutoff ) options->significance_cutoff = SIGNIFICANCE_CUTOFF;
        /* rvet set as default */
        if ( !options->number_of_methods) {
                options->number_of_methods = 1;
                options->method = emalloc (options->number_of_methods*sizeof(int));
                if ( !options->method) return 1;
                options->method_name = chmatrix (options->number_of_methods, SHORTSTRING);
                if ( !options->method_name) return 1;
                options->method[0] = RVET;
                sprintf (options->method_name[0], "%s", "rvet");
        }

        if (  would_use_extern && !options->scorename[0]) {
                fprintf ( stderr, "To use external scoring method, \"insc\" ");
                fprintf ( stderr, "keyword must be given in %s file,\n", filename);
                fprintf ( stderr, "followed by the external score file name.\n");
                return 1;
        }

        int number_of_methods = options->number_of_methods;


        if ( !would_use_extern && options->scorename[0] ) {
                number_of_methods++;
        }
        if (  options->queried_node ) {
                number_of_methods++;
        }

        if ( number_of_methods !=  options->number_of_methods ) {
                int *aux_int     = emalloc (number_of_methods*sizeof(int));
                if ( !aux_int ) return 1;
                char ** aux_name = chmatrix (number_of_methods, SHORTSTRING);
                if ( !aux_name ) return 1;

                memcpy (aux_int, options->method,
                        options->number_of_methods*sizeof(int));
                memcpy (aux_name[0], options->method_name[0],
                        options->number_of_methods*SHORTSTRING*sizeof(char) );

                number_of_methods = options->number_of_methods;
                if (  !would_use_extern && options->scorename[0] ) {
                        aux_int[number_of_methods] = EXTERN;
                        sprintf (aux_name[number_of_methods], "%s", "r4s");
                        number_of_methods++;
                }
                if (  options->queried_node ) {
                        aux_int[number_of_methods] = DET;
                        sprintf (aux_name[number_of_methods], "%s", "det");
                        number_of_methods++;
                }

                free (options->method);
                free (options->method_name[0]);
                free (options->method_name);


                options->method = aux_int;
                options->method_name = aux_name;
                options->number_of_methods = number_of_methods;
        }

        return 0;

}

/****************************************************************/
int find_method ( char *token, int * method, char * method_name) {


        if (!strcmp ( token, "entr" )) {
                *method = ENTROPY;
                sprintf (method_name, "%s", "entr");

        } else if (!strcmp ( token, "valdar" )) {
                *method = VALDAR;
                sprintf (method_name, "%s", "valdar");

        } else if (!strcmp ( token, "rvet" )) {
                *method = RVET;
                sprintf (method_name, "%s", "rvet");

        } else if (!strcmp ( token, "rvn" )) {
                *method = RVET_NORM;
                sprintf (method_name, "%s", "rvn");

        } else if (!strcmp ( token, "maj" )) {
                *method = MAJORITY;
                sprintf (method_name, "%s", "majf");

        } else if (!strcmp ( token, "pheno" )) {
                *method = PHENO;
                sprintf (method_name, "%s", "pheno");

        } else if (!strcmp ( token, "ivet" )) {
                *method = IVET;
                sprintf (method_name, "%s", "ivet");

        } else if (!strcmp ( token, "carb" )) {
                *method = CARBONE;
                sprintf (method_name, "%s", "carb");

        } else if (!strcmp ( token, "entr_s" )) {
                *method = ENTR_W_SIM;
                sprintf (method_name, "%s", "entr_s");

        } else if (!strcmp ( token, "rv_s" )) {
                *method = RV_W_SIM;
                sprintf (method_name, "%s", "rvs");

        } else if (!strcmp ( token, "maj_s" )) {
                *method = MAJORITY_W_SIM;
                sprintf (method_name, "%s", "majf_s");

        } else if (!strcmp ( token, "ivet_s" )) {
                *method = IV_W_SIM;
                sprintf (method_name, "%s", "ivet_s");

        } else if (!strcmp ( token, "phen_h" )) {
                *method = PHENO_HYBRID;
                sprintf (method_name, "%s", "phen_h");

        } else if (!strcmp ( token, "phen_n" )) {
                *method = PHENO_NODES;
                sprintf (method_name, "%s", "phen_n");
        } else if  (  !strncmp (token, "r4s", 4)  ) {
                *method = EXTERN;

        } else {
                /* method not recognized */
                fprintf ( stderr, "Method \"%s\" not recognized.\n", token);
                return 1;
        }

        return 0;

}
/****************************************************************************/
int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
             char * fmt, char * warnstr) {

        fprintf ( fptr, "\tError on line %3d:     %s", line_ctr, line);
        fprintf ( fptr, fmt, warnstr);
        return 0;
}
