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


typedef enum {UNK, GCG, FASTA} AfaFormat;

int allocate_alignment_space (Alignment * alignment, Options * options, int number_of_seqs, int almt_length);
int get_name_from_fasta_hdr  (char * line, char * name);
int housekeeping_and_sanity_checking (Options * options, Alignment * alignment);
int read_gcg   (Options * options, Alignment * alignment, FILE * fptr);
int read_fasta (Options * options, Alignment * alignment, FILE * fptr);


/*******************************************************/
int read_alignment ( Options * options, Alignment * alignment){

        FILE * fptr        = NULL;
        AfaFormat format   = UNK;
        char line[BUFFLEN] = {'\0'};
        int retval         = 0;


        /***************************************************/
        /* open file                                       */
        fptr = efopen ( options->almtname, "r");
        if ( !fptr ) return 1;

        /***************************************************/
        /* try to guess the input format                   */
        while(fgets(line, BUFFLEN, fptr)!=NULL) {
                if ( line[0]=='>') {
                        format = FASTA;
                        break;

                } else if ( strstr(line, "MSF:") ) {
                        format = GCG;
                        break;
                }
        }
        if ( format == UNK ) {
                fprintf ( stderr, "Alignment format in %s not recognized. Is the format gcg or fasta?\n",
                          options->almtname);
                return 1;
        }

        /***************************************************/
        /*  read                                           */
        rewind(fptr);
        if (format == GCG) {
                retval = read_gcg   (options, alignment, fptr);
        } else {
                retval = read_fasta (options, alignment, fptr);
        }
        fclose(fptr);
        if (retval) return retval;


        retval = housekeeping_and_sanity_checking (options, alignment);
        if (retval) return retval;

        return 0;
}


/********************************************************************************************/
int read_gcg ( Options * options, Alignment * alignment, FILE * fptr){

        char line[BUFFLEN];
        int number_of_seqs, almt_length, ctr;
        int * seq_pos, pos, pdbseq_pos_ctr,pdbseq_found;
        int * pos_ctr;
        char * seq_ptr;
        char ** sequence;
        char ** name;
        char curr_name[BUFFLEN];

        memset (alignment, 0, sizeof(Alignment) );

        /* find the alignment length info */
        almt_length = 0;
        while(fgets(line, BUFFLEN, fptr)!=NULL) {
                if ( strstr(line, "MSF:" ) ) {
                        sscanf (line, "%*s %d", &almt_length);
                        break;
                }
        }
        if ( almt_length ) {
                //printf ( "Alignment length in %s is %d.\n", options->almtname, almt_length);
        } else {
                fprintf ( stderr, "Alignment length info not found in %s. Is the format gcg?\n",
                          options->almtname );
                return 1;
        }

        /* determine the number of sequences */
        number_of_seqs = 0;
        pdbseq_found   = 0;
        while(fgets(line, BUFFLEN, fptr)!=NULL) {
                if ( !strncmp (line, "//", 2) ) break;
                if ( strstr(line, "Name:" ) ) {
                        number_of_seqs++;
                        sscanf (line, "%*s %s", curr_name);
                        if ( options->pdbseq_name[0] && !strcmp(curr_name, options->pdbseq_name) )
                                pdbseq_found = 1;
                }
        }
        if ( !number_of_seqs ) {
                fprintf ( stderr, "No sequences found in %s. Is the format gcg?\n",
                          options->almtname);
                return 1;
        }

        if ( options->pdbname[0] ) {
                if ( pdbseq_found ) {
                        if (options->skip_pdbseq) number_of_seqs--; /* we will store the pdbseq separately */
                } else {
                        fprintf ( stderr, "Sequence corresponding to the structure (%s) not found in %s.\n",
                                  options->pdbseq_name,  options->almtname);
                        return 1;

                }
        }

        /****************************************************/
        allocate_alignment_space (alignment, options, number_of_seqs, almt_length);
        seq_pos = (int *) emalloc ( number_of_seqs*sizeof(int));
        if ( !seq_pos ) return 1;
        sequence = alignment->sequence; /* alias */
        name     = alignment->name;


        /***************************************************/
        /* read in the names                               */
        rewind(fptr);
        ctr = 0;
        pdbseq_pos_ctr = 0;
        while(fgets(line, BUFFLEN, fptr)!=NULL) {
                if (!strncmp (line, "//", 2) ) break;
                if ( strstr(line, "Name:" ) ) {
                        sscanf (line, "%*s %s", curr_name);
                        if (options->pdbname[0]  &&
                            options->skip_pdbseq &&  !strcmp (curr_name, options->pdbseq_name) )
                                continue;
                        sprintf (name[ctr], "%s", curr_name);
                        ctr++;
                }
        }

        /***************************************************/
        /* read in the sequences */
        while(fgets(line, BUFFLEN, fptr)!=NULL) {
                if ( isspace (line[0] ) ) continue;
                sscanf (line, "%s", curr_name);
                if (  options->pdbname[0] && options->skip_pdbseq  && !strcmp (options->pdbseq_name, curr_name)  ) {
                        seq_ptr = alignment->pdbseq;
                        pos_ctr = &pdbseq_pos_ctr;
                } else {
                        ctr = 0;
                        while (  ctr <number_of_seqs &&  strcmp (name[ctr], curr_name) ) ctr++;
                        if ( ctr >= number_of_seqs ) {
                                fprintf ( stderr, "The name %s not found in the header of %s.\n",
                                          curr_name,  options->almtname);
                                return 1;
                        }
                        seq_ptr = sequence [ctr];
                        pos_ctr = seq_pos + ctr;
                }
                pos = 0;
                while ( !isspace(line[pos]) ) pos++; /* skip the name */
                while  (line[pos] != '\n' && pos < BUFFLEN) {
                        if ( !isspace(line[pos] ) ) {
                                if ( *pos_ctr == almt_length) {
                                        fprintf (stderr,
                                                 "Sequence %s is longer  than the alignment (%d positions).\n",
                                                 curr_name, almt_length);
                                        return 1;
                                }
                                seq_ptr [ *pos_ctr ] = line[pos];
                                (*pos_ctr)++;
                        }
                        pos++;
                }
        }

        /* return values */
        alignment->number_of_seqs = number_of_seqs;
        alignment->length         = almt_length;
        alignment->refseq_name    = options->refseq_name;
        alignment->no_refseqs     = options->no_refseqs;

        /* free */
        free (seq_pos);

        return 0;
}


/********************************************************************************************/
int read_fasta ( Options * options, Alignment * alignment, FILE * fptr){

        char line[BUFFLEN];
        int number_of_seqs, almt_length, ctr;
        int pos, pdbseq_found;
        int seq_length;
        int line_pos, seq_pos;
        char * seq_ptr;
        char ** sequence;
        char ** name;
        char curr_name[BUFFLEN];

        memset (alignment, 0, sizeof(Alignment) );

        /* determine the number of sequences */
        number_of_seqs = 0;
        almt_length    = 0;
        pdbseq_found   = 0;
        seq_length     = 0;
        while(fgets(line, BUFFLEN, fptr)!=NULL) {

                if ( line[0]== '>' ) {
                        /*process the previous seq*/
                        if (seq_length && !almt_length) {
                                almt_length = seq_length;
                        } else if (seq_length != almt_length) {
                                fprintf (stderr, "Sequence length mismatch for %s: ", curr_name);
                                fprintf (stderr, "(first sequence length: %d, this sequence: %d).\n", seq_length, almt_length);
                                return 1;
                        }
                        number_of_seqs++;
                        if ( get_name_from_fasta_hdr (line, curr_name) ) return 1;

                        if ( options->pdbseq_name[0] && !strcmp(curr_name, options->pdbseq_name) )
                                pdbseq_found = 1;

                        seq_length = 0;

                } else {
                        for (pos=0; pos < BUFFLEN; pos++) {
                                if (line[pos] =='\n') break;
                                if (!isspace(line[pos]) ) seq_length++;
                        }
                }
        }
        if ( !number_of_seqs ) {
                fprintf ( stderr, "No sequences found in %s. Is the format fasta?\n",
                          options->almtname);
                return 1;
        }


        if ( options->pdbname[0] ) {
                if ( pdbseq_found ) {
                        if (options->skip_pdbseq) number_of_seqs--; /* we will store the pdbseq separately */
                } else {
                        fprintf ( stderr, "Sequence corresponding to the structure (%s) not found in %s.\n",
                                  options->pdbseq_name,  options->almtname);
                        return 1;
                }
        }

        /******************************************************************************/
        allocate_alignment_space (alignment, options, number_of_seqs, almt_length);
        sequence = alignment->sequence; /* alias */
        name     = alignment->name;

        /****************************/
        /* read in                  */
        rewind(fptr);
        ctr     = 0;
        seq_ptr = NULL;
        seq_pos = 0;
        while(fgets(line, BUFFLEN, fptr)!=NULL) {

                if ( line[0]== '>' ) {

                        if ( get_name_from_fasta_hdr (line, curr_name) ) return 1;

                        if (  options->pdbname[0] && options->skip_pdbseq  && !strcmp (options->pdbseq_name, curr_name)  ) {
                                seq_ptr = alignment->pdbseq;
                        } else {
                                sprintf (name[ctr], "%s",  curr_name);
                                seq_ptr = sequence [ctr];
                                ctr++;
                        }
                        seq_pos = 0;

                } else {
                        for (line_pos = 0; line[line_pos] != '\n' && line_pos < BUFFLEN; line_pos++) {
                                if (  isspace(line[line_pos] ) ) continue;
                                seq_ptr [ seq_pos ] = line[line_pos];
                                seq_pos++;
                        }
                }
        }

        /* return values */
        alignment->number_of_seqs = number_of_seqs;
        alignment->length         = almt_length;
        alignment->refseq_name    = options->refseq_name;
        alignment->no_refseqs     = options->no_refseqs;

        return 0;
}



/************************************************************************************************************/
int allocate_alignment_space (Alignment * alignment, Options * options, int number_of_seqs, int almt_length) {

        //printf ( "Number of sequences in %s is %d.\n",  options->almtname, number_of_seqs);

        /* allocate */
        alignment->sequence = chmatrix (number_of_seqs, almt_length);
        if ( !alignment->sequence ) return 1;

        if (  options->pdbname[0] && options->skip_pdbseq ) {
                alignment->pdbseq     =  emalloc (almt_length*sizeof(char));
                if ( !alignment->pdbseq) return 1;
        }

        alignment->name     = chmatrix (number_of_seqs, ALMT_NAME_LENGTH);
        if ( !alignment->name ) return 1;

        alignment->refseq   = emalloc (options->no_refseqs*sizeof(char*));
        if ( !alignment->refseq ) return 1;

        /* this seems a bit out of place here */
        alignment->node     = emalloc (almt_length*sizeof(Node*));
        if ( !alignment->node ) return 1;

        alignment->no_refseqs = options->no_refseqs;

        return 0;
}


/************************************************************************************************************/
int get_name_from_fasta_hdr (char * line, char * name ) {

        char token[MAX_TOK][MEDSTRING] = {{'\0'}};
        int retval, max_token;
        char comment_char;


        /* the first chracter is '>' */
        retval = tokenize ( token, &max_token, line+1, comment_char= '!' );
        switch ( retval ) {
        case  TOK_TOOMNY:
                fprintf ( stderr, "%s\n%s\n", line, "Too many tokens.");
                return 1;
        case TOK_TOOLONG:
                fprintf ( stderr, "%s\n%s\n", line, "Token too long.");
                return 1;
        }
        if ( max_token < 0 ) return 1;

        sprintf(name, "%s", token[0]);


        return 0;
}

/*******************************************************************************************/
int seq_cleanup (char * seq, int length){

        if ( !seq) return 0; /* not my job */
        int pos;
        for (pos=0; pos<length; pos++) {
                /* --> turn to uppercase */
                if ((seq[pos]>=97)&&(seq[pos]<=122)) {seq[pos] -= 32;}
                /* turn dash  to dot */
                if ( seq[pos]==45 )                  {seq[pos]  = 46;}
                /* turn tweedle to dot */
                if ( seq[pos]==126)                  {seq[pos]  = 46;}
                /* turn X to dot */
                if ( seq[pos]==88)                   {seq[pos]  = 46;}

        }
        return 0;
}


/*******************************************************************************************/
int housekeeping_and_sanity_checking (Options * options, Alignment * alignment) {

        int ctr, refseq_ctr;
        int refseq_found = 0;

        /* referece sequences - names and pointers to the sequence          */
        if ( options->no_refseqs) {

                alignment->refseq_name = chmatrix (options->no_refseqs*sizeof(char*), ALMT_NAME_LENGTH);
                if ( !alignment->refseq_name) return 1;
                for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++) {
                        sprintf ( alignment->refseq_name[refseq_ctr], "%s", options->refseq_name[refseq_ctr]);
                }
                /* we'll asign these below */
                alignment->refseq = emalloc (options->no_refseqs*sizeof(char*));
                if ( !alignment->refseq) return 1;

                /* do we have all refseqs that were anounced in the cmd file? use refseq_found to keep track*/
                for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  {
                        refseq_found = 0;
                        for (ctr=0; ctr < alignment->number_of_seqs; ctr++ ) {
                                if ( !strcmp(alignment->name[ctr], options->refseq_name[refseq_ctr]) ) {
                                        refseq_found = 1;
                                        break;
                                }
                        }
                        if ( !refseq_found) {
                                fprintf (stderr, "Reference sequence %s not found in %s\n",
                                         options->refseq_name[refseq_ctr], options->almtname);
                                return 1;
                        }
                }
                /* pointers to reference sequences        */
                for (ctr=0; ctr < alignment->number_of_seqs; ctr++ ) {
                        for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  {
                                if (!strcmp ( alignment->name[ctr], alignment->refseq_name[refseq_ctr]) ) {
                                        alignment->refseq[refseq_ctr] = alignment->sequence[ctr];
                                        break;
                                }
                        }
                }
        }

        /* check for duplicate names                       */
        for (ctr = 0; ctr <alignment->number_of_seqs; ctr++) {
                int ctr2;
                for (ctr2 = ctr+1; ctr2 <alignment->number_of_seqs; ctr2++ ) {
                        if ( !strcmp (alignment->name[ctr], alignment->name[ctr2]) ) {
                                fprintf ( stderr, "Duplicate sequence names  found in %s:  names %d and %d: %s and %s\n",
                                          options->almtname, ctr, ctr2, alignment->name[ctr], alignment->name[ctr2]);
                                return 1;
                        }
                }
        }

        /* length check                                    */
        for (ctr=0; ctr < alignment->number_of_seqs; ctr++ ) {
                if (strlen(alignment->sequence[ctr]) < alignment->length ) {
                        fprintf (stderr,
                                 "Sequence %s is shorter (%d positions) than the alignment (%d positions).\n",
                                 alignment->name[ctr], (int)strlen(alignment->sequence[ctr]), alignment->length  );
                        return 1;
                }
        }


        /* pointer to the sequence that corresponds to the structure */
        if (  options->pdbname[0] && !options->skip_pdbseq ) {
                for (ctr=0; ctr < alignment->number_of_seqs; ctr++ ) {
                        if (!strcmp ( alignment->name[ctr], options->pdbseq_name) ) {
                                alignment->pdbseq = alignment->sequence[ctr];
                                break;
                        }
                }
        }

        /* cleanup the sequences:  */
        for (ctr=0; ctr < alignment->number_of_seqs; ctr++ ) {
                seq_cleanup (alignment->sequence[ctr],alignment->length);
        }

        if (options->pdbname[0]) seq_cleanup (alignment->pdbseq, alignment->length);

        return 0;
}
