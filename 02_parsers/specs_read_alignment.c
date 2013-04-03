# include "specs.h"


typedef enum {UNK, GCG, FASTA} AfaFormat;

int allocate_alignment_space (Alignment * alignment, Options * options, int number_of_seqs, int almt_length);
int get_name_from_fasta_hdr  (char * line, char * name);
int read_gcg   (Options * options, Alignment * alignment);
int read_fasta (Options * options, Alignment * alignment);


/*******************************************************/
int read_alignment ( Options * options, Alignment * alignment){

    FILE * fptr = NULL;
    AfaFormat format = UNK;
    
    /***************************************************/
    /*  some housekeeping                              */
    if ( options->no_refseqs) {
	alignment->refseq_name = chmatrix (options->no_refseqs*sizeof(char*), ALMT_NAME_LENGTH);
	if ( ! alignment->refseq_name) return 1;
	for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++) {
	    sprintf ( alignment->refseq_name[refseq_ctr], "%s", options->refseq_name[refseq_ctr]);
	}
	/* we'll asign these below */
	alignment->refseq = emalloc (options->no_refseqs*sizeof(char*));
	if ( ! alignment->refseq) return 1;

	refseq_found = emalloc (options->no_refseqs*sizeof(int));
	if ( !refseq_found) return 1;
    }
    
    /***************************************************/
    /* open file                                       */
    fptr = efopen ( options->almtname, "r");
    if ( !fptr ) return 1;

    /***************************************************/
    /* try to guess the input format                  */
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( line[0]=='>') {
	    format = FASTA;
	    break;
	    
	} else if ( strstr(line, "MSF:") ){
	    format = GCG;
	    break;
	} 
    }
    if ( format == UNK ) {
	fprintf ( stderr, "Alignment format not recognized. Is the format gcg or fasta?\n",
		  options->almtname);
	return 1;
    }

    /***************************************************/
    /*  read                                           */
    if (format == GCG) {
	read_gcg   (options, alignment);
    } else {
	read_fasta (options, alignment);
    }
    
    fclose(fptr);

    /***************************************************/
    /* sanity check                                    */
    for (ctr=0; ctr < number_of_seqs; ctr++ ) {
	if ( seq_pos[ctr] >  almt_length ) {
	    fprintf (stderr,
		     "Sequence %s is longer (%d positions) than the alignment (%d positions).\n",
		     name[ctr],  seq_pos[ctr], almt_length);
	    return 1;
	} else if ( seq_pos[ctr] <  almt_length ) {
	    fprintf (stderr,
		     "Sequence %s is shorter (%d positions) than the alignment (%d positions).\n",
		     name[ctr],  seq_pos[ctr], almt_length);
	    return 1;
	}
    }

    /***************************************************/
    /* more housekeeping                               */
    /*  reference and pdb sequence handling            */
    for (ctr=0; ctr < number_of_seqs; ctr++ ) {
	for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  {
	    if (! strcmp ( name[ctr], alignment->refseq_name[refseq_ctr]) ) {
		alignment->refseq[refseq_ctr] = sequence[ctr];
		break;
	    }
	}
    }

     if (  options->pdbname[0] && !options->skip_pdbseq ) {
	for (ctr=0; ctr < number_of_seqs; ctr++ ) {
	    if (! strcmp ( name[ctr], options->pdbseq_name) ) {
		alignment->pdbseq = sequence[ctr];
		break;
	    }
	}
    }
    
    /* return values */
    alignment->number_of_seqs = number_of_seqs;
    alignment->length         = almt_length;
    alignment->sequence       = sequence;
    alignment->name           = name;
    return 0;
}


/********************************************************************************************/
int read_gcg ( Options * options, Alignment * alignment){
    
    char line[BUFFLEN];
    int  number_of_seqs, almt_length, ctr, ctr2, refseq_ctr;
    int * seq_pos, pos, pdbseq_pos_ctr,pdbseq_found;
    int * refseq_found = NULL;
    int * pos_ctr;
    char * seq_ptr;
    char ** sequence;
    char ** name;
    char curr_name[BUFFLEN];
     
    memset (alignment, 0, sizeof(Alignment) );

    /* find the alignment length info */
    almt_length = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( strstr(line, "MSF:" ) ){
	    sscanf (line, "%*s %d", &almt_length);
	    break;
	}
    }
    if ( almt_length ) {
	/* printf ( "Alignment length in %s is %d.\n", cwname, almt_length); */
    } else {
	fprintf ( stderr, "Alignment length info not found in %s. Is the format gcg?\n",
		  options->almtname );
	return 1;
    }

    /* determine the number of sequences */
    number_of_seqs = 0;
    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++) refseq_found[refseq_ctr] = 0;
    pdbseq_found = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( ! strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) {
	    number_of_seqs++;
	    sscanf (line, "%*s %s", curr_name);
	    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  {
		if ( !strcmp(curr_name, alignment->refseq_name[refseq_ctr]) ) {
		    refseq_found[refseq_ctr] = 1;
		}
	    }
	    if ( options->pdbseq_name[0] && !strcmp(curr_name, options->pdbseq_name) )
		pdbseq_found = 1;
	}
    }
    if ( !number_of_seqs ) {
	fprintf ( stderr, "No sequences found in %s. Is the format gcg?\n",
		  options->almtname);
	return 1;
    }
    

    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  {
	if ( !  refseq_found[refseq_ctr] ) {
	    fprintf ( stderr, "Refseq  %s not found in %s.\n",
		      alignment->refseq_name[refseq_ctr],  options->almtname);
	    return 1;
	}
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

    /****************************/
    allocate_alignment_space (alignment, options, number_of_seqs, almt_length);


    /****************************/
    /* read in                  */
    rewind(fptr);
    ctr = 0;
    pdbseq_pos_ctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if (!  strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) {
	    sscanf (line, "%*s %s", curr_name);
	    if (options->pdbname[0]  &&
		options->skip_pdbseq &&  !strcmp (curr_name, options->pdbseq_name) )
		continue;
	    sprintf (name[ctr], "%s", curr_name);		
	    ctr ++;
	}
    }

    
    /* check for duplicate names */
    for (ctr = 0; ctr <number_of_seqs;  ctr++) {
	for (ctr2 = ctr+1; ctr2 <number_of_seqs;  ctr2++ ) {
	    if ( ! strcmp (name[ctr], name[ctr2]) ) {
		fprintf ( stderr, "Duplicate names  found in the header of %s:  %s (names %d and %d: %s and %s)\n",
			  options->almtname, name[ctr], ctr, ctr2, name[ctr], name[ctr2]);
		return 1;
	    }
	}
    }
    
    /* read in the sequences */
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( isspace (line[0] ) ) continue;
	sscanf (line, "%s", curr_name);
	if (  options->pdbname[0] && options->skip_pdbseq  && ! strcmp (options->pdbseq_name, curr_name)  ) {
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
	while ( ! isspace(line[pos]) ) pos++;
	while  (line[pos] != '\n' && pos < BUFFLEN) {
	    if ( !  isspace(line[pos] ) ){
                /* --> turn to uppercase */
		if ((line[pos]>=97)&&(line[pos]<=122)) {line[pos] -= 32;}
		/* turn dash to dot */
		if ( line[pos]==45 )                   {line[pos]  = 46;} 
		/* turn tweedle to dot */
		if ( line[pos]==126)                   {line[pos]  = 46;} 
		/* turn X to dot */
		if ( line[pos]==88)                    {line[pos]  = 46;} 
		seq_ptr [ *pos_ctr ] = line[pos];
		(*pos_ctr)++;
	    }
	    pos ++;
	}
    }
 	

    /* free */
    free (seq_pos);
    free (refseq_found);
    
    return 0;
}
/********************************************************************************************/
int read_fasta ( Options * options, Alignment * alignment){
    
    char line[BUFFLEN];
    int  number_of_seqs, almt_length, ctr, ctr2, refseq_ctr;
    int * seq_pos, pos, pdbseq_pos_ctr,pdbseq_found;
    int * refseq_found = NULL;
    int  pos_ctr;
    char * seq_ptr;
    char ** sequence;
    char ** name;
    char curr_name[BUFFLEN];
     
    memset (alignment, 0, sizeof(Alignment) );
    

    /* determine the number of sequences */
    number_of_seqs =  0;
    almt_length    =  0;
    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++) refseq_found[refseq_ctr] = 0;
    pdbseq_found = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){

	if ( line[0]== '>' ) {
	    
	    /*process the previous seq*/ 
	    if (!almt_length) {
		almt_length = seq_length;
	    } else if (seq_length != almt_length) {
		fprintf (stderr, "Sequence length mismatch for %s: ", curr_name);
		fprintf (stderr, "(alignment length: %d, this sequence: %d).\n", seq_length, almt_length);
		return 1;
	    }
	    number_of_seqs++;
	    if get_name_from_fasta_hdr (line, curr_name) return 1;
	    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  {
		if ( !strcmp(curr_name, alignment->refseq_name[refseq_ctr]) ) {
		    refseq_found[refseq_ctr] = 1;
		}
	    }
	    if ( options->pdbseq_name[0] && !strcmp(curr_name, options->pdbseq_name) )
		pdbseq_found = 1;
	    seq_length = 0
		
        } else {
	    for (pos=0; pos < BUFFLEN; pos++) {
		if  (!isspace(line[pos]) ) seq_length ++;
	    }
	    
	}
    }
    if ( !number_of_seqs ) {
	fprintf ( stderr, "No sequences found in %s. Is the format gcg?\n",
		  options->almtname);
	return 1;
    }
    

    for (refseq_ctr=0; refseq_ctr<options->no_refseqs; refseq_ctr++)  {
	if ( !  refseq_found[refseq_ctr] ) {
	    fprintf ( stderr, "Refseq  %s not found in %s.\n",
		      alignment->refseq_name[refseq_ctr],  options->almtname);
	    return 1;
	}
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
    
    /****************************/
    allocate_alignment_space (alignment, options, number_of_seqs, almt_length);

    
    /****************************/
    /* read in                  */
    rewind(fptr);

    while(fgets(line, BUFFLEN, fptr)!=NULL){

	if ( line[0]== '>' ) {
	    if get_name_from_fasta_hdr (line, curr_name) return 1;

	    if (  options->pdbname[0] && options->skip_pdbseq  && ! strcmp (options->pdbseq_name, curr_name)  ) {
		seq_ptr = alignment->pdbseq;
	    } else {
		ctr = 0;
		while (  ctr <number_of_seqs &&  strcmp (name[ctr], curr_name) ) ctr++;
		if ( ctr >= number_of_seqs ) {
		    fprintf ( stderr, "The name %s not found in the header of %s.\n",
			      curr_name,  options->almtname);
		    return 1;
		}
		seq_ptr = sequence [ctr];
	    }
	    
	} else {
	    pos_ctr = 0;
	    pos     = 0;
	    while ( ! isspace(line[pos]) ) pos++;
	    while  (line[pos] != '\n' && pos < BUFFLEN) {
		if ( !  isspace(line[pos] ) ){
		    /* --> turn to uppercase */
		    if ((line[pos]>=97)&&(line[pos]<=122)) {line[pos] -= 32;}
		    /* turn dash to dot */
		    if ( line[pos]==45 )                   {line[pos]  = 46;} 
		    /* turn tweedle to dot */
		    if ( line[pos]==126)                   {line[pos]  = 46;} 
		    /* turn X to dot */
		    if ( line[pos]==88)                    {line[pos]  = 46;} 
		    seq_ptr [ pos_ctr ] = line[pos];
		    pos_ctr++;
		}
		pos ++;
	    }

	}
    }

    
    
    /* check for duplicate names */
    for (ctr = 0; ctr <number_of_seqs;  ctr++) {
	for (ctr2 = ctr+1; ctr2 <number_of_seqs;  ctr2++ ) {
	    if ( ! strcmp (name[ctr], name[ctr2]) ) {
		fprintf ( stderr, "Duplicate names  found in the header of %s:  %s (names %d and %d: %s and %s)\n",
			  options->almtname, name[ctr], ctr, ctr2, name[ctr], name[ctr2]);
		return 1;
	    }
	}
    }

    /* free */
    free (seq_pos);
    free (refseq_found);
    
    return 0;
}



/************************************************************************************************************/
int allocate_alignment_space (Alignment * alignment, Options * options, int number_of_seqs, int almt_length) {
    
    printf ( "Number of sequences in %s is %d.\n",  options->almtname, number_of_seqs);
    
    /* allocate */
    sequence = chmatrix (number_of_seqs, almt_length);
    if ( !sequence ) return 1;


    if (  options->pdbname[0] && options->skip_pdbseq ) {
	alignment->pdbseq     =  emalloc (almt_length*sizeof(char));
	if ( !alignment->pdbseq) return 1;
    }
    
    alignment->seq_dist = dmatrix (number_of_seqs, number_of_seqs);
    if ( !alignment->seq_dist ) return 1;

    alignment->node = emalloc (almt_length*sizeof(Node*));
    if ( !alignment->node ) return 1;
    
    name     = chmatrix (number_of_seqs, ALMT_NAME_LENGTH);
    if ( !name ) return 1;

    seq_pos = (int *) emalloc ( number_of_seqs*sizeof(int));
    if ( !seq_pos ) return 1;


    return 0;
}


/************************************************************************************************************/
int get_name_from_fasta_hdr (char * line, char * name ) {
    
    char token[MAX_TOK][MEDSTRING] = {{'\0'}};
    int retval, max_token;
    cha comment_char;
    
    retval = tokenize ( token, &max_token, line, comment_char= '!' );
    switch ( retval ) {
    case  TOK_TOOMNY:
	errmsg ( log, line_ctr, line, "\t\t %s\n", "Too many tokens.");
	fclose (log);
	return 1;
    case TOK_TOOLONG:
	errmsg ( log, line_ctr, line, "\t\t %s\n", "Token too long.");
	fclose (log);
	return 1;
    }
    if ( max_token < 0 ) return 1;

    speintf(name, "%s", token[0]);

    
    return 0;
}
