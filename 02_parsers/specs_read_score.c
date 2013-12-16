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



/* this is geared towards r4s format*/
int read_score_almt ( char *filename, Alignment * alignment, double *score ) {
    
    FILE *fptr = NULL;
    char line[LONGSTRING];
    double score_val;
    int ctr, line_ctr, aux_ctr, *renumber;
    char aa, *refseq, *aux_seq;
    
    refseq = alignment->refseq[0];
    
    /* need to do some hacking bcs r4s skips the X's if they appear in the input */
    if ( ! (aux_seq  = (char *) emalloc ( alignment->length *sizeof(char))) ) return 1;
    if ( ! (renumber = (int *)  emalloc ( alignment->length *sizeof(int))) ) return 1;
     
    aux_ctr = 0;
    for(ctr=0; ctr < alignment->length; ctr++ ) {
	if ( refseq[ctr] != 'X' ) {
	    aux_seq[aux_ctr] = refseq[ctr];
	    renumber[aux_ctr] = ctr;
	    aux_ctr ++;
	}
    }
    fptr   = efopen ( filename, "r" );
    if (! fptr ) return 1;
    memset ( line, 0, LONGSTRING);
    ctr = 0;
    line_ctr = 0;
    while(fgets(line,LONGSTRING,fptr)!=NULL){
	line_ctr++;
	if ( ! isdigit (line[4]) ) continue;
	if ( line[0] == '#' ) continue;
	sscanf (line, " %d  %c  %lf \n", &ctr, &aa,  &score_val); /* notice the r4s format */
	if ( ctr >  alignment->length ) {
	    fprintf ( stderr,
		      "Error reading in the score: sequence in the score file too long, line %d.\n", line_ctr);
	    return 1;
	}
	ctr --;
	if ( aa != aux_seq[ctr] ) {
	    fprintf ( stderr,
		      "Error reading in the score: score file/referene sequence type mismatch at pos %d,"
		      " line %d (alignment: %c  score file: %c)\n",
		      ctr+1, line_ctr, aux_seq[ctr], aa );
	    return 1;
	}
	score[renumber[ctr]]      = score_val;
   }
    fclose (fptr);

    free (aux_seq);
    free ( renumber);

# if 0
    /* further hacking: in r4s some conserved residues are more conserved the others - force them to
       have the same score */
    double  * entr, min_score;
    int entropy ( Alignment * alignment, double *score);
    if ( ! (entr = (double *) emalloc ( alignment->length *sizeof(double))) ) return 1;
    entropy ( alignment, entr);
    min_score = 10000;
    for(ctr=0; ctr < alignment->length; ctr++ ) {
	if ( min_score > score[ctr] ) min_score = score[ctr];
    }
    for(ctr=0; ctr < alignment->length; ctr++ ) {
	if ( ! entr[ctr] ) score[ctr] = min_score;
    }

    free (entr);
# endif
    return 0;
}

/********************************************************************************/
int read_score ( char *filename, Protein * protein, int * score2prot, double *score ) {
    
    FILE *fptr = NULL;
    char line[LONGSTRING];
    char pdb_id[PDB_ATOM_RES_NO+2];
    double score_val;
    int res_ctr, ctr;
    
   
    fptr   = efopen ( filename, "r" );
    if (! fptr ) return 1;
    memset ( line, 0, LONGSTRING);
    ctr = 0;
    while(fgets(line,LONGSTRING,fptr)!=NULL){
	sscanf (line, " %s %lf \n", pdb_id,  &score_val);
	for ( res_ctr=0; res_ctr < protein->length; res_ctr++ ) {
	    if ( ! strncmp( protein->sequence[res_ctr].pdb_id, pdb_id, PDB_ATOM_RES_NO_LEN+2) ) {
		score2prot[ctr] = res_ctr;
		score[ctr]      = score_val;
		ctr ++;
	    }
	}
    }
    fclose (fptr);
    return 0;
}
