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
#include <errno.h>
#include <fcntl.h>

#include <sys/types.h>


void * emalloc(int	size)
{
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %u bytes", size);
	return NULL;
    }

    return ptr;
}



FILE * efopen(char * name, char * mode)
{

    FILE * fp;


    if ((fp = fopen(name, mode)) == NULL) {
	fprintf (stderr,  
	      "Cannot open \"%s\" for \"%s\"\n", name, mode);
	return NULL;
    }

    return fp;

}




/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] 
 */
char **chmatrix(int rows, int columns){
    char **m;
    int i;
        /* allocate pointers to rows */
    m=(char **) malloc(rows*sizeof(char*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(char *) calloc( rows*columns, sizeof(char));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}

int **intmatrix(int rows, int columns){
    int **m;
    int i;
        /* allocate pointers to rows */
    m=(int **) malloc(rows*sizeof(int*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(int *) calloc( rows*columns, sizeof(int));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}

double **dmatrix(int rows, int columns){
    double **m;
    int i;
        /* allocate pointers to rows */
    m=(double **) malloc(rows*sizeof(double*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(double *) calloc( rows*columns, sizeof(double));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}



/* free a  matrix  */
void free_matrix(void **m)
{
    free(m[0]);
    free(m);
}


/* free a  matrix  - save some warnings*/
void free_cmatrix(char **m)
{
    free(m[0]);
    free(m);
}
void free_imatrix(int **m)
{
    free(m[0]);
    free(m);
}
void free_dmatrix(double **m)
{
    free(m[0]);
    free(m);
}
void free_d3matrix(double ***m)
{
    free(m[0][0]);
    free(m[0]);
    free(m);
}
void free_strmatrix(char ***m)
{
    free(m[0][0]);
    free(m[0]);
    free(m);
}



/********************************************************************************************/


/* sort array according to the score in the other */
/* I couldn't declare pos_cmp within array_qsort  bcs it   crashed on mac */

double * score_array;

int pos_cmp (const void * a0, const void * b0) {
    
    int * a= (int*) a0;
    int * b= (int*) b0;
    if ( score_array[*a] > score_array[*b]) {
	return 1;
    }
    if ( score_array[*a] < score_array[*b]) {
	return -1;
    }
    return 0;
}

int array_qsort (int * sorted_pos, double * sa, int sequence_length ) {
    /* position comparison function */
    score_array = sa;

    qsort (sorted_pos, sequence_length, sizeof(int), pos_cmp);

    return 0;
}
/********************************************************************************************/
int tokenize ( char token[MAX_TOK][MEDSTRING], int * max_token, char * line , char comment_char) {
    /* assumes the tokens to be no bigger than MEDSTRING */ 
    
    char * chrptr, *last; 
    int current_token, current_char = 0;
    int reading;
   
    memset (token[0], 0, MAX_TOK*MEDSTRING*sizeof(char)); 
    chrptr = line;
    last   = chrptr + strlen (line);
    current_token = -1;
    current_char  =  0;
    reading = 0;
    while ( chrptr <= last) {
	if ( *chrptr == comment_char ) break;
	if ( *chrptr == '\n' ) break;
	if ( *chrptr && ! isspace(*chrptr) ) {
	    if ( ! reading ) {
		reading = 1;
		current_char = 0;
		current_token++;
		if ( current_token >= MAX_TOK ) {
		    return TOK_TOOMNY; /* defined in possum_utils.h */
		}
	    }
	    if ( current_char >= MEDSTRING ) {
		return TOK_TOOLONG;
	    }
	    token[current_token][current_char] = *chrptr;
	    current_char++;
	} else {
	    if ( reading ) {
		reading = 0;
	    }
	}
	chrptr++;
    }
    *max_token = current_token;

    return 0;
    
}

/**********************************************************/
/* get rid of spaces in a string */
int  string_clean ( char* string, int length) {
    int ctr;
    for (ctr = 0; ctr < length; ctr ++) {
	if ( isspace (string[ctr]) ) string[ctr] = '\0';
    }
    ctr=0;
    while ( !string[ctr]) ctr++;
    if ( ctr ) {
	memmove (string, string+ctr, length-ctr);
	memset ( string+length-1-ctr, 0, ctr);
    }

    return 0;
}
