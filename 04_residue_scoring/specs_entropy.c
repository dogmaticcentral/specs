# include "specs.h"

int entropy ( Alignment * alignment, int * similar_to,  int normalize, double *score){

    int col, seq, ctr, aa;
    int freq[ASCII]; /* ASCII == 128,  ASCII size  */
    double p, entropy;
    
    for (col = 0; col < alignment->length; col++){
	/* find frequencies */
	int col_length = 0;
	memset (freq, 0, ASCII*sizeof(int));
	for (seq=0; seq< alignment->number_of_seqs; seq++){
	    aa = (int) alignment->sequence[seq][col];
	    if ( similar_to) aa= similar_to[aa];
	    freq [aa] ++;
	    if (normalize && aa != '.') col_length ++;
	}
	/* find entropy */
	entropy = 0.0;
	for ( ctr=0; ctr < 128; ctr++) {
	    if ( freq[ctr] ) {
		if (normalize){
		    p = (double)freq[ctr]/col_length;
		} else {
		    p = (double)freq[ctr]/alignment->number_of_seqs;
		}
		entropy -= p*log(p);
	    }
	}
	score[col] = entropy;
    }

    return 0;
}
