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
