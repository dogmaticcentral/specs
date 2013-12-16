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

int determine_adj_matrix ( int ** adj_matrix, Residue * sequence, int no_res, double cutoff_dist){
    
    int resctr1, resctr2, atomctr1, atomctr2, next, noc;
    double dist, aux, cutoff_sq, x1, y1, z1;
    Atom * atomptr1, * atomptr2;

    if ( ! adj_matrix ) {
	fprintf (stderr,
		 "Error:  determine_adj_matrix() expects allocated space for the adj matrix.\n");
	return 1;
    } 

    cutoff_sq = cutoff_dist*cutoff_dist; 

    memset ( 	adj_matrix[0], 0, no_res*no_res*sizeof(int));
    /*diagonal*/
    for ( resctr1=0; resctr1<no_res; resctr1++) {
	adj_matrix [resctr1][resctr1] = 1;
    }
    noc = 0;
    for ( resctr1=0; resctr1<no_res; resctr1++) {
	for ( resctr2=resctr1+1; resctr2<no_res; resctr2++) {


	    next = 0;
	    for (atomctr1=0; atomctr1 <sequence[resctr1].no_atoms && !next; atomctr1++ ) {
		atomptr1 = sequence[resctr1].atom + atomctr1;
		x1 = atomptr1->x;
		y1 = atomptr1->y;
		z1 = atomptr1->z;
		
		for (atomctr2=0; atomctr2 <sequence[resctr2].no_atoms && !next; atomctr2++ ) {
		    
		    atomptr2 = sequence[resctr2].atom + atomctr2;
	
		    dist = 0.0;
		    aux = x1 -  atomptr2->x;
		    dist += aux*aux;
		    aux = y1 -  atomptr2->y;
		    dist += aux*aux;
		    aux = z1 -  atomptr2->z;
		    dist += aux*aux;
		    if ( dist < cutoff_sq  ) {
			adj_matrix [resctr1][resctr2] =  1;
			adj_matrix [resctr2][resctr1] =  1;
			next = 1;
			noc++;
		    }
		} /* end atomctr2 loop */
	    } /* end atomctr1 loop */
	    
	} 
    }

#if 0
    int contctr = 0;
    for ( resctr1=0; resctr1<no_res; resctr1++) {
	for ( resctr2= resctr1+1; resctr2<no_res; resctr2++) {
	    if ( adj_matrix [resctr1][resctr2] ) {
		//printf ("%4s %4s  %1d \n", sequence[resctr1].pdb_id,
		//sequence[resctr2].pdb_id, adj_matrix [resctr1][resctr2]);
		contctr ++;
	    }
	}
	//printf ("\n");
    }
	printf (" number of residues:  %4d               number  of contacts: %4d \n\n", no_res, contctr);
	exit(0);
# endif    
    return 0;
}

/*************************************************************************************************/
/*************************************************************************************************/
/*************************************************************************************************/
