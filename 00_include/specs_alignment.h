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
# define ALMT_NAME_LENGTH 150

typedef struct{
    int number_of_seqs;
    int length;
    char ** sequence;
    char ** name;
    int * seq_gaps;
    int * column_gaps;
    Node **node;
    double **seq_dist;
    int  **aligned_sites;
    int  **identical_sites;
    int  **similar_sites;
    char * pdbseq;
    char ** refseq;
    char **refseq_name;
    
    int  no_refseqs;
    int * protected_position;
    int number_of_protected_positions;
    int refseq_gaps;

    int * number_of_exon_bdries; /* in each reference sequence */
    char ** exon_bdry;

}  Alignment;
