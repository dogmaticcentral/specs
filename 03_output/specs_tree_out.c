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

int output_tree (char * base_filename, Node *root) {
    FILE * fptr;
    char filename[BUFFLEN];
    int print_tree (FILE * fptr,Node * node );
    
    sprintf (filename, "%s.nhx", base_filename);
    fptr = efopen (filename, "w");
    if ( !fptr) return 1;
    
    print_tree (fptr, root);
    fprintf (fptr, "\n");
    
    fclose (fptr);

    return 0;

}
