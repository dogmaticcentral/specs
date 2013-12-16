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


int set_similarity ( Options * options, int * similar_to ) {

    int i;

    for(i=0;i<ASCII;i++){
	similar_to[i] = i;
    }
# if 1
    similar_to['I'] = 'V';
    similar_to['L'] = 'V';
    similar_to['S'] = 'T';
    similar_to['D'] = 'E';
    similar_to['K'] = 'R';
    similar_to['Q'] = 'N';
    similar_to['.'] = '.';

    similar_to['A'] = 'V';
    similar_to['M'] = 'V';
    similar_to['G'] = 'V';
    similar_to['F'] = 'Y';
    similar_to['H'] = 'R';
# else

    similar_to['A'] = 'V';
    similar_to['I'] = 'V';
    similar_to['L'] = 'V';
    similar_to['G'] = 'V';
    similar_to['S'] = 'V';
    similar_to['T'] = 'V';
    similar_to['C'] = 'V';
    similar_to['M'] = 'V';
    
# endif   

    return 0;
}
