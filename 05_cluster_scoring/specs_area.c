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

# if 0
double area_over_coverage (int * int_cvg, double * value, int no_of_res ) {
     
     double step_size = 1.0/no_of_res;
     double area;
     int cvg_ctr;
	
     area  = 0;
     for (cvg_ctr=0; cvg_ctr < no_of_res && int_cvg[cvg_ctr]; cvg_ctr ++) {
	 area += value[cvg_ctr];
     }
     
     area *= step_size;
     return area;
		
}
# endif

double area_over_coverage (int * int_cvg, double * value, int no_of_res)  {
     
     double step_size = 1.0/no_of_res;
     double area;
     int cvg_ctr, res_ctr;
     int empty_bins;
     double  prev,  smaller;
	
     area  = 0;
     cvg_ctr = 0;
     empty_bins = 0;
     prev = 0;
     
     for (res_ctr=1; res_ctr <= no_of_res; res_ctr ++) {
	 
	 if ( res_ctr <  int_cvg[cvg_ctr] ) {
	     empty_bins ++;
	     continue;
	 }
	 smaller = ( prev < value[cvg_ctr]) ? prev : value[cvg_ctr];
	 area += empty_bins*smaller + value[cvg_ctr];

	 empty_bins = 0;
	 cvg_ctr ++;
	 prev = value[cvg_ctr];
     }
     area *= step_size;
     return area;
		
 }

