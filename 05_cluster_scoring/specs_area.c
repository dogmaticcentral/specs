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

