# include "specs.h"

int output_clustering (char * base_filename,  int * int_cvg, double *clustering_score, int length){

    int ctr;
    FILE * fptr;
    char filename[BUFFLEN];
    double area_over_coverage (int * int_cvg, double * value, int no_of_res ) ;
    
    sprintf (filename, "%s.clustering", base_filename);

    fptr = efopen (filename, "w");
    if ( !fptr) return 1;
    
    fprintf ( fptr, "%%%8s%8s%12s\n", "cum res", "cvg ", "z  ");
  
    for (ctr=0; ctr < length &&int_cvg[ctr] ; ctr++ ) {
	fprintf ( fptr, "%8d%8.3lf%12.3le\n",
		  int_cvg[ctr], (double)int_cvg[ctr]/length, clustering_score[ctr]);
	
    }

    fprintf ( fptr, "\n%%  area:  %8.3lf\n", area_over_coverage(int_cvg, clustering_score, length) );
    fclose (fptr);
    
    return 0;
}
