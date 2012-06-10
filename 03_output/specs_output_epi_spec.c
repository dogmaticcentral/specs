# include "specs.h"

int  output_epi_spec ( Options * options, Protein * protein, int ** int_cvg,
		       double **sensitivity, double **specificity,int surface){
    int L = protein->length;
    int ctr, cvg_ctr, score_ctr;
    double current_sens, current_spec, cvg;
    double ** new_sens, ** new_spec;
    FILE * fptr;
    char filename[BUFFLEN];
    int number_of_scores = options->number_of_methods;
    
    if (! (new_sens = dmatrix(number_of_scores, L)))  exit (1);
    if (! (new_spec = dmatrix(number_of_scores, L)))  exit (1);
    
    if ( surface ) {
	int resctr;
	L = 0;
	for ( resctr=0; resctr< protein->length; resctr++) {
	    L += protein->sequence[resctr].solvent_accessible;
	}
    } else {
	L = protein->length;
    }

    /* to compare different scoring methods,
       to each coverage for which there is no score, ascribe the previous value */
    for ( score_ctr=0; score_ctr < number_of_scores; score_ctr++) {
	cvg_ctr = 0;
	current_sens = 0.0;
	current_spec = 1.0;
	for ( ctr=0; ctr<L && int_cvg[score_ctr][cvg_ctr]; ctr++) {
	    if ( ctr+1 >= int_cvg[score_ctr][cvg_ctr] ) {
		current_sens = sensitivity[score_ctr][cvg_ctr];
		current_spec = specificity[score_ctr][cvg_ctr];
		cvg_ctr ++;
	    }
	    new_sens[score_ctr][ctr] = current_sens;
	    new_spec[score_ctr][ctr] = current_spec;
	}
    }
    
    sprintf (filename, "%s.overlap", options->outname);

    fptr = efopen (filename, "w");
    if ( !fptr) return 1;
    fprintf (fptr, "%%%8s", "cvg");
    for ( score_ctr=0; score_ctr<number_of_scores; score_ctr++) {
	fprintf (fptr, " %16s", options->method_name[score_ctr]);
    }
    fprintf (fptr, "\n");
    for ( ctr=0; ctr<L; ctr++) {
	cvg = (double)(ctr+1)/L;
	fprintf (fptr, " %8.3lf", cvg);
	for ( score_ctr=0; score_ctr < number_of_scores; score_ctr++) {
	    fprintf (fptr, " %8.3lf%8.3lf", new_spec[score_ctr][ctr], new_sens[score_ctr][ctr]);
	}
	fprintf ( fptr, "\n");
   }
    fclose(fptr);
    
    free_matrix ((void*)new_sens);
    free_matrix ((void*)new_spec);

    return 0;
}
