# include "specs.h"

int read_epitope (Options * options, Protein * protein, int *epitope){

    char buf [LONGSTRING];
    char aux [PDB_ATOM_RES_NO_LEN+2];
    char res_type;
    int  found, total, res_ctr;
    FILE * fptr;
   
    fptr   = efopen ( options->epi_file_name, "r" );
    if (! fptr ) return 1;
    
    memset ( epitope, 0, protein->length*sizeof(int));
    total = 0;
    while ( ! feof ( fptr) ) {
 	if ( fgets (buf, LONGSTRING, fptr ) && strlen(buf) > 1 ) {
	    if ( buf[0] == '#' ) continue;
	    if (! (sscanf (buf, "%s %c", aux, &res_type)  ) ) {
		fprintf (stderr, "Error reading %s.\n", options->epi_file_name);
		fclose (fptr);
		return 1;
	    }
	    found = 0;
	    for ( res_ctr=0; res_ctr < protein->length; res_ctr++) {
		if ( ! strcmp ( protein->sequence[res_ctr].pdb_id, aux) ) {
		    if ( res_type !=  protein->sequence[res_ctr].res_type_short ) {
			fprintf (stderr, "While reading %s: Type mismatch for res id %s (pdb:%c, here:%c) .\n",
				 options->epi_file_name, aux, protein->sequence[res_ctr].res_type_short, res_type);
			fclose (fptr);
			return 1;
		    }
		    found = 1;
		    epitope[res_ctr]= 1;
		    total++;
		}
	    }
	    if ( ! found ) {
		fprintf (stderr, "While reading %s: Residue with the id %s not found.\n",
			 options->epi_file_name, aux);
		fclose (fptr);
		return 1;
	    }
	}
    }
    
    if ( !total ) {
	fprintf (stderr, "No positions matching pdb found in %s.\n", options->epi_file_name);
	return 1;
    }

    fclose (fptr);
    return 0;
}
