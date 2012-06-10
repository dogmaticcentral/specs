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
