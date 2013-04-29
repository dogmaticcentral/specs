# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>

# define DEFAULT_MAX_ID 0.98
# define DEFAULT_MIN_ID 0.30
# define DEFAULT_MIN_QRY_ID 0.30
# define BUFFLEN  150
# define ALMT_NAME_LENGTH 50

# define WINDOW_LENGTH 20

typedef struct{
    int number_of_seqs;
    int length;
    char ** sequence;
    char ** name;
    int * seq_gaps;
    int * column_gaps;
    int * marked;
}  Alignment;

/************************************************/
int main ( int argc, char * argv[]) {

    char infile[BUFFLEN]= {'\0'}, outfile[BUFFLEN]= {'\0'};
    char query[ALMT_NAME_LENGTH];
    double max_id, min_id, min_qry_id;
    Alignment alignment;
    int read_afa (char *filenname, Alignment * alignment);
    int compare_seqs (Alignment * alignment, char *query, double max_id,
		      double min_id, double min_qry_id);
    int write_afa (char *filename, Alignment * alignment);
   
    if ( argc < 4 ) {
	printf (
	    "Usage: %s  <infile>  <outfile>  <query> [<max_id> [<min_id> ][<min_id_qry>] ].\n",
	    argv[0]);
	exit (1);
    }

    memset ( &alignment, 0, sizeof(Alignment) );
    
    sprintf (infile, "%s", argv[1]);
    sprintf (outfile,"%s", argv[2]);
    max_id = DEFAULT_MAX_ID;
    min_id = DEFAULT_MIN_ID;
    min_qry_id = DEFAULT_MIN_ID;
    if (argc > 3 ) sprintf (query, "%s", argv[3]);
    if (argc > 4 ) max_id = atof(argv[4]);
    if (argc > 5 ) min_id = atof(argv[5]);
    if (argc > 6 ) min_qry_id = atof(argv[5]);


    if ( read_afa( infile, &alignment) ) exit (1);
    compare_seqs (&alignment, query, max_id, min_id,  min_qry_id );
    if ( write_afa ( outfile, &alignment) ) exit (1);

    
     
    return 0;
}



/************************************************/
int compare_seqs (Alignment * alignment, char *query, double max_id, double min_id, double min_qry_id){
    
    int seq_ctr_1,  seq_ctr_2, ctr;
    int is_query;
    int window_length;
    char * seq1,  *seq2;
    double pid, avg_pid, stdev, * pid_w_query = NULL;
    int count_gaps (Alignment * alignment);
    int window_dist (char * seq1, char *seq2, int seq_length,
		     int window_length, double *avg_pid_ptr, double * stdev_ptr) ;

    if ( count_gaps(alignment ) ) exit(1);

    if ( ! (  pid_w_query =  calloc (alignment->number_of_seqs, sizeof(double) ) ) ) {
	fprintf ( stderr, "Error: calloc error in compare_seqs().\n");
	exit (1);
    }
    
    is_query = -1;
    for (seq_ctr_1=0; seq_ctr_1 < alignment->number_of_seqs; seq_ctr_1++ ) {
	if ( ! strncmp (alignment->name[seq_ctr_1], query, strlen(query) ) &&
	     ! strncmp (alignment->name[seq_ctr_1], query, strlen(alignment->name[seq_ctr_1]) ) ) {
	    is_query = seq_ctr_1;
	    break;
	}
    }


    if ( is_query < 0 ) {
	fprintf ( stderr, "Error: query %s not found\n", query);
	exit (1);
    }
    
    /* round 1: find distance to query, and drop all with dist to query < min_id_query */ 
    /* also can remove too similar, to make it faster later */
    seq_ctr_1 = is_query;
    seq1 = alignment->sequence[seq_ctr_1];
    for (seq_ctr_2=0; seq_ctr_2 < alignment->number_of_seqs; seq_ctr_2++ ) {
	if ( seq_ctr_2 == is_query ) continue;
	if ( alignment->marked[seq_ctr_2] ) continue;
	seq2 = alignment->sequence[seq_ctr_2];
	window_dist (seq1, seq2,  alignment->length, window_length=alignment->length, &pid, &stdev);
	pid_w_query[seq_ctr_2] = pid;
	if ( pid  > max_id ) {
	    printf ( " near query:  %s  %s  %8.2lf  %8.2lf \n", alignment->name[seq_ctr_1],
		     alignment->name[seq_ctr_2], pid, max_id);
	    alignment->marked[seq_ctr_2] = 1;
	    continue;
	}
	window_dist (seq1, seq2,  alignment->length, WINDOW_LENGTH, &avg_pid, &stdev);
	if ( avg_pid < min_qry_id ) {
	    printf ( " far from query: %10s  %10s  %8.2lf  %8.2lf  %8.2lf \n",
		     alignment->name[seq_ctr_1],  alignment->name[seq_ctr_2],
		     pid, avg_pid, min_qry_id);
	    alignment->marked[seq_ctr_2]  = 1;
	}
	
    }
 

    /* round 2: for each pair if dist < min_id, take out the sequence more distant from the query */ 
	
    for (seq_ctr_1=0; seq_ctr_1 < alignment->number_of_seqs; seq_ctr_1++ ) {
	if ( seq_ctr_1 == is_query ) continue;
	if ( alignment->marked[seq_ctr_1] ) continue;
	seq1 = alignment->sequence[seq_ctr_1];
	
	for (seq_ctr_2=seq_ctr_1+1; seq_ctr_2 < alignment->number_of_seqs; seq_ctr_2++ ) {
	    if ( seq_ctr_2 == is_query ) continue;
	    if ( alignment->marked[seq_ctr_2] ) continue;
	    seq2 = alignment->sequence[seq_ctr_2];
	    
		
	    window_dist (seq1, seq2,  alignment->length, window_length=alignment->length, &pid, &stdev);
	    window_dist (seq1, seq2,  alignment->length, WINDOW_LENGTH, &avg_pid, &stdev);
		     
	    if ( avg_pid > max_id ) {
		printf ( "near: %10s  %10s  %8.2lf  %8.2lf  %8.2lf \n",
			 alignment->name[seq_ctr_1],  alignment->name[seq_ctr_2],
			 pid, avg_pid, max_id);
		if ( pid_w_query[seq_ctr_2] < pid_w_query[seq_ctr_1] ) {
		    alignment->marked[seq_ctr_1] = 1;
		} else {
		    alignment->marked[seq_ctr_2]  = 1;
		}
		
	    } else if ( avg_pid < min_id ) {
		printf ( " far: %10s  %10s  %8.2lf  %8.2lf  %8.2lf \n",
			 alignment->name[seq_ctr_1],  alignment->name[seq_ctr_2],
			 pid, avg_pid, min_id);
		if ( pid_w_query[seq_ctr_2] > pid_w_query[seq_ctr_1] ) {
		    alignment->marked[seq_ctr_1] = 1;
		} else {
		    alignment->marked[seq_ctr_2]  = 1;
		}
	    } 
	    
	    
	}
	
    }

    ctr = 0;
    for (seq_ctr_1=0; seq_ctr_1 < alignment->number_of_seqs; seq_ctr_1++ ) {
	ctr += alignment->marked[seq_ctr_1];
    }
    printf ( " removed %d sequences   remaining %d\n", ctr,  alignment->number_of_seqs-ctr); 
    
    
    free( pid_w_query);
    
    return 0;
}

/************************************************/


int window_dist (char * seq1, char *seq2, int seq_length,int window_length,
		 double *avg_pid_ptr, double * stdev_ptr) {

    double pid, avg_pid, avg_pid_sq,  stdev;
    int begin, pos, common_length, identical_length;
    int window_ctr, all_gaps;

    avg_pid    = 0.0;
    avg_pid_sq = 0.0;
    window_ctr = 0;
    for (begin=0; begin <= seq_length - window_length; begin++) {
	common_length = identical_length = 0;
	all_gaps = 0;
	for (pos=begin; pos < begin+window_length; pos++ ) {
	    if ( seq1[pos] == '.' || seq2[pos] == '.' ) {
		all_gaps += ( seq1[pos] == seq2[pos] );
		continue;
	    }
	    common_length ++;
	    if ( seq1[pos] ==  seq2[pos] ) identical_length++;
	}
	if ( all_gaps == window_length ) continue;
	pid = common_length ? (double)identical_length/common_length:0.0;
	avg_pid += pid;
	avg_pid_sq += pid*pid;
	window_ctr ++;
    }
    if ( window_ctr ) {
	avg_pid /= window_ctr;
	avg_pid_sq /= window_ctr;
	stdev = avg_pid_sq - avg_pid*avg_pid;
	stdev = (stdev>0) ? sqrt(stdev):0.0;
    } else {
	stdev = 0.0;
    }

    *avg_pid_ptr = avg_pid;
    *stdev_ptr   = stdev;
    
    return 0;
}

/************************************************/
int count_gaps (Alignment * alignment) {

    int seq_ctr, pos_ctr;
    void * emalloc(int	size);
    if ( ! alignment->seq_gaps  ) {
	alignment->seq_gaps    = (int *) emalloc (alignment->number_of_seqs*sizeof(int));
	if (!alignment->seq_gaps) return 1;
    }
    if ( ! alignment->column_gaps) {
	alignment->column_gaps = (int *) emalloc (alignment->length*sizeof(int));
	if (!alignment->column_gaps) return 1;
    }

    memset ( alignment->seq_gaps, 0, (alignment->number_of_seqs*sizeof(int)) );
    memset ( alignment->column_gaps, 0, (alignment->length*sizeof(int)) );
		
    
    for ( seq_ctr=0; seq_ctr<alignment->number_of_seqs; seq_ctr++ ) {
	if ( alignment->marked[seq_ctr] ) continue;
	for ( pos_ctr=0; pos_ctr<alignment->length; pos_ctr++) {
	    if ( alignment->sequence[seq_ctr][pos_ctr] == '.' ) {
		alignment->column_gaps[pos_ctr] ++;
		alignment->seq_gaps[seq_ctr] ++;
	    }
	}
    }
    return 0;
}

/************************************************/
int write_afa (char *filename, Alignment * alignment) {

    FILE * fptr = NULL;
    int seq_ctr, pos_ctr, c, pruned_no_seqs;
    FILE * efopen(char * name, char * mode);
    /* open file */
    fptr = efopen ( filename, "w");
    if ( !fptr ) return 1;

    count_gaps (alignment);
    
    pruned_no_seqs = 0;
    for (seq_ctr=0; seq_ctr < alignment->number_of_seqs; seq_ctr++ ) {
	if ( alignment->marked[seq_ctr] ) continue;
	pruned_no_seqs++;
    }
    
    for (seq_ctr=0; seq_ctr < alignment->number_of_seqs; seq_ctr++ ) {
	if ( alignment->marked[seq_ctr] ) continue;
	fprintf (fptr,  ">%s\n", alignment->name[seq_ctr]);
	
	c = 0; /* counter for newline insertion */ 
	for (pos_ctr=0; pos_ctr < alignment->length; pos_ctr++) {
	    
	    if ( alignment->column_gaps[pos_ctr] == pruned_no_seqs ) continue;
	    c++;  if ( c && !(c%50) ) fprintf (fptr, "\n");
	    fprintf ( fptr, "%c", alignment->sequence[seq_ctr][pos_ctr] );
	}
	fprintf (fptr, "\n");
    }
 
    fclose (fptr);
    return 0;
}

/************************************************/
int read_afa (char *filename, Alignment * alignment) {

    FILE * fptr = NULL;
    char line[BUFFLEN];
    int number_of_seqs, almt_length, ctr;
    int reading, seq_ctr;
    int seq_pos, begin_name;
    int * marked = NULL;
    char ** sequence;
    char ** name;
    char **chmatrix(int rows, int columns);
    FILE * efopen(char * name, char * mode);
    void * emalloc(int	size);
    /* open file */
    fptr = efopen ( filename, "r");
    if ( !fptr ) return 1;
    
    /* find the alignment length  and number of seqs info */
    almt_length = 0;
    number_of_seqs = 0;
    reading = 0;
    ctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( line[0] == '>' ){
	    number_of_seqs++;
	    reading = (! almt_length);
	} else if ( reading ) {
	    ctr =0;
	    while  ( ctr < BUFFLEN && line[ctr] != '\n' )  {
		if ( !isspace (line[ctr]) ) almt_length ++;
		ctr ++;
	    }
	}
    }
    printf ( "number of seqs = %d\n",   number_of_seqs);
    printf ( "almt_length = %d\n", almt_length);
    
     /* allocate */
    sequence = chmatrix (number_of_seqs, almt_length);
    if ( !sequence ) return 1;
    name     = chmatrix (number_of_seqs, ALMT_NAME_LENGTH);
    if ( !name ) return 1;
    marked = (int * )emalloc ( number_of_seqs*sizeof(int));
    if ( !marked ) return 1;
    
    
    /* read in */
    rewind(fptr);
    seq_ctr = -1;
    seq_pos = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( line[0] == '>' ){
	    seq_pos = 0;
	    seq_ctr++;
	    
	    /* chomp*/
	    ctr = 1;
	    while ( ctr < BUFFLEN && line[ctr] != '\n' ) ctr++;
	    if ( ctr < BUFFLEN ) line[ctr] = '\0';
	    
	    ctr = 1;
	    while ( isspace (line[ctr]) ) ctr++;
	    begin_name = ctr;
	    while ( !isspace (line[ctr]) ) ctr++;
	    line[ctr] = '\0';
	    /* make sure the name is not too long */
	    if ( strlen ( &line[ctr]) > ALMT_NAME_LENGTH ) {
		line[ALMT_NAME_LENGTH+ctr-1] = '\0';
	    }
	    sprintf ( name[seq_ctr], "%s", &line[begin_name]);
	} else  {
	    ctr =0;
	    while  ( ctr < BUFFLEN && line[ctr] != '\n' )  {
		if ( !isspace (line[ctr]) ) {
		    if ( seq_pos >= almt_length ) {
			fprintf (stderr, "Error: sequence %s longer than the first sequence.\n",
				 name[seq_ctr]);
			fprintf (stderr, "[Is %s really an alignment (afa),", filename);
			fprintf (stderr, " or just a sequence file (fasta)?].\n");
			exit (1);
		    }
		    if (line[ctr] == '-' ) {
			sequence[seq_ctr][seq_pos] = '.';
		    } else {
			sequence[seq_ctr][seq_pos] = line[ctr];
		    }
		    seq_pos ++;
		}
		ctr ++;
	    }
	}
    }

    
   
    /* return values */
    alignment->number_of_seqs = number_of_seqs;
    alignment->length         = almt_length;
    alignment->sequence       = sequence;
    alignment->name           = name;
    alignment->marked           = marked;

    
    fclose (fptr);
    
    return 0;
}


/************************************************/
char **chmatrix(int rows, int columns){
    char **m;
    int i;
        /* allocate pointers to rows */
    m=(char **) malloc(rows*sizeof(char*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(char *) calloc( rows*columns, sizeof(char));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}

/************************************************/
FILE * efopen(char * name, char * mode)
{

    FILE * fp;


    if ((fp = fopen(name, mode)) == NULL) {
	fprintf (stderr,  
	      "Cannot open \"%s\" for \"%s\"\n", name, mode);
	return NULL;
    }

    return fp;

}



/************************************************/
void * emalloc(int	size)
{
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %u bytes", size);
	return NULL;
    }

    return ptr;
}


