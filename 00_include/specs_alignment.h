# define ALMT_NAME_LENGTH 30

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
    char * refseq;
    char refseq_name[ALMT_NAME_LENGTH];
    int  refseq_gaps;
}  Alignment;
