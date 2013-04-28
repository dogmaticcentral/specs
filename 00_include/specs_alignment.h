# define ALMT_NAME_LENGTH 150

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
    char ** refseq;
    char **refseq_name;
    int  no_refseqs;
    int * protected_position;
    int number_of_protected_positions;
    int refseq_gaps;
}  Alignment;
