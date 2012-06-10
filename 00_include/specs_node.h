# define NODE_PRINT( fptr, node )                      \
    fprintf ( fptr, "id:%3d     type:%3d      value:%10p    parent:%10p    left:%10p  right:%10p \n",\
	    node ->id,  node ->type, node, node ->parent,  node ->left,  node ->right);\

# define NAME_LENGTH 30
# define ASCII 128 /* number of characters */

# define ROOT  1
# define INNER 2
# define LEAF  4
# define CENTRAL 8 /*for handling unrooted trees*/  

# define UP    1
# define DOWN  2
# define LEFT  4
# define RIGHT 8

typedef struct Node {

    struct Node *left, *right, *parent;
    int id;
    int type;
    int number_of_leaves;  /* in the corresponding subtree */
    int marked;
    int  bin; /*similarity bin the node belongs to*/
    double dist_to_parent, dist_to_left, dist_to_right;
    double avg_sim;
    char * seq;
    char * name;
    char consensus;
    int consensus_length;
    int population[ASCII]; /* population of amino acid types */
    double entropy;
    double entropy_complement;

} Node;

typedef struct Tree{

    Node *root;
    Node *leaf; /* this is actually the beginning of the node storage array */
    Node  ***group_root,  ***group_root_sim; /* pointers to nodes which form groups at given rank */
    int size; /* leaves + inner nodes */
    int no_of_leaves;

} Tree;




# define NODE_PRINT( fptr, node )                      \
    fprintf ( fptr, "id:%3d     type:%3d      value:%10p    parent:%10p    left:%10p  right:%10p \n",\
	    node ->id,  node ->type, node, node ->parent,  node ->left,  node ->right);\

Node ***node_matrix(int rows, int columns);
double  node_distance ( double **seq_dist, Node* node1, Node* node2 );
