# include "specs.h"


int set_similarity ( Options * options, int * similar_to ) {

    int i;

    for(i=0;i<ASCII;i++){
	similar_to[i] = i;
    }
# if 1
    similar_to['I'] = 'V';
    similar_to['L'] = 'V';
    similar_to['S'] = 'T';
    similar_to['D'] = 'E';
    similar_to['K'] = 'R';
    similar_to['Q'] = 'N';
    similar_to['.'] = '.';

    similar_to['A'] = 'V';
    similar_to['M'] = 'V';
    similar_to['G'] = 'V';
    similar_to['F'] = 'Y';
    similar_to['H'] = 'R';
# else

    similar_to['A'] = 'V';
    similar_to['I'] = 'V';
    similar_to['L'] = 'V';
    similar_to['G'] = 'V';
    similar_to['S'] = 'V';
    similar_to['T'] = 'V';
    similar_to['C'] = 'V';
    similar_to['M'] = 'V';
    
# endif   

    return 0;
}
