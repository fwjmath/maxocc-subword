/******************************************************************************

    swmain - a research tool for minimizing maximal subword occurrences
    
    Copyright (C) 2024 Wenjie Fang <fwjmath@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "swcnt.hpp"
#include "swexhaust.hpp"
#include "swutils.hpp"
#include "swmeta.hpp"

/*
Four modes of operations:

1. Finding words with minimal maxocc of subwords using hinted exhaustive search
2. Metaheuristic search to obtain reasonable hint
3. Histogram of maxocc of subwords
4. Compute (an underestimation of) the best maxocc of "increments" of a given word by inserting a bit
*/

int main(int argc, char** argv){
    binom_precompute();
    
    time_t mytime = time(NULL);
    printf("%s", ctime(&mytime));
    if(argc <= 1){
        printf("Needs at least the number of bits\n");
        return 0;
    }
    int n = atoi(argv[1]);
    if(n <= 0 || n > 64){
        printf("Invalid argument, the number of bits is between 1 and 64\n");
    }
    
    u64 hint = 0;
    bool computed = false;
    if(argc >= 3){
        if(strcmp(argv[2], "histo") == 0){
            printf("Producing histogram for maxocc with %d bits.\n", n);
            histo_subword(n);
            computed = true;
        }else if(strcmp(argv[2], "meta") == 0){
            if(argc <= 5){
                printf("Insufficient arguments.\n");
                printf("Needs exhaustive search radius, max sampling number.");
            }
            printf("Metaheuristic search for hint with with %d bits.\n", n);
            mixed_descent(n, atoi(argv[3]), atoi(argv[4]));
            computed = true;
        }else if(strcmp(argv[2], "word") == 0){
            if(argc < 4 || strlen(argv[3]) != n){
                printf("Invalid argument for the word.\n");
            }else{
                printf("Maxocc for single word.\n");
                compute_maxfreq_subword(argv[3]);
            }
            computed = true;
        }else if(strcmp(argv[2], "insert") == 0){
            if(argc < 4 || strlen(argv[3]) != n - 1){
                printf("Invalid argument for the word.\n");
            }else{
                printf("Heuristic by inserting a bit.\n");
                insert_heuristic(argv[3]);
            }
            computed = true; 
        }else{
           hint = atoi(argv[2]);
        }
    }
    if(!computed){
        if(hint == 0){
            printf("Invalid hint. Need an over-estimation of minimal maxocc.\n");
            printf("Using default hint 2^n.\n");
            hint = 1;
            hint <<= n;
        }
        hinted_search(n, hint);
    }
    mytime = time(NULL);
    printf("%s", ctime(&mytime)); 
    return 0;
}
