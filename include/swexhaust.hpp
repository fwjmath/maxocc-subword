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

#ifndef __SWEXHAUST__
#define __SWEXHAUST__

#include "swcnt.hpp"
#include "fibogen.hpp"
#include <map>

// number of threads used in parallel mode, should always be a power of 2
#ifndef THREAD_COUNT
#define THREAD_COUNT 4
#endif

typedef std::map<u64, u64> Histogram;

// thread information for parallelism
typedef struct {
    int n;
    int thread_id;
    u64 record;
    Rec_occ* minrec;
} Thread_info;

// Returns the words with the minimal value of most frequence occurrences
// Each word may have several subwords reaching the same number of occurrences
// And we may have several words with the same numbers
Rec_occ min_maxfreq_subword_hinted(int n, u64 record);

// The same as the function above, but only for some of the words
// Used for the parallel version
void* min_maxfreq_subword_hinted_parallel(void* info);

// compute the maxfreq for subwords in a given word. Used in metaheuristics.
Rec_sw maxfreq_subword_hinted_fast(Word w, u64 record);

// compute most frequent subwords for a single given word. Used in computing for a single word.
Rec_sw maxfreq_subword_single(Word w, u64 record);

// returns histogram of max subword occurrences
Histogram maxfreq_subword_histo(int n);

// compute max freq subword of some lenghts, for metaheuristics.
u64 maxfreq_subword_fast(Word w);

// for metaheuristics
Rec_sw maxfreq_subword_hinted_fast(Word w, u64 record);

#endif