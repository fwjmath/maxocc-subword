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

#ifndef __SWUTILS__
#define __SWUTILS__

#include <stdio.h>
#include <string.h>
#include "swcnt.hpp"
#include "swexhaust.hpp"

// build a word according to a 0-1 string
Word build_word_str(const char* str, bool is_subword);

// exhaustive search for minimal subword entropy, using a hint
// using subword with large number of occurrences from the last word as a hint
void hinted_search(int n, u64 hint);

// exhaustive search for minimal subword entropy, using a hint, branch and bound
void pruned_search(int n, u64 hint);

// build a histogram for subword occurrences
void histo_subword(int n);

// print record
void print_record(Rec_sw* minrec);

#endif