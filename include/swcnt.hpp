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

#ifndef __SWCNT__
#define __SWCNT__

#include <stdint.h>
#include <stdio.h>
#include <tuple>
//#include <unordered_map>
#include <tbb/concurrent_unordered_map.h>
#include <vector>
#include <bit>
#include <utility>

#define MAXLEN 64
#define MAX_CACHE_RUN 8

typedef uint64_t u64;

typedef std::pair<u64, u64> u64pair;

typedef int Runtab[MAXLEN];

struct pairhash {
    public:
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const {
        return std::hash<T>()(x.first) * 23 + std::hash<U>()(x.second);
    }
};

// typedef std::unordered_map<u64pair, u64, pairhash> Cache; // Cache for subword counting

typedef tbb::concurrent_unordered_map<u64pair, u64, pairhash> Cache; // Cache for subword counting

typedef struct{
    u64 bits;
    int* run;  // run lengths from left to right, leftmost run always contains 0
               // should not use this field outside the called function
               // as most of the time it will be local variables
    int runcnt;
    int len;   // number of bits
} Word;

// record for subword occurrences of a given word
typedef struct {
    Word word;
    std::vector<Word> subwords;
    u64 occ;
} Rec_sw;

typedef struct {
    std::vector<Rec_sw> recs;
    u64 occ;
} Rec_occ;

// precompute binomial coefficients
void binom_precompute();

// return the precomputed binomial coefficients
u64 binomial(int i, int j);

// build the struct Word with an external array
Word build_word(u64 wordbin, int len, int* tab);

// debug purpose
void print_word(Word* word);

// printing the word, does not use the field "run"
void print_word_bin(Word word);

// debug purpose
bool is_equal_word(Word* word1, Word* word2);

// increment the word as binary representation by 1
// returns false if cannot continue
bool increment_word(Word* word);

// increment the word as binary representation by 2
// returns false if cannot continue
bool increment_word_2(Word* word);

// add a bit at the end
void add_bit(Word* word, int bit);

// remove the bit at the end
void remove_bit(Word* word);

// returns the number of subword occurrences
u64 subword_cnt(Word word, Word subword);

#endif
