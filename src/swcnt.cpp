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

#include "swcnt.hpp"

u64 binom[MAXLEN][MAXLEN];
Cache swcnt[MAXLEN][MAXLEN]; // indices are #bits
u64 parallel_mask = -1; // -1 by default means serial mode. 0 to set means parallel mode

// precompute the table
void binom_precompute(){
    // clean the table
    for(int i = 0; i < MAXLEN; i++)
        for(int j = 0; j < MAXLEN; j++)
            binom[i][j] = 0;
    // compute
    binom[0][0] = 1;
    for(int i = 1; i < MAXLEN; i++){
        binom[i][i] = binom[i][0] = 1;
        for(int j = 1; j < i; j++){
            binom[i][j] = binom[i - 1][j - 1] + binom[i - 1][j];
        }
    }
    return;
}

// extraction of precomputed results
u64 binomial(int i, int j){
    return binom[i][j];
}

// build struct word with external array
Word build_word(u64 wordbin, int len, Runtab tab){
    Word w;
    w.bits = wordbin;
    w.len = len;
    w.run = tab;
    // first reverse the bits so that it is easier later
    u64 rev = 0;
    for(int i = 0; i < len; i++){
        rev <<= 1;
        rev += wordbin & 1;
        wordbin >>= 1;
    }
    // then we count the runs using intrinsics
    if(rev & 1) rev = ~rev; // always keeping rightmost run 0
    int unread = len;
    for(w.runcnt = 0; ; w.runcnt++) { // from right to left in runs here
        int zcnt = std::countr_zero(rev);
        if(zcnt >= unread){ // all zero
            w.run[w.runcnt] = unread;
            w.runcnt++;
            break;
        }
        w.run[w.runcnt] = zcnt;
        unread -= zcnt;
        rev >>= zcnt;
        rev = ~rev; // exchange 0 and 1, keep rightmost run 0
    }
    return w;
}

// print the word, does not use the field "run"
void print_word_bin(Word word){
    // convert binary representation
    char str[65];
    u64 bits = word.bits;
    int len = word.len;
    for(int i = len - 1; i >= 0; i--){
        str[i] = '0' + (bits & 1);
        bits >>= 1;
    }
    str[len] = 0;
    printf("%s\n", str);
    return;
}

// debug purpose
void print_word(Word* word){
    // convert binary representation
    char str[65];
    u64 bits = word->bits;
    int len = word->len;
    for(int i = len - 1; i >= 0; i--){
        str[i] = '0' + (bits & 1);
        bits >>= 1;
    }
    str[len] = 0;
    printf("Binary: %s\nRuns: ", str);
    for(int i = 0; i < word->runcnt; i++) printf("%d ", word->run[i]);
    printf("\n");
    return;
}

// debug purpose
bool is_equal_word(Word* word1, Word* word2){
    if(!(word1->bits == word2->bits 
      && word1->len  == word2->len
      && word1->runcnt  == word2->runcnt
    )) return false;
    for(int i = 0; i < word1->runcnt; i++){
        if(word1->run[i] != word2->run[i]) return false;
    }
    return true;
}

// increment a word, we suppose that the word here is not a segment
bool increment_word(Word* word){
    word->bits++;
    if(word->runcnt & 1){
        // odd number of segments, thus ending with 0
        if(word->run[word->runcnt - 1] > 1){
            // more than one 0, new run of 1 of size 1
            word->run[word->runcnt - 1]--;
            word->run[word->runcnt] = 1;
            word->runcnt++;
        } else {
            // only one 0, fusion with previous run of 1
            word->runcnt--;
            word->run[word->runcnt - 1]++;
        }
    } else {
        // even number of segments, thus ending with 1
        if(word->run[word->runcnt - 2] > 1){ // at least two runs
            // more than one 0 in the segment before, the last run is now zero
            // new run of size 1
            word->run[word->runcnt] = word->run[word->runcnt - 1];
            word->run[word->runcnt - 1] = 1;
            word->run[word->runcnt - 2]--;
            word->runcnt++;
        } else {
            // only one 0 in the segment before, fusion with previous run
            // test ending condition: 01111...111
            if(word->runcnt == 2) return false;
            // at least 4 runs
            word->run[word->runcnt - 3]++;
            word->run[word->runcnt - 2] = word->run[word->runcnt - 1];
            word->runcnt--;
        }
    }
    return true;
}

// increment a word, assuming starting and ending with 0, and we add 2
bool increment_word_2(Word* word){
    if(!increment_word(word)) return false;
    return increment_word(word);
}

void add_bit(Word* word, int bit){
    // normalize
    bit &= 1;
    // add the bit and extend the length
    word->bits <<= 1;
    word->bits += bit;
    word->len++;
    // adjust the runs
    if((word->runcnt & 1) == bit){ // different letter
        word->run[word->runcnt] = 1;
        word->runcnt++;
    } else {
        word->run[word->runcnt - 1]++;
    }
    return;
}

void remove_bit(Word* word){
    word->bits >>= 1;
    word->len--;
    word->run[word->runcnt - 1]--;
    if(word->run[word->runcnt - 1] == 0) word->runcnt--;
    return;
}

// compute the index of the first run of w after the leftmost possible
// occurrences of sw in w
inline static int count_run_idx(int* wrun, int wruncnt, int* swrun, int swruncnt){
    int idx = 0;
    int sidx = 0;
    int curswrun = swrun[sidx];
    while(sidx < swruncnt){
        if(idx >= wruncnt) return -1;
        if(wrun[idx] < curswrun){
            curswrun -= wrun[idx];
            idx += 2;
        } else {
            sidx++;
            idx++;
            curswrun = swrun[sidx];
        }
    }
    return idx;
}

// compute the index of the first run of w after the leftmost possible
// occurrences of sw in w
inline static int count_run_idx_rev(int* wrun, int wruncnt, int* swrun, int swruncnt){
    int idx = wruncnt - 1;
    int sidx = swruncnt - 1;
    int curswrun = swrun[sidx];
    while(sidx >= 0){
        if(idx < 0) return -1;
        if(wrun[idx] < curswrun){
            curswrun -= wrun[idx];
            idx -= 2;
        } else {
            sidx--;
            idx--;
            curswrun = swrun[sidx];
        }
    }
    return idx;
}

// cut the word at index run, and return the first part
static inline Word cut_word_front(Word w, int run){
    Word neww = w;
    int accu = 0;
    for(int i = 0; i < run; i++) accu += w.run[i];
    neww.bits >>= w.len - accu; // only on the left that remains
    // run table pointer does not change
    neww.runcnt = run;
    neww.len = accu;
    return neww;
}

// cut the word at index run, and return the second part
static inline Word cut_word_back(Word w, int run){
    Word neww = w;
    int accu = 0;
    for(int i = run; i < w.runcnt; i++) accu += w.run[i];
    neww.bits &= (1ULL << accu) - 1; // only the right that remains
    neww.run += run;
    neww.runcnt = w.runcnt - run;
    neww.len = accu;
    return neww;
}

// assuming w and sw starts with the same letter, and end also the same
static u64 subword_cnt_raw(Word w, Word sw, int orig_wlen){
    if(sw.runcnt == 0) return 1; // empty subword
    if(w.runcnt < sw.runcnt) return 0; // not enough run
    u64 accu = 0;
    // cut the subword into two
    int mididx = sw.runcnt / 2;
    int midseg = sw.run[mididx];
    // compute the runs of w that sw may span, get left and right index for wl
    int lidx = count_run_idx(w.run, w.runcnt, sw.run, mididx);
    int ridx = count_run_idx_rev(w.run, w.runcnt, 
                                 sw.run + mididx + 1, sw.runcnt - mididx - 1);
    if(lidx > ridx || lidx < 0 || ridx >= w.runcnt){
        // debug info
        /* 
        printf("Left and right index check failed: %d, %d\n", lidx, ridx);
        print_word(&w);
        print_word(&sw);
        */
        return 0; // not possible
    }
    // lookup
    if(w.runcnt < MAX_CACHE_RUN){
        auto search = swcnt[w.len][sw.len].find(u64pair(w.bits, sw.bits));
        if(search != swcnt[w.len][sw.len].end()) return search->second;
    }
    // cut in the middle and recursion (divide and conquer)
    // we look at where the middle segment could span, 
    // then cut the word and the subword into have and do recursion
    Word swfront = cut_word_front(sw, mididx);
    Word swback = cut_word_back(sw, mididx + 1);
    if(lidx == ridx){
        accu = binom[w.run[lidx]][midseg]; // middle span
        accu *= subword_cnt_raw(cut_word_front(w, lidx), swfront, orig_wlen);
        accu *= subword_cnt_raw(cut_word_back(w, lidx + 1), swback, orig_wlen); 
    } else {
        for(int k = lidx; k < ridx + 2; k += 2){
            int wsegtotal = 0;
            for(int l = k; l < ridx + 2; l += 2){
                /*
                int wsegtotal = 0;
                for(int i = k; i < l + 2; i += 2) wsegtotal += w.run[i];
                */
                wsegtotal += w.run[l];
                int wsegin = wsegtotal - w.run[k] - w.run[l];
                int64_t mult = 0;
                // inclusion-exclusion principle for the middle span
                // if k == l, only the first term will be non-zero
                mult += binom[wsegtotal][midseg];
                mult -= binom[wsegin + w.run[k]][midseg];
                mult -= binom[wsegin + w.run[l]][midseg];
                if(wsegin >= 0) mult += binom[wsegin][midseg];
                if(mult > 0){
                    mult *= subword_cnt_raw(cut_word_front(w, k), swfront, orig_wlen);
                    mult *= subword_cnt_raw(cut_word_back(w, l + 1), swback, orig_wlen);
                    accu += mult;
                }
            }
        }
    }
    // restriction on length to limit memory usage and control for modification for parallelism
    if(w.runcnt < MAX_CACHE_RUN && w.len <= orig_wlen){
        swcnt[w.len][sw.len].insert({u64pair(w.bits, sw.bits), accu});
    }
    // debug info
    /*
    printf("Subcase result: %lu, %d, %d\n", accu, lidx, ridx);
    print_word(&w);
    print_word(&sw);
    // print_word(&swfront);
    // print_word(&swback);
    */
    return accu;
}

// count subword occurrences
u64 subword_cnt(Word word, Word subword){
    // get the words with the same tail
    if(((word.bits >> (word.len - 1)) & 1) != (subword.bits >> (subword.len - 1)) & 1){
        word.len -= word.run[0];
        word.run++;
        word.runcnt--;
    }
    // get the words with the same head
    if((word.bits & 1) != (subword.bits & 1)){
        word.len -= word.run[word.runcnt - 1];
        word.bits >>= word.run[word.runcnt - 1];
        word.runcnt--;
    }
    return subword_cnt_raw(word, subword, word.len & parallel_mask);
}

// set to be used in a parallel way, i.e., forbidding insertions in cache
void set_parallel_mode(){
    parallel_mask = 0;
    return;
}