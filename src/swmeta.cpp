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

#include "swmeta.hpp"

// initialize combination represented by bits, all on the left
static inline u64 init_comb(int n, int k){
    return ((1ULL << k) - 1) << (n - k);
}

// get the next combination using a well-known algorithm (the same as in Python)
static inline bool next_comb(int n, int k, u64* cur){
    // get the position of the rightmost 1
    int rpos = std::countr_zero(*cur);
    // if we can move it to the right, then do it
    if(rpos > 0){
        *cur -= (1ULL << (rpos - 1));
        return true;
    }
    // otherwise, move the whole pack of 1, first count their numbers
    int ocnt = std::countr_zero(~(*cur));
    if(ocnt == k) return false;
    // find the next 1, move it left and attach the pack
    u64 ncur = *cur; 
    ncur -= (1ULL << ocnt) - 1;
    int npos = std::countr_zero(ncur);
    *cur = ncur - (1ULL << (npos - 1 - ocnt));
    return true;
}

// debug purpose
void check_comb(int n, int k){
    u64 comb = init_comb(n, k);
    u64 oldcomb = comb;
    u64 cnt = 1;
    while(next_comb(n, k, &comb)) {
        if(oldcomb <= comb){
            printf("Error at %lu\n", oldcomb);
            break;
        }
        cnt++;
    }
    printf("n: %d, k: %d, total: %lu\n", n, k, cnt);
    return; 
}

// search for the whole neighborhood, hinted by record and previous subwords
static inline Rec_sw local_search(Word w, int k, u64 record){
    Rec_sw minrec = {w, std::vector<Word>(), record};
    int n = w.len;
    int tab[MAXLEN];
    bool flag = true;
    Runtab swruns;
    while(flag) {
        u64 comb = init_comb(n - 1, k);
        u64 recbits = w.bits;
        flag = false;
        do {
            Word curw = build_word(w.bits ^ comb, n, swruns);
            Rec_sw maxrec = maxfreq_subword_hinted_fast(curw, record);
            if(maxrec.occ < minrec.occ){
                minrec = maxrec;
                record = minrec.occ;
                flag = true;
                break;
            }
        } while(next_comb(n - 1, k, &comb));
    }
    minrec.word = build_word(minrec.word.bits, n, swruns);
    return minrec;
}

static inline Rec_sw local_search_full(Word w, int k, u64 record){
    Rec_sw minrec = {w, std::vector<Word>(), record};
    for(int kk = 1; kk <= k; kk++){
        Rec_sw maxrec = local_search(w, kk, record);
        if(maxrec.occ < minrec.occ){
            minrec = maxrec;
            record = minrec.occ;
            break;
        }
    }
    return minrec;
}

// generate a random word with given length
static inline Word random_word(int n){
    u64 bits = 0;
    Runtab swruns;
    for(int i = 0; i < n - 1; i++){
        bits <<= 1;
        bits += (rand() >> 5) & 1;
    }
    return build_word(bits, n, swruns);
}

// metaheuristic, mixing iteratively
void mixed_descent(int n, int maxk, u64 maxiter){
    time_t mytime = time(NULL);
    srand(mytime);
    printf("Starting with n = %d, maxk = %d, maxiter = %lu, %s",
           n, maxk, maxiter, ctime(&mytime));
    // initial record, dummy hint
    Word w = random_word(n);
    Rec_sw currec = local_search_full(w, maxk, maxfreq_subword_fast(w));
    Rec_sw bestrec = currec;
    Runtab swruns;
    mytime = time(NULL);
    printf("%s", ctime(&mytime));
    print_record(&bestrec);
    // random flips
    int flipcnt = maxk + 2;
    u64 itercnt = 0;
    bool flag = true;
    while(true){
        itercnt++;
        if(flag){
            // stagnant, flip more bits
            if(itercnt >= maxiter){
                itercnt = 0;
                flipcnt++;
                printf("Current flipcnt: %d\n", flipcnt);
                if(flipcnt * 3 > n) return;
            }
        } else { // we move, so restart more conservatively
            flipcnt = maxk + 2;
        }
        // flip random bits
        u64 bits = bestrec.word.bits;
        for(int i = 0; i < n - 1; i++){
            if(rand() / (float) RAND_MAX * (n - 1) < flipcnt) bits ^= 1ULL << i;
        }
        currec = local_search_full(build_word(bits, n, swruns), 
                                   maxk, bestrec.occ);
        if(currec.occ < bestrec.occ){
            bestrec = currec;
            mytime = time(NULL);
            printf("%s", ctime(&mytime));
            print_record(&bestrec);
            flag = false;
        } else {
            flag = true;
        }
    }
}