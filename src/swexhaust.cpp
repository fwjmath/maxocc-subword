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

#include "swexhaust.hpp"
#include <time.h>

// compute max frequence subword with given length, for histogram, no speed up
static inline u64 maxfreq_subword_len(Word w, int k){
    u64 maxocc = 0;
    Runtab wruns;
    Word sw = build_word(w.bits & 1, k, wruns);
    do {
        u64 occ = subword_cnt(w, sw);
        if(occ >= maxocc){
            maxocc = occ;
        }
    } while(increment_word_2(&sw));
    return maxocc;
}

// compute max freq subword of all lenghts, for histogram
u64 maxfreq_subword(Word w){
    u64 maxocc = 0;
    // check different lengths
    for(int curk = 2; curk < w.len - 1; curk++){
        u64 occ = maxfreq_subword_len(w, curk);
        if(occ >= maxocc){
            maxocc = occ;
        }
    }
    return maxocc;
}

// compute max freq subword up to some length, for starting metaheuristics
u64 maxfreq_subword_fast(Word w){
    u64 maxocc = 0;
    // check different lengths
    for(int curk = w.len / 4; curk < w.len / 2; curk++){
        u64 occ = maxfreq_subword_len(w, curk);
        if(occ >= maxocc){
            maxocc = occ;
        }
    }
    return maxocc;
}

// test whether it is the smallest in its equivalence class
static inline int sym_mult(u64 bits, int len){
    // clean up
    bits &= (1ULL << len) - 1;
    // reverse the bits
    u64 rev = 0;
    u64 orig = bits;
    for(int i = 0; i < len; i++){
        rev <<= 1;
        rev += orig & 1;
        orig >>= 1;
    }
    if(bits & 1){ // no chance for rev, as it starts with 1, better switch 0 and 1
        rev = ~rev;
        rev &= (1 << len) - 1;
    }
    if(bits < rev)
        return 2;
    else if(bits == rev)
        return 1;
    return 0; 
}

// returns histogram of max subword occurrences
Histogram maxfreq_subword_histo(int n){
    // construct the word
    Runtab wordruns;
    Word w = build_word(0, n, wordruns);
    // initialize the histogram
    Histogram histo = Histogram(); 
    do {
        int mult = sym_mult(w.bits, n);
        if(mult == 0) continue; // only test primitive ones
        u64 freq = maxfreq_subword(w);
        if(!histo.contains(freq)){
            histo[freq] = 0;
        }
        histo[freq] += mult;
    } while(increment_word(&w));
    return histo;
}

// compute max frequence subword with given length
// exhaustive, but stops once we find a subword breaking record
// as we will be taking maximum for a given word, when one subwordsuch is found,
// no need to test further as it will not improve the record
static void maxfreq_subword_len_hinted(Rec_sw* maxrec, int k, u64 record){
    // first check: are there enough subwords occurrences?
    // need to check if the record is a real one or just the max
    // TODO: can we improve this?
    if(record != (1ULL << maxrec->word.len) && binomial(maxrec->word.len, k) < record) return;
    // initialization
    maxrec->occ = 0;
    maxrec->subwords.clear();
    Runtab swruns;
    Word w = maxrec->word;
    Word sw = build_word(w.bits & 1, k, swruns);
    // the loop
    do {
        u64 occ = subword_cnt(w, sw);
        if(occ >= maxrec->occ){
            if(occ > maxrec->occ) maxrec->subwords.clear();
            maxrec->occ = occ;
            maxrec->subwords.push_back(sw);
            if(occ > record) break;
        }
    } while(increment_word_2(&sw));
    return;
}

static Rec_sw maxfreq_subword_hinted(Word w, u64 record, Word* lastsw){
    Rec_sw maxrec = {w, std::vector<Word>(), 1};
    Rec_sw maxrec_len = {w, std::vector<Word>(), 1};
    int lastsw_len = lastsw->len;
    u64 lastsw_bits = lastsw->bits;
    Runtab swruntab;
    // filter with heuristics
    // if one of the following constructed subword give something bigger than
    // the record, then we can stop
    u64 newbits[3];
    // first filter: replace the last bit
    newbits[0] = lastsw_bits & (~(u64)1) + (w.bits & 1);
    // second filter: add the last bit
    newbits[1] = (lastsw_bits << 1) + (w.bits & 1);
    // third filter: remains the same
    newbits[2] = lastsw_bits;
    for(int i = 0; i < 3; i++){
        int newlen = (i == 1 ? lastsw_len + 1 : lastsw_len); 
        Word newsw = build_word(newbits[i], newlen, swruntab);
        u64 filter_occ = subword_cnt(w, newsw);
        if(filter_occ > record){
            maxrec.subwords.push_back(newsw);
            maxrec.occ = filter_occ;
            return maxrec;
        }
    }
    // another filter: flip a bit
    for(int i = 1; i < lastsw_len - 1; i++){
        u64 modsw = lastsw_bits ^ (1ull << i);
        Word newsw = build_word(modsw, lastsw_len, swruntab);
        u64 filter_occ = subword_cnt(w, newsw);
        if(filter_occ > record){
            maxrec.subwords.push_back(newsw);
            maxrec.occ = filter_occ;
            return maxrec;
        }
    }
    // yet another filter: flip two bits
    for(int i = 1; i < lastsw_len - 2; i++){
        for(int j = i + 1; j < lastsw_len - 1; j++){
            u64 modsw = lastsw_bits ^ (1ull << i) ^ (1ull << j);
            Word newsw = build_word(modsw, lastsw_len, swruntab);
            u64 filter_occ = subword_cnt(w, newsw);
            if(filter_occ > record){
                maxrec.subwords.push_back(newsw);
                maxrec.occ = filter_occ;
                return maxrec;
            }
        }
    }
    // again another filter: words with run length only 1 and 2
    // needs more test to see if it leads to speedup for larger n
    // speeds up for n=37, about 15%, so promoted to regular usage
    Fibo_state fbst;
    if(lastsw_len >= 3 && fibogen_init(lastsw_len, &fbst)){
        u64 bits = 0;
        while(true){
            bool contd = fibogen_next(&bits, &fbst);
            Word newsw = build_word(bits, lastsw_len, swruntab);
            u64 filter_occ = subword_cnt(w, newsw);
            if(filter_occ > record){
                maxrec.subwords.push_back(newsw);
                maxrec.occ = filter_occ;
                // we update here because it may change a lot
                *lastsw = maxrec.subwords[0];
                return maxrec;
            }
            if(!contd) break;
        }
    }
    // check different lengths with most probable order
    int curk = lastsw_len;
    int curdev = 0;
    while(true){
        maxfreq_subword_len_hinted(&maxrec_len, curk, record);
        if(maxrec_len.occ >= maxrec.occ){
            if(maxrec_len.occ > maxrec.occ) maxrec.subwords.clear();
            maxrec.occ = maxrec_len.occ;
            maxrec.subwords.insert(maxrec.subwords.end(),
                                   maxrec_len.subwords.begin(),
                                   maxrec_len.subwords.end());
            if(maxrec.occ > record){
                *lastsw = maxrec.subwords[0];
                break;
            }
        }
        do {
            if(curdev < 0) 
                curdev = -curdev;
            else 
                curdev = -curdev - 1;
            curk = lastsw_len + curdev;
        } while (curk < 2);
        if(curk == w.len - 1) break;
    }
    return maxrec;
}

// compute most frequent subwords for a single given word
Rec_sw maxfreq_subword_single(Word w, u64 record){
    Word lastsw = {0, NULL, 1, 2};
    return maxfreq_subword_hinted(w, record, &lastsw);
}

// for metaheuristics
Rec_sw maxfreq_subword_hinted_fast(Word w, u64 record){
    Rec_sw maxrec = {w, std::vector<Word>(), 1};
    Rec_sw maxrec_len = {w, std::vector<Word>(), 1};
    // check different lengths with most probable order
    for(int curk = w.len / 4; curk < w.len / 2; curk++){
        maxfreq_subword_len_hinted(&maxrec_len, curk, record);
        if(maxrec_len.occ >= maxrec.occ){
            if(maxrec_len.occ > maxrec.occ) maxrec.subwords.clear();
            maxrec.occ = maxrec_len.occ;
            maxrec.subwords.insert(maxrec.subwords.end(),
                                   maxrec_len.subwords.begin(),
                                   maxrec_len.subwords.end());
            if(maxrec.occ > record){
                break;
            }
        }
    }
    return maxrec;
}

// test whether it is the smallest in its equivalence class
static inline bool is_primitive(u64 bits, int len){
    return sym_mult(bits, len) > 0;
}

static inline void update_minrec(Rec_occ* minrec, Rec_sw maxrec){
    if(minrec->occ >= maxrec.occ){
        if(minrec->occ > maxrec.occ){
            minrec->recs.clear();
        }
        minrec->occ = maxrec.occ;
        minrec->recs.push_back(maxrec);
    }
    return;
}

// exhaustive search with a hint
Rec_occ min_maxfreq_subword_hinted(int n, u64 record){
    // construct the word
    Runtab wruns;
    Word w = build_word(0, n, wruns);
    Word lastsw = {0, NULL, 1, 2};
    // initialize the records
    Rec_occ minrec;
    minrec.occ = record;
    minrec.recs = std::vector<Rec_sw>();
    do {
        if(!is_primitive(w.bits, n)) continue; // only test primitive ones
        update_minrec(&minrec, maxfreq_subword_hinted(w, record, &lastsw));
        record = minrec.occ;
    } while(increment_word(&w));
    return minrec;
}

// exhaustive search with a hint, parallel version
void* min_maxfreq_subword_hinted_parallel(void* info){
    // get information
    Thread_info tinfo = *((Thread_info*) info);
    int n = tinfo.n;
    int tid = tinfo.thread_id;
    u64 record = tinfo.record;
    // construct the word
    Runtab wruns;
    Word w = build_word(0, n, wruns);
    int segstart = n >> 1;
    // initialize the records
    tinfo.minrec->occ = record;
    tinfo.minrec->recs = std::vector<Rec_sw>();
    Word lastsw = {0, NULL, 1, 2};
    bool has_next = true; // there are still words to check
    do {
        if(((w.bits >> segstart) & (THREAD_COUNT - 1)) == tid
            && is_primitive(w.bits, n)){ // only test primitive ones
            update_minrec(tinfo.minrec, maxfreq_subword_hinted(w, record, &lastsw));
            record = tinfo.minrec->occ;
        }
    } while(increment_word(&w));
    // measure the time
    time_t mytime = time(NULL);
    printf("Thread %d finished at %s", tid, ctime(&mytime));
    return NULL;
}