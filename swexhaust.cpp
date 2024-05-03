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

// compute max frequence subword with given length, for histogram, no speed up
static inline u64 maxfreq_subword_len(Word w, int k){
    u64 maxocc = 0;
    Word sw = build_word(w.bits & 1, k, true);
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
    Word w = build_word(0, n, false);
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
    maxrec->occ = 0;
    maxrec->subwords.clear();
    Word w = maxrec->word;
    Word sw = build_word(w.bits & 1, k, true);
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

Rec_sw maxfreq_subword_hinted(Word w, u64 record){
    static Word lastsw = {0, NULL, 1, 2};
    Rec_sw maxrec = {w, std::vector<Word>(), 1};
    Rec_sw maxrec_len = {w, std::vector<Word>(), 1};
    int lastsw_len = lastsw.len;
    u64 lastsw_bits = lastsw.bits;
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
        Word newsw = build_word(newbits[i], newlen, true);
        u64 filter_occ = subword_cnt(w, newsw);
        if(filter_occ > record){
            maxrec.subwords.push_back(newsw);
            maxrec.occ = filter_occ;
            return maxrec;
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
                lastsw = maxrec.subwords[0];
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
    Word w = build_word(0, n, false);
    // initialize the records
    Rec_occ minrec;
    minrec.occ = record;
    minrec.recs = std::vector<Rec_sw>();
    do {
        if(!is_primitive(w.bits, n)) continue; // only test primitive ones
        update_minrec(&minrec, maxfreq_subword_hinted(w, record));
        record = minrec.occ;
    } while(increment_word(&w));
    return minrec;
}

// extend the candidate word in pruned exhaustive search
static void min_maxfreq_prune(int n, Word* pw, Rec_occ* minrec){
    if(pw->len == n && is_primitive(pw->bits, n)){ // we have a full word
        update_minrec(minrec, maxfreq_subword_hinted(*pw, minrec->occ));
        return; 
    }
    // we don't have a full word, try to extend
    // add bits
    for(int bit = 0; bit < 2; bit++){
        add_bit(pw, bit);
        // printf("Add bit %d\n", bit);
        if(maxfreq_subword_hinted(*pw, minrec->occ).occ <= minrec->occ){
            min_maxfreq_prune(n, pw, minrec);
        }
        remove_bit(pw);
    }
    return;
}

// exhaustive search by constructing partial words and pruning those with more
// than the hinted number of subword occurrences.
Rec_occ min_maxfreq_subword_pruned(int n, u64 record){
    int localruns[MAXLEN];
    int len = n >> 1;
    Rec_occ minrec = {std::vector<Rec_sw>(), record};
    Word pw = build_word_ext(0, len + 1, localruns);
    do {
        min_maxfreq_prune(n, &pw, &minrec);
    } while(increment_word(&pw));
    return minrec;
}