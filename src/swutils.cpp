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

#include "swutils.hpp"

// build a word according to a 0-1 string
Word build_word_str(const char* str, Runtab wruns){
    int n = strlen(str);
    u64 bits = 0;
    // assuming 0-1 string
    for(int i = 0; i < n; i++){
        bits <<= 1;
        bits += str[i] - '0';
    }
    return build_word(bits, n, wruns);
}

// print record
void print_record(Rec_sw* minrec){
    print_word_bin(minrec->word);
    for(auto sw : minrec->subwords){
        printf("Subword: ");
        print_word_bin(sw);
    }
    return;
}

// exhaustive search for minimal subword entropy, using a hint
// using subword with large number of occurrences from the last word as a hint
void hinted_search(int n, u64 hint){
    Rec_occ minrec = min_maxfreq_subword_hinted(n, hint);
    printf("%d bits, hint %lu, found %lu\n", n, hint, minrec.occ);
    for(auto rec : minrec.recs){
        print_record(&rec);
    }
    return;
}

// exhaustive search for minimal subword entropy, using a hint, branch and bound
void pruned_search(int n, u64 hint){
    Rec_occ minrec = min_maxfreq_subword_pruned(n, hint);
    printf("%d bits, hint %lu, found %lu\n", n, hint, minrec.occ);
    for(auto rec : minrec.recs){
        print_record(&rec);
    }
    return;
}

// build a histogram for subword occurrences
void histo_subword(int n){
    Histogram histo = maxfreq_subword_histo(n);
    printf("Maximal subword occurrences histogram for %d bits\n{\n", n);
    for(const auto& [freq, cnt] : histo){
        printf("%lu: %lu\n", freq, cnt);
    }
    printf("}\n");
    return;
}

// compute the most frequent subwords of a given word
void compute_maxfreq_subword(char* wstr){
    int n = strlen(wstr);
    Runtab wruns;
    Rec_sw minrec = maxfreq_subword_hinted(build_word_str(wstr, wruns), 1ul << n);
    printf("Word %s, maxocc %lu\n", wstr, minrec.occ);
    print_record(&minrec);
    return;
}

// adds a letter somewhere in a hinted word (previous record), using incomplete
// computation
void insert_heuristic(char* wstr){
    int n = strlen(wstr);
    Runtab wruns;
    Word oldw = build_word_str(wstr, wruns);
    u64 wbits = oldw.bits;
    u64 recw = wbits;
    u64 recocc = maxfreq_subword_fast(oldw) << 1;
    // insert a new bit at each possible way
    for(int i = 0; i < n; i++){
        for(int bit = 0; bit < 2; bit++){
            // insert a new bit at position i
            u64 newbits = wbits >> i;
            newbits <<= 1;
            newbits += bit;
            newbits <<= i;
            newbits += wbits & ((1ul << i) - 1);
            Word w = build_word(newbits, n + 1, wruns);
            u64 swocc = maxfreq_subword_fast(w);
            if(swocc < recocc){
                recocc = swocc;
                recw = newbits;
            }
        }
    }
    oldw = build_word(recw, n + 1, wruns);
    print_word_bin(oldw);
    printf("Maxocc (fast): %lu\n", recocc);
    return;
}
