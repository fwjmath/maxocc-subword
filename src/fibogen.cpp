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

#include "fibogen.hpp"

static void fill_runs2(int remains, Fibo_state* state){
    u64 lastbit = state->curbits & 1; // in fact, lastbit contains two bits...
    lastbit += lastbit << 1;
    while(remains >= 2){
        lastbit = 3 - lastbit; // take the opposite
        state->curbits <<= 2;
        state->curbits += lastbit;
        remains -= 2;
    }
    if(remains == 1){
        lastbit = 3 - lastbit;
        state->curbits <<= 1;
        state->curbits += (lastbit & 1);
    }
    return;
}

// initilization, returns whether it is successful
bool fibogen_init(int n, Fibo_state* state){
    if(n < 2 || n > 64) return false; // only doing with 64 bits
    state->bitcnt = n;
    // we start with a word fill with runs of length 2
    state->curbits = 0;
    fill_runs2(state->bitcnt - 2, state);
    return true;
}

// return the next word in a pointer, and returns whether there is a next word
// we do it in the order of first doing run length 2, then run length 1
bool fibogen_next(u64* bits, Fibo_state* state){
    *bits = state->curbits;
    // determine the current run length
    int lastrun = 2 - ((state->curbits ^ (state->curbits >> 1)) & 1);
    // first case: run of size 2, we can replace it by two runs of size 1
    if(lastrun == 2){
        state->curbits ^= 1; // change the last bit
        return true;
    }
    // second case: run of size 1
    // we frist remove all runs of size 1, then fill with runs of size 2
    int remains = 0;
    while(((state->curbits ^ (state->curbits >> 1)) & 1)){
        state->curbits >>= 1;
        remains++;
    }
    // we stop when it is filled with 1
    // note that we started with first bit 0, so we will only reach bitcnt - 1
    if(remains == state->bitcnt - 1) return false;
    // now we are sur to have a run of length 2, and we reduce its length
    state->curbits >>= 1;
    remains++;
    // fill with runs of length 2 till the end
    fill_runs2(remains, state);
    return true;
}