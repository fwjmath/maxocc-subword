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

#ifndef __FIBOGEN__
#define __FIBOGEN__

#include <stdint.h>

typedef uint64_t u64;

// intialize the generator of words containing only runs of length 1 or 2
// they are counted by the Fibonacci numbers, thus the name of this module
// returns whether it has been properly initialized
bool fibogen_init(int n);

// returns the next word in a pointer, and returns whether there is a next word
bool fibogen_next(u64* bits);

#endif