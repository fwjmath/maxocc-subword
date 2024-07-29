# maxocc-subword

By Wenjie Fang, released under GPLv3.

This is the research code for my work on minimizing maximal occurrences of binary words of a given length. It is not exactly well-written, but it seems to work.

Given a binary word `w`, we say that `v` of length `k` is a **subword** of `w` if there are positions `p1` < `p2` < ... < `pk` such that reading `w` at positions from `p1` to `pk` gives the word `v`, and those positions form an **occurrence** of `v` in `w`. The **maximal subword occurrence** of `w` is given by the number of occurrences of subwords with the highest number of occurrences.

As it is easy to find words with large maximal subword occurrences, we try instead to minimize this quantity. It is thus the goal of the program here, which comes in several modes:

1. Exhaustive search on single thread (with possibly a hint)

2. Exhaustive search on multithread (with possibly a hint)

3. Metaheuristic search, to obtain a hint

4. Letter insertion in a given word, to obtain a hint (proposed by Jim Fill)

5. Histogram, to study the landscape of the metaheuristic search and some statistics of maximal subword occurrences.

## Algorithm

The core of the algorithm is a function that compute the number of occurrences of a given subword in a given word. In this algorithm, the unit of computation is "runs" instead of bits, and we perform a dichotomy to reduce the problem to smaller subproblems. When the subproblems become smaller than a threshold, we store the result in a cache to accelerate future computation.

Then, for the hinted exhaustive search, a hint (upper bound of minimal maxocc) is given to stop the search once we know that a certain word cannot produce better result. Some heuristics are also used to find likely frequent subwords, which can certify that the current word cannot produce results lower than the hint, thus can be ignored. Such heuristics greatly accelerate the computation.

To obtain a reasonable hint, we have a meta-heuristic search that combines iterative deepening exhaustive local search and stochastic jumps.

## Different branches

There are three branches of this project, each offering slightly different functionnality.

1. Branch `main` contains multithread functionnality using simple `pthread`, with precomputed lookup table.

2. Branch `tbb-mt` contains multithread functionnality using `pthread` and Thread Building Blocks (TBB), with dynamic lookup table as hashmap from TBB.

3. Branch `single` does not contain multithread functionnality, but is slightly more optimized than the multithread version on a single core.

## Compilation

On a machine with `gcc` and `make`, compile with

```
make
```

## Execution

Suppose that we want to compute words with 19 bits that minimizes the maximal subword occurrences. It suffices to run:

```
./swmain 19
```

If you want to accelerate it with a hint, just add the hint after it:

```
./swmain 19 1000
```

Note that if the hint is smaller than the real number, the program will return nothing, as every word would be pruned.

For multithread functionnality, it suffices to run:

```
./swmain 19 mt 1000
```

The last hint (here `1000`) is optional. The number of thread is fixed upon compilation, with the variable `THREADCNT`. The default value is 4, meaning that 4 threads are used at the same time. There is no workload balancing, so some thread may exit before others. Thus, this program does not fully use all the cores for a certain among of time at the end of the computation.

To obtain a reasonable hint, we may run meta-heursitic search with:

```
./swmain 19 meta 2 10
```

Here, the first parameter after "meta" is the size of the neighborhood that we will be performing exhaustive local search for the minimal element, and the second parameter is the number of stochastic jumps we will be performing when the local search fails to improve the result.

Another way to obtain a reasonable hint is to use Jim Fill's proposed heuristic that inserting a letter in a word reaching minimal maxocc usually gives a word with quite small maxocc, which can be used as a hint. For this, we may run:

```
./swmain 19 insert 011100100101110001
```

The word given as the last argument is of length one smaller, and it gives a hint for the length given in the first argument.

For histogram, we may run:

```
./swmain 19 histo
```

It produces in stdout the histogram of maximal subword occurrences of words with 19 bits, in the format of Python dictionary.

## Extra

We provide a Sagemath notebook `maxocc-periodic-gf.ipynb` that computes the generating function of occurrences of periodic subwords with a given period in periodic words with another given period.
