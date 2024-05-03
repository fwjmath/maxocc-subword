# maxocc-subword

By Wenjie Fang, released under GPLv3.

This is the research code for my work on minimizing maximal occurrences of binary words of a given length. It is not exactly well-written, but it seems to work.

Given a binary word `w`, we say that `v` of length `k` is a **subword** of `w` if there are positions `p1` < `p2` < ... < `pk` such that reading `w` at positions from `p1` to `pk` gives the word `v`, and those positions form an **occurrence** of `v` in `w`. The **maximal subword occurrence** of `w` is given by the number of occurrences of subwords with the highest number of occurrences.

As it is easy to find words with large maximal subword occurrences, we try instead to minimize this quantity. It is thus the goal of the program here, which comes in 3 modes:

1. Exhaustive search (with possibly a hint)

2. Metaheuristic search, to obtain a hint

3. Histogram, to study the landscape of the metaheuristic search and some statistics of maximal subword occurrences.

## Algorithm

The core of the algorithm is a function that compute the number of occurrences of a given subword in a given word. In this algorithm, the unit of computation is "runs" instead of bits, and we perform a dichotomy to reduce the problem to smaller subproblems. When the subproblems become smaller than a threshold, we store the result in a cache to accelerate future computation.

Then, for the hinted exhaustive search, a hint (upper bound of minimal maxocc) is given to stop the search once we know that a certain word cannot produce better result. Some heuristics are also used to find likely frequent subwords, which can certify that the current word cannot produce results lower than the hint, thus can be ignored. Such heuristics greatly accelerate the computation.

To obtain a reasonable hint, we have a meta-heuristic search that combines iterative deepening exhaustive local search and stochastic jumps.

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

To obtain a reasonable hint, we may run meta-heursitic search with:

```
./swmain 19 meta 2 10
```

Here, the first parameter after "meta" is the size of the neighborhood that we will be performing exhaustive local search for the minimal element, and the second parameter is the number of stochastic jumps we will be performing when the local search fails to improve the result.

For histogram, we may run:

```
./swmain 19 histo
```

It produces in stdout the histogram of maximal subword occurrences of words with 19 bits, in the format of Python dictionary.
