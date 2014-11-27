lda-bump-cpp
============

Latent Dirichlet allocation (LDA) with bumping variational inference.

Implements three versions of LDA
* Coordinate ascent (mean-field)
* Stochastic variational inference (https://github.com/Blei-Lab/onlineldavb)
* Bumping variational inference [1]

[1] Alp Kucukelbir and David M Blei. _Profile predictive inference._ arXiv
preprint arXiv:1411.0292, 2014.


Requirements
------------

`lda-bump-cpp` is written in C++11. It requires a modern compiler. It also
depends on Eigen 3, Boost, and CMake. It uses docopt (provided).

* Eigen 3: http://eigen.tuxfamily.org/
* Boost: http://www.boost.org/
* CMake: http://www.cmake.org/

Refer to platform-specific instructions for installation.
(I recommend homebrew on Mac OS X.)


Instructions to Build and Run
-----------------------------

The driver (main) program runs all three algorithms.

`cmake .`
`make driver`

A toy dataset of arXiv abstracts is provided.

Example
```
./driver --topics=5
         --vocabulary=data/arxiv-vocab.dat
         --datatr=data/arxiv-train-1k.dat
         --datatest=data/arxiv-test-1k.dat
```

For more help, run
```
./driver -h
```


Vocabulary Data Format
----------------------

A text file with each word (`[term_1]` through `[term_N]`) on a separate line.


Corpus Data Format
------------------

A text file where each line is of the form (the LDA-C format):

`[M] [term_1]:[count] [term_2]:[count] ...  [term_N]:[count]`

where `[M]` is the number of unique terms in the document, and the
`[count]` associated with each term is how many times that term appeared
in the document.


Visualizing the Output
----------------------

A python script visualizes the topics.
(Modified from https://github.com/Blei-Lab/onlineldavb)

```
./printtopics.py data/arxiv-vocab.dat
                 results/Thu_Nov_27_10-45-09_2014/lambda_coord_ascent.dat

./printtopics.py data/arxiv-vocab.dat
                 results/Thu_Nov_27_10-45-09_2014/lambda_svi.dat

./printtopics.py data/arxiv-vocab.dat
                 results/Thu_Nov_27_10-45-09_2014/lambda_bumping.dat
```

























