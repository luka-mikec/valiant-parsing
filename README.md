# Valiant's parsing algorithm

... is, to my knowledge, asymptotically the fastest known parsing algorithm.

Algorithm goes as follows:
 1. Convert grammar to Chomsky normal form.
 2. Build a Cocke–Younger–Kasami-like parsing matrix. Define semiring operations on sets of nonterminals (*clauses*).
 3. Unlike *CYK*, Compute closure using Valiant's divide-and-conquer algorithm.
 4. Use binary (boolean) matrix  multiplication. *Not implemented yet.*


### Things inside
 - *matrix* and *matrix_view* generic structures. Both expect certain operations to be defined on used element type. Inside: matrix.h
 - *grammar* structure, representing a context-free grammar. Can convert itself to CNF. Doesn't support epsilons/blanks. Can be loaded from a file (see **sample_input.txt**), Inside: grammar.h, grammar.cpp.
 - bunch of functions implementing Valiant's algorithm. See main() for usage.

### Compilation

```sh
$ git clone https://github.com/luka-mikec/valiant-parsing.git
$ cd valiant-parsing
$ g++ *.cpp -std=c++11 -o valiant
```

### Running
```sh
$ ./valiant
```

