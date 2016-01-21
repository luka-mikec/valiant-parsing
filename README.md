# Valiant's parsing algorithm

... is the asymptotically fastest known parsing algorithm.

Algorithm goes as follows:
 1. Convert grammar to Chomsky normal form.
 2. Build a Cocke–Younger–Kasami-like (*CYK*) parsing matrix. Define semiring operations on sets of nonterminals (*clauses*).
 3. Unlike *CYK*, compute closure using Valiant's divide-and-conquer algorithm.
 4. Use binary (boolean) matrix  multiplication.


## Things inside
 - *matrix* and *matrix_view* generic structures. Both expect certain operations to be defined on used element type. Inside: matrix.h.
 - *grammar* structure, representing a context-free grammar. Can convert itself to CNF. Doesn't support epsilons/blanks. Can be loaded from a file (see **sample_input.txt**), Inside: grammar.h, grammar.cpp.
 - bunch of functions implementing Valiant's algorithm. See main() for usage. Inside: main.cpp.


## Compilation
On Unix-based systems:
```sh
$ git clone https://github.com/luka-mikec/valiant-parsing.git
$ cd valiant-parsing
$ g++ *.cpp -std=c++11 -o valiant
```

On Windows, use any modern C++ compiler with "-std=c++11" or equivalent compiler flag.

## Running
### Examples
Evaluating (0+1) and (1*(0+0)) in BA:
```sh
$. echo "(0+1) (1*(0+0))" | ./valiant -g boolean_algebra.txt
```
Search for closed algebraic expressions of BA evaluating to 1, up to length 5, by iterating all words over alphabet:
```sh
$ ./valiant -g boolean_algebra.txt -language 5 -skip_interactive
```

### Full syntax
```sh
$ ./valiant [-g <grammar_file_address>] [-language <inclusive_maximal_length>] [-skip_interactive] [-bmm] [-show_grammar] [-show_table]
```
Arguments:
 1. *-g <grammar_file_address>* specifies grammar. See below for details on this file's format. Default value: boolean_algebra.txt
 2. *-language <inclusive_maximal_length>* will enumerate all words over alphabet of up to given length (inclusive).
 3. *-skip_interactive* will skip waiting for STDIN. Without using this parameter, waiting can still be skipped later by pressing Ctrl + D while in interactive mode.
 4. *-bmm* will use binary matrix multiplication. This won't improve the algorithm in any way, but won't make it asymptotically worse either. It **will** make the calculation **slower** though. Implemented only because Valiant used it in the article.
 5. *-show_grammar* will display loaded grammar before and after conversion to Chomsky normal form.
 6. *-show_table* will display final (CYK-based) parsing matrix.

Only the first argument is mandatory.

### Usage
Strings to test can be entered at runtime. Use *Ctrl + D* to skip interactive mode. Use *Ctrl + C* to exit.
Otherwise, write all the strings in a file, separated by whitespace.

Output is a list of 0's and 1's, *i*th number representing whether *i*th input string belongs to language.
Output is given one character per line.

## Grammar file
### Examples
See boolean_algebra.txt for a larger example. New lines are optional.


Grammar (V, T, P, S) = ({S, A, B}, {a, b}, {S -> AB, A -> AB | aB, B .> b}, S) is given by
```tex
a b
#
S A B
#
S -> AB
A -> AB | aB
B -> b
```

### Full syntax
Grammar file consists of:
 1. Whitespace-separated list of terminal symbols. Terminal symbol is any non-whitespace ASCII byte, followed by
 2. the symbol #, a whitespace and then:
 3. Whitespace-separated list of nonterminal symbols. Nonterminal symbol is any non-whitespace ASCII byte that wasn't already used. The list should probably include *S*, the root nonterminal symbol. Then:
 4. the symbol #, a whitespace and then:
 5. Whitespace-separated list of productions, where production is:
     1. A left-hand-side, which is any previosly declared nonterminal symbol.
     2. A right-hand-side, which is one of the following:
         1. A chain of symbols, each of which is either a terminal or a nonterminal.
         2. A right-hand-side, followed by whitespace, symbol |, whitespace, then another right-hand-side.


