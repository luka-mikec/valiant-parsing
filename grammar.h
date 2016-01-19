#ifndef VALIANT_PARSING_GRAMMAR_H
#define VALIANT_PARSING_GRAMMAR_H

#include <algorithm>
#include <iostream>
#include <list>
#include <set>
#include <stdexcept>
#include <string>
#include <valarray>
#include <vector>

using std::begin;
using std::end;
using std::equal;
using std::istream;
using std::list;
using std::ostream;
using std::runtime_error;
using std::set;
using std::string;
using std::to_string;
using std::valarray;
using std::vector;

struct symbol
{
  string identifier;
  enum class types {terminal, nonterminal} type;

  symbol();
  symbol(string _identifier, types _type = types::nonterminal);
  symbol(const char* _identifier, types _type = types::nonterminal);

  bool operator==(const symbol& other) const;

  operator string() const;

  operator char() const;

  bool operator <(const symbol & other) const;
};

ostream& operator<<(ostream& ostr, const symbol & n);


// union of symbols
struct clause
{
  clause() {}
  clause(const symbol & nont) { disjuncts = {nont};}
  set<symbol> disjuncts;

  bool is_sym_set(symbol s) { return disjuncts.find(s) != disjuncts.end(); }
  void set_sym(symbol s) {  disjuncts.insert(s); }
};

ostream& operator<<(ostream& ostr, const clause& c);
clause operator+(const clause &a, const clause &b);
clause operator+=(clause &a, const clause &b);


struct production
{
  symbol left;
  valarray<symbol> right;
  production(symbol l, valarray<symbol> r);
  symbol& operator[](int i);
  const symbol& operator[](int i) const;
  int size() const;

  bool operator==(const production& other) const;
};


struct grammar
{

  typedef symbol::types types;

  vector<symbol> nonterminals;
  vector<symbol> terminals;
  vector<production> prods;

  void add_terminal(string s);
  void add_nonterminal(string s);
  void add_production(symbol lhs, valarray<symbol> rhs);

  types which_is_it(string s);

  vector<production> productions_of(symbol s);

  int baking_seed = 0;
  symbol bake_nonterminal(string name = "");

  // returns {G | s -*-> G} \ {s}
  set<symbol> unit_closure(symbol sym);

  void unit_elimination();

  void to_cfg();

  void read_from_stream(istream& stream);
};

ostream& operator<<(ostream& ostr, const grammar& g);


#endif //VALIANT_PARSING_GRAMMAR_H
