#include "grammar.h"


symbol::symbol() : identifier("NOT INITIALIZED."), type(types::nonterminal) {}

symbol::symbol(string _identifier, types _type)
  : identifier(_identifier), type(_type) {}

symbol::symbol(const char* _identifier, types _type)
  : identifier(_identifier), type(_type) {}

bool symbol::operator==(const symbol& other) const
{
  return type == other.type && identifier == other.identifier;
}

symbol::operator string() const
{
  return identifier;
}
symbol::operator char() const
{
  return identifier[0];
}

bool symbol::operator <(const symbol & other) const
{
  return identifier < other.identifier;
}


ostream& operator<<(ostream& ostr, const symbol & n)
{
  ostr << string(n);
  return ostr;
}

ostream& operator<<(ostream& ostr, const clause& c)
{
  ostr << "[";
  //cout << c.disjuncts.size() << endl;
  for (auto& nont: c.disjuncts)
    ostr << nont << " ";
  ostr << "]";
  return ostr;
}

clause operator+(const clause &a, const clause &b)
{
  clause c = a;
  c.disjuncts.insert(b.disjuncts.begin(), b.disjuncts.end());
  return c;
}

clause operator+=(clause &a, const clause &b)
{
  a.disjuncts.insert(b.disjuncts.begin(), b.disjuncts.end());
  return a;
}


production::production(symbol l, valarray<symbol> r)
  : left(l), right(r)
{
  if (left.type == symbol::types::terminal)
    throw runtime_error("lhs must be nonterminal (" + left.identifier + ")");
  if (right.size() == 0)
    throw runtime_error("rhs can't be epsilon (" + left.identifier + ")");
}

symbol& production::operator[](int i)
{
  return right[i];
}

const symbol& production::operator[](int i) const
{
  return right[i];
}

int production::size() const
{
  return right.size();
}

bool production::operator==(const production& other) const
{
  return left == other.left && equal(begin(right), end(right), begin(other.right));
}



void grammar::add_terminal(string s)
{
  terminals.push_back(symbol(s, types::terminal));
}

void grammar::add_nonterminal(string s)
{
  nonterminals.push_back(symbol(s, types::nonterminal));
}

void grammar::add_production(symbol lhs, valarray<symbol> rhs)
{
  prods.push_back(production(lhs, rhs));
}

symbol::types grammar::which_is_it(string s) {
  for (auto& s:nonterminals) {if (s.type == types::nonterminal) return types::nonterminal;}
  for (auto& s:terminals) {if (s.type == types::terminal) return types::terminal;}
  throw runtime_error("unknown symbol '" + s + "'");
}

vector<production> grammar::productions_of(symbol s)
{
  vector<production> res;
  for (auto&prod: prods)
  {
    if (prod.left == s)
    {
      res.push_back(prod);
    }
  }
  return res;
}

symbol grammar::bake_nonterminal(string name)
{
  string s = name == "" ? "__" + to_string(baking_seed++) : name;
  add_nonterminal(s);
  return nonterminals.back();
}

// returns {G | s -*-> G} \ {s}
set<symbol> grammar::unit_closure(symbol sym)
{
  set<symbol> visited;
  list<symbol> qu = {sym};
  set<symbol> res;

  while (!qu.empty())
  {
    auto s = qu.front();
    visited.insert(s);
    res.insert(s);

    auto producibles = productions_of(s);

    for (production &prod : producibles)
      if (prod.size() == 1 && prod[0].type == types::nonterminal)
      if (visited.find(prod[0]) == visited.end())
        qu.push_back(prod[0]);

    qu.pop_front();
  }

  res.erase(sym);
  return res;
}

void grammar::unit_elimination()
{
  vector<set<symbol>> partition; // induced by unit_closure, also not really a partition :D
  vector<production> new_productions;

  for (auto s : nonterminals)
    partition.push_back(unit_closure(s));

  // exchange and removal step
  for (int k = 0; k < nonterminals.size(); ++k)
  {
    symbol nt = nonterminals[k];
    const set<symbol>& ps = partition[k];

    for (int i = 0; i < prods.size(); ++i)
    {
      production &prod = prods[i];
      if (ps.find(prod.left) != ps.end())
        new_productions.push_back(production(nt, prod.right));
      if (prod.right.size() == 1 && prod.right[0].type == types::nonterminal)
      {
        prods.erase(prods.begin() + i--);
      }
    }
  }
  for (auto& prod : new_productions)
    prods.push_back(prod);
}

void grammar::to_cfg()
{
  set<symbol> terminals_set;

  // deleting X -> x
  for (auto& term : terminals)
  {
    symbol box = bake_nonterminal("_" + term.identifier);

    for (auto& prod : prods)
      for (int i = 0; i < prod.size(); ++i)
        if (prod[i].identifier == term.identifier)
          prod[i] = box;

    add_production(box, {term});
    terminals_set.insert(term);
  }

  // deleting X -> Y
  unit_elimination();

  // deleting X -> QYZ
  for (int i = 0; i < prods.size(); ++i)
  {
    production& prod = prods[i];
    if (prod.size() >= 3)
    {
      vector<production> new_prods;
      symbol parent_nonterminal = prod.left;
      if (prod.size() % 2 == 0)
      {
        for (int i = 0; i < prod.size() / 2; ++i)
        {
          if (i < prod.size() / 2 - 1)
          {
            symbol x0 = bake_nonterminal("_" + prod[2 * i].identifier + prod[2 * i + 1].identifier);
            symbol x1 = bake_nonterminal("_" + prod[2 * i + 2].identifier + prod[2 * i + 3].identifier+"...");

            new_prods.push_back(production(parent_nonterminal, {x0, x1} ));
            new_prods.push_back(production(x0, {prod[2 * i], prod[2 * i + 1]} ));
            parent_nonterminal = x1;
          }
          else
          {
            new_prods.push_back(production(parent_nonterminal, {prod[2 * i], prod[2 * i + 1]}));
          }
        }
      }
      else
      {
        for (int i = 0; i < prod.size() / 2; ++i)
        {
          if (i < prod.size() / 2 - 1)
          {
            symbol x0 = bake_nonterminal("_" + prod[2 * i].identifier + prod[2 * i + 1].identifier);
            symbol x1 = bake_nonterminal("_" + prod[2 * i + 2].identifier + prod[2 * i + 3].identifier+"...");

            new_prods.push_back(production(parent_nonterminal, {x0, x1} ));
            new_prods.push_back(production(x0, {prod[2 * i], prod[2 * i + 1]} ));
            parent_nonterminal = x1;
          }
          else
          {
            symbol x0 = bake_nonterminal("_" + prod[2 * i].identifier + prod[2 * i + 1].identifier);

            new_prods.push_back(production(parent_nonterminal, {x0, prod[2 * i + 2]} ));
            new_prods.push_back(production(x0, {prod[2 * i], prod[2 * i + 1]} ));
          }
        }
      }

      for (auto& p : new_prods)
        prods.push_back(p);

      prods.erase(prods.begin() + i--);
    }
  }
}

void grammar::read_from_stream(istream& stream)
{
  string str;

  bool exp_term = true;
  while (stream >> str)
  {
    if (exp_term)
    if (str == "#")
      exp_term = false;
    else
      add_terminal(str);
    else
    if (str == "#")
      break;
    else
      add_nonterminal(str);
  }

  int state = 0; // 0 wait lhs, 1 expect ->, 2 rhs, 3 expect | or lhs
  string lhs_name;
  while (stream >> str)
  {
    switch (state)
    {
      case 0:
        lhs_name = str;
        ++state;
        break;
      case 1:
        if (stream && (str != "->")) throw runtime_error("parse err, -> expected, got " + str);
        ++state;
        break;
      case 2:
      {
        valarray<symbol> syms(str.size());
        for (int i = 0; i < str.size(); ++i)
        {
          string sym = str.substr(i, 1);
          types t = which_is_it(sym);
          syms[i] = symbol(sym, t);
        }
        add_production(symbol(lhs_name), syms);
        ++state;
      }
        break;
      case 3:
        if (stream && (str != "|")){
          lhs_name = str;
          state = 1;
        }
        else if (str == "|"){
          state = 2;
        }
        break;
    }
  }
  //to_cfg();
}

ostream& operator<<(ostream& ostr, const grammar& g)
{
  ostr << "terminals: ";
  for (auto& s : g.terminals)
    ostr << s << " ";
  ostr << "\nnonterminals: ";
  for (auto& s : g.nonterminals)
    ostr << s << " ";
  ostr << "\nproductions: ";
  for (auto& s : g.prods)
  {
    ostr << "\n  " << s.left << " -> ";
    for (int i = 0; i < s.size(); ++i)
      ostr << s[i] << " ";
  }
  return ostr;
}
