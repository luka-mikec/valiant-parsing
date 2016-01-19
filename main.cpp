#include <iostream>
#include <valarray>
#include <functional>
#include <list>
#include <set>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <string>

#include "matrix.h"
#include "grammar.h"

using namespace std;



map<symbol, set<pair<symbol, symbol>>> inside_productions; // A -> AA
map<symbol, set<symbol>> outside_productions;  // A -> a
set<symbol> terminals, nonterminals;
map<symbol, int> nt_idxs;
matrix<clause> backtrack; // if P -> AB, then P \in [A][B]

template<>
matrix<clause> matrix_view<clause>::operator* (const matrix_view<clause>& o)
{
  int h = nonterminals.size();
  matrix<bool> *a = new matrix<bool>[h];
  matrix<bool> *b = new matrix<bool>[h];

  auto sym = nonterminals.begin();
  for (int i = 0; i < h; ++i, ++sym)
  {
    a[i] = matrix<bool>(rows(), cols());
    b[i] = matrix<bool>(o.rows(), o.cols());

    for (int j = 0; j < rows(); ++j)
      for (int k = 0; k < cols(); ++k)
        a[i][j][k] = (*this)(j, k).disjuncts.find(*sym) != (*this)(j, k).disjuncts.end();

    for (int j = 0; j < o.rows(); ++j)
      for (int k = 0; k < o.cols(); ++k)
        b[i][j][k] = o(j, k).disjuncts.find(*sym) != o(j, k).disjuncts.end();

    //cout << i << endl << a[i] << endl << b[i] << endl;
  }

  matrix<bool> *c = new matrix<bool>[h * h];

  for (int i = 0; i < h*h; ++i)
  {
    c[i] = matrix_view<bool>(a[i / h]) * b[i % h];
    // cout << i << endl << c[i];
  }

  matrix<clause> res(rows(), o.cols());
  for (int j = 0; j < rows(); ++j)
    for (int k = 0; k < o.cols(); ++k)
      for (int i = 0; i < h * h; ++i)
        if (c[i][j][k])
          res[j][k] += backtrack[i / h][i % h];

  return res;
  //return mul_op(*this, o);
}


clause producible_by(symbol a, symbol b)
{
  clause res;
  auto prod_pair = make_pair(a, b);

  for (auto& p : inside_productions)
    for (auto &rhs : p.second)
      if (rhs == prod_pair)
        res.disjuncts.insert(p.first);

  return res;
}

clause producible_by(symbol a)
{
  clause res;

  for (auto& p : outside_productions)
    for (auto &rhs : p.second)
      if (rhs == a)
        res.disjuncts.insert(p.first);

  return res;
}

clause operator*(symbol a, symbol b)
{
  return producible_by(a, b);
}


clause operator*(const clause &a, const clause &b)
{
  clause c;
  for (auto& aa : a.disjuncts)
    for (auto& bb : b.disjuncts)
    {
      clause r = aa * bb;
      c.disjuncts.insert(r.disjuncts.begin(), r.disjuncts.end());
    }
  return c;
}


typedef matrix<clause> valiant;
typedef matrix_view<clause> valiant_view;


valiant setup_parse_mat(string item)
{
  int min_size = item.size() + 1;
  int actual_size = 1;
  while (actual_size < min_size)
    actual_size *= 2;

  valiant _m(actual_size, actual_size);
  valiant_view whole(_m);
  for (int i = 0; i < item.size(); ++i)
    whole(i, /*actual_size - min_size +*/ i + 1 ) = producible_by(symbol(item.substr(i, 1), symbol::types::terminal));


  return _m;
}


void closure(valiant_view &p);

// pretpostavlja da su odgovarajuće podmatrice inplace zatvorene
void lemma(valiant_view &b, int r)
{
  int n = b.cols(); // r = n / 2;
  if (n == 1)
    return;
  if (n == 2)
  {
    valiant c = b * b;
    valiant_view whole_c(c);
    b += whole_c;
    return;
  }
  if (2*r < n)
    throw runtime_error("lemma: r < n/2 :( is matrix 2^n?");

  valiant_view upper_left(b, 0, r, 0, r),
    bottom_right(b, n-r, n, n-r, n),
    upper_middle(b, 0, n-r, n-r, r),
    right_middle(b, n-r, r, r, n),
    top_right(b, 0, n - r, r, n);

  valiant c = upper_middle * right_middle;
  //cout << "?" << upper_middle << endl << "?" << right_middle << endl << "?" << c << endl;
  valiant messy_d = b.to_matrix();
  valiant_view top_right_messy_d(messy_d, 0, n - r, r, n);
  top_right_messy_d += c;

  valiant d = messy_d.without_range(n-r, r, n-r, r);
  valiant_view whole_d(d);
  closure(whole_d);

  // todo: d+

  top_right += valiant_view(d, 0, n-r, n-r, 2*(n-r));

  //cout << "$" << c;

}

void p3(valiant_view &p);
void p4(valiant_view &p);

void p2(valiant_view &p)
{
  int n = p.rows();
  if (n < 4)
    return;

  valiant_view middle(p, n / 4, 3 * n / 4, n / 4, 3 * n / 4);
  p2(middle);

  valiant_view upper_left(p, 0, 3*n/4, 0, 3*n/4),
    bottom_right(p, n/4, n, n/4, n);
  p3(upper_left);
  //cout  << "$"<< bottom_right;
  p3(bottom_right);

  p4(p);
}

void p3(valiant_view &p)
{
  int n = p.rows();
  if (n < 3)
    return;

  lemma(p, 2*n/3);
}

void p4(valiant_view &p)
{
  int n = p.rows();
  if (n < 4)
    return;

  lemma(p, 3*n/4);
}

void closure(valiant_view &p)
{
  int n = p.rows();
  if (n == 1)
    return;
  if (n % 2 == 1)
    throw runtime_error("closure of non 2^n matrix");

  valiant_view upper_left(p, 0, n/2, 0, n/2),
    bottom_right(p, n/2, n, n/2, n);
  closure(upper_left);
  closure(bottom_right);

  p2(p);
}




void use_grammar(const grammar &h)
{
  grammar g = h;
  g.to_cfg();

  terminals.insert(g.terminals.begin(), g.terminals.end());
  nonterminals.insert(g.nonterminals.begin(), g.nonterminals.end());

  int ns = nonterminals.size();
  auto n_itr = nonterminals.begin();
  for (int i = 0; i < ns; ++i, ++n_itr)
    nt_idxs[*n_itr] = i;

  backtrack.init(ns, ns);

  for (auto& prod : g.prods)
  {
    if (prod.size() == 1)
      outside_productions[prod.left].insert(prod.right[0]);
    else
    {
      inside_productions[prod.left].insert(make_pair(prod.right[0], prod.right[1]));
      backtrack[nt_idxs[prod.right[0]]][nt_idxs[prod.right[1]]].set_sym(prod.left);
    }
  }


}

int main()
{
  ifstream ist("sample_input.txt");

  grammar g;
  g.read_from_stream(ist);
  use_grammar(g);

  /*valiant test = setup_parse_mat("ab");
  valiant_view base(test);
  valiant res = base * base;
  cout << res << endl << res;

  return 0;*/

  std::function<set<string>(string)> words_gen;

  for (int len = 1; len < 5; ++len)
  {
    words_gen = [&g, &words_gen, len](string prefix)->set<string> {
      set<string> res;
      for (int i = 0; i < g.terminals.size(); ++i)
      {
        string w = prefix + (string)g.terminals[i];
        if (w.size() < len)
        {
          set<string> more = words_gen(w);
          res.insert(more.begin(), more.end());
        }
        else
          res.insert(w);
      }
      return res;
    };
    set<string> exprs = words_gen("");
    for (string expr : exprs)
    {
      valiant test = setup_parse_mat(expr); // setup_parse_mat("((p+q)*((q+p)+r))");
      valiant_view base(test);

      //cout << base << endl;

      closure(base);

      if (base(0, len).disjuncts.find(symbol("S")) != base(0, len).disjuncts.end())
        cout << expr << "\t" << base(0, len) << endl;
    }
  }



  return 0;
}