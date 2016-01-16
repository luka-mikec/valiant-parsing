#include <iostream>
#include <valarray>
#include <functional>
#include <list>
#include <set>
#include <fstream>
#include <stdexcept>

#include <map>
#include <utility>


#include <string>

using namespace std;

struct symbol
{
  string identifier;
  enum class types {terminal, nonterminal} type;

  symbol(string _identifier, types _type = types::nonterminal)
    : identifier(_identifier), type(_type) {}
  symbol(const char* _identifier, types _type = types::nonterminal)
    : identifier(_identifier), type(_type) {}

  operator string() const {return identifier; }
  operator char() const {return identifier[0];}
  bool operator <(const symbol & other) const { return identifier < other.identifier; }
};



struct clause
{
  clause() {}
  clause(const symbol & nont) { disjuncts = {nont};}
  set<symbol> disjuncts;
};



ostream& operator<<(ostream& ostr, const symbol & n)
{
  ostr << string(n);
  return ostr;
}

ostream& operator<<(ostream& ostr, const clause& c)
{
  ostr << "[";
  for (auto& nont: c.disjuncts)
    ostr << nont << " ";
  ostr << "]";
  return ostr;
}

map<symbol, set<pair<symbol, symbol>>> inside_productions; // A -> AA
map<symbol, set<symbol>> outside_productions;  // A -> a



// expects A -> a, or A -> AA format, all one-letter
void read_from_stream(istream& stream)
{
  string str;

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
        if (isupper(str[0]))
          inside_productions[lhs_name].insert(make_pair(symbol(str.substr(0, 1)), symbol(str.substr(1, 1))));
        else
          outside_productions[lhs_name].insert(symbol(str, symbol::types::terminal));
        ++state;
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

template <typename T> struct matrix_view;

template <typename T>
struct matrix
{
  int rows() const
  {
    return cells.size();
  }

  int cols() const
  {
    return cells[0].size();
  }

  void init(int r, int s)
  {
    cells.resize(r, valarray<T>(s));
  }

  matrix() {}
  matrix(int r, int s) { init(r, s); }

  valarray<T>& operator[](int i)
  {
    return cells[i];
  }

  const valarray<T>& operator[](int i) const
  {
    return cells[i];
  }

  matrix without_range(int dr_beg, int dr_end, int dc_beg, int dc_end)
  {
    int newrows = rows() - (dr_end - dr_beg);
    int newcols = cols() - (dc_end - dc_beg);

    matrix new_cells;
    new_cells.cells.resize(newrows);
    
    for (int i = 0; i < newrows; ++i)
      new_cells[i] = valarray<T>(newcols);

    int citr = 0;
    for (int i = 0; i < dc_beg; ++i)
    {
      for (int j = 0; j < newrows; ++j)
        new_cells[j][citr] = cells[j + (j < dr_beg ? 0 : (dr_end - dr_beg))][i];
      ++ citr;
    }
    for (int i = dc_end; i < cols(); ++i)
    {
      for (int j = 0; j < newrows; ++j)
        new_cells[j][citr] = cells[j + (j < dr_beg ? 0 : (dr_end - dr_beg))][i];
      ++ citr;
    }

    return new_cells;
  }

  valarray<valarray<T>> cells;

  //operator matrix_view<T>() const;

};

typedef matrix<clause> valiant;

template <typename T>
class matrix_view
{
  matrix<T> &m;
  int r_beg, r_end, c_beg, c_end;

public:
  matrix_view(matrix<T>& m, int r_beg, int r_end, int c_beg, int c_end)
    : r_beg(r_beg), r_end(r_end), c_beg(c_beg), c_end(c_end), m(m)
  {

  }

  matrix_view(matrix_view& m, int r_beg, int r_end, int c_beg, int c_end)
    : r_beg(r_beg + m.r_beg), r_end(r_end + m.r_beg), c_beg(c_beg + m.c_beg), c_end(c_end + m.c_beg), m(m.m)
  {

  }

  matrix_view(matrix<T>& m)
    : r_beg(0), r_end(m.rows()), c_beg(0), c_end(m.cols()), m(m)
  {

  }


  const T& operator()(int r, int s) const
  {
    return m.cells[r_beg + r][c_beg + s];
  }

  T& operator()(int r, int s)
  {
    return m.cells[r_beg + r][c_beg + s];
  }

  int rows() const
  {
    return r_end - r_beg;
  }

  int cols() const
  {
    return c_end - c_beg;
  }

  static matrix<T> slow_mul(const matrix_view<T>& a, const matrix_view<T>& b)
  {
    matrix<T> c(a.rows(), b.cols());

    for (int i = 0; i < a.rows(); ++i)
      for (int j = 0; j < b.cols(); ++j)
        for (int k = 0; k < a.cols(); ++k)
        {
          if ((i == 0) && (j == 3))
          {
            //cout << "a";
          }
          c[i][j] += a(i, k) * b(k, j);
        }

    return c;
  }

  matrix<T> operator* (const matrix_view<T>& o)
  {
    return mul_op(*this, o);
  }

  matrix_view operator += (const matrix_view<T>& other)
  {
    for (int i = 0; i < other.rows(); ++i)
      for (int j = 0; j < other.cols(); ++j)
      {
        (*this)(i, j) += other(i, j);
      }
    return *this;
  }

  function<matrix<T>(const matrix_view&, const matrix_view&)> mul_op = slow_mul;

  void emplace(matrix_view<T> other, int at_r, int at_c)
  {
    for (int i = 0; i < other.rows(); ++i)
      for (int j = 0; j < other.cols(); ++j)
        (*this)(at_r + i, at_c + j) = other(i, j);
  }

  matrix<T> to_matrix()
  {
    matrix<T> r(rows(), cols());
    for (int i = 0; i < rows(); ++i)
      for (int j = 0; j < cols(); ++j)
        r[i][j] = (*this)(i, j);
    return r;
  }


};

typedef matrix_view<clause> valiant_view;

/*template <typename T>
matrix<T>::operator matrix_view<T>() const
{
  return matrix_view<T>(*this);
}*/

template <class T>
ostream& operator<< (ostream& ostr, const matrix_view<T> &mat)
{
  for (int i = 0; i < mat.rows(); ++i, ostr << endl)
    for (int j = 0; j < mat.cols(); ++j)
      ostr << mat(i,j) << "\t\t\t";
  return ostr;
}

template <class T>
ostream& operator<< (ostream& ostr, const matrix<T> &_mat)
{
  matrix<T> __mat = _mat; // mat_view doesn't guarantee constness
  matrix_view<T> mat(__mat);
  ostr << mat;
  return ostr;
}

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


int main()
{
  ifstream ist("/home/luka/bin/valiant_parser/prop_log.txt");
  //ifstream ist("/home/luka/bin/valiant_parser/aaa.txt");
  read_from_stream(ist);

  valiant test = setup_parse_mat("(p+(r*(p+s)))"); // setup_parse_mat("((p+q)*((q+p)+r))");
  valiant_view base(test);

  //cout << base;

  closure(base);

  cout << base;

  //valiant potencija = test;
  //valiant_view whole(potencija);

  //cout << v.without_range(2,3,1,2);

  /*for (int i = 0; i < 25; ++i)
  {
    cout << whole << endl;

    valiant tmp = whole * whole;
    valiant_view whole_tmp(tmp);
    whole += whole_tmp;
  }*/

  return 0;
}