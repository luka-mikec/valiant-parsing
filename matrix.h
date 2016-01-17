#ifndef VALIANT_PARSING_MATRIX_H
#define VALIANT_PARSING_MATRIX_H

#include <algorithm>
#include <iostream>
#include <functional>
#include <stdexcept>
#include <valarray>

using std::endl;
using std::istream;
using std::function;
using std::ostream;
using std::runtime_error;
using std::valarray;

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
};


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

#endif //VALIANT_PARSING_MATRIX_H
