/* **********************************
    Quantum Heisenber Model

    a Monte Carlo Simulation

    File: matrix.cpp

    Copyright ¨ 2001 Niklas Karlsson
   ********************************** */

using namespace std;

#include <iostream>
#include <math.h>
#include "matrix.h"

Matrix::Matrix()
{
    matrix = NULL;
}

// Create n*m matrix
Matrix::Matrix(int m, int n)
{
    matrix = new MATRIXDATA(m, n);
    Zero();
}

// Create a square n*n matrix
Matrix::Matrix(int m)
{
    matrix = new MATRIXDATA(m, m);
    Zero();
}

// Copy constructor
Matrix::Matrix(const Matrix& m)
{
    matrix = new MATRIXDATA(m.matrix->n, m.matrix->m);
    matrix->diagonal = m.matrix->diagonal;
    for (int i = 0; i < (m.matrix->n*m.matrix->m); i++)
        matrix->elements[i] = m.matrix->elements[i];
}

Matrix::~Matrix()
{
    if (matrix != NULL) 
        delete(matrix);
}

void Matrix::Zero()
{
    int size = matrix->n*matrix->m;
    for (int i = 0; i < size; i++)
        matrix->elements[i] = 0;
}

void Matrix::Identity()
{
    // It's only possible to have square identity matrices
    if (matrix->n == matrix->m)
        for (int i = 0; i < matrix->n; i++)
            matrix->elements[i*matrix->m + i] = 1;
}

double Matrix::Trace()
{
    double trace = 0;

    // Can only calculate the trace of square matrices
    if (matrix->n == matrix->m)
        for (int i = 0; i < matrix->n; i++)
            trace += matrix->elements[i*matrix->m + i];
 
    return trace; 
}

void Matrix::Transpose()
{
    Matrix tmp(matrix->m, matrix->n);

    for (int i = 0; i < matrix->n; i++)
        for (int j = 0; j < matrix->m; j++)
            tmp(j,i) = (*this)(i,j);

    MATRIXDATA* tmpmatrix = matrix;
    matrix = tmp.matrix;
    tmp.matrix = tmpmatrix;
}

void Matrix::Diagonalize()
{
#define ROTATE(a,i,j,k,l) g = a(i, j);h = a(k, l);a(i, j) = g - s*(h + g*tau);a(k, l) = h + s*(g - h*tau);
    double *b, *d, *z;
    double tresh, theta, tau, t, sm, s, h, g, c;
    int i, j, ip, iq;

    // Can only diagonalize symmetric (square) matrices
    if (matrix->n == matrix->m) {
        int n = matrix->n;
  
        Matrix v(matrix->n);
        v.Identity();

        b = new double[n + 1];
        d = new double[n + 1];
        z = new double[n + 1];

        for (ip = 1; ip <= n; ip++) {
            b[ip] = d[ip] = (*this)(ip, ip);
            z[ip] = 0.0;
        }

        for (i = 1; i <= 50; i++) {
            sm =0.0;
            for (ip = 1; ip <= (n - 1); ip++) {
                for (iq = ip+1; iq <= n; iq++)
                    sm += fabs((*this)(ip, iq));
            }

            if (sm == 0.0) {
                matrix->diagonal = true;
                Zero();
                Identity();
                for (ip = 1; ip <= n; ip++)
                    (*this)(ip, ip) = d[ip];

                delete[] b;
                delete[] d;
                delete[] z;

                return;
            }

            if (i < 4)
                tresh = 0.2*(sm/(n*n)); 
            else
                tresh = 0.0;

            for (ip = 1; ip <= (n - 1); ip++) {
                for (iq = ip+1; iq <= n; iq++) {
                    g = 100.0*fabs((*this)(ip, iq));
                    if ((i > 4) && ((double)(fabs(d[ip])+g) == (double)fabs(d[ip])) &&
                        ((double)(fabs(d[iq])+g) == (double)fabs(d[iq])))
                        (*this)(ip, iq) = 0.0;
                    else if (fabs((*this)(ip, iq)) > tresh) {
                        h = d[iq] - d[ip];
                        if ((double)(fabs(h) + g) == (double)(fabs(h)))
                            t = (*this)(ip, iq)/h;
                        else {
                            theta = 0.5*h/((*this)(ip, iq));
                            t = 1.0/(fabs(theta)+sqrt(1.0 + theta*theta));
                            if (theta < 0.0)
                                t = -t;
                        }
                        c = 1.0/sqrt(1.0 + t*t);
                        s = t*c;
                        tau = s/(1.0 + c);
                        h = t*(*this)(ip, iq);
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= h;
                        d[iq] += h;
                        (*this)(ip, iq) = 0.0;

                        for (j = 1; j <= (ip - 1); j++) {
                            ROTATE((*this),j,ip,j,iq);
                        }
                        for (j = ip+1; j <= (iq - 1); j++) {
                            ROTATE((*this),ip,j,j,iq);
                        }
                        for (j = iq+1; j <= n; j++) {
                            ROTATE((*this),ip,j,iq,j);
                        }
                        for (j = 1; j <= n; j++) {
                            ROTATE(v,j,ip,j,iq);
                        }
                    }
                }
            }

            for (ip = 1; ip <= n; ip++) {
                b[ip] += z[ip];
                d[ip] = b[ip];
                z[ip] = 0.0;
            }
        }
    }
}

void Matrix::Exp()
{
    int m = matrix->m;
    int n = matrix->n;

    if (matrix->diagonal == true) {
        for (int i = 1; i <= m; i++)
            (*this)(i, i) = exp((*this)(i, i));
    }
}

std::ostream& operator << (std::ostream& os, const Matrix& m)
{
    for (int i = 0; i < m.matrix->m; i++) {
        for (int j = 0; j < m.matrix->n; j++)
            os << m(i+1, j+1) << "\t";
        os << endl;
    }
    return os;
}

// row i, column j
double& Matrix::operator () (int i, int j)
{
    if ((i >= 1 && i <= matrix->m) && (j >= 1 && j <= matrix->n))
        return matrix->elements[(i-1)*matrix->m + (j-1)];
}

const double Matrix::operator () (int i, int j) const
{
    if ((i >= 1 && i <= matrix->m) && (j >= 1 && j <= matrix->n))
        return matrix->elements[(i-1)*matrix->m + (j-1)];
}

Matrix& Matrix::operator = (const Matrix& m)
{
    if (matrix != NULL)
        delete(matrix);

    matrix = new MATRIXDATA(m.matrix->n, m.matrix->m);
    for (int i = 0; i < (m.matrix->n*m.matrix->m); i++)
        matrix->elements[i] = m.matrix->elements[i];

    return *this;
}

Matrix& Matrix::operator += (const Matrix& m)
{
    if ((m.matrix->n == matrix->n) && (m.matrix->m == matrix->m))
        for (int i = 0; i < (matrix->n*matrix->m); i++)
            matrix->elements[i] += m.matrix->elements[i];
 
    return *this;
}

Matrix& Matrix::operator -= (const Matrix& m)
{
    if ((m.matrix->n == matrix->n) && (m.matrix->m == matrix->m))
        for (int i = 0; i < (matrix->n*matrix->m); i++)
            matrix->elements[i] -= m.matrix->elements[i];
 
    return *this;
}

Matrix& Matrix::operator *= (const Matrix& m)
{
    if (matrix->m == m.matrix->n)
        *this = *this * m;
 
    return *this;
}

Matrix& Matrix::operator *= (double s)
{
    for (int i = 0; i < (matrix->n*matrix->m); i++)
        matrix->elements[i] *= s;
 
    return *this;
}

Matrix& Matrix::operator /= (double s)
{
    for (int i = 0; i < (matrix->n*matrix->m); i++)
        matrix->elements[i] /= s;
 
    return *this;
}

Matrix Matrix::operator + ()
{
    return *this;
}

Matrix Matrix::operator - ()
{
    Matrix ret(matrix->n, matrix->m);
 
    for (int i = 1; i <= matrix->n; i++)
        for (int j = 1; j <= matrix->m; j++)
            ret(i,j) = -(*this)(i,j);

    return ret;
}

Matrix operator + (const Matrix& m1, const Matrix& m2)
{
    if ((m1.matrix->n == m2.matrix->n) && (m1.matrix->m == m2.matrix->m)) {
        Matrix ret(m1.matrix->n, m1.matrix->m);

        for (int i = 1; i <= m1.matrix->n; i++)
            for (int j = 1; j <= m1.matrix->m; j++)
                ret(i,j) = m1(i,j) + m2(i,j);

        return ret;
    } else
        return m1;
}

Matrix operator - (const Matrix& m1, const Matrix& m2)
{
    if ((m1.matrix->n == m2.matrix->n) && (m1.matrix->m == m2.matrix->m)) {
        Matrix ret(m1.matrix->n, m1.matrix->m);

        for (int i = 1; i <= m1.matrix->n; i++)
            for (int j = 1; j <= m1.matrix->m; j++)
                ret(i,j) = m1(i,j) - m2(i,j);

        return ret;
    } else
        return m1;
}

Matrix operator * (const Matrix& m1, const Matrix& m2)
{
    if (m1.matrix->m == m2.matrix->n) {
        Matrix ret(m1.matrix->n, m2.matrix->m);

        for (int i = 1; i <= m1.matrix->n; i++)
            for (int j = 1; j <= m1.matrix->m; j++) {
                ret(i,j) = 0;
                for (int k = 1; k <= m1.matrix->m; k++)
                    ret(i,j) += m1(i,k) * m2(k,j);
            }
        return ret;
    } else
        return m1;
}

Matrix operator * (const Matrix& m, const double& s)
{
    Matrix ret(m.matrix->n, m.matrix->m);
 
    for (int i = 1; i <= m.matrix->n; i++)
        for (int j = 1; j <= m.matrix->m; j++)
            ret(i,j) = m(i,j) * s;

    return ret;
}

Matrix operator * (const double& s, const Matrix& m)
{
    return (m * s);
}

