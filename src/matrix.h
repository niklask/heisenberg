/* **********************************
    Quantum Heisenber Model

    a Monte Carlo Simulation
    
    File: matrix.h
    
    Copyright ¨ 2001 Niklas Karlsson
   ********************************** */

#ifndef _Matrix_H_
#define _Matrix_H_

class Matrix {

    // Public methods (and constructors)
  public:
    Matrix();
    Matrix(int m, int n);
    Matrix(int m);
    Matrix(const Matrix& m);
    virtual ~Matrix();
  
    virtual void Identity();
    virtual void Zero();
    virtual double Trace();
    virtual void Transpose();
    virtual void Diagonalize();
    virtual void Exp(); 

    // Insertion and extraction operators
    friend std::ostream& operator << (std::ostream& os, const Matrix& m); 
    //friend istream& operator >> (istream& is, Matrix& m); 

    // Access operators
    double& operator () (int i, int j);
    const double operator () (int i, int j) const;

    // Assignment operators
    Matrix& operator = (const Matrix& m);
    Matrix& operator += (const Matrix& m);
    Matrix& operator -= (const Matrix& m);
    Matrix& operator *= (const Matrix& m);
    Matrix& operator *= (double s);
    Matrix& operator /= (double s);

    // Unary operators
    Matrix operator + ();
    Matrix operator - ();
  
    // Binary operators
    friend Matrix operator + (const Matrix& m1, const Matrix& m2);
    friend Matrix operator - (const Matrix& m1, const Matrix& m2);
  
    friend Matrix operator * (const Matrix& m1, const Matrix&m2);
    friend Matrix operator * (const Matrix& m, const double& s);
    friend Matrix operator * (const double& s, const Matrix& m);
   
    // Private fields
  private:
    typedef struct MATRIXDATA {
        double* elements;
        int m, n;
        bool diagonal;
        
        MATRIXDATA(int m, int n) : m (m), n (n) {
            diagonal = false;
            elements = new double[m*n];
        }
        
        ~MATRIXDATA(){
            delete[] elements;
        }
    };

    // The actual matrix
    MATRIXDATA* matrix;
};

#endif
