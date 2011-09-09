/* **********************************
  Quantum Heisenber Model

  a Monte Carlo Simulation

  File: statevector.h

  Copyright ¨ 2001 Niklas Karlsson
  ********************************** */

#ifndef _StateVector_H_
#define _StateVector_H_

class StateVector {

    // Public methods (and constructors)
  public:
    StateVector();
    StateVector(int size);
    StateVector(int size, int s);
    StateVector(int size, bool random);
    StateVector(const StateVector& s);
    virtual ~StateVector();
  
    virtual void null();

    virtual int getSize();

    virtual StateVector& Sz(int n);
    virtual StateVector& Sp(int n);
    virtual StateVector& Sm(int n);
    virtual void flip(int i);

    // Insertion and extraction operators
    friend ostream& operator << (ostream& os, const StateVector& s); 
    //friend istream& operator >> (istream& is, Matrix& m); 
 
    // Access operators
    int& operator () (int i);
    const int operator () (int i) const;

    // Assignment operators
    StateVector& operator = (const StateVector& s);
  
    // Logical operators
    friend bool operator == (const StateVector& s1, const StateVector& s2);
    friend bool operator != (const StateVector& s1, const StateVector& s2);

    // Binary operators
    friend double operator | (const StateVector& s1, const StateVector& s2);

    // Private fields
  private:

    typedef struct STATEVECTOR {
        double c;
        int* states;
        int size;
        
        STATEVECTOR(int size) : size (size) {
            c = 1;
            states = new int[size];
        }
        
        ~STATEVECTOR(){
            delete[] states;
        }
    };

    STATEVECTOR* vector;
};

#endif
