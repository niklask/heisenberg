/* **********************************
    Quantum Heisenber Model

    a Monte Carlo Simulation

    File: hamiltonian.cpp

    Copyright ¨ 2001 Niklas Karlsson
  ********************************** */

using namespace std;

#include <iostream>
#include <math.h>

#include "matrix.h"
#include "hamiltonian.h"
#include "statevector.h"

// ---------------------------------------------------------------------------------
//  Hamiltonian::Hamiltonian(int size)
//
// Create the Hamiltonian of given size 
// ---------------------------------------------------------------------------------
Hamiltonian::Hamiltonian(int size) : Matrix((int)(exp(size*log(2.))) + 1)
{
    int msize = (int)(exp(size*log(2.))) + 1;

    for (int m = 0; m < msize; m++)
        for (int n = 0; n < msize; n++)
            (*this)(m+1, n+1) = Hmn(m, n, size);
}

// ---------------------------------------------------------------------------------
//  Hamiltonian::Hmn(int m, int n, int size)
//
// Create matrix element mn of the Hamiltonian
// ---------------------------------------------------------------------------------
double Hamiltonian::Hmn(int m, int n, int size)
{
    StateVector sl(size, m);
    StateVector sr(size, n);
    StateVector tmp;

    double H = 0;

    for (int i = 0; i < size; i++) {
        int j = (i + 1) % size;

        tmp = sr;
        tmp.Sz(i); tmp.Sz(j);
        H += (sl | tmp);

        tmp = sr;
        tmp.Sp(i); tmp.Sm(j);
        H += 0.5 * (sl | tmp);

        tmp = sr;
        tmp.Sp(j); tmp.Sm(i);
        H += 0.5 * (sl | tmp);
    }

    return H;
}
