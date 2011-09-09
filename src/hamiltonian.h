/* **********************************
    Quantum Heisenber Model

    a Monte Carlo Simulation

    File: hamiltonian.h

    Copyright ¨ 2001 Niklas Karlsson
  ********************************** */

#ifndef _Hamiltonian_H_
#define _Hamiltonian_H_

class Hamiltonian : public Matrix {

  public:
    Hamiltonian(int size);

  private:
    virtual double Hmn(int m, int n, int size);

};

#endif
