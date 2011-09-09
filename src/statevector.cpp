/* **********************************
    Quantum Heisenber Model

    a Monte Carlo Simulation
  
    File: statevector.cpp

    Copyright © 2001 Niklas Karlsson
   ********************************** */

using namespace std;

#include <iostream>

#include "random.h"
#include "statevector.h"

StateVector::StateVector()
{
    vector = NULL;
}

StateVector::StateVector(int size)
{
    vector = new STATEVECTOR(size);
    null();
}

StateVector::StateVector(int size, bool random)
{
    vector = new STATEVECTOR(size);
    null();

    if (random == true)
        for (int i = 0; i < size; i++)
            if (ran() < 0.5)
                flip(i);
}

StateVector::StateVector(int size, int s)
{
    vector = new STATEVECTOR(size);
    null();

    double r = s;
    for (int i = 0; i < vector->size; i++) {
        if (r < 1)
            vector->states[vector->size - (i+1)] = 0;
        else {
            vector->states[vector->size - (i+1)] = ((int)(r) % 2);
            r = (r / 2) - ((int)(r) % 2)*0.5;
        }
    }
}

// Copy constructor
StateVector::StateVector(const StateVector& s)
{
    vector = new STATEVECTOR(s.vector->size);
    vector->c = s.vector->c;
    for (int i = 0; i < (s.vector->size); i++)
        vector->states[i] = s.vector->states[i];
}

StateVector::~StateVector()
{
    if (vector != NULL)
        delete(vector);
}

void StateVector::null()
{
    vector->c = 1;
    for (int i = 0; i < vector->size; i++)
        vector->states[i] = 0;
}

int StateVector::getSize()
{
    return vector->size;
}

StateVector& StateVector::Sz(int n)
{
    if (vector->states[n] == 1) {
        // Operating on spin down
        vector->c *= -0.5;
    } else if (vector->states[n] == 0) {
        // Operating on spin up
        vector->c *= 0.5;
    }
    return *this;
}

StateVector& StateVector::Sp(int n)
{
    if (vector->states[n] == 0) {
        // Operating on spin up
        vector->c *= 0;
    } else if (vector->states[n] == 1) {
        // Operating on spin down
        vector->states[n] = 0;
    }
    return *this;
}

StateVector& StateVector::Sm(int n)
{
    if (vector->states[n] == 1) {
        // Operating on spin down
        vector->c *= 0;
    } else if (vector->states[n] == 0) {
        // Operating on spin up
        vector->states[n] = 1;
    }
    return *this;
}

void StateVector::flip(int i)
{
    // Flip the spin at site i
    if (vector->states[i] == 1) {
        vector->states[i] = 0;
    } else {
        vector->states[i] = 1;
    }
}

ostream& operator << (ostream& os, const StateVector& s)
{
    os << s.vector->c << "*[";
    for (int i = 0; i < s.vector->size; i++)
        os << s.vector->states[i];
    os << "]";
    return os;
}

int& StateVector::operator () (int i)
{
    return vector->states[i];
}

const int StateVector::operator () (int i) const
{
    return vector->states[i];
}

StateVector& StateVector::operator = (const StateVector& s)
{
    if (vector != NULL)
        delete(vector);

    vector = new STATEVECTOR(s.vector->size);
    for (int i = 0; i < (s.vector->size); i++)
        vector->states[i] = s.vector->states[i];

    return *this;
}

bool operator == (const StateVector& s1, const StateVector& s2)
{
    for (int i = 0; i < s1.vector->size; i++) {
        if (s1.vector->states[i] != s2.vector->states[i])
            return false;
    }

    return true;
}

bool operator != (const StateVector& s1, const StateVector& s2)
{
    if (s1 == s2)
        return false;
    else
        return true;
}

double operator | (const StateVector& s1, const StateVector& s2)
{
    if (s1 == s2)
        return (s1.vector->c * s2.vector->c);
    else
        return 0;
}
