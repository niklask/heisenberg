/* **********************************
    Quantum Heisenberg Model

    a Monte Carlo Simulation
   
    File: spindisplay.cpp

    Copyright © 2001 Niklas Karlsson
   ********************************** */

using namespace std;

#include <iostream>

//#include "up.xpm"
//#include "down.xpm"

//#include <stdlib.h>
//#include <unistd.h>

// include the Qt headers
#include <qapplication.h>
#include <qcanvas.h>

#include "statevector.h"
#include "operatorlist.h"
#include "spindisplay.h"

SpinDisplay::SpinDisplay(StateVector* stateVector, OperatorList* operatorList)
{
    this->stateVector = stateVector;
    this->operatorList = operatorList;

    InitDisplay();        
}

SpinDisplay::~SpinDisplay()
{
}

int SpinDisplay::InitDisplay()
{
    int spins = stateVector->getSize();
    int length = operatorList->getSize();

    int width = spins * 20 + (spins - 1) * 10 + 40;
    int height = length * 20 + (length - 1) * 10 + 40;
    return 0;
}

