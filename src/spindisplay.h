/* **********************************
   Quantum Heisenber Model

   a Monte Carlo Simulation

   File: spindisplay.h

   Copyright ¨ 2001 Niklas Karlsson
  ********************************** */

#ifndef _SpinDisplay_H_
#define _SpinDisplay_H_

class SpinDisplay {
    
  // Public methods (and constructors)
  public:
    SpinDisplay(StateVector* stateVector, OperatorList* operatorList);
    virtual ~SpinDisplay();

    virtual int InitDisplay();
    //virtual void Update();

  private:
    StateVector* stateVector;
    OperatorList* operatorList;    
};

#endif
