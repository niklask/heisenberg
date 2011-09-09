/* **********************************
   Quantum Heisenber Model

   a Monte Carlo Simulation

   File: operatorlist.h

   Copyright ¨ 2001 Niklas Karlsson
  ********************************** */

#ifndef _OperatorList_H_
#define _OperatorList_H_

typedef struct OPERATOR {
    int op;
    int spin;
};

class OperatorList {

    // Public methods (and constructors)
  public:
    OperatorList();
    OperatorList(int size);
    OperatorList(const OperatorList& o);
    virtual ~OperatorList();
  
    virtual void unit();

    virtual int getSize();
    virtual int countDOD();
    virtual int countDiagonal();
    virtual int countOffDiagonal();

    virtual void diagonalUpdate(StateVector& sv, double beta);
    virtual void loopUpdate(StateVector& sv, double beta);
    virtual void changeOp(int i);

    virtual void display(StateVector sv);
    virtual bool check(StateVector sv);
 
    // Insertion and extraction operators
    friend ostream& operator << (ostream& os, const OperatorList& o); 
 
    // Access operators
    OPERATOR& operator () (int i);
    const OPERATOR& operator () (int i) const;

    // Private fields
  private:

    typedef struct OPERATORLIST {
        //int n;  // n is the no of diagonal and off-diagonal operators in the list
        OPERATOR* operators;
        int size;

        OPERATORLIST(int size) : size (size) {
            //n = 0;
            operators = new OPERATOR[size];
        }

        ~OPERATORLIST(){
            delete[] operators;
        }
    };

    OPERATORLIST* list;
};

#endif
