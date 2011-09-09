/* **********************************
    Quantum Heisenberg Model

    a Monte Carlo Simulation
   
    File: operatorlist.cpp

    Copyright © 2001 Niklas Karlsson
   ********************************** */

// Available operators:
//   0 - unit operator
//   1 - diagonal operator
//   2 - off-diagonal operator

using namespace std;

#include <iostream>
#include <vector>

#include "random.h"
#include "statevector.h"
#include "operatorlist.h"

// -------------------------------------------------------------------------------------------------
//  OperatorList::OperatorList()
//
//  Use: default constructor
//  Input: none
//  Output: none
// -------------------------------------------------------------------------------------------------
OperatorList::OperatorList()
{
    list = NULL;
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::OperatorList(int size)
//
//  Use: main constructor
//  Input: int size - the size of the operator list
//  Output: none
// -------------------------------------------------------------------------------------------------
OperatorList::OperatorList(int size)
{
    list = new OPERATORLIST(size);
    unit();    
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::OperatorList(const OperatorList& o)
//
//  Use: copy constructor
//  Input: const OperatorList& o - the operator list to copy
//  Output: none
// -------------------------------------------------------------------------------------------------
OperatorList::OperatorList(const OperatorList& o)
{
    list = new OPERATORLIST(o.list->size);
    for (int i = 0; i < (o.list->size); i++)
        list->operators[i] = o.list->operators[i];
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::~OperatorList()
//
//  Use: default destructor
//  Input: none
//  Output: none
// -------------------------------------------------------------------------------------------------
OperatorList::~OperatorList()
{
    if (list != NULL)
        delete(list);
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::unit()
//
//  Use: Sets all operators to unit operators
//  Input: none
//  Output: none
// -------------------------------------------------------------------------------------------------
void OperatorList::unit()
{
    for (int i = 0; i < list->size; i++) {
        list->operators[i].op = 0;
        // unit operators don't act on spins, therefore set spin to a neg. value
        list->operators[i].spin = -1;
    }
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::getSize()
//
//  Use: Returns the size if the operator list
//  Input: none
//  Output: int - the size as integer
// -------------------------------------------------------------------------------------------------
int OperatorList::getSize()
{
    return list->size;
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::countDiagonal()
//
//  Use: Count the number of diagonal operators
//  Input: none
//  Output: int - the number of diagonal operators
// -------------------------------------------------------------------------------------------------
int OperatorList::countDiagonal()
{
    int n = 0;
    for (int i = 0; i < list->size; i++) 
        if (list->operators[i].op == 1)
            n++;

    return n;
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::countDiagonal()
//
//  Use: Count the number of off-diagonal operators
//  Input: none
//  Output: int - the number of off-diagonal operators
// -------------------------------------------------------------------------------------------------
int OperatorList::countOffDiagonal()
{
    int n = 0;
    for (int i = 0; i < list->size; i++) 
        if (list->operators[i].op == 2)
            n++;

    return n;
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::countDOD()
//
//  Use: Count the number of diagonal AND off-diagonal operators
//  Input: none
//  Output: int - the number of diagonal AND off-diagonal operators
// -------------------------------------------------------------------------------------------------
int OperatorList::countDOD()
{
    int n = 0;
    for (int i = 0; i < list->size; i++) 
        if (list->operators[i].op > 0)
            n++;

    return n;
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::diagonalMove(StateVector& sv, double beta)
//
//  Use: Performs a diagonal update on the operator list
//  Input: StateVector& sv - the spin configuration
//         double beta - the inverse temperature
//  Output: none
// -------------------------------------------------------------------------------------------------
void OperatorList::diagonalUpdate(StateVector& sv, double beta)
{
    int N = sv.getSize();       // Number of spins
    int L = getSize();          // Operator list length

    for (int i = 0; i < L; i++) {
        int n = countDOD();   // Number of diagonal and off-diagonal operators
        int op = list->operators[i].op;
        if (op == 0) {
            // unit operator
            // Calculate probability to change operator
            double P = (double)(N*beta/(2*(L-n)));
            // Generate a random spin that this operator acts on
            int spin1 = (int)(ran()*N);
            // and its adjacent spin
            int spin2 = (spin1 + 1) % N;
            
            // Change operator to diagonal only if the spins are _not_ parallel...
            if (sv(spin1) != sv(spin2)) {
                // ...and with probability P
                double r = ran();

                if (r < P) {
                    //cout << "Change..." << endl;
                    list->operators[i].op = 1;
                    list->operators[i].spin = spin1;
                }
            }
        } else if (op == 1) { 
            // diagonal operator
            // Calculate probability to change to unit operator
            double P = (double)(2*(L-n+1)/(N*beta));
            double r = ran();
            
            // Get the spins this off-diagonal operator acts on
            int spin1 = list->operators[i].spin;
            int spin2 = (spin1 + 1) % N;
            
            if (r < P) {
                list->operators[i].op = 0;    
                // unit ops cannot act on spins
                list->operators[i].spin = -1;    
            }
        } else if (op == 2) { 
            // off-diagonal operator
            // Get the spins this off-diagonal operator acts on
            int spin1 = list->operators[i].spin;
            int spin2 = (spin1 + 1) % N;
            
            // Flip both spins
            sv.flip(spin1);
            sv.flip(spin2);                
        }
    }
} 

// -------------------------------------------------------------------------------------------------
//  OperatorList::loopUpdate(StateVector& sv, double beta)
//
//  Use: Performs a loop cluster update on the operator list
//  Input: StateVector& sv - the spin configuration
//         double beta - the inverse temperature
//  Output: none
// -------------------------------------------------------------------------------------------------
void OperatorList::loopUpdate(StateVector& sv, double beta)
{
    int N = sv.getSize();     // Number of spins
    int L = getSize();        // Operator list length

    int noDOD = countDOD();
    if (noDOD > 0) {
        // j is the operator loop index
        int j;
        // Select a random non-unit operator...
        do {
            j = (int)(ran()*L);
            if (j > (L - 1)) 
              j = L - 1;
        } while (list->operators[j].op == 0); 

        // i is the state loop index
        int i = j;

        // direction = -1 : backward direction
        // direction = 1 : forward direction
        int direction = 1;
        int loops = 0;
        int spin = list->operators[j].spin;
        int startop = j;
        int startstate = i;
        int startspin = spin;
        int spin1 = spin;
        int spin2 = (spin1 + 1) % N;

        // change this operator from diagonal to off-diagonal and vice versa...
        changeOp(j);

        // this is the actual loop update
        do {
            if (loops > 0) {
                // Break the loop when we return to the starting operator
                if ((i == startstate) &&(j == startop))
                    if (spin == startspin) {
                        changeOp(j);
                        break;
                    } else if ((spin - 1) == startspin)
                        break;
            }

            j += direction;
            // we use periodic boundary conditions
            if (j < 0) {
                j = L - 1;
                sv.flip(spin);
            }
            if (j > (L - 1)) {
                j = 0;
                sv.flip(spin);
            }

            int op = list->operators[j].op;
            spin1 = list->operators[j].spin;
            spin2 = (spin1 + 1) % N;

            if (op == 0) {
                // in case of a unit operator just update state index
                i += direction;
                if (i < 0) i = L - 1;
                if (i > (L-1)) i = 0;
            } else {
                // diagonal or off-diagonal operator
                if ((spin != spin1) && (spin != spin2)) {
                    i += direction;
                    if (i < 0) i = L - 1;
                    if (i > (L-1)) i = 0;
                } else if (spin == spin1) {
                    spin++;
                    if (spin > (N - 1))
                        spin = 0;
                    direction *= -1;
                    changeOp(j);
                } else if (spin == spin2) {
                    spin--;
                    if (spin < 0)
                        spin = N - 1;
                    direction *= -1;
                    changeOp(j);
                }
            }

            loops = 1;
        } while(true);
    }
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::changeOp(int i)
//
//  Use: Change an operator from diagonal to off-diagonal and vice versa
//  Input: int i - operator to change
//  Output: none
// -------------------------------------------------------------------------------------------------
void OperatorList::changeOp(int i)
{
    if (list->operators[i].op == 1)
        list->operators[i].op = 2;
    else if (list->operators[i].op == 2)
        list->operators[i].op = 1;
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::display(StateVector sv)
//
//  Use: Output a text visualisation of the operator list
//  Input: StateVector sv - the spin configuration
//  Output: none
// -------------------------------------------------------------------------------------------------
void OperatorList::display(StateVector sv)
{
    int N = sv.getSize();
    int L = getSize();

    const char* spins[2] = {"U", "D"};
    const char* ops[3] = {" ", "-", "="};

    for (int i = 0; i < L; i++) {
        cout << " ";
        for (int j = 0; j < N; j++)
            cout << spins[sv(j)] << " ";
        cout  << " (" << i << ")" << endl;

        int op = list->operators[i].op;
        if (op == 0) 
            cout << endl;
        else {
            int spin1 = list->operators[i].spin;
            int spin2 = spin1 + 1;
            if (spin2 > (N - 1))
                spin2 = 0;
            
            if (spin2 == 0) {
                cout << " " << ops[op];
                for (int k = 0; k < (2*N - 3); k++)
                    cout << " ";
                cout << ops[op] << endl;
            } else {
                for (int k = 0; k < (2*spin1 + 1); k++)
                    cout << " ";            
                for (int k = 0; k < 3; k++)
                    cout << ops[op];
                cout << endl;
            }

            // flip spins if off diagonal operator
            if (op == 2) {
                sv.flip(spin1);
                sv.flip(spin2);
            }            
        }
    }
    cout << " ";
    for (int j = 0; j < N; j++)
        cout << spins[sv(j)] << " ";
    cout << endl << endl << endl;
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::check(StateVector sv)
//
//  Use: Check the operator list for consitency
//  Input: StateVector sv - the spin configuration
//  Output: bool - true if ok, false if not
// -------------------------------------------------------------------------------------------------
bool OperatorList::check(StateVector sv)
{
    // Check that no diagonal or off-diagonal operator acts
    // on two parallel spins
    int N = sv.getSize();
    int L = getSize();

    for (int i = 0; i < L; i++) {
        int op = list->operators[i].op;
        int spin1 = list->operators[i].spin;
        int spin2 = spin1 + 1;
        if (spin2 > (N - 1))
            spin2 = 0;

        if ((op > 0) && (sv(spin1) == sv(spin2))) {
            cout << "Error: Operating on parallel spins!!" << endl;
            cout << "operator " << op << " on " << spin1 << "," << spin2;
            cout << " in " << sv << endl;
            cout << (*this) << endl;
            return false;
        }

        if (op == 2) {
            sv.flip(spin1);
            sv.flip(spin2);
        }
    }

    return true;
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::operator () (int i)
//
//  Use: Access operator to get a certain operator from the list
//  Input: int i - the operator index
//  Output: OPERATOR& - the operator with spin
// -------------------------------------------------------------------------------------------------
OPERATOR& OperatorList::operator () (int i)
{
    return list->operators[i];
}

// -------------------------------------------------------------------------------------------------
//  OperatorList::operator () (int i) const
//
//  Use: Access operator to get a certain operator from the list
//  Input: int i - the operator index
//  Output: const OPERATOR& - the operator with spin
// -------------------------------------------------------------------------------------------------
const OPERATOR& OperatorList::operator () (int i) const
{
    return list->operators[i];
}

// -------------------------------------------------------------------------------------------------
//  operator << (ostream& os, const OperatorList& o)
//
//  Use: Creates output for the operator list
//  Input: ostream& os - stream to output to
//         const OperatorList& o - the operator list to output
//  Output: ostream& - the stream
// -------------------------------------------------------------------------------------------------
ostream& operator << (ostream& os, const OperatorList& o)
{
    for (int i = 0; i < o.list->size; i++) {
        os << "[" << o.list->operators[i].op << ",";
        os << o.list->operators[i].spin << "] ";
    }
    return os;
}
