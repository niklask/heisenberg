/* **********************************
    Quantum Heisenber Model

    a Monte Carlo Simulation
    
    File: fileinput.h

    Copyright ¨ 2001 Niklas Karlsson
   ********************************** */

#ifndef _FileInput_H_
#define _FileInput_H_

#include <string>
#include <map>

using namespace std;

typedef map<string, string> INPUTDATA;

class FileInput {
  
    // Public methods (and constructors)
  public:
    FileInput() {};
    FileInput(const char* filename);
    virtual ~FileInput();
  
    virtual int parse(const char* filename);
    virtual double getDoubleValue(string variable);
    virtual int getIntValue(string variable);
    virtual bool getBooleanValue(string variable);
    virtual string getStringValue(string variable);
  
    // Private fields
  private:
    // Input variables read
    INPUTDATA input;
};

#endif 
