/* **********************************
    Quantum Heisenber Model

    a Monte Carlo Simulation

    File: fileinput.cpp

    Copyright ¨ 2001 Niklas Karlsson
   ********************************** */

#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <iostream>

using namespace std;

#include "fileinput.h"

FileInput::FileInput(const char* filename)
{
    if (parse(filename) != 0)
        cout << "Parse error in file " << filename << endl;
}

FileInput::~FileInput()
{
}

int FileInput::parse(const char* filename)
{
    string parameter, value, buffer;

    ifstream file(filename);
    
    getline(file, buffer);
    while (!file.eof()) {
        // only parse line if is not empty
        if (!buffer.empty()) {
            // parse only if first character is not a comment sign
            if ((buffer.substr(0, 1)).compare(";") != 0) {
                // create a stream on the input buffer
                istringstream ist(buffer.c_str());
                // get first token as parameter name
                ist >> parameter;
                ist >> value;
                // is there an equal sign after the parameter name
                if (value.compare("=") == 0) {
                    // if so, get the value as third token
                    ist >> value;
                    // add parameter to map<>
                    input[parameter] = value;
                } else {
                    // parse error
                    return 1;
                }
            }
        }
        getline(file, buffer);
    }

    file.close();
    return 0;
}

double FileInput::getDoubleValue(string variable)
{
    return atof(((string)input[variable]).c_str());
}

int FileInput::getIntValue(string variable)
{
    return (int)(atof(((string)input[variable]).c_str()));
}

bool FileInput::getBooleanValue(string variable)
{
    if (((string)input[variable]).compare("true") == 0)
        return true;
    else if (((string)input[variable]).compare("false") == 0)
        return false;

    return false;
}

string FileInput::getStringValue(string variable)
{
    return input[variable];
}
