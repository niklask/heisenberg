using namespace std;

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include <math.h>

#include "fileinput.h"
#include "statevector.h"
#include "operatorlist.h"
#include "random.h"

int no_loopMoves;

string int2str(int i)
{
    stringstream stream;
    stream << i;
    return stream.str();
}

void mcStep(StateVector& sv, OperatorList& ol, double beta)
{
    ol.diagonalUpdate(sv, beta);

    if (ol.check(sv) != true)
        exit(0);

    for (int i = 0; i < no_loopMoves; i++) {
        ol.loopUpdate(sv, beta);

        if (ol.check(sv) != true)
            exit(0);
    }

}

double simulation(int no_spins, int spin_conf, int warmup_steps, int no_mcSteps,
                  int list_size, double beta)
{
    StateVector* sv;

    // Create the spin configuration...
    if (spin_conf == 0) {
        // if spin_conf is set to zero we generate a random spin configuration
        sv = new StateVector(no_spins, true);
    } else {
        sv = new StateVector(no_spins, spin_conf);
    }

    // and an operator list with only unit operators
    OperatorList ol(list_size);

    // some warmup steps
    if (warmup_steps > 0) {
        for (int i = 0; i < warmup_steps; i++) {
            mcStep(*sv, ol, beta);
        }
    }

    double energy = 0;
    if (no_mcSteps > 0) {
        // Do MC steps and calculate energy after each step
        for (int i = 0; i < no_mcSteps; i++) {
            mcStep(*sv, ol, beta);
            int n = ol.countDOD();
            energy += (double)(-n/(beta*no_spins) + 0.25);
        }
    }

    delete sv;
    return (double)(energy/no_mcSteps);
}

int main(int argc, char* argv[])
{
    // load the configuration file
    FileInput input("qmc.input");

    // read the parameters
    double start_temp = input.getDoubleValue("start_temp");
    double end_temp = input.getDoubleValue("end_temp");
    double temp_step = input.getDoubleValue("temp_step");

    int list_size = input.getIntValue("list_size");
    int no_spins = input.getIntValue("no_spins");
    int spin_conf = input.getIntValue("spin_conf");

    no_loopMoves = input.getIntValue("no_loopMoves");

    int warmup_steps = input.getIntValue("warmup_steps");
    int no_mcSteps = input.getIntValue("no_mcSteps");

    // seed the random number generator
    long long int seed = -(long long int)input.getDoubleValue("seed");
    seed_ran(seed);

    cout << setprecision(10);
    for (double t = start_temp; t <= end_temp; t += temp_step) {
        // calculate the inverse temperature
        double beta = 1/t;

        if (t == 1 || t == 2 || t == 3 || t == 4)
            no_mcSteps *= 10;
        // perform one simulation and echo the result
        double energy = simulation(no_spins, spin_conf, warmup_steps, 
                                   no_mcSteps, list_size, beta);
        cout << t << "\t" << energy << endl;
    }
}
