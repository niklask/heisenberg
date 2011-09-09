#include <fstream>
#include <string>
#include <sstream>

#include "matrix.h"
#include "hamiltonian.h"
#include "fileinput.h"

string int2str(int i)
{
    stringstream stream;
    stream << i;
    return stream.str();
}
 
//  calculate solution for a given system
void calculate(int size, double start_temp, double end_temp,
               double temp_step, const char* outputfile)
{

    ofstream file(outputfile);

    Hamiltonian H = Hamiltonian(size);
    H.Diagonalize();
    for (double t = start_temp; t <= end_temp; t += temp_step) {
        double beta = (1/t);
        Matrix eH = H;
        eH *= (-beta);
        eH.Exp();
        double energy = (H*eH).Trace() / eH.Trace();
        energy = energy / size;
        file << t << "\t" << energy << endl;
    }
    
    file.close();
}

int main(int argc, char* argv[])
{
    const char* inputfile = "diag.input";
    FileInput input(inputfile);

    string file_prefix = input.getStringValue("file_prefix");
    int no_runs = input.getIntValue("no_runs");
    int Nspin_start = input.getIntValue("Nspin_start");
    int spin_step = input.getIntValue("spin_step");
    double start_temp = input.getDoubleValue("start_temp");
    double end_temp = input.getDoubleValue("end_temp");
    double temp_step = input.getDoubleValue("temp_step");

    for (int i = 0; i < no_runs; i++) {
        int size = i * spin_step + Nspin_start;
        string outputfile = file_prefix + "-N" + int2str(size);

        calculate(size, start_temp, end_temp, temp_step, outputfile.c_str());
    }   
}
