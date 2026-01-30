#include "im.h"
#include<vector>


using std::vector;

/* 
This project will implement a Monte-Carlo-type simulator for the 2D Ising model of magnetisation.
We thus need to numerically implement a class containing Ising magnets at each grid site,
perhaps this is simplest by implementing first a class called Spin that we can then place at 
each site? Or is this computationally inefficient? 
The only thing that class would contain would be a value +/- 1, so it's probably superfluous.
Go straight to class Ising.

*/


class Ising {
    private:
        int L;
        int N;
        double J;
        double H;
        double T;
        vector<vector<double>> grid;

    public:
        Ising(int L, double J, double H, double T) : L{L}, J{J}, H{H}, T{T} 
        {
            N = L*L;

        }
        
};