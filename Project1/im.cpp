#include "im.h"
#include<vector>
#include<array>
#include<iostream>
#include<random>
#include<utility>

using std::vector, std::array;

/* 
This project will implement a Monte-Carlo-type simulator for the 2D Ising model of magnetisation.
We thus need to numerically implement a class containing Ising magnets at each grid site,
perhaps this is simplest by implementing first a class called Spin that we can then place at 
each site? Or is this computationally inefficient? 
The only thing that class would contain would be a value +/- 1, so it's probably superfluous.
Go straight to class Ising.

*/

template<std::size_t L>
class MClattice{
    private:
        array<array<double, L>, L> lat;

    public:
        double& operator()(std::size_t i, std::size_t j) {
            return lat[i][j];
        }

        const double& operator()(std::size_t i, std::size_t j) const {
            return lat[i][j];
        }

        template<std::size_t M>
        friend std::ostream& operator<<(std::ostream& os, const MClattice<M>& lat);
};

template<std::size_t L>
std::ostream& operator<<(std::ostream& os, const MClattice<L>& lat) {
    for (std::size_t i = 0; i < L; ++i) {
        for (std::size_t j = 0; j < L; ++j) {
            os << lat(i,j) << " ";
        }
        os << "\n";
    }
    return os;
}

template<std::size_t L>
class Ising {
    private:
        size_t N;
        double J;
        double H;
        double T;
        double mu = 1;
        double kB = 1;
        MClattice<L> grid;
        array<array<double, 2>, 5> boltzmannFactors; // for our model, we can precompute Boltzmann factors 
                                           // and save much comp time in exponential functions
        
        
        std::mt19937 rng{42};
        std::uniform_int_distribution<int> U_disc{0, L-1};
        std::uniform_real_distribution<double> U{0, 1};

    public:
        Ising(double J, double H, double T) : J{J}, H{H}, T{T} 
        {
            N = L*L;
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    double rv = U(rng);
                    if (rv < 0.5) {
                        grid(i,j) = -1;
                    } else {
                        grid(i,j) = 1; 
                    }
                }
            }
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 2; ++j) {
                    // CALCULATE BOLTZMANN FACTORS HERE
                    int K = i-2;
                    int m = 2*j-1;
                    int k = 2*K;
                    double deltaE = 2*(m*k*J + m*H);
                    boltzmannFactors[i][j] = std::exp(-deltaE/(kB*T)); 
                }
            }
        }
        template<std::size_t M> // overload << operator to print Ising lattice
        friend std::ostream& operator<<(std::ostream& os, const Ising<M>& lattice);

        double getTotalEnergy() {
            /* Naive implementation for initial setup */
            double fieldEnergy = 0;
            double nnEnergy = 0;
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    fieldEnergy -= mu*H*grid(i,j); 
                    nnEnergy -= J*(grid((i-1 + L)%L, j)*grid(i,j) + grid((i+1)%L, j)*grid(i,j) 
                                + grid(i, (j-1 + L)%L)*grid(i,j) + grid(i, (j+1)%L)*grid(i,j)); 
                                /* Correct? Inefficient at least. */
                }
            }
            return fieldEnergy + nnEnergy;
        }

        std::pair<int, int> getLocalEnergyDifference(int i, int j) {
            /* Function to calculate local energy change when flipping 1 spin */ 
            /*double deltaE = 2*grid(i,j)*(J*(grid((i-1 + L)%L, j) + grid((i+1)%L, j)) 
                                    + grid(i, (j-1+L)%L) + grid(i, (j+1)%L)
                                    + H);*/
            int K = grid((i-1 + L)%L, j) + grid((i+1)%L, j) + grid(i, (j-1+L)%L) + grid(i, (j+1)%L);
            int m = grid(i,j);
            return {K, m};
        }

        double getMagnetisation() {
            double magnetisation = 0;
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    magnetisation += grid(i,j); 
                }
            }
            return magnetisation / N;
        }

        double calcBoltzmannFactor(double E) {
            double bF = std::exp(-E/(kB*T));
            return bF;
        }

        double getBoltzmannFactor(int K, int m) {
            /* Function to extract the correct Boltzmann factor from array given 
            a certain value of K = {-4, -2, 0, 2, 4} and m = {-1, 1}. */
            double bF = boltzmannFactors[(K+4)/2][(m+1)/2];
            return bF;
        }

        void metropolisStep() {
            int x = U_disc(rng);
            int y = U_disc(rng);
            
            grid(x,y) *= -1; // spin has been flipped, meaning it has to be flipped back if not accepted
            auto deltaE = getLocalEnergyDifference(x, y); // verify that this has correct sign (calculate)
            int K = deltaE.first;
            int m = deltaE.second;
            double bF = getBoltzmannFactor(K, m);
            if (U(rng) < bF) {
                return; // if flip is accepted, do nothing
            } else { grid(x,y) *= -1; } // if the flip is not accepted, flip back
        }
};

template<std::size_t L> //overload << operator for printing 
std::ostream& operator<<(std::ostream& os, const Ising<L>& lattice) {
    os << lattice.grid;
    return os;
} 

template<std::size_t L>
class Simulation {
    /* A class to implement all simulator-based functions, i.e. the class is responsible for 
    performing simulations given a set of instructions and a certain Ising input. */
    private:
        Ising<L> syst; 
        size_t N = L*L;
        int sweeps = 0;

        std::mt19937 rng{42};
        std::uniform_int_distribution<int> U_disc{0, L-1};
        std::uniform_real_distribution<double> U{0, 1};

    public:
        Simulation(Ising<L> system) : syst(system) {

        }

        void metropolisSweep() {
            /* Function to implement a Metropolis algorithm sweep, that is (strictly speaking)
                a series of spins flips such that all spins have eventually been flipped. 
                However, we simplify and say that one sweep is N attempted flips (to avoid
                having to check for the condition that at least one has been successful per site). */
            for (int i = 0; i < N; ++i) {
                syst.metropolisStep();
            }

            sweeps++;
        }
        void printLattice() {
            std::cout << syst << std::endl;
        }


};










int main() {
    // std::cout << (0-1 + 40)%40 << std::endl;
    // std::cout << (0+1 + 40)%40 << std::endl;


    Ising<10> printTest(1, 1, 1); 
    std::cout << printTest <<std::endl;
    std::cout << printTest.getMagnetisation() << std::endl;
    Simulation<10> testSim(printTest);
    // testSim.testFunc();
    testSim.metropolisSweep();
    testSim.testFunc();
    // std::cout << printTest.getMagnetisation() << std::endl;

    return 0;
}