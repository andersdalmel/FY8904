#include "im.h"
#include<vector>
#include<array>
#include<iostream>
#include<random>
#include<utility>
#include<fstream>

using std::vector, std::array;

/* 
This project will implement a Monte-Carlo-type simulator for the 2D Ising model of magnetisation.
We thus need to numerically implement a class containing Ising magnets at each grid site,
perhaps this is simplest by implementing first a class called Spin that we can then place at 
each site? Or is this computationally inefficient? 
The only thing that class would contain would be a value +/- 1, so it's probably superfluous.
Go straight to class Ising.
 
Maybe it would be cool to simulate enough samples to approximate M(T)? 
*/

template<std::size_t L>
class MClattice{
    private:
        array<array<int8_t, L>, L> lat;

    public:
        // Allow indexing MClattice
        int8_t& operator()(std::size_t i, std::size_t j) {
            return lat[i][j];
        }

        const int8_t& operator()(std::size_t i, std::size_t j) const {
            return lat[i][j];
        }

        // Make MClattice iterable
        auto begin() { return lat.begin(); }
        auto end()   { return lat.end(); }

        auto begin() const { return lat.begin(); }
        auto end()   const { return lat.end(); }

        // Overload << for printing
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
        double NM;

        int flipSuccess;

        std::mt19937 rng;
        std::uniform_int_distribution<int> U_disc{0, L-1};
        std::uniform_real_distribution<double> U{0, 1};

    public:
        Ising(double J, double H, double T, unsigned int seed) : J{J}, H{H}, T{T}, rng{seed} {
            
            N = L*L;
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    double rv = U(rng);
                    if (rv < 0.5) {
                        grid(i,j) = -1;
                        NM -= 1;
                    } else {
                        grid(i,j) = 1; 
                        NM += 1;
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



        /* ------------------- */
        /* --- calculators --- */
        /* ------------------- */

        double calculateTotalEnergy() {
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

        std::pair<int, int> calculateLocalEnergyDifference(int i, int j) {
            /* Function to calculate local energy change when flipping 1 spin */ 
            /*double deltaE = 2*grid(i,j)*(J*(grid((i-1 + L)%L, j) + grid((i+1)%L, j)) 
                                    + grid(i, (j-1+L)%L) + grid(i, (j+1)%L)
                                    + H);*/
            int K = grid((i-1 + L)%L, j) + grid((i+1)%L, j) + grid(i, (j-1+L)%L) + grid(i, (j+1)%L);
            int m = grid(i,j);
            return {K, m};
        }

        double calculateMagnetisation() {
            double magnetisation = 0;
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    magnetisation += grid(i,j); 
                }
            }
            return magnetisation / N;
        }

        double calculateBoltzmannFactor(double E) {
            double bF = std::exp(-E/(kB*T));
            return bF;
        }
                
        /* ------------------- */
        /* getters and setters */
        /* ------------------- */
        
        double getBoltzmannFactor(int K, int m) {
            /* Function to extract the correct Boltzmann factor from array given 
            a certain value of K = {-4, -2, 0, 2, 4} and m = {-1, 1}. */
            double bF = boltzmannFactors[(K+4)/2][(m+1)/2];
            return bF;
        }


        double getNM() {
            return NM;
        }

        int getFlipSuccess() {
            return flipSuccess;
        }


                
        /* ------------------- */
        /* Monte Carlo methods */
        /* ------------------- */

        void metropolisStep() {
            int x = U_disc(rng);
            int y = U_disc(rng);
            
            auto deltaE = calculateLocalEnergyDifference(x, y); // verify that this has correct sign (calculate)
            int K = deltaE.first;
            int m = deltaE.second;
            double bF = getBoltzmannFactor(K, m);
            if (U(rng) < bF) {
                grid(x,y) *= -1;
                NM += 2*grid(x,y);
                ++flipSuccess;
                return; // if flip is accepted, perform flip, update total magnetisation, note success.
            } 
        }

        void saveLatticeSnapshot(const std::string& filename, bool append=true) {
            std::ofstream out(filename, append ? std::ios::binary | std::ios::app : std::ios::binary);
            if (!out) {
                throw std::runtime_error("Cannot open file " + filename);
            }
            
            // Write row by row
            for (const auto& row : grid) {
                out.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(int8_t));
            }

            out.close();
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
        vector<double> magnetisation;
        
        double totalEnergy;
        double finalMagnetisation;

        std::mt19937 rng;
        std::uniform_int_distribution<int> U_disc{0, L-1};
        std::uniform_real_distribution<double> U{0, 1}; 

    public:
        Simulation(Ising<L> system, unsigned int seed) : syst(system), rng{seed} { }

        void metropolisSweep() {
            /* Function to implement a Metropolis algorithm sweep, i.e. (strictly speaking)
                a series of spins flips such that all spins have eventually been flipped. 
                However, we simplify and say that one sweep is N attempted flips (to avoid
                having to check for the condition that at least one has been successful per site). */
            for (int i = 0; i < N; ++i) {
                syst.metropolisStep();
            }
            sweeps++;
        }

        vector<double> getMagnetisation() {
            return magnetisation;
        }

        double getTotalEnergy() {
            return totalEnergy;
        }

        double getFinalMagnetisation() {
            return finalMagnetisation;
        }

        void printLattice() {
            std::cout << syst << std::endl;
        }

        void saveMagnetisation(const std::string& filename) {
            std::ofstream file(filename);
            if (!file) {
                throw std::runtime_error("Cannot open file");
            }
            file << "time,magnetisation\n";
            for (int i = 0; i < magnetisation.size(); ++i) {
                file << i*100 << "," << magnetisation[i] << "\n";
            }
        } 

        void simulation(int n, int saveRate = 0) {
            /* Function to run a MC simulation of n sweeps. */
            if (saveRate > 0) {
                for (int i = 0; i < n; ++i) {
                    metropolisSweep(); // only save every 100th value to reduce correlation & memory usage 
                    if (i % 100 == 0) { // this should also implement a equilibrium assertion
                        double NM = syst.getNM();
                        magnetisation.push_back(NM/N);
                        std::cout << "At sweep no. " << sweeps << std::endl;
                    }
                    if (i % saveRate == 0) {
                        syst.saveLatticeSnapshot("latticeSnaps.bin", true);
                    }
                }
                totalEnergy = syst.calculateTotalEnergy();
                int flipSuccesses = syst.getFlipSuccess();
                double successPercentage = static_cast<double>(flipSuccesses)/(n*N);
                std::cout << "DONE! \nPercentage of spin flips successful: " << successPercentage << std::endl;
            } else {
                for (int i = 0; i < n; ++i) {
                    metropolisSweep(); // only save every 100th value to reduce correlation & memory usage 
                    if (i % 100 == 0) { // this should also implement a equilibrium assertion
                        double NM = syst.getNM();
                        magnetisation.push_back(NM/N);
                        std::cout << "At sweep no. " << sweeps << std::endl;
                    }
                }
                totalEnergy = syst.calculateTotalEnergy();
                finalMagnetisation = syst.getNM() / N;
                int flipSuccesses = syst.getFlipSuccess();
                double successPercentage = static_cast<double>(flipSuccesses)/(n*N);
                std::cout << "DONE! \nPercentage of spin flips successful: " << successPercentage << std::endl;
            }
        }
        

};



template<std::size_t L> 
class Ensemble {
    private:
        int S;
        double J;
        double H;
        vector<double> T;
        unsigned int seed;
        int simLength;
        double energy;          // running average of system energy at the end of simulations
        double magnetisation;   // running average of system magnetisation at the end of simulations

    public:
        Ensemble(int samples, double J, double H, vector<double> T, unsigned int seed, int simLength) : S{samples}, J{J}, H{H}, T{T}, 
                                            seed{seed}, simLength{simLength} { }

        void runSimulations() {
            for (double temp : T) {    
                for (int s = 0; s < S; ++s) {
                    Ising<L> syst(J, H, T, seed + i);
                    Simulation<L> sim(syst, seed + i);
                    sim.simulation(simLength, 0); // change saverate of snaps manually if snaps are desired
                    energy = energy*(i - 1)/i + sim.getTotalEnergy() / i; // calc. running avg
                    magnetisation = magnetisation*(i - 1)/i + sim.getFinalMagnetisation() / i; // calc. running avg
                }
            }
        }



};





int main() {
    // ISING(J, H, T)
    Ising<40> testIsing(1, 0, 2.1, 42); 
    Simulation<40> testSim(testIsing, 42);
    testSim.simulation(10000, 1);
    testSim.saveMagnetisation("testMagnetisation.csv");
   
    return 0;
}