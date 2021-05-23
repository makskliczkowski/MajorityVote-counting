#pragma once
#include "General.h"

#include  <vector>
#include <random>
#include <fstream>
#include <numeric>
#include <string>
#include <iostream>
#include <execution>
using namespace std;

inline std::random_device rd{};
inline std::mt19937 eng{rd()};

/*struct random_num {
    // Hold RNG state as a member variable
    std::mt19937 eng;
    random_num()
	{
        eng = std::mt19937{static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count()) };
	}

    double randZero_One_Uni()
    {
        return std::uniform_real_distribution<double>(0, 1)(this->eng);
    }
    int randInt_Uni(int L)
    {
        return std::uniform_int_distribution<>(1, L)(eng);
    }
};*/


double randZero_One_Uni();
int randInt_Uni(int L);
int myModuloEuclidean(int a, int b);

double proba(int nei_sum, double p, double q, int spin);
void simulationAverages(int N, int num_nei, int mcs, int av_num, double pstart,double pmax, int pnum, std::vector<double> qs,std::string type = "square", bool plot_in_time = false);
void simulationAveragesNonInfty(int Lx,int Ly, int mcs, int av_num, double pstart,double pmax, int pnum, std::vector<double> qs,std::string type = "square",std::string savename="simu");
std::vector<short> make_ran_lattice(int N);

