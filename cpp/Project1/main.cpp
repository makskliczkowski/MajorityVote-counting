#include "simu.h"


int main()
{
	//std::vector<double> qs = {0.00,0.05, 0.1,7.0/30 - 0.01,0.3};
	std::vector<double> qs = { 0.03,0.01};
	//simulationAverages(10000, 6, 100, 200, 0,1.0, 40, qs,"Results\\Triangular\\triangular");
	//simulationAverages(10000, 4, 80, 100, 0,0.5, 20, qs,"Results\\Square\\square");
	std::vector<int> Ls_square = {70};
	std::vector<int> Ls_triangle = {70};
	for(auto L:Ls_triangle)
	{
		simulationAveragesNonInfty(L, L, 150000, 6000, 0.03,0.32, 160, qs,"triangular", "Results\\Non_infty\\Triangular\\triangular");

	}
	for(auto L:Ls_square)
	{
		simulationAveragesNonInfty(L, L, 150000, 6000, 0.03,0.32, 160, qs,"square","Results\\Non_infty\\Square\\square");
	}
	return 0;
}