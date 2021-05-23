#include "simu.h"


int main()
{
	//std::vector<double> qs = {0.00,0.05, 0.1,7.0/30 - 0.01,0.3};
	std::vector<double> qs = {0.00,0.05, 0.1, 0.02 ,0.3};
	//simulationAverages(10000, 6, 100, 200, 0,1.0, 40, qs,"Results\\Triangular\\triangular");
	//simulationAverages(10000, 4, 80, 100, 0,0.5, 20, qs,"Results\\Square\\square");
	simulationAveragesNonInfty(25, 25, 800, 500, 0.005,0.35, 30, qs,"triangular", "Results\\Non_infty\\Triangular\\triangular");
	simulationAveragesNonInfty(25, 25, 800, 500, 0.005,0.35, 30, qs,"square","Results\\Non_infty\\Square\\square");
	return 0;
}