#include <iostream>
#include "MTwisterFunctions.h"
#include "simulation_engine.h"

int main()
{
	int i;

	simulation_engine my_sim;
	my_sim.display_parameters();

	std::cout << "\n\nRAND_MAX = " << RAND_MAX << "\n";

	double drnd1, drnd2;
	std::cout << "\nUniformly Distributed:\t";
	drnd1 = genrand();
	std::cout << drnd1 << "\n";

	std::cout << "\nBivariate Normally Distributed:\t";
	randbivnorm(5, 2, 0.8, drnd1, drnd2);
	std::cout << drnd1 << "\t" << drnd2 << "\n";

	std::cin >> i;
	return 0;
}


