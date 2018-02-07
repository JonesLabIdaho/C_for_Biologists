#include <iostream>
#include "MTwisterFunctions.h"
#include "simulation_engine.h"

int main()
{
	simulation_engine my_sim;
	my_sim.display_parameters();

	bool initialization_success;
	initialization_success = my_sim.initialize_population();
	if (!initialization_success)
	{
		std::cout << "\nSimulation Initialization Failure!\n";
		return 0;
	}

	my_sim.output_adults();

	char end_it;
	std::cout << "\n\nEnter any character to exit...";
	std::cin >> end_it;
	return 0;
}

