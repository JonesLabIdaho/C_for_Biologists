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

	std::cout << "\nAdults:";
	my_sim.output_adults();
	my_sim.polygynous_mating();
	my_sim.mutation();
	my_sim.natural_selection();
	std::cout << "\n\nProgeny:";
	my_sim.output_progeny();

	char end_it;
	std::cout << "\n\nEnter any character to exit...";
	std::cin >> end_it;
	return 0;
}

