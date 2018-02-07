#include <iostream>
#include "MTwisterFunctions.h"
#include "simulation_engine.h"

int main()
{
	int generations;

	simulation_engine my_sim;
	my_sim.display_parameters();

	bool initialization_success;
	initialization_success = my_sim.initialize_population();
	if (!initialization_success)
	{
		std::cout << "\nSimulation Initialization Failure!\n";
		return 0;
	}

	int ig;
	for (ig = 0; ig < my_sim.getNumberOfInitialGenerations(); ig++)
	{
		my_sim.monogamy();
		my_sim.mutation();
		my_sim.selection(true);
		my_sim.population_regulation();
	}

	// Turn the lifecycle halfway, so the
	// experimental generations start with progeny
	my_sim.calculate_values_progeny();
	my_sim.monogamy();
	my_sim.mutation();

	for (generations = 0; generations < my_sim.getNumberOfGenerations(); generations++)
	{
		if (!my_sim.is_extinct())
		{
			my_sim.selection(false);
			my_sim.calculate_values_progeny();
			my_sim.population_regulation();
			my_sim.monogamy();
			my_sim.mutation();
			my_sim.calculate_values_adults();
			my_sim.save_population_variables(generations);
		}

	}

	char end_it;
	std::cout << "\n\nEnter any character to exit...";
	std::cin >> end_it;
	return 0;
}

