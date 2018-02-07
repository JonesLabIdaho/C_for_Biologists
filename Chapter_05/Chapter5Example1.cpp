#include <iostream>
#include "simulation_engine.h"

int main()
{
	simulation_engine my_sim;

	my_sim.set_dad("Fred", 44, 210, 69);
	my_sim.set_mom("Wilma", 44, 150, 64);
	my_sim.set_kid(0, "Joe", 4, 35, 40);
	my_sim.set_kid(1, "Marge", 6, 44, 46);
	my_sim.set_kid(2, "Kelly", 8, 65, 50);

	my_sim.output_dad();
	std::cout << "\n";

	my_sim.output_mom();
	std::cout << "\n";

	for (size_t i = 0; i < 3; i++)
	{
		my_sim.output_kid(i);
		std::cout << "\n";
	}

	char end_it;
	std::cout << "\n\nEnter any character to exit...";
	std::cin >> end_it;

	return 0;
}
