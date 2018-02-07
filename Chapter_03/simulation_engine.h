#pragma once
#include <iostream>
#include "MTwisterFunctions.h"
#include <time.h>

class simulation_engine
{
private:
	int NumberOfGenerations;
	int PopulationSize;
	double MutationalCorrelation, SelectionalCorrelation;

public:
	simulation_engine() // Constructor: Initialize parameters and variables here
	{
		srand(static_cast<int>(time(NULL)));
		sgenrand(rand());

		NumberOfGenerations = 10;
		PopulationSize = 100;
		MutationalCorrelation = 0.5;
		SelectionalCorrelation = 0.2;
	}

	void display_parameters()
	{
		std::cout << "\nNumber of Generations:  \t" << NumberOfGenerations;
		std::cout << "\nPopulation Size:        \t" << PopulationSize;
		std::cout << "\nMutational Correlation: \t" << MutationalCorrelation;
		std::cout << "\nSelectional Correlation:\t" << SelectionalCorrelation;
	}

};
