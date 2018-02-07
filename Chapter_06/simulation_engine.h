#pragma once
#include "MTwisterFunctions.h"
#include <time.h>
#include <string>
#include <iostream>

class individual
{
public:
	double Allele1trait0[50];
	double Allele2trait0[50];
	double Allele1trait1[50];
	double Allele2trait1[50];
	double Allele1both[50][2];
	double Allele2both[50][2];
	double Genotype[2];
	double Phenotype[2];
	bool Female;
};

class simulation_engine
{
private:
	individual *adult;
	individual *progeny;
	int NadultMax, NprogMax;
	int NumberOfGenerations;
	int PopulationSize;
	double MutationalCorrelation, SelectionalCorrelation;
	int NumLociTrait0, NumLociTrait1; // The number of non-pleiotropic loci for each trait
	int NumLociBoth; // The number of pleiotropic loci (affecting both traits)
	int Fecundity; // Number of offspring each female can produce
	double MutationRatePerLocus; // The per-locus mutation rate
	int CarryingCapacity; // The maximum adult population size
	double MutationalVariance[2];
	double SelectionStrength[2];

public:
	simulation_engine() // Constructor: Initialize parameters and variables here
	{
		srand(static_cast<int>(time(NULL)));
		sgenrand(rand());

		// Set the maximum number of adults and progeny in the population
		// Allocate memory for the pointers for the adults and progeny
		NadultMax = 2000;
		NprogMax = 10000;
		adult = new individual[NadultMax];
		progeny = new individual[NprogMax];

		// Initialize the Parameters

		// Demographic Parameters
		NumberOfGenerations = 10;
		PopulationSize = 100;
		CarryingCapacity = PopulationSize;
		Fecundity = 4;

		// Genetic Parameters
		NumLociTrait0 = 10;
		NumLociTrait1 = 10;
		NumLociBoth = 10;

		// Mutational Parameters
		MutationalVariance[0] = 0.05;
		MutationalVariance[1] = 0.05;
		MutationalCorrelation = 0.5;
		MutationRatePerLocus = 0.0002;

		// Selection Parameters
		SelectionStrength[0] = 49;
		SelectionStrength[1] = 49;
		SelectionalCorrelation = 0.2;
	}

	~simulation_engine()
	{
		delete[] adult;
		delete[] progeny;
	}

	void display_parameters()
	{
		std::cout << "\nNumber of Generations:  \t" << NumberOfGenerations;
		std::cout << "\nPopulation Size:        \t" << PopulationSize;
		std::cout << "\nMutational Correlation: \t" << MutationalCorrelation;
		std::cout << "\nSelectional Correlation:\t" << SelectionalCorrelation;
	}
};
