#pragma once
#include <iostream>
#include "MTwisterFunctions.h"
#include <time.h>
#include <string>

class individual
{
public:
	std::string name;
	int age;
	double weight;
	double height;
};


class simulation_engine
{
private:
	individual mom;
	individual dad;
	individual kid[3];

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

	void set_dad(std::string dad_name, int dad_age, double dad_weight, double dad_height)
	{
		dad.name = dad_name;
		dad.age = dad_age;
		dad.weight = dad_weight;
		dad.height = dad_height;
	}

	void set_mom(std::string mom_name, int mom_age, double mom_weight, double mom_height)
	{
		mom.name = mom_name;
		mom.age = mom_age;
		mom.weight = mom_weight;
		mom.height = mom_height;
	}

	void set_individual(individual &ind, std::string ind_name, int ind_age,
		double ind_weight, double ind_height)
	{
		ind.name = ind_name;
		ind.age = ind_age;
		ind.weight = ind_weight;
		ind.height = ind_height;
	}

	void output_mom()
	{
		std::cout << "\nMom:";
		std::cout << "\nName:\t" << mom.name;
		std::cout << "\nAge:\t" << mom.age;
		std::cout << "\nHeight:\t" << mom.height;
		std::cout << "\nWeight:\t" << mom.weight;
	}

	void output_dad()
	{
		std::cout << "\nDad:";
		std::cout << "\nName:\t" << dad.name;
		std::cout << "\nAge:\t" << dad.age;
		std::cout << "\nHeight:\t" << dad.height;
		std::cout << "\nWeight:\t" << dad.weight;
	}

	void output_individual(individual &ind)
	{
		std::cout << "\nName:\t" << ind.name;
		std::cout << "\nAge:\t" << ind.age;
		std::cout << "\nHeight:\t" << ind.height;
		std::cout << "\nWeight:\t" << ind.weight;
	}

	void set_kid(int whichkid, std::string kid_name, int kid_age, double kid_weight,
		double kid_height)
	{
		set_individual(kid[whichkid], kid_name, kid_age, kid_weight, kid_height);
	}

	void output_kid(int whichkid)
	{
		std::cout << "\nKid " << whichkid + 1;
		output_individual(kid[whichkid]);
	}



};
