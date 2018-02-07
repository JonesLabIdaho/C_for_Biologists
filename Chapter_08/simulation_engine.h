#pragma once
#include "MTwisterFunctions.h"
#include <time.h>
#include <string>
#include <iostream>
#include <iomanip>

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

	void calculate_genotypic_values(int n_loci_0, int n_loci_1, int n_loci_both)
	{
		int i;
		Genotype[0] = 0;
		for (i = 0; i < n_loci_0; i++)
			Genotype[0] = Genotype[0] + Allele1trait0[i] + Allele2trait0[i];
		for (i = 0; i < n_loci_both; i++)
			Genotype[0] = Genotype[0] + Allele1both[i][0] + Allele2both[i][0];

		Genotype[1] = 0;
		for (i = 0; i < n_loci_1; i++)
			Genotype[1] = Genotype[1] + Allele1trait1[i] + Allele2trait1[i];
		for (i = 0; i < n_loci_both; i++)
			Genotype[1] = Genotype[1] + Allele1both[i][1] + Allele2both[i][1];
	}

	void calculate_phenotype(double env_st_dev_0, double env_st_dev_1)
	{
		Phenotype[0] = Genotype[0] + randnorm(0, env_st_dev_0);
		Phenotype[1] = Genotype[1] + randnorm(0, env_st_dev_1);
	}

	void set_sex()
	{
		if (genrand() < 0.5)
			Female = true;
		else
			Female = false;
	}

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
	double EnvironmentalVariance[2];
	double EnvironmentalStDev[2];
	int MaxMatingEncounters;
	int Nprogeny;


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
		PopulationSize = 6;
		CarryingCapacity = PopulationSize;
		Fecundity = 4;
		MaxMatingEncounters = 500;

		// Genetic Parameters
		NumLociTrait0 = 1;
		NumLociTrait1 = 1;
		NumLociBoth = 1;
		EnvironmentalVariance[0] = 1;
		EnvironmentalVariance[1] = 1;
		EnvironmentalStDev[0] = sqrt(EnvironmentalVariance[0]);
		EnvironmentalStDev[1] = sqrt(EnvironmentalVariance[1]);

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
		std::cout << "Parameter_Values:\n";
		// Demographic Parameters
		std::cout << "Demographic_Parameters:\n";
		std::cout << "No_Generations:   \t" << NumberOfGenerations << "\n";
		std::cout << "Initial_Pop_Size: \t" << PopulationSize << "\n";
		std::cout << "Carrying_Capacity:\t" << CarryingCapacity << "\n";
		std::cout << "Female_Fecundity: \t" << Fecundity << "\n";
		std::cout << "Max_Mating_Enc:   \t" << MaxMatingEncounters << "\n";

		// Quantitative Genetic Parameters
		std::cout << "Quantitative_Genetic_Parameters:\n";
		std::cout << "No_Loci_Trait0:   \t" << NumLociTrait0 << "\n";
		std::cout << "No_Loci_Trait1:   \t" << NumLociTrait1 << "\n";
		std::cout << "No_Loci_Pleiotrop:\t" << NumLociBoth << "\n";
		std::cout << "Env_Variance_Trt0:\t" << EnvironmentalVariance[0] << "\n";
		std::cout << "Env_Variance_Trt1:\t" << EnvironmentalVariance[1] << "\n";


		// Mutational Parameters
		std::cout << "Mutational_Parameters:\n";
		std::cout << "Mut_Var_Trait0:   \t" << MutationalVariance[0] << "\n";
		std::cout << "Mut_Var_Trait1:   \t" << MutationalVariance[1] << "\n";
		std::cout << "Mut_Correlation:  \t" << MutationalCorrelation << "\n";
		std::cout << "Mutation_Rate:    \t" << MutationRatePerLocus << "\n";

		// Selection Parameters
		std::cout << "Selection_Parameters:\n";
		std::cout << "Omega_Trait0:     \t" << SelectionStrength[0] << "\n";
		std::cout << "Omega_Trait1:     \t" << SelectionStrength[1] << "\n";
		std::cout << "Selection_Corr:   \t" << SelectionalCorrelation << "\n";
	}

	bool initialize_population()
	{
		bool good_initialization = true;
		int i, j;
		double dIP1, allelic_std_dev_0, allelic_std_dev_1;

		// Make sure the population size is not larger than the size of the adult array
		if (PopulationSize > NadultMax)
		{
			PopulationSize = NadultMax;
			good_initialization = false;
		}

		// Set the starting allelic values for the adults in the population
		dIP1 = NumLociTrait0 + NumLociBoth;
		if (dIP1 > 0)
			allelic_std_dev_0 = sqrt(MutationalVariance[0]) / dIP1;
		else
			allelic_std_dev_0 = 0;

		dIP1 = NumLociTrait1 + NumLociBoth;
		if (dIP1 > 0)
			allelic_std_dev_1 = sqrt(MutationalVariance[1]) / dIP1;
		else
			allelic_std_dev_1 = 0;

		for (i = 0; i < PopulationSize; i++)
		{
			// Set allelic values to zero for loci affecting trait 0
			for (j = 0; j < NumLociTrait0; j++)
			{
				adult[i].Allele1trait0[j] = randnorm(0, allelic_std_dev_0);
				adult[i].Allele2trait0[j] = randnorm(0, allelic_std_dev_0);
			} // end of j

			// Set allelic values to zero for loci affecting trait 1
			for (j = 0; j < NumLociTrait1; j++)
			{
				adult[i].Allele1trait1[j] = randnorm(0, allelic_std_dev_1);
				adult[i].Allele2trait1[j] = randnorm(0, allelic_std_dev_1);
			} // end of j

			  // Set allelic values to zero for pleiotropic loci
			for (j = 0; j < NumLociBoth; j++)
			{
				adult[i].Allele1both[j][0] = randnorm(0, allelic_std_dev_0);
				adult[i].Allele1both[j][1] = randnorm(0, allelic_std_dev_1);
				adult[i].Allele2both[j][0] = randnorm(0, allelic_std_dev_0);
				adult[i].Allele2both[j][1] = randnorm(0, allelic_std_dev_1);
			} // end of j

			adult[i].calculate_genotypic_values(NumLociTrait0, NumLociTrait1, NumLociBoth);
			adult[i].calculate_phenotype(EnvironmentalStDev[0], EnvironmentalStDev[1]);
			adult[i].set_sex();
		} // end of i

		return good_initialization;
	}

	void output_adults()
	{
		std::cout << "\n\nID\tSex\tTrait\tGeno\tPheno\t";
		int i, j, k;
		for (j = 0; j < NumLociTrait0; j++)
			std::cout << "Locus_" << j << "Trait0\t";
		for (j = 0; j < NumLociTrait1; j++)
			std::cout << "Locus_" << j << "Trait1\t";
		for (j = 0; j < NumLociBoth; j++)
			std::cout << "Locus_" << j << "Pleio\t";

		// Output two rows for each adult
		// The first row will be for trait 0
		// The second row will be for trait 1

		for (i = 0; i < PopulationSize; i++)
		{
			for (k = 0; k < 2; k++)
			{
				std::cout << "\n" << i << "\t";
				if (adult[i].Female)
					std::cout << "female\t";
				else
					std::cout << "male\t";
				std::cout << k << "\t";
				std::cout << std::setprecision(3) << std::fixed << adult[i].Genotype[k] << "\t";
				std::cout << std::setprecision(3) << std::fixed << adult[i].Phenotype[k] << "\t";
				if (k == 0)
				{
					for (j = 0; j < NumLociTrait0; j++)
						std::cout << std::setprecision(3) << std::fixed << adult[i].Allele1trait0[j]
						<< "/" << adult[i].Allele2trait0[j] << "\t";
					for (j = 0; j < NumLociTrait1; j++)
						std::cout << "         \t";
				} // end of if (k == 0)
				else
				{
					for (j = 0; j < NumLociTrait0; j++)
						std::cout << "         \t";
					for (j = 0; j < NumLociTrait1; j++)
						std::cout << std::setprecision(3) << std::fixed << adult[i].Allele1trait1[j]
						<< "/" << adult[i].Allele2trait1[j] << "\t";
				} // end of else
				for (j = 0; j < NumLociBoth; j++)
				{
					std::cout << std::setprecision(3) << std::fixed << adult[i].Allele1both[j][k]
						<< "/" << adult[i].Allele2both[j][k] << "\t";
				}
			} // end of k loop
		} // end of i loop

	}

	void polygynous_mating()
	{
		// This function implements strict polygyny.
		// Under this mating system, each female mates once, but
		// each male can mate an unlimited number of times.
		// Females choose males at random.

		int i, j, m;
		int iPC;
		bool mate_found;
		int mateID, counter, rnum;
		double dRnum;

		// Check to make sure at least one male is present in the population
		bool males_present = false;
		for (i = 0; i < PopulationSize; i++)
		{
			if (!adult[i].Female)
				males_present = true;
		}

		iPC = 0;
		for (i = 0; i < PopulationSize; i++)
		{
			if (adult[i].Female && males_present)
			{

				// Find a mate for this female
				// Mating is random, so any male will do

				mate_found = false;
				counter = 0;
				while (counter < MaxMatingEncounters && !mate_found)
				{
					rnum = randnum(PopulationSize);
					if (!adult[rnum].Female)
					{
						mateID = rnum;
						mate_found = true;
					}
					counter++;
				} // end of while

				// If a mate is found, produce progeny
				if (mate_found)
				{
					for (m = 0; m < Fecundity; m++)
					{
						if (iPC >= NprogMax)
							iPC = NprogMax - 1;

						// First, let's take care of Mendelian assortment of the trait0 loci
						for (j = 0; j < NumLociTrait0; j++)
						{
							// The progeny needs one maternal allele and one paternal allele.
							// The maternal allele will be from adult[i] (the mother), and
							// we determine which allele with basically a coinflip. The
							// function genrand() produces a value in the range [0,1), so
							// if this number is less than 0.5 we choose one allele. Otherwise,
							// we choose the other allele.

							dRnum = genrand();
							if (dRnum < 0.5)
								progeny[iPC].Allele1trait0[j] = adult[i].Allele1trait0[j];
							else
								progeny[iPC].Allele1trait0[j] = adult[i].Allele2trait0[j];

							// The procedure is the same for the father, adult[mateID].
							dRnum = genrand();
							if (dRnum < 0.5)
								progeny[iPC].Allele2trait0[j] = adult[mateID].Allele1trait0[j];
							else
								progeny[iPC].Allele2trait0[j] = adult[mateID].Allele2trait0[j];

						} // end of j loop

						// Second, take care of Mendelian assortment for the trait1 loci.
						// The procedure is essentially identical to that for trait0.
						for (j = 0; j < NumLociTrait1; j++)
						{
							dRnum = genrand();
							if (dRnum < 0.5)
								progeny[iPC].Allele1trait1[j] = adult[i].Allele1trait1[j];
							else
								progeny[iPC].Allele1trait1[j] = adult[i].Allele2trait1[j];

							dRnum = genrand();
							if (dRnum < 0.5)
								progeny[iPC].Allele2trait1[j] = adult[mateID].Allele1trait1[j];
							else
								progeny[iPC].Allele2trait1[j] = adult[mateID].Allele2trait1[j];
						} // end of j loop

						// The pleiotropic loci are more complicated. We have to be sure to
						// keep the allelic effects (on the two traits) together for each allele.
						// Otherwise, the loci would not behave realistically like actual
						// pleiotropic loci.
						for (j = 0; j < NumLociBoth; j++)
						{
							dRnum = genrand();
							if (dRnum < 0.5)
							{
								progeny[iPC].Allele1both[j][0] = adult[i].Allele1both[j][0];
								progeny[iPC].Allele1both[j][1] = adult[i].Allele1both[j][1];
							}
							else
							{
								progeny[iPC].Allele1both[j][0] = adult[i].Allele2both[j][0];
								progeny[iPC].Allele1both[j][1] = adult[i].Allele2both[j][1];
							}
							dRnum = genrand();
							if (dRnum < 0.5)
							{
								progeny[iPC].Allele2both[j][0] = adult[mateID].Allele1both[j][0];
								progeny[iPC].Allele2both[j][1] = adult[mateID].Allele1both[j][1];
							}
							else
							{
								progeny[iPC].Allele2both[j][0] = adult[mateID].Allele2both[j][0];
								progeny[iPC].Allele2both[j][1] = adult[mateID].Allele2both[j][1];
							}
						} // end of j loop

						progeny[iPC].set_sex();
						progeny[iPC].calculate_genotypic_values(NumLociTrait0, NumLociTrait1, NumLociBoth);
						progeny[iPC].calculate_phenotype(EnvironmentalStDev[0], EnvironmentalStDev[1]);
						iPC++;

					} // end of m loop
				} // end of if (mate_found)
			} // end of if (adult[i].Female && males_present)
		} // end of i loop

		Nprogeny = iPC;

	}

	void output_progeny()
	{
		std::cout << "\n\nID\tSex\tTrait\tGeno\tPheno\t";

		int i, j, k;
		for (j = 0; j < NumLociTrait0; j++)
			std::cout << "Locus_" << j << "Trait0\t";
		for (j = 0; j < NumLociTrait1; j++)
			std::cout << "Locus_" << j << "Trait1\t";
		for (j = 0; j < NumLociBoth; j++)
			std::cout << "Locus_" << j << "Pleio\t";

		// Output two rows for each offspring
		// The first row will be for trait 0
		// The second row will be for trait 1

		for (i = 0; i < Nprogeny; i++)
		{
			for (k = 0; k < 2; k++)
			{
				std::cout << "\n" << i << "\t";
				if (progeny[i].Female)
					std::cout << "female\t";
				else
					std::cout << "male\t";
				std::cout << k << "\t";
				std::cout << std::setprecision(3) << std::fixed << progeny[i].Genotype[k] << "\t";
				std::cout << std::setprecision(3) << std::fixed << progeny[i].Phenotype[k] << "\t";
				if (k == 0)
				{
					for (j = 0; j < NumLociTrait0; j++)
						std::cout << std::setprecision(3) << std::fixed << progeny[i].Allele1trait0[j]
						<< "/" << progeny[i].Allele2trait0[j] << "\t";
					for (j = 0; j < NumLociTrait1; j++)
						std::cout << "         \t";
				} // end of if (k == 0)
				else
				{
					for (j = 0; j < NumLociTrait0; j++)
						std::cout << "         \t";
					for (j = 0; j < NumLociTrait1; j++)
						std::cout << std::setprecision(3) << std::fixed << progeny[i].Allele1trait1[j]
						<< "/" << progeny[i].Allele2trait1[j] << "\t";
				} // end of else

				for (j = 0; j < NumLociBoth; j++)
				{
					std::cout << std::setprecision(3) << std::fixed << progeny[i].Allele1both[j][k]
						<< "/" << progeny[i].Allele2both[j][k] << "\t";
				}
			} // end of k loop
		} // end of i loop
	}


};
