#pragma once
#include "MTwisterFunctions.h"
#include <time.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

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
	bool Alive;

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
	double GaussianPreferenceVariance;
	bool PopulationExtinct;

	// Variables corresponding to population-level summary statistics:
	double phenotypic_mean[2], genotypic_mean[2];
	double phenotypic_variance[2], genotypic_variance[2];
	double phenotypic_covariance, genotypic_covariance;
	double phenotypic_correlation, genotypic_correlation;


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
		PopulationSize = 10;
		CarryingCapacity = PopulationSize;
		Fecundity = 4;
		PopulationExtinct = false;

		// Mating Parameters
		MaxMatingEncounters = 500;
		GaussianPreferenceVariance = 100;

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
		MutationRatePerLocus = 0.05;

		// Selection Parameters
		SelectionStrength[0] = 9;
		SelectionStrength[1] = 9;
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

		// Mating Parameters
		std::cout << "Mating_Parameters:\n";
		std::cout << "Max_Mating_Enc:   \t" << MaxMatingEncounters << "\n";
		std::cout << "Gaussian_Pref_Var:\t" << GaussianPreferenceVariance << "\n";

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

		std::cout << "Alive";

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

				if (progeny[i].Alive)
					std::cout << "Alive";
				else
					std::cout << "Dead";

			} // end of k loop
		} // end of i loop
	}

	void mutation()
	{
		int i;
		double dRnd1, dRnd2;
		int iRnd1;
		double mut_rate_per_ind, total_number_loci;
		int i_total_number_loci;
		int mutated_locus;
		double mutational_effect[2];
		double mutational_std_dev[2];

		mutational_std_dev[0] = sqrt(MutationalVariance[0]);
		mutational_std_dev[1] = sqrt(MutationalVariance[1]);
		i_total_number_loci = NumLociTrait0 + NumLociTrait1 + NumLociBoth;
		total_number_loci = i_total_number_loci;
		mut_rate_per_ind = 2.0 * total_number_loci * MutationRatePerLocus;

		for (i = 0; i < Nprogeny; i++)
		{
			dRnd1 = genrand();
			if (genrand() < mut_rate_per_ind)
			{
				iRnd1 = randnum(i_total_number_loci);
				if (iRnd1 < NumLociTrait0)
				{
					mutated_locus = iRnd1;
					dRnd2 = genrand();
					mutational_effect[0] = randnorm(0, mutational_std_dev[0]);
					if (dRnd2 < 0.5)
						progeny[i].Allele1trait0[mutated_locus] = progeny[i].Allele1trait0[mutated_locus]
						+ mutational_effect[0];
					else
						progeny[i].Allele2trait0[mutated_locus] = progeny[i].Allele2trait0[mutated_locus]
						+ mutational_effect[0];
				} // end of if -- mutation affecting trait 0 loci
				if (iRnd1 >= NumLociTrait0 && iRnd1 < (NumLociTrait0 + NumLociTrait1))
				{
					mutated_locus = iRnd1 - NumLociTrait0;
					dRnd2 = genrand();
					mutational_effect[1] = randnorm(0, mutational_std_dev[1]);
					if (dRnd2 < 0.5)
						progeny[i].Allele1trait1[mutated_locus] = progeny[i].Allele1trait1[mutated_locus]
						+ mutational_effect[1];
					else
						progeny[i].Allele2trait1[mutated_locus] = progeny[i].Allele2trait1[mutated_locus]
						+ mutational_effect[1];
				} // end of if -- mutation affecting trait 1 loci
				if (iRnd1 >= (NumLociTrait0 + NumLociTrait1))
				{
					mutated_locus = iRnd1 - (NumLociTrait0 + NumLociTrait1);
					dRnd2 = genrand();
					randbivnorm(mutational_std_dev[0], mutational_std_dev[1], MutationalCorrelation,
						mutational_effect[0], mutational_effect[1]);
					if (dRnd2 < 0.5)
					{
						progeny[i].Allele1both[mutated_locus][0] =
							progeny[i].Allele1both[mutated_locus][0] + mutational_effect[0];
						progeny[i].Allele1both[mutated_locus][1] =
							progeny[i].Allele1both[mutated_locus][1] + mutational_effect[1];
					}
					else
					{
						progeny[i].Allele2both[mutated_locus][0] =
							progeny[i].Allele2both[mutated_locus][0] + mutational_effect[0];
						progeny[i].Allele2both[mutated_locus][1] =
							progeny[i].Allele2both[mutated_locus][1] + mutational_effect[1];
					}
				} // end of if -- mutation affecting pleiotropic loci

				progeny[i].calculate_genotypic_values(NumLociTrait0, NumLociTrait1, NumLociBoth);
				progeny[i].calculate_phenotype(EnvironmentalStDev[0], EnvironmentalStDev[1]);

			} // end of if

		} // end of i loop

	}

	void natural_selection()
	{
		int i;
		double survival_prob, dRnum1, optimum[2];
		double dSu, dSv, dSc;
		double SSsqrt[2];
		dSc = 2 * (1 - SelectionalCorrelation * SelectionalCorrelation);
		SSsqrt[0] = sqrt(SelectionStrength[0]);
		SSsqrt[1] = sqrt(SelectionStrength[1]);


		optimum[0] = 0;
		optimum[1] = 0;

		for (i = 0; i < Nprogeny; i++)
		{
			// Univariate Gaussian selection on trait 0:
			//survival_prob = exp(-1.0 * (progeny[i].Phenotype[0] -
			//	optimum[0])*(progeny[i].Phenotype[0] - optimum[0]) /
			//	(2 * SelectionStrength[0]));

			// Bivariate selection on traits 0 and 1:
			dSu = (progeny[i].Phenotype[0] - optimum[0]) / SSsqrt[0];
			dSv = (progeny[i].Phenotype[1] - optimum[1]) / SSsqrt[1];
			survival_prob = exp((2 * SelectionalCorrelation*dSu*dSv - dSu * dSu - dSv * dSv) / dSc);

			dRnum1 = genrand();
			if (dRnum1 < survival_prob)
				progeny[i].Alive = true;
			else
				progeny[i].Alive = false;
		} // end of i loop

	}

	void population_regulation()
	{
		int i;
		double carrying_capacity_unfilled, progeny_left, keep_prob;
		int number_adults_chosen;
		double drnd1;

		progeny_left = 0;
		for (i = 0; i < Nprogeny; i++)
			if (progeny[i].Alive)
				progeny_left++;

		carrying_capacity_unfilled = CarryingCapacity;
		number_adults_chosen = 0;
		for (i = 0; i < Nprogeny; i++)
		{
			if (progeny[i].Alive)
			{
				keep_prob = carrying_capacity_unfilled / progeny_left;
				drnd1 = genrand();
				if (drnd1 < keep_prob)
				{
					adult[number_adults_chosen] = progeny[i];
					carrying_capacity_unfilled = carrying_capacity_unfilled - 1;
					number_adults_chosen++;
				}
				progeny_left = progeny_left - 1;
			} // end of if (progeny[i].Alive)
		} // end of i
		PopulationSize = number_adults_chosen;
	}

	void calculate_values_progeny()
	{
		int i;
		double dNP = Nprogeny;
		for (i = 0; i < 2; i++)
		{
			phenotypic_mean[i] = 0;
			genotypic_mean[i] = 0;
			phenotypic_variance[i] = 0;
			genotypic_variance[i] = 0;
		}
		phenotypic_covariance = 0;
		genotypic_covariance = 0;

		if (dNP > 0)
		{
			// Calculate the Means
			for (i = 0; i < Nprogeny; i++)
			{
				phenotypic_mean[0] = phenotypic_mean[0] + progeny[i].Phenotype[0];
				phenotypic_mean[1] = phenotypic_mean[1] + progeny[i].Phenotype[1];
				genotypic_mean[0] = genotypic_mean[0] + progeny[i].Genotype[0];
				genotypic_mean[1] = genotypic_mean[1] + progeny[i].Genotype[1];
			}
			phenotypic_mean[0] = phenotypic_mean[0] / dNP;
			phenotypic_mean[1] = phenotypic_mean[1] / dNP;
			genotypic_mean[0] = genotypic_mean[0] / dNP;
			genotypic_mean[1] = genotypic_mean[1] / dNP;

			// Calculate Variances and Covariances
			for (i = 0; i < Nprogeny; i++)
			{
				phenotypic_variance[0] = phenotypic_variance[0]
					+ (progeny[i].Phenotype[0] - phenotypic_mean[0])
					* (progeny[i].Phenotype[0] - phenotypic_mean[0]);
				phenotypic_variance[1] = phenotypic_variance[1]
					+ (progeny[i].Phenotype[1] - phenotypic_mean[1])
					* (progeny[i].Phenotype[1] - phenotypic_mean[1]);
				phenotypic_covariance = phenotypic_covariance
					+ (progeny[i].Phenotype[0] - phenotypic_mean[0])
					* (progeny[i].Phenotype[1] - phenotypic_mean[1]);

				genotypic_variance[0] = genotypic_variance[0]
					+ (progeny[i].Genotype[0] - genotypic_mean[0])
					* (progeny[i].Genotype[0] - genotypic_mean[0]);
				genotypic_variance[1] = genotypic_variance[1]
					+ (progeny[i].Genotype[1] - genotypic_mean[1])
					* (progeny[i].Genotype[1] - genotypic_mean[1]);
				genotypic_covariance = genotypic_covariance
					+ (progeny[i].Genotype[0] - genotypic_mean[0])
					* (progeny[i].Genotype[1] - genotypic_mean[1]);
			} // end of i
			phenotypic_variance[0] = phenotypic_variance[0] / dNP;
			phenotypic_variance[1] = phenotypic_variance[1] / dNP;
			genotypic_variance[0] = genotypic_variance[0] / dNP;
			genotypic_variance[1] = genotypic_variance[1] / dNP;
			phenotypic_covariance = phenotypic_covariance / dNP;
			genotypic_covariance = genotypic_covariance / dNP;
		}

		// Calculate the phenotypic and genotypic correlations
		double dtemp;
		dtemp = sqrt(phenotypic_variance[0] * phenotypic_variance[1]);
		if (dtemp > 0)
			phenotypic_correlation = phenotypic_covariance / dtemp;
		else
			phenotypic_correlation = 0;

		dtemp = sqrt(genotypic_variance[0] * genotypic_variance[1]);
		if (dtemp > 0)
			genotypic_correlation = genotypic_covariance / dtemp;
		else
			genotypic_correlation = 0;

	}

	void output_population_variables(int generation)
	{
		if (generation == 0) // Output the header in generation zero
		{
			std::cout << "\nGen\tzbar0\tzbar1\tP00\tP11\tP12\tr(P)\tgbar0\tgbar1\tG00\tG11\tG01\tr(G)";
		}

		std::cout << "\n" << generation;
		std::cout << std::setprecision(3) << std::fixed;
		std::cout << "\t" << phenotypic_mean[0];
		std::cout << "\t" << phenotypic_mean[1];
		std::cout << "\t" << phenotypic_variance[0];
		std::cout << "\t" << phenotypic_variance[1];
		std::cout << "\t" << phenotypic_covariance;
		std::cout << "\t" << phenotypic_correlation;
		std::cout << "\t" << genotypic_mean[0];
		std::cout << "\t" << genotypic_mean[1];
		std::cout << "\t" << genotypic_variance[0];
		std::cout << "\t" << genotypic_variance[1];
		std::cout << "\t" << genotypic_covariance;
		std::cout << "\t" << genotypic_correlation;
	}

	int getNumberOfGenerations()
	{
		return NumberOfGenerations;
	}

	void save_population_variables(int generation)
	{
		std::ofstream outfile;
		if (generation == 0) // Output the header in generation zero
		{
			outfile.open("output.csv");
			outfile << "Gen,N,zbar0,zbar1,P00,P11,P12,r(P),gbar0,gbar1,G00,G11,G01,r(G)";
			outfile.close();
		}

		outfile.open("output.csv", std::fstream::app);
		outfile << "\n" << generation;
		outfile << "," << PopulationSize;
		outfile << "," << phenotypic_mean[0];
		outfile << "," << phenotypic_mean[1];
		outfile << "," << phenotypic_variance[0];
		outfile << "," << phenotypic_variance[1];
		outfile << "," << phenotypic_covariance;
		outfile << "," << phenotypic_correlation;
		outfile << "," << genotypic_mean[0];
		outfile << "," << genotypic_mean[1];
		outfile << "," << genotypic_variance[0];
		outfile << "," << genotypic_variance[1];
		outfile << "," << genotypic_covariance;
		outfile << "," << genotypic_correlation;
		outfile.close();
	}

	void gaussian_mating()
	{
		// This function implements gaussian mate choice in a
		// polygynous mating system. It should be used instead
		// of other mating functions (like polygynous_mating).
		// In this function, each female mates at most once, and
		// each male can mate an unlimited number of times, so
		// the mating system is polygynous.

		// Females choose males based on their trait values.
		// Mating preferences are Gaussian in the sense that each female
		// has an ideal preferred male phenotype and her preferences
		// fall off as males depart from her preferred phenotype.
		// The drop in mating probability as the male departs from
		// the preferred phenotype is modeled as a Gaussian-shaped
		// function.

		// Absolute or relative preferences? This function uses relative
		// preferences, and they are relative to the phenotypic mean of 
		// the males. Consequently, if a female's preference trait has
		// a value of 0.72, then that particular female's ideal male
		// has a trait value 0.72 units greater than the male mean trait
		// value. If she encountered such a male, she would mate with
		// probability 1. Her probability of mating would drop off for
		// males with trait values larger or smaller than the preferred
		// value.

		int i, j, m;
		int iPC;
		bool mate_found;
		int mateID, counter, rnum;
		double dRnum;
		double mate_prob;
		double dNmales;
		double mean_male_trait0;
		double ideal;

		// Check to make sure at least one male is present in the population.
		// Also calculate the mean of the male ornament trait (trait 0). 

		bool males_present = false;
		dNmales = 0;
		mean_male_trait0 = 0;
		for (i = 0; i < PopulationSize; i++)
		{
			if (!adult[i].Female)
			{
				males_present = true;
				dNmales++;
				mean_male_trait0 = mean_male_trait0 + adult[i].Phenotype[0];
			}
		}

		if (dNmales > 0)
		{
			mean_male_trait0 = mean_male_trait0 / dNmales;
		}
		else
		{
			mean_male_trait0 = 0;
			PopulationExtinct = true;
		}

		iPC = 0;
		for (i = 0; i < PopulationSize; i++)
		{
			if (adult[i].Female && males_present)
			{
				// Find a mate for this female
				// She only gets MaxMatingEncounters tries to find
				// someone. We also include encounters with females
				// in the count, so it will be harder to find a mate
				// when the sex ratio is extremely female-biased.
				// If she doesn't find a mate in the allotted number
				// of tries, then she produces no progeny.

				// For each female, we first have to determine her
				// ideal mate phenotype. Trait 1 is the female preference
				// trait and it describes the deviation from the male
				// mean of her preferred mate at Trait 0, the ornament. 

				ideal = adult[i].Phenotype[1] + mean_male_trait0;

				mate_found = false;
				counter = 0;
				while (counter < MaxMatingEncounters && !mate_found)
				{
					rnum = randnum(PopulationSize);
					if (!adult[rnum].Female)
					{
						// Calculate the focal female's (adult[i]) probability 
						// of mating with this random male (adult[rnum]). The
						// probability is calculated using the Gaussian 
						// probability density function without the normalization
						// term.
						mate_prob = exp(-1 * (adult[rnum].Phenotype[0] - ideal)*
							(adult[rnum].Phenotype[0] - ideal)
							/ (2 * GaussianPreferenceVariance));

						// Generate a random number (roll the dice!)
						dRnum = genrand();

						if (dRnum < mate_prob)
						{
							mateID = rnum;
							mate_found = true;
						}
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
							// if this number is less than 0.5 we choose one allele. 
							// Otherwise, we choose the other allele.

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
						  // pleiotropic loci. [Note that this part isn’t properly indented in
						  // this text document – it should be indented a bit more in your actual
						  // .cpp file. Keep in mind, however, that indentation is just cosmetic.
						  // The C++ compiler ignores all spaces and tabs.]
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

	bool is_extinct()
	{
		return PopulationExtinct;
	}

};
