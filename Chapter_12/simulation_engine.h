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
	double MatingSuccess;

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
	int NumberOfInitialGenerations;
	double InitialSelectionStrength[2], InitialSelectionalCorrelation;
	double InitialOptimum[2];
	double Optimum[2];
	bool InitialSelectionSexLimited, ExperimentalSelectionSexLimited;

	// Variables corresponding to population-level summary statistics:
	double phenotypic_mean[2], genotypic_mean[2];
	double phenotypic_variance[2], genotypic_variance[2];
	double phenotypic_covariance, genotypic_covariance;
	double phenotypic_correlation, genotypic_correlation;
	double EigenValue[2];
	double EigenVector1[2];
	double EigenVector2[2];
	double LeadAngle, Sigma, Epsilon;
	double sel_diff_trt_0, sel_diff_trt_1;

	double PrevEval[2], PrevAngle, PrevSigma;
	double PrevEpsilon, PrevG00, PrevG11, PrevG01, PrevRg;
	double cPrevEval[2], cPrevAngle, cPrevSigma;
	double cPrevEpsilon, cPrevG00, cPrevG11, cPrevG01, cPrevRg;

	double SexRatio, Im, If, MdiffMales, MdiffFemales;

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

		// Initial Generations Parameters
		NumberOfInitialGenerations = 1000;
		InitialSelectionStrength[0] = 49;
		InitialSelectionStrength[1] = 49;
		InitialSelectionalCorrelation = 0;
		InitialOptimum[0] = 0;
		InitialOptimum[1] = 0;

		// Demographic Parameters
		NumberOfGenerations = 2000;
		PopulationSize = 500;
		CarryingCapacity = PopulationSize;
		Fecundity = 4;
		PopulationExtinct = false;

		// Mating Parameters
		MaxMatingEncounters = 50;
		GaussianPreferenceVariance = 9;

		// Genetic Parameters
		NumLociTrait0 = 25;
		NumLociTrait1 = 25;
		NumLociBoth = 25;
		EnvironmentalVariance[0] = 1;
		EnvironmentalVariance[1] = 1;
		EnvironmentalStDev[0] = sqrt(EnvironmentalVariance[0]);
		EnvironmentalStDev[1] = sqrt(EnvironmentalVariance[1]);

		// Mutational Parameters
		MutationalVariance[0] = 0.05;
		MutationalVariance[1] = 0.05;
		MutationalCorrelation = 0;
		MutationRatePerLocus = 0.0002;

		// Selection Parameters
		SelectionStrength[0] = 49;
		SelectionStrength[1] = 49;
		SelectionalCorrelation = 0;
		Optimum[0] = 0;
		Optimum[1] = 0;
		InitialSelectionSexLimited = true;
		ExperimentalSelectionSexLimited = true;

		// Initialize previous generation variables to zero
		PrevEval[0] = 0;
		PrevEval[1] = 0;
		PrevAngle = 0;
		PrevSigma = 0;
		PrevEpsilon = 0;
		PrevG00 = 0;
		PrevG11 = 0;
		PrevG01 = 0;
		PrevRg = 0;

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
		std::cout << "Optimum_Trait0    \t" << Optimum[0] << "\n";
		std::cout << "Optimum_Trait1    \t" << Optimum[1] << "\n";
		if (ExperimentalSelectionSexLimited)
			std::cout << "Sex_Limited_Sel:  \ttrue\n";
		else
			std::cout << "Sex_Limited_Sel:  \tfalse\n";

		// Initial Generations
		std::cout << "Initial_Generations_Parameters:\n";
		std::cout << "No_Initial_Gens:  \t" << NumberOfInitialGenerations << "\n";
		std::cout << "Init_Omega_Trait0:\t" << InitialSelectionStrength[0] << "\n";
		std::cout << "Init_Omega_Trait1:\t" << InitialSelectionStrength[1] << "\n";
		std::cout << "Init_Sel_Corr:    \t" << InitialSelectionalCorrelation << "\n";
		std::cout << "Init_Opt_Trt0     \t" << InitialOptimum[0] << "\n";
		std::cout << "Init_Opt_Trt1     \t" << InitialOptimum[1] << "\n";
		if (InitialSelectionSexLimited)
			std::cout << "Init_Sex_Lim_Sel: \ttrue\n";
		else
			std::cout << "Init_Sex_Lim_Sel: \tfalse\n";
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
					adult[i].MatingSuccess++;
					adult[mateID].MatingSuccess++;
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

			progeny[i].MatingSuccess = 0;

		} // end of i loop

	}

	void natural_selection(bool init_gens)
	{
		int i;
		double survival_prob, dRnum1;
		double dSu, dSv, dSc;
		double SSsqrt[2];
		double local_sel_str[2], local_sel_corr, local_opt[2];

		if (init_gens)
		{
			for (i = 0; i < 2; i++)
			{
				local_sel_str[i] = InitialSelectionStrength[i];
				local_opt[i] = InitialOptimum[i];
			}
			local_sel_corr = InitialSelectionalCorrelation;
		}
		else
		{
			for (i = 0; i < 2; i++)
			{
				local_sel_str[i] = SelectionStrength[i];
				local_opt[i] = Optimum[i];
			}
			local_sel_corr = SelectionalCorrelation;
		}

		dSc = 2 * (1 - local_sel_corr * local_sel_corr);
		SSsqrt[0] = sqrt(local_sel_str[0]);
		SSsqrt[1] = sqrt(local_sel_str[1]);

		for (i = 0; i < Nprogeny; i++)
		{
			//Bivariate selection on traits 0 and 1:
			if (SSsqrt[0] > 0 && SSsqrt[1] > 0)
			{
				dSu = (progeny[i].Phenotype[0] - local_opt[0]) / SSsqrt[0];
				dSv = (progeny[i].Phenotype[1] - local_opt[1]) / SSsqrt[1];
				survival_prob = exp((2 * local_sel_corr*dSu*dSv - dSu * dSu - dSv * dSv) / dSc);
			}

			// Selection on trait 0 only
			if (SSsqrt[0] > 0 && SSsqrt[1] == 0)
			{
				survival_prob = exp(-1.0 *
					(progeny[i].Phenotype[0] - local_opt[0])*
					(progeny[i].Phenotype[0] - local_opt[0])
					/ (2 * local_sel_str[0]));
			}

			// Selection on trait 1 only
			if (SSsqrt[0] == 0 && SSsqrt[1] > 0)
			{
				survival_prob = exp(-1.0 *
					(progeny[i].Phenotype[1] - local_opt[1])*
					(progeny[i].Phenotype[1] - local_opt[1])
					/ (2 * local_sel_str[1]));
			}

			// No selection on either trait
			if (SSsqrt[0] == 0 && SSsqrt[1] == 0)
			{
				survival_prob = 1;
			}

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

		// Check to make sure the variances are both non-zero before
		// attempting to calculate the G-matrix.

		if (genotypic_variance[0] > 0 && genotypic_variance[1] > 0)
		{

			// Calculate the eigenvectors and eigenvalues of the G-matrix

			double a, b, c, d;
			double ev[2];
			a = genotypic_variance[0];
			b = genotypic_covariance;
			c = genotypic_covariance;
			d = genotypic_variance[1];
			const double m_pi = 3.14159265358979323846;

			// Use the formula for a 2x2 matrix to calculate eigenvalues
			if ((a + d)*(a + d) > 4 * (a*d - b * c))
			{
				ev[0] = ((a + d) + sqrt((a + d)*(a + d) - 4 * (a*d - b * c))) / 2;
				ev[1] = ((a + d) - sqrt((a + d)*(a + d) - 4 * (a*d - b * c))) / 2;
			}
			else
			{
				ev[0] = 0;
				ev[1] = 0;
			}

			// Make sure EigenValue[0] is the larger of the two
			if (ev[0] > ev[1])
			{
				EigenValue[0] = ev[0];
				EigenValue[1] = ev[1];
			}
			else
			{
				EigenValue[0] = ev[1];
				EigenValue[1] = ev[0];
			}

			// Calculate eigenvectors using 2x2 matrix formulae
			if (c != 0)
			{
				EigenVector1[0] = (EigenValue[0] - d) / sqrt((EigenValue[0] - d)*(EigenValue[0] - d) + c * c);
				EigenVector1[1] = c / sqrt((EigenValue[0] - d)*(EigenValue[0] - d) + c * c);
				EigenVector2[0] = (EigenValue[1] - d) / sqrt((EigenValue[1] - d)*(EigenValue[1] - d) + c * c);
				EigenVector2[1] = c / sqrt((EigenValue[1] - d)*(EigenValue[1] - d) + c * c);
			}
			if (c == 0)
			{
				if (d > a)
				{
					EigenVector1[0] = 0;
					EigenVector1[1] = 1;
					EigenVector2[0] = 1;
					EigenVector2[1] = 0;
				}
				else
				{
					EigenVector1[0] = 1;
					EigenVector1[1] = 0;
					EigenVector2[0] = 0;
					EigenVector2[1] = 1;
				}
			}
			// Calculate the angle of the leading eigenvector in degrees
			if (EigenVector1[0] != 0)
			{
				LeadAngle = atan(EigenVector1[1] / EigenVector1[0])*(180 / m_pi);
			}
			else
			{
				LeadAngle = 90;
			}

			Sigma = EigenValue[0] + EigenValue[1];
			Epsilon = EigenValue[1] / EigenValue[0];
		}
		else // deal with the case where one or both genetic variances are zero
		{
			// If either variance is zero, then the covariance is also zero.
			// If the covariance is zero, then the eigenvalues are the same as the variances.
			EigenValue[0] = genotypic_variance[0];
			EigenValue[1] = genotypic_variance[1];
			Sigma = EigenValue[0] + EigenValue[1];
			Epsilon = 0; // If either variance is zero, so is epsilon

						 // if only one variance is zero, we will define the angle to 
						 // point along the axis of the trait with the non-zero variance.

			EigenVector1[0] = 0;
			EigenVector1[1] = 0;
			EigenVector2[0] = 0;
			EigenVector2[1] = 0;
			LeadAngle = 0;

			if (EigenValue[0] > 0)
			{
				EigenVector1[0] = 1;
				EigenVector1[1] = 0;
				EigenVector2[0] = 0;
				EigenVector2[1] = 1;
			}
			if (EigenValue[1] > 0)
			{
				EigenVector1[0] = 0;
				EigenVector1[1] = 1;
				EigenVector2[0] = 1;
				EigenVector2[1] = 0;
				LeadAngle = 90;
			}
		}

		// Calculate selection differentials on trait 0 and trait 1
		// These selection differentials represent viability 
		// selection, because they include only progeny survival
		// during the juvenile phase of the lifecycle.

		// The selection differential is the covariance between
		// trait values and fitness. Here, fitness is based on 
		// whether or not the individual survived viability selection.

		// We use some of the variables calculated above:
		// phenotypic_mean[2], dNP (number of progeny).

		double mean_fitness;
		if (dNP > 0)
		{
			// Tally the mean fitness of the offspring
			mean_fitness = 0;
			for (i = 0; i < Nprogeny; i++)
			{
				if (progeny[i].Alive)
				{
					mean_fitness++;
				}
			}
			mean_fitness = mean_fitness / dNP;
		}

		if (dNP > 0 && mean_fitness > 0 && mean_fitness < 1)
		{
			sel_diff_trt_0 = 0;
			sel_diff_trt_1 = 0;
			for (i = 0; i < Nprogeny; i++)
			{
				if (progeny[i].Alive) // Alive and dead have different fitnesses (1 vs. 0)
				{
					sel_diff_trt_0 = sel_diff_trt_0 + (progeny[i].Phenotype[0] 
						- phenotypic_mean[0]) * (1.0 / mean_fitness - 1.0);
					sel_diff_trt_1 = sel_diff_trt_1 + (progeny[i].Phenotype[1] 
						- phenotypic_mean[1]) * (1.0 / mean_fitness - 1.0);
				}
				else
				{
					sel_diff_trt_0 = sel_diff_trt_0 + (progeny[i].Phenotype[0] 
						- phenotypic_mean[0]) * (-1.0);
					sel_diff_trt_1 = sel_diff_trt_1 + (progeny[i].Phenotype[1] 
						- phenotypic_mean[1]) * (-1.0);
				}
			}
			sel_diff_trt_0 = sel_diff_trt_0 / dNP;
			sel_diff_trt_1 = sel_diff_trt_1 / dNP;
		}
		else
		{
			sel_diff_trt_0 = 0;
			sel_diff_trt_1 = 0;
		}

		cPrevEval[0] = fabs(PrevEval[0] - EigenValue[0]);
		cPrevEval[1] = fabs(PrevEval[1] - EigenValue[1]);
		cPrevAngle = fabs(PrevAngle - LeadAngle);
		cPrevSigma = fabs(PrevSigma - Sigma);
		cPrevEpsilon = fabs(PrevEpsilon - Epsilon);
		cPrevG00 = fabs(PrevG00 - genotypic_variance[0]);
		cPrevG11 = fabs(PrevG11 - genotypic_variance[1]);
		cPrevG01 = fabs(PrevG01 - genotypic_covariance);
		cPrevRg = fabs(PrevRg - genotypic_correlation);

		// Ensure that we are calculating the smallest 
		// possible change in eigenvector angle.
		if (cPrevAngle > 90)
			cPrevAngle = 180 - cPrevAngle;

		PrevEval[0] = EigenValue[0];
		PrevEval[1] = EigenValue[1];
		PrevAngle = LeadAngle;
		PrevSigma = Sigma;
		PrevEpsilon = Epsilon;
		PrevG00 = genotypic_variance[0];
		PrevG11 = genotypic_variance[1];
		PrevG01 = genotypic_covariance;
		PrevRg = genotypic_correlation;

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
			outfile << ",Lambda1,Lambda2,EvecX,EvecY,Angle,Size,Eccen,Strt0,Strt1";
			outfile << ",~Lmbd1,~Lmbd2,~Ang,~Size,~Ecc,~G00,~G11,~G01,~r(g)";
			outfile << ",ASR,Im,If,MdifM,MdifF";
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
		outfile << "," << EigenValue[0];
		outfile << "," << EigenValue[1];
		outfile << "," << EigenVector1[0];
		outfile << "," << EigenVector1[1];
		outfile << "," << LeadAngle;
		outfile << "," << Sigma;
		outfile << "," << Epsilon;
		outfile << "," << sel_diff_trt_0;
		outfile << "," << sel_diff_trt_1;
		outfile << "," << cPrevEval[0];
		outfile << "," << cPrevEval[1];
		outfile << "," << cPrevAngle;
		outfile << "," << cPrevSigma;
		outfile << "," << cPrevEpsilon;
		outfile << "," << cPrevG00;
		outfile << "," << cPrevG11;
		outfile << "," << cPrevG01;
		outfile << "," << cPrevRg;
		outfile << "," << SexRatio;
		outfile << "," << Im;
		outfile << "," << If;
		outfile << "," << MdiffMales;
		outfile << "," << MdiffFemales;
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
		int iNmales;
		int iMaleIndex;
		int *iMaleList = new int[PopulationSize];

		// Check to make sure at least one male is present in the population.
		// Also calculate the mean of the male ornament trait (trait 0). 

		bool males_present = false;
		iNmales = 0;
		mean_male_trait0 = 0;
		for (i = 0; i < PopulationSize; i++)
		{
			if (!adult[i].Female)
			{
				males_present = true;
				iMaleList[iNmales] = i;
				iNmales++;
				mean_male_trait0 = mean_male_trait0 + adult[i].Phenotype[0];
			}
		}
		dNmales = iNmales;

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
				// someone. 

				// She only encounters males.

				// If she doesn't find a mate in the allotted number
				// of tries, then she produces no progeny.

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
					rnum = randnum(iNmales);
					iMaleIndex = iMaleList[rnum];

					if (!adult[rnum].Female)
					{
						// Calculate the focal female's (adult[i]) probability 
						// of mating with this random male (adult[iMaleIndex]). The
						// probability is calculated using the Gaussian 
						// probability density function without the normalization
						// term.
						if (GaussianPreferenceVariance > 0)
						{
							mate_prob = exp(-1 * (adult[iMaleIndex].Phenotype[0] - ideal)*
								(adult[iMaleIndex].Phenotype[0] - ideal)
								/ (2 * GaussianPreferenceVariance));
						}
						else
						{
							mate_prob = 1;
						}

						// Generate a random number (roll the dice!)
						dRnum = genrand();

						if (dRnum < mate_prob)
						{
							mateID = iMaleIndex;
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
					adult[i].MatingSuccess++;
					adult[mateID].MatingSuccess++;
				} // end of if (mate_found)
			} // end of if (adult[i].Female && males_present)
		} // end of i loop
		Nprogeny = iPC;
		delete[] iMaleList;
	}

	bool is_extinct()
	{
		return PopulationExtinct;
	}

	void sex_limited_selection(bool init_gens)
	{
		int i;
		double survival_prob, dRnum1, optimum[2];
		double local_sel_str[2];

		if (init_gens)
		{
			local_sel_str[0] = InitialSelectionStrength[0];
			local_sel_str[1] = InitialSelectionStrength[1];
			optimum[0] = InitialOptimum[0];
			optimum[1] = InitialOptimum[1];
		}
		else
		{
			local_sel_str[0] = SelectionStrength[0];
			local_sel_str[1] = SelectionStrength[1];
			optimum[0] = Optimum[0];
			optimum[1] = Optimum[1];
		}

		for (i = 0; i < Nprogeny; i++)
		{
			// Calculate survival probabilities for males and females separately
			if (!progeny[i].Female)
			{
				if (local_sel_str[0] > 0)
				{
					survival_prob = exp(-1.0 *
						(progeny[i].Phenotype[0] - optimum[0])*
						(progeny[i].Phenotype[0] - optimum[0])
						/ (2 * local_sel_str[0]));
				}
				else
				{
					survival_prob = 1; // SelectionStrength = 0 means no selection.
				}
			}
			else
			{
				if (local_sel_str[1] > 0)
				{
					survival_prob = exp(-1.0 *
						(progeny[i].Phenotype[1] - optimum[1])*
						(progeny[i].Phenotype[1] - optimum[1])
						/ (2 * local_sel_str[1]));
				}
				else
				{
					survival_prob = 1; // SelectionStrength = 0 means no selection.
				}
			}

			// Generate a random number to see if the individual survives
			dRnum1 = genrand();
			if (dRnum1 < survival_prob)
				progeny[i].Alive = true;
			else
				progeny[i].Alive = false;
		} // end of i loop
	}

	int getNumberOfInitialGenerations()
	{
		return NumberOfInitialGenerations;
	}

	void calculate_values_adults()
	{
		int i;
		double MeanMSmales, MeanMSfemales, MeanTrait0males, MeanTrait1females;
		double VarMSmales, VarMSfemales, Nma, Nfe;

		Nma = 0;
		Nfe = 0;
		MeanMSmales = 0;
		MeanMSfemales = 0;
		MeanTrait0males = 0;
		MeanTrait1females = 0;

		// Calculate mean mating success and mean trait values
		for (i = 0; i < PopulationSize; i++)
		{
			if (adult[i].Female)
			{
				Nfe++;
				MeanMSfemales = MeanMSfemales + adult[i].MatingSuccess;
				MeanTrait1females = MeanTrait1females + adult[i].Phenotype[1];
			}
			else
			{
				Nma++;
				MeanMSmales = MeanMSmales + adult[i].MatingSuccess;
				MeanTrait0males = MeanTrait0males + adult[i].Phenotype[0];
			}
		} // end of i

		if (Nma > 0)
		{
			MeanMSmales = MeanMSmales / Nma;
			MeanTrait0males = MeanTrait0males / Nma;
		}
		else
		{
			MeanMSmales = 0;
			MeanTrait0males = 0;
		}

		if (Nfe > 0)
		{
			MeanMSfemales = MeanMSfemales / Nfe;
			MeanTrait1females = MeanTrait1females / Nfe;
		}
		else {
			MeanMSfemales = 0;
			MeanTrait1females = 0;
		}

		// Calculate the Sex Ratio.
		if (Nma + Nfe > 0)
			SexRatio = Nma / (Nma + Nfe);
		else
			SexRatio = 0;

		// Calculate opportunities for sexual selection
		// and mating differentials

		VarMSmales = 0;
		VarMSfemales = 0;
		MdiffMales = 0;
		MdiffFemales = 0;

		for (i = 0; i < PopulationSize; i++)
		{
			if (adult[i].Female)
			{
				VarMSfemales = VarMSfemales + (adult[i].MatingSuccess - MeanMSfemales)
					*(adult[i].MatingSuccess - MeanMSfemales);
				if (MeanMSfemales > 0)
					MdiffFemales = MdiffFemales + (adult[i].Phenotype[1] - MeanTrait1females)
					*(adult[i].MatingSuccess / MeanMSfemales - 1);
			}
			else
			{
				VarMSmales = VarMSmales + (adult[i].MatingSuccess - MeanMSmales)
					*(adult[i].MatingSuccess - MeanMSmales);
				if (MeanMSmales > 0)
					MdiffMales = MdiffMales + (adult[i].Phenotype[0] - MeanTrait0males)
					*(adult[i].MatingSuccess / MeanMSmales - 1);
			}
		} // end of i

		if (Nfe > 0)
		{
			VarMSfemales = VarMSfemales / Nfe;
			MdiffFemales = MdiffFemales / Nfe;
		}
		else
		{
			VarMSfemales = 0;
			MdiffFemales = 0;
		}

		if (Nma > 0)
		{
			VarMSmales = VarMSmales / Nma;
			MdiffMales = MdiffMales / Nma;
		}
		else
		{
			VarMSmales = 0;
			MdiffMales = 0;
		}

		if (MeanMSmales > 0)
			Im = VarMSmales / (MeanMSmales*MeanMSmales);
		else
			Im = 0;

		if (MeanMSfemales > 0)
			If = VarMSfemales / (MeanMSfemales*MeanMSfemales);
		else
			If = 0;


	}

	void selection(bool initial_generations)
	{
		if (initial_generations)
		{
			if (InitialSelectionSexLimited)
				sex_limited_selection(true);
			else
				natural_selection(true);
		}
		else
		{
			if (ExperimentalSelectionSexLimited)
				sex_limited_selection(false);
			else
				natural_selection(false);
		}
	}


};
