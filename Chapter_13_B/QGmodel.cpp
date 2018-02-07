#include <iostream>
#include "MTwisterFunctions.h"
#include "simulation_engine.h"

int main()
{
	int generations, ig, rep, number_of_replicates;
	number_of_replicates = 20;

	std::string outfile_name, extended_outfile_name;
	outfile_name = "outfile";

	std::string parameter_filename;
	parameter_filename = outfile_name + "_parameters.csv";

	simulation_engine *my_sim = new simulation_engine[number_of_replicates];
	mean_recorder *mean_array = new mean_recorder[number_of_replicates];

	my_sim[0].save_parameter_values(parameter_filename);

	for (rep = 0; rep < number_of_replicates; rep++)
	{
		
		extended_outfile_name = outfile_name + "_run_" + std::to_string(rep + 1) + ".csv";

		bool initialization_success;
		initialization_success = my_sim[rep].initialize_population();
		if (!initialization_success)
		{
			std::cout << "\nSimulation Initialization Failure!\n";
			return 0;
		}

		for (ig = 0; ig < my_sim[rep].getNumberOfInitialGenerations(); ig++)
		{
			my_sim[rep].monogamy();
			my_sim[rep].mutation();
			my_sim[rep].selection(true);
			my_sim[rep].population_regulation();
		}

		// Turn the lifecycle halfway, so the
		// experimental generations start with progeny
		my_sim[rep].calculate_values_progeny();
		my_sim[rep].monogamy();
		my_sim[rep].mutation();

		for (generations = 0; generations < my_sim[rep].getNumberOfGenerations(); generations++)
		{
			if (!my_sim[rep].is_extinct())
			{
				my_sim[rep].selection(false);
				my_sim[rep].calculate_values_progeny();
				my_sim[rep].population_regulation();
				my_sim[rep].monogamy();
				my_sim[rep].mutation();
				my_sim[rep].calculate_values_adults();
				my_sim[rep].store_variables_in_memory();
			}
		}

		my_sim[rep].save_stored_variables(extended_outfile_name);

		my_sim[rep].calc_run_means();
		my_sim[rep].append_means_to_output_file(extended_outfile_name);

		mean_array[rep] = my_sim[rep].report_means();
	}

	int i, j, iNmeans;
	double dNreps = number_of_replicates;
	std::vector<double> mean_of_means;
	std::vector<double> std_dev_of_means;

	iNmeans = static_cast<int>(mean_array[0].mean_list.size());
	for (i = 0; i < iNmeans; i++)
	{
		mean_of_means.push_back(0);
		std_dev_of_means.push_back(0);
	}

	// Calculate the means of the means
	for (i = 0; i < number_of_replicates; i++)
	{
		for (j = 0; j < iNmeans; j++)
		{
			mean_of_means[j] = mean_of_means[j] + mean_array[i].mean_list[j];
		}
	}

	for (j = 0; j < iNmeans; j++)
		mean_of_means[j] = mean_of_means[j] / dNreps;

	// Calculate the standard deviations of the means
	for (i = 0; i < number_of_replicates; i++)
	{
		for (j = 0; j < iNmeans; j++)
		{
			std_dev_of_means[j] = std_dev_of_means[j] +
				(mean_array[i].mean_list[j] - mean_of_means[j]) *
				(mean_array[i].mean_list[j] - mean_of_means[j]);
		}
	}

	for (j = 0; j < iNmeans; j++)
	{
		std_dev_of_means[j] = std_dev_of_means[j] / dNreps;
		std_dev_of_means[j] = sqrt(std_dev_of_means[j]);
	}

	// Save to a file
	std::string summary_filename;
	summary_filename = outfile_name + "_summary.csv";
	std::ofstream outfile;
	outfile.open(summary_filename);
	outfile << "Rep,N,zbar0,zbar1,P00,P11,P12,r(P),gbar0,gbar1,G00,G11,G01,r(G)";
	outfile << ",Lambda1,Lambda2,EvecX,EvecY,Angle,Size,Eccen,Strt0,Strt1";
	outfile << ",~Lmbd1,~Lmbd2,~Ang,~Size,~Ecc,~G00,~G11,~G01,~r(g)";
	outfile << ",ASR,Im,If,MdifM,MdifF";

	for (i = 0; i < number_of_replicates; i++)
	{
		outfile << "\nRep" << i + 1;
		for (j = 0; j < iNmeans; j++)
		{
			outfile << "," << mean_array[i].mean_list[j];
		}
	}

	outfile << "\nMean_of_means";
	for (j = 0; j < iNmeans; j++)
	{
		outfile << "," << mean_of_means[j];
	}

	outfile << "\nStdev_of_means";
	for (j = 0; j < iNmeans; j++)
	{
		outfile << "," << std_dev_of_means[j];
	}

	outfile.close();

	delete[] my_sim;
	delete[] mean_array;

	char end_it;
	std::cout << "\n\nEnter any character to exit...";
	std::cin >> end_it;
	return 0;
}


