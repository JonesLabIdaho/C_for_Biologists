#include <iostream>
#include <iomanip>
#include <fstream>

int main()
{
	double Pstart, Qstart, Waa, Wab, Wbb;
	int NumberOfGenerations;
	int i;
	double p, q;
	double Pprime, Qprime;
	double Wbar;

	std::ofstream out_file;
	out_file.open("output.csv");

	// Set parameters and starting conditions here 

	Pstart = 0.5;
	Qstart = 1 - Pstart;
	Waa = 1.2;
	Wab = 1.0;
	Wbb = 0.8;
	NumberOfGenerations = 20;

	// Single-locus, constant viability selection 
	// in a diploid population with two alleles.

	p = Pstart;
	q = Qstart;

	std::cout << "Gen\tp\tq\n";
	std::cout << "0\t" << Pstart << "\t" << Qstart;

	out_file << "Gen,p,q\n";
	out_file << "0," << Pstart << "," << Qstart;

	for (i = 0; i < NumberOfGenerations; i++)
	{
		Wbar = p*p*Waa + 2 * p*q*Wab + q*q*Wbb;
		Pprime = p*(p*Waa + q*Wab) / Wbar;
		Qprime = q*(p*Wab + q*Wbb) / Wbar;

		std::cout << std::setprecision(4) << "\n" << i + 1 << "\t" << Pprime << "\t" << Qprime;
		out_file << "\n" << i + 1 << "," << Pprime << "," << Qprime;

		p = Pprime;
		q = Qprime;

	} // end of i loop

	out_file.close();

	char end_it;
	std::cout << "\n\nEnter any character to exit...";
	std::cin >> end_it;

	return 0;
}
