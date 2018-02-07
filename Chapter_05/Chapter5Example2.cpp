
#include <iostream>
#include <string>
#include <time.h>

int main()
{
	std::string name[6];
	size_t i;
	int area_code[6];
	int prefix[6];
	int number[6];


	name[0] = "Mike";
	name[1] = "Sarah";
	name[2] = "Fred";
	name[3] = "Carrie";
	name[4] = "Ted";
	name[5] = "Julie";

	std::srand(time(NULL));
	for (i = 0; i < 6; i++)
	{
		area_code[i] = rand() % 900 + 100;
		prefix[i] = rand() % 900 + 100;
		number[i] = rand() % 9000 + 1000;
	}

	std::cout << "Name\tNumber\n";
	for (i = 0; i < 6; i++)
	{
		std::cout << name[i] << "\t(" << area_code[i] << ") " << prefix[i] << "-" << number[i] << "\n";
	}

	char end_it;
	std::cout << "\n\nEnter any character to exit...";
	std::cin >> end_it;

	return 0;
}
