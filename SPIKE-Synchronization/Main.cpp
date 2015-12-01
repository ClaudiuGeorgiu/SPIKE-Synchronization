#include "SPIKESynchronization.h"
#include <iostream>
#include <vector>

using namespace std;

int main(void)
{
	vector<int> input1 = { 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1 };
	vector<int> input2 = { 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1 };

	for (int n = 0; n < input1.size(); ++n)
	{
		cout << input1[n] << "   ";
	}
	cout << "\n";

	for (int n = 0; n < input2.size(); ++n)
	{
		cout << input2[n] << "   ";
	}
	cout << "\n";

	SPIKESynchronization* spike = new SPIKESynchronization();

	vector<int> result = spike->CoincidenceVector(input1, input2);
	for (int n = 0; n < result.size(); ++n)
	{
		cout << result[n] << "   ";
	}
	cout << "\n";

	if (spike != NULL)
		delete spike;

	getchar();

	return 0;
}