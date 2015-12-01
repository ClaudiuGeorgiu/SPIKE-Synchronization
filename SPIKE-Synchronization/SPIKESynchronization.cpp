#include "SPIKESynchronization.h"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

SPIKESynchronization::SPIKESynchronization()
{
}

SPIKESynchronization::~SPIKESynchronization()
{
}

int SPIKESynchronization::getPreviousSpikeIndex(vector<int> input, int index)
{
	if (index < 1 || index >= input.size())
		return -1;

	for (int n = index - 1; n >= 0; --n)
	{
		if (input[n] == 1)
			return n;
	}

	return -1;
}

int SPIKESynchronization::getNextSpikeIndex(vector<int> input, int index)
{
	if (index < 0 || index >= input.size() - 1)
		return -1;

	for (int n = index + 1; n < input.size(); ++n)
	{
		if (input[n] == 1)
			return n;
	}

	return -1;
}

double SPIKESynchronization::getTau(vector<int> input1, vector<int> input2, int index1, int index2)
{
	vector<int> temp = vector<int>();

	int nextSpike1 = getNextSpikeIndex(input1, index1);
	if (nextSpike1 != -1)
		temp.push_back(nextSpike1 - index1); // v_i (1)

	int prevSpike1 = getPreviousSpikeIndex(input1, index1);
	if (prevSpike1 != -1)
		temp.push_back(index1 - prevSpike1); // v_(i-1) (1)

	int nextSpike2 = getNextSpikeIndex(input2, index2);
	if (nextSpike2 != -1)
		temp.push_back(nextSpike2 - index2); // v_j (2)

	int prevSpike2 = getPreviousSpikeIndex(input2, index2);
	if (prevSpike2 != -1)
		temp.push_back(index2 - prevSpike2); // v_(j-1) (2)

	return 0.5 * *min_element(begin(temp), end(temp));
}

vector<int> SPIKESynchronization::CoincidenceVector(vector<int> input1, vector<int> input2)
{
	// Check if the two inputs have the same size.
	if (input1.size() != input2.size())
	{
		cout << "Error: the inputs don't have the same size.";
		return vector<int>();
	}

	int trainSize = input1.size();
	vector<int> coincidenceVector = vector<int>(trainSize);

	for (int i = 0; i < trainSize; ++i)
	{
		// If this is not a spike, continue to next element in the vector.
		if (input1[i] == 0)
			continue;

		int jMin = trainSize;
		int minDistance = trainSize;

		for (int j = 0; j < trainSize; ++j)
		{
			// If this is not a spike, continue to next element in the vector.
			if (input2[j] == 0)
				continue;

			// Get the closest spike j (2nd spike train) to the current
			// spike i (1st spike train).
			if (abs(i - j) < minDistance)
			{
				jMin = j;
				minDistance = abs(i - j);
			}
		}

		if (minDistance < getTau(input1, input2, i, jMin))
			coincidenceVector[i] = 1;
		else coincidenceVector[i] = 0;
	}

	return coincidenceVector;
}
