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
    vector<int> temp;

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
    // The spike trains in input can have different sizes,
    // and the coincidence train will have the same size
    // as the longest of the two trains in input.
    int trainSize1 = input1.size();
    int trainSize2 = input2.size();
    int longestTrain = std::max(trainSize1, trainSize2);

    vector<int> coincidenceVector(longestTrain, -1);

    for (int i = 0; i < trainSize1; ++i)
    {
        // If this is not a spike, continue to next element in the vector.
        if (input1[i] != 1)
            continue;

        // The maximum index inside the second train is equal to train's size - 1.
        int jMin = trainSize2 - 1;

        // We set the highest value for the distance between two spikes. This value
        // will be changed if there exist spikes with closer distance.
        int minDistance = longestTrain;

        for (int j = 0; j < trainSize2; ++j)
        {
            // If this is not a spike, continue to next element in the vector.
            if (input2[j] != 1)
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

vector<int> SPIKESynchronization::MergeCoincidences(vector<int> input1, vector<int> input2)
{
    // Check if the two coincidence vectors have the same size.
    if (input1.size() != input2.size())
    {
        cout << "Error: the coincidence vectors don't have the same size.";
        return vector<int>();
    }

    // The coincidence vectors in input must have the same size,
    // so the train size is equal for both.
    int trainSize = input1.size();

    vector<int> mergedCoincidence(trainSize, -1);

    // Put 1 for coincident spikes.
    // Put 0 for non coincident spikes.
    // Put -1 for non-spikes.
    for (int i = 0; i < trainSize; ++i)
    {
        if (input1[i] == 1 || input2[i] == 1)
            mergedCoincidence[i] = 1;
        else if (input1[i] == 0 || input2[i] == 0)
            mergedCoincidence[i] = 0;
    }

    return mergedCoincidence;
}

double SPIKESynchronization::SYNCValue(vector<int> coincidenceVector)
{
    double syncValue = 0;
    double totalSpikes = 0;

    for (int i = 0; i < coincidenceVector.size(); ++i)
    {
        if (coincidenceVector[i] == 0)
        {
            ++totalSpikes;
        }
        else if (coincidenceVector[i] == 1)
        {
            ++syncValue;
            ++totalSpikes;
        }
    }

    return syncValue / totalSpikes;
}

double SPIKESynchronization::SYNCDistance(vector<int> coincidenceVector)
{
    return 1 - SYNCValue(coincidenceVector);
}

vector<vector<double>> SPIKESynchronization::CoincidenceVectorMultivariate(vector<vector<int>> input)
{
    vector<vector<int>> coincidenceVectorPairs;

    for (int i = 0; i < input.size(); ++i)
    {
        for (int j = i + 1; j < input.size(); ++j)
        {
            coincidenceVectorPairs.push_back(CoincidenceVector(input[i], input[j]));
        }
    }

    int multivariateCoeff = input.size() - 1;
    vector<vector<double>> coincidenceVectorMultivariate;

    for (int i = 0; i < coincidenceVectorPairs.size(); ++i)
    {
        coincidenceVectorMultivariate.push_back(vector<double>(coincidenceVectorPairs[i].size(), -1));
    }

    for (int i = 0; i < coincidenceVectorPairs.size(); ++i)
    {
        for (int j = 0; j < coincidenceVectorPairs.size(); ++j)
        {
            if (i == j)
                continue;

            for (int n = 0; n < coincidenceVectorPairs[i].size(); ++n)
            {
                if (coincidenceVectorPairs[j][n] != -1)
                {
                    if (coincidenceVectorMultivariate[i][n] == -1)
                        coincidenceVectorMultivariate[i][n] = coincidenceVectorPairs[j][n];
                    else coincidenceVectorMultivariate[i][n] += coincidenceVectorPairs[j][n];
                }
            }
        }
    }

    for (int i = 0; i < coincidenceVectorMultivariate.size(); ++i)
    {
        for (int n = 0; n < coincidenceVectorMultivariate[i].size(); ++n)
        {
            if (coincidenceVectorMultivariate[i][n] != -1)
                coincidenceVectorMultivariate[i][n] /= multivariateCoeff;
        }
    }

    return coincidenceVectorMultivariate;
}

vector<double> SPIKESynchronization::MergeCoincidencesMultivariate(vector<vector<double>> input)
{
    int trainSize = 0;

    // The final coincidence vector will have the length of the longest
    // coincidence vector obtained from a pair of spike trains.
    for (int n = 0; n < input.size(); ++n)
    {
        if (input[n].size() > trainSize)
            trainSize = input[n].size();
    }

    vector<double> mergedCoincidenceMultivariate(trainSize, -1);

    for (int i = 0; i < trainSize; ++i)
    {
        for (int j = 0; j < input.size(); ++j)
        {
            if (input[j][i] != -1 && input[j][i] > mergedCoincidenceMultivariate[i])
                mergedCoincidenceMultivariate[i] = input[j][i];
        }
    }

    return mergedCoincidenceMultivariate;
}

double SPIKESynchronization::SYNCValueMultivariate(vector<vector<double>> input)
{
    double syncValue = 0;
    double totalSpikes = 0;

    for (int i = 0; i < input.size(); ++i)
    {
        for (int n = 0; n < input[i].size(); ++n)
        {
            if (input[i][n] == 0)
            {
                ++totalSpikes;
            }
            else if (input[i][n] > 0)
            {
                ++syncValue;
                ++totalSpikes;
            }
        }
    }

    return syncValue / totalSpikes;
}

double SPIKESynchronization::SYNCDistanceMultivariate(vector<vector<double>> input)
{
    return 1 - SYNCValueMultivariate(input);
}
