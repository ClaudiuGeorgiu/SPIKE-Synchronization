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

int SPIKESynchronization::getPreviousSpikeIndex(vector<int> inputTrain, int index)
{
    // Check the validity of the provided index.
    if (index < 1 || index >= inputTrain.size())
        return -1;

    // Search for the closest previous spike in the input train, starting from
    // the provided index and going backwards.
    for (int n = index - 1; n >= 0; --n)
    {
        if (inputTrain[n] == 1)
            return n;
    }

    // If we are here, no previous spike was found.
    return -1;
}

int SPIKESynchronization::getNextSpikeIndex(vector<int> inputTrain, int index)
{
    // Check the validity of the provided index.
    if (index < 0 || index >= inputTrain.size() - 1)
        return -1;

    // Search for the closest next spike in the input train, starting from
    // the provided index and going forward.
    for (int n = index + 1; n < inputTrain.size(); ++n)
    {
        if (inputTrain[n] == 1)
            return n;
    }

    // If we are here, no next spike was found.
    return -1;
}

double SPIKESynchronization::getTau(vector<int> inputTrain1, vector<int> inputTrain2, int index1, int index2)
{
    // A temporary vector to store the inter-spike intervals.
    vector<int> temp;

    int nextSpike1 = getNextSpikeIndex(inputTrain1, index1);
    if (nextSpike1 != -1)
        temp.push_back(nextSpike1 - index1); // v_i (1)

    int prevSpike1 = getPreviousSpikeIndex(inputTrain1, index1);
    if (prevSpike1 != -1)
        temp.push_back(index1 - prevSpike1); // v_(i-1) (1)

    int nextSpike2 = getNextSpikeIndex(inputTrain2, index2);
    if (nextSpike2 != -1)
        temp.push_back(nextSpike2 - index2); // v_j (2)

    int prevSpike2 = getPreviousSpikeIndex(inputTrain2, index2);
    if (prevSpike2 != -1)
        temp.push_back(index2 - prevSpike2); // v_(j-1) (2)

    // If there are no inter-spike intervals.
    if (temp.size() == 0)
        return 0;

    // Take the minimum inter-spike interval multiplied by 1/2, as described in the paper.
    return 0.5 * (*min_element(begin(temp), end(temp)));
}

vector<int> SPIKESynchronization::CoincidenceVectorPair(vector<int> inputTrain1, vector<int> inputTrain2)
{
    // The spike trains in input can have different sizes,
    // and the coincidence train will have the same size
    // as the longest of the two trains in input.
    int trainSize1 = inputTrain1.size();
    int trainSize2 = inputTrain2.size();
    int longestTrain = std::max(trainSize1, trainSize2);

    vector<int> coincidenceVector(longestTrain, -1);

    for (int i = 0; i < trainSize1; ++i)
    {
        // If this is not a spike, continue to next element in the vector.
        if (inputTrain1[i] != 1)
            continue;

        // The maximum index inside the second train is equal to train's size - 1.
        int jMin = trainSize2 - 1;

        // We set the highest value for the distance between two spikes. This value
        // will be changed if there exist spikes with smaller distance.
        int minDistance = longestTrain;

        for (int j = 0; j < trainSize2; ++j)
        {
            // If this is not a spike, continue to next element in the vector.
            if (inputTrain2[j] != 1)
                continue;

            // Get the closest spike j (2nd spike train) to the current
            // spike i (1st spike train).
            if (abs(i - j) < minDistance)
            {
                jMin = j;
                minDistance = abs(i - j);
            }
        }

        // If the distance between the closest spikes is smaller than the coincidence
        // windows, this is a coincidence.
        if (minDistance < getTau(inputTrain1, inputTrain2, i, jMin))
            coincidenceVector[i] = 1;
        else coincidenceVector[i] = 0;
    }

    return coincidenceVector;
}

vector<int> SPIKESynchronization::MergeCoincidencesPair(vector<int> coincidenceVector1, vector<int> coincidenceVector2)
{
    // Check if the two coincidence vectors have the same size.
    if (coincidenceVector1.size() != coincidenceVector2.size())
    {
        cout << "Error: the coincidence vectors don't have the same size.";
        return vector<int>();
    }

    // The coincidence vectors in input must have the same size,
    // so the train size is equal for both.
    int trainSize = coincidenceVector1.size();

    vector<int> mergedCoincidence(trainSize, -1);

    // Put 1 for coincident spikes.
    // Put 0 for non coincident spikes.
    // Put -1 for non-spikes.
    for (int i = 0; i < trainSize; ++i)
    {
        if (coincidenceVector1[i] == 1 || coincidenceVector2[i] == 1)
            mergedCoincidence[i] = 1;
        else if (coincidenceVector1[i] == 0 || coincidenceVector2[i] == 0)
            mergedCoincidence[i] = 0;
    }

    return mergedCoincidence;
}

vector<vector<double>> SPIKESynchronization::CoincidenceVectorMultivariate(vector<vector<int>> inputTrains)
{
    // Contains coincidence vectors for pairs of spike trains.
    vector<vector<vector<int>>> coincidenceVectorPairs;

    // Generate the coincidence vector for all the pairs of input trains.
    for (int i = 0; i < inputTrains.size(); ++i)
    {
        // Contains the pairs of the i-th input train.
        coincidenceVectorPairs.push_back(vector<vector<int>>());

        for (int j = 0; j < inputTrains.size(); ++j)
        {
            if (i != j)
                coincidenceVectorPairs[i].push_back(CoincidenceVectorPair(inputTrains[i], inputTrains[j]));
        }
    }

    int multivariateCoeff = inputTrains.size() - 1; // N - 1

    // In the multivariate case we can have double numbers.
    vector<vector<double>> coincidenceVectorMultivariate;

    for (int i = 0; i < coincidenceVectorPairs.size(); ++i)
    {
        int maxCoincidenceSize = 0;

        for (int j = 0; j < coincidenceVectorPairs[i].size(); ++j)
        {
            if (coincidenceVectorPairs[i][j].size() > maxCoincidenceSize)
                maxCoincidenceSize = coincidenceVectorPairs[i][j].size();
        }

        // Create an average coincidence vector for each input train.
        coincidenceVectorMultivariate.push_back(vector<double>(maxCoincidenceSize, -1));
    }

    for (int h = 0; h < coincidenceVectorPairs.size(); ++h)
    {
        for (int i = 0; i < coincidenceVectorPairs[h].size(); ++i)
        {
            for (int j = 0; j < coincidenceVectorPairs[h].size(); ++j)
            {
                // Compute the total coincidence counter for each spike in every spike train.
                for (int n = 0; n < coincidenceVectorPairs[h][j].size(); ++n)
                {
                    if (j == 0)
                        coincidenceVectorMultivariate[h][n] = coincidenceVectorPairs[h][j][n];

                    // Skip non-spikes.
                    else if (coincidenceVectorPairs[h][j][n] != -1)
                    {
                        if (coincidenceVectorMultivariate[h][n] == -1)
                            coincidenceVectorMultivariate[h][n] = coincidenceVectorPairs[h][j][n];
                        else coincidenceVectorMultivariate[h][n] += coincidenceVectorPairs[h][j][n];
                    }
                }
            }
        }
    }

    for (int i = 0; i < coincidenceVectorMultivariate.size(); ++i)
    {
        for (int n = 0; n < coincidenceVectorMultivariate[i].size(); ++n)
        {
            // Compute the average coincidence counter for each spike in every spike train.
            if (coincidenceVectorMultivariate[i][n] != -1)
                coincidenceVectorMultivariate[i][n] /= multivariateCoeff;
        }
    }

    return coincidenceVectorMultivariate;
}

vector<double> SPIKESynchronization::MergeCoincidencesMultivariate(vector<vector<double>> coincidenceVectors)
{
    int trainSize = 0;

    // The final coincidence vector will have the length of the longest
    // coincidence vector obtained from a pair of spike trains.
    for (int n = 0; n < coincidenceVectors.size(); ++n)
    {
        if (coincidenceVectors[n].size() > trainSize)
            trainSize = coincidenceVectors[n].size();
    }

    vector<double> mergedCoincidenceMultivariate(trainSize, -1);

    for (int i = 0; i < trainSize; ++i)
    {
        for (int j = 0; j < coincidenceVectors.size(); ++j)
        {
            if (coincidenceVectors[j][i] != -1 && coincidenceVectors[j][i] > mergedCoincidenceMultivariate[i])
                mergedCoincidenceMultivariate[i] = coincidenceVectors[j][i];
        }
    }

    return mergedCoincidenceMultivariate;
}

double SPIKESynchronization::SYNCValue(vector<double> coincidenceProfile)
{
    double syncValue = 0;
    double totalSpikes = 0;

    for (int i = 0; i < coincidenceProfile.size(); ++i)
    {
        if (coincidenceProfile[i] == 0)
        {
            ++totalSpikes;
        }
        else if (coincidenceProfile[i] > 0)
        {
            syncValue += coincidenceProfile[i];
            ++totalSpikes;
        }
    }

    if (totalSpikes == 0)
        return 0;

    return syncValue / totalSpikes;
}

double SPIKESynchronization::SYNCDistance(vector<double> coincidenceProfile)
{
    return 1 - SYNCValue(coincidenceProfile);
}
