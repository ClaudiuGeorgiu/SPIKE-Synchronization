#include "SPIKESynchronization.h"

// Libraries used for the data types in input.
#include <vector>
#include <map>

// Used to get the maximum value of a double number.
#include <cfloat>

// Used to get the minimum element in a vector.
#include <algorithm>

using namespace std;

SPIKESynchronization::SPIKESynchronization()
{
}

SPIKESynchronization::~SPIKESynchronization()
{
}



/*******************************************************************************************************************************/
/* Used only for vector inputs containing 1 where spikes occur, -1 otherwise.                                                  */
/*******************************************************************************************************************************/

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
    if (temp.empty())
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
        double minDistance = DBL_MAX;

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
        // window, this is a coincidence.
        if (minDistance < getTau(inputTrain1, inputTrain2, i, jMin))
            coincidenceVector[i] = 1;
        else coincidenceVector[i] = 0;
    }

    return coincidenceVector;
}

vector<vector<double>> SPIKESynchronization::CoincidenceVectorMultivariate(vector<vector<int>> inputTrainsVector)
{
    // Contains coincidence vectors for pairs of spike trains.
    vector<vector<vector<int>>> coincidenceVectorPairs;

    // Generate the coincidence vector for all the pairs of input trains.
    for (int i = 0; i < inputTrainsVector.size(); ++i)
    {
        // Contains the pairs of the i-th input train.
        coincidenceVectorPairs.push_back(vector<vector<int>>());

        for (int j = 0; j < inputTrainsVector.size(); ++j)
        {
            if (i != j)
                coincidenceVectorPairs[i].push_back(CoincidenceVectorPair(inputTrainsVector[i], inputTrainsVector[j]));
        }
    }

    int multivariateCoeff = inputTrainsVector.size() - 1; // N - 1

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

vector<double> SPIKESynchronization::MergeCoincidencesMultivariate(vector<vector<double>> coincidenceVectorsVector)
{
    int trainSize = 0;

    // The final coincidence vector will have the length of the longest
    // coincidence vector obtained from a pair of spike trains.
    for (int n = 0; n < coincidenceVectorsVector.size(); ++n)
    {
        if (coincidenceVectorsVector[n].size() > trainSize)
            trainSize = coincidenceVectorsVector[n].size();
    }

    vector<double> mergedCoincidenceMultivariate(trainSize, -1);

    for (int i = 0; i < trainSize; ++i)
    {
        for (int j = 0; j < coincidenceVectorsVector.size(); ++j)
        {
            if (coincidenceVectorsVector[j][i] != -1 && coincidenceVectorsVector[j][i] > mergedCoincidenceMultivariate[i])
                mergedCoincidenceMultivariate[i] = coincidenceVectorsVector[j][i];
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

/*******************************************************************************************************************************/



/*******************************************************************************************************************************/
/* Used only for vector inputs containing the times at which the spikes occur.                                                 */
/*******************************************************************************************************************************/

int SPIKESynchronization::getPreviousSpikeIndex(vector<double> inputTrain, int index)
{
    // Check the validity of the provided index.
    if (index < 1 || index >= inputTrain.size())
        return -1;

    return index - 1;
}

int SPIKESynchronization::getNextSpikeIndex(vector<double> inputTrain, int index)
{
    // Check the validity of the provided index.
    if (index < 0 || index >= inputTrain.size() - 1)
        return -1;

    return index + 1;
}

double SPIKESynchronization::getTau(vector<double> inputTrain1, vector<double> inputTrain2, int index1, int index2)
{
    // A temporary vector to store the inter-spike intervals.
    vector<double> temp;

    int nextSpike1 = getNextSpikeIndex(inputTrain1, index1);
    if (nextSpike1 != -1)
        temp.push_back(inputTrain1[nextSpike1] - inputTrain1[index1]); // v_i (1)

    int prevSpike1 = getPreviousSpikeIndex(inputTrain1, index1);
    if (prevSpike1 != -1)
        temp.push_back(inputTrain1[index1] - inputTrain1[prevSpike1]); // v_(i-1) (1)

    int nextSpike2 = getNextSpikeIndex(inputTrain2, index2);
    if (nextSpike2 != -1)
        temp.push_back(inputTrain2[nextSpike2] - inputTrain2[index2]); // v_j (2)

    int prevSpike2 = getPreviousSpikeIndex(inputTrain2, index2);
    if (prevSpike2 != -1)
        temp.push_back(inputTrain2[index2] - inputTrain2[prevSpike2]); // v_(j-1) (2)

    // If there are no inter-spike intervals.
    if (temp.empty())
        return 0;

    // Take the minimum inter-spike interval multiplied by 1/2, as described in the paper.
    return 0.5 * (*min_element(begin(temp), end(temp)));
}

map<double, int> SPIKESynchronization::CoincidenceVectorPair(vector<double> inputTrain1, vector<double> inputTrain2)
{
    int trainSize1 = inputTrain1.size();
    int trainSize2 = inputTrain2.size();

    map<double, int> coincidenceVector;

    // Merge the times contained in the inputs for the coincidence vector.
    for (int j = 0; j < trainSize1; ++j)
    {
        coincidenceVector[inputTrain1[j]] = 0;
    }
    for (int j = 0; j < trainSize2; ++j)
    {
        coincidenceVector[inputTrain2[j]] = 0;
    }

    for (int i = 0; i < trainSize1; ++i)
    {
        // If the input train is ordered, the furthest spike is the last in the train.
        int jMin = inputTrain2.size() - 1;

        // We set the highest value for the distance between two spikes. This value
        // will be changed if there exist spikes with smaller distance.
        double minDistance = DBL_MAX;

        for (int j = 0; j < trainSize2; ++j)
        {
            // Get the closest spike j (2nd spike train) to the current
            // spike i (1st spike train).
            if (abs(inputTrain1[i] - inputTrain2[j]) < minDistance)
            {
                jMin = j;
                minDistance = abs(inputTrain1[i] - inputTrain2[j]);
            }
        }

        // If the distance between the closest spikes is smaller than the coincidence
        // window, this is a coincidence.
        if (minDistance < getTau(inputTrain1, inputTrain2, i, jMin))
            coincidenceVector[inputTrain1[i]] = 1;
        else coincidenceVector[inputTrain1[i]] = 0;
    }

    return coincidenceVector;
}

vector<map<double, double>> SPIKESynchronization::CoincidenceVectorMultivariate(vector<vector<double>> inputTrainsTime)
{
    // Contains coincidence vectors for pairs of spike trains.
    vector<vector<map<double, int>>> coincidenceVectorPairs;

    // Generate the coincidence vector for all the pairs of input trains.
    for (int i = 0; i < inputTrainsTime.size(); ++i)
    {
        // Contains the pairs of the i-th input train.
        coincidenceVectorPairs.push_back(vector<map<double, int>>());

        for (int j = 0; j < inputTrainsTime.size(); ++j)
        {
            if (i != j)
                coincidenceVectorPairs[i].push_back(CoincidenceVectorPair(inputTrainsTime[i], inputTrainsTime[j]));
        }
    }

    int multivariateCoeff = inputTrainsTime.size() - 1; // N - 1
 
    vector<map<double, double>> coincidenceVectorMultivariate;
 
    for (int h = 0; h < coincidenceVectorPairs.size(); ++h)
    {
        coincidenceVectorMultivariate.push_back(map<double, double>());

        for (int i = 0; i < coincidenceVectorPairs[h].size(); ++i)
        {
            // Compute the total coincidence counter for each spike in every spike train.
            for (int j = 0; j < coincidenceVectorPairs[h].size(); ++j)
            {
                if (j == 0)
                {
                    for (auto const &timeSpikePair : coincidenceVectorPairs[h][j])
                    {
                        coincidenceVectorMultivariate[h][timeSpikePair.first] = timeSpikePair.second;
                    }
                }
                else
                {
                    for (auto const &timeSpikePair : coincidenceVectorPairs[h][j])
                    {
                        coincidenceVectorMultivariate[h][timeSpikePair.first] += timeSpikePair.second;
                    }
                }
            }
        }
    }

    for (int i = 0; i < coincidenceVectorMultivariate.size(); ++i)
    {
        for (auto const &timeSpikePair : coincidenceVectorMultivariate[i])
        {
            coincidenceVectorMultivariate[i][timeSpikePair.first] /= multivariateCoeff;
        }
    }

    return coincidenceVectorMultivariate;
}

map<double, double> SPIKESynchronization::MergeCoincidencesMultivariate(vector<map<double, double>> coincidenceVectorsTime)
{
    map<double, double> mergedCoincidenceMultivariate;

    for (auto const &timeSpikePair1 : coincidenceVectorsTime[0])
    {
        // The coincidence vector of the first input is taken as a starting point for the final coincidence.
        mergedCoincidenceMultivariate[timeSpikePair1.first] = timeSpikePair1.second;

        // Iterate over all the coincidences at a given time.
        for (int j = 1; j < coincidenceVectorsTime.size(); ++j)
        {
            if (coincidenceVectorsTime[j].count(timeSpikePair1.first) > 0)
            {
                // Take the coincidence with the highest value.
                if (coincidenceVectorsTime[j][timeSpikePair1.first] > timeSpikePair1.second)
                    mergedCoincidenceMultivariate[timeSpikePair1.first] = coincidenceVectorsTime[j][timeSpikePair1.first];
            }
        }
    }

    return mergedCoincidenceMultivariate;
}

double SPIKESynchronization::SYNCValue(map<double, double> coincidenceProfile)
{
    double syncValue = 0;
    double totalSpikes = 0;

    for (auto const &timeSpikePair : coincidenceProfile)
    {
        if (timeSpikePair.second == 0)
        {
            ++totalSpikes;
        }
        else if (timeSpikePair.second > 0)
        {
            syncValue += timeSpikePair.second;
            ++totalSpikes;
        }
    }

    if (totalSpikes == 0)
        return 0;

    return syncValue / totalSpikes;
}

double SPIKESynchronization::SYNCDistance(map<double, double> coincidenceProfile)
{
    return 1 - SYNCValue(coincidenceProfile);
}

/*******************************************************************************************************************************/