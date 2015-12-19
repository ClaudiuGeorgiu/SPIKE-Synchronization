#include <vector>
#include <map>

#ifndef SPIKESYNCHRONIZATION_H
#define SPIKESYNCHRONIZATION_H

class SPIKESynchronization
{
    protected:

        // Get the index of the previous spike in the input vector, starting from the provided index.
        // Return -1 if no valid index was found.
        int getPreviousSpikeIndex(std::vector<int> inputTrain, int index);
        int getPreviousSpikeIndex(std::vector<double> inputTrain, int index);

        // Get the index of the next spike in the input vector, starting from the provided index.
        // Return -1 if no valid index was found.
        int getNextSpikeIndex(std::vector<int> inputTrain, int index);
        int getNextSpikeIndex(std::vector<double> inputTrain, int index);
        
        // Get the coincidence window from the inputs and the indices, as described in the paper.
        double getTau(std::vector<int> inputTrain1, std::vector<int> inputTrain2, int index1, int index2);
        double getTau(std::vector<double> inputTrain1, std::vector<double> inputTrain2, int index1, int index2);

        // Get the vector containing the coincidence indices for a pair of spike trains.
        std::vector<int> CoincidenceVectorPair(std::vector<int> inputTrain1, std::vector<int> inputTrain2);
        std::map<double, int> CoincidenceVectorPair(std::vector<double> inputTrain1, std::vector<double> inputTrain2);

    public:

        SPIKESynchronization();
        virtual ~SPIKESynchronization();

        // Get a list of vectors containing the coincidence indices for the pairs of spike trains in input.
        std::vector<std::vector<double>> CoincidenceVectorMultivariate(std::vector<std::vector<int>> inputTrainsVector);
        std::vector<std::map<double, double>> CoincidenceVectorMultivariate(std::vector<std::vector<double>> inputTrainsTime);

        // Get the SPIKE-Synchronization profile by merging all the coincidence vectors of all the spike trains.
        std::vector<double> MergeCoincidencesMultivariate(std::vector<std::vector<double>> coincidenceVectorsVector);
        std::map<double, double> MergeCoincidencesMultivariate(std::vector<std::map<double, double>> coincidenceVectorsTime);

        double SYNCValue(std::vector<double> coincidenceProfile);
        double SYNCValue(std::map<double, double> coincidenceProfile);

        double SYNCDistance(std::vector<double> coincidenceProfile);
        double SYNCDistance(std::map<double, double> coincidenceProfile);
};

#endif