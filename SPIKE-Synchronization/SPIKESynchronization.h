#include <vector>

#ifndef SPIKESYNCHRONIZATION_H
#define SPIKESYNCHRONIZATION_H

class SPIKESynchronization
{
    public:
        SPIKESynchronization();
        virtual ~SPIKESynchronization();

        int getPreviousSpikeIndex(std::vector<int> input, int index);
        int getNextSpikeIndex(std::vector<int> input, int index);
        double getTau(std::vector<int> input1, std::vector<int> input2, int index1, int index2);
        std::vector<int> CoincidenceVector(std::vector<int> input1, std::vector<int> input2);
        std::vector<int> MergeCoincidences(std::vector<int> input1, std::vector<int> input2);
        double SYNCValue(std::vector<int> coincidenceVector);
        double SYNCDistance(std::vector<int> coincidenceVector);
};

#endif