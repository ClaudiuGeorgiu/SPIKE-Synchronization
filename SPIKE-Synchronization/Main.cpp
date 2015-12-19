#include "SPIKESynchronization.h"
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

// Used for printing aligned numbers in the console.
#include <iomanip>

using namespace std;

// Print a vector to the console.
void printVector(vector<int>);
void printVector(vector<double>);

int main(void)
{
    SPIKESynchronization* spike = new SPIKESynchronization();



    /*************************************************************************************************************************************/
    /* Compute the SPIKE-Synchronization profile for vector inputs containing 1 where spikes occur, -1 otherwise.                        */
    /*************************************************************************************************************************************/

    // Some sample spike trains as vectors. 
    //  1 - spike
    // -1 - non-spike
    vector<vector<int>> inputTrainsVector =
    {
        { 1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, -1, -1, 1 },
        { -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1 },
        { 1, 1, -1, -1, -1, 1, -1, 1, -1, -1, 1, -1, -1, 1, 1, 1, -1 }
    };

    cout << "Input trains:\n";

    // Print the sample spike trains.
    for (int n = 0; n < inputTrainsVector.size(); ++n)
    {
        printVector(inputTrainsVector[n]);
    }
    cout << "\n";

    vector<double> synchronizationProfile = spike->MergeCoincidencesMultivariate(spike->CoincidenceVectorMultivariate(inputTrainsVector));

    cout << "SPIKE-Synchronization profile:\n";
    printVector(synchronizationProfile);

    cout << "\nSYNC value: " << spike->SYNCValue(synchronizationProfile);
    cout << "\nSYNC distance: " << spike->SYNCDistance(synchronizationProfile);

    /*************************************************************************************************************************************/



    /*************************************************************************************************************************************/
    /* Compute the SPIKE-Synchronization profile for vector inputs containing only the times at which the spikes occur.                  */
    /*************************************************************************************************************************************/

    // Some sample spike trains (containing the times at which the spikes occur). 
    vector<vector<double>> inputTrainsTime =
    {
        { 0.0, 1.0, 1.5, 3.0, 5.0, 6.5, 7.0, 7.5, 8.0, 10.0, 100.0 },
        { 0.2, 0.4, 1.2, 3.1, 5.0, 6.5, 6.9, 7.9, 8.0, 11.0, 102.0 },
        { 0.0, 1.0, 1.5, 3.0, 5.0, 6.5, 7.0, 7.5, 8.0, 10.0, 100.0 }
    };

    cout << "\n\n\n\n\n\nInput trains times:\n";

    // Print the sample spike trains.
    for (int n = 0; n < inputTrainsTime.size(); ++n)
    {
        printVector(inputTrainsTime[n]);
    }
    cout << "\n";

    map<double, double> synchronizationProfileTime = spike->MergeCoincidencesMultivariate(spike->CoincidenceVectorMultivariate(inputTrainsTime));

 
    cout << "SPIKE-Synchronization profile:\n";
    for (auto const &timeSpikePair : synchronizationProfileTime)
    {
        cout << right << setw(10) << timeSpikePair.second;
    }
    cout << "\n";
    for (auto const &timeSpikePair : synchronizationProfileTime)
    {
        cout << right << setw(10) << timeSpikePair.first;
    }
    cout << "\n";

    cout << "\nSYNC value: " << spike->SYNCValue(synchronizationProfileTime);
    cout << "\nSYNC distance: " << spike->SYNCDistance(synchronizationProfileTime);

    /*************************************************************************************************************************************/



    if (spike != NULL)
        delete spike;

    // Keep showing the console until the user hits a button.
    getchar();

    return 0;
}

void printVector(vector<int> input)
{
    for (int n = 0; n < input.size(); ++n)
    {
        cout << right << setw(10) << input[n];
    }
    cout << "\n";
}

void printVector(vector<double> input)
{
    for (int n = 0; n < input.size(); ++n)
    {
        cout  << right << setw(10) << input[n];
    }
    cout << "\n";
}