#include "SPIKESynchronization.h"
#include <iostream>
#include <vector>
#include <algorithm>

// Used for printing the numbers in the console aligned.
#include <iomanip>

using namespace std;

// Print a vector to the console.
void printVector(vector<int>);
void printVector(vector<double>);

int main(void)
{
    // Some sample spike trains. 
    //  1 - spike
    // -1 - non-spike
    vector<vector<int>> inputTrains =
    {
        { 1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, -1, -1, 1 },
        { -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1 },
        { 1, 1, -1, -1, -1, 1, -1, 1, -1, -1, 1, -1, -1, 1, 1, 1, -1 }
    };

    // Print the sample spike trains.
    for (int n = 0; n < inputTrains.size(); ++n)
    {
        printVector(inputTrains[n]);
    }
    cout << "\n";

    SPIKESynchronization* spike = new SPIKESynchronization();

    vector<double> synchronizationProfile = spike->MergeCoincidencesMultivariate(spike->CoincidenceVectorMultivariate(inputTrains));

    // Print the synchronization profile.
    printVector(synchronizationProfile);

    cout << "\nSYNC value: " << spike->SYNCValue(synchronizationProfile);
    cout << "\nSYNC distance: " << spike->SYNCDistance(synchronizationProfile);

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
        cout << std::right << std::setw(7) << input[n];
    }
    cout << "\n";
}

void printVector(vector<double> input)
{
    for (int n = 0; n < input.size(); ++n)
    {
        cout << std::setprecision(2) << std::right << std::setw(7) << input[n];
    }
    cout << "\n";
}