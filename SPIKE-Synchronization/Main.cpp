#include "SPIKESynchronization.h"
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

void printVector(vector<int>);
void printVector(vector<double>);

int main(void)
{
    vector<int> input1 = { 1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, -1, -1, 1 };
    vector<int> input2 = { -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1 };

    printVector(input1);
    printVector(input2);
    cout << "\n";

    SPIKESynchronization* spike = new SPIKESynchronization();

//     vector<int> coincidence = spike->MergeCoincidences(spike->CoincidenceVector(input1, input2),
//                                                        spike->CoincidenceVector(input2, input1));
// 
//     printVector(coincidence);
// 
//     cout << "\n" << spike->SYNCValue(coincidence);
//     cout << "\n" << spike->SYNCDistance(coincidence);

    

    vector<vector<int>> multivariateInput;
    multivariateInput.push_back(input1);
    multivariateInput.push_back(input2);

    auto temp = spike->CoincidenceVectorMultivariate(multivariateInput);
    auto coincidenceMultivariate = spike->MergeCoincidencesMultivariate(temp);
    printVector(coincidenceMultivariate);
    cout << "\n" << spike->SYNCValueMultivariate(temp);
    cout << "\n" << spike->SYNCDistanceMultivariate(temp);

    if (spike != NULL)
        delete spike;

    getchar();

    return 0;
}

void printVector(vector<int> input)
{
    for (int n = 0; n < input.size(); ++n)
    {
        cout << std::right << std::setw(5) << input[n];
    }
    cout << "\n";
}

void printVector(vector<double> input)
{
    for (int n = 0; n < input.size(); ++n)
    {
        cout << std::right << std::setw(5) << input[n];
    }
    cout << "\n";
}