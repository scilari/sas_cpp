#include <iostream>
#include <vector>
#include "src/SystematicAliasSampler.h"


int main() {
    std::cout << "Hello Systematic Alias Sampler!" << std::endl;

    // Setting up possible values and the corresponding probabilities (point mass function)
    std::vector<int> values = {0, 1, 2, 3, 4};
    std::vector<double> pmf = {0.1, 0.3, 0.4, 0.1, 0.1};

    // Constructing the sampler
    scilari::SystematicAliasSampler<int> sampler(values, pmf);

    // Sampling some values and printing the frequencies
    std::vector<int> sampleCounts = {10, 33, 100, 333, 1000, 3333, 10000};

    for(int sampleCount : sampleCounts){
        std::cout << "Sampling " << sampleCount << " values..." << std::endl;

        auto samples = sampler.Sample(sampleCount);

        // counting the occurrences
        std::vector<int> occurrences = {0, 0, 0, 0, 0};
        for(int s: samples){
            occurrences[s]++;
        }

        // printing the empirical ratios
        for(int i = 0; i < occurrences.size(); ++i){
            double ratio = occurrences[i]/ static_cast<double>(sampleCount);
            double expected = pmf[i];
            std::cout << "Expected ratio of " << i << ": " << expected << " empirical ratio: " << ratio << std::endl;
        }

        std::cout << std::endl;
    }

    return 0;
}