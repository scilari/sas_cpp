#include <iostream>
#include <vector>
#include "src/SystematicAliasSampler.h"


int main() {
    std::vector<int> values = {1, 2, 3, 4};
    std::vector<double> pmf = {0.1, 0.3, 0.4, 0.2};

    scilari::SystematicAliasSampler<int> sampler(values, pmf);

    auto samples = sampler.Sample(100);
    int f1 = 0, f2 = 0, f3 = 0, f4 = 0;
    for(int i : samples){
        if(i == 1) f1++;
        if(i == 2) f2++;
        if(i == 3) f3++;
        if(i == 4) f4++;
    }

    double size = (double) samples.size();
    double q1 = f1/size;
    double q2 = f2/size;
    double q3 = f3/size;
    double q4 = f4/size;

    std::cout << std::to_string(q1) << ": " << std::to_string(pmf[0])  << std::endl;
    std::cout << std::to_string(q2) << ": " << std::to_string(pmf[1])  << std::endl;
    std::cout << std::to_string(q3) << ": " << std::to_string(pmf[2])  << std::endl;
    std::cout << std::to_string(q4) << ": " << std::to_string(pmf[3])  << std::endl;

    return 0;
}