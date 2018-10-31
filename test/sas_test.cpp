//
// Created by ilari on 29/10/2018.
//
#include "gtest/gtest.h"
#include "../src/AliasTable.h"
#include "../src/SystematicAliasSampler.h"

using namespace scilari;

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

class AliasTableTest : public :: testing::Test{
public:
    std::vector<double> pmf = {0.05, 0.05, 0.1, 0.2, 0.3, 0.2, 0.1};
    AliasTable alias = AliasTable(pmf);
    size_t n = pmf.size();
    double epsilon = 1e-6;
};

class AliasSamplerTest : public AliasTableTest{
public:
    std::vector<int> values = {0, 1, 2, 3, 4, 5, 6};
    SystematicAliasSampler<int> sampler = SystematicAliasSampler<int>(values, pmf);

    std::vector<int> sampleCounts = std::vector<int>(n);
};

TEST_F(AliasTableTest, TableSumsToOriginal){
    std::vector<double> prob(n);
    for(int i = 0; i < n; ++i){
        prob[i] += (1.0 - alias.aliasProbabilities[i])/n;
        prob[alias.aliasIndices[i]] += alias.aliasProbabilities[i]/n;
    }

    for(int i = 0; i < n; ++i){
        EXPECT_NEAR(prob[i], pmf[i], epsilon);
    }
}

TEST_F(AliasSamplerTest, EmpiricalDistributionMatchesOriginal){
    double tolerance = 1e-2;
    int sampleCount = 1000;
    auto samples = sampler.Sample(sampleCount);

    for(int s : samples){
        sampleCounts[s]++;
    }

    for(int v : values){
        double ratio = sampleCounts[v]/ static_cast<double>(sampleCount);
        double expected = pmf[v];
        EXPECT_NEAR(ratio, expected, tolerance);
    }
}

TEST_F(AliasSamplerTest, SystematicMatchesBetterThanIid){
    int sampleCount = 1000;

    std::vector<int> samplesIid;
    while(samplesIid.size() < sampleCount){
        samplesIid.push_back(sampler.Sample());
    }

    auto samples = sampler.Sample(sampleCount);

    for(int s : samples){
        sampleCounts[s]++;
    }

    std::vector<int> sampleCountsIid = std::vector<int>(n);
    for(int s : samplesIid){
        sampleCountsIid[s]++;
    }

    // There might be some random fails in individual values due to randomness
    int failCount = 0;
    double maxDiff = 0;
    double maxDiffIid = 0;
    for(int v : values){
        double ratio = sampleCounts[v]/ static_cast<double>(sampleCount);
        double ratioIid = sampleCountsIid[v] / static_cast<double>(sampleCount);

        double expected = pmf[v];

        double diff = abs(ratio - expected);
        double diffIid = abs(ratioIid - expected);

        if(diff > diffIid) failCount++;

        maxDiff = fmax(maxDiff, diff);
        maxDiffIid = fmax(maxDiffIid, diffIid);
    }

    EXPECT_LE(failCount, 2);
    EXPECT_LE(maxDiff, maxDiffIid);
}

