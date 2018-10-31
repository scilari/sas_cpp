//
// Created by ilari vallivaara on 25/10/2018.
//

#ifndef SAS_CPP_SYSTEMATICALIASSAMPLER_H
#define SAS_CPP_SYSTEMATICALIASSAMPLER_H

#include <cassert>
#include <chrono>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include "AliasTable.h"

namespace scilari{

    /**
     * Systematic Alias Sampler described fully in the corresponding paper (TOMS, ACM, Vallivaara et. al). Allows
     * efficient batch sampling, and produces empirical distributions with good fit.
     * @tparam E Element type of the values returned by the sampler.
     */
    template <class E> class SystematicAliasSampler {

    public:
        /**
         * Constructs the sampler
         * @param values Values corresponding to the probability masses (pmf)
         * @param pmf Probability masses corresponding to the values.
         */
        SystematicAliasSampler(const std::vector<E> &values, const std::vector<double> &pmf)
                : values_(values), aliasedValues_(values.size()),
                  binCount_(values.size()) {

            assert(values.size() == pmf.size());
            InitRNGs();
            InitializeAliasTable(NormalizedProbabilities(pmf));
        }

        /**
         * Returns a single i.i.d. sample.
         * @return I.i.d. sample.
         */
        E Sample() { return Sample(RandomInt(), RandomDouble()); }

        /**
         * Returns a non-i.i.d. batch of samples from the pmf given in constructor.
         * @param sampleCount Number of samples
         * @return Vector containing the samples (partially ordered).
         */
        std::vector<E> Sample(int sampleCount) {
            std::vector<E> v(sampleCount);
            return Sample(sampleCount, v);
        }

        /**
         * Returns a non-i.i.d. batch of samples from the pmf given in constructor.
         * @param sampleCount Number of samples.
         * @param output Pre-initialized vector of size sampleCount used as the return value.
         * @return Vector containing the samples (partially ordered).
         */
        std::vector<E> Sample(int sampleCount, std::vector<E> &output) { return SampleSystematic(sampleCount, output, 0); }


    private:
        unsigned int binCount_;
        std::vector<E> values_;
        std::vector<E> aliasedValues_;
        std::vector<double> aliasProbabilities_;

        // constants used as parameters
        constexpr static double divEpsilon_ = 0.07;
        constexpr static int minBatchSize_ = 16;
        constexpr static int minRecurSize_ = 4 * minBatchSize_;
        constexpr static int batchSplitNumerator = 7;
        constexpr static int batchSplitDenominator = 13;

        // RNGs
        std::mt19937 randomGen;
        std::uniform_real_distribution<double> randomDouble_;
        std::uniform_int_distribution<int> randomInt_;

        std::vector<E> SampleSystematic(int sampleCount, std::vector<E> &samples,
                                        int fillFrom) {
            if (sampleCount > minBatchSize_ &&
                IsDivisibilityProblem(binCount_, sampleCount)) {
                int splitIndex =
                        sampleCount <= minRecurSize_
                        ? minBatchSize_
                        : sampleCount * batchSplitNumerator / batchSplitDenominator;

                SampleSystematic(sampleCount - splitIndex, samples, fillFrom);
                return SampleSystematic(splitIndex, samples,
                                        fillFrom + sampleCount - splitIndex);
            } else {
                double step = binCount_ / static_cast<double>(sampleCount);
                double r = step * (1.0 - RandomDouble()); // (0, step]

                // going backwards in order to prevent index out of bounds error
                // (this is fine, because sample(-0, -epsilon) = sample(0, 0))
                double x = binCount_ - r;
                double i = fillFrom;
                int fillTo = fillFrom + sampleCount;

                while (i < fillTo) {
                    double ri = floor(x);
                    double rf = x - ri;
                    samples[i] = Sample(ri, rf);
                    x -= step;
                    ++i;
                }

                return samples;
            }
        }

        E Sample(int intPart, double fracPart) {
            if (fracPart > aliasProbabilities_[intPart])
                return values_[intPart];
            else
                return aliasedValues_[intPart];
        }

        double RandomDouble() { return randomDouble_(randomGen); }

        int RandomInt() { return randomInt_(randomGen); }

        bool IsDivisibilityProblem(int binCount, int sampleCount) {
            return AlmostDivides(binCount, sampleCount, divEpsilon_) ||
                   AlmostDivides(binCount * 4, sampleCount, divEpsilon_) ||
                   AlmostDivides(binCount * 5, sampleCount, divEpsilon_) ||
                   AlmostDivides(binCount * 6, sampleCount, divEpsilon_);
        }

        bool AlmostDivides(double x, double y, double eps) {
            double q = x / y;
            double d = DistanceFromInt(q);
            return d < eps;
        }

        double DistanceFromInt(double x) {
            double y = x - std::floor(x);
            return std::min(y, 1.0 - y);
        }

        void InitRNGs() {
            auto seed =
                    std::chrono::high_resolution_clock::now().time_since_epoch().count();
            randomGen = std::mt19937(seed);
            randomDouble_ = std::uniform_real_distribution<>(0.0, 1.0);
            randomInt_ = std::uniform_int_distribution<>(0, binCount_ - 1);
        }

        std::vector<double> NormalizedProbabilities(const std::vector<double> &v) {
            double sum = 0.0;
            for (double x : v) {
                if (x < 0.0) {
                    throw std::invalid_argument("Probabilities cannot be negative.");
                }
                sum += x;
            }

            if (sum == 0.0) {
                throw std::invalid_argument("Probabilities cannot be all zero.");
            }

            double invSum = 1.0 / sum;

            std::vector<double> normalized;
            normalized.reserve(v.size());
            for (double x : v) {
                normalized.push_back(x * invSum);
            }

            return normalized;
        }

        void InitializeAliasTable(const std::vector<double> &pmf) {
            auto alias = AliasTable(pmf);
            // populate aliasedValues using alias indices
            for (int i = 0; i < binCount_; ++i) {
                int ix = alias.aliasIndices[i];
                aliasedValues_[i] = values_[ix];
            }

            aliasProbabilities_ = std::move(alias.aliasProbabilities);
        }
    };

} // namespace scilari

#endif // SAS_CPP_SYSTEMATICALIASSAMPLER_H
