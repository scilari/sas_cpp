//
// Created by ilari on 30/10/2018.
//

#ifndef SAS_CPP_ALIASTABLE_H
#define SAS_CPP_ALIASTABLE_H

#include <vector>
/**
  * Kronmal and Peterson 1979 implementation of alias method by Walker (1974, 1977) as described in
  * "An Analysis of the Alias Method for Discrete Random-Variate Generation" by Smith and Jacobson (2005)
  *
  * After initialization with the point mass function (pmf), the object acts as a wrapper for the aliased
  * value indices and their corresponding probabilities (ratio of the bin).
  */
namespace scilari{

    class AliasTable {
    public:
        /// Constructs alias table by using the given point mass function (pmf). Pmf does not need to be normalized.
        explicit AliasTable(const std::vector<double>& pmf);

        /// Indices of elements corresponding to the probability mass in aliasProbabilities
        std::vector<int> aliasIndices;

        /// Probabilities corresponding to the aliased values (represented by the indices)
        std::vector<double> aliasProbabilities;
    };

} // namespace scilari

#endif //SAS_CPP_ALIASTABLE_H
