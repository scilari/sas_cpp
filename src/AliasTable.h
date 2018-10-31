//
// Created by ilari on 30/10/2018.
//

#ifndef SAS_CPP_ALIASTABLE_H
#define SAS_CPP_ALIASTABLE_H


#include <vector>

namespace scilari{

    class AliasTable {
    public:
        AliasTable(const std::vector<double>& pmf);

        std::vector<int> aliasIndices;
        std::vector<double> aliasProbabilities;
    };
} // namespace scilari




#endif //SAS_CPP_ALIASTABLE_H
