//
// Created by ilari on 30/10/2018.
//

#include "AliasTable.h"

#include <stack>

using namespace scilari;

AliasTable::AliasTable(const std::vector<double> &pmf) {
    size_t binCount_ = pmf.size();
    // creating a scaled vector where 1.0 corresponds to original 1.0/n (even
    // distribution value)
    std::vector<double> q;
    q.reserve(binCount_);
    for (double p : pmf) {
        q.push_back(p * binCount_);
    }

    // initialize alias indices with some sane index values
    std::vector<int> a(binCount_);
    for (int i = 0; i < binCount_; ++i) {
        a[i] = i;
    }

    std::vector<double> s(binCount_); // working copy of alias probabilities

    // creating two stacks for indices of smaller and larger elements
    // (w.r.t. 1.0)
    std::stack<int> g, h;
    for (int i = 0; i < binCount_; ++i) {
        double p = q[i];
        if (p >= 1.0)
            g.push(i);
        else
            h.push(i);
    }

    // Compute alias indices and probabilities
    // checking both because of possible numerical inaccuracies large may be
    // emptied first
    while (!h.empty() && !g.empty()) {
        // selecting a large and small bin
        int j = h.top();
        int k = g.top();

        h.pop(); // removing the small bin index as it will be filled

        // moving probability mass from large to small bin to fill it up, set
        // probability, and update the indices
        a[j] = k;
        s[j] = 1.0 - q[j];
        q[k] = (q[k] + q[j]) - 1.0; // numerically more stable

        // move large bin index to small bin stack if needed
        if (q[k] < 1.0) {
            g.pop();
            h.push(k);
        }
    }

    aliasIndices = std::move(a);
    aliasProbabilities = std::move(s);

}


