#pragma once

#include "KeywordTree.hpp"
#include "SuffixNode.hpp"

class SuffixTree {
    SuffixNode root;
public:
    SuffixTree(const KeywordTree& kt);
    const std::string search(const std::string& str) const;
};
