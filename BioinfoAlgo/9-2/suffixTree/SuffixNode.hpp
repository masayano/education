#pragma once

#include "typedef.hpp"

#include "KeywordNode.hpp"

class SuffixNode {
    std::string suffix;
    std::vector<POSITION>   indexes;
    std::vector<SuffixNode> childNodes;
public:
    SuffixNode(const KeywordNode& kn);
    const SuffixNode& search(
            const SuffixNode& root,
            const std::string& str) const;
    const std::string toString() const;
    /* getter */
    const std::string& getSuffix() const;
};
