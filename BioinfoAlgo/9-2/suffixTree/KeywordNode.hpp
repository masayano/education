#pragma once

#include "typedef.hpp"

#include <string>
#include <vector>

class KeywordNode {
    const char keyword;
    std::vector<POSITION>    indexes;
    std::vector<KeywordNode> childNodes;
public:
    KeywordNode(const char chr);
    void regString(
            const std::string& str,
            const POSITION&    index);
    const KeywordNode& search(
            const KeywordNode& root,
            const std::string& str) const;
    const std::string toString() const;
    /* getter */
    char getKeyword() const;
    const std::vector<POSITION>&    getIndexes   () const;
    const std::vector<KeywordNode>& getChildNodes() const;
};
