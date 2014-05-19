#pragma once

#include "typedef.hpp"

#include <string>
#include "KeywordNode.hpp"

class KeywordTree {
    KeywordNode root;
    STR_NUM     strCount;
public:
    KeywordTree();
    void regString(const std::string& str);
    const std::string search(const std::string& str) const;
    /* getter */
    const KeywordNode& getRoot() const;
};
