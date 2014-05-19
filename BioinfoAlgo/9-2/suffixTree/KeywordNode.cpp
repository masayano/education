#include "KeywordNode.hpp"

KeywordNode::KeywordNode(const char chr) : keyword(chr) {}

NODE_INDEX getChildIdx(
        const std::vector<KeywordNode>& childNodes,
        const char chr) {
    const int length = childNodes.size();
    for(int i = 0; i < length; ++i) {
        if(childNodes[i].getKeyword() == chr) {
            return i;
        }
    }
    return NOT_FOUND;
}

const POSITION nextIndex(const POSITION& index) {
    auto next = index;
    next.first += 1;
    return next;
}

void KeywordNode::regString(
        const std::string& str,
        const POSITION& index) {
    if(!str.empty()) {
        indexes.emplace_back(index);
        const auto chr = str[0];
        auto nodeIdx = getChildIdx(childNodes, chr);
        if(nodeIdx == NOT_FOUND) {
            childNodes.emplace_back(chr);
            nodeIdx = childNodes.size()-1;
        }
        childNodes[nodeIdx].regString(str.substr(1), nextIndex(index));
    }
}

#include <iostream>

const KeywordNode& KeywordNode::search(
        const KeywordNode& root,
        const std::string& str) const {
    std::cout << keyword << std::endl;
    if(str.empty()) {
        return *this;
    } else {
        const auto chr = str[0];
        const auto nodeIdx = getChildIdx(childNodes, chr);
        if(nodeIdx == NOT_FOUND) {
            return root;
        } else {
            return childNodes[nodeIdx].search(root, str.substr(1));
        }
    }
}


#include <sstream>

const std::string KeywordNode::toString() const {
    std::stringstream ss;
    for(const auto index : indexes) {
        ss << "idx: " << index.first << " of seq:" << index.second << std::endl;
    }
    ss << std::endl;
    return ss.str();
}

    /* getter */

char KeywordNode::getKeyword() const { return keyword; }
const std::vector<POSITION>&    KeywordNode::getIndexes   () const { return indexes; }
const std::vector<KeywordNode>& KeywordNode::getChildNodes() const { return childNodes; }
