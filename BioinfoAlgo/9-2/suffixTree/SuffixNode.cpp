#include "SuffixNode.hpp"

std::pair<std::string, std::vector<KeywordNode> > makeSuffix(const KeywordNode& kn) {
    std::string suffix = "";
    auto* pTempNode   = &kn;
    auto* pChildNodes = &(pTempNode->getChildNodes());
    while(1) {
        suffix += pTempNode->getKeyword();
        if(pChildNodes->size() == 1) {
            pTempNode  = &((*pChildNodes)[0]);
            pChildNodes = &(pTempNode->getChildNodes());
        } else {
            break;
        }
    }
    return std::make_pair(suffix, *pChildNodes);
}

SuffixNode::SuffixNode(const KeywordNode& kn) {
    indexes = kn.getIndexes();
    const auto suffix_childNodes = makeSuffix(kn);
    suffix = suffix_childNodes.first;
    for(const auto& keywordNode : suffix_childNodes.second) {
        childNodes.emplace_back(keywordNode);
    }
}

#include <iostream>


int getMatchLength(
        const std::string& str,
        const std::string& suffix) {
    int end = 0;
    const auto length = std::min(str.size(), suffix.size());
    for(auto i = 0U; i < length; ++i) {
        if(str[i] == suffix[i]) {
            end++;
        } else {
            break;
        }
    }
    return end;
}

NODE_INDEX getChildIdx(
        const std::vector<SuffixNode>& childNodes,
        const std::string& nextStr) {
    const int length = childNodes.size();
    for(int i = 0; i < length; ++i) {
        if(nextStr[0] == childNodes[i].getSuffix()[0]) {
            return i;
        }
    }
    return NOT_FOUND;
}

const SuffixNode& SuffixNode::search(
        const SuffixNode& root,
        const std::string& str) const {
    std::cout << suffix << std::endl;
    if(str.empty() || (suffix.find(str) == 0)) {
        return *this;
    } else {
        const auto matchLength = getMatchLength(str, suffix);
        const auto nextStr = str.substr(matchLength);
        const auto nodeIdx = getChildIdx(childNodes, nextStr);
        if(nodeIdx == NOT_FOUND) {
            return root;
        } else {
            return childNodes[nodeIdx].search(root, nextStr);
        }
    }
}

#include <sstream>

const std::string SuffixNode::toString() const {
    std::stringstream ss;
    for(const auto& index : indexes) {
        ss << "idx: " << index.first << " of seq:" << index.second << std::endl;
    }
    ss << std::endl;
    return ss.str();
}

    /* getter */

const std::string& SuffixNode::getSuffix() const {
    return suffix;
}
