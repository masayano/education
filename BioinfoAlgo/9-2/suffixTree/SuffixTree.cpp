#include "SuffixTree.hpp"

SuffixTree::SuffixTree(const KeywordTree& kt) : root(kt.getRoot()) {}

const std::string SuffixTree::search(const std::string& str) const {
    return root.search(root, "$"+str).toString();
}
