#include "KeywordTree.hpp"

KeywordTree::KeywordTree() : root('$') {
    strCount = 0;
}

void KeywordTree::regString(const std::string& str) {
    const int length = str.size();
    for(int i = 0; i < length; ++i) {
        root.regString(str.substr(i), std::make_pair(i, strCount));
    }
    strCount += 1;
}

const std::string KeywordTree::search(const std::string& str) const {
    return root.search(root, str).toString();
}

    /* getter */

const KeywordNode& KeywordTree::getRoot() const {
    return root;
}
