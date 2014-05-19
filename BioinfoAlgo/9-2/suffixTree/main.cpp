#include "KeywordTree.hpp"
#include "SuffixTree.hpp"

#include <iostream>
#include <string>

int main() {
    KeywordTree kt;
    kt.regString("Albert Einstein (/ˈælbərt ˈaɪnstaɪn/; German: [ˈalbɐt ˈaɪnʃtaɪn] ( listen); 14 March 1879 – 18 April 1955) was a German-born theoretical physicist.");
    kt.regString("ETH Zürich (German: Eidgenössische Technische Hochschule Zürich) is an engineering, science, technology, mathematics and management university in the city of Zürich, Switzerland.");
    kt.regString("Zürich or Zurich (German: Zürich [ˈt͡syːrɪç]; Swiss German: Züri [ˈt͡syɾi]) is the largest city in Switzerland and the capital of the canton of Zürich.");
    SuffixTree st(kt);
    const std::string query = "German:";
    std::cout << "[KeywordTree]" << std::endl << "\"" << query << "\" is in" << std::endl << kt.search(query);
    std::cout << "[SuffixTree]"  << std::endl << "\"" << query << "\" is in" << std::endl << st.search(query);
    return 0;
}
