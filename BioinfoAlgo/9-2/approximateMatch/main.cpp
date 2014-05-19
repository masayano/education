#include <string>
#include <vector>

std::vector<int> approximateMatch(
        const std::string& reference,
        const std::string& query,
        const int capability) {
    std::vector<int> output;
    const auto rLen   = reference.size();
    const auto qLen   = query    .size();
    const auto length = rLen - qLen + 1;
    for(auto i = 0U; i < length; ++i) {
        int temp = 0;
        for(auto j = 0U; j < qLen; ++j) {
            if(reference[i+j] == query[j]) {
            } else {
                ++temp;
            }
        }
        if(temp <= capability) {
            output.push_back(i);
        }
    }
    return output;
}

#include <iostream>

int main() {
    const std::string reference = "My name is Alice.";
    const std::string query     = "Akice";
    const int capability = 1;
    const auto output = approximateMatch(reference, query, capability);
    for(const auto idx : output) {
        std::cout << idx << std::endl;
    }
    return 0;
}
