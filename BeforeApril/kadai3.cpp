// 課題3: λファージゲノムをEcoRI,HindIII,BamHI,NotIでそれぞれ切った時の断片の長さを列挙せよ。
//        また、全部で切ったときの断片の長さも列挙せよ。
//        ラムダファージのゲノム配列は以下のURLからダウンロードせよ。
//        ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Enterobacteria_phage_lambda_uid14204/NC_001416.fna

#include <algorithm> // std::sortなどの様々な関数を使用可能にする
#include <fstream>
#include <iostream>
#include <map>       // 連想配列 std::map を使用可能にする
#include <sstream>   // 様々な型をつなげていくことが出来る文字列ストリームを使用可能にする
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>

const char* indent = "    ";

std::string readFasta(const char* fileName) {
    std::ifstream ifs(fileName); 
    std::stringstream seqStream("");
    std::string buf;
    while(ifs && std::getline(ifs, buf)) {
        // ファイル中に一個しか配列が含まれていないと分かっているので
        // このような単純な処理で済ませている
        if(buf[0] != '>') { seqStream << buf; }
    }
    return seqStream.str();
}

std::map<std::string, std::string> readRestrictionEnzymeData(const char* fileName) {
    std::ifstream ifs(fileName);
    std::string buf;
    std::map<std::string, std::string> restrictionEnzymeMap;
    while(ifs && std::getline(ifs, buf)) {
        typedef boost::char_separator<char> char_separator;
        typedef boost::tokenizer<char_separator> tokenizer;
        char_separator sep("\t", "", boost::keep_empty_tokens);
        tokenizer tok(buf, sep);
        auto it = tok.begin();
        const auto& enzymeName   = *it;
        ++it;
        const auto& recognizeSeq = *it;
        // std::map への要素の追加は insert(std::make_pair(key, val)) で行う
        restrictionEnzymeMap.insert(std::make_pair(enzymeName, recognizeSeq));
    }
    return restrictionEnzymeMap;
}

void printRestrictEnzymeMap(const std::map<std::string, std::string>& restrictionEnzymeMap) {
    for(const auto& restrictPair : restrictionEnzymeMap) {
        std::cout << indent << restrictPair.first << ": " << restrictPair.second << std::endl;
    }
}

std::vector<int> getCutIdxArray(
        const std::string& seq,
        const std::string& query) {
    std::vector<int> idxArray;
    auto idx = seq.find(query);
    const auto halfLength = query.size() / 2;
    while(idx != std::string::npos) {
        // デバッグに使ったコードはできるだけ残しておく
        //std::cout << idx << std::endl;
        idxArray.push_back(idx + halfLength);
        idx = seq.find(query, idx+1);
    }
    return idxArray;
}

std::map<std::string, std::vector<int> > restrictSequence(
        const std::string& seq,
        const std::map<std::string, std::string>& restrictionEnzymeMap) {
    std::map<std::string, std::vector<int> > restrictedIdxMap;
    // std::map の全要素を洗うためにはイテレータを使う
    // const 付きのクラスのイテレータは const_iterator などでなければならない
    for(const auto& restrictPair : restrictionEnzymeMap) {
        const auto& enzymeName   = restrictPair.first;  // std::pair<key, val> を参照するときは
        const auto& recognizeSeq = restrictPair.second; // key は first, val は second メンバを参照する
        restrictedIdxMap.insert(std::make_pair(enzymeName, getCutIdxArray(seq, recognizeSeq)));
    }
    return restrictedIdxMap;
}

void printRestrictedIdxArray(const std::vector<int>& cutIdxArray) {
    for(const auto& cutIdx : cutIdxArray) {
        std::cout << cutIdx << " ";
    }
    std::cout << std::endl;
}

void printRestrictedIdxMap(const std::map<std::string, std::vector<int> >& restrictedIdxMap) {
    for(const auto& restrictedIdxPair : restrictedIdxMap) {
        std::cout << indent << restrictedIdxPair.first << ": ";
        const auto& cutIdxArray = restrictedIdxPair.second;
        printRestrictedIdxArray(cutIdxArray);
    }
}

std::vector<int> vectorCompaction(const std::vector<int>& in) {
    auto tmp = in;
    std::sort(begin(tmp), end(tmp));
    auto it = begin(tmp);
    std::vector<int> out;
    out.push_back(*it);
    for(++it; it != end(tmp); ++it) {
        const int val = (*it);
        if(out.back() != val) { out.push_back(val); }
    }
    return out;
}

std::vector<int> mergeArray(const std::map<std::string, std::vector<int> >& arrayMap) {
    auto it = arrayMap.begin();
    auto tmp = (*it).second;
    //std::cout << (*it).first << std::endl;
    for(++it; it != end(arrayMap); ++it) {
        const auto& newArray = (*it).second;
        tmp.insert(end(tmp), begin(newArray), end(newArray));
        //std::cout << tmp.size() << std::endl;
    }
    //std::cout << tmp.size() << std::endl;
    const auto merged = vectorCompaction(tmp);
    return merged;
}

void getFragmentLengthArray(
        const int seqLength,
        const std::vector<int>& cutIdxArray,
        const int index,
        const std::vector<bool>& flagArray,
        std::vector<int>& fragmentLengthArray) {
    const int length = cutIdxArray.size();
    if(index < length) {
        auto flagArray_false = flagArray;
        auto flagArray_true  = flagArray;
        flagArray_false.push_back(false);
        flagArray_true .push_back(true);
        getFragmentLengthArray(
                seqLength,
                cutIdxArray,
                index + 1,
                flagArray_false,
                fragmentLengthArray);
        getFragmentLengthArray(
                seqLength,
                cutIdxArray,
                index + 1,
                flagArray_true,
                fragmentLengthArray);
    } else {
        std::vector<int> tmpCutIdxArray;
        tmpCutIdxArray.push_back(0);
        for(int i = 0; i < length; ++i) {
            if(flagArray[i]) {
                //std::cout << cutIdxArray[i] << std::endl;
                tmpCutIdxArray.push_back(cutIdxArray[i]);
            }
        }
        const auto fragmentNum = tmpCutIdxArray.size();
        tmpCutIdxArray.push_back(seqLength);
        for(auto i = 1U; i < fragmentNum; ++i) {
            const auto fragmentLength = tmpCutIdxArray[i] - tmpCutIdxArray[i-1];
            fragmentLengthArray.push_back(fragmentLength);
        }
    }
}

std::vector<int> getFragmentLengthArray(
        const int seqLength,
        const std::vector<int>& cutIdxArray) {
    std::vector<int>  tmpFragmentLengthArray;
    std::vector<bool> flagArray;
    getFragmentLengthArray(
            seqLength,
            cutIdxArray,
            0,
            flagArray,
            tmpFragmentLengthArray);
    const auto length = tmpFragmentLengthArray.size();
    //std::cout << length << std::endl;
    if(length > 0) {
        const auto fragmentLengthArray = vectorCompaction(tmpFragmentLengthArray);
        return fragmentLengthArray;
    } else {
        return tmpFragmentLengthArray;
    }
}

std::map<std::string, std::vector<int> > getFragmentLengthArray(
        const int seqLength,
        const std::map<std::string, std::vector<int> >& cutIdxArrayMap) {
    std::map<std::string, std::vector<int> > fragmentLengthArrayMap;
    for(const auto& cutIdxArrayPair : cutIdxArrayMap) {
        const auto& enzymeName          = cutIdxArrayPair.first;
        std::cout << "Calculating fragment length of enzyme \"" << enzymeName << "\"...";
        const auto  fragmentLengthArray = getFragmentLengthArray(seqLength, cutIdxArrayPair.second);
        std::cout << "finished." << std::endl;
        fragmentLengthArrayMap.insert(std::make_pair(enzymeName, fragmentLengthArray));
    }
    return fragmentLengthArrayMap;
}

void printFragmentLengthArray(const std::map<std::string, std::vector<int> >& fragmentLengthArrayMap) {
    for(const auto& fragmentLengthArrayPair : fragmentLengthArrayMap) {
        const auto& enzymeName          = fragmentLengthArrayPair.first;
        std::cout << "Printing fragment length of enzyme \"" << enzymeName << "\":" << std::endl;
        const auto& fragmentLengthArray = fragmentLengthArrayPair.second;
        const auto fragmentTypeNum = fragmentLengthArray.size();
        const auto max             = fragmentTypeNum - 1;
        std::cout << indent;
        for(auto i = 0U; i < fragmentTypeNum; ++i) {
            std::cout << fragmentLengthArray[i];
            if(i != max) {
                std::cout << ", ";
            } else {
                std::cout << std::endl;
            }
        }
    }
}

void printFragmentLengthArray(const std::vector<int>& fragmentLengthArray_merged) {
    const auto fragmentTypeNum = fragmentLengthArray_merged.size();
    const auto max             = fragmentTypeNum - 1;
    for(auto i = 0U; i < fragmentTypeNum; ++i) {
        std::cout << fragmentLengthArray_merged[i];
        if(i != max) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
}

int main() {
    // 実行ファイルと同じ場所にあらかじめ NC_001416.fna を置いておくこと
    std::cout << "Reading fasta file...";
    const auto seq = readFasta("NC_001416.fna");
    std::cout << "finished." << std::endl;
    std::cout << "   Sequence length: " << seq.size() << std::endl;

    // 実行ファイルと同じ場所にあらかじめ kadai3.dat を置いておくこと
    // kadai3 のフォーマットは "[酵素名]<タブ>[認識配列]" である
    std::cout << "Reading enzyme data...";
    const auto restrictionEnzymeMap = readRestrictionEnzymeData("kadai3.dat");
    std::cout << "finished." << std::endl;
    printRestrictEnzymeMap(restrictionEnzymeMap);

    std::cout << "Calculating restrict index...";
    const auto restrictedIdxMap = restrictSequence(seq, restrictionEnzymeMap);
    std::cout << "finished." << std::endl;
    printRestrictedIdxMap(restrictedIdxMap);

    std::cout << "Merging restrict index...";
    const auto cutIdxArray_merged = mergeArray(restrictedIdxMap);
    std::cout << "finished." << std::endl;
    std::cout << indent << "Merged: ";
    printRestrictedIdxArray(cutIdxArray_merged);

    const auto seqLength = seq.size();
    const auto fragmentLengthArrayMap = getFragmentLengthArray(seqLength, restrictedIdxMap);

    std::cout << "Calculating fragment length of all enzyme...";
    const auto fragmentLengthArray_merged = getFragmentLengthArray(seqLength, cutIdxArray_merged);
    std::cout << "finished." << std::endl;

    printFragmentLengthArray(fragmentLengthArrayMap);

    std::cout << "Printing fragment length of all enzyme: " << std::endl;
    std::cout << indent;
    printFragmentLengthArray(fragmentLengthArray_merged);

    return 0;
}
