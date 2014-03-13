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
        tokenizer::iterator it = tok.begin();
        const std::string enzymeName   = *it;
        ++it;
        const std::string recognizeSeq = *it;
        // std::map への要素の追加は insert(std::make_pair(key, val)) で行う
        restrictionEnzymeMap.insert(std::make_pair(enzymeName, recognizeSeq));
    }
    return restrictionEnzymeMap;
}

std::vector<int> getCutIdxArray(
        const std::string& seq,
        const std::string& query) {
    std::vector<int> idxArray;
    int idx = seq.find(query);
    const int halfLength = query.size() / 2;
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
    for    (std::map<std::string, std::string>::const_iterator it = restrictionEnzymeMap.begin();
            it != restrictionEnzymeMap.end();
            ++it) {
        const std::string& enzymeName   = (*it).first;  // std::pair<key, val> を参照するときは
        const std::string& recognizeSeq = (*it).second; // key は first, val は second メンバを参照する
        restrictedIdxMap.insert(std::make_pair(enzymeName, getCutIdxArray(seq, recognizeSeq)));
    }
    return restrictedIdxMap;
}

std::vector<int> vectorCompaction(const std::vector<int>& in) {
    std::vector<int> tmp = in;
    std::sort(tmp.begin(), tmp.end());
    std::vector<int>::const_iterator it = tmp.begin();
    std::vector<int> out;
    out.push_back(*it);
    for(++it; it != tmp.end(); ++it) {
        const int val = (*it);
        if(out.back() != val) { out.push_back(val); }
    }
    return out;
}

std::vector<int> mergeArray(const std::map<std::string, std::vector<int> >& arrayMap) {
    std::map<std::string, std::vector<int> >::const_iterator it = arrayMap.begin();
    std::vector<int> tmp = (*it).second;
    //std::cout << (*it).first << std::endl;
    for(++it; it != arrayMap.end(); ++it) {
        const std::vector<int>& newArray = (*it).second;
        tmp.insert(tmp.end(), newArray.begin(), newArray.end());
        //std::cout << tmp.size() << std::endl;
    }
    //std::cout << tmp.size() << std::endl;
    std::vector<int> merged = vectorCompaction(tmp);
    return merged;
}

std::vector<int> getFragmentLengthArray(
        const int seqLength,
        const std::vector<int>& cutIdxArray) {
    std::vector<int> fragmentLengthArray;
    // TODO: どう書こうここ？
    return fragmentLengthArray;
}

std::map<std::string, std::vector<int> > getFragmentLengthArray(
        const int seqLength,
        const std::map<std::string, std::vector<int> >& cutIdxArrayMap) {
    std::map<std::string, std::vector<int> > fragmentLengthArrayMap;
    for    (std::map<std::string, std::vector<int> >::const_iterator it = cutIdxArrayMap.begin();
            it != cutIdxArrayMap.end();
            ++it) {
        const std::string      enzymeName          = (*it).first;
        std::cout << "Calculating fragment length of enzyme \"" << enzymeName << "\"...";
        const std::vector<int> fragmentLengthArray = getFragmentLengthArray(seqLength, (*it).second);
        std::cout << "finished." << std::endl;
        fragmentLengthArrayMap.insert(std::make_pair(enzymeName, fragmentLengthArray));
    }
    return fragmentLengthArrayMap;
}

void printFragmentLengthArray(const std::map<std::string, std::vector<int> >& fragmentLengthArrayMap) {
    for    (std::map<std::string, std::vector<int> >::const_iterator it = fragmentLengthArrayMap.begin();
            it != fragmentLengthArrayMap.end();
            ++it) {
        const std::string      enzymeName          = (*it).first;
        std::cout << "Printing fragment length of enzyme \"" << enzymeName << "\":" << std::endl;
        const std::vector<int> fragmentLengthArray = (*it).second;
        const int fragmentTypeNum = fragmentLengthArray.size();
        for(int i = 0; i < fragmentTypeNum; ++i) {
            std::cout << fragmentLengthArray[i] << std::endl;
        }
    }
}

void printFragmentLengthArray(const std::vector<int>& fragmentLengthArray_merged) {
    const int fragmentTypeNum = fragmentLengthArray_merged.size();
    for(int i = 0; i < fragmentTypeNum; ++i) {
        std::cout << fragmentLengthArray_merged[i] << std::endl;
    }
}

int main() {
    // 実行ファイルと同じ場所にあらかじめ NC_001416.fna を置いておくこと
    std::cout << "Reading fasta file...";
    const std::string seq = readFasta("NC_001416.fna");
    std::cout << "finished." << std::endl;

    // 実行ファイルと同じ場所にあらかじめ kadai3.dat を置いておくこと
    // kadai3 のフォーマットは "[酵素名]<タブ>[認識配列]" である
    std::cout << "Reading enzyme data...";
    const std::map<std::string, std::string> restrictionEnzymeMap
            = readRestrictionEnzymeData("kadai3.dat");
    std::cout << "finished." << std::endl;

    std::cout << "Calculating restrict index...";
    const std::map<std::string, std::vector<int> > restrictedIdxMap
            = restrictSequence(seq, restrictionEnzymeMap);
    std::cout << "finished." << std::endl;

    std::cout << "Merging restrict index...";
    const std::vector<int> cutIdxArray_merged = mergeArray(restrictedIdxMap);
    std::cout << "finished." << std::endl;

    const int seqLength = seq.size();
    const std::map<std::string, std::vector<int> > fragmentLengthArrayMap
            = getFragmentLengthArray(seqLength, restrictedIdxMap);

    std::cout << "Calculating fragment length of all enzyme...";
    const std::vector<int> fragmentLengthArray_merged
            = getFragmentLengthArray(seqLength, cutIdxArray_merged);
    std::cout << "finished." << std::endl;

    printFragmentLengthArray(fragmentLengthArrayMap);

    std::cout << "Printing fragment length of all enzyme: " << std::endl;
    printFragmentLengthArray(fragmentLengthArray_merged);

    return 0;
}
