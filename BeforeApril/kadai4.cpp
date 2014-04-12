// 課題4: 添付のファイルはKEGG DBに登録されているEnzymeのリスト（をランダムにシャッフルしたもの）である。
//        このリストを、EC番号順にソートせよ。

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

class EnzymeData { // 酵素データをまとめて扱うためにクラスを定義
    std::vector<int> num;  // ec番号
    std::string      data; // ec番号でない部分
public:
    EnzymeData(std::string str) { // クラス名と同じ名前の返り値がない関数をコンストラクタという。
                                  // クラスの実体(インスタンス)が生成されるときに呼ばれる。
        typedef boost::char_separator<char> char_separator;
        typedef boost::tokenizer<char_separator> tokenizer;
        char_separator sep1("\t", "", boost::keep_empty_tokens);
        tokenizer tok1(str, sep1);
        tokenizer::iterator it1 = begin(tok1);
        const auto ecString = *it1;
        ++it1;
        data = *it1;
        try {
            char_separator sep2(".", "", boost::keep_empty_tokens);
            const auto numString = ecString.substr(3); // ecString の 4文字目以降を見る
            tokenizer tok2(numString, sep2);
            for (const auto& token : tok2) {
                // std::cout << ecString << " > " << (*it2).size() << std::endl;
                num.push_back(boost::lexical_cast<int>(token));
            }
        } catch(boost::bad_lexical_cast) {
            std::cout << "error at ec string: " << ecString << std::endl;
        }
    }
    std::string print() const {
        std::stringstream ss("");
        ss << "ec:";
        const auto length = num.size();
        const auto last   = length - 1;
        for(auto i = 0U; i < length; ++i) {
            ss << num[i];
            if(i != last) { ss << "."; }
        }
        ss << "\t";
        ss << data;
        return ss.str();
    }
    // クラス内部のデータは関数を通じて操作する(推奨)
    const std::vector<int>& getNum () const { return num;  }
    const std::string&      getData() const { return data; }
    // "<", ">", "==" の演算子を定義することでクラス間に 大小 の概念を与えることができる
    // つまり ソート可能 になる
    bool operator<(const EnzymeData& e) const {
        const auto& n = e.getNum();
        const auto length = num.size();
        for(auto i = 0U; i < length; ++i) {
            const auto lValue = num[i];
            const auto rValue = n[i];
            if     (lValue < rValue) { return true; }
            else if(lValue > rValue) { break; }
        }
        return false;
    }
}; // クラス定義の最後はセミコロンをつける

int main() {
    // "ec:1.2.3.4<タブ>bra bra bra..."
    // みたいなフォーマットのファイルを用意して同じ場所に置いておく
    std::ifstream ifs("enzyme_randomized.list");
    std::vector<EnzymeData> enzymeArray;
    std::string buf;
    while(ifs && std::getline(ifs, buf)) {
        EnzymeData e(buf);        // 酵素データクラスの実体(インスタンス)を文字列を素に作成
        enzymeArray.push_back(e); // 酵素データを配列に登録
    }
    std::sort(begin(enzymeArray), end(enzymeArray)); // ソート
    std::ofstream ofs("enzyme_randomized.list.sort");  // 書き出し先のファイルを開く
    for(const auto& enzyme : enzymeArray) {
        ofs << enzyme.print() << std::endl; // 書き出し
    }
    return 0;
}
