// 課題4: 添付のファイルはKEGG DBに登録されているEnzymeのリスト（をランダムにシャッフルしたもの）である。このリストを、EC番号順にソートせよ。

#include <algorithm> // std::sortなどの様々な関数を使用可能にする
#include <fstream>
#include <sstream>   // 文字をベタベタつなげるときは stringstream がよい
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
        tokenizer::iterator it1 = tok1.begin();
        const std::string ecString = *it1;
        ++it1;
        data = *it1;
        try {
            char_separator sep2(".", "", boost::keep_empty_tokens);
            const std::string numString = ecString.substr(3); // ecString の 4文字目以降を見る
            tokenizer tok2(numString, sep2);
            for (tokenizer::iterator it2 = tok2.begin(); it2 != tok2.end(); ++it2) {
                // std::cout << ecString << " > " << (*it2).size() << std::endl;
                // デバッグに使ったコードはできるだけ残しておく
                num.push_back(boost::lexical_cast<int>(*it2));
            }
        } catch(boost::bad_lexical_cast) {
            std::cout << "error at ec string: " << ecString << std::endl;
        }
    }
    std::string print() const {
        std::stringstream ss("");
        ss << "ec:";
        for(int i = 0; i < num.size(); ++i) {
            ss << num[i];
            if(i != (num.size() - 1)) {
                ss << ".";
            }
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
        const std::vector<int>& n = e.getNum();
        if     (num[0] <  n[0]) { return true; } 
        else if(num[0] == n[0]) {
            if     (num[1] <  n[1]) { return true; }
            else if(num[1] == n[1]) {
                if     (num[2] <  n[2]) { return true; }
                else if(num[2] == n[2]) {
                    if(num[3] < n[3]) { return true; } } } }
        return false;
    }
    bool operator>(const EnzymeData& e) const {
        const std::vector<int>& n = e.getNum();
        if     (num[0] >  n[0]) { return true; } 
        else if(num[0] == n[0]) {
            if     (num[1] >  n[1]) { return true; }
            else if(num[1] == n[1]) {
                if     (num[2] >  n[2]) { return true; }
                else if(num[2] == n[2]) {
                    if(num[3] > n[3]) { return true; } } } }
        return false;
    }
    bool operator==(const EnzymeData& e) const {
        const std::vector<int>& n = e.getNum();
        if((num[0]==n[0])&&(num[1]==n[1])&&(num[2]==n[2])&&(num[3]==n[3])) { return true; }
        return false;
    }
}; // クラス定義の最後はセミコロンをつける

int main() {
    std::ifstream ifs("enzyme_randomized.list");
    std::vector<EnzymeData> enzymeArray;
    std::string buf;
    while(ifs && std::getline(ifs, buf)) {
        EnzymeData e(buf);        // 酵素データクラスの実体(インスタンス)を文字列を素に作成
        enzymeArray.push_back(e); // 酵素データを配列に登録
    }
    std::sort(enzymeArray.begin(), enzymeArray.end()); // ソート
    std::ofstream ofs("enzyme_randomized.list.sort");  // 書き出し先のファイルを開く
    for(std::vector<EnzymeData>::iterator it = enzymeArray.begin(); it != enzymeArray.end(); ++it) {
        ofs << (*it).print() << std::endl; // 書き出し
    }
    return 0;
}
