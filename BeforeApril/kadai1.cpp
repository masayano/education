// 課題1: 以下のURLからgenbank形式のファイルをダウンロードして、含まれる遺伝子数と遺伝子コード領域の総塩基数をカウントせよ。

#include <fstream>  // ファイル入出力 std::ifstream, std::ofstream を使えるようにする
#include <iostream> // 標準入出力     std::cout,     std::cin      を使えるようにする
#include <string>   // 文字列クラス   std::string                  を使えるようにする
#include <vector>   // 配列管理クラス std::vector                  を使えるようにする

// C++用の便利ライブラリ boost を使う
// ubuntu 12.04 では "sudo apt-get install libboost-all-dev" でインストールできる (2014.3.13現在)
#include <boost/tokenizer.hpp>    // 区切り文字認識を可能にする
#include <boost/lexical_cast.hpp> // 任意の型変換を可能にする

int main() {
    // 事前に
    // ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__DH10B_uid58979/NC_010473.gff
    // をダウンロードして
    // プログラムと同じ場所に置いておく
    std::ifstream ifs("NC_010473.gff");    // ファイルを開く
    std::string   buf;                     // 一行ずつ読み取っていく際に、行を格納するための文字列
    int geneNum    = 0; // 遺伝子数
    int codeLength = 0; // コード領域の長さ
    while(ifs && std::getline(ifs, buf)) { // ファイルを読めなくなるまで一行ずつ読む
        if(buf[0] != '#') {                                     // 行の先頭が # でない場合に以下の処理を行う
            std::vector<std::string> line;                      // 区切られた部分ごとに記憶するための文字列の配列を準備
            typedef boost::char_separator<char> char_separator; // 型が長いので略称をつける
            typedef boost::tokenizer<char_separator> tokenizer; // 同上
            char_separator sep("\t", "", boost::keep_empty_tokens); // タブ文字 "\t" を区切り文字にして、空のエントリ "" も残す
            tokenizer tokens(buf, sep);                             // 区切ったら
            for (tokenizer::iterator it = tokens.begin(); it != tokens.end(); ++it) { // 先頭から一個ずつ iterator でアクセスする
                line.push_back(*it); // 文字列の配列に新しく文字列を追加
            }
            const std::string& type = line[2]; // 場所の種別が記録されてるのは先頭から3つ目の文字列
            if(type == "gene") {
                ++geneNum;
            } else if(type == "CDS") {
                try {
                    const int start = boost::lexical_cast<int>(line[3]); // 開始位置
                    const int end   = boost::lexical_cast<int>(line[4]); // 終了位置
                    codeLength += (end - start + 1);                     // コードされている長さを加算
                } catch(boost::bad_lexical_cast) {                 // 数値に位置が変換できなかった場合
                    std::cout << "error at: " << buf << std::endl; // どんな行でエラーがあったか報告
                }
            }
        }
    }
    std::cout << "gene num   : " << geneNum    << std::endl; // 出力
    std::cout << "code length: " << codeLength << std::endl; // 同上
    return 0;
}
