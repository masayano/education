// 課題2: 任意のアミノ酸配列に対応する塩基配列を列挙するプログラムを作成せよ。

#include <iostream>
#include <string>
#include <vector>

// ASCIIコードというルールでは
// char型8bit 256通り のうち
// 文字は正の数 128通り で表現される
const int charVariationNum = 128;

std::vector<std::vector<std::string> > makeAminoAcid2dna() {
    std::vector<std::vector<std::string> > matrix(charVariationNum); // 配列の1次元目の要素数を設定

    matrix[static_cast<int>('K')].push_back("AAA");
    matrix[static_cast<int>('N')].push_back("AAC");
    matrix[static_cast<int>('K')].push_back("AAG");
    matrix[static_cast<int>('N')].push_back("AAT");
    matrix[static_cast<int>('T')].push_back("ACA");
    matrix[static_cast<int>('T')].push_back("ACC");
    matrix[static_cast<int>('T')].push_back("ACG");
    matrix[static_cast<int>('T')].push_back("ACT");

    matrix[static_cast<int>('R')].push_back("AGA");
    matrix[static_cast<int>('S')].push_back("AGC");
    matrix[static_cast<int>('R')].push_back("AGG");
    matrix[static_cast<int>('S')].push_back("AGT");
    matrix[static_cast<int>('I')].push_back("ATA");
    matrix[static_cast<int>('I')].push_back("ATC");
    matrix[static_cast<int>('M')].push_back("ATG");
    matrix[static_cast<int>('I')].push_back("ATT");

    matrix[static_cast<int>('Q')].push_back("CAA");
    matrix[static_cast<int>('H')].push_back("CAC");
    matrix[static_cast<int>('Q')].push_back("CAG");
    matrix[static_cast<int>('H')].push_back("CAT");
    matrix[static_cast<int>('P')].push_back("CCA");
    matrix[static_cast<int>('P')].push_back("CCC");
    matrix[static_cast<int>('P')].push_back("CCG");
    matrix[static_cast<int>('P')].push_back("CCT");

    matrix[static_cast<int>('R')].push_back("CGA");
    matrix[static_cast<int>('R')].push_back("CGC");
    matrix[static_cast<int>('R')].push_back("CGG");
    matrix[static_cast<int>('R')].push_back("CGT");
    matrix[static_cast<int>('L')].push_back("CTA");
    matrix[static_cast<int>('L')].push_back("CTC");
    matrix[static_cast<int>('L')].push_back("CTG");
    matrix[static_cast<int>('L')].push_back("CTT");


    matrix[static_cast<int>('E')].push_back("GAA");
    matrix[static_cast<int>('D')].push_back("GAC");
    matrix[static_cast<int>('E')].push_back("GAG");
    matrix[static_cast<int>('D')].push_back("GAT");
    matrix[static_cast<int>('A')].push_back("GCA");
    matrix[static_cast<int>('A')].push_back("GCC");
    matrix[static_cast<int>('A')].push_back("GCG");
    matrix[static_cast<int>('A')].push_back("GCT");

    matrix[static_cast<int>('G')].push_back("GGA");
    matrix[static_cast<int>('G')].push_back("GGC");
    matrix[static_cast<int>('G')].push_back("GGG");
    matrix[static_cast<int>('G')].push_back("GGT");
    matrix[static_cast<int>('V')].push_back("GTA");
    matrix[static_cast<int>('V')].push_back("GTC");
    matrix[static_cast<int>('V')].push_back("GTG");
    matrix[static_cast<int>('V')].push_back("GTT");

    matrix[static_cast<int>('*')].push_back("TAA");
    matrix[static_cast<int>('Y')].push_back("TAC");
    matrix[static_cast<int>('*')].push_back("TAG");
    matrix[static_cast<int>('Y')].push_back("TAT");
    matrix[static_cast<int>('S')].push_back("TCA");
    matrix[static_cast<int>('S')].push_back("TCC");
    matrix[static_cast<int>('S')].push_back("TCG");
    matrix[static_cast<int>('S')].push_back("TCT");

    matrix[static_cast<int>('W')].push_back("TGA");
    matrix[static_cast<int>('*')].push_back("TGC");
    matrix[static_cast<int>('C')].push_back("TGG");
    matrix[static_cast<int>('C')].push_back("TGT");
    matrix[static_cast<int>('L')].push_back("TTA");
    matrix[static_cast<int>('F')].push_back("TTC");
    matrix[static_cast<int>('L')].push_back("TTG");
    matrix[static_cast<int>('F')].push_back("TTT");

    return matrix;
}

void printDNA(
        const std::vector<std::vector<std::string> >& aminoAcid2dna, // すでに存在するクラスを関数に渡すときは
        const std::string& aminoAcidArray,                           // ポインタか参照で渡すこと
        const std::size_t index,                                     // (時間とメモリの節約のため)
        const std::string& head) {
    if(index < aminoAcidArray.size()) {        // 受け取ったアミノ酸配列の終端まで来ていない場合に文字列を伸ばす処理を実行
        const auto aa = aminoAcidArray[index]; // 文字を受け取る
        const auto& list = aminoAcid2dna[static_cast<int>(aa)]; // 候補となるコドンのリストを取得
        for(const auto& codon : codonArray) {
            const auto newHead = head + " " + codon; // 見やすくするためコドン間に空白を入れておく
            // 再帰はたくさんしすぎてはいけないが
            // アミノ酸配列の長さはたかがしれているので
            // ここは簡単のためやる
            printDNA(
                    aminoAcid2dna,
                    aminoAcidArray,
                    index + 1,
                    newHead);
        }
    } else {
        // 終端まで来た場合は終了コドンをつけて出力
        std::cout << head << " TAA" << std::endl;
        std::cout << head << " TAG" << std::endl;
        std::cout << head << " TGC" << std::endl;
    }
}

int main() {
    // アミノ酸を添字に入力すると
    // 候補のコドンリストを 文字列の配列 として返してくれる
    // 要素が 文字列 である 2次元配列 を作成
    // (本当はもっと高速な方法があるが、理解しやすいコードにするため妥協)
    const auto aminoAcid2dna = makeAminoAcid2dna();

    std::string aminoAcidArray;
    std::cout << "Please input amino acid array: ";
    std::cin >> aminoAcidArray; // ユーザーからのキー入力を記録

    std::string head = "";
    printDNA(aminoAcid2dna, aminoAcidArray, 0U, head);

    return 0;
}
