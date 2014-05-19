# マクロ定義部
CC       = g++
OBJS     = main.o KeywordTree.o KeywordNode.o SuffixTree.o SuffixNode.o
CPPFLAGS = -Wall -std=c++11

# 生成規則部
suffixTree: $(OBJS)
	$(CC) -o $@ $(OBJS)

.c.o:
	$(CC) $(CPPFLAGS) -c $<

main.o: main.cpp
KeywordTree.o: KeywordTree.cpp
KeywordNode.o: KeywordNode.cpp
SuffixTree.o: SuffixTree.cpp
SuffixNode.o: SuffixNode.cpp

clean:
	rm *.o *~ suffixTree
