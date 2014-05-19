# マクロ定義部
CC       = g++
OBJS     = main.o
CPPFLAGS = -Wall -std=c++11

# 生成規則部
approximateMatch: $(OBJS)
	$(CC) -o $@ $(OBJS)

.c.o:
	$(CC) $(CPPFLAGS) -c $<

main.o: main.cpp

clean:
	rm *.o *~ approximateMatch
