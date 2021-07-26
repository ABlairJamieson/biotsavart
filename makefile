

all: biotsavart.exe

biotsavart.exe: biotsavart.o
	g++ -std=c++11 -Wall -O2 biotsavart.o `root-config --libs` -o biotsavart.exe

biotsavart.o: biotsavart.C t2kstyle.h 
	g++ -std=c++11 -Wall -O2 -I`root-config --incdir` -c biotsavart.C

clean:
	rm *.o *.exe

.PHONY: all clean
