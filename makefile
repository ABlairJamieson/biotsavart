

all: hyperkbs2.exe

hyperkbs2.exe: hyperkbs2.o MCMC.o 
	g++ -std=c++11 -Wall -O2 hyperkbs2.o MCMC.o `root-config --libs` -o hyperkbs2.exe

hyperkbs2.o: hyperkbs2.C t2kstyle.h MCMCFcn.h
	g++ -std=c++11 -Wall -O2 -I`root-config --incdir` -c hyperkbs2.C

MCMC.o: MCMC.C MCMCFcn.h
	g++ -std=c++11 -Wall -O2 -I`root-config --incdir` -c MCMC.C

clean:
	rm *.o *.exe

.PHONY: all clean
