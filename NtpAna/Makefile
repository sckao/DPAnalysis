CC = g++
ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs --glibs`
CXXFALGS = -g -fPIC -Wall -O -ansi -D_GNU_SOURCE -O2
CXX = -m64

OBJECTS = main.o MathFunctions.o AnaInput.o hJetTime.o PhotonAna.o timeVsAmpliCorrector.o

ana.exe : $(OBJECTS)
	$(CC) $(CXXFLAGS) $(CXX) $(ROOTFLAG)  $(OBJECTS) -lMinuit -o ana.exe

main.o : main.cc
	$(CC) -c $(CXX)  main.cc $(ROOTFLAG)
MathFunctions.o : MathFunctions.cc
	$(CC) -c $(CXX)  MathFunctions.cc $(ROOTFLAG)
AnaInput.o : AnaInput.cc
	$(CC) -c $(CXX)  AnaInput.cc $(ROOTFLAG)
hJetTime.o : hJetTime.cc
	$(CC) -c $(CXX) hJetTime.cc $(ROOTFLAG)
PhotonAna.o : PhotonAna.cc
	$(CC) -c $(CXX) PhotonAna.cc $(ROOTFLAG)
timeVsAmpliCorrector.o : timeVsAmpliCorrector.cc
	$(CC) -c $(CXX)  timeVsAmpliCorrector.cc $(ROOTFLAG)

clean : 
	rm -rf *.o ana.exe
