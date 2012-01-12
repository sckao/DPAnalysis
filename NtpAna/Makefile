CC=g++ -g -fPIC -Wall -O -ansi -D_GNU_SOURCE -g -O2 -m64
ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags`
LIB=`${ROOTSYS}/bin/root-config --libs`
GLIB=`${ROOTSYS}/bin/root-config --glibs`

OBJECTS=MathFunctions.o AnaInput.o hJetTime.o Selection.o PhotonAna.o timeVsAmpliCorrector.o Sync.o

all:ana.exe

MathFunctions.o : MathFunctions.cc
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)
        
AnaInput.o : AnaInput.cc 
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

hJetTime.o : hJetTime.cc 
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

Selection.o : Selection.cc 
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

timeVsAmpliCorrector.o : timeVsAmpliCorrector.cc 
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

PhotonAna.o : PhotonAna.cc 
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

Sync.o : Sync.cc 
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

ana.exe : main.cc $(OBJECTS)
	$(CC) -o $@ $< $(OBJECTS) $(ROOTFLAG) $(LIB) $(GLIB)

clean : 
	rm -rf *.o ana.exe
