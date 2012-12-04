CC=g++ -g -fPIC -Wall -O -ansi -D_GNU_SOURCE -g -O2 -m64
ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags`
LDFLAGS = `${ROOTSYS}/bin/root-config  --ldflags --glibs` -lRooFit -lRooFitCore -lRooStats -lMinuit
LIB=`${ROOTSYS}/bin/root-config --libs`

OBJECTS=AnaInput.o hDraw.o StatsLib.o StatsTool.o 

all:ana.exe

AnaInput.o : AnaInput.cc AnaInput.h
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

hDraw.o : hDraw.cc hDraw.h
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

StatsLib.o : StatsLib.cc StatsLib.h
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

StatsTool.o : StatsTool.cc StatsTool.h
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB) $(LDFLAGS)

ana.exe : main.cc $(OBJECTS)
	$(CC) -o $@ $< $(OBJECTS) $(ROOTFLAG) $(LIB) $(LDFLAGS)

clean : 
	rm -rf *.o ana.exe
