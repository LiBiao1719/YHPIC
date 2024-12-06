
CXX = mpicxx
DEFINES = -w -Wunused-result -DDEBUG  -DNOSCALE -DHIGHORDERCLOUD -DPARALLEL -DPARALLEL_IO  -DMPICH_SKIP_MPICXX 
INCLUDES = -I/usr/local/hdf5-1.8.18/include 
LIBS = -L /usr/local/hdf5-1.8.18/lib -lhdf5 -lz -lm

CXXFLAGS = -O3 $(DEFINES)

SOURCES = \
fieldinfo.cpp laser.cpp main.cpp parparam.cpp readfile.cpp mytime.cpp\
cellinfo.cpp diagnostic.cpp ostack.cpp cell.cpp diagnostic_trace.cpp\

HEADERS = \
fieldinfo.h laser.h matrix.h parparam.h readfile.h otypen.h vector.h species.h mytime.h\
cell.h particle.h cellinfo.h diagnostic.h ostack.h diagnostic_trace.h\

OBJECTS = \
fieldinfo.o laser.o main.o parparam.o readfile.o mytime.o\
cellinfo.o diagnostic.o ostack.o cell.o diagnostic_trace.o\

oopic3d: $(SOURCES) $(HEADERS)
	$(CXX) -o lared-p3d $(CXXFLAGS) $(INCLUDES) $(SOURCES) $(LIBS)
	mv lared-p3d ../run/
.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< 
clean:
	rm *.o ../run/lared-p3d

