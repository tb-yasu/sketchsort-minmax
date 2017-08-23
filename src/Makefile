CXX = g++

CXXFLAGS = -O3 -g -DNDEBUG -Wno-deprecated -ffast-math -pipe -I./ 
%CXXFLAGS = -O3 -g -DNDEBUG -Wno-deprecated -pedantic -ansi -finline-functions -foptimize-sibling-calls -Wcast-qual -Wwrite-strings -Wsign-promo -Wcast-align -Wno-long-long -fexpensive-optimizations -funroll-all-loops -ffast-math -fomit-frame-pointer -pipe -I./ -L./ 
%CXXFLAGS = -O3 -DNDEBUG -W -Wall -Wno-deprecated -pedantic -ansi -finline-functions -foptimize-sibling-calls -Wcast-qual -Wwrite-strings -Wsign-promo -Wcast-align -Wno-long-long -fexpensive-optimizations -funroll-all-loops -ffast-math -fomit-frame-pointer -pipe -I./ -L./
%CXXFLAGS = -O3 -DNDEBUG -W -Wall -Wno-deprecated
LINKFLAGS = -lm -lboost_system -L./

SRCS = \
	Main.cpp \
	SketchSort.cpp \

OBJS = $(SRCS:%.cpp=%.o)

all: sketchsort-minmax

sketchsort-minmax: $(OBJS)
	$(CXX) $(OTHERFLAGS) $(OBJS) $(LINKFLAGS) -o sketchsort-minmax
debug:
	make all CXXFLAGS="-ggdb -W -Wall -pedantic"

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) -c $<

clean:
	rm -f sketchsort-minmax *.o *~
