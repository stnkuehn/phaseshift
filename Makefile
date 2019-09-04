# Dr. Steffen Kuehn / 2019

# target
BIN = phaseshift
CFILES = $(wildcard *.cpp)
OBJ = $(addprefix obj/,$(notdir $(CFILES:.cpp=.o)))
STYLEFIX_FILES = $(wildcard *.cpp) $(wildcard *.h)
	
# compiler
CPP = $(TOOLCHAIN)-g++

# flags
LIBS = -lm -lrfftw -lfftw -L/opt/picoscope/lib -lps2000 -lsndfile -lgsl -lgslcblas
FLAGS = $(INCS) -Wall -std=c++11
INCS = -I/opt/picoscope/include/libps2000-2.1

ifdef DEBUG
	FLAGS += -g2
else
	FLAGS += -O3
endif

# create all new
all: clean target 

$(BIN): $(OBJ)
	$(CPP) $(OBJ) -o $(BIN) $(LIBS)

obj/%.o: %.cpp
	@mkdir -p obj
	$(CPP) $(FLAGS) -c -o $@ $<
	
target: $(BIN)
	
stylefix:
	astyle --options=none --lineend=windows --style=allman --indent=force-tab=4 --break-blocks --indent-switches --pad-oper --pad-header --unpad-paren $(STYLEFIX_FILES)

# clean up all
clean:
	rm -rf obj
	rm -rf *.orig
	rm -f $(BIN)




