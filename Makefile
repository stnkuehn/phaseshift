# Dr. Steffen Kuehn / 2019

# target
CFILES = $(wildcard *.cpp)
OBJ = $(addprefix obj/,$(notdir $(CFILES:.cpp=.o)))
STYLEFIX_FILES = $(wildcard *.cpp) $(wildcard *.h)
	
# compiler
CPPLINUX = g++
BIN = phaseshift

# flags
LIBS = -lm -lrfftw -lfftw -L/opt/picoscope/lib -lps2000 -lps3000a -lsndfile -lgsl -lgslcblas
FLAGS = $(INCS) -Wall -std=c++11
INCS += -I/opt/picoscope/include/libps2000-2.1
INCS += -I/opt/picoscope/include/libps3000a-1.1

ifdef DEBUG
	FLAGS += -g2
else
	FLAGS += -O3
endif

# create all new
all: clean target 

$(BIN): $(OBJ)
	$(CPPLINUX) $(OBJ) -o $(BIN) $(LIBS)

obj/%.o: %.cpp
	@mkdir -p obj
	$(CPPLINUX) $(FLAGS) -c -o $@ $<
	
target: $(BIN)
	
stylefix:
	astyle --options=none --lineend=windows --style=allman --indent=force-tab=4 --break-blocks --indent-switches --pad-oper --pad-header --unpad-paren $(STYLEFIX_FILES)

# clean up all
clean:
	rm -rf obj
	rm -rf *.orig
	rm -f $(BIN)




