# Compiler
CXX = g++-14
CXXFLAGS = -std=c++17 -O3 -fopenmp -Wall -Wextra -pedantic -Werror #-O0 -g#:デバッグ情報を追加, -O0:最適化を完全に無効化（変数や関数の情報が省略されなくなる）

# Target
TARGET = test
SRCS = vicsek_heiretsu_version0.cpp #main.cpp dist.cpp getNeighboringCell.cpp initialize_particles.cpp Cell_Map.cpp update_particles_and_order_parameter.cpp

# Build target
all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Clean up
clean:
	rm -f $(TARGET)
	rm -f *.o

# Run the program
run: $(TARGET)
	./$(TARGET)

