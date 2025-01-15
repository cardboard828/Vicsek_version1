// getNeighboringCell.cpp
#include "getNeighboringCell.h"

std::vector<std::array<int, 2>> getNeighboringCell(int cell_x, int cell_y, int L) {
    std::vector<std::array<int, 2>> neighbors;
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            int neighbor_x = (cell_x + dx + L) % L;
            int neighbor_y = (cell_y + dy + L) % L;
            neighbors.push_back({neighbor_x, neighbor_y});
        }
    }
    return neighbors;
}
