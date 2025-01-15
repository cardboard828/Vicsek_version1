#include "Cell_Map.h"
#include <cmath>

void createParticleMap(const std::vector<Particle>& particles, int N, int L,
                       std::vector<std::vector<std::vector<Particle>>>& cells) {
    // 1x1セルのマップを作成
    cells.resize(L, std::vector<std::vector<Particle>>(L));

    // 各パーティクルに対して、セルの位置を計算し、そのセルに格納
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        int cell_x = static_cast<int>(std::floor(particles[i].x)) % L;
        int cell_y = static_cast<int>(std::floor(particles[i].y)) % L;
        cells[cell_x][cell_y].push_back(particles[i]);
    }
}

