#ifndef CELL_MAP_H
#define CELL_MAP_H

#include <vector>
#include "Particle.h"

// 関数宣言
void createParticleMap(const std::vector<Particle>& particles, int N, int L,
                       std::vector<std::vector<std::vector<Particle>>>& cells);

#endif

