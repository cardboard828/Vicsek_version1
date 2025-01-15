// get_orientations.h
#ifndef GET_ORIENTATIONS_H
#define GET_ORIENTATIONS_H

#include <vector>
#include <array>
#include "Particle.h"  // Particle構造体を含むヘッダーファイル

void get_orientations(
    std::vector<Particle>& particles, 
    std::vector<std::vector<std::vector<Particle>>>& cells,
    std::vector<std::array<double, 2>>& new_orientation_2,
    int L
);

#endif // GET_ORIENTATIONS_H
