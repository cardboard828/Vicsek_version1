#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <fstream>
#include <time.h>
#include <omp.h> 
#include <chrono>
#include "Particle.h"
#include "dist.h"
#include "getNeighboringCell.h"
#include "Cell_Map.h"
#include "initialize_particles.h"
#include "get_orientations.h"
#include "update_particles_and_order_parameter.h"


int main()
{
    // System parameters
    const int L = 256;
    const double rho = 2.0;
    const int N = 131072;              //L * L * rho;
    const double v0 = 0.5;
    // const double eta = 0.4;
    // Radius of interaction is 1
    // Simulation parameters
    const int nsteps = 10000;
    const double dt = 1;

    const std::string output_folder = "output";
    
    
    //eta-ruup
    for (int k=0; k<50; k++)
    {
        double eta=0.45+0.001*k;
        double sum=0;
        double order_parameter_ev=0;
        //whiteノイズの定義
        std::uniform_real_distribution<> white_noise(-eta * M_PI, eta * M_PI);
        // ファイル名をetaに基づいて生成
        std::string filename = output_folder + "/eta=0.45+0.001*" + std::to_string(k) + ".txt";
        // ファイルを開く
        std::ofstream file(filename);

        // Output the parameters
        file << "#L = " << L << '\n';
        file << "#rho = " << rho << '\n';
        file << "#N = " << N << '\n';
        file << "#v0 = " << v0 << '\n';
        file << "#nsteps = " << nsteps << '\n';
        file << "#dt = " << dt << '\n';
        
        // Random number generator
        std::random_device seed_gen;
        std::mt19937 engine(seed_gen());

        Particle particles[N];

        // Initialize particles
        initialize_particles(particles, N, L, engine);
        
        for (int t = 0; t < nsteps; t++)
        {
            std::cout << "t = " << t << std::endl;

            std::vector<double> order_parameter_x(N, 0);
            std::vector<double> order_parameter_y(N, 0);
            double order_parameter=0;

            // 1x1セルの粒子マップを作成する
            std::vector<Particle> particles(N); 
            std::vector<std::vector<std::vector<Particle>>> cells;

            createParticleMap(particles, N, L, cells);  // 関数を呼び出す

            // an vector with the shape shape to store the new orientation
            std::vector<std::array<double, 2>> new_orientation_2(N);
            #pragma omp parallel for
            for (int i = 0; i < N; i++) 
            {
                new_orientation_2[i][0] = 0; // x座標
                new_orientation_2[i][1] = 0; // y座標
            }

            get_orientations(particles, cells, new_orientation_2, L);

            update_particles_and_order_parameter(
                particles, 
                new_orientation_2, 
                order_parameter_x, 
                order_parameter_y, 
                white_noise, 
                engine, 
                N, 
                v0, 
                dt, 
                L,
                order_parameter
            );
            std::cout << "秩序変数: " << order_parameter << std::endl;
            sum+=order_parameter;
        }
        order_parameter_ev=sum/nsteps;

        file << "#order_parameter " << " " << "#eta"<< '\n';
        file << order_parameter_ev << " " << eta << '\n';     
        file.close();
    }
    return 0;
}
