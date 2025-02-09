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

int main()
{
    std::chrono::system_clock::time_point  start, end; // 型は auto で可
    start = std::chrono::system_clock::now(); // 計測開始時間

    // System parameters
    const int L = 256;
    const double rho = 2.0;
    const int N = 131072;              //L * L * rho;
    const double v0 = 0.5;
    // const double eta = 0.4;
    // Radius of interaction is 1

    // Simulation parameters
    const int nsteps = 100000;
    const double dt = 1;

    const std::string output_folder = "ver0_output_variable=eta";
    
    //eta-loop
    for (int k=0; k<10; k++)
    {
        double eta=0.4757-0.0001*k;
        double sum=0;
        double order_parameter_ev=0;
        // ファイル名をetaに基づいて生成
        std::string filename = output_folder +"/ver0_output_eta=" + std::to_string(eta) + ".txt";
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

        // Initialize particles
        Particle particles[N];
        std::uniform_real_distribution<> init_pos(0, L);
        std::uniform_real_distribution<> dist_angle(0, 2 * M_PI);
        std::uniform_real_distribution<> white_noise(-eta * M_PI, eta * M_PI);

        #pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
            particles[i].x = init_pos(engine);
            particles[i].y = init_pos(engine);
            particles[i].orientation = dist_angle(engine);
        }
        
        for (int t = 0; t < nsteps; t++)
        {
            std::cout << "t = " << t << std::endl;

            int num_threads = omp_get_num_threads();
            //一旦32threadまで対応
            double order_parameter_x[32]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double order_parameter_y[32]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            double parameter_x=0;
            double parameter_y=0;
            double order_parameter=0;

            // Create a map of 1x1 cells to stored particles pointers
            std::vector<std::vector<std::vector<Particle>>> cells(L, std::vector<std::vector<Particle>>(L));
            for (int i = 0; i < N; i++)
            {
                int cell_x = particles[i].x;
                int cell_y = particles[i].y;
                cells[cell_x][cell_y].push_back(particles[i]);
            }

            // an vector with the shape shape to store the new orientation
            std::vector<std::array<double, 2>> new_orientation_2(N);
            #pragma omp parallel for
            for (int i = 0; i < N; i++) 
            {
                new_orientation_2[i][0] = 0; // x座標
                new_orientation_2[i][1] = 0; // y座標
            }

            // loop over all cells
            #pragma omp parallel for
            for (int i = 0; i < N; i++)
            {
                    int cell_x = floor(particles[i].x);
                    int cell_y = floor(particles[i].y);
                    // Get neighboring cells
                    std::vector<std::array<int,2>> neighbors = getNeighboringCell(cell_x, cell_y, L);
                        // Loop over all neighboring cells
                        for (size_t j = 0; j < neighbors.size(); j++)
                        {
                            int neighbor_x = neighbors[j][0];
                            int neighbor_y = neighbors[j][1];

                            // Loop over all particles in the neighboring cell
                            for (size_t k = 0; k < cells[neighbor_x][neighbor_y].size(); k++)
                            {
                                // Calculate distance between particles
                                double d = dist(particles[i], cells[neighbor_x][neighbor_y][k]);
                                if (d < 1)
                                {
                                    // Add angle to the sum
                                    new_orientation_2[i][0]+= cos(cells[neighbor_x][neighbor_y][k].orientation);
                                    new_orientation_2[i][1]+= sin(cells[neighbor_x][neighbor_y][k].orientation);
                                }
                            }
                        }
            }

            //並列処理に使用したスレッド数をゲット
            // Update particle orientation and position
                #pragma omp parallel
                {
                    int thread_id = omp_get_thread_num();
                    #pragma omp for 
                    for (int i = 0; i < N; i++)
                    {
                        particles[i].orientation = atan2(new_orientation_2[i][1], new_orientation_2[i][0]);
                        // Add white noise to the orientation
                        particles[i].orientation += white_noise(engine);
                        particles[i].x += v0 * cos(particles[i].orientation) * dt;
                        particles[i].y += v0 * sin(particles[i].orientation) * dt;
                        if(particles[i].x < 0)
                        {
                            particles[i].x += L;
                        }
                        else if(particles[i].x >= L)
                        {
                            particles[i].x -= L;
                        }
                        if(particles[i].y < 0)
                        {
                            particles[i].y += L;
                        }
                        else if(particles[i].y >= L)
                        {
                            particles[i].y -= L;
                        }
                        order_parameter_x[thread_id]+=cos(particles[i].orientation);
                        order_parameter_y[thread_id]+=sin(particles[i].orientation);
                    }
                }
                std::cout << "thread数: " << num_threads << std::endl;
                for (int i=0; i< num_threads; i++)
                {
                    parameter_x+=order_parameter_x[i];
                    parameter_y+=order_parameter_y[i];
                }
                order_parameter=sqrt(pow(parameter_x,2)+pow(parameter_y,2))/N;
                std::cout << "秩序変数: " << order_parameter << std::endl;
                sum+=order_parameter;
        }
        order_parameter_ev=sum/nsteps;
        end = std::chrono::system_clock::now();  // 計測終了時間
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count(); //処理に要した時間をミリ秒に変換

        file << "処理時間: " << elapsed << " ミリ秒" << '\n';
        file << "#order_parameter " << " " << "#eta"<< '\n';
        file << order_parameter_ev << " " << eta << '\n';     
        file.close();
    }
    return 0;
}
