// SystemParams.h
#ifndef SYSTEM_PARAMS_H
#define SYSTEM_PARAMS_H

struct SystemParams {
    int L;
    double rho;
    int N;
    double v0;
    int nsteps;
    double dt;

    SystemParams(int L, double rho, int N, double v0, int nsteps, double dt)
        : L(L), rho(rho), N(N), v0(v0), nsteps(nsteps), dt(dt) {}
};

#endif // SYSTEM_PARAMS_H
