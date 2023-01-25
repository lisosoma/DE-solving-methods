#include <stdio.h>
#include <math.h>

double a = 80., b = 24., c = 0., d = -5.; // PARAMS OF ODE


void GetLocalMatrix(double* MFE, int pow, double L){
    if(pow==1){
        MFE[0] = -a/L - b/2 + c*L/3;
        MFE[1] =  a/L + b/2 + c*L/6;
        MFE[2] =  a/L - b/2 + c*L/6;
        MFE[3] = -a/L + b/2 + c*L/3;
    }
    else{
        MFE[0] = -a*7/3/L  - b/2   + c*2*L/15;
        MFE[1] =  a*8/3/L  + b*2/3 + c*L/15;
        MFE[2] = -a/3/L    - b/6   - c*L/30;
        MFE[3] =  a*8/3/L  - b*2/3 + c*L/15;
        MFE[4] = -a*16/3/L - b*0   + c*8*L/15;
        MFE[5] =  a*8/3/L  + b*2/3 + c*L/15;
        MFE[6] = -a/3/L    + b/6   - c*L/30;
        MFE[7] =  a*8/3/L  - b*2/3 + c*L/15;
        MFE[8] = -a*7/3/L  + b/2   + c*2*L/15;
    }
}


void MakeEmpty(double* F, int N){
    for (int i=0; i<N; i++)
        F[i] = 0;
}


void GetGlobalMatrix(double* M, int N, int pow, double L){
    int i, j, k, ii, p;
    double MFE[(pow+1)*(pow+1)];
    GetLocalMatrix(MFE, pow, L);
    for (i=0; i<N; i+=pow){
        for (ii=i, k=0; k<pow+1; k++, ii++){
            for (j=i, p=0; p<pow+1; j++, p++){
                M[ii*(N+1) + j] += MFE[k*(pow+1) + p];
            }
        }
    }

    M[0] =              1; // Dirichlet boundary conditions
    M[1] =              0; // Dirichlet boundary conditions
    M[2] =              0; // Dirichlet boundary conditions
    M[(N+1)*(N+1)-1] =  1; // Dirichlet boundary conditions
    M[(N+1)*(N+1)-2] =  0; // Dirichlet boundary conditions
    M[(N+1)*(N+1)-3] =  0; // Dirichlet boundary conditions
}


void GetGlobalVector(double* F, int N, int pow, double L)
{
    if(pow==1){
        for (int i=1; i<N; i++)
            F[i] = -d*2*L/2;
    }
    else {
        for (int i=1; i<N; i+=2) {
            F[i] +=   -d*L/6;
            F[i+1] += -d*L*2/3;
            F[i+2] += -d*L/6;
        }
    }
    F[0] = 12.;  // Dirichlet boundary conditions
    F[N] = -13.; // Dirichlet boundary conditions
}


void Gauss(double* U, double* M, double* F, int N){
    double sum;
    for (int k=0; k<N; k++) {
        for (int i=k+1; i<N; i++) {
            sum = M[i*N + k] / M[k*N + k];
            for (int j=k; j<N; j++) {
                M[i*N + j] -= M[k*N + j] * sum;
            }
            F[i] -= F[k] * sum;
        }
    }

    U[N] = F[N] / M[(N)*(N)];
    for (int i=N-1; i>=0; i--){
        sum=0;
        for (int j=i+1; j<N; j++){
            sum += M[i*N + j] * U[j];
        }
        U[i] = (F[i] - sum) / M[i*N + i];
    }
    U[0] = 12;
}


double AnalyticSolution(double x){
    return -17.74383392 + 53.43758282 * exp(-3. * x / 10.) + 5. * x / 24.;
}


void to_csv(double* U, int N, double x){
    double j; int i;
    FILE *file = fopen("pow1_N40.csv", "w");
    for (i=0, j=2; i<N; i++, j+=x){
        fprintf(file, "%f,%f,%f\n", AnalyticSolution(j), U[i], j);
        printf("%f,%f,%f\n", AnalyticSolution(j), U[i], j); // DEBUG :)
    }
}

int main(int argc, char* argv[]){
    int pow = 1, N = 40; // PARAMS OF FEA
    int NN;               // Size of grid

    if (pow==1) NN = N;
    else NN = 2*N;

    double M[(NN+1)*(NN+1)]; // Global Matrix
    double F[NN+1];          // Global Vector
    double U[NN+1];          // Solution
    double l_all = 8;        // x in [2, 10]
    double L = l_all / N;

    MakeEmpty(M, (NN+1)*(NN+1));
    MakeEmpty(F, (NN+1));
    GetGlobalMatrix(M, NN, pow, L);
    GetGlobalVector(F, NN, pow, L);
    Gauss(U, M, F, NN+1);

    to_csv(U, NN+1, L/pow);
}
