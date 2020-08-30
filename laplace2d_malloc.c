/*
 *  Copyright 2012 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/*
 * 
 *  The original code has been modified to use CUDA Unified memory and dynamic arrays.
 *
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "timer.h"

#define NN 4096
#define NM 4096

int main(int argc, char** argv)
{
    const int n = NN;
    const int m = NM;
    const int iter_max = 1000;
    double aux1;
    double **A,**Anew;
    int j;
    
    const double tol = 1.0e-6;
    double error     = 1.0;
    
    
    //Matrices A and Anew are allocated
    A = malloc(n * sizeof(double *));
    for(j = 0; j < n; j++)
        A[j] = malloc(m * sizeof(double));
    
    Anew = malloc(n * sizeof(double *));
    for(j = 0; j < n; j++)
        Anew[j] = malloc(m * sizeof(double));
    
    //Matrices A and Anew are initiallized
    for (int j = 0; j < n; j++)
    {
        A[j][0]    = 1.0;
        Anew[j][0] = 1.0;
        for (int k = 1; k< m; k++)
        {
            A[j][k]    = 0.0;
            Anew[j][k] = 0.0;               
        }
    }
    
    printf("A   [1][1] %lf \n", A[1][1]);
   
   
    #pragma acc enter data create(A[0:n][0:m],Anew[0:n][0:m])  //allocate memory in the device for both matrices
    #pragma acc update device(A[0:n][0:m])                 //initialize A in the device using the data from the host, because A has already been defined in the host.
    
    
    
    
    printf("Jacobi relaxation Calculation: %d x %d mesh\n", n, m);
    
    StartTimer();
    int iter = 0;
    
    while ( error > tol && iter < iter_max )
    {
        error = 0.0;
#pragma acc parallel loop present(A[0:n][0:m],Anew[0:n][0:m])  // "present" means that the data within brackets has already been allocated in the device. The reduction clause may not be neccessary...
        for( int j = 1; j < n-1; j++)
        {
              #pragma acc loop //this is optional. Como A es una matriz definida con un malloc, si no ponemos esto, no se paraleliza el bucle interior
            for( int i = 1; i < m-1; i++ )
            {
                Anew[j][i] = 0.25 * ( A[j][i+1] + A[j][i-1]
                                    + A[j-1][i] + A[j+1][i]);
                error = fmax( error, fabs(Anew[j][i] - A[j][i]));
            }
        }

#pragma acc parallel loop present(A[0:n][0:m],Anew[0:n][0:m]) private(aux1) //if we have auxiliary loop variables, we define them as private
        for( int j = 1; j < n-1; j++)
        {
              #pragma acc loop //this is optional. Como A es una matriz definida con un malloc, si no ponemos esto, no se paraleliza el bucle interior
            for( int i = 1; i < m-1; i++ )
            {
                A[j][i] = Anew[j][i];  
                aux1=(float)i+j;
                A[i][j]=A[i][j];//+aux1*0.000001;
            }
        }

        if(iter % 100 == 0) printf("%5d, %0.6f\n", iter, error);
        
        iter++;
    }
    
    
    printf("If we don't do <update self>, we have no value: A[1][1]=%lf \n", A[1][1]);
    
    #pragma acc update self(A[0:n][0:m])    //this allows to get the data from the device to the host in order to print it
    printf("When doing it: A[1][1]=%lf \n", A[1][1]);
    
    

    
    
    double runtime = GetTimer();
 
    printf(" total: %f s\n", runtime / 1000);
}
