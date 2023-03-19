// Copyright 2023 Grigorios Piperagkas. All rights reserved.
// Use of this source code is governed by a BSD-3-clause
// license that can be found in the LICENSE file.
 /*
//////////////////////////////////////////////////////////////////////////////

 Helper functions for Serial Multidirectional Search proposed by VJ Torczon(1989).
 
 Author:
 G.S. Piperagkas
 Date:
 18/03/2023
//////////////////////////////////////////////////////////////////////////////
*/

/*
//////////////////////////////////////////////////////////////////////////////
//
//CHECK BOUNDS FUNCTION
//
//////////////////////////////////////////////////////////////////////////////
*/
void check_bounds(problem prob, double * simplex, int i){
    //check for bounds and replace the simplex points into box
    int j,check;
    check=1;
        for (j=0;j<prob.dim;j++){
            if (simplex[i*prob.dim+j]>prob.Xmax[j]){
                simplex[i*prob.dim+j]=prob.Xmax[j];
                check=0;
            }
            else if (simplex[i*prob.dim+j]<prob.Xmin[j]){
                simplex[i*prob.dim+j]=prob.Xmin[j];
                check=0;
            }
        }
}

/*
//////////////////////////////////////////////////////////////////////////////
//
//RESTART SIMPLEX FUNCTION
//
//////////////////////////////////////////////////////////////////////////////
*/
void restart(problem prob, double *simplex,double *Fsimplex,double Fvec[],int fevs){
    int i,j;
    double vector[prob.dim];
    for (i=1;i<prob.dim;i++){
        for (j=0;j<prob.dim;j++){
            simplex[i*prob.dim+j]=prob.Xmin[j] + drand48()*(prob.Xmax[j]-prob.Xmin[j]);
        }
        check_bounds(prob,simplex,i);
        Fvec[i]=Objective(vector,prob);
        fevs++;
        Fsimplex[i]=Fvec[i];
    }
}





//////////////////////////////////////////////////////////////////////////////
// SORT 2 ARRAYS STRAIGHT INSERTION
//
// Sort 2 arrays in terms of the first arrays values...ascending order
// O(n^2)
// First array dim: N*M //// second array dim: N (sort N items)
/////////////////////////////////////////////////////////////////////////////

void pik2srt (int N,int M, double array[] ,double arr[])
{

// INPUT ARGUMENTS
// ---------------
// N, M   : Dimension of ARR.
// arr : Array for sorting.
// arr2: 2dim array for swaps

    int i, j, k;
    double a;
    double b[M];
    // SORTING
    // --------
    for (j=1; j<N; j++) // Pick out each element in turn
    {
        a = arr[j];
        for ( k=0; k<M;k++ ){
            b[k]= array[M*j+k];
        }
        i = j-1;
        while ((i>=0) && (arr[i]>a))    // Find where to insert it
        {
            arr[i+1] = arr[i];  // Move element i at position i+1
            for (k=0;k<M;k++){
                array[M*(i+1)+k]=array[M*i+k];
            }
            i--;
        }
        arr[i+1] = a;           // Insert it
        for (k=0;k<M;k++){
            array[M*(i+1)+k]=b[k];
        }
    }
}
