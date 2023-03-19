// Copyright 2023 Grigorios Piperagkas. All rights reserved.
// Use of this source code is governed by a BSD-3-clause
// license that can be found in the LICENSE file.
 /*
//////////////////////////////////////////////////////////////////////////////

 Main program for Serial Multidirectional Search proposed by VJ Torczon(1989).
 
 Author:
 G.S. Piperagkas
 Date:
 18/03/2023
//////////////////////////////////////////////////////////////////////////////
*/

#include<iostream>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<ctime>
#include<iomanip>
#include<cstdio>

using namespace std;

// Algorithm Parameters
struct ALGO_par {
    double Pexp;
    double Pcontr;
};

// Problem parameters
struct problem {
    int dim;
    int bench;
    int nfpar;
    long maxfevs;
    long maxiters;
    double e;
    double* Xmin;
    double* Xmax;
    double* Fparams;
    int* Fiparams;
    
};


// Additional source files

#include "SOURCE_objective.cpp"
#include "SOURCE_helper.cpp"


/*
 //////////////////////////////////////////////////////////////////////////////
 
 MAIN PROGRAM
 
 //////////////////////////////////////////////////////////////////////////////
 */

int main(int argc, char * argv[])
{
    
    int i,j,iterations,fevs;
    double Finob;
    
    //Bounds for benchmarking functions
    double xmini[]={-100,-30,-5.12,-600,-20,-8,-10,-10,-4,-10,-10,-100,-65.536,-500};
    double xmaxi[]={100,30,5.12,600,30,8,10,10,5,10,10,100,65.536,500};
    
    
    //Define problem and algorithms structures
    ALGO_par algo;
    problem prob;
    
    prob.maxiters=1000;
    prob.maxfevs=1000000;
    algo.Pexp=2;
    algo.Pcontr=0.5;
    prob.nfpar=0;
    
    //-------------------------
    // DEFINE ALL PROBLEM DATA
    //-------------------------
    
    prob.dim = 100;
    prob.bench = 0; //select benchmarking function
    
    double Fvec[prob.dim+1];
    for (i=0;i<prob.dim;i++) Fvec[i]=0;
    
    
    //-------------------------
    // Define dynamic vectors
    //-------------------------
    
    double * simplex=NULL;
    simplex = new double[(prob.dim+1)*prob.dim];
    double * Fsimplex=NULL;
    Fsimplex = new double[prob.dim+1];
    double * refl=NULL;
    refl = new double[(prob.dim+1)*prob.dim];
    double * Frefl=NULL;
    Frefl = new double[prob.dim+1];
    double * expa=NULL;
    expa = new double[(prob.dim+1)*prob.dim];
    double * Fexpa=NULL;
    Fexpa = new double[prob.dim];
    double * contr=NULL;
    contr= new double[(prob.dim+1)*prob.dim];
    double * Fcontr=NULL;
    Fcontr = new double[prob.dim];
    
    double vector[prob.dim];
    
    prob.Xmin= new double[prob.dim];
    prob.Xmax= new double[prob.dim];
    
    for (i=0;i<prob.dim;i++){ //set bounds correctly from data
        prob.Xmin[i]=xmini[prob.bench];
        prob.Xmax[i]=xmaxi[prob.bench];
    }
    
    //-------------------------------------
    // Initialize parameters
    //-------------------------------------
    
    
    // Random number generator
    srand48(time(NULL));
    
    //============================================================
    //
    //Initialize simplex randomly changing one dimension at a time
    //
    //============================================================
    fevs=0;
    
    for (i=0;i<(prob.dim+1);i++){
        for (j=0;j<prob.dim;j++){
            simplex[i*prob.dim+j]=prob.Xmin[j] + drand48()*(prob.Xmax[j]-prob.Xmin[j]);
            vector[j]=simplex[i*prob.dim+j];
        }
        check_bounds(prob,simplex,i);
        Fvec[i]=Objective(vector,prob);
        fevs++;
        Fsimplex[i]=Fvec[i];
    }
    
    //initial sorting of values
    pik2srt(prob.dim+1,prob.dim,simplex,Fsimplex);
    for (i=0;i<prob.dim+1;i++){
        fprintf(stderr,"Sorted initial values of simplex: %f \n",Fsimplex[i]);
    }
    
    //=============================================
    //
    // LOOP UNTIL STOPPING CRITERION
    //
    //=============================================
    
    
    iterations=0;
    
    while (iterations<prob.maxiters){
        
        //if the simplex is stuck to local minima, restart while keeping best point
        if ((fabs(Fsimplex[0]-Fsimplex[prob.dim])<1)){
            restart(prob,simplex,Fsimplex,Fvec,fevs);
            cout << "RESTART!" << "\n";
        }
        
        //===================
        //
        //    Reflection
        //
        //===================
        i=0;
        while (i<prob.dim+1){
            for (j=0;j<prob.dim;j++){
                refl[i*prob.dim+j]=simplex[j]-(simplex[i*prob.dim+j]-simplex[j]);
            }
            check_bounds(prob,refl,i);
            for (j=0;j<prob.dim;j++){
                vector[j]=refl[i*prob.dim+j];
            }
            Fvec[i]=Objective(vector,prob);
            fevs++;
            Frefl[i]=Fvec[i];
            i++;
        }
        //sorting for reflection!
        pik2srt(prob.dim+1,prob.dim,refl,Frefl);
        
        //end of reflection
        
        //===========================
        //
        //    Expansion
        //
        //===========================
        
        if ((Frefl[0]<Fsimplex[0])){
            i=0;
            while ((i<prob.dim+1)){
                for (j=0;j<prob.dim;j++){
                    expa[i*prob.dim+j]=simplex[j]-algo.Pexp*(simplex[i*prob.dim+j]-simplex[j]);
                }
                check_bounds(prob,expa,i);
                for (j=0;j<prob.dim;j++){
                    vector[j]=expa[i*prob.dim+j];
                }
                Fvec[i]=Objective(vector,prob);
                fevs++;
                Fexpa[i]=Fvec[i];
                i++;
            }
            //sorting for expansion
            pik2srt(prob.dim+1,prob.dim,expa,Fexpa);
            //end of expansion
            
            if ((Fexpa[0]<Frefl[0])){
                for (i=0;i<prob.dim+1;i++){
                    for (j=0;j<prob.dim;j++){
                        simplex[i*prob.dim+j]=expa[i*prob.dim+j];
                    }
                    Fsimplex[i]=Fexpa[i];
                }
            }else{
                for (i=0;i<prob.dim+1;i++){
                    for (j=0;j<prob.dim;j++){
                        simplex[i*prob.dim+j]=refl[i*prob.dim+j];
                    }
                    Fsimplex[i]=Frefl[i];
                }
            }
        }else{
            //===================
            //
            //    Contraction
            //
            //===================
            for (i=0;i<prob.dim+1;i++){
                for (j=0;j<prob.dim;j++){
                    contr[i*prob.dim+j]=simplex[j]-algo.Pcontr*(simplex[i*prob.dim+j]-simplex[j]);
                }
                for (j=0;j<prob.dim;j++){
                    vector[j]=contr[i*prob.dim+j];
                }
                Fvec[i]=Objective(vector,prob);
                fevs++;
                Fcontr[i]=Fvec[i];
            }
            //sorting for contraction
            pik2srt(prob.dim+1,prob.dim,contr,Fcontr);
            //end of contraction
            for (i=0;i<prob.dim+1;i++){
                for (j=0;j<prob.dim;j++){
                    simplex[i*prob.dim+j]=contr[i*prob.dim+j];
                }
                Fsimplex[i]=Fcontr[i];
            }
        }
        iterations++;
        fprintf(stderr,"iterations= %d \n",iterations);
        fprintf(stderr,"best Obj value found: %14.7f \n", Fsimplex[0]);
        fprintf(stderr,"Worst Ob value found: %14.7f \n",Fsimplex[prob.dim]);
    }//end of while
    
    pik2srt(prob.dim+1,prob.dim,simplex,Fsimplex);
    for (i=0;i<prob.dim+1;i++){
        fprintf(stderr,"Final sorted value: %f \n",Fsimplex[i]);
    }
    
    fprintf(stderr,"Function evaluations: %d \n",fevs);
    for (i=0;i<prob.dim;i++) vector[i]=simplex[i];
    Finob=Objective(vector,prob);
    cout << "Final objective value: "<< Finob <<"\n";
    cout << "Decision vector: \n" ;
    for (i=0;i<prob.dim;i++){
        cout << vector[i] << "\n";
    }
}//end of main...

