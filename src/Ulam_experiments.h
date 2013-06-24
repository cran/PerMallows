//
//  Ulam_experiments.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 16/06/14.
//  Copyright (c) 2014 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Ulam_experiments__
#define __perms_mallows__Ulam_experiments__

#include <iostream>
#include "Ulam_disk.h"


class Ulam_experiments;
using namespace std;
class Ulam_experiments: public Ulam{
//class Ulam_experiments: public Ulam_disk{
public:
//    Ulam_experiments(int n):Ulam_disk(n){
    Ulam_experiments(int n):Ulam(n){
        step_ = 0 ;
        min_dist_ = 0;
        calls_ = 0;
        sum_dist_ = NULL;
    };
    ~Ulam_experiments(){
        if (sum_dist_ != NULL )
            delete [] sum_dist_;
    }
    void gibbs_sampling_experiments( int m, double *theta, int model, int **samples, double *error, double * time) ;
    void distances_sampling_experiments(int m, double theta, int**samples, double *error, double * time) ;
    
    int step_, min_dist_, num_perms_bucket_, calls_;
    long double *sum_dist_ ;
    int set_median(int m, int **samples, int *sigma_0);
    long double get_likelihood(int m, int **samples, int model, int *sigma_0);
    
    void learning_experiments_approx_1(int m, double *theta, int ** samples );
    void learning_experiments_approx_2(int**samples, int m, int model,  int*time,double * likeli);
    long double get_likelihood_given_theta(int m, int** samples, int model, int * sigma_0, double * theta);
};

#endif /* defined(__perms_mallows__Ulam_experiments__) */
