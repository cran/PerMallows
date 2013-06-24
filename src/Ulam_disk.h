//
//  Ulam_disk.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 31/01/14.
//  Copyright (c) 2014 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Ulam_disk__
#define __perms_mallows__Ulam_disk__

#define GEN_NEXT  0 //ok, print and continue
#define GEN_TERM  1 //ok, terminate
#define GEN_EMPTY 2 //ok, print EMPTY SET and continue
#define GEN_ERROR 3 //an error occured, print an error message and terminate

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#include <vector>
#include <set>
#include <algorithm>
#include "Ulam.h"
#include "Newton_raphson.h"
#include "Ferrers_diagram.h"
#include "Exponential_model.h"


class Ulam_disk;
using namespace std;
class Ulam_disk: public Ulam{
protected:
/*    
    int     gen_part_init(unsigned char *vector, const unsigned char n, unsigned char *k);
    int     gen_part_next(unsigned char *vector, unsigned char *k, int bound);
 */
    void    read_permus_per_dist();
    void    read_permus_per_shape(int d, int bound, int * shape, int*shape_len);
    void    distances_sampling(int m, double theta, int**samples);
    void    generate_permu_with_given_LIS(int lis, int*sigma);
    void    estimate_theta(int m, int *sigma_0, int **samples, double *theta);
    long double get_likelihood(int m, int** samples, int model, int * sigma_0);
    

public:
    //int n_;
    char str_base_path [500] ;//= "/Users/eki/Desktop/aux_files/";
    Ulam_disk(int n) : Ulam(n){
        strcpy(str_base_path, "./");///Users/eki/Library/Developer/Xcode/DerivedData/perms_mallows-fxsohiuulqqerjgnrdpitiwnjkin/Build/Products/Debug/permus_per_shape_5_0
    }

    ~Ulam_disk(){};
    
    void    save_counts_to_file();
    double  expectation(double theta);
    /*************          ABSTRACT        ***************/
    //// local
    long double num_permus_at_distance(int d) {
        read_permus_per_dist();
        return num_permus_per_dist_[ d ];
    }
    /* same as the parent class Ulam        
     long double get_likelihood(int m, int** samples, int model, int * sigma_0) ;
     int     maximum_distance();
    int     distance(int*s1, int*s2);
    double  probability(int*s, int*s_0, double*theta) ;
    void    random_sample_at_dist(int dist, int m, int **samples) ;
    void    distances_sampling(int m, double theta, int**samples) ;
    void    estimate_consensus_approx (int model, int m, int **samples, int *sigma_0) ;
     int     perm2dist_decomp_vector(int*sigma, int*vec );
     void    dist_decomp_vector2perm(int* vec, int* sigma);
     
     void    multistage_sampling(int m, double*theta, int**samples);
     void    gibbs_sampling(int m, double *theta, int model, int **samples);
     
     double  estimate_consensus_exact  (int model, int m, int **samples, int*sigma_0_ini, int *sigma_0) ;*/
    
    /*************          ABSTRACT   NOT SUPPORTED     ***************/




};



#endif /* defined(__perms_mallows__Ulam_disk__) */


















