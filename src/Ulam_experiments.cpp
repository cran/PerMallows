//
//  Ulam_experiments.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 16/06/14.
//  Copyright (c) 2014 Ekhine Irurozki. All rights reserved.
//

#include "Ulam_experiments.h"
#include "Generic.h"

void Ulam_experiments::gibbs_sampling_experiments(int m, double *theta, int model, int **samples, double *error, double *time){

    //read_permus_per_dist();//OJO WARNING - si es de clase
    fill_shapes_of_n();
    
    int     burning_period_samples = n_*n_;
    int     *sigma = new int[ n_ ];
    int     *sigma_prime = new int[ n_ ];
    Generic  gen;
    gen.generate_random_permutation( n_ , 1, sigma);
    int     dist_sigma = distance(sigma);
    //exp
    double acum_time = 0 ;
    gen.ini_chrono();
    int cont =0 , s_prev = 0;
    //    int s_prev;
    double  d_acumul_all_avg =  0,  d_acumul_last = 0, d_expected;
    //fin exp
    
    
    //for(int sample = 0 ; sample < m + burning_period_samples ; sample ++){ /// 2
     for ( int sample = 0 ; cont < 15 ; sample ++){////2
        int a;
        int b ;
        int dist_p = 0 ;
        do {
            a = rand() % ( n_ );
            b = rand() % ( n_ );
        }while (a == b);
        gen.insert_at(sigma, n_ , a , b , sigma_prime);
        bool make_swap = false;
        if( dist_sigma > ( dist_p =distance(sigma_prime) ) )  make_swap = true;
        else{
            double rand_double = (double)rand()/RAND_MAX;
            if(rand_double < exp(-theta[0])) make_swap = true;
        }
        if(make_swap){
            for(int i = 0  ; i < n_ ; i ++) sigma[ i ] = sigma_prime[ i ];
            dist_sigma = dist_p;
        }
        if(sample>=burning_period_samples){
            //samples[sample-burning_period_samples]=new int[ n_ ];
            //for(int i = 0  ; i < n_ ; i ++)   samples[ sample - burning_period_samples ][ i ] = sigma[ i ];
            d_acumul_last += dist_sigma;
        }
        
        //// 2
        int s = sample - burning_period_samples;
        if ((s+1) % 1000 == 0 ) {//s+1 permutations
            acum_time += gen.end_chrono();
            if ( acum_time >= 1 * 1000 ){
                if (cont == 0 ) time[cont] = acum_time;
                else time[cont] = time[cont-1]+ acum_time;
                acum_time = 0 ;
                
                d_acumul_all_avg = (double)((d_acumul_all_avg * ( s_prev ) ) + d_acumul_last )/(s+1+s_prev);
                //cout<<"s_g "<<s<<endl;
                s_prev += (s+1);
                sample -= s;
                d_expected = expectation(theta[0]);//
                d_acumul_last = 0;
                error[ cont ] = 0;
                error[cont ] = abs(d_expected-d_acumul_all_avg);
                cont++;
                error[ cont ] = (s+1+s_prev);//all but the last one are overwritte
            }
            gen.ini_chrono();
        }
        //// 1
        /*int s = sample - burning_period_samples;
        if (sample>=burning_period_samples && (s+1) % BUCKET_EXP_HAM == 0 ) {//s+1 permutations
            acum_time += gen.end_chrono();
            if (cont == 0 ) time[cont] = acum_time;
            else time[cont] = time[cont-1]+ acum_time;
            acum_time = 0 ;
            
            d_acumul_all_avg = (d_acumul_all_avg * (s + 1 - BUCKET_EXP_HAM)  + d_acumul_last ) /(s+1);
            
            d_expected = expectation(theta[0]);//
            d_acumul_last = 0;
            error[cont ] = abs(d_expected-d_acumul_all_avg);
            cont++;
            gen.ini_chrono();
        }*/
    }
    delete [] sigma_prime;
}

void Ulam_experiments::distances_sampling_experiments(int m, double theta, int **samples, double *error, double *time){
    double  distance_acum = 0, rand_distance = 0;
    double  *proba_acumul        = new double[ n_ ];
    int     target_distance;
    int * sigma_test = new int[n_];
    
    //read_permus_per_dist();//OJO WARNING - si es de clase
    fill_shapes_of_n();

    //exp
    double acum_time = 0 ;
    Generic gen;
    gen.ini_chrono();
    int cont =0 , s_prev = 0;
    //    int s_prev;
    double  d_acumul_all_avg =  0,  d_acumul_last = 0, d_expected;
    //fin exp

    //calculate_probas_at_each_distance( theta, proba, bound_l);
    proba_acumul[ 0 ] = 1; // exp(-theta*d) = exp (-theta *0)
    for (int i = 1 ; i < n_ ; i++)//acumulate the number of permus at each distance
        proba_acumul[i] = num_permus_per_dist_[i] * exp ( -theta * i ) + proba_acumul[ i - 1 ];
    //for (int s = 0 ; s < m ; s ++){
    for (int s = 0 ; cont < 15 ; s ++) {//halt (x-axis) time
        rand_distance = (double) rand() / (double)(RAND_MAX) * proba_acumul[ n_ - 1 ];
        target_distance = 0;
        while(proba_acumul[ target_distance ] <= rand_distance) target_distance++;
        //samples[ i ] = new int[ n_ ];
        //cout<<"ulam distance greene_niejenhuis_wilf "<<target_distance<<endl;
        generate_permu_with_given_LIS( n_ - target_distance, sigma_test);
        
        d_acumul_last += target_distance;
        distance_acum += (double)target_distance;
        ////2
        if ((s+1) % 100 == 0 ) {//s+1 permutations
            acum_time += gen.end_chrono();
            if ( acum_time >= 1 * 1000 ){
                if (cont == 0 ) time[cont] = acum_time;
                else time[cont] = time[cont-1]+ acum_time;
                acum_time = 0 ;
                d_acumul_all_avg = (d_acumul_all_avg * (s_prev)  + d_acumul_last ) /(s+1+s_prev);
                //h_acumul_all_avg[ i ] = (double)((h_acumul_all_avg[ i ]* (s_prev ) ) + h_acumul_last[ i ] )/(s+1+s_prev);
                
                //cout<<"s_d "<<s<<endl;
                s_prev += (s+1);
                s = 0 ;
                d_expected = expectation(theta);//
                d_acumul_last = 0;
                error[ cont ] = 0;
                error[cont ] = abs(d_expected-d_acumul_all_avg);
                cont++;
                error[ cont ] = (s+1+s_prev);//all but the last one are overwritten
            }
            gen.ini_chrono();
        }
        ////1
        /*if ((s+1) % BUCKET_EXP_HAM == 0 ) {//s+1 permutations
            acum_time += gen.end_chrono();
            if (cont == 0 ) time[cont] = acum_time;
            else time[cont] = time[cont-1]+ acum_time;
            acum_time = 0 ;
            
            d_acumul_all_avg = (d_acumul_all_avg * (s + 1 - BUCKET_EXP_HAM)  + d_acumul_last ) /(s+1);
            
            d_expected = expectation(theta);//
            d_acumul_last = 0;
            error[cont ] = abs(d_expected-d_acumul_all_avg);
            cont++;
            gen.ini_chrono();
        }*/
    }
    //cout<<"end"<<endl;
    //cout<<"average distance (dist_sampling) "<<distance_acum/m<<endl;
    delete [] proba_acumul;

}


void Ulam_experiments::learning_experiments_approx_1(int m, double *theta, int ** samples ){
    distances_sampling(m , theta[0], samples);
}

void Ulam_experiments::learning_experiments_approx_2(int**samples, int m, int model,  int*time,double * likeli){
    
    //Newton_raphson nr(n_);
    Generic gen;
    int     * sigma_0 = NULL;// not used new int[ n_ ];
    
    gen.ini_chrono();
    set_median(m , samples , sigma_0);
    *likeli = get_likelihood(m, samples, model, sigma_0);
    *time = gen.end_chrono();
    
    //delete [] sigma_0;
}

long double Ulam_experiments::get_likelihood(int m, int **samples, int model, int *sigma_0){
    Newton_raphson newton(n_);
    long double likelihood, psi = 0;
    int     dist_avg = 0 ;
    double  theta;
    
    //read_permus_per_dist();//OJO WARNING - si es de clase
    fill_shapes_of_n();
    
    dist_avg = min_dist_;// distance_to_sample(samples, m, sigma_0);
    theta    = newton.Newton_raphson_method( (double)dist_avg/m, -1.001,ULAM_DISTANCE, NULL, MALLOWS_MODEL, num_permus_per_dist_);
    for (int i = 0 ; i < n_ ; i ++ ) psi += num_permus_per_dist_[ i ] * exp (-theta * i );
    likelihood = - dist_avg * theta - m* log ( psi );
    return likelihood;
    
}

int Ulam_experiments::set_median(int m, int **samples, int *sigma_0){
    Generic gen;
    //gen.ini_chrono();
    if ( sum_dist_ == NULL ){
        num_perms_bucket_ = m;
        sum_dist_ = new long double [num_perms_bucket_ * BUCKET_EXP_HAM_LER ];
        for (int i = 0 ; i < num_perms_bucket_ * BUCKET_EXP_HAM_LER; i ++) sum_dist_[ i ] = 0;
    }
    for(int i = 0; i < m - num_perms_bucket_ ; i++)
		for(int j = m - num_perms_bucket_ ; j < m  ; j++){
			int d = distance(samples[ i ], samples[ j ]);
			sum_dist_[ i ] += d;
            sum_dist_[ j ] += d;
            calls_++;
		}
	for(int i = m - num_perms_bucket_ ; i < m ; i++)
		for(int j = i + 1 ; j < m ; j++){
			int d = distance(samples[ i ], samples[ j ]);
			sum_dist_[ i ] += d;
			sum_dist_[ j ] += d;
            calls_++;
		}
	
    min_dist_ = sum_dist_[ 0 ];
	int min_pos = 0 ;
	for(int i = 1; i < m ; i++) if ( min_dist_ > sum_dist_[ i ] ) {min_dist_ = sum_dist_[ i ]; min_pos = i;}
	//for(int i = 0; i < n_ ; i++) sigma_0[ i ] = samples[ min_pos][ i ];

    //double time = gen.end_chrono();
    //cout<<(m*m -m)/2<<" "<<calls_<<" "<<time<<endl;
    
    step_ ++;
    return min_dist_;
}


long double Ulam_experiments::get_likelihood_given_theta(int m, int **samples, int model, int *sigma_0, double *theta){
    Newton_raphson newton(n_);
    long double likelihood, psi = 0;
    int     dist_avg = 0 ;
    
    //read_permus_per_dist();//OJO WARNING - si es de clase
    fill_shapes_of_n();
    
    dist_avg = distance_to_sample(samples, m, sigma_0);
    for (int i = 0 ; i < n_ ; i ++ ) psi += num_permus_per_dist_[ i ] * exp (-theta[0] * i );
    likelihood = - dist_avg * theta[0] - m* log ( psi );
    return likelihood;
    
}



