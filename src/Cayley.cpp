//
//  Cayley.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#include "Cayley.h"
#include "Generic.h"
#include "Newton_raphson.h"
#include "Ulam.h"
#include <cmath>
#include <cfloat>

double  Cayley::expectation(double theta) {
        double sum=0;
        for(int j = 1 ; j < n_ ; j++){
            double ex = exp(theta);
            double denom = j+ex;
            sum += (double)j/denom;
            //sum += (double) j /(double)(exp(theta)+(double)j);
        }
        return (double)(sum );
    }
    
void    Cayley::expectation(double *theta, double *expect) {
  for (int i = 0 ; i < n_ - 1; i ++){
    int j = i + 1;
    expect[ i ] = (n_-j) / ((n_-j)+exp(theta[i]));
  }
  expect[ n_ - 1 ] = 0;
}

double Cayley::probability(int *s, int *s_0, double *theta){
    int *x = new int[ n_ ],*comp = new int[ n_ ], *inv = new int[ n_ ];
    for(int i= 0;i < n_; i ++ ) inv[s_0[ i ]-1] =i+1;
    double proba=1;
    for (int i = 0 ; i < n_; i ++ ) comp[ i ] =s[inv[ i ]-1]; 
    double *psi = new double[ n_ ];
    perm2dist_decomp_vector(comp, x);
    calculate_psi(theta, psi);
    for(int i = 0 ; i < n_ - 1 ; i ++ ) proba *= exp(-theta[ i ]*x[ i ])/psi[ i ];
    delete []psi;
    delete []comp;
    delete []inv;
    delete []x;
    return proba;
}

double Cayley::get_theta_log_likelihood(int m, int *x_acumul, int *x_acumul_variation, double *theta_estim){
    double likeli= 0;
    double *psi;
    psi = new  double[ n_ ];
    theta_estim[n_-1] = 0;
    for(int i = 0;i < n_ - 1; i ++ ){
        int x_i = x_acumul[ i ];
        if(x_acumul_variation != NULL)x_i += x_acumul_variation[ i ];
        int j = i+1;
        if(x_i == 0)x_i=1;
        if(x_i == m)  x_i = m-1;
        double xav = (double) x_i/m;
        if(xav != 0){
            theta_estim[ i ] =log(n_-j)-log(xav / (1-xav));
            //if(theta_estim[ i ] < 0)theta_estim[ i ] = 0;
            psi[ i ] =1+(n_-j)*exp(-theta_estim[ i ]);
            likeli += x_i * theta_estim[ i ] + m * log(psi[ i ]);
            if ( likeli != likeli)
                likeli= 0;
        }else theta_estim[ i ] = 0;
    }
    delete [] psi;
    return (-1 * likeli);
}
void Cayley::get_x_lower_bound_freq(int m, int ** samples_inv_freq, int ini_pos, int *min_bound_x){
    int *freq = new int[ n_ ]; for (int i = 0 ; i < n_; i ++ )freq[ i ] = 0;
    int maxFreq= 0;
    for(int j=ini_pos;j<n_-1;j++){
        for(int s= 0;s<n_;s++){
            freq[ s ] += samples_inv_freq[ j ][ s ];
            if(freq[ s ]> maxFreq)
                maxFreq=freq[ s ];
        }
        min_bound_x[ j ] = m - maxFreq;
        if(min_bound_x[ j ]<0)min_bound_x[ j ] = 0;
    }
    delete []freq;
}
void Cayley::get_x_lower_bound(int m, int ** sample, int ini_pos, int *x_min_bound){
    int *freq = new int[ n_ ]; for (int i = 0 ; i < n_; i ++ )freq[ i ] = 0;
    int max_freq= 0;
    for(int j=ini_pos;j<n_-1;j++){
        for(int s= 0;s<m;s++){
            freq[sample[ s ][ j ]-1]++;
            if(freq[sample[ s ][ j ]-1]> max_freq)
                max_freq=freq[sample[ s ][ j ]-1];
        }
        x_min_bound[ j ] =m-max_freq;
        if(x_min_bound[ j ]<0)x_min_bound[ j ] = 0;
    }
    //  for (int i = 0 ; i < n; i ++ ) cout<<xMinBound[ i ]<<" ";cout<<" xbound, pos: "<<iniPos<<endl;
    delete []freq;
}

/* returns 2 vectors. item cycle_items[ i ] belongs to cycle cycle_indices[ i ]*/
int Cayley::get_cycles(int *sigma, int *cycle_items, int *cycle_indices){
    bool*visited = new bool[ n_ ];
    for (int i = 0 ; i < n_; i ++ )visited[ i ] =false;
    int item_index, cont=0, cycle_index= 0;
    while(cont<n_) {
        item_index= 0;
        while(visited[item_index])item_index++;
        while (!visited[item_index]) {
            visited[item_index] = true;
            cycle_items[cont] =item_index+1;
            cycle_indices[cont] =cycle_index;
            item_index=sigma[item_index]-1;
            cont++;
        }
        cycle_index++;
    }
    delete [] visited;
    return cycle_index;
}

void Cayley::random_sample_at_dist(int d, int m, int **samples){
    //samples must be initialized before as samples = new int*[ m ];
    for (int i = 0 ; i < m; i ++ ){
        samples[ i ] =generate_permu_with_k_cycles(n_, n_ - d );
    }
}
double Cayley::calculate_psi(double *theta, double *psi_vector){
    double psi = 1;
    for (int i = 0 ; i <  n_ - 1; i ++ ){
        int j = i + 1;
        psi_vector[ i ] = 1 + ( n_ - j ) *(double) exp(-theta[ i ]);
        psi *= psi_vector[ i ];
    }
    return psi;
}

void Cayley::dist_decomp_vector2perm(int* vec, int* sigma) {
    //x has length n
    //forward generation of the permtation
    int *sigma_inv = new int [ n_ ];
    int item ;//from 0..n-1 !!!!
    for (int i = 0 ; i < n_ ; i++) sigma_inv [ i ] = -1;
    for (int pos = 0 ; pos < n_ ; pos ++) {
        if (vec[ pos ] == 0){// an item <= pos can close the cycle
            item = pos ;
            while ( sigma_inv [ item ] != -1 || ! item_closes_cycle(pos, item, sigma, sigma_inv)) item --;
        }else{
            item = - 1 ;//from 0..n-1 !!!!
            int random = rand() % (n_ - pos - 1);
            while ( random  >= 0 ){
                item++;
                if (sigma_inv [ item  ] == -1 &&  ! item_closes_cycle(pos , item, sigma, sigma_inv)  )     random --;
            }
        }
        sigma_inv[ item ] = pos + FIRST_ITEM;
        sigma[ pos ] = item + FIRST_ITEM;
    }
    delete [] sigma_inv;
}

bool Cayley::item_closes_cycle(int pos, int item, int *sigma, int *sigma_inv){
    while (item < pos ) item = sigma[ item ] - FIRST_ITEM;
    if ( item == pos ) return true;
    return false;
}
/*void Cayley::get_permutation_from_X_from_ini(int *x, int *sigma){
 int *sigma_inv = new int [ n_ ];
 int item ;
 for (int i = 0 ; i < n_ ; i++) sigma_inv [ i ] = -1;
 for (int pos = 0 ; pos < n_ ; pos ++) {
 if (x[ pos ] == 0){// an item <= pos can close the cycle
 item = pos + FIRST_ITEM;
 while ( sigma_inv [ item - FIRST_ITEM ] != -1 || ! item_closes_cycle(pos, item, sigma, sigma_inv)) item --;
 }else{
 item = FIRST_ITEM - 1 ;
 int random = rand() % (n_ - pos - 1);
 while ( random  >= 0 ){
 item++;
 if (sigma_inv [ item - FIRST_ITEM ] == -1 &&  ! item_closes_cycle(pos , item, sigma, sigma_inv)  )     random --;
 }
 }
 sigma_inv[ item - FIRST_ITEM ] = pos ;
 sigma[ pos ] = item ;
 }
 
 delete [] sigma_inv;
 }
 
 bool Cayley::item_closes_cycle(int pos, int item, int *sigma, int *sigma_inv){
 item -= FIRST_ITEM;
 while (item < pos ) item = sigma[ item ] - FIRST_ITEM;
 if ( item == pos ) return true;
 return false;
 }
*/

void Cayley::x_vector_to_permutation_backwards_old(int *x, int *permu){
    int distance= 0; for(int i = 0; i < n_ ; i++)if(x[ i ] == 1)distance++;
    // cyles[][] guarda un ciclo en cada fila. El prmer numero de cada fila
     // cycles[ i ][0] es un item que tiene x_j=0 y es el mayor del ciclo
     
     // para cada elem con x_j=1 se asigma aleatoriameente a cada alguno de los ciclos
     // con la restriccion de que pertenezca a un ciclo cuyo maximo item
     //
     //
    bool mute=true;
    int *cycle = new int[ n_ ];
    Generic gen;
    //for (int i = 0 ; i < n; i ++ ) cout<< x[ i ]<<" ";cout <<" x"<<endl;
    //v2
    //int len1=distance;//num posiciones con x[ i ] ==1
    int len0=n_-distance;//num posiciones con x[ i ] ==0
    
    int *cont = new int[len0];
    int **cycles = new int*[len0];
    for (int i = 0 ; i < len0; i ++ )cycles[ i ]  = new int[distance+1];
    for (int i = 0 ; i < len0; i ++ ){
        for(int j= 0;j<distance+1;j++)cycles[ i ][ j ] =-1;
        cont[ i ] = 0;
    }
    int iniCont= 0;
    //en ccycles van los indices de 0..n-1
    //al construir la permu se ponen lo indices de 1..n
    for (int i = 0 ; i < n_; i ++ )
        if(x[ i ] ==0){
            cycles[iniCont][cont[iniCont]++] =i;//se meten los items desde 0..n-1
            iniCont++;
        }
    
    //if(!mute) for (int i = 0 ; i < len0; i ++ )cout<<"cycles first "<<cycles[ i ][0]<<endl;
    
    double prob=1, partialProb= 0;
    for(int xind=n_-1;xind>= 0;xind--){
        if(x[xind] ==1){
            int cycleNum = get_most_prob_cycle(xind, cycles, len0, cont);
            prob *= partialProb;
            cycles[cycleNum][cont[cycleNum]++] =xind;//se meten los items desde 0..n-1
        }
    }
    // double permsInCyle=1;
    
    for (int i = 0 ; i < len0; i ++ ){
        int cycleLength = cont[ i ];
        int *rand = new int[cycleLength];
        gen.generate_random_permutation(cycleLength, 0, rand);//rand:0..n-1
        cycle=cycles[ i ];
        permu[cycle[rand[cycleLength-1]]] = cycle[rand[0]] + 1; //+1 <= permu[ i ] \in {1..n}
        for(int j= 0;j<cycleLength-1;j++)
            permu[cycle[rand[ j ]]] = cycle[rand[j+1]] + 1;////+1 <=permu[ i ] \in {1..n}
        delete []cycle;
    }
}
long double Cayley::num_permus_at_distance(int d){
    return stirling_matrix_[n_ ][ n_ - d];
}
long double Cayley::count_permus_with_cycles(int d){
    return stirling_matrix_[n_ ][ n_ - d];
}
void Cayley::print_stirling_matrix(){
    Generic gen;
    gen.print_long_double_matrix(stirling_matrix_, n_ + 1, n_ + 1);
}

int Cayley::distance(int * s, int * t){
    int *comp = new int[ n_ ], *sigma_inv = new int[ n_ ];
    Generic gen;
    for(int j = 0 ; j < n_ ; j++) sigma_inv[t[ j ]-1] =j + 1;
    for(int i = 0 ; i < n_ ; i++) comp[ i ] = s [ sigma_inv [ i ] - 1 ];
    int dist = perm2dist_decomp_vector(comp, NULL);
    delete [] sigma_inv;
    delete [] comp;
    return dist;
}

int Cayley::perm2dist_decomp_vector(int*sigma, int*vec ){
    //x is a vector of length n
    //also updates the x vector if it isnot null
    if(vec!=NULL)for (int i = 0 ; i < n_; i ++ )vec[ i ] =1;
    //for (int i = 0 ; i < n; i ++ )cout<<sigma[ i ]<<" ";cout<<endl;
    int num_cycles=0, num_visited=0, item= 0;
    bool*visited = new bool[ n_ ];
    for (int i = 0 ; i < n_; i ++ )visited[ i ] =false;
    while(num_visited < n_ ){
        item = num_cycles;
        while ( visited[ item ]) item++;
        num_cycles++;
        int max_item_in_cycle= 0;
        do{
            if ( item > max_item_in_cycle ) max_item_in_cycle = item;
            visited[ item ] =true;
            num_visited++;
            item = sigma[ item ]-1;
        }while ( !visited[ item ] );
        if(vec != NULL )vec[ max_item_in_cycle ] = 0;
    }
    delete [] visited;
    return (n_ - num_cycles );
}


void Cayley::distances_sampling(int m, double theta, int **samples) {
    int     target_dist = 0;
    double  rand_val;
    long double *acumul = new long double[ n_ ];//+1
    acumul[0] =exp(-theta  * 0) * stirling_matrix_[ n_ ][ n_ ];//stirling1(n_,n_);
    for(int dista=1;dista< n_ ;dista++)
        acumul[dista] =acumul[dista-1] +  exp(-theta  * dista) * stirling_matrix_[ n_ ][ n_ - dista ];//stirling1(n_,n_-dista);
    for (int i = 0 ; i < m; i ++ ){
        target_dist = 0;
        rand_val = (double) acumul[ n_ - 1 ] * (double) rand() / RAND_MAX;
        while(acumul[target_dist] <= rand_val) target_dist++;
        int *sigma=generate_permu_with_k_cycles(n_,(n_-target_dist));
        samples[ i ] =sigma;
    }
    delete []acumul;
}



void Cayley::multistage_sampling(int m, double *theta, int **samples){
    double *psi = new double[n_-1];
    int *x = new int[ n_ ];
    
    calculate_psi(theta, psi);
    for(int samp = 0 ; samp < m ; samp++){
        for(int i= 0;i < n_ - 1; i ++ ){
            double probaX_j = (double) 1/psi[ i ];//(double) exp(-theta[ i ])/psi[ i ]; //
            if(((double)rand() / RAND_MAX) < probaX_j) x[ i ] = 0;
            else x[ i ] =1;
        }
        x[n_-1] = 0;
        int *sigma = new int[ n_ ];
        dist_decomp_vector2perm(x, sigma);
        samples[samp] = sigma;
    }
    delete [] x;
    delete [] psi;
}


bool Cayley::same_cycle(int i, int j, int *sigma){
    int index = sigma[ i ] - 1;
    while(index != i && index != j)index = sigma[index]-1;
    if(index == j) return true;
    return false;
}

void Cayley::get_max_item_in_future_cycles(int *sigma, int i, int j, int *max_i, int *max_j){
    
    int pos=sigma[ i ]-1;
    *max_j = pos;

    while(pos != j ){
        pos  = sigma[pos]-1;
        if(pos > *max_j ) *max_j = pos;
    }
    pos=sigma[ j ]-1;
    *max_i = pos;
    while(pos != i){
        pos  = sigma[pos]-1;
        if(pos > *max_i ) *max_i = pos;
    }
}

void Cayley::get_max_item_in_current_cycle(int *sigma, int i, int *max_i ){
    int pos=sigma[ i ]-1;
    *max_i = pos;
    while(pos != i ){
        pos  = sigma[pos]-1;
        if(pos > *max_i ) *max_i = pos;
    }
}


void Cayley::gibbs_sampling(int m, double *theta, int model, int **samples) {
    int burning_period_samples = n_*log(n_);
    int *sigma = new int[ n_ ];
    Generic  gen;
    gen.generate_random_permutation(n_, 1, sigma);
    
    for(int sample= 0;sample<m+burning_period_samples;sample++){
        int i,j, max_i=-1, max_j=-1, min;
        do{
            i = rand() % n_;
            j = rand() % n_;
        }while(i == j);
        bool make_swap = false;
        if(  same_cycle(i, j, sigma) )  make_swap=true;
        else{
            double rand_double = (double)rand()/RAND_MAX;
            if(model == MALLOWS_MODEL){
                if(rand_double < exp(-theta[0])) make_swap = true;
            }else{
                get_max_item_in_current_cycle(sigma, i, &max_i);
                get_max_item_in_current_cycle(sigma, j, &max_j);
                //get_max_item_in_future_cycles(sigma, i, j, &max_i, &max_j);
                min = (max_i < max_j) ? max_i : max_j;
                if(rand_double < exp(-theta[min])) make_swap = true;
            }
        }
        if(make_swap){
            int aux=sigma[ i ];
            sigma[ i ] =sigma[ j ];
            sigma[ j ] =aux;
        }
        //gen.print_int_vector(sigma, n_);
        if(sample>=burning_period_samples){
            samples[sample-burning_period_samples]  = new int[ n_ ];
            for (int i = 0 ; i < n_; i ++ )   samples[sample-burning_period_samples][ i ] =sigma[ i ];
        }
    }
}


int Cayley::get_most_prob_cycle(int ind, int **cycles, int len, int *leng_cycles){
    //returns the cycle index
    Generic gen;
    int cont = 0;
    while(ind>cycles[cont][0]) cont++;
    //each cycles from cycle[cont] to cycle[len-1] has a particuler
    //probability of being selected as the terget
    //this prob is
    float*proba_cycles = new float[len];
    float*proba_cyclesNoAcum = new float[len];
    for (int i = 0 ; i < len; i ++ ) proba_cycles[ i ] = 0;
    float acumul= 0;
    for(int i=cont ; i < len; i ++ ){
        for(int j=cont;j<len;j++){
            //if(j==i)proba_cycles[ i ] *= (float)gen.factorial(lengCycles[ j ]);
            //else proba_cycles[ i ] *= (float) gen.factorial(lengCycles[ j ]-1);
            float multi;
            if(j==i) multi= (float)gen.factorial(leng_cycles[ j ]);
            else  multi= (float) gen.factorial(leng_cycles[ j ]-1);
            if(proba_cycles[ i ] ==0) proba_cycles[ i ] =multi;
            else proba_cycles[ i ] *= multi;
        }
        acumul += proba_cycles[ i ];
    }
    for (int i = 0 ; i < len; i ++ )proba_cyclesNoAcum[ i ] =(float)proba_cycles[ i ]/acumul;
    for(int i=cont ; i < len; i ++ )
        if(i!=0)proba_cycles[ i ] = (float)proba_cycles[ i ]/acumul + (float)proba_cycles[i-1];
        else proba_cycles[ i ] = (float)proba_cycles[ i ]/acumul ;
    float ran = (float)rand()/(RAND_MAX);
    delete []proba_cycles;
    for(int i=cont ; i < len; i ++ )
        if(proba_cycles[ i ]>ran){
            //*partialProb=proba_cyclesNoAcum[ i ];
            return i;
        }
    
    return len-1;
}
int * Cayley::generate_permu_with_k_cycles(int n, int k){//generateUARPermuWithKCyles
    int *sigma = new int[n];
    if(n==k) for (int i = 0 ; i < n; i ++ )sigma[ i ] =i+1;
    else
        if(k==1){
            Generic gen; int *cycle = new int[n];
            gen.generate_random_permutation(n, 0, cycle);
            for (int i = 0 ; i < n-1; i ++ ) sigma[cycle[ i ]] = cycle[i+1]+1;
            sigma[cycle[n-1]] =cycle[0]+1;
            delete []cycle;
        }else{
            long double  ran1 =((long double)rand() / ((double )RAND_MAX +1 ) );
            int ran2=-1;
            if(ran1 < (long double)(stirling_matrix_[ n - 1 ][ k - 1 ] / stirling_matrix_[ n ][ k ])){//stirling1(n-1, k-1) / stirling1(n,k))){////n is in a cycle by itself
                /* long double ran1 = rand() % (long double)stirling1(n,k),  ran2=-1; if(ran1 < stirling1(n-1, k-1)){ */
                int *small=generate_permu_with_k_cycles(n-1, k-1);
                for (int i = 0 ; i < n-1; i ++ )sigma[ i ] = small[ i ];
                sigma[n-1] =n;
                delete []small;
            }else{
                int *small=generate_permu_with_k_cycles(n-1, k);
                for (int i = 0 ; i < n-1; i ++ )sigma[ i ] =small[ i ];
                ran2 = rand() % (n-1) ;//0..n-2
                sigma[n-1] =sigma[ran2];
                sigma[ran2] =n;
                delete []small;
            }
        }
    return sigma;
}

long double Cayley::get_likelihood(int m, int** samples, int model, int * sigma_0){
    double  psi;
    int     dist    = 0 ;
    double  * theta = new double[ n_ ];
    long double loglikelihood;

    if(model == MALLOWS_MODEL){
        Newton_raphson newton(n_);
        double  *psi_vec = new double [n_];
        dist = distance_to_sample(samples, m, sigma_0);
        theta[ 0 ] = newton.Newton_raphson_method((double)dist/m, -10.001,CAYLEY_DISTANCE, MALLOWS_MODEL, NULL, NULL);
         for (int i = 1 ; i < n_ -1; i++) theta[ i ] = theta[ 0 ];
         psi = calculate_psi(theta, psi_vec);
         delete [] psi_vec;
         loglikelihood = - theta[0] * dist - m * log (psi);
    }else{
        int *x = new int[ n_ ], *x_acumul = new int[ n_ ], *inv = new int[ n_ ], *comp = new int[ n_ ];
        for(int i = 0 ; i < n_ ; i++) x_acumul[ i ] = 0;
        for(int i = 0 ; i < n_ ; i++) inv [ sigma_0[ i ] - 1 ] = i + 1;
        for (int i = 0 ; i < m; i ++ ){
            for(int j = 0 ; j < n_ ; j ++) comp[ j ] = samples[ i ] [ inv [ j ] - 1 ];
            perm2dist_decomp_vector(comp, x);
            for(int j = 0 ; j < n_ - 1 ; j++) x_acumul[ j ] += x[ j ];
        }
        loglikelihood = get_theta_log_likelihood(m, x_acumul, NULL, theta);//theta is an array of length n
        delete [] x;
        delete [] x_acumul;
        delete [] inv;
        delete [] comp;
    }
    delete [] theta;
    return loglikelihood;
}


 void Cayley::estimate_theta(int m, int *sigma_0, int **samples, int model, double *theta){
     if(model == MALLOWS_MODEL){
        int dist = distance_to_sample(samples, m, sigma_0);
        Newton_raphson newton(n_);
        theta[0] = newton.Newton_raphson_method((double)dist/m, -10.001,CAYLEY_DISTANCE, MALLOWS_MODEL, NULL, NULL);
    }else{
        int *x = new int[ n_ ], *x_acumul = new int[ n_ ], *inv = new int[ n_ ], *comp = new int[ n_ ];
        for(int i = 0 ; i < n_ ; i++) x_acumul[ i ] = 0;
        for(int i = 0 ; i < n_ ; i++) inv [ sigma_0[ i ] - 1 ] = i + 1;
        for (int i = 0 ; i < m; i ++ ){
            for(int j = 0 ; j < n_ ; j ++) comp[ j ] = samples[ i ] [ inv [ j ] - 1 ];
            perm2dist_decomp_vector(comp, x);
            for(int j = 0 ; j < n_ - 1 ; j++) x_acumul[ j ] += x[ j ];
        }
        get_theta_log_likelihood(m, x_acumul, NULL, theta);//theta is an array of length n
        delete [] x;
        delete [] x_acumul;
        delete [] inv;
        delete [] comp;
    }
}

int Cayley::distance_to_sample(int **samples, int m, int *sigma){
    int distance= 0;
    int *comp = new int[ n_ ], *sigma_inv = new int[ n_ ];
    for(int j = 0 ; j < n_ ; j ++) sigma_inv[sigma[ j ] - 1 ] = j + 1;
    for(int s = 0 ; s < m ; s ++){
        for(int i = 0 ; i < n_ ; i ++) comp[ i ] = samples[ s ][ sigma_inv [ i ] - 1 ];
        distance += perm2dist_decomp_vector(comp, NULL);
    }
    delete []sigma_inv;
    delete []comp;
    return distance ;
}

long double Cayley::count_permus_by_x(int *x){
    int d = 0, c_len = 0;
    for(int i = 0 ; i < n_ ; i++) if(x[ i ] == 1 ) d++;
    c_len = n_ - d;
    long double *c = new long double[c_len];
    for ( int i = 0 ; i < c_len ; i ++) c[ i ] = 0 ;
    int *k= new int[ n_ ], k_in=c_len;
    for( int i = n_ - 1 ; i >= 0 ; i --){// first top cycle item with which can be joind¡ed in the cycle
        if ( x[ i ] == 0 )k_in --;
        k[ i ] = k_in;
    }
    long double l = count_permus_by_x_core(0, x, c, c_len, k);
    delete [] c;
    delete [] k;
    return l;
}

long double Cayley::count_permus_by_x_core(int index, int *x, long double *c, int c_len, int *k){
    Generic gen;
    long double l = 0, p = 1;
    if( index == n_ ) {
        for(int i = 0 ; i < c_len ; i ++) p *= gen.factorial(c[ i ]);
        return p ;
    }
    if (x[index] == 0) return count_permus_by_x_core(index + 1, x, c , c_len, k);
    for(int i = k[index] ; i < c_len ; i ++){
        long double *c_i = new long double[c_len];
        for(int j = 0 ; j < c_len ; j ++) c_i[ j ] = c[ j ];
        c_i[ i ] ++;
        l += count_permus_by_x_core(index+1, x, c_i, c_len, k);
        delete [] c_i;
    }
    return l;
}

double Cayley::estimate_consensus_exact_gmm(int m, int **samples, int*sigma_0_ini, int *sigma_0){
    int     ** samples_inv = new int*[ m ];
    int     ** samples_inv_freq = new int*[ n_ ];
    int     * x_acum = new int[ n_ ];
    int     * sigma_0_aux = new int[ n_ ];
    int     * sigma_0_inv_aux = new int[ n_ ];
    double  likelihood = -DBL_MAX;

    
    Generic gen;
    for (int i = 0 ; i < m; i ++ ){
        samples_inv[ i ]  = new int[ n_ ];
        gen.invert(n_, samples[ i ], samples_inv[ i ]);
    }
    for (int i = 0 ; i < n_; i ++ ){ //ini
        samples_inv_freq[ i ]  = new int[ n_ ];
        for(int j= 0;j<n_;j++)samples_inv_freq[ i ][ j ] = 0;
        x_acum[ i ] = 0;
    }
    for(int i = 0 ; i < m; i ++ )for(int j= 0;j<n_;j++)samples_inv_freq[ j ][samples_inv[ i ][ j ]-1]++;
    for (int i = 0 ; i < n_; i ++ ) {sigma_0_aux[ i ] = -1 ; sigma_0_inv_aux[ i ] = -1 ;}
    
    //bounds
    estimate_consensus_approx(GENERALIZED_MALLOWS_MODEL , m, samples, sigma_0);
    likelihood = get_likelihood(m , samples , GENERALIZED_MALLOWS_MODEL , sigma_0);
    
    if ( sigma_0_ini != NULL ){
        //obtain the likelihod of the proposed sigma_0_ini.
        double like_ini = get_likelihood(m , samples , GENERALIZED_MALLOWS_MODEL , sigma_0_ini);
        //check both solutions: the proposed sigma_0_ini vs. the approx. Discard the worse
        if (like_ini > likelihood) {
            likelihood = like_ini;
            for (int i = 0 ; i < n_ ; i++) sigma_0[ i ] = sigma_0_ini[ i ];
        }
    }
    
    
    double visited_nodes = estimate_consensus_exact_gmm_core(m, 0 , samples, samples_inv, samples_inv_freq,
                                                   x_acum, sigma_0_aux, sigma_0_inv_aux, 0, sigma_0, &likelihood);
    //delete
    for (int i = 0 ; i < m; i ++ ) delete []samples_inv[ i ];
    delete [] samples_inv;
    for (int i = 0 ; i < n_; i ++ ) delete []samples_inv_freq[ i ];
    delete [] samples_inv_freq;
    delete [] x_acum;
    delete [] sigma_0_aux;
    delete [] sigma_0_inv_aux;
    return visited_nodes;
}

double Cayley::estimate_consensus_exact_gmm_core(int m, int pos, int ** samples, int ** samples_inv, int **samples_inv_freq,
                                       int * x_acum, int * current_sigma, int * current_sigma_inv, double current_likeli_bound,
                                       int * best_sigma, double * best_likeli){
    if(pos == n_ && current_likeli_bound >= (*best_likeli)){
        for (int i = 0 ; i < n_; i ++ )best_sigma[ i ] =current_sigma[ i ];
        (*best_likeli)=current_likeli_bound;
        return 1;
    }
    double visited_nodes= 0;
    bool trace=false;
    int *x_acum_var = new int[ n_ ];
    double *theta_estim = new double[ n_ ];
    for(int it= 0;it<n_;it++){//sigmaInv(pos)=it
        if(current_sigma[it] ==-1){
            int xIncr= 0;
            current_sigma_inv[pos] =it+1;
            current_sigma[it] =pos+1;
            int *pos_swaps = new int[ m ];
            for(int s= 0;s<m;s++){
                pos_swaps[ s ] =-1;
                if(samples[ s ][it] != current_sigma[it]){
                    int x= samples[ s ][it];
                    int y= samples_inv[ s ][pos]-1;
                    samples[ s ][it] =pos+1;
                    samples[ s ][y] =x;
                    samples_inv[ s ][pos] =it+1;
                    samples_inv[ s ][x-1] =y+1;
                    pos_swaps[ s ] =y;
                    xIncr++;
                    samples_inv_freq[pos][y]--;
                    samples_inv_freq[x-1][it]--;
                    samples_inv_freq[pos][it]++;
                    samples_inv_freq[x-1][y]++;
                }
            }
            
            int *xbound = new int[ n_ ];
            for (int i = 0 ; i < n_; i ++ )xbound[ i ] = 0;
            get_x_lower_bound_freq(m, samples_inv_freq, pos+1, xbound);
            
            for (int i = 0 ; i < n_; i ++ )x_acum_var[ i ] =x_acum[ i ];
            x_acum_var[pos] += xIncr;
            double likeliBound=get_theta_log_likelihood(m, x_acum_var, xbound, theta_estim);
            delete []xbound;
            
            if(likeliBound >= (*best_likeli) )
                visited_nodes += estimate_consensus_exact_gmm_core(m, pos+1, samples, samples_inv, samples_inv_freq, x_acum_var, current_sigma, current_sigma_inv, likeliBound, best_sigma, best_likeli);
            //else {cout<<"bounded at "<<pos<<" - ";for (int i = 0 ; i < n; i ++ )cout<<current_sigma[ i ]<<" ";cout<<endl;}
            current_sigma_inv[pos] =-1;
            current_sigma[it] =-1;
            for (int s = 0 ; s < m ; s ++){
                if(pos_swaps[ s ] != -1){
                    int y=pos_swaps[ s ];
                    int x=samples[ s ][y];
                    samples[ s ][it] =x;
                    samples[ s ][y] =pos+1;
                    samples_inv[ s ][pos] =y+1;
                    samples_inv[ s ][x-1] =it+1;
                    pos_swaps[ s ] =-1;
                    samples_inv_freq[pos][y]++;
                    samples_inv_freq[x-1][it]++;
                    samples_inv_freq[pos][it]--;
                    samples_inv_freq[x-1][y]--;
                }
            }
            delete []pos_swaps;
        }
    }
    delete [] theta_estim;
    delete [] x_acum_var;
    return visited_nodes+1;
}

double Cayley::estimate_consensus_exact_mm(int m, int **samples, int*sigma_0_ini, int *sigma_0){
    int **samples_inv = new int*[ m ];
    int *x_acum = new int[ n_ ];
    int *sigma_0_aux = new int[ n_ ];
    int *sigma_0_inv_aux = new int[ n_ ];
    Generic gen;
    for (int i = 0 ; i < m; i ++ ){
        samples_inv[ i ]  = new int[ n_ ];
        gen.invert(n_, samples[ i ], samples_inv[ i ]);
    }
    double best_distance = (n_ - 1) * m;//maximum
    for (int i = 0 ; i < n_; i ++ ) {
        sigma_0_aux[ i ] =-1;
        sigma_0_inv_aux[ i ] =-1;
        x_acum[ i ] = 0;
    }
    //bounds
    estimate_consensus_approx(MALLOWS_MODEL , m, samples, sigma_0);
    best_distance = distance_to_sample(samples , m , sigma_0);
    
    if ( sigma_0_ini != NULL ){
        //obtain the likelihod of the proposed sigma_0_ini.
        double dist_ini = distance_to_sample(samples , m , sigma_0_ini);
        //check both solutions: the proposed sigma_0_ini vs. the approx. Discard the worse
        if (dist_ini < best_distance) {
            best_distance = dist_ini;
            for (int i = 0 ; i < n_ ; i++) sigma_0[ i ] = sigma_0_ini[ i ];
        }
    }

    double visited_nodes = estimate_consensus_exact_mm_core(m, 0 , samples, samples_inv, x_acum, sigma_0_aux, sigma_0_inv_aux, 0, sigma_0, &best_distance);
    for (int i = 0 ; i < n_; i ++ ) delete []samples_inv[ i ];
    delete [] samples_inv;
    delete [] x_acum;
    delete [] sigma_0_aux;
    delete [] sigma_0_inv_aux;
    return visited_nodes;
}
double Cayley::estimate_consensus_exact_mm_core(int m, int pos, int ** samples, int ** samples_inv, int * x_acum,  int * current_sigma,
                              int * current_sigma_inv, double current_dist_bound, int *best_sigma, double *best_dist){
    if(pos == n_ && current_dist_bound <= (*best_dist)){
        for (int i = 0 ; i < n_; i ++ )best_sigma[ i ] =current_sigma[ i ];
        (*best_dist)=current_dist_bound;
        return 1;
    }
    double visited_nodes= 0;
    bool trace=false, enc=false;
    int *x_acum_var = new int[ n_ ], *candVec = new int[ n_ ];
    int cand=n_;
    int *freq = new int[ n_ ]; for (int i = 0 ; i < n_; i ++ )freq[ i ] = 0;
    for(int s = 0;(!enc && s<m );s++){
        if((freq[samples_inv[ s ][pos]-1]++) > (m/2)){
            candVec[0] =samples_inv[ s ][pos]-1;
            enc=true;
            cand=1;
        }
    }
    if(!enc){
        for(int it = 0 ; it < n_ ; it ++) candVec[ it ] = it;
        cand=n_;
    }
    for(int index = 0;index < cand ; index++){//sigmaInv(pos)=it
        int it = candVec[index];
        if(current_sigma[it] ==-1){
            int x_incr= 0;
            current_sigma_inv[pos] =it+1;
            current_sigma[it] =pos+1;
            int *pos_swaps = new int[ m ];
            //for (int i = 0 ; i < m; i ++ )posSwaps[ i ] =-1;
            for(int s= 0;s<m;s++){
                pos_swaps[ s ] =-1;
                if(samples[ s ][it] != current_sigma[it]){
                    int x = samples[ s ][it];
                    int y = samples_inv[ s ][pos]-1;
                    samples[ s ][it] = pos+1;
                    samples[ s ][y] = x;
                    samples_inv[ s ][pos] = it + 1;
                    samples_inv[ s ][x-1] = y + 1;
                    pos_swaps[ s ] =y;
                    x_incr++;
                }
            }
            
            double distance_bound = 0;
            int *xbound = new int[ n_ ];for (int i = 0 ; i < n_; i ++ )xbound[ i ] = 0;
            get_x_lower_bound(m, samples_inv, pos+1, xbound);
            for(int i=pos+1 ; i < n_; i ++ ) distance_bound+=xbound[ i ];
            delete []xbound;
            for (int i = 0 ; i < n_; i ++ ){
                x_acum_var[ i ] =x_acum[ i ];
                distance_bound += x_acum_var[ i ];
            }
            x_acum_var[pos] += x_incr;
            distance_bound += x_incr;
            
            if(distance_bound <= (*best_dist)) visited_nodes+=estimate_consensus_exact_mm_core(m, pos+1, samples, samples_inv, x_acum_var, current_sigma, current_sigma_inv, distance_bound, best_sigma, best_dist);
            //else {cout<<"bounded at "<<pos<<" - ";for (int i = 0 ; i < n; i ++ )cout<<current_sigma[ i ]<<" ";cout<<endl;}
            current_sigma_inv[pos] =-1;
            current_sigma[it] =-1;
            for(int s= 0;s<m;s++){
                if(pos_swaps[ s ] != -1){
                    int y = pos_swaps[ s ];
                    int x = samples[ s ][y];
                    samples[ s ][it] =x;
                    samples[ s ][y] =pos+1;
                    samples_inv[ s ][pos] =y+1;
                    samples_inv[ s ][x-1] =it+1;
                    pos_swaps[ s ] =-1;
                }
            }
            delete [] pos_swaps;
        }
    }
    
    delete [] x_acum_var;
    delete [] candVec;
    delete [] freq;
    return visited_nodes + 1 ;
}

void Cayley::estimate_consensus_approx(int model, int m, int **samples, int *sigma_0){
    int ** samples_inv = new int*[ m ];
    int ** samples_copy = new int*[ m ];
    for (int i = 0 ; i < m; i++) {
        samples_copy[ i ] = new int[ n_ ];
        samples_inv [ i ] = new int[ n_ ];
        for(int j = 0 ; j < n_ ; j ++) {samples_inv[ i ][ samples[ i ][ j ] - 1 ] = j + 1;samples_copy[ i ][ j ] = samples[ i ][ j ];}
    }
    double bl_data;
    double *best_likeli = &bl_data;
    if( model == MALLOWS_MODEL )
        estimate_consensus_approx_mm(m, samples_copy, samples_inv, sigma_0, best_likeli);
    else
        estimate_consensus_approx_gmm(m, samples_copy, samples_inv, sigma_0, best_likeli);
    variable_neighborhood_search(m, samples, sigma_0, model , best_likeli);
    for (int i = 0 ; i < m; i++) {
        delete [] samples_inv[ i ];
        delete [] samples_copy[ i ];
    }
    delete [] samples_copy;
    delete [] samples_inv;
}

void Cayley::estimate_consensus_approx_gmm(int m, int **samples_copy, int **samples_inv, int *sigma_0, double *best_likeli){
    // pre: sigma, sigmaInv =-1
    int *freq = new int[ n_ ];     
    int *x_acum = new int[ n_ ];   for (int i = 0; i < n_; i++) x_acum[ i ] = 0;
    int * sigma_0_inv = new int[ n_ ];
    
    for (int i = 0; i < n_; i++){ sigma_0[ i ] = -1;sigma_0_inv[ i ] = -1;}

    for (int item = 0 ; item < n_ ; item ++){
        //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n; i ++ )cout<<samples[ s ][ i ]<<" ";cout<<"samples  trace  1"<<endl;}
        for (int i = 0; i < n_; i++) freq[ i ] = 0;
        int max_freq = 0;
        int pj = -1; // 0..n-1
        //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n; i ++ )cout<<samples[ s ][ i ]<<" ";cout<<"samples  trace  2"<<endl;}
        for (int s = 0;  s < m; s++){
            freq[ samples_inv[ s ][item] - 1] ++;
            if ((freq[ samples_inv[ s ][item] - 1]) > max_freq){
                max_freq = freq[ samples_inv[ s ][item] - 1];
                pj = samples_inv[ s ][item] - 1;
            }
        }
        //for (int i = 0 ; i < n; i ++ ){for(int j= 0;j<n;j++)cout<<freq[ i ][ j ]<<" ";cout<<endl;}
        sigma_0_inv[item] = pj + 1;
        sigma_0 [pj] = item + 1;
        for (int s = 0; s < m; s++) {
            if (samples_copy[ s ][pj] != item + 1) {//swap
                int x = samples_copy[ s ][pj];
                int y = samples_inv[ s ][item] - 1;
                samples_copy[ s ][pj] = item + 1;
                samples_copy[ s ][y] = x;
                samples_inv [ s ][item] = pj + 1;
                samples_inv [ s ][x - 1] = y + 1;
                x_acum[ item ]++;
            }
        }
    }
    
    double *theta = new double[ n_ ];
    *best_likeli = get_theta_log_likelihood(m, x_acum, NULL, theta);
    delete [] theta;
    delete [] x_acum;
    delete [] sigma_0_inv;
    delete [] freq;
}

void Cayley::estimate_consensus_approx_mm(int m, int **samples_copy, int **samples_inv, int *sigma_0, double *best_distance){
    // pre: sigma, sigmaInv =-1
    int distance_increase = 0, remaining = n_;
    int **freq = new int*[ n_ ];for (int i = 0; i < n_; i++) freq[ i ] = new int[ n_ ];
    
    for (int i = 0; i < n_; i++) sigma_0[ i ] = -1;
    do{ //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n; i ++ )cout<<samples[ s ][ i ]<<" ";cout<<"samples  trace  1"<<endl;}
        for (int i = 0; i < n_; i++)
            for (int j = 0; j < n_; j++)freq[ i ][ j ] = 0;
        int max_freq = 0;
        int pi = -1, pj = -1; // 0..n-1
        bool dirty_items = true ; //si hay alguno q aparezca >m/2 veces?
        
        //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n; i ++ )cout<<samples[ s ][ i ]<<" ";cout<<"samples  trace  2"<<endl;}
        Generic gen;
        //gen.print_int_matrix(samples_copy, m, n_);
        for (int s = 0; dirty_items && s < m; s++)
            for (int i = 0; ( dirty_items && i < n_ ) ; i++){
                //gen.print_int_matrix(freq, n_, n_);
                if (sigma_0[ i ] == -1) (freq[ i ][samples_copy[ s ][ i ] - 1])++;
                if ((freq[ i ][samples_copy[ s ][ i ] - 1]) > max_freq){
                    max_freq = (freq[ i ][samples_copy[ s ][ i ] - 1]);
                    pi = i;
                    pj = samples_copy[ s ][ i ] - 1;
                    if (max_freq > m/2 ) dirty_items = false ;
                }
            }
        //for (int i = 0 ; i < n; i ++ ){for (int j = 0 ; j < n;j++)cout<<freq[ i ][ j ]<<" ";cout<<endl;}
        sigma_0[pi] = pj + 1;
        for (int s = 0; s < m; s++) {
            if (samples_copy[ s ][pi] != pj + 1) {//////////////swap
                //swapSamplesDo(s, pi, pj, samples, samples_inv, posSwap);
                int x = samples_copy[ s ][pi];
                int y = samples_inv[ s ][pj] - 1;
                samples_copy[ s ][pi] = pj + 1;
                samples_copy[ s ][y] = x;
                samples_inv[ s ][pj] = pi + 1;
                samples_inv[ s ][x - 1] = y + 1;
                distance_increase++;
            }
        }
        remaining--;
    }while(remaining > 0) ;
    (*best_distance) = distance_increase;
    for (int i = 0 ; i < n_ ; i ++) delete [] freq[ i ];
    delete [] freq;
}

void Cayley::variable_neighborhood_search(int m, int **samples, int *sigma, int model, double *f_eval){
    bool improve;
    do{
        double f_eval_ini = (*f_eval);
        //cout<<" distance ini"<<(*f_eval)<<endl;
        if( model == MALLOWS_MODEL ) local_search_swap_mm  (m, samples, sigma, f_eval);
        else local_search_swap_gmm  (m, samples, sigma, f_eval);
        //cout<<" distance new swap "<<(*f_eval)<<endl;
        local_search_insert(m, samples, sigma, model, f_eval);
        //cout<<" distance new ins "<<(*f_eval)<<endl;
        improve=false ;
        if( model == MALLOWS_MODEL  && (f_eval_ini) > *f_eval ) improve = true;
        if( model == GENERALIZED_MALLOWS_MODEL  && (f_eval_ini) < *f_eval ) improve = true;
    }while(improve);
}

void Cayley::local_search_swap_mm(int m, int **samples, int *sigma_0, double *f_eval){
    int ** samples_comp = new int*[ m ];
    for (int s = 0 ; s < m; s++)  samples_comp[ s ] = new int[ n_ ];
    int *sigma_0_inv = new int[ n_ ];
    int *cycle_items = new int [n_ ], *cycle_index = new int[ n_];
    int **same_cycle = new int*[ n_ ];
    for (int i = 0 ; i < n_ ; i ++){ same_cycle[ i ] = new int[ n_ ]; for (int j = 0 ; j < n_ ; j++) same_cycle[ i][ j ] = 0;}
    Generic gen;
    int index_i,index_j, distance_variation = 0;
    for(int i = 0 ; i < n_ ; i ++) sigma_0_inv[ sigma_0 [ i ] - 1 ] = i + 1;
    bool improve;
    do {//iterating use just sigma_0_inv
        int max_freq = 0;
        for (int i = 0 ; i < n_ ; i ++) for (int j = 0 ; j < n_ ; j++) same_cycle[ i][ j ] = 0;
        for (int s = 0 ; s < m; s++) {
            for(int j = 0 ; j < n_ ; j ++ ) samples_comp[ s ][ j ] = samples[ s ][ sigma_0_inv [ j ] - 1 ];
            get_cycles(samples_comp[ s ], cycle_items, cycle_index);
            //cout<<"Cycles ";gen.print_int_vector(cycle_items, n_); gen.print_int_vector(cycle_index, n_);
            for ( int i = 0 ; i < n_ ; i ++){
                for (int j = i + 1; j < n_ && (cycle_index[ i ] == cycle_index[ j ]) ; j++){
                    int max, min; // for a triangular matrix
                    if (cycle_index [ i ] > cycle_index[ j ]){ max = cycle_items[ i ] - 1; min  = cycle_items[ j ]-1;}
                    else {min  = cycle_items[ i ]-1; max  = cycle_items[ j ]-1;}
                    same_cycle[ min ][ max ] ++;
                    if (max_freq < same_cycle[ min ][ max ] ) {
                        max_freq = same_cycle[ min ][ max ];
                        index_i = min;
                        index_j = max;
                    }
                    //if (max_freq > m)//ERROR
                    //{ cout<<"Cycles ";gen.print_int_vector(cycle_items, n_); gen.print_int_vector(cycle_index, n_);gen.print_int_matrix(same_cycle, n_, n_);}
                }
            }
        }
        //there are m oermus, max_freq of them are going to decrease dist in 1. (m-max_freq) are going to incresse
        distance_variation = m - 2 * max_freq;
        
        improve = false;
        if( distance_variation < 0 ) {
            improve = true;
            int aux = sigma_0_inv[index_i];
            sigma_0_inv[index_i] = sigma_0_inv[index_j];
            sigma_0_inv[index_j] = aux;
            
            (*f_eval) += distance_variation;
            for(int i = 0 ; i < n_ ; i ++) sigma_0[ sigma_0_inv[ i ] - 1 ] = i + 1;
            
        }
        
        //cout<<"trace in local search swap, distance variation "<<distance_variation<<" mxa freq "<<max_freq<<endl;
    } while (improve );
    
    delete [] sigma_0_inv;
    delete [] cycle_index;
    delete [] cycle_items;
    for (int i = 0 ; i < n_ ; i++) delete [] same_cycle[ i ];
    delete [] same_cycle;
    for (int i = 0 ; i < m ; i++)  delete [] samples_comp[ i ];
    delete [] samples_comp;
}

void Cayley::local_search_swap_gmm(int m, int **samples, int *sigma_0, double *f_eval){
    int ** samples_comp = new int*[ m ];
    for (int s = 0 ; s < m; s++)  samples_comp[ s ] = new int[ n_ ];
    int *sigma_0_inv = new int[ n_ ];
    int *x_acum = new int [n_ ], *x_acum_var = new int[ n_], *x = new int[ n_], *best_x_acum = new int[ n_ ], *x_test = new int[ n_ ];
    double *theta = new double [ n_ ];
    Generic gen;
    int index_i,index_j;
    for(int i = 0 ; i < n_ ; i ++) {sigma_0_inv[ sigma_0 [ i ] - 1 ] = i + 1;x_acum[ i ] = 0;x_acum_var[ i ] = 0;}
    bool improve;
    for (int s = 0 ; s < m; s++){
        for(int j = 0 ; j < n_ ; j ++ )
            samples_comp[ s ][ j ] = samples[ s ][ sigma_0_inv [ j ] - 1 ];
        perm2dist_decomp_vector(samples_comp[ s ], x);
        for(int j = 0 ; j < n_ ; j ++ ) x_acum[ j ] += x[ j ];
    }

    do {//iterating use just sigma_0_inv
        double best_like = 0;
        for (int iter_i = 0; iter_i < n_ - 1 ; iter_i++) {
            for (int iter_j = iter_i + 1 ; iter_j < n_ ; iter_j++) {
                //test swap i j on samples_comp
                for(int k = 0 ; k < n_ ; k ++) {x_acum_var[ k ] = 0;}
                for (int s = 0 ; s < m ; s++) {
                    int max_i, max_j;
                    //cout<<"sample "; Generic gen; gen.print_int_vector(samples_comp[ s ], n_);
                    if (same_cycle(iter_i, iter_j, samples_comp[ s ])){
                        get_max_item_in_future_cycles(samples_comp[ s ], iter_i, iter_j, &max_i , &max_j);
                        x_acum_var[(max_i < max_j ? max_i : max_j )] --;
                        //cout<<" resta en "<<(max_i < max_j ? max_i : max_j )<<endl;
                    }else {
                        get_max_item_in_current_cycle(samples_comp[ s ], iter_i, &max_i );
                        get_max_item_in_current_cycle(samples_comp[ s ], iter_j, &max_j );
                        x_acum_var[(max_i < max_j ? max_i : max_j )] ++;
                        //cout<<" suma en "<<(max_i < max_j ? max_i : max_j )<<endl;
                    }
                }
                //////test
                double like = get_theta_log_likelihood(m, x_acum, x_acum_var, theta);
                if ( like > best_like || best_like == 0){
                    for(int i = 0 ; i < n_ ; i ++) best_x_acum[ i ] = x_acum_var[ i ];
                    best_like = like;
                    index_i = iter_i;
                    index_j = iter_j;
                }
            }
        }
        //there are m oermus, max_freq of them are going to decrease dist in 1. (m-max_freq) are going to incresse
        improve = false;
        

        if( best_like > *f_eval ) {
            improve = true;
            (*f_eval) = best_like;
            int aux = sigma_0_inv[index_i];
            sigma_0_inv[index_i] = sigma_0_inv[index_j];
            sigma_0_inv[index_j] = aux;
            for(int i = 0 ; i < n_ ; i ++) sigma_0[ sigma_0_inv[ i ] - 1 ] = i + 1;
            for (int s = 0 ; s < m ; s++) {
                int aux = samples_comp[ s ][index_i];
                samples_comp [ s ][index_i] = samples_comp[ s ][index_j];
                samples_comp [ s ][index_j] = aux;
            }
            for(int j = 0 ; j < n_ ; j ++ ) x_acum[ j ] += best_x_acum[ j ];
            //test a ver si son iguales
        }
        //cout<<"trace in local search swap, distance variation "<<distance_variation<<" mxa freq "<<max_freq<<endl;
    } while (improve );


    delete [] theta;
    delete [] best_x_acum;
    delete [] sigma_0_inv;
    delete [] x_acum;
    delete [] x_acum_var;
    delete [] x;
    delete [] x_test;
    for (int i = 0 ; i < m ; i++)  delete [] samples_comp[ i ];
    delete [] samples_comp;
}

void Cayley::local_search_insert (int m, int **samples, int *sigma_0, int model, double *f_eval){
    int *best_sol = new int[ n_ ], *next_sol = new int[ n_ ], *x_acumul = new int[ n_ ], *x = new int[ n_ ];
    double * theta = new double[ n_ ];
    double best_eval= 0;
    bool better;
    Generic gen;
    do {
        better = false;
        best_eval = 0;
        for(int i = 0 ; i < n_ ; i ++){
            for ( int j = 0 ; j < n_ ; j++){
                if(i != j){
                    gen.insert_at (sigma_0, n_, i, j, next_sol);
                    if(model == MALLOWS_MODEL){
                        int dist = distance_to_sample(samples, m , next_sol);
                        if ( dist < best_eval || best_eval == 0){
                            best_eval = (double )dist;
                            for (int i = 0 ; i < n_ ; i ++) best_sol [ i ] = next_sol[ i ];
                        }
                    }else {
                        for(int index = 0 ; index < n_ ; index ++ ) x_acumul[ index ] = 0;
                        for (int s = 0 ; s < m; s++){
                            perm2dist_decomp_vector(samples[ s ], x);
                            for(int index = 0 ; index < n_ ; index ++ ) x_acumul[index] += x[index];
                        }
                        double likeli = get_theta_log_likelihood(m,x_acumul, NULL, theta);
                        if (likeli > best_eval || best_eval == 0) {
                            best_eval = likeli;
                            for (int i = 0 ; i < n_ ; i ++) best_sol [ i ] = next_sol[ i ];
                        }
                    }
                }
            }
        }
        
        if( (model == MALLOWS_MODEL && best_eval < *f_eval ) || ( model == GENERALIZED_MALLOWS_MODEL && best_eval > *f_eval ) ){
            *f_eval = best_eval;
            for (int i = 0 ; i < n_ ; i ++) sigma_0 [ i ] = best_sol[ i ];
            better = true;
        }
    } while (better);
    delete [] next_sol;
    delete [] theta;
    delete [] x;
    delete [] x_acumul;
    delete [] best_sol;
}
