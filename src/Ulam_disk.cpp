//
//  Ulam_disk.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 31/01/14.
//  Copyright (c) 2014 Ekhine Irurozki. All rights reserved.
//


#include "Ulam_disk.h"

void    Ulam_disk::read_permus_per_dist(){
    if ( num_permus_per_dist_[ 0 ] == 0 ){//not initialized
        //num_permus_per_dist_num_permus_per_dist_ = new long double[ n_ ];
        char integer_string[5];
        sprintf(integer_string, "%d", n_);
        char str_permus_per_dist[600];
        strcpy(str_permus_per_dist, str_base_path);
        strcat(str_permus_per_dist, "permus_per_dist_");
        strcat(str_permus_per_dist, integer_string);
        ifstream  file;
        file.open(str_permus_per_dist);
        if(!((file ))){
            //cout << "Cannot read input file. Have you generated the files? "<<endl;
            return;
        }
        for (int i = 0 ; i < n_ ; i ++)
            file >> num_permus_per_dist_[ i ];
        file.close();
    }
}

void    Ulam_disk::read_permus_per_shape(int d, int bound, int * shape, int*shape_len){
    string line;
    long double acum;
    bool end = false;
    char integer_string[5];
    
    sprintf(integer_string, "%d", n_);
    char str_permus_per_shape[500] ;// = "/Users/eki/Dropbox/permus/prj/perms_mallows/permus_per_shape_";
    strcpy (str_permus_per_shape, str_base_path);
    strcat(str_permus_per_shape, "permus_per_shape_");
    strcat(str_permus_per_shape, integer_string); 
    strcat(str_permus_per_shape, "_");
    sprintf(integer_string, "%d", d);
    strcat(str_permus_per_shape, integer_string);
    
    ifstream  file;
    file.open(str_permus_per_shape);
    if(!((file ))){
        //cout  << "Cannot read input file - permus_per_dist "<<endl;
        //exit(1);
        return;
    }
    do{
        file >> acum ;//read the first int: acumul number of permutations at distance d
        //cout<<acum<<endl;
        if ( acum < bound )//if not reached, discars theline
            file.ignore(n_ * 2 , '\n');
        else {
            getline(file, line);
            end = true;
        }
    }while ( ! end );
    file.close();
    
    stringstream tokenizer( line );//otherwise read the associated partition from the rest of the line
    int lambda_i;
    while( tokenizer >> lambda_i ){
        //cout << "Token: " << lambda_i << endl;
        shape[ (*shape_len)++ ] =lambda_i;
    }
}

void Ulam_disk::save_counts_to_file(){
    /* adapted from ZS1: Zoghbi, Antoine and Stojmenovic, Ivan: Fast Algorithms for Generating Integer Partitions. International Journal of Computer Mathematics, 70, 1998, 319-332.
     generaes partitions in  anti-lexicographic order */
    ofstream file_permus_per_dist;
//    ofstream file_shapes_per_dist;
    ofstream file_permus_per_shape;
    
    char integer_string[5];
    sprintf(integer_string, "%d", n_);
    char str_permus_per_shape[500] ;//= "/Users/eki/Dropbox/permus/prj/perms_mallows/permus_per_shape_"; // un file e estos por cada n y distance
    strcpy(str_permus_per_shape, str_base_path);
    strcat(str_permus_per_shape, "permus_per_shape_");
    
    char str_permus_per_shape_n_d[600];
    char str_permus_per_dist[600] ;//= "/Users/eki/Dropbox/permus/prj/perms_mallows/permus_per_dist_";
    strcpy(str_permus_per_dist, str_base_path);
    strcat(str_permus_per_dist, "permus_per_dist_");

    strcat(str_permus_per_shape, integer_string);
    strcat(str_permus_per_dist, integer_string);
    
    unsigned char k;                  //length of figures
    unsigned char *vector     = NULL; //where the current figure is stored
    int           gen_result;         //return value of generation functions
    int           prev_distance = -1 , dist ;
    int           part_len;
    long double     permus_per_shape = 0;
    vector = (unsigned char *)malloc(sizeof(unsigned char) * n_ );
    if(vector == NULL)    {
        //fprintf(stderr, "error: insufficient memory\n");
        //exit(EXIT_FAILURE);
        return;
    }
    for (int i = 0 ; i < n_ ;i++) num_permus_per_dist_[ i ] = 0 ;
    
    Ferrers_diagram*f;
    gen_result = gen_part_init(vector, n_, &k);
    while(gen_result == GEN_NEXT ) {
        part_len = (int)k;
        int*part = new int[part_len];//DO NOT delete, member of Ferrers_diagram
        for(int i = 0 ; i < part_len; i++) part[i]=(int)vector[i];
        dist = part[ 0 ];
        
        f = new Ferrers_diagram(n_, part , part_len);
        f->calculate_hook_length(facts_[ n_ ]);
        dist = f->get_resulting_distance();
        num_permus_per_dist_[ dist ] += f->get_num_permus();
        if ( dist != prev_distance ){
            permus_per_shape = 0;
            //cout<<"Generating shape at distance "<<dist<<endl;
            if (file_permus_per_shape.is_open() )
                file_permus_per_shape.close();
            strcpy(str_permus_per_shape_n_d, str_permus_per_shape);
            strcat(str_permus_per_shape_n_d, "_");
            sprintf(integer_string, "%d", dist );
            strcat(str_permus_per_shape_n_d, integer_string);
            file_permus_per_shape.open(str_permus_per_shape_n_d);
        }
        permus_per_shape += f->get_num_permus();
        file_permus_per_shape << permus_per_shape <<" ";
        for ( int i =  0; i < part_len ; i++ )
            file_permus_per_shape<<part[ i ]<<" ";
        file_permus_per_shape<<endl;
        prev_distance = dist;
        gen_result = gen_part_next(vector, &k, 0);
        delete f ;
    }
    free(vector);
    
    file_permus_per_shape.close();
    file_permus_per_dist.open (str_permus_per_dist);
    for (int i = 0 ; i < n_; i++ ){
        file_permus_per_dist<< num_permus_per_dist_[ i ]<<endl ;
    }
    //delete [] num_permus_per_dist_;
    file_permus_per_dist.close();
}

double Ulam_disk::expectation(double theta){
    read_permus_per_dist();//OJO WARNING - si es de clase
    long double numer = 0, denom = 0;
    for (int d = 0 ; d < n_ - 1; d++){
        long double aux = num_permus_per_dist_[ d ] * exp(-theta *d ) ;
        numer += aux * d;
        denom += aux;
    }
    return (double)numer / denom;
}

void Ulam_disk::distances_sampling(int m, double theta, int **samples){
    double  rand_distance = 0;
    double  *proba        = new double[ n_ ];
    int     target_distance;
    
    read_permus_per_dist();//OJO WARNING - si es de clase
    proba[ 0 ] = 1; // exp(-theta*d) = exp (-theta *0)
    for (int i = 1 ; i < n_ ; i++)//acumulate the number of permus at each distance
        proba[i] = num_permus_per_dist_[i] * exp ( -theta * i ) + proba[ i - 1 ];
    for (int i = 0 ; i < m ; i ++){
        rand_distance = (double) rand() / (double)(RAND_MAX) * proba[ n_ - 1 ];
        target_distance = 0;
        while(proba[ target_distance ] <= rand_distance) target_distance++;
        samples[ i ] = new int[ n_ ];
        //cout<<"ulam distance greene_niejenhuis_wilf "<<target_distance<<endl;
        generate_permu_with_given_LIS( n_ - target_distance, samples[i]);
        //distance_acum += (double)target_distance;
    }
    //cout<<"average distance (dist_sampling) "<<distance_acum/m<<endl;
    delete [] proba;
    
}


void Ulam_disk::generate_permu_with_given_LIS (int l, int *sigma){
    int     d = n_ - l;
    int     to_insert;
    int     col, row, aux, new_col, new_row;
    int     *col_index = new int[n_], *row_index = new int[n_];
    
    read_permus_per_dist();//OJO WARNING - si es de clase
    double target_permus = (double)rand() / (double)(RAND_MAX) * num_permus_per_dist_[ d ];
    int * shape = new int[ n_ ] , shape_len = 0;
    read_permus_per_shape(d , target_permus, shape, &shape_len);
    
    int* shape1 = new int [ shape_len ];//member of f1
    int* shape2 = new int [ shape_len ];
    for (int i = 0 ; i < shape_len ; i ++){
        shape1[ i ] = shape[ i ];
        shape2[ i ] = shape[ i ];
    }
    Ferrers_diagram* f1 = new Ferrers_diagram(n_   , shape1 , shape_len);
    Ferrers_diagram* f2 = new Ferrers_diagram(n_   , shape2 , shape_len);
    
    f1->random_SYT();
    f2->random_SYT();
    int ** tableau1 = f1->get_syt();//
    int ** tableau2 = f2->get_syt();
    
    for (int i = 0 ; i < f2->get_ferrers_shape_length() ; i++){
        for (int j =  0 ; j < f2->get_ferrers_shape()[i] ; j++) {
            row_index[ tableau2[ i ][ j ] - 1 ] = i;
            col_index[ tableau2[ i ][ j ] - 1 ] = j;
        }
    }
    
    for (int index = n_ - 1 ; index >= 0 ; index --){
        col = col_index[ index ];
        row = row_index[ index ];
        to_insert = tableau1[ row ][ col ];
        while (row != 0) {
            new_col=0, new_row = row - 1;
            while (f1->get_ferrers_shape()[new_row] > new_col+1
                   && tableau1[new_row][new_col + 1 ] < to_insert)
                new_col++;
            aux = tableau1[new_row][new_col];
            tableau1[new_row][new_col] = to_insert;
            to_insert = aux;
            row = new_row;
            col = new_col;
        }
        sigma[index ] = to_insert;
        tableau1[row_index[ index ]][ col_index[ index ]] = n_ + 1;
        //gen.print_int_vector(sigma, n_);
    }
    delete [] col_index;
    delete [] row_index;
    delete [] shape;
    delete f1;
    delete f2;
}


long double Ulam_disk::get_likelihood(int m, int **samples, int model, int *sigma_0){
    Newton_raphson newton(n_);
    long double likelihood, psi = 0;
    int     dist_avg = 0 ;
    double  theta;
    
    read_permus_per_dist();//OJO WARNING - si es de clase
    dist_avg = distance_to_sample(samples, m, sigma_0);
    theta    = newton.Newton_raphson_method( (double)dist_avg/m, -1.001,ULAM_DISTANCE, NULL, MALLOWS_MODEL, num_permus_per_dist_);
    for (int i = 0 ; i < n_ ; i ++ ) psi += num_permus_per_dist_[ i ] * exp (-theta * i );
    likelihood = - dist_avg * theta - m* log ( psi );
    return likelihood;

}

void Ulam_disk::estimate_theta(int m, int *sigma_0, int **samples, double *theta){
    Newton_raphson newton(n_);
    read_permus_per_dist();//OJO WARNING - si es de clase
    int dist_avg = distance_to_sample(samples, m, sigma_0);
    *theta = newton.Newton_raphson_method((double)dist_avg/m, -1.001,ULAM_DISTANCE, NULL, MALLOWS_MODEL, num_permus_per_dist_);
}















