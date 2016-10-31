//
//  main.cpp
//  Class tutorial
//
//  Created by joss briody on 01/10/2014.
//  Copyright (c) 2014 jossbriody. All rights reserved.
//

#include <iostream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <iostream>
#include <vector>
#include <math.h>


std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);
std::uniform_real_distribution<> dis2(0, 100);

double random_num = dis(gen); // only one random num per run
double random_num2 = dis2(gen);
const double T = 300;

class lattice
{
friend class lattice_energy;
friend class metropolis;
private:
    int height;
    int width;
    std::vector<int> spins;
public:
    void setwidth( int wid ) { width = wid; }
    void setheight( int hei ) { height = hei; }
    int getsize() { return width * height; }
    void initialise_lattice( std::vector<int> &init_spins ) {
        int current = 0;
        for(auto &idx : init_spins){
            if( round( dis( gen) ) == 0 )
                idx =  -1;
            else
                idx = 1; // rounding up may cause issue?
            ++current;
        }
        spins = init_spins;
    }
    void setspins( int lattice_idx, int new_spin ) { spins[ lattice_idx ] = new_spin; }
    
    std::vector<int> printspins() {
        return spins;
    }
    lattice operator + (const lattice& l )
    {
        lattice Lattice;
        Lattice.width = this->width + l.width;
        Lattice.height = this->height + l.height;
        return Lattice;
    }
};

class lattice_energy   // make this the child of the class lattice - use a virtual function
{
friend class lattice;
friend class metropolis;
public:
    lattice_energy( lattice );
    double sum_self_energy( lattice );
    double calculate_self_energy( int, lattice );
    double sum_interaction_energy( lattice  );
    double calculate_interaction_energy( int, lattice );
    double printenergy() { return energy; }
    static double J;
    static double Moment;
    static double field_strength;
private:
    double interaction_energy;
    double self_energy;
    double energy;
};

double lattice_energy::sum_interaction_energy( lattice ensemble ) {//std::vector<int> spins ) {
    double tot_interactionE = 0.0; // in units of J
    int i = 0;
    while( i < ensemble.spins.size() )
    {
        tot_interactionE += calculate_interaction_energy( i, ensemble );
        i++;
    };
    return (- J * tot_interactionE ) / 2;
}
double lattice_energy::calculate_interaction_energy( int i,  lattice ensemble ) {
    int x_idx = i % ensemble.width;
    int y_idx = i % ensemble.height;
    double interactionE;
    if( x_idx == 0 || x_idx == ensemble.width - 1 || y_idx == 0 || y_idx == ensemble.height - 1 ) // continuous b.c ( i.e surface of a sphere)
    {
        if( y_idx != 0 && y_idx != ensemble.height - 1 )             //x edge
        {
            x_idx == 0 ? interactionE = ( ( ensemble.spins[ i ] * ensemble.spins[ i + 1 ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i + ( ensemble.width - 1 ) ] )
                                          + ( ensemble.spins[ i ] * ensemble.spins[ i + ensemble.width ] )  + ( ensemble.spins[ i ] * ensemble.spins[ i - ensemble.width ] ) )
            :
            interactionE += ( ( ensemble.spins[ i ] * ensemble.spins[ i - ( ensemble.width - 1 ) ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i - 1 ] )
                             + ( ensemble.spins[ i ] * ensemble.spins[ i + ensemble.width ] )  + ( ensemble.spins[ i ] * ensemble.spins[ i - ensemble.width ] ) );
            
        }
        else if( x_idx != 0 && x_idx != ensemble.width - 1 )         //y edge
        {
            y_idx == 0 ? interactionE = ( ( ensemble.spins[ i ] * ensemble.spins[ i + 1 ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i - 1 ] )
                                          + ( ensemble.spins[ i ] * ensemble.spins[ i + ensemble.width ] )  + ( ensemble.spins[ i ] * ensemble.spins[ i + ( ensemble.height - 1 ) * ensemble.width ] ) )
            :
            interactionE += ( ( ensemble.spins[ i ] * ensemble.spins[ i + 1 ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i - 1 ] )
                             + ( ensemble.spins[ i ] * ensemble.spins[ i - ensemble.width ] )  + ( ensemble.spins[ i ] * ensemble.spins[ i - ( ( ensemble.height - 1 ) * ensemble.width ) ] ) );
        }
        else if( x_idx == 0 && y_idx == 0)                          // top-left corner
        {
            interactionE = ( ensemble.spins[ i ] * ensemble.spins[ i + 1 ] ) + ( ensemble.spins[ i ] *  ensemble.spins[ i + ensemble.width ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i + ( ensemble.width - 1 ) ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i + ( ensemble.height - 1 ) * ensemble.width ] ) ;
            
        }
        else if( x_idx == 0 && y_idx == ensemble.height - 1 )       //bottom-left corner
        {
            interactionE = ( ensemble.spins[ i ] * ensemble.spins[ i + 1 ] ) + ( ensemble.spins[ i ] *  ensemble.spins[ i - ( ( ensemble.height - 1 ) * ensemble.width ) ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i + ( ensemble.width - 1 ) ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i - ensemble.width ] ) ;
        }
        else if( x_idx == ensemble.width - 1 && y_idx == 0)         //top-right corner
        {
            interactionE = ( ensemble.spins[ i ] * ensemble.spins[ i - 1 ] ) + ( ensemble.spins[ i ] *  ensemble.spins[ i + ( ( ensemble.height - 1 ) * ensemble.width ) ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i - ( ensemble.width - 1 ) ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i + ensemble.width ] ) ;
        }
        else // This is the case: ( x_idx == ensemble.width - 1 && y_idx == ensemble.height - 1 )  ie. bottom-right corner
        {
            interactionE = ( ensemble.spins[ i ] * ensemble.spins[ i - 1 ] ) + ( ensemble.spins[ i ] *  ensemble.spins[ i - ( ( ensemble.height - 1 ) * ensemble.width ) ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i - ( ensemble.width - 1 ) ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i - ensemble.width ] ) ;
        }
    }
    else // in the middle of the lattice
    {
        interactionE = ( ensemble.spins[ i ] * ensemble.spins[ i + 1 ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i - 1 ] ) + ( ensemble.spins[ i ] * ensemble.spins[ i + ensemble.width ] )  +
        ( ensemble.spins[ i ] * ensemble.spins[ i - ensemble.width ] );
    }
    return interactionE;
    //std::cout<< x_idx << "\t" << y_idx << std::endl;
}

double lattice_energy::sum_self_energy( lattice ensemble )
{
    std::vector<int> spins = ensemble.spins;
    double sum_self_energy = 0.0;
    int i = 0;
    while( i < spins.size() )
    {
        sum_self_energy += calculate_self_energy( i, ensemble );
        i++;
    }
    return sum_self_energy;
}

double lattice_energy::calculate_self_energy( int i, lattice ensemble )
{
    std::vector<int> spins = ensemble.spins;
    return ( - Moment * field_strength * spins[ i ] );
}

lattice_energy::lattice_energy( lattice ensemble ) : interaction_energy( sum_interaction_energy( ensemble ) ), self_energy( sum_self_energy( ensemble ) ) {
    energy = interaction_energy + self_energy; }

double lattice_energy::J = 1.0;
double lattice_energy::Moment = 1.0;
double lattice_energy::field_strength = 1.0;

class metropolis {
    friend class lattice;
    friend class lattice_energy;
private:
    std::vector<int> temp_spin;
    public:
    metropolis( lattice, lattice_energy );
    double calculate_dE( lattice ensemble, lattice_energy energy ) {
        int spin_idx = round( dis2( gen ) );
        double dE = - 2 * energy.calculate_interaction_energy( spin_idx , ensemble);
        return dE;
    }
    void flip_spin( double dE, int spin_idx, std::vector<int> & );
    std::vector<int> iterate_lattice( lattice ensemble, lattice_energy ); //std::vector<int> &spin
    void update_spin( lattice & );
};

metropolis::metropolis( lattice ensemble, lattice_energy energy) : temp_spin( iterate_lattice( ensemble, energy ) ) { }

void metropolis::flip_spin( double dE, int spin_idx, std::vector<int> &spins) {
    double boltzman_num = dis(gen);
    if( dE < 0 )
        spins[ spin_idx ] *= - 1;
    else if( exp( - dE / T) > boltzman_num )
        spins[ spin_idx ] *= - 1; //in units of Kb
}

std::vector<int> metropolis::iterate_lattice( lattice ensemble, lattice_energy energy  ) {
    int count = 0;
    std::vector<int> spins = ensemble.spins;
    while( count < spins.size() )
    {
        double dE = calculate_dE( ensemble, energy );
        flip_spin( dE, count, spins);
        count ++;
    }
    return spins;
}

void metropolis::update_spin( lattice &ensemble ) {
//    for( auto &s : ensemble.spins )
//        std::cout << s << std::endl;
    ensemble.spins = temp_spin;
//    std:: cout << "new spins" << std::endl;
//    for( auto &s2 : ensemble.spins )
//        std::cout << s2 << std::endl;
}

struct historic_states {
    static int time_steps;
    std::vector<double> Energies;
    std::vector< std::vector< int > >  Spin_states;
};

int historic_states::time_steps = 10;

int main()
{
    historic_states States;
    
    lattice lattice1;
    lattice1.setwidth( 10 );
    lattice1.setheight( 10 );
    
    std::vector<int> init_spins( lattice1.getsize(), 0 );
    std::vector<int> new_spins = init_spins;
    
    lattice1.initialise_lattice( init_spins) ;      //initialise spins randomly up or down
    
    for( int i = 0; i < States.time_steps; i++ )
    {
        lattice_energy energy( lattice1 );
        States.Spin_states.push_back( lattice1.printspins() );
        States.Energies.push_back( energy.printenergy() );
        metropolis evolve( lattice1, energy );
        evolve.update_spin( lattice1 );
    }
        
    return 0;
}

