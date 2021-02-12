#include <iostream>
#include <fstream>
#include <random>
#include <eigen3/Eigen/Dense>
#include <cmath>

struct kuramoto{
    const unsigned N;
    double dt;
    double K;
    Eigen::VectorXd ts;
    Eigen::Array<double, -1, -1> ks;
    Eigen::Array<double, -1, -1> omegas;
    Eigen::Array<double, -1, -1> thetas;
    Eigen::Array<double, -1, -1> theta_prime;
    std::mt19937 rnd;
    std::normal_distribution<double> normal_dist{0., 1. };
    std::uniform_real_distribution<double> uniform_dist{0., 1.};


// INITIALIZATION FUNCTIONS
    kuramoto(unsigned n, double TT, unsigned size_ts) :
        N(n), omegas(n, 1), thetas(n, 1), theta_prime(n,1), ks(1, 1){
            std::random_device r;
            std::seed_seq seed2{r(), r(), r(), r(), r(), r(), r(), r()};
            rnd.seed(seed2); 

            ts.setLinSpaced(size_ts,0,TT);
            dt = ts(1) - ts(0);
            ks(0) = 1;
            K = 1;

            for(int i = 0; i < thetas.rows(); i++){
                thetas(i) = uniform_dist(rnd);
            }
        }

    inline double get_gaussian(){
      return normal_dist(rnd);
    }
    
    inline double get_uniform(){
      return uniform_dist(rnd);
    }

    void gaussian_omegas(){
        for(int i = 0; i < thetas.rows(); i++){
            omegas(i) = get_gaussian();
        }
    }

    void uniform_omegas(){
        for(int i = 0; i < thetas.rows(); i++){
            omegas(i) = get_uniform();
        }
    }

    void ks_linspace(unsigned size_ks, double k0, double k1){

        Eigen::Array<double, -1, -1> newks(size_ks,1);
        Eigen::ArrayXf v = Eigen::ArrayXf::LinSpaced(size_ks, k0, k1);
            for(int i = 0; i < v.rows(); i++){
                newks(i) = v[i];
            }
        ks = newks;
    } 
//

    void evolve(){
        update_theta_prime();
        thetas += dt * theta_prime;
    }

    void update_theta_prime(){
        for(int i = 0; i < theta_prime.rows(); i++){
            theta_prime(i) = 0;
            for(int j = 0; j < theta_prime.rows(); j++){
                theta_prime(i) += sin(thetas(j)-thetas(i)) ;
            }
        }
        theta_prime /= (K / N);
        theta_prime += omegas;
    }

};

int main(){


// exercise 1
    kuramoto kura(10, 100, 10001);
    kura.gaussian_omegas();
    kura.ks_linspace(26, 0, 5);
    std::cout << kura.thetas;
    kura.evolve();
    std::cout << kura.thetas;

}

