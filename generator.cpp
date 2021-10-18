#include "generator.h"
#include <random>
#include <ctime>
#include <iostream>
#include <mutex>

std::mt19937_64 engine(std::time(nullptr));
std::uniform_real_distribution<double> distribution{0, 1};
std::ofstream crds ("C:\\Qt\\progects\\En_dep\\cardS.txt");
std::ofstream crds_ ("C:\\Qt\\progects\\En_dep\\cardS1.txt");
std::ofstream neas ("C:\\Qt\\progects\\En_dep\\NEAS.txt");
constexpr double pi{3.14159265358979323846};

std::vector<double> all2;

void Generator::randomcp(){
    double gamma1{}, gamma2{};
    gamma1 = distribution(engine);
    gamma2 = distribution(engine);
    x_core = radius*std::sqrt(gamma1)*std::cos(gamma2*2*pi);
    y_core = radius*std::sqrt(gamma1)*std::sin(gamma2*2*pi);
}

double Generator::gamma_dep_inside(double &E_Gev){
    if(E_Gev < 0.01) return 0.057*std::pow(pie, -0.2);
    else if (E_Gev >= 0.01) return 0.025*std::pow(pie, 0.16);
}

void Generator::cherenok(double &distance, double &pie, int&pid, int& i){
    if((distance < 0.1 && pie > 10 && pid <= 3) || (distance < 0.1 && pie > 100 && pid >=5)) en_deposit[i] += 2;
}

double Generator::electrons_dep_inside(double &pie){
    if(pie < 10 && pie > 1) return 3e-4*std::pow(pie, 3.5);
    else if(pie >=10) return std::pow((pie/10), 0.1);
    else if(pie < 1) return 3e-4*E_Gev;
}

double Generator::neutrons_dep_inside(double &E_Gev){
    if(E_Gev < 0.1 && E_Gev > 0.01) return 2.3*std::pow((E_Gev/0.1), 0.67);
    else if(E_Gev >= 0.1 && E_Gev < 100) return 2.3*std::pow((E_Gev/0.1), 0.2);
    else if(E_Gev >= 100) return 23.3*std::sqrt(E_Gev/1000);
}

double Generator::protons_dep_inside(double &E_Gev){
    if(E_Gev > 0.045 && E_Gev < 0.2) return 0.5/E_Gev;
    else if(E_Gev >= 0.2 && E_Gev < 100) return 2.3*std::pow((E_Gev/0.1), 0.2);
    else if(E_Gev >= 100) return 23.3*std::sqrt(E_Gev/1000);
}

double Generator::gamma_dep_outside(double &E_Gev){
    if(E_Gev < 0.0215) return 300 * std::pow(E_Gev, 4);
    else if (E_Gev < 0.1 && E_Gev >= 0.0215) return 3e-6*std::pow(E_Gev, -0.8);
    else if (E_Gev < 0.3 && E_Gev >= 0.1) return 0.0025*std::pow(E_Gev, 2);
    else if(E_Gev >= 0.3) return 4e-4*std::sqrt(E_Gev);
}

void Generator::process(){
    if (pid > 65) return;
    E_Gev = pie/1000;
    X_shift = pix/100 + x_core;
    Y_shift = piy/100 + y_core;

    if(std::sqrt(X_shift*X_shift+Y_shift*Y_shift) < 8) Ne+=1;                               // full Ne w/o gammas within R<Rm

    if(pid > 6 && std::sqrt(X_shift*X_shift+Y_shift*Y_shift) < 1) Number_hadrons_core+=1;   // hadrons inside r=1m : the core

    if (std::abs(X_shift) < 20 && std::abs(Y_shift) < 20){
        if (pid ==13) Number_neutrons_circle+=1;
        if(pid > 6) Number_hadrons_circle+=1;
    }
    else return;                                                        // around dets in CARPET, R=20 m

    double distance{}, eps{};
    for(int i = 0; i < Number_det; ++i)
    {
        eps = 0;
        distance = std::sqrt((X_shift - detx[i])*(X_shift - detx[i]) + (Y_shift - dety[i])*(Y_shift - dety[i]));
        if(distance <= 0.35){
            if(pid == 1){                                                   // gammas
                eps = gamma_dep_inside(E_Gev);
                en_dep(eps, distance, i);
            }
            if(pid == 2 || pid == 3){                                       // e+-
                eps = electrons_dep_inside(pie);
                cherenok(distance,pie,pid,i);
                en_dep(eps, distance, i);
            }
            if((pid == 5 || pid == 6) && E_Gev > 0.1){                          // muons
                eps = std::pow((E_Gev/0.3), 0.1);
                cherenok(distance, pie, pid, i);
                if(E_Gev > 1) muon_det[i] += 1;
                Number_muons+=1;
                en_dep(eps, distance, i);
            }
            if(pid == 13){                                                  // neutrons
                eps = neutrons_dep_inside(E_Gev);
                en_dep(eps, distance, i);
            }
            if(pid == 14){                                                  // protons
                eps = protons_dep_inside(E_Gev);
                en_dep(eps, distance, i);
            }
        }
        else en_dep(eps, distance, i);
    }
}

void Generator::en_dep(double &eps, double &distance, int &i){
    double gamma{};
    gamma = distribution(engine);
    if(gamma <= 0.8) en_deposit[i] += eps*1.25;
    double pg{}, p{};
    int N{};
    if(distance<=40){
        if(pid > 7 || (pid == 1 && pie > 10)){
            if(pid == 1) pg = gamma_dep_outside(E_Gev);                            // gammas inside 20 m around det.
            else pg = 0.0125*std::pow(E_Gev, 0.45);                                // for hadrons from GEANT
            p = std::exp(-distance/0.45)*pg*0.635;                                 // due to pulse selection efficiency
            if(p < 5){
                std::poisson_distribution<int> dis(p);
                N = dis(engine);
            }
            else N = round(0.5+p);
            integ_par[i] += N;                                                            // No of n (atmospheric n added)
            Number_neutrons+=N;
        }
    }
}

void Generator::output(double Ns, double ECRTOT, double THETACR, double PHICR){
    int *trigged_det = new int[4];
    double *number_par_adc = new double[64];
    double *en_sum = new double[4];
    int  M1{}, M2{}, N{}, M{};
    double p{}, lg_Ne{};

    for(int j = 0; j < 4; ++j){
        for(int i = 0; i < 16; ++i){
            p = en_deposit[i+16*j]/par_per_channel;
            if(p < 10){
                std::poisson_distribution<int> dis(p);
                N = dis(engine);
            }
            else N = round(0.5+p);

            en_deposit[i+16*j] = number_par_adc[i+16*j] = N;
            if(en_deposit[i+16*j] >= threshold) trigged_det[j]+=1;

            en_sum[j] += number_par_adc[i+16*j];
            if(integ_par[i+16*j] > 999) integ_par[i+16*j] = 999;
            if(number_par_adc[i+16*j] > 99999) number_par_adc[i+16*j] = 99999;
        }
        if(trigged_det[j] >= krat) ++M1;
        if(en_sum[j] >= 300) ++M2;
  }
        if(M1 >=1) M = 1;
        if(M2 >= 1) M += 2;
        if(Number_neutrons >=10) M += 4;
        if(M1 > 1 || M2 > 1) M += 8;

        if(M < 1){
            Number_neutrons = Number_muons = Number_hadrons_circle = M = Ne = 0;
            for(int i = 0; i < Number_det; ++i) en_deposit[i] = number_par_adc[i] = integ_par[i] = muon_det[i] = 0;
        }
        else{
        ++Number_event;

        if(Ne > 2) lg_Ne = std::log10(Ne);
        else lg_Ne = 1;

        if(Number_neutrons > 2500) Number_neutrons = 2500;

        std::cout << "Ns = " << Ns << " evs = " << Number_event << " M = " << M << " n = " << Number_neutrons << " mu= " << Number_muons
                  << " E= " << ECRTOT/1e3 << ' ' << "Number_hadrons_circle= " << Number_hadrons_circle << '\n';

        //crds_ << x_core << "," << y_core << "," << Ns  << "," << Number_neutrons  << "," << Number_muons << "," << M << ",";
        for(int j = 0; j < 64; ++j) crds_ << number_par_adc[j] << "," << integ_par[j] << "," << muon_det[j] << ",";

        all.insert(all.end(), {x_core, y_core, Ns, static_cast<double>(Number_neutrons), static_cast<double>(Number_muons), static_cast<double>(M)});
        for(int j = 0; j < Number_det; ++j){
            all.push_back(number_par_adc[j]);
            all.push_back(integ_par[j]);
            all.push_back(muon_det[j]);
        }
        all.insert(all.end(), {lg_Ne, ECRTOT/1e3, THETACR, PHICR, static_cast<double>(Number_hadrons_circle)});

        //crds_ << lg_Ne << "," << ECRTOT/1e3 << "," << THETACR << "," << PHICR << "," << Number_hadrons_circle << '\n';

        //neas << "N_EAS = " << Ns << ' ' << "evs = " << N_event << ' ' << "am0 = " << am0 << " Nd = " << Ndet << '\n';
        }
        Number_neutrons = Number_muons = Number_hadrons_circle = M = Ne = 0;
        for(int i = 0; i < Number_det; ++i) en_deposit[i] = number_par_adc[i] = integ_par[i] = muon_det[i] = 0;
}

void Generator::print_all(){
    for(int i = 0; i < all.size(); ++i){
        if((i+1)%203!=0) crds << all[i] << ",";
        if((i+1)%203==0){
            crds << all[i] <<'\n';
        }
    }
}
