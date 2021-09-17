#include "generator.h"
#include <random>
#include <ctime>
#include <iostream>

std::mt19937_64 engine(std::time(nullptr));
std::uniform_real_distribution<double> distribution{0, 1};

void Generator::randomcp(){
    double gamma1{}, gamma2{};
    gamma1 = distribution(engine);
    gamma2 = distribution(engine);
    x0 = r0*(2*gamma1 - 1.);
    y0 = r0*(2*gamma2 - 1.);
//  x0 = 31.4;                 // для теста -37.2	-18.5
//  y0 = 19.4;                 // 31.4  19.4
}

int Generator::poisson(double &p, int &N){
    double gamma{}, TR{1};
    N = 0;
    gamma = distribution(engine);
    TR = gamma;
    while(TR - std::exp(-p) > 0){
        gamma = distribution(engine);
        TR*=gamma;
        N+=1;
    }
    return N;
}

void Generator::process(){

    if (pid > 65) return;
    E = pie/1000;
    X = pix/100 + x0;
    Y = piy/100 + y0;
    RH = std::sqrt(X*X+Y*Y);

    if(RH < 8) En+=1;                               // full Ne w/o gammas within R<Rm
    EL = E;                                         // ??? это в if или нет?

    if(pid > 6 && RH < 1) Nh_core+=1;                    // hadrons inside r=1m : the core

    if (std::abs(X) < 20 && std::abs(Y) < 20){
        if (pid ==13) Nn_sq+=1;
        if(pid > 6) Nn_sq+=1;
    }
    else return;                                    // around dets in CARPET, R=20 m

    double rdx{}, eps{};
    for(int i = 0; i < 64; ++i)
    {
        rdx = std::sqrt((X - detx[i])*(X - detx[i]) + (Y - dety[i])*(Y - dety[i]));
        if(rdx <= 0.35){
            if(pid == 1){                                                   // gammas
                if(E < 0.01) eps = 0.057*std::pow(pie, -0.2);
                else if (E >= 0.01) eps = 0.025*std::pow(pie, 0.16);
                en_dep(eps, rdx, i);
            }
            if(pid == 2 || pid == 3){                                       // e+-
                if(pie < 10 && pie > 1) eps = 3e-4*std::pow(pie, 3.5);
                else if(pie >=10) eps = std::pow((pie/10), 0.1);            // expectation for e in mip
                else if(pie < 1) eps = 3e-4*E;
                if(rdx < 0.1 && pie > 10) ed[i] += 2;                       // Cherenok in glass = 5 prtcls
                en_dep(eps, rdx, i);
            }
            if((pid == 5 || pid == 6) && E > 0.1){                          // muons
                eps = std::pow((E/0.3), 0.1);
                if(rdx < 0.1 && pie > 100) ed[i] += 2;                      // Cherenok in glass = 5 prtcls
                if(E > 1) nedmu[i] += 1;
                Nm+=1;
                en_dep(eps, rdx, i);
            }
            if(pid == 13){                                                  // neutrons
                if(E < 0.1 && E > 0.01) eps = 2.3*std::pow((E/0.1), 0.67);
                else if(E >= 0.1 && E < 100) eps = 2.3*std::pow((E/0.1), 0.2);
                else if(E >= 100) eps = 23.3*std::sqrt(E/1000);             //sqrt(pie)
                en_dep(eps, rdx, i);
            }
            if(pid == 14){                                                  // protons  ?? ge or eq
                if(E > 0.045 && E < 0.2) eps = 0.5/E;
                else if(E >= 0.2 && E < 100) eps = 2.3*std::pow((E/0.1), 0.2);
                else if(E >= 100) eps = 23.3*std::sqrt(E/1000);
                en_dep(eps, rdx, i);
            }
            if(E>=0.01 && pid < 13 && pid >=8){                             // hadrons
                eps = 1;                                                    // if no interaction
            }
        }
    }
}

void Generator::en_dep(double &eps, double &rdx, int &i){
    double gamma{};
    gamma = distribution(engine);
    if(gamma <= 0.8) ed[i] += eps*1.25;
    double pg{}, p{};
    if((pid == 1 && pie > 10)){                                             // gammas inside 20 m around det.
       if(EL < 0.0215) pg = 300 * std::pow(EL, 4);                          // for gammas 1<E<25 MeV
       else if (EL < 0.1 && EL >= 0.0215) pg = 3e-6*std::pow(EL, -0.8);     // for gammas 25 < E < 100 MeV
       else if (EL < 0.3 && EL >= 0.1) pg = 0.0025*std::pow(EL, 2);         // for gammas 100<E<300 MeV
       else if(EL >= 0.3) pg = 4e-4*std::sqrt(EL);                          // E > 300 MeV
    }
    else if(pid > 7) pg = 0.0125*std::pow(EL, 0.45);                        // for hadrons from GEANT
    p = std::exp(-rdx/0.45)*pg*0.635;                                       // due to pulse selection efficiency
    int N{};
    if(p < 5)  N = poisson(p, N);                                           // Пуассон
    else N = round(0.5+p);
    ien[i] += N;                                                            // No of n (atmospheric n added)
    Nn+=N;
}

void Generator::output(double Ns, double ECRTOT, double THETACR, double PHICR){
    int *n16 = new int[4];
    double *jd = new double[64];
    double *esum = new double[4];

    std::ofstream crds ("C:\\Qt\\progects\\En_dep\\cardS.csv");
    std::ofstream neas ("C:\\Qt\\progects\\En_dep\\NEAS.csv");

    int M{}, M1{}, M2{}, N{};
    double p{}, sEn{};
    for(int j = 0; j < 4; ++j){
        for(int i = 0; i < 16; ++i){
        p = ed[i+16*j]/rp;
        if(p < 10) N = poisson(p, N);
        else N = round(0.5+p);

        ed[i+16*j] = jd[i+16*j] = N;
        if(ed[i+16*j] >= am0){
            n16[j]+=1;
        }
        esum[j] += jd[i+16*j];
        if(ien[i+16*j] > 999) ien[i+16*j] = 999;
        if(jd[i+16*j] > 99999) jd[i+16*j] = 99999;
    }
        if(n16[j] >= krat) ++M1;
        if(esum[j] >= 300) ++M2;
  }
        if(M1 >=1) M = 1;
        if(M2 >= 1) M += 2;
        if(Nn >=10) M += 4;
        if(M1 > 1 || M2 > 1) M += 8;

        if(M < 1){
            Nn = Nm = Nn_sq = M = En = 0;
            for(int i = 0; i < Ndet; ++i) ed[i] = jd[i] = ien[i] = nedmu[i] = 0;
        }
        ++N_event;
        if(En > 2) sEn = std::log10(En);
        else sEn = 1;

        if(Nn > 2500) Nn = 2500;

        std::cout << "Ns = " << Ns << " evs = " << N_event << " M = " << M << " n = " << Nn << " mu= " << Nm << " E= " << ECRTOT/1e3 << ' '
        << "Nn_sq= " << Nn_sq << '\n';

        crds << "x0" << "," << "y0" << "," << "Ns " << "," << "Nn " << "," << "Nm " << "," << "M " << "," << "jd " << "," << "ien "
             << "," << "nedmu " << "," << "sEn " << "," << "ECR " << "," << "THETACR " << "," << "PHICR " << "," << "Nn_sq " << '\n';

        for(int i = 0; i < 64; ++i){
           if(i==0){
           crds << x0 << "," << y0 << "," << Ns  << "," << Nn  << "," << Nm << "," << M << "," << jd[i] << "," << ien[i]
           << "," << nedmu[i] << "," << sEn << "," << ECRTOT/1e3 << "," << THETACR << "," << PHICR << "," << Nn_sq << '\n';
           }
           else{
               crds << ' ' << "," << ' '  << "," << ' '  << "," << ' '  << "," << ' ' << "," << ' ' << "," << jd[i] << "," << ien[i]
               << "," << nedmu[i] << "," << ' ' << "," << ' ' << "," << ' ' << "," << ' ' << "," << ' ' << '\n';
           }
        }

        neas << "N_EAS = " << Ns << ' ' << "evs = " << N_event << ' ' << "am0 = " << am0 << " Nd = " << Ndet << '\n';
}
