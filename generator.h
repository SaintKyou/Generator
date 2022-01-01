#ifndef GENERATOR_H
#define GENERATOR_H

#include <vector>
#include <fstream>
#include <tr1/random>
#include <time.h>
#include <thread>

class Generator
{
private:
//    double *detx = new double[64];
//    double *dety = new double[64];
    double detx[16];
    double dety[16];
    int krat, Number_det;                 // 2, 64
    float par_per_channel, threshold, radius;               //  1.46, 6.0 2.44   !   2.44*3/5=1.46        ! particles per channel of ADC in MSC
    // Corsika
    int pid, pit;
    double pie, pix, piy;
    std::vector<std::vector<int>> times_det;
public:
//    Generator()
//    {
//        radius = 60;
//        krat = 2;
//        Number_det = 64;
//        par_per_channel = 1.46;
//        threshold = 6.0;
//    }
    std::tr1::random_device rd;
    typedef std::tr1::mt19937 MyRng;
    MyRng gen;
    Generator(double _radius, int _krat, int _Number_det, float _par_per_channel, float _threshold, unsigned int _seed): times_det(16, std::vector<int>()), gen(rd())
    {
    gen.seed(static_cast<std::tr1::mt19937::result_type>(_seed));
    //gen.seed(static_cast<std::tr1::mt19937::result_type>(::time(NULL)));
    //gen.seed(static_cast<std::tr1::mt19937::result_type>(std::hash<std::thread::id>()(std::this_thread::get_id())));
    radius = _radius;
    krat = _krat;
    Number_det = _Number_det;
    par_per_channel = _par_per_channel;
    threshold = _threshold;
    std::ifstream conf("C:\\Qt\\progects\\e_dep_clust\\test\\config_rect.txt");
    for(int i = 0; i < 16; ++i) conf >> detx[i] >> dety[i];
/*    for(int i = 0; i < 16; ++i)
    {
        detx[i] = (-4 + 0.5 + i%4)*5;
        detx[i+16] = detx[i] + 20;
        detx[i+32] = detx[i];
        detx[i+48] = detx[i+16];

        dety[i] = (4 - 0.5 - i/4)*5;
        dety[i+16] = dety[i];
        dety[i+32] = dety[i] - 20;
        dety[i+48] = dety[i+32];
    }*/

    Number_muons = Number_hadrons_circle = Number_hadrons_core = Number_neutrons_circle = Number_neutrons = 0;
    E_Gev = Y_shift = X_shift = Number_event = Ne = 0;
    for(int i = 0; i < 16; ++i) en_deposit[i] = muon_det[i] = integ_par[i] = /*k_index[i] = time_trig[i]=*/ 0;
    //fill_coord();
    }
    ~Generator() {}
    double x_core, y_core;
    int Number_muons, Number_hadrons_circle, Number_hadrons_core, Number_neutrons_circle, Number_neutrons /*, M*/;
    double E_Gev, Y_shift, X_shift, Number_event, Ne /*, lg_Ne{}*/;
//    double *en_deposit = new double[64];
//    double *muon_det = new double[64];
//    double *integ_par = new double[64];
    double en_deposit[16];
    double muon_det[16];
    double integ_par[16];

    //double *number_par_adc = new double[64];

//    std::vector <int> t1;
//    int time_trig[16];
//    int time_det1 = 0;
//    int t[16][5];
//    int k_index[16];

    std::vector <double> all;
    void set(int *pid1, double *pie1, double *pix1, double *piy1, int *pit1){
        pid = *pid1;
        pie = *pie1;
        pix = *pix1;
        piy = *piy1;
        pit = *pit1;
    }
    void set_seed(int);
    double uniform(int, int/*, std::tr1::mt19937*/);
    int poisson(double);
    void print_all(/*std::string*/);
    double gamma_dep_inside(double&);
    double gamma_dep_outside(double&);
    double electrons_dep_inside(double&);
    double muons_dep_inside(double&);
    double neutrons_dep_inside(double&);
    double protons_dep_inside(double&);
    void cherenok(double &, double&, int&, int&);

    void randomcp();
    void process();
    void en_dep(double&, double&, int&);
    void output(double Ns, double ECRTOT, double THETACR, double PHICR);
    void print(double Ns, double ECRTOT, double THETACR, double PHICR);
};

#endif // GENERATOR_H

