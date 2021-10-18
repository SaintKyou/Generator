#ifndef GENERATOR_H
#define GENERATOR_H

#include <vector>
#include <fstream>

class Generator
{
private:
    double *detx = new double[64];
    double *dety = new double[64];
    void fill_coord(){                                  // заполнение координат детекторов
        for(int i = 0; i < 16; ++i)
        {
            detx[i] = (-4 + 0.5 + i%4)*5;
            detx[i+16] = detx[i] + 20;
            detx[i+32] = detx[i];
            detx[i+48] = detx[i+16];

            dety[i] = (4 - 0.5 - i/4)*5;
            dety[i+16] = dety[i];
            dety[i+32] = dety[i] - 20;
            dety[i+48] = dety[i+32];
        }
    }
    int radius{}, krat{}, Number_det{};                 //60, 2, 64
    float par_per_channel{}, threshold{};               //  1.46, 6.0 2.44   !   2.44*3/5=1.46        ! particles per channel of ADC in MSC
    // Corsika
    int pid{};
    double pie{}, pix{}, piy{};
public:
    Generator()
    {
        radius = 60;
        krat = 2;
        Number_det = 64;
        par_per_channel = 1.46;
        threshold = 6.0;
    }
    Generator(int _radius, int _krat, int _Number_det, float _par_per_channel, float _threshold)
    {
    radius = _radius;
    krat = _krat;
    Number_det = _Number_det;
    par_per_channel = _par_per_channel;
    threshold = _threshold;
    fill_coord();
    }
    ~Generator() {}
    double x_core{}, y_core{};
    int Number_muons{}, Number_hadrons_circle{}, Number_hadrons_core{}, Number_neutrons_circle{}, Number_neutrons{};
    double E_Gev{}, Y_shift{}, X_shift{}, Number_event{}, Ne{};
    double *en_deposit = new double[64];
    double *muon_det = new double[64];
    double *integ_par = new double[64];

    //double *number_par_adc = new double[64];

    std::vector <double> all;
    void set(int *pid1, double *pie1, double *pix1, double *piy1){
        pid = *pid1;
        pie = *pie1;
        pix = *pix1;
        piy = *piy1;
    }
    void print_all();
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
