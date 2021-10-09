#ifndef GENERATOR_H
#define GENERATOR_H

#include <vector>
#include <fstream>

class Generator
{
private:
    double *detx = new double[64];
    double *dety = new double[64];
    void fill_coord(){                                      // заполнение координат детекторов
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
    int r0{}, krat{}, Ndet{};         //60, 2, 64
    float rp{}, am0{};               //  1.46, 6.0 rp 2.44   !   2.44*3/5=1.46        ! particles per channel of ADC in MSC
    // Corsika
    int pid{};
    double pie{}, pix{}, piy{};
public:
    Generator()
    {
        r0 = 60;
        krat = 2;
        Ndet = 64;
        rp = 1.46;
        am0 = 6.0;
    }
    Generator(int _r0, int _krat, int _Ndet, float _rp, float _am0)
    {
    r0 = _r0;
    krat = _krat;
    Ndet = _Ndet;
    rp = _rp;
    am0 = _am0;
    fill_coord();
    }
    ~Generator() {}
    double x0{}, y0{};
    int Nm{}, Nh_sq{}, Nh_core{}, Nn_sq{}, Nn{}, mask{}; // mask - M1 + M2 + 1 bit (Nn>10)
    double E{}, Y{}, X{}, N_event{}, En{};
    double *ed = new double[64];
    double *nedmu = new double[64];
    double *ien = new double[64];
    std::vector <double> all;
    void set(int *pid1, double *pie1, double *pix1, double *piy1){
        pid = *pid1;
        pie = *pie1;
        pix = *pix1;
        piy = *piy1;
    }
    void print_all();
    int poisson(double &, int&);                        // Пуассон

    void randomcp();
    void process();
    void en_dep(double&, double&, int&);
    void output(double Ns, double ECRTOT, double THETACR, double PHICR);
};

#endif // GENERATOR_H
