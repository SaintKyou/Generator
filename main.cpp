#include <iostream>
#include <future>
#include <thread>
#include <random>
#include <ctime>
#include <cmath>
#include "generator.h"

using namespace std;

void loop(unsigned int start, unsigned int end)
{
    for(unsigned int k = start; k < end; ++k){
//        for(unsigned int i = 0; i < AllEvent.size(); ++i){
//            ++ShowerSimCounter;
            //randomcp_();
//            for(unsigned int j = 0; j<AllEvent[i].size(); ++j){
              //  process_(&AllEvent[i][j].id,&AllEvent[i][j].p,&AllEvent[i][j].x,&AllEvent[i][j].y);
            }
            //output_(&ShowerSimCounter,&ShowerEnergy[i], &ShowerTheta[i], &ShowerPhi[i]);
        }
//   }
//}

int main ()
{
  std::ifstream f1("C:\\Users\\ASUS Zephyrus\\Desktop\\test\\cors_out.txt");
  constexpr int N = 897117;
  int* id = new int[N];
  double* p = new double[N];
  double* x = new double[N];
  double* y = new double[N];
  Generator par(60, 2, 64, 1.46, 6.0);

  //x0,y0,Ns,Nn,Mm,M,(jd(i),ien(i),nedmu(i),i=1,64),& sEn,ECR,THETACR,PHICR,Mh0

  int throws{4255}, ShowerCounter{};
//  std::cout << "Бросков: " << '\n';
//  std::cin >> throws;
//  std::cout << "Бросков " << throws << std::endl;
  for(int j = 0; j < throws; ++j){
    ++ShowerCounter;
    par.randomcp();
    for(int i = 0; i < N; ++i){
        f1 >> id[i] >> p[i] >> x[i] >> y[i];
        par.set(&id[i], &p[i], &x[i], &y[i]);
        par.process();
    }
    par.output(ShowerCounter, 500000, 0, 0);
  }

  //std::cout << par.x0 << ' ' << par.y0 << endl;
  //for(int i = 0; i < 64; ++i) std::cout << par.ed[i] << ',';
  //std::cout << '\n';

  //std::cout << par.Nn << ' ' << par.Mm << ' ' << par.ed[0] << endl;;
  //par.output(1, 500000, 0, 0);


  return 0;
}
