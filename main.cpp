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

int main (int argc, char **argv)
{
  std::ifstream f1("C:\\Users\\ASUS Zephyrus\\Desktop\\test\\2\\cors_out.txt");
  constexpr int N = 897117;
  int* id = new int[N];
  double* p = new double[N];
  double* x = new double[N];
  double* y = new double[N];
  Generator par(60, 2, 64, 1.46, 6.0);

  //x0,y0,Ns,Nn,Mm,M,(jd(i),ien(i),nedmu(i),i=1,64),& sEn,ECR,THETACR,PHICR,Mh0
  par.randomcp();
  for(int i = 0; i < N; ++i){
      f1 >> id[i] >> p[i] >> x[i] >> y[i];

      par.set(&id[i], &p[i], &x[i], &y[i]);
      par.process();
      //par.output()

  }
  std::cout << par.x0 << ' ' << par.y0 << endl;
  //for(int i = 0; i < 64; ++i) std::cout << par.ed[i] << ',';
  //std::cout << '\n';

  //std::cout << par.Nn << ' ' << par.Mm << ' ' << par.ed[0] << endl;;
  par.output(1, 500000, 0, 0);



  //std::cout << id[1] << ' ' << p[1] << ' ' << x[1] << ' ' << y[1] << endl;


  return 0;
}
