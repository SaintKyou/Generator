#include <iostream>
#include <future>
#include <thread>
#include <random>
#include <ctime>
#include <cmath>
#include "generator.h"

#include <future>
#include <memory>
#include <list>

int ShowerCounter{};
constexpr int N = 897117;
int throws{100};
Generator par(60, 2, 64, 1.46, 6.0);
int* id = new int[N];
double* p = new double[N];
double* x = new double[N];
double* y = new double[N];
std::ifstream f1("C:\\Users\\ASUS Zephyrus\\Desktop\\test\\cors_out.txt");

std::vector<double> all;

void loop(unsigned int start, unsigned int end)
{
    for(unsigned int k = start; k < end; ++k){
        ++ShowerCounter;
        par.randomcp();
        for(int i = 0; i < N; ++i){
            f1 >> id[i] >> p[i] >> x[i] >> y[i];
            par.set(&id[i], &p[i], &x[i], &y[i]);
            par.process();
        }
        par.output(ShowerCounter, 500000, 0, 0);
        if(k == end-1) par.print_all(); /*all.insert(std::end(all), std::begin(par.all), std::end(par.all));*/  //std::cout << par.all.size() << '\n';
        }
}

void loop2(Generator par2)
{
    for(unsigned int k = 0; k < 500; ++k){
        ++ShowerCounter;
        //std::cout << ShowerCounter << '\n';
        par2.randomcp();
        for(int i = 0; i < N; ++i){
            //f1 >> id[i] >> p[i] >> x[i] >> y[i];
            par2.set(&id[i], &p[i], &x[i], &y[i]);
            par2.process();
        }
        par2.output(k, 500000, 0, 0);
        if(k == 500-1) par2.print_all(); //std::cout << par.all.size() << '\n';
        }
}

//std::ofstream crdss ("C:\\Qt\\progects\\En_dep\\cardS.txt");

//void print()
//{
//    for(int i = 0; i < all.size(); ++i){
//        crdss << all[i] << ",";
//        if((i+1)%203==0){
//            crdss << '\n';
//        }
//    }
//}


std::vector<std::thread> thread_pool;
std::list<std::future<void>> taskList;

int main ()
{
  for(int i = 0; i < N; ++i) f1 >> id[i] >> p[i] >> x[i] >> y[i];
  int step = throws/10; // std::thread::hardware_concurrency()


    std::thread t1(loop2, par);
    std::thread t2(loop2, par);
    std::thread t3(loop2, par);
    std::thread t4(loop2, par);
    std::thread t5(loop2, par);
    std::thread t6(loop2, par);
    std::thread t7(loop2, par);
    std::thread t8(loop2, par);
    std::thread t9(loop2, par);
    std::thread t10(loop2, par);
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();
    t7.join();
    t8.join();
    t9.join();
    t10.join();

   // print();

//    for(int i = 0; i < 10; ++i){
//        thread_pool.emplace_back(loop2, par, i);
//    }

//    for(auto& entry: thread_pool)
//        entry.join();



//e    loop(0, throws);

    //std::cout << par.all.size() << '\n'; //~5000*(6+ 64*3)   9*36*3


//  int count{};
//  std::vector<std::future<void>> futures;
//  for(auto i = 0; i < 10; ++i){
//  futures.push_back(std::async(loop, count, count*(i+1)*step));
//  count+=step;
//  }



//  std::thread t1(loop, 0, 5000);
//  std::thread t2(loop, 0, 500);
//  t1.join();
//  t2.join();

//      int count{}, coun2{step};
//      std::future<void> moveThis;

//      for(auto i = 0; i < 10; ++i){
//      taskList.push_back(std::async(loop, count, coun2));
//      std::cout << count << ' ' << coun2 << '\n';
//      count +=step;
//      coun2+=step;
//      }

  //x0,y0,Ns,Nn,Mm,M,(jd(i),ien(i),nedmu(i),i=1,64),& sEn,ECR,THETACR,PHICR,Mh0


//  for(int j = 0; j < throws; ++j){
//    ++ShowerCounter;
//    par.randomcp();
//    for(int i = 0; i < N; ++i){
//        f1 >> id[i] >> p[i] >> x[i] >> y[i];
//        par.set(&id[i], &p[i], &x[i], &y[i]);
//        par.process();
//    }
//    par.output(ShowerCounter, 500000, 0, 0);
//  }


//    for(unsigned i = 0; i < 10; ++i){
//        thread_pool.emplace_back(loop, count, coun2);
//        //std::cout << coun2 << ' ' << count*16 << '\n';
//        //++count;
//        //coun2 = (count-1)*16;
//        count+=step;
//        coun2+=step;
//    }

//    for(auto& entry: thread_pool)
//        entry.join();



  return 0;
}
