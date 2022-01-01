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
int* id = new int[N];
double* p = new double[N];
double* x = new double[N];
double* y = new double[N];
int* t = new int[N];
//std::ifstream f2("C:\\Users\\ASUS Zephyrus\\Desktop\\test\\cors_out.txt");
std::ifstream f2("C:\\Users\\ASUS Zephyrus\\Desktop\\test\\str.txt");

std::vector<double> all;
#include <mutex>
std::mutex g_lock;

//void loop2(Generator par2, int number)
//{
//    for(unsigned int k = 0; k < number; ++k){
//        //++ShowerCounter;
//        par2.randomcp();
//        for(int i = 0; i < N; ++i){
//            //f1 >> id[i] >> p[i] >> x[i] >> y[i];
//            par2.set(&id[i], &p[i], &x[i], &y[i]);
//            par2.process();
//            }
//        par2.output(k, 500000, 0, 0);
//        }
//    par2.print_all();
//}

void loop3(int number, unsigned int seed/*, std::string f*/)
{
    Generator par2(20, 2, 16, 1.46, 6.0, seed);
    //par2.set_seed(0);
    for(unsigned int k = 0; k < number; ++k){
        //++ShowerCounter;
        par2.randomcp();
        for(int i = 0; i < N; ++i){
            //f1 >> id[i] >> p[i] >> x[i] >> y[i];
            par2.set(&id[i], &p[i], &x[i], &y[i], &t[i]);
            par2.process();
            }
        par2.output(k, 500000, 0, 0);
        }
    g_lock.lock();
    par2.print_all(/*f*/);
    g_lock.unlock();
}

std::vector<std::thread> thread_pool;
//std::list<std::future<void>> taskList;

int main ()
{
  for(int i = 0; i < N; ++i) f2 >> id[i] >> p[i] >> x[i] >> y[i] >> t[i];
//  for(int i = 0; i < 1000; ++i) xy >> xx[i] >> yy[i];
//  Generator par[16];
  int throws{5000};
  int number_of_threads = 8;
  int step = throws/number_of_threads; // std::thread::hardware_concurrency()

//  srand((unsigned int)time(0));

  for(int i = 0; i < number_of_threads; ++i){
      if(i==number_of_threads-1) step = throws - (number_of_threads-1)*step;
      thread_pool.emplace_back(loop3, step, i/*, "cardS.txt"*/);

  }

  for(int i = 0; i < number_of_threads; ++i){
      std::thread& entry = thread_pool[i];
      entry.join();
  }

//  std::cout << ::time(NULL) << '\n';

  return 0;
}
