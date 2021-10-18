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

constexpr int N = 897117;
int* id = new int[N];
double* p = new double[N];
double* x = new double[N];
double* y = new double[N];
std::ifstream f1("C:\\Users\\ASUS Zephyrus\\Desktop\\test\\cors_out.txt");

void loop3(int number)
{
    Generator par2(60, 2, 64, 1.46, 6.0);
    for(unsigned int k = 0; k < number; ++k){
        //++ShowerCounter;
        par2.randomcp();
        for(int i = 0; i < N; ++i){
            //f1 >> id[i] >> p[i] >> x[i] >> y[i];
            par2.set(&id[i], &p[i], &x[i], &y[i]);
            par2.process();
            }
        par2.output(k, 500000, 0, 0);
        }
    par2.print_all();
}

std::vector<std::thread> thread_pool;
//std::list<std::future<void>> taskList;

int main ()
{
  for(int i = 0; i < N; ++i) f1 >> id[i] >> p[i] >> x[i] >> y[i];
  int throws{5000};
  int number_of_threads = 16;
  int step = throws/number_of_threads; // std::thread::hardware_concurrency()

  for(int i = 0; i < number_of_threads; ++i){
      if(i==number_of_threads-1) step = throws - (number_of_threads-1)*step;
      thread_pool.emplace_back(loop3, step);

  }

  for(auto& entry: thread_pool)
      entry.join();

//loop3(throws);    // 1 thread

  return 0;
}
