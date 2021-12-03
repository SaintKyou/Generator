#include <iostream>
#include <thread>
#include <ctime>
#include <cmath>
#include <mutex>
#include "generator.h"

int ShowerCounter{};
constexpr int N = 897117;
int* id = new int[N];
double* p = new double[N];
double* x = new double[N];
double* y = new double[N];
std::ifstream f2("C:\\Users\\ASUS Zephyrus\\Desktop\\test\\cors_out.txt");

std::mutex g_lock;

void loop3(int number)
{
    Generator par2(60, 1, 16, 1.46, 6.0);
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
    g_lock.lock();
    par2.print_all();
    g_lock.unlock();
}

std::vector<std::thread> thread_pool;

int main ()
{
  for(int i = 0; i < N; ++i) f2 >> id[i] >> p[i] >> x[i] >> y[i];
  int throws{1000};
  int number_of_threads = 8;
  int step = throws/number_of_threads; // std::thread::hardware_concurrency()

  for(int i = 0; i < number_of_threads; ++i){
      if(i==number_of_threads-1) step = throws - (number_of_threads-1)*step;
      thread_pool.emplace_back(loop3, step);

  }

  for(int i = 0; i < number_of_threads; ++i){
      std::thread& entry = thread_pool[i];
      entry.join();
  }

  return 0;
}
