#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <unordered_set>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <map>
#include <stack>
#include <deque>
#include <Eigen/Dense>
#include <random>
#include <numeric>
#include <memory>
#include <cmath>
#include <chrono>
#include <random>
#include <array>
using namespace Eigen;

class Data
{
public:
  long long sampled = 0;
  long long sieved = 0;
  long long compared = 0;
  std::map<int, int> L_norms;
  std::map<int, int> S_norms;
  std::chrono::milliseconds timeTaken;
};

/**
 * @brief sample independently samples each coordinate of a vector of length n
 * from a discrete gaussian with parameter s
 *
 * @param n is the dimension of the vector outputted
 * @param s is the discrete gaussian parameter
 * @return a vector of length n where each coordinate is sampled independently
 * from the discrete gaussian with parameter s
 */
Matrix<long, Dynamic, 1> sample(int n, int s)
{
  std::random_device rd;
  std::mt19937_64 mt(rd());

  std::uniform_int_distribution<int> disc_gaussian_range_dist(-10 * s, 10 * s);
  Matrix<long, Dynamic, 1> v(n);

  for (int j = 0; j < n; j++)
  {
    bool accept = false;
    int x;
    double p;
    while (!accept)
    {
      x = disc_gaussian_range_dist(mt);
      p = exp(-M_PI * pow(x, 2) / pow(s, 2));
      std::bernoulli_distribution acceptance_coin(p);
      accept = acceptance_coin(mt);
    }
    v(j) = x;
  }
  return v;
}

/**
 * @brief gauss_reduce finds a collision in the integer lattice by sieving a
 * lattice vector v with the vectors in the pairwise reduced set of vectors L.
 * It then adds short combinations of vectors in L with v to S to resieve and
 * accordingly deletes their longer counterparts from L.
 *
 * @param v is the vector being sieved
 * @param L is the set of pairwise reduced vectors we use to sieve v
 * @param S is the set of vectors we add to if the reduced v has short
 * combinations with vectors in L
 * @return v after being sieved by vectors in L
 */
Matrix<long, Dynamic, 1> gauss_reduce(Matrix<long, Dynamic, 1> v, std::vector<Matrix<long, Dynamic, 1>> *L, std::vector<Matrix<long, Dynamic, 1>> *S, Data *output)
{
  bool updated_v = false;
  long norm_v = v.dot(v);
  do
  {
    updated_v = false;
    for (auto it_w = L->begin(); it_w != L->end();)
    {
      output->compared++;
      long norm_w = it_w->dot(*it_w);
      Matrix<long, Dynamic, 1> diff = v - *it_w;
      long diff_norm = diff.dot(diff);
      if (norm_w <= norm_v && diff_norm < norm_v)
      {
        float proj_v_on_w = ((float)(v.dot(*it_w))) / norm_w;
        int c = (int)round(proj_v_on_w);
        v = v - c * (*it_w);
        norm_v = v.dot(v);
        updated_v = true;
      }
      it_w++;
    }
  } while (updated_v == true);

  for (auto it_w = L->begin(); it_w != L->end();)
  {
    output->compared++;
    long norm_w = it_w->dot(*it_w);
    Matrix<long, Dynamic, 1> diff = *it_w - v;
    long diff_norm = diff.dot(diff);
    if (norm_w > norm_v && diff_norm < norm_w)
    {
      it_w = L->erase(it_w);
      S->push_back(diff);
    }
    else
    {
      it_w++;
    }
  }

  return v;
}

/**
 * @brief gauss_sieve feeds gauss_reduce with vectors either sampled from the
 * discrete gaussian or taken from Sin order to find a collision in the integer
 * lattice. In conjunction, it collects data on the number of vectors sampled,
 * number of vectors sieved. Finally, at a collision it records the distribution
 * of the squares of the norms of vectors in L and S.
 *
 * @param n is the dimension of the discrete gaussian we sample from
 * @param s is the discrete gaussian parameter
 * @param max_collisions is the number of collisions before we terminate our sieving
 * @return a tuple containing number of vectors sampled, number of vectors
 * sieved, the distribution of the squares of the norms vectors in L, and the
 * distribution of the squares of the norms of vectors in S.
 */
Data gauss_sieve(int n, int s, int max_collisions)
{
  Data output;
  std::vector<Matrix<long, Dynamic, 1>> L;
  std::vector<Matrix<long, Dynamic, 1>> S;

  int collisions = 0;
  std::chrono::milliseconds startTimeMS = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::system_clock::now().time_since_epoch());
  while (collisions < max_collisions)
  {
    Matrix<long, Dynamic, 1> v(n);
    output.sieved++;
    if (S.empty())
    {
      output.sampled++;
      v = sample(n, s);
    }
    else
    {
      v = S.back();
      S.pop_back();
    }

    v = gauss_reduce(v, &L, &S, &output);
    L.push_back(v);
    if (v.dot(v) == 0)
    {
      collisions++;
    }
  }

  std::map<int, int> Lnorms;
  for (int i = 0; i < L.size(); i++)
  {
    int curr = L[i].dot(L[i]);
    auto search = Lnorms.find(curr);
    if (search != Lnorms.end())
    {
      (search->second)++;
    }
    else
    {
      Lnorms.insert({curr, 1});
    }
  }
  std::map<int, int> Snorms;
  for (int i = 0; i < S.size(); i++)
  {
    int curr = S[i].dot(S[i]);
    auto search = Snorms.find(curr);
    if (search != Snorms.end())
    {
      (search->second)++;
    }
    else
    {
      Snorms.insert({curr, 1});
    }
  }
  output.L_norms = Lnorms;
  output.S_norms = Snorms;
  std::chrono::milliseconds finishTimeMS = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::system_clock::now().time_since_epoch());

  output.timeTaken = finishTimeMS - startTimeMS;
  return output;
}

int main(int argc, char **argv)
{
  if (argc != 2)
  {
    std::cout << "Bad Argument" << std::endl;
    return -1;
  }
  std::ofstream file;
  std::string filename = (argv[1]);
  file.open(filename);
  std::array<int, 3> s = {10, 100, 1000};
  std::array<int, 3> d = {16, 32, 64};
  for (auto it_d = d.begin(); it_d != d.end(); it_d++)
  {
    for (auto it_s = s.begin(); it_s != s.end(); it_s++)
    {
      file << "Parameters: dimension = " << *it_d << " s = " << *it_s << "\n";
      Data data = gauss_sieve(*it_d, *it_s, 1);
      file << "Time Taken:" << std::to_string((double)data.timeTaken.count() / 1000.0) << "\n";
      file << "Number of Sampled Vectors:" << data.sampled << "\n";
      file << "Number of Sieved Vectors:" << data.sieved << "\n";
      file << "Number of Compared Vectors:" << data.compared << "\n";
      file << "Norm Distribution in L:"
           << "\n";
      for (const auto &x : data.L_norms)
        file << std::setw(10) << x.first << " | " << x.second << "\n";
      file << "Norm Distribution in S:"
           << "\n";
      for (const auto &x : data.S_norms)
        file << std::setw(10) << x.first << " | " << x.second << "\n";
      file << "\n";
      file.flush();
    }
  }

  file.close();
  return 0;
}