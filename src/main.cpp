#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include "readData.h"

#if defined(_WIN32) || defined(_WIN64)
    double constante = 1000.0;
#elif defined(__linux__) || defined(__unix__)
    double constante = 1000000.0;
#endif

#define IMAX 50
#define MAXDB 99999999999.99

using namespace std;

double ** matrizAdj; // matriz de adjacencia
int dimension; // quantidade total de vertices

template <typename A, typename B>
void zip(
    const vector<A> &a,
    const vector<B> &b,
    vector<pair<A,B> > &zipped)
{
    for(size_t i=0; i<a.size(); ++i)
    {
        zipped.push_back(make_pair(a[i], b[i]));
    }
}

// Escrever o primeiro e segundo elemento do par do
// vetor zippado em a e b. (Assumindo que
// os vetores tem o mesmo tamanho)
template <typename A, typename B>
void unzip(
    const vector<pair<A, B> > &zipped,
    vector<A> &a,
    vector<B> &b)
{
    for(size_t i=0; i<a.size(); i++)
    {
        a[i] = zipped[i].first;
        b[i] = zipped[i].second;
    }
}

void printData();
void printSolution(vector<int>& sol, double cost);
void construction(vector<int>& sol, double* cost, float alfa);
void rvnd(vector<int>& sol, double* cost);
void Swap(vector<int>& sol, double* cost);
void opt2(vector<int>& sol, double* cost);
void reinsertion(vector<int>& sol, double* cost);
void reinsertion2(vector<int>& sol, double* cost);
void reinsertion3(vector<int>& sol, double* cost);
void dbridge(vector<int>& sol, double* cost);

int main(int argc, char** argv)
{

  readData(argc, argv, &dimension, &matrizAdj);
  //printData();

  float alfa;
  int iterIls = 0;
  int Iils = (dimension >= 150) ? dimension / 2 : dimension;

  time_t end_t;
  time_t begin_t;

  //Mean* sol_mean = createStruct(argv[1]); //saveMean

  double current_cost;
  double better_cost;
  double best_cost = MAXDB;
  vector<int> current_solution;
  vector<int> better_solution;
  vector<int> best_solution;

  srand(time(NULL));

  //for(int j = 0; j < 10; j++) //10 iterations to calculate mean
  //{
    //GILS-RVND
    //begin_t = clock();
    for(int i = 0; i < IMAX; i++)
    {
      iterIls = 0;
      alfa = rand() % 101;
      alfa = alfa/100.0;

      current_solution.clear();
      construction(current_solution, &current_cost, alfa); //solucao inicial

      better_cost = current_cost;
      better_solution = current_solution;

      while(iterIls < Iils)
      {
        rvnd(current_solution, &current_cost); //busca por otimo local

        if(current_cost < better_cost)
        {
          iterIls = 0;
          better_cost = current_cost;
          better_solution = current_solution;
        }
        else
        {
          current_solution = better_solution;
          current_cost = better_cost;
        }

        dbridge(current_solution, &current_cost); //perturbacao

        iterIls++;
      }

      if(better_cost < best_cost)
      {
        best_cost = better_cost;
        best_solution = better_solution;
      }
    }
    //end_t = clock();
    //GILS-RVND

    //addInfo(sol_mean, (end_t - begin_t) / constante, best_cost); //saveMean
    
    //printf("\nTIME: %.2lf\n", (end_t - begin_t) / constante);
    //printSolution(best_solution, best_cost);

    //current_cost = 0;
    //better_cost = 0;
    //best_cost = MAXDB;
    //current_solution.clear();
    //better_solution.clear();
    //best_solution.clear();
  //}

  //calculateMean(sol_mean); //saveMean
  //writeFile(sol_mean, argv[1]); //saveMean]

  return 0;
}

void construction(vector<int>& sol, double* cost, float alfa)
{
  sol.push_back(1);
  sol.push_back(1);

  vector<int> candidates;
  for(int i = 1; i < dimension; i++)
  {  
    candidates.push_back(i+1);
  }

  int index = 0;
  double custo_aux = 0;
  vector<double> deltas;
  vector<int> noInserido;
  vector<int> arestaRemovida;

  while(candidates.size() != 0)
  {
    for(int i = index, j = i+1; i < sol.size()-1; i++, j++)
    {
      custo_aux = matrizAdj[sol[i]][sol[j]];

      for(int k = 0; k < candidates.size(); k++)
      {
        deltas.push_back(matrizAdj[sol[i]][candidates[k]] + matrizAdj[candidates[k]][sol[j]] - custo_aux);
        noInserido.push_back(k);
        arestaRemovida.push_back(i);
      }

      index = i+1;
    }

    vector<pair<double,int> > zipped;
    vector<pair<double,int> > zipped2;
    zip(deltas, noInserido, zipped);
    zip(deltas, arestaRemovida, zipped2);

    sort(zipped.begin(), zipped.end()); //ordenar os nos inseridos pelos deltas
    sort(zipped2.begin(), zipped2.end()); //ordenar as arestas removidas pelos deltas

    unzip(zipped, deltas, noInserido);
    //zipped.clear();
    unzip(zipped2, deltas, arestaRemovida);
    //zipped2.clear();


    int aux = (int)(floor(alfa*candidates.size()));
    int r = 0;

    if(aux != 0)
    {
      r = rand() % aux;
    }
    
    sol.insert(sol.begin() + (arestaRemovida[r]+1), candidates[noInserido[r]]);
    candidates.erase(candidates.begin() + noInserido[r]);

    deltas.clear();
    noInserido.clear();
    arestaRemovida.clear();
  }

  *cost = 0;

  for(int i = 0; i < sol.size()-1; i++)
  {
    *cost = matrizAdj[sol[i]][sol[i+1]] + *cost;
  }
}

void rvnd(vector<int>& sol, double* cost)
{
  int r = 0;
  double cheaper_cost = *cost;

  vector <int> better_sol = sol;
  vector <int> NL;
  for(int i = 0; i < 5; i++)
  {
    NL.push_back(i);
  }

  while(NL.size() != 0)
  {

    r = rand() % NL.size();

    switch(NL[r])
    {
      case 0:
        Swap(better_sol, &cheaper_cost);
        break;

      case 1:
        opt2(better_sol, &cheaper_cost);
        break;

      case 2:
        reinsertion(better_sol, &cheaper_cost);
        break;

      case 3:
        reinsertion2(better_sol, &cheaper_cost);
        break;

      case 4:
        reinsertion3(better_sol, &cheaper_cost);
        break;
    }

    if(cheaper_cost < *cost)
    {
      sol = better_sol;
      *cost = cheaper_cost;

      NL.clear();
      for(int i = 0; i < 5; i++)
      {
        NL.push_back(i);
      }
    }
    else 
    {
      better_sol = sol;
      cheaper_cost = *cost;
      NL.erase(NL.begin() + r);
    } 
  }
}

void Swap(vector<int>& sol, double* cost)
{
  int best_i;
  int best_j;
  double decrease;
  double best_decrease = 0;
  bool improved = false;

  int solSize = sol.size();
  double delta1;
  double delta2;
  
  /*for(int i = 1; i < solSize - 2; i++)
  {

    delta1 = matrizAdj[sol[i-1]][sol[i]];
    delta2 = matrizAdj[sol[i]][sol[i+1]];

    for(int j = i + 1; j < solSize - 1; j++)
    {

      if(j == i + 1)
      {
        decrease = matrizAdj[sol[i-1]][sol[j]] +  matrizAdj[sol[i]][sol[j+1]] -  delta1 -  matrizAdj[sol[j]][sol[j+1]];
      }
      else
      {
        decrease = matrizAdj[sol[i-1]][sol[j]] + matrizAdj[sol[j]][sol[i+1]] + matrizAdj[sol[j-1]][sol[i]] + matrizAdj[sol[i]][sol[j+1]] - delta1 - delta2 - matrizAdj[sol[j-1]][sol[j]] - matrizAdj[sol[j]][sol[j+1]];
      }

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }
  }*/

  for(int i = 1; i < solSize - 2; i++)
  {

    delta1 = matrizAdj[sol[i-1]][sol[i]];
    delta2 = matrizAdj[sol[i]][sol[i+1]];

    //primeiro testa swap adjacente ao no
    int k = i+1;
    decrease = matrizAdj[sol[i-1]][sol[k]] +  matrizAdj[sol[i]][sol[k+1]] -  delta1 -  matrizAdj[sol[k]][sol[k+1]];
    
    if(decrease < best_decrease)
    {
      best_i = i;
      best_j = k;
      improved = true;
      best_decrease = decrease;
    }

    for(int j = i + 2; j < solSize - 1; j++) //depois testa todas as outras possibilidades de swap
    {

      decrease = matrizAdj[sol[i-1]][sol[j]] + matrizAdj[sol[j]][sol[i+1]] + matrizAdj[sol[j-1]][sol[i]] + matrizAdj[sol[i]][sol[j+1]] - delta1 - delta2 - matrizAdj[sol[j-1]][sol[j]] - matrizAdj[sol[j]][sol[j+1]];

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }
  }

  if(improved)
  {
    //printf("best_i = %d\nbest_j = %d\ndecrease of = %lf\n", sol[best_i], sol[best_j], best_decrease);
    *cost = *cost + best_decrease;
    swap(sol[best_i], sol[best_j]);
  }

}

void opt2(vector<int>& sol, double* cost)
{
  int best_i;
  int best_j;
  double decrease;
  double best_decrease = 0;
  bool improved = false;

  int solSize = sol.size();
  double delta1;

  for(int i = 1; i < solSize - 4; i++)
  {
    
    delta1 = matrizAdj[sol[i-1]][sol[i]];

    for(int j = i + 3; j < solSize - 1; j++)
    {
      /*if(i == 1 && j == solSize - 1) //nao inverter o caminho todo
        continue;*/

      decrease = matrizAdj[sol[i-1]][sol[j]] + matrizAdj[sol[i]][sol[j+1]] - delta1 - matrizAdj[sol[j]][sol[j+1]];

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }
  }

  if(improved)
  {
    //printf("best_i = %d\nbest_j = %d\ndecrease of = %lf\n", sol[best_i], sol[best_j], best_decrease);
    *cost = *cost + best_decrease;
    reverse(sol.begin() + best_i, sol.begin() + best_j + 1);
  }
}

void reinsertion(vector<int>& sol, double* cost)
{
  int best_i;
  int best_j;
  double decrease;
  double best_decrease = 0;
  bool improved = false;

  int solSize = sol.size();
  double delta1;
  double delta2;
  double delta3;

  /*for(int i = 1; i < solSize - 1; i++)
  {

    delta1 = matrizAdj[sol[i-1]][sol[i+1]];
    delta2 = matrizAdj[sol[i-1]][sol[i]];
    delta3 = matrizAdj[sol[i]][sol[i+1]];

    for(int j = 1; j < solSize - 1; j++)
    {
      if(j >= (i-1) && j <= (i+1))
        continue;

      if(i < j)
      {
        decrease = delta1 + matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i]][sol[j+1]] - delta2 - delta3 - matrizAdj[sol[j]][sol[j+1]];
      }

      if(i > j)
      {
        decrease = matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i]][sol[j+1]] + delta1 - matrizAdj[sol[j]][sol[j+1]] - delta2 - delta3;
      }

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }
  }*/

  for(int i = 1; i < solSize - 1; i++)
  {

    delta1 = matrizAdj[sol[i-1]][sol[i+1]];
    delta2 = matrizAdj[sol[i-1]][sol[i]];
    delta3 = matrizAdj[sol[i]][sol[i+1]];

    for(int j = i+2; j < solSize - 1; j++) //primeiro todas as possibilidades com j > i
    {

      decrease = delta1 + matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i]][sol[j+1]] - delta2 - delta3 - matrizAdj[sol[j]][sol[j+1]];

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }

    for(int j = 1; j <= i-2; j++) //depois todas as possibilidades com i > j
    {

      decrease = matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i]][sol[j+1]] + delta1 - matrizAdj[sol[j]][sol[j+1]] - delta2 - delta3;

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }
  }

  if(improved)
  {
    //printf("best_i = %d\nbest_j = %d\ndecrease of = %lf\n", sol[best_i], sol[best_j], best_decrease);
    *cost = *cost + best_decrease;

    int node = sol[best_i];

    if(best_i > best_j)
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 1);
      sol.insert(sol.begin() + best_j + 1, node);
    }
    else
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 1);
      sol.insert(sol.begin() + best_j, node);
    }
  }

}

void reinsertion2(vector<int>& sol, double* cost)
{
  int best_i;
  int best_j;
  double decrease;
  double best_decrease = 0;
  bool improved = false;

  int solSize = sol.size();
  double delta1;
  double delta2;
  double delta3;

  /*for(int i = 1; i < solSize - 2; i++)
  {

    delta1 = matrizAdj[sol[i-1]][sol[i+2]];
    delta2 = matrizAdj[sol[i-1]][sol[i]];
    delta3 = matrizAdj[sol[i+1]][sol[i+2]];

    for(int j = 1; j < solSize - 1; j++)
    {
      if(j >= (i-1) && j <= (i+1))
        continue;

      if(i < j)
      {
        decrease = delta1 + matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i+1]][sol[j+1]] - delta2 - delta3 - matrizAdj[sol[j]][sol[j+1]];
      }

      if(i > j)
      {
        decrease = matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i+1]][sol[j+1]] + delta1 - matrizAdj[sol[j]][sol[j+1]] - delta2 - delta3;
      }

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }
  }*/

  for(int i = 1; i < solSize - 2; i++)
  {

    delta1 = matrizAdj[sol[i-1]][sol[i+2]];
    delta2 = matrizAdj[sol[i-1]][sol[i]];
    delta3 = matrizAdj[sol[i+1]][sol[i+2]];

    for(int j = i+2; j < solSize - 1; j++) //primeiro todas as possibilidades com j > i
    {

      decrease = delta1 + matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i+1]][sol[j+1]] - delta2 - delta3 - matrizAdj[sol[j]][sol[j+1]];

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }

    for(int j = 1; j <= i-2; j++) //depois todas as possibilidades com i > j
    {

      decrease = matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i+1]][sol[j+1]] + delta1 - matrizAdj[sol[j]][sol[j+1]] - delta2 - delta3;

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }
  }

  if(improved)
  {
    //printf("best_i = %d\nbest_j = %d\ndecrease of = %lf\n", sol[best_i], sol[best_j], best_decrease);
    *cost = *cost + best_decrease;

    vector<int> subsequence(sol.begin() + best_i, sol.begin() + best_i + 2);

    if(best_i > best_j)
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 2);
      sol.insert(sol.begin() + best_j + 1, subsequence.begin(), subsequence.end());
    }
    else
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 2);
      sol.insert(sol.begin() + best_j - 1, subsequence.begin(), subsequence.end());
    }
  }

}

void reinsertion3(vector<int>& sol, double* cost)
{
  int best_i;
  int best_j;
  double decrease;
  double best_decrease = 0;
  bool improved = false;

  int solSize = sol.size();
  double delta1;
  double delta2;
  double delta3;

  /*for(int i = 1; i < solSize - 3; i++)
  {

    delta1 = matrizAdj[sol[i-1]][sol[i+3]];
    delta2 = matrizAdj[sol[i-1]][sol[i]];
    delta3 = matrizAdj[sol[i+2]][sol[i+3]];

    for(int j = 1; j < solSize - 1; j++)
    {
      if(j >= (i-1) && j <= (i+2))
        continue;

      if(i < j)
      {
        decrease = delta1 + matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i+2]][sol[j+1]] - delta2 - delta3 - matrizAdj[sol[j]][sol[j+1]];
      }

      if(i > j)
      {
        decrease = matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i+2]][sol[j+1]] + delta1 - matrizAdj[sol[j]][sol[j+1]] - delta2 - delta3;
      }

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }
  }*/

  for(int i = 1; i < solSize - 3; i++)
  {

    delta1 = matrizAdj[sol[i-1]][sol[i+3]];
    delta2 = matrizAdj[sol[i-1]][sol[i]];
    delta3 = matrizAdj[sol[i+2]][sol[i+3]];

    for(int j = i+3; j < solSize - 1; j++) //primeiro todas as possibilidades com j > i
    {

      decrease = delta1 + matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i+2]][sol[j+1]] - delta2 - delta3 - matrizAdj[sol[j]][sol[j+1]];

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }

    for(int j = 1; j <= i-2; j++) //depois todas as possibilidades com i > j
    {

      decrease = matrizAdj[sol[j]][sol[i]] + matrizAdj[sol[i+2]][sol[j+1]] + delta1 - matrizAdj[sol[j]][sol[j+1]] - delta2 - delta3;

      if(decrease < best_decrease)
      {
        best_i = i;
        best_j = j;
        improved = true;
        best_decrease = decrease;
      }
    }
  }

  if(improved)
  {
    //printf("best_i = %d\nbest_j = %d\ndecrease of = %lf\n", sol[best_i], sol[best_j], best_decrease);
    *cost = *cost + best_decrease;

    vector<int> subsequence(sol.begin() + best_i, sol.begin() + best_i + 3);

    if(best_i > best_j)
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 3);
      sol.insert(sol.begin() + best_j + 1, subsequence.begin(), subsequence.end());
    }
    else
    {
      sol.erase(sol.begin() + best_i, sol.begin() + best_i + 3);
      sol.insert(sol.begin() + best_j - 2, subsequence.begin(), subsequence.end());
    }
  }

}

void dbridge(vector<int>& sol, double* cost)
{
  int tam1 = 0; //tam subsequencia 1
  int tam2 = 0; //tam subsequencia 2
  int index1 = 0; //indice de onde comeca subsequencia 1
  int index2 = 0; //indice de onde comeca subsequencia 2

  //tam das subsequencias entre 2 e |V|/10
  if(dimension/10 <= 2)
  {
    tam1 = 2;
    tam2 = 2;
  }
  else
  {
    tam1 = rand() % ((dimension/10)-1) + 2;
    tam2 = rand() % ((dimension/10)-1) + 2;
  }
  
  //indice de inicio da primeira subsequencia dentro de um intervalo valido
  index1 = rand() % ((sol.size()-2)-(tam1-1)) + 1; //rand() % (sol.size()-2) + 1 da um indice dentro vetor excluindo origem e destino

  //procura indice de inicio da segunda subsequencia dentro de um intervalo valido
  while(1)
  {
    index2 = rand() % ((sol.size()-2)-(tam2-1)) + 1; //rand() % (sol.size()-2) + 1 da um indice dentro vetor excluindo origem e destino

    if((index1 < index2) && (index1+tam1 > index2)) //subsequencias com intersecao
      continue;

    if((index2 < index1) && (index2+tam2 > index1)) //subsequencias com intersecao
      continue;

    if(index1 == index2)
      continue;
    
    break;
  }

  if(index2 < index1)
  {
    int aux_index = index1;
    index1 = index2;
    index2 = aux_index;

    int aux_tam = tam1;
    tam1 = tam2;
    tam2 = aux_tam;
  }

  if(index1 + tam1 == index2)
  {
    *cost = *cost + (matrizAdj[sol[index1-1]][sol[index2]] + matrizAdj[sol[index2+tam2-1]][sol[index1]] + matrizAdj[sol[index1+tam1-1]][sol[index2+tam2]] - matrizAdj[sol[index1-1]][sol[index1]] - matrizAdj[sol[index1+tam1-1]][sol[index1+tam1]] - matrizAdj[sol[index2+tam2-1]][sol[index2+tam2]]);
  }
  else
  {
    *cost = *cost + (matrizAdj[sol[index1-1]][sol[index2]] + matrizAdj[sol[index2+tam2-1]][sol[index1+tam1]] + matrizAdj[sol[index2-1]][sol[index1]] + matrizAdj[sol[index1+tam1-1]][sol[index2+tam2]] - matrizAdj[sol[index1-1]][sol[index1]] - matrizAdj[sol[index1+tam1-1]][sol[index1+tam1]] - matrizAdj[sol[index2-1]][sol[index2]] - matrizAdj[sol[index2+tam2-1]][sol[index2+tam2]]);
  }

  vector<int> subsequence1(sol.begin() + index1, sol.begin() + index1 + tam1);
  vector<int> subsequence2(sol.begin() + index2, sol.begin() + index2 + tam2);

  //printf("best_i = %d\ntam_i = %d\nbest_j = %d\ntam_j = %d\nnew cost = %lf\n", sol[index1], tam1, sol[index2], tam2, *cost);
  
  sol.erase(sol.begin() + index1, sol.begin() + index1 + tam1); //remove a primeira subsequencia
  sol.insert(sol.begin() + index1, subsequence2.begin(), subsequence2.end()); //insere a segunda no lugar

  index2 += tam2 - tam1; //indice da segunda subsequencia deslocado

  sol.erase(sol.begin() + index2, sol.begin() + index2 + tam2); //remove a segunda subsequencia
  sol.insert(sol.begin() + index2, subsequence1.begin(), subsequence1.end()); //insere a primeira no lugar
}

void printSolution(vector<int>& sol, double cost)
{
  for(int i = 0; i < sol.size(); i++)
  {
    cout << sol[i] << " ";
  }
  cout << "cost = " << cost << endl;

  cost = 0;
  for(int i = 0; i < sol.size()-1; i++)
  {
    cost += matrizAdj[sol[i]][sol[i+1]];
  }
  cout << "calculated cost = " << cost << "\n" << endl;
}

void printData()
{
  cout << "dimension: " << dimension << endl;
  for (int i = 1; i <= dimension; i++)
  {
    for (int j = 1; j <= dimension; j++)
    {
      cout << matrizAdj[i][j] << " ";
    }
    cout << endl;
  }
}