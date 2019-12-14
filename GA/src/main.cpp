#include "readData.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>    


using namespace std;
using namespace std::chrono; 

double ** matrizAdj;
int dimension; // number of vertices in the problem

struct solution{
  vector<int> cities;
  double distance; //delta after insert the node k in the edge {i,j}
  bool operator<(const solution& a) const
  {
      return distance < a.distance;
  }
};

/*This function prints all the elements of a list*/
void printList(vector<int> &list){
  for (auto i = list.begin(); i != list.end(); ++i) 
    cout << *i << " "; 
}
/*This function calculates the distance given a solution of TSP problem.*/
double getDistance(vector<int> &listaCidade){
  double distance = 0;

  for(int i = 0; i < listaCidade.size()- 1; i++){
      distance += matrizAdj[listaCidade.at(i)][listaCidade.at(i+1)];
  }
  return distance; 
}
vector<solution> randomConstruction(int populationSize){

  vector<solution> population;
  solution chromosome;
  vector<int> candidateNodes; //List of candidate nodes


  // Creating a list with the candidate nodes starting from the city number 2.
  for (int i = 1; i <= dimension; i++) candidateNodes.push_back(i); candidateNodes.push_back(1);


  for(int i = 0; i < populationSize;i++){
    random_shuffle(candidateNodes.begin()+1, candidateNodes.end()-1);
    chromosome.cities = candidateNodes;
    chromosome.distance = getDistance(chromosome.cities);
    population.push_back(chromosome);
  }

  return population;
}

void PMX(vector<vector<int> > &parents, vector<solution> &childrens){
  vector<int> subsequence;
  std::vector<int>::iterator it;
  solution kid;
  int index, n, indexTest;
  bool found = false;

  //Selecting the part of DNA that will be moved to the kids
  index = (rand() % (parents[0].size()/2)) + 1;
  subsequence.push_back(index);
  index = (rand() % ((parents[0].size()-1)/2)) + parents[0].size()/2;
  if (index == subsequence[0]) ++index;
  subsequence.push_back(index);

  //Generating the children
  vector<int> children1(parents[0].size(),0);
  vector<int> children2(parents[1].size(),0);
  for(int i = subsequence[0]; i < subsequence[1]+1; i++){
    children1[i] = parents[1][i];
    children2[i] = parents[0][i];
  }
  
  //Copying the elements from the left
  for(int i = 0; i < subsequence[0];i++){
    found = false;
    for(int j = subsequence[0]; j <= subsequence[1]; j++){
      if(parents[0][i] == children1[j]){
        found = true;
        break;
      }
    }
    if(!found) children1[i] = parents[0][i];
    
    found = false;
    for(int j = subsequence[0]; j <= subsequence[1]; j++){
      if(parents[1][i] == children2[j]){
        found = true;
        break;
      }
    }
    if(!found) children2[i] = parents[1][i];
  }
  
  //Copying the elements from the left
  for(int i = subsequence[1]+1; i < parents[0].size();i++){
    found = false;
    for(int j = subsequence[0]; j <= subsequence[1]; j++){
      if(parents[0][i] == children1[j]){
        found = true;
        break;
      }
    }
    if(!found) children1[i] = parents[0][i];

    found = false;
    for(int j = subsequence[0]; j <= subsequence[1]; j++){
      if(parents[1][i] == children2[j]){
        found = true;
        break;
      }
    }
    if(!found) children2[i] = parents[1][i];
  }

  //Improving the solution for the problem generation
  n = count(children1.begin(), children1.end(), 0);
  for (int i = 0; i < n; i++)
  {
    it = find (children1.begin(), children1.end(), 0);
    index = distance(children1.begin(), it);
    indexTest = index;
    while (it != children1.end())
    {
      it = find (children1.begin(), children1.end(), parents[1][indexTest]);
      if(it != children1.end()){
        it = find (parents[0].begin(), parents[0].end(), parents[1][indexTest]);
        indexTest = distance(parents[0].begin(), it);
      }
      else{
        children1[index] = parents[1][indexTest];
      }
    }
  }

  n = count(children2.begin(), children2.end(), 0);
  for (int i = 0; i < n; i++)
  {
    it = find (children2.begin(), children2.end(), 0);
    index = distance(children2.begin(), it);
    indexTest = index;
    while (it != children2.end())
    {
      it = find (children2.begin(), children2.end(), parents[0][indexTest]);
      if(it != children2.end()){
        it = find (parents[1].begin(), parents[1].end(), parents[0][indexTest]);
        indexTest = distance(parents[1].begin(), it);
      }
      else{
        children2[index] = parents[0][indexTest];
      }
    }
  }

  kid.cities = children1; 
  kid.distance = getDistance(children1);
  childrens.push_back(kid);
  kid.cities = children2; 
  kid.distance = getDistance(children2);
  childrens.push_back(kid);
}

void crossover(vector<solution> &population, int dimension){

  bool chromosomeUsed[dimension] = {false};
  vector<vector<int> > parents;
  vector<solution> kids;
  vector<solution> newGeneration;

  //Generating children
  for(int i = 0, n; i < dimension; i = i +2){
    for(int j = 0; j < 2; j++){
      do{
        n = rand()%dimension;// x will be the 4 digit number  
      }while(chromosomeUsed[n]);
      parents.push_back(population[n].cities);
      chromosomeUsed[n] = true;
    }
    //Reproducing
    PMX(parents, kids);
    parents.clear();
  }
  //Seleting best chromosomes
  newGeneration.reserve( population.size() + kids.size() ); // preallocate memory
  newGeneration.insert( newGeneration.end(), population.begin(), population.end() );
  newGeneration.insert( newGeneration.end(), kids.begin(), kids.end() );
  sort(newGeneration.begin(), newGeneration.end());
  newGeneration.erase(newGeneration.begin()+dimension, newGeneration.end());
}

void mutation(vector<solution> & pop, double rate, int dimension){

  /* Variables that control the randomness of the solution*/
  std::default_random_engine generator (time(NULL));
  std::discrete_distribution<int> distribution {1-rate,rate};
  int indexs[2];
  for(int i = 0; i < pop.size(); i++){
    if(distribution(generator)){
      indexs[0] = rand() % dimension;
      if(indexs[0] == 0) indexs[0] = 1;
      indexs[1] = rand() % dimension;
      if(indexs[1] == dimension) indexs[1] = dimension-1;
      swap(pop[i].cities[indexs[0]], pop[i].cities[indexs[1]]);
      pop[i].distance = getDistance(pop[i].cities);
    }    
  }
}

int main(int argc, char** argv) {

  unsigned int populationSize = 50;
  int numberOfGenerations = 20;
  double mutationRate = 0.01;
  vector<solution> population;
  solution bestRoute;

  /* Generate seed randomly*/
  srand (time(NULL));
  readData(argc, argv, &dimension, &matrizAdj);

  population = randomConstruction(populationSize);

  for (int generation = 1; generation <= numberOfGenerations; generation++)
  {
    crossover(population, populationSize);
    mutation(population, mutationRate, dimension);
  }

  sort(population.begin(), population.end());
  bestRoute = population.front();
  printList(bestRoute.cities); cout << "\t" << bestRoute.distance << endl;

  return 0;  
}