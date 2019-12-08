#include "readData.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>


using namespace std;
using namespace std::chrono; 

double ** matrizAdj;
int dimension; // number of vertices in the problem
/*This array is responsable to choose the alpha variable to control the randomness level in the construction step*/
vector<double> R = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10
                    , 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21
                    , 0.21, 0.22, 0.23, 0.24, 0.25};
//TODO change to IMAX                    
int IMAX = 1;  //number of iteration
/*This struct helps the construction step*/
struct InsertionInfo{
  int insertedNode; // node k to be inserted
  int removedEdge; //edge {i,j} where the node k will be inserted
  double cost; //delta after insert the node k in the edge {i,j}

  bool operator() (InsertionInfo i,InsertionInfo j) { return (i.cost<j.cost);}
};

struct solution{
  vector<int> cities;
  double distance; //delta after insert the node k in the edge {i,j}
};

/*This function print all the adjacency matrix*/
void printData() {
  cout << "dimension: " << dimension << endl;
  for (size_t i = 1; i <= dimension; i++) {
    for (size_t j = 1; j <= dimension; j++) {
      cout << matrizAdj[i][j] << " ";
    }
    cout << endl;
  }
}
/*This function prints all the elements of a list*/
void printList(vector<int> &list){
  for (auto i = list.begin(); i != list.end(); ++i) 
    cout << *i << " "; 
  cout << endl;
}
/*This function calculates the distance given a solution of TSP problem.*/
double getDistance(vector<int> &listaCidade){
  double distance = 0;

  for(int i = 0; i < listaCidade.size()- 1; i++){
      distance += matrizAdj[listaCidade.at(i)][listaCidade.at(i+1)];
  }
  return distance; 
}
/*This function builds a initial solution of TSP problem using the cheapest insertion heuristic.
The solution returned by this function is done using some level of randomness given by alpha.*/
vector<int>  Construction(double &alpha){ 

  vector<int>  s = {1,1}; //Initial Solution
  vector<int> candidateNodes; //List of candidate nodes
  /* Variables that control the randomness of the solution*/
  std::default_random_engine generator (time(NULL));
  std::discrete_distribution<int> distribution {1-alpha,alpha};
  int choosed, number; 
  //Compute the best insertion in the solution
  double cost = 0, lowerCost = 100000000;
  vector <InsertionInfo> bestCandidates; //Vector that stores the best candidates in case of tie
  InsertionInfo candidate;

  // Creating a list with the candidate nodes starting from the city number 2.
  for (int i = 2; i <= dimension; i++) 
    candidateNodes.push_back(i); 

  
  //Calculate the best node to insertion in the entire candidate nodes
  while(candidateNodes.size()){
    lowerCost = 100000000;
    for(int i = 0; i < s.size()-1;i++){
      for(int j = 0; j < candidateNodes.size();j++){
        cost = matrizAdj[s[i]][candidateNodes[j]] + matrizAdj[candidateNodes[j]][s[i+1]] - matrizAdj[s[i]][s[i+1]];
        if(cost < lowerCost){
          candidate.cost = cost;
          candidate.insertedNode = j;
          candidate.removedEdge = i;
          bestCandidates.clear();
          bestCandidates.insert(bestCandidates.begin(), candidate);
          lowerCost = cost;
        }else if (cost == lowerCost)
        {
          candidate.cost = cost;
          candidate.insertedNode = j;
          candidate.removedEdge = i;
          bestCandidates.push_back(candidate);
        }
      }
    }
    //Check if the next node inserted will be the cheapest insertion or a random node.
    number = distribution(generator);
    if(!number){
      choosed = rand() % bestCandidates.size();
      s.insert(s.begin()+1+bestCandidates[choosed].removedEdge , candidateNodes[bestCandidates[choosed].insertedNode]);
      candidateNodes.erase(candidateNodes.begin() + bestCandidates[choosed].insertedNode);
    } 
    else{
      choosed = rand()%candidateNodes.size();
      s.insert(s.begin()+1 + (rand()%(s.size()-1)), candidateNodes[choosed]);
      candidateNodes.erase(candidateNodes.begin() + choosed);  
    } 
  }

  return s;
}

/*This function calculates all the possibles 'swap' between the nodes and chooses the best one (Descresing the total distance).*/
void Swap(solution &s){

  double delta = 0, deltaMinimium = 0, deltaFixed = 0;
  int firstNode = 0, secondNode = 0, length = 0;

    length = s.cities.size();
    for(int i = 1; i < length-2; i++){
      deltaFixed = - matrizAdj[s.cities[i-1]][s.cities[i]] - matrizAdj[s.cities[i]][s.cities[i+1]];   
      for(int j = i+1; j < length-1; j++){  
        if(j == i+1){
          delta = deltaFixed - matrizAdj[s.cities[j]][s.cities[j+1]]
                    + matrizAdj[s.cities[i-1]][s.cities[j]] + matrizAdj[s.cities[j]][s.cities[i]] + matrizAdj[s.cities[i]][s.cities[j+1]];
          if(deltaMinimium > delta){
            deltaMinimium = delta;
            firstNode = i; secondNode = j;  
          }     
        } 
        else{
          delta = deltaFixed - matrizAdj[s.cities[j-1]][s.cities[j]] - matrizAdj[s.cities[j]][s.cities[j+1]]
                  + matrizAdj[s.cities[i-1]][s.cities[j]] + matrizAdj[s.cities[j]][s.cities[i+1]] + matrizAdj[s.cities[j-1]][s.cities[i]] + matrizAdj[s.cities[i]][s.cities[j+1]];
          if(deltaMinimium > delta){
            deltaMinimium = delta;
            firstNode = i; secondNode = j;  
          }      
        } 
      }
    }
    //Check if exists an improvement
    if(deltaMinimium < 0){
      s.distance += deltaMinimium;
      swap(s.cities[firstNode], s.cities[secondNode]);
    }
}
/*This function calculates all the reverse subsenquences in the solution and chooses the best one (Descresing the total distance).*/
void _2opt(solution &s){

  double delta = 0, deltaMinimium = 0, deltaFixed = 0;
  int firstNode = 0, secondNode = 0, length = 0;
  vector<int> subsequence;

    length = s.cities.size();
    for(int i = 1 ; i < length-2; i++){
      deltaFixed = - matrizAdj[s.cities[i-1]][s.cities[i]];
      for(int j = i +1; j < length-1; j++){
        delta = deltaFixed - matrizAdj[s.cities[j]][s.cities[j+1]]
                  + matrizAdj[s.cities[i-1]][s.cities[j]] + matrizAdj[s.cities[i]][s.cities[j+1]];
        if(deltaMinimium > delta){
          deltaMinimium = delta;
          firstNode = i; secondNode = j;  
        }
      }
    }
    //Check if exists an improvement
    if(deltaMinimium < 0){
      //creting new reverse sequence 
      for(int i = secondNode; i >= firstNode; i--) subsequence.push_back(s.cities[i]);
      //replacing the sequence with the old one
      s.cities.insert(s.cities.begin()+ firstNode , subsequence.begin(),subsequence.end());
      s.cities.erase(s.cities.begin() + firstNode + subsequence.size()  , s.cities.begin() + firstNode + 2*subsequence.size());
      //updating the new distance after moving the sequence
      s.distance +=  deltaMinimium;
    }
}

/*This function calculates the best insertion of a subsequence of size 'k' in the solution and chooses the best one (Descresing the total distance).*/
void orkOpt(solution &s, int k){

  vector<int> subsequence;
  double delta = 0, deltaMinimium = 0, deltaFixed = 0;
  int firstNode = 0, secondNode = 0, length = 0;

    length = s.cities.size();
    for(int i = 1; i < length -k; i++){
      deltaFixed = - matrizAdj[s.cities[i-1]][s.cities[i]] - matrizAdj[s.cities[i+k-1]][s.cities[i+k]] + matrizAdj[s.cities[i-1]][s.cities[i+k]];
      //cout << deltaFixed << endl;
      for(int j = 0; j < length -1; j++){
        if(j >= i-1 && j < i+k) continue;
        else{
          delta = deltaFixed + matrizAdj[s.cities[j]][s.cities[i]] + matrizAdj[s.cities[i+k-1]][s.cities[j+1]] - matrizAdj[s.cities[j]][s.cities[j+1]];
          if(deltaMinimium > delta){
            deltaMinimium = delta;
            firstNode = i; secondNode = j;  
          }
        }
      }
    }
    
    //If an improvement exists
    if(deltaMinimium < 0){
      //Creating the subsequence that will be moved  
      for(int j = 0; j < k; j++ ){
        subsequence.push_back(s.cities[j+firstNode]);
      }
      s.cities.erase(s.cities.begin()+firstNode, s.cities.begin()+firstNode+k);  
      if(firstNode < secondNode){
        s.cities.insert(s.cities.begin()+secondNode+1-k, subsequence.begin(), subsequence.end());
      }else{
        s.cities.insert(s.cities.begin()+secondNode+1, subsequence.begin(), subsequence.end());  
      }
      s.distance += deltaMinimium;
    } 
}
 /*This function pertubs the solution (chooses two subsequences with random size with no matching and swaps both).*/
solution Pertub(solution s){

  vector<int> subsequence;
  std::vector<int>::iterator it;
  int index;

  //Choose the first index
  index = (rand() % (s.cities.size()-2)) + 1;
  //Choose the next three indexs randomly
  for(int i = 0; i <= 3; i++){
    do{
      index = (rand() % (s.cities.size()-2)) + 1;
      it = find (subsequence.begin(), subsequence.end(), index);
    }while(it != subsequence.end());
    subsequence.push_back(index);
  } 
  std::sort (subsequence.begin(), subsequence.end());
  if(subsequence[1] - subsequence[0] > 10) subsequence[1] = subsequence[0]+10;
  if(subsequence[3] - subsequence[2] > 10) subsequence[3] = subsequence[2]+10;

  //Swap the subvectors sequences
  vector<int> subvector1 = std::vector<int>(s.cities.begin()+subsequence[0],s.cities.begin()+subsequence[1]+1);
  vector<int> subvector2 = std::vector<int>(s.cities.begin()+subsequence[2],s.cities.begin()+subsequence[3]+1);
  s.cities.insert(s.cities.begin()+ subsequence[2] , subvector1.begin(),subvector1.end());
  s.cities.erase(s.cities.begin() + subsequence[2] + subvector1.size(), s.cities.begin() + subsequence[3]+1  + subvector1.size());
  s.cities.erase(s.cities.begin() + subsequence[0], s.cities.begin() + subsequence[1]+1);
  s.cities.insert(s.cities.begin()+ subsequence[0] , subvector2.begin(),subvector2.end());

  return s;
}

/*This function implements the Random VND and return the minimum distance found.*/
void RVND(solution &s){
  
  vector<int> NL;
  solution s_;
  //set some values:
  for (int i=1; i<=5; ++i) NL.push_back(i);
  //using built-in random generator:
  random_shuffle(NL.begin(), NL.end());

  s_ = s;

  while(!NL.empty()){
    switch (NL.front()){ 
      case 1:
        //cout << "Swap" << endl;
        Swap(s_);
        break;
      case 2:
        //cout << "2-opt" << endl;
        _2opt(s_);
        break;
      case 3:
        //cout << "Reinsertion" << endl;
        orkOpt(s_, 1);
        break;
      case 4:
        //cout << "Or-2-opt" << endl;
        orkOpt(s_, 2);
        break;
      case 5:
        //cout << "Or-3-opt" << endl;
        orkOpt(s_, 3);
        break;
    }
    if(s_.distance < s.distance){
      s = s_;
      NL.clear();
      for (int i=1; i<=5; ++i) NL.push_back(i);
      random_shuffle(NL.begin(), NL.end());
    }else{
      NL.erase(NL.begin());
    }
  }
}



int main(int argc, char** argv) {

  double alpha;
  solution  s;
  solution  s_;
  int f = 10000000;
  solution bestRoute;
  

  /* Generate seed randomly*/
  srand (time(NULL));
  readData(argc, argv, &dimension, &matrizAdj);

  int IILS = 0 ;
  if(dimension >= 150){
    IILS = dimension/2;
  }else{
    IILS = dimension;
  }
  
  
  // Get starting timepoint 
  auto start = high_resolution_clock::now(); 
  //printData();
  for(int i = 0; i < IMAX; i++){
    alpha = R[rand()%R.size()];
    auto inicio = high_resolution_clock::now(); 
    s.cities = Construction(alpha);
    s.distance = getDistance(s.cities);
    auto fim = high_resolution_clock::now(); 
    auto duracao = duration_cast<microseconds>(fim - inicio); 
    //cout << "Construção: " << duracao.count()/1000000.0 << endl;
    s_ = s;
    for(int iterILS = 0; iterILS < IILS; iterILS++){
      RVND(s);
      if(s.distance < s_.distance){
        s_ = s;
        iterILS = 0;
      }
      s = Pertub(s_);
      s.distance = getDistance(s.cities);
    }
    if(s_.distance < f){
      bestRoute = s_;
      f = s_.distance;
    }
  }
  // Get ending timepoint 
  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 

  //cout << "Best solution:" << endl;
  //printList(bestRoute);
  cout << f << ", " << duration.count()/1000000.0 << endl;
  //cout << getDistance(bestRoute);  


  return 0;  
}