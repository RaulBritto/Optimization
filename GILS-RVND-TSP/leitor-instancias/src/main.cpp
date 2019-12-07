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
int IMAX = 50;  //number of iteration
/*This struct helps the construction step*/
struct InsertionInfo{
  int insertedNode; // node k to be inserted
  int removedEdge; //edge {i,j} where the node k will be inserted
  double cost; //delta after insert the node k in the edge {i,j}

  bool operator() (InsertionInfo i,InsertionInfo j) { return (i.cost<j.cost);}
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

  vector<int>  s = {1,1};
  vector<int> candidateNodes;
  int initialSubtourSize = 2;

  // Creating a list with the candidate nodes starting from the city number 2.
  for (int i = 2; i <= dimension; i++) 
    candidateNodes.push_back(i); 

  // Creating an initial Subtour choosed randomly.
  for(int i = 0; i < initialSubtourSize; i++){
    int j = (rand()) % candidateNodes.size();
    s.insert(s.begin()+1, candidateNodes[j]);
    candidateNodes.erase(candidateNodes.begin() + j);
  }

  vector <InsertionInfo> insertionCost((s.size()-1)*candidateNodes.size());
  /* Variables that control the randomness of the solution*/
  std::default_random_engine generator (time(NULL));
  std::discrete_distribution<int> distribution {1-alpha,alpha};
  int choosed, number; 

  // Compute the insertion cost for all candidate Nodes
  while(candidateNodes.size()){
    for(int i = 1, j = i+1, l = 0 ; i <= s.size() - 1; i++, j++){
      for(auto k : candidateNodes){
        insertionCost[l].cost = matrizAdj[s[i-1]][k] + matrizAdj[s[j-1]][k] - matrizAdj[s[i-1]][s[j-1]];
        insertionCost[l].insertedNode = k;
        insertionCost[l].removedEdge = i ;
        l++;
      }
    }
    //Sort the vector InsertionCost by ascending order
    std::sort(insertionCost.begin(), insertionCost.end(), InsertionInfo());

    //Check if the next node inserted will be the cheapest insertion or a random node.
    number = distribution(generator);
    if(!number) choosed = 0;
    else choosed = rand()%insertionCost.size();

    //Insert the node choosed and remove from candidate nodes
    s.insert(s.begin()+insertionCost[choosed].removedEdge, insertionCost[choosed].insertedNode);
    candidateNodes.erase(remove(candidateNodes.begin(), candidateNodes.end(), insertionCost[choosed].insertedNode), candidateNodes.end());
    insertionCost.resize((s.size()-1)*candidateNodes.size());
  } 

  return s;
}

/*This function calculates all the possibles 'swap' between the nodes and chooses the best one (Descresing the total distance).*/
int Swap(vector<int> &s, double distancia){

  double delta = 0, deltaMinimium = 0, deltaFixed = 0;
  int firstNode = 0, secondNode = 0, length = 0;

    length = s.size();
    for(int i = 1, n = 1; i < length-2; i++){
      deltaFixed = - matrizAdj[s[i-1]][s[i]] - matrizAdj[s[i]][s[i+1]];   
      for(int j = i+1; j < length-1; j++){  
        if(j == i+1){
          delta = deltaFixed - matrizAdj[s[j]][s[j+1]]
                    + matrizAdj[s[i-1]][s[j]] + matrizAdj[s[j]][s[i]] + matrizAdj[s[i]][s[j+1]];
          if(deltaMinimium > delta){
            deltaMinimium = delta;
            firstNode = i; secondNode = j;  
          }     
        } 
        else{
          delta = deltaFixed - matrizAdj[s[j-1]][s[j]] - matrizAdj[s[j]][s[j+1]]
                  + matrizAdj[s[i-1]][s[j]] + matrizAdj[s[j]][s[i+1]] + matrizAdj[s[j-1]][s[i]] + matrizAdj[s[i]][s[j+1]];
          if(deltaMinimium > delta){
            deltaMinimium = delta;
            firstNode = i; secondNode = j;  
          }      
        } 
      }
    }
    //Check if exists an improvement
    if(deltaMinimium < 0){
      distancia = distancia + deltaMinimium;
      swap(s[firstNode], s[secondNode]);
    }
      return distancia;      
}
/*This function calculates all the reverse subsenquences in the solution and chooses the best one (Descresing the total distance).*/
int _2opt(vector<int> &s, double distancia){

  double delta = 0, deltaMinimium = 0, deltaFixed = 0;
  int firstNode = 0, secondNode = 0, length = 0;
  vector<int> subsequence;

    length = s.size();
    for(int i = 1 ; i < length-2; i++){
      deltaFixed = - matrizAdj[s[i-1]][s[i]];
      for(int j = i +1; j < length-1; j++){
        delta = deltaFixed - matrizAdj[s[j]][s[j+1]]
                  + matrizAdj[s[i-1]][s[j]] + matrizAdj[s[i]][s[j+1]];
        if(deltaMinimium > delta){
          deltaMinimium = delta;
          firstNode = i; secondNode = j;  
        }
      }
    }
    //Check if exists an improvement
    if(deltaMinimium < 0){
      //creting new reverse sequence 
      for(int i = secondNode; i >= firstNode; i--) subsequence.push_back(s[i]);
      //replacing the sequence with the old one
      s.insert(s.begin()+ firstNode , subsequence.begin(),subsequence.end());
      s.erase(s.begin() + firstNode + subsequence.size()  , s.begin() + firstNode + 2*subsequence.size());
      //updating the new distance after moving the sequence
      distancia = distancia + deltaMinimium;
      subsequence.clear();
    }

  return distancia; 
}

/*This function calculates the best insertion of a subsequence of size 'k' in the solution and chooses the best one (Descresing the total distance).*/
int orkOpt(vector<int> &s, double distancia, int k){

  vector<int> subsequence;
  double delta = 0, deltaMinimium = 0, deltaFixed = 0;
  int firstNode = 0, secondNode = 0, length = 0;

    length = s.size();
    for(int i = 1; i < length -k; i++){
      deltaFixed = - matrizAdj[s[i-1]][s[i]] - matrizAdj[s[i+k-1]][s[i+k]] + matrizAdj[s[i-1]][s[i+k]];
      //cout << deltaFixed << endl;
      for(int j = 0; j < length -1; j++){
        if(j >= i-1 && j < i+k) continue;
        else{
          delta = deltaFixed + matrizAdj[s[j]][s[i]] + matrizAdj[s[i+k-1]][s[j+1]] - matrizAdj[s[j]][s[j+1]];
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
        subsequence.push_back(s[j+firstNode]);
      }
      s.erase(s.begin()+firstNode, s.begin()+firstNode+k);  
      if(firstNode < secondNode){
        s.insert(s.begin()+secondNode+1-k, subsequence.begin(), subsequence.end());
      }else{
        s.insert(s.begin()+secondNode+1, subsequence.begin(), subsequence.end());  
      }
      distancia = distancia + deltaMinimium;

    } 

  return distancia; 
}
 /*This function pertubs the solution (chooses two subsequences with random size with no matching and swaps both).*/
void Pertub(vector<int> &s){

  vector<int> subsequence;
  std::vector<int>::iterator it;
  int index;

  //printList(s);

  //Choose the first index
  index = (rand() % (s.size()-2)) + 1;
  //Choose the next three indexs randomly
  for(int i = 0; i <= 3; i++){
    do{
      index = (rand() % (s.size()-2)) + 1;
      it = find (subsequence.begin(), subsequence.end(), index);
    }while(it != subsequence.end());
    subsequence.push_back(index);
  } 
  std::sort (subsequence.begin(), subsequence.end());
  
  if(subsequence[1] - subsequence[0] > 10) subsequence[1] = subsequence[1]+10;
  if(subsequence[3] - subsequence[2] > 10) subsequence[3] = subsequence[2]+10;
  //printList(subsequence);
  
  //Swap the subvectors sequences
  vector<int> subvector1 = std::vector<int>(s.begin()+subsequence[0],s.begin()+subsequence[1]+1);
  vector<int> subvector2 = std::vector<int>(s.begin()+subsequence[2],s.begin()+subsequence[3]+1);
  s.insert(s.begin()+ subsequence[2] , subvector1.begin(),subvector1.end());
  s.erase(s.begin() + subsequence[2] + subvector1.size(), s.begin() + subsequence[3]+1 +  subvector1.size());
  s.erase(s.begin() + subsequence[0], s.begin() + subsequence[1]+1);
  s.insert(s.begin()+ subsequence[0] , subvector2.begin(),subvector2.end());

  //printList(s);
  
}

/*This function implements the Random VND and return the minimum distance found.*/
int RVND(vector<int> &s){
  
  vector<int> NL;
  vector<int> s_;
  int distance = 0, distance_ = 0, choosed = 0;
  //set some values:
  for (int i=1; i<=5; ++i) NL.push_back(i);
  //using built-in random generator:
  random_shuffle(NL.begin(), NL.end());

  s_ = s;
  distance = getDistance(s);
  distance_ = distance;

  while(!NL.empty()){
    switch (NL.front()){ 
      case 1:
        //cout << "Swap" << endl;
        distance_ = Swap(s_, distance_);
        break;
      case 2:
        //cout << "2-opt" << endl;
        distance_ = _2opt(s_, distance_);
        break;
      case 3:
        //cout << "Reinsertion" << endl;
        distance_ = orkOpt(s_, distance_, 1);
        break;
      case 4:
        //cout << "Or-2-opt" << endl;
        distance_ = orkOpt(s_, distance_, 2);
        break;
      case 5:
        //cout << "Or-3-opt" << endl;
        distance_ = orkOpt(s_, distance_, 3);
        break;
    }
    if(distance_ < distance){
      s = s_;
      distance = distance_;
      NL.clear();
      for (int i=1; i<=5; ++i) NL.push_back(i);
      random_shuffle(NL.begin(), NL.end());
    }else{
      NL.erase(NL.begin());
    }
  }

  return distance;
}



int main(int argc, char** argv) {

  // Get starting timepoint 
  auto start = high_resolution_clock::now(); 

  double alpha;
  vector<int>  s;
  int distance = 0, distance_ = 0;
  vector<int>  s_;
  int f = 10000000;
  vector<int> bestRoute;

  srand (time(NULL));
 
  readData(argc, argv, &dimension, &matrizAdj);
  //TODO put dimension
  int IILS = 0 ;
  if(dimension >= 150){
    IILS = dimension/2;
  }else{
    IILS = dimension;
  }
  
  
  //printData();
  for(int i = 0; i < IMAX; i++){
    alpha = R[rand()%R.size()];
    s = Construction(alpha);
    s_ = s;
    distance_ = getDistance(s_);

    for(int iterILS = 0; iterILS < IILS; iterILS++){
      distance = RVND(s);
      if(distance < distance_){
        s_ = s;
        distance_ = distance;
        iterILS = 0;
      }
      Pertub(s_);
    }
    if(distance_ < f){
      bestRoute = s_;
      f = distance_;
    }
  }
  // Get ending timepoint 
  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 

  //cout << "Best solution:" << endl;
  //printList(bestRoute);
  cout << f << " " << duration.count()/1000000.0 << endl;
  

  return 0;  
}