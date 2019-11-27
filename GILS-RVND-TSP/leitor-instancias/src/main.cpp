#include "readData.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>


using namespace std;

double ** matrizAdj;
int dimension; // quantidade total de vertices
vector<double> R = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10
                    , 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21
                    , 0.21, 0.22, 0.23, 0.24, 0.25};
int IMAX = 10;

struct InsertionInfo{
  int noInserido; // no k a ser inserido
  int arestaRemovida; //aresta {i,j} onde o no k sera inserido
  double custo; //delta ao inserir k na aresta {i,j}

  bool operator() (InsertionInfo i,InsertionInfo j) { return (i.custo<j.custo);}
};

void printData() {
  cout << "dimension: " << dimension << endl;
  for (size_t i = 1; i <= dimension; i++) {
    for (size_t j = 1; j <= dimension; j++) {
      cout << matrizAdj[i][j] << " ";
    }
    cout << endl;
  }
}

void printList(vector<int> &list){
  for (auto i = list.begin(); i != list.end(); ++i) 
    cout << *i << " "; 
  cout << endl;
}

double getDistance(vector<int> &listaCidade){
  double distance = 0;

  for(int i = 0; i < listaCidade.size()- 1; i++){
      distance += matrizAdj[listaCidade.at(i)][listaCidade.at(i+1)];
  }
  //cout << "distancia " << distance << endl;
  return distance; 
}

vector<int>  Construction(double &alpha){ 

  vector<int>  s = {1,1};
  vector<int> listaDeCandidatos;
  int tamanhoSubtourInicial = 2;

  // Criando lista de cidades candidatas (indices) a partir da cidade 2
  for (int i = 2; i <= dimension; i++) 
    listaDeCandidatos.push_back(i); 

  //Criando subtour com as cidades candidatas
  for(int i = 0; i < tamanhoSubtourInicial; i++){
    int j = (rand()) % listaDeCandidatos.size();
    s.insert(s.begin()+1, listaDeCandidatos[j]);
    listaDeCandidatos.erase(listaDeCandidatos.begin() + j);
  }

  vector <InsertionInfo> custoInsercao((s.size()-1)*listaDeCandidatos.size());
  /* Variavéis de controle da aleatorieadade*/
  std::default_random_engine generator (time(NULL));
  std::discrete_distribution<int> distribution {1-alpha,alpha};
  int choosed, number; 

  while(listaDeCandidatos.size()){
    for(int i = 1, j = i+1, l = 0 ; i <= s.size() - 1; i++, j++){
      for(auto k : listaDeCandidatos){
        custoInsercao[l].custo = matrizAdj[s[i-1]][k] + matrizAdj[s[j-1]][k] - matrizAdj[s[i-1]][s[j-1]];
        custoInsercao[l].noInserido = k;
        custoInsercao[l].arestaRemovida = i ;
        l++;
      }
    }

    std::sort(custoInsercao.begin(), custoInsercao.end(), InsertionInfo());

    number = distribution(generator);
    if(!number) choosed = 0;
    else choosed = rand()%custoInsercao.size();

    s.insert(s.begin()+custoInsercao[choosed].arestaRemovida, custoInsercao[choosed].noInserido);
    listaDeCandidatos.erase(remove(listaDeCandidatos.begin(), listaDeCandidatos.end(), custoInsercao[choosed].noInserido), listaDeCandidatos.end());
    custoInsercao.resize((s.size()-1)*listaDeCandidatos.size());
  } 

  return s;
}

int Swap(vector<int> &s, double distancia){

  double delta = 0, deltaMinimium = 0;
  int firstNode = 0, secondNode = 0;
  
    for(int n = 0, i = 1; n < (s.size()-1)*(s.size()-2)/2; n++,i++){
      for(int j = i+1; j < s.size()-1; j++){  
        if(j == i+1){
          delta = - matrizAdj[s[i-1]][s[i]] - matrizAdj[s[i]][s[i+1]] - matrizAdj[s[j]][s[j+1]]
                    + matrizAdj[s[i-1]][s[j]] + matrizAdj[s[j]][s[i]] + matrizAdj[s[i]][s[j+1]];
          if(deltaMinimium > delta){
            deltaMinimium = delta;
            firstNode = i; secondNode = j;  
          }     
        } 
        else{
          delta = - matrizAdj[s[i-1]][s[i]] - matrizAdj[s[i]][s[i+1]] - matrizAdj[s[j-1]][s[j]] - matrizAdj[s[j]][s[j+1]]
                  + matrizAdj[s[i-1]][s[j]] + matrizAdj[s[j]][s[i+1]] + matrizAdj[s[j-1]][s[i]] + matrizAdj[s[i]][s[j+1]];
          if(deltaMinimium > delta){
            deltaMinimium = delta;
            firstNode = i; secondNode = j;  
          }      
        } 
      }
    }
    //Se houve melhora
    if(deltaMinimium < 0){
      distancia = distancia + deltaMinimium;
      swap(s[firstNode], s[secondNode]);
    }
      return distancia;      
}

int _2opt(vector<int> &s, double distancia){

  double delta = 0, deltaMinimium = 0;
  int firstNode = 0, secondNode = 0;
  vector<int> subsequence;

    for(int i = 1 ; i < s.size()-2; i++){
      for(int j = i +1; j < s.size()-1; j++){
        delta = - matrizAdj[s[i-1]][s[i]] - matrizAdj[s[j]][s[j+1]]
                  + matrizAdj[s[i-1]][s[j]] + matrizAdj[s[i]][s[j+1]];
        if(deltaMinimium > delta){
          deltaMinimium = delta;
          firstNode = i; secondNode = j;  
        }
      }
    }
    //Se houve melhora
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

int orkOpt(vector<int> &s, double distancia, int k){

  vector<int> subsequence;
  vector<int> s_;
  double delta = 0, deltaMinimium = 0;

    for(int i = 1; i < s.size() -k; i++){
      //Creating the subsequence that will be moved  
      for(int j = 0; j < k; j++ ){
        subsequence.push_back(s[i+j]);
      }
      //Removing the subsequence from the solution s
      s.erase(s.begin()+i, s.begin()+i+subsequence.size());
      //To test the new possibilities
      for(int j = 1; j <= s.size()-1;j++){
        //this move generates the same sequence of cities
        if(i == j) continue;
        delta = - matrizAdj[s[i-1]][subsequence[0]] - matrizAdj[subsequence[k-1]][s[i]] - matrizAdj[s[j-1]][s[j]]
                + matrizAdj[s[i-1]][s[i]] + matrizAdj[s[j-1]][subsequence[0]] + matrizAdj[subsequence[k-1]][s[j]];
        if(deltaMinimium > delta){
          deltaMinimium = delta;
          s.insert(s.begin()+j, subsequence.begin(), subsequence.end());
          s_ = s;
          s.erase(s.begin()+j, s.begin()+j+subsequence.size());
        }
      }
      s.insert(s.begin()+i , subsequence.begin(), subsequence.end());
      subsequence.clear();
    }
    //If an improvement exists
    if(deltaMinimium < 0){
      distancia = distancia + deltaMinimium;
      s = s_;
    } 

  return distancia; 
}

vector<int> Pertub(vector<int> s){

  vector<int> subsequence;
  std::vector<int>::iterator it;
  int index;

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
  
  //Swap the subvectors sequences
  vector<int> subvector1 = std::vector<int>(s.begin()+subsequence[0],s.begin()+subsequence[1]+1);
  vector<int> subvector2 = std::vector<int>(s.begin()+subsequence[2],s.begin()+subsequence[3]+1);
  s.insert(s.begin()+ subsequence[2] , subvector1.begin(),subvector1.end());
  s.erase(s.begin() + subsequence[2] + subvector1.size(), s.begin() + subsequence[3]+1 + + subvector1.size());
  s.erase(s.begin() + subsequence[0], s.begin() + subsequence[1]+1);
  s.insert(s.begin()+ subsequence[0] , subvector2.begin(),subvector2.end());
  
  return s;
}

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

  double alpha;
  vector<int>  s;
  int distance = 0, distance_ = 0;
  vector<int>  s_;
  int f = 10000000;
  vector<int> bestRoute;

  /* Variavéis de controle da aleatorieadade*/
  srand (time(NULL));
 
  readData(argc, argv, &dimension, &matrizAdj);
  //TODO: ALTERAR 1 para dimension
  int IILS = min(dimension,100); 
  
  //printData();
  for(int i = 0; i < IMAX; i++){
    alpha = R[rand()%R.size()];
    s = Construction(alpha);
    s_ = s;
    distance_ = getDistance(s_);

    //cout << "Solução inicial: ";
    //printList(s_);
    //cout << "Menor distancia: " << distance_ << endl;

    for(int iterILS = 0; iterILS < IILS; iterILS++){
      distance = RVND(s);
      //cout << iterILS << ")" << endl;
      //printList(s);
      //cout << "Menor distancia: " << distance_ << endl;
      //cout << "distance apos otimizacao: " << distance << endl;
      if(distance < distance_){
        s_ = s;
        distance_ = distance;
        iterILS = 0;
      }
      s = Pertub(s_);
      distance = getDistance(s);
      //cout << "Pertubando:" << endl;
      //printList(s);
      //cout << "distance: " << distance << endl;
    }
    if(distance_ < f){
      bestRoute = s_;
      f = distance_;
    }
    //cout << i << ")" << getDistance(s_) << endl;
  }

  cout << "Melhor solucao:" << endl;
  printList(bestRoute);
  cout << f << endl;
    
  return 0;  
}