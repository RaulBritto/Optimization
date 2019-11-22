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
int IMAX = 1;

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

void printCities(vector<int> &list){
  for (auto i = list.begin(); i != list.end(); ++i) 
    cout << *i << " "; 
  cout << endl;
}

double getDistance(vector<int> &listaCidade){
  double distance = 0;

  for(int i = 0; i < listaCidade.size()- 1; i++){
      distance += matrizAdj[listaCidade.at(i)][listaCidade.at(i+1)];
  }
  cout << "distancia " << distance << endl;
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
  
  do{
    delta = 0;
    deltaMinimium = 0;

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
      //cout << "[" <<s.at(i) << "," << s.at(j) << "]" << "  delta = " << delta <<endl;
      //cout  <<  i << "," << j << " " << distancia  <<endl;
      //if(deltaMinimium > delta) deltaMinimium = delta;
      }
    }
    //Se houve melhora
    if(deltaMinimium < 0){
      cout << "Menor delta: " << deltaMinimium << " [Nós: " << s[firstNode] << " " << s[secondNode] << "]" <<endl;
      distancia = distancia + deltaMinimium;
      swap(s[firstNode], s[secondNode]);
      printCities(s);
      cout << " NewDistance " << distancia << endl;
    }
  }while(deltaMinimium < 0);  
      return distancia;      
}

int RVND(vector<int> &s){
  
  vector<int> NL;
  int distance = 0;
  //set some values:
  for (int i=1; i<=5; ++i) NL.push_back(i);
  //using built-in random generator:
  random_shuffle(NL.begin(), NL.end());

  distance = getDistance(s);

  while(!NL.empty()){
    switch (NL.front()){ 
      case 1:
        //cout << "Swap" << endl;
        Swap(s, distance);
        break;
      case 2:
        //cout << "2-opt" << endl;
        break;
      case 3:
        //cout << "Reinsertion" << endl;
        break;
      case 4:
        //cout << "Or-2-opt" << endl;
        break;
      case 5:
        //cout << "Or-3-opt" << endl;
        break;
    }
    NL.erase(NL.begin());
  }

  return distance;
}



int main(int argc, char** argv) {

  double alpha;
  vector<int>  s;
  int distance = 0;
  vector<int>  s_;
  int distance_ = 0;

  //TODO: ALTERAR 1 para dimension
  int IILS = min(1,100); 
  
  /* Variavéis de controle da aleatorieadade*/
  srand (time(NULL));
 
  readData(argc, argv, &dimension, &matrizAdj);
  //printData();
  for(int i = 0; i < IMAX; i++){
    alpha = R[rand()%R.size()];
    s = Construction(alpha);
    s_ = s;

    for(int interILS = 0; interILS < IILS; interILS++){
      printCities(s);
      distance = RVND(s);
      //if(interILS == 9) interILS = 0;
    }
  }
    
  return 0;  
}