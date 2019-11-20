#include "readData.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <time.h>      
#include <map>   
#include <random>
#include <tuple>

using namespace std;

double ** matrizAdj;
int dimension; // quantidade total de vertices
int tamanhoSubtourInicial = 2;
double alpha = 0.15;
int IMAX = 1;
int IILS = 1;
/*
int matrizAdj[7][7] =  {
                      {0,0,0,0,0,0,0},
                      {0,0,2,1,4,9,1},
                      {0,2,0,5,9,7,2},
                      {0,1,5,0,3,8,6},
                      {0,4,9,3,0,2,5},
                      {0,9,7,8,2,0,2},
                      {0,1,2,6,5,2,0}
    }; // matriz de adjacencia
*/


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

double calculaDistancia(vector<int> listaCidade){
  double distancia = 0;

  for(int i = 0; i < listaCidade.size()- 1; i++){
      distancia += matrizAdj[listaCidade.at(i)][listaCidade.at(i+1)];
  }

  cout << "distancia " << distancia << endl;

  return distancia; 
}

void imprimeCidade(vector<int> lista){
    for (auto i = lista.begin(); i != lista.end(); ++i) 
        cout << *i << " "; 
}

vector<int>&  Construction(vector<int> &s, vector<int> &listaDeCandidatos) 
{ 
    vector <InsertionInfo> custoInsercao((s.size()-1)*listaDeCandidatos.size());
    /* Variavéis de controle da aleatorieadade*/
    srand (time(NULL));
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

      //for(int i = 0; i < custoInsercao.size(); i++)
      //  cout << custoInsercao.at(i).noInserido << ' ' <<  custoInsercao.at(i).arestaRemovida << ' ' << custoInsercao.at(i).custo << endl;
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

    for(int n = 0, i = 0; n < (s.size()-1)*(s.size()-2)/2; n++,i++){
      for(int j = i+1; j < s.size()-1; j++){
        if(i==0){
      /*  
        //primeiro nó e o consectuivo
        if(j== 1) {
          delta = - matrizAdj[s[i]][s[j]] - matrizAdj[s[j]][s[j+1]] - matrizAdj[s.size()-1][i] 
                  + matrizAdj[s[j]][s[i]] + matrizAdj[s[i]][s[j+1]] + matrizAdj[s.size()-1][s[i]];
          deltaMinimium = delta;
          firstNode = i; secondNode = j;  
        }
        //primeiro nó e penúltimo  
        else if (j== s.size()-2){
          delta = - matrizAdj[s[i]][s[i+1]] - matrizAdj[s[j-1]][s[j]] - matrizAdj[s[j]][s[j+1]]
                  + matrizAdj[s[j]][s[i+1]] + matrizAdj[s[j-1]][s[i]] + matrizAdj[s[i]][s[j]];
          if(deltaMinimium > delta){
            deltaMinimium = delta;
            firstNode = i; secondNode = j;  
          } 
        } 
        else {
          //primeiro nó e os nós do meio
          delta = - matrizAdj[s[i]][s[i+1]] - matrizAdj[s.size()-1][s[i]] - matrizAdj[s[j-1]][s[j]] - matrizAdj[s[j]][s[j+1]]
                  + matrizAdj[s[j]][s[i+1]] + matrizAdj[s[j-1]][s[i]] + matrizAdj[s[i]][s[j+1]] + matrizAdj[s[s.size()-2]][s[j]];
          if(deltaMinimium > delta){
            deltaMinimium = delta;
            firstNode = i; secondNode = j;  
          } 
        }
      */
        }
        else{
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

      
      //cout << "[" <<s.at(i) << "," << s.at(j) << "]" << "  delta = " << delta <<endl;
      //cout  <<  i << "," << j << " " << distancia  <<endl;
      //if(deltaMinimium > delta) deltaMinimium = delta;

    }
  }
    //Se houve melhora
    if(deltaMinimium < 0){
      cout << "Menor delta: " << deltaMinimium << " [Nós: " << s[firstNode] << " " << s[secondNode] << "]" <<endl;
      distancia = distancia + deltaMinimium;
      if(!firstNode){
        swap(s[firstNode], s[secondNode]);
        s[s.size()-1] = s[0];
      } 
      else{
        swap(s[firstNode], s[secondNode]);
      }
      imprimeCidade(s);
      cout << " NewDistance " << distancia << endl;
    }
  }while(deltaMinimium < 0);  
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

int main(int argc, char** argv) {

  vector<int>  s = {1,1};
  vector<int> s_;  
  vector<int> listaDeCandidatos;
  
  double distance = 0, _distance = 0, distanceMin = 10000000000;

  /* Variavéis de controle da aleatorieadade*/
  srand (time(NULL));
 
  readData(argc, argv, &dimension, &matrizAdj);
  //printData();
  for(int iteration = 0; iteration < IMAX; iteration++){
    // Criando lista de cidades candidatas (indices) a partir da cidade 2
    for (int i = 2; i <= dimension; i++) 
      listaDeCandidatos.push_back(i);   

    //Criando subtour com as cidades candidatas
    for(int i = 0; i < tamanhoSubtourInicial; i++){
      int j = (rand()) % listaDeCandidatos.size();
      s.insert(s.begin()+1, listaDeCandidatos[j]);
      listaDeCandidatos.erase(listaDeCandidatos.begin() + j);
    }

    cout << "s: " ;
    imprimeCidade(s);
    cout << endl;
    //calculaDistancia(s);
    
    s = Construction(s, listaDeCandidatos);
    s_ = s;
    cout << "solucao inicial: ";
    imprimeCidade (s);
    cout << endl;
    distance = calculaDistancia(s);
    

    for(int iterILS = 0; iterILS < IILS; iterILS++){
      _distance =  Swap(s,distance); 
      if(_distance < distance){
        s_ = s;
        iterILS = 0;
      }
      s = Pertub(s_); 
    }
    

      cout << endl <<"solucao final rodada : "  << endl;    
      imprimeCidade (s);
      cout << endl;
      distance = calculaDistancia(s);   




  }

  cout << "Menor distancia = " << _distance << endl;
    return 0;  
}