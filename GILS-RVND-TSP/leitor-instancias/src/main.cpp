#include "readData.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <time.h>      
#include <map>   
#include <random>

using namespace std;

//double ** matrizAdj;
int dimension = 6; // quantidade total de vertices
double alpha = 0.5;

int matrizAdj[7][7] =  {
                      {0,0,0,0,0,0,0},
                      {0,0,2,1,4,9,1},
                      {0,2,0,5,9,7,2},
                      {0,1,5,0,3,8,6},
                      {0,4,9,3,0,2,5},
                      {0,9,7,8,2,0,2},
                      {0,1,2,6,5,2,0}
    }; // matriz de adjacencia



struct InsertionInfo{
  int noInserido; // no k a s e r i n s e r i d o
  int arestaRemovida; // a r e s t a { i , j } onde o no k s e r a i n s e r i d o
  double custo; // d e l t a ao i n s e r i r k na a r e s t a { i , j }

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

int calculaDistancia(vector<int> listaCidade){
  int distancia = 0;

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


int main(int argc, char** argv) {

  vector<int>  s = {1,1};  
  vector<int> listaDeCandidatos;
  int tamanhoSubtourInicial = 2;

  /* Variavéis de controle da aleatorieadade*/
  srand (time(NULL));
  std::default_random_engine generator (time(NULL));
  std::discrete_distribution<int> distribution {1-alpha,alpha};
  int choosed, number;

 
  //readData(argc, argv, &dimension, &matrizAdj);
  printData();

  // Criando lista de cidades candidatas (indices) a partir da cidade 2
  for (int i = 2; i <= dimension; i++) 
    listaDeCandidatos.push_back(i); 

  //Criando subtour com as cidades candidatas
  for(int i = 0; i < tamanhoSubtourInicial; i++){
    int j = (0*rand()) % listaDeCandidatos.size();
    s.insert(s.begin()+1, listaDeCandidatos[j]);
    listaDeCandidatos.erase(listaDeCandidatos.begin() + j);
  }

  cout << "s: " ;
  imprimeCidade(s);
  cout << endl;
  //calculaDistancia(s);


  vector <InsertionInfo> custoInsercao((s.size()-1)*listaDeCandidatos.size());

  while(listaDeCandidatos.size()){
    for(int i = 1, j = i+1, l = 0 ; i <= s.size() - 1; i++, j++){
      for(auto k : listaDeCandidatos){
        custoInsercao[l].custo = matrizAdj[s[i-1]][k] + matrizAdj[s[j-1]][k] - matrizAdj[s[i-1]][s[j-1]];
        custoInsercao[l].noInserido = k;
        custoInsercao[l].arestaRemovida = i ;
        l++;
      }
    }

    //for(int i = 0; i < custoInsercao.size(); i++)
    //  cout << custoInsercao.at(i).noInserido << ' ' <<  custoInsercao.at(i).arestaRemovida << ' ' << custoInsercao.at(i).custo << endl;
    std::sort(custoInsercao.begin(), custoInsercao.end(), InsertionInfo());

    number = distribution(generator);
    if(!number) choosed = 0;
    else choosed = rand()%custoInsercao.size();
    
    s.insert(s.begin()+custoInsercao[choosed].arestaRemovida, custoInsercao[choosed].noInserido);
    listaDeCandidatos.erase(remove(listaDeCandidatos.begin(), listaDeCandidatos.end(), custoInsercao[choosed].noInserido), listaDeCandidatos.end());
    custoInsercao.resize((s.size()-1)*listaDeCandidatos.size());
  }

  

  cout << "s: " ;
  imprimeCidade(s);
  cout << endl;
  calculaDistancia(s);

    
    return 0;  
}