#include "readData.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>

using namespace std;

double ** matrizAdj;
int dimension ; // quantidade total de vertices
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
 
  readData(argc, argv, &dimension, &matrizAdj);
  //printData();

  
  // Criando lista de cidades candidatas (indices) a partir da cidade 2
  for (int i = 2; i <= dimension; i++) 
    listaDeCandidatos.push_back(i); 

  

  //Criando subtour com as cidades candidatas
  for(int i = 0; i < tamanhoSubtourInicial; i++){
    int j = rand() % listaDeCandidatos.size();
    s.insert(s.begin()+1, listaDeCandidatos[j]);
    listaDeCandidatos.erase(listaDeCandidatos.begin() + j);
  }

/*  
  s.insert(s.begin()+1, 3);
  listaDeCandidatos.erase(remove(listaDeCandidatos.begin(), listaDeCandidatos.end(), 3), listaDeCandidatos.end());
  s.insert(s.begin()+2, 2);
  listaDeCandidatos.erase(remove(listaDeCandidatos.begin(), listaDeCandidatos.end(), 2), listaDeCandidatos.end());
 */ 

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
    s.insert(s.begin()+custoInsercao[0].arestaRemovida, custoInsercao[0].noInserido);
    listaDeCandidatos.erase(remove(listaDeCandidatos.begin(), listaDeCandidatos.end(), custoInsercao[0].noInserido), listaDeCandidatos.end());
    custoInsercao.resize((s.size()-1)*listaDeCandidatos.size());
  }

  cout << "s: " ;
  imprimeCidade(s);
  cout << endl;
  calculaDistancia(s);

    
    return 0;  
}