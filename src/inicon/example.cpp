#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <vector>
#include "splitandfit.h"
#include "trial-functions.h"

using namespace std;

/*
 * type = 0 : D dimensional cubic domain
 * type = 1 : 2 dimensional triangular domain
 * type = 2 : 3 dimensional tetraedric domain
 */
 
/*
 * Exemplo do uso da Biblioteca de particionamento de domínio.
 * 
 * Os inputs são:
 * 1)Uma função, a ser usada no molde da gsl, 
 *   como as presentes no arquivo "trial-functions.h" ou abaixo tkgauss
 *   Esta função seria posta como membro de uma gsl_monte_function, 
 *   que tem que estar declarada junto com o código.
 * 
 * 2)Uma variavel do tipo wparams, que servirá para compatibilizar a 
 *   função com o domínio de integração. Caso haja algum parâmetro
 *   setável dentro da função, utilize o membro p de wparams para
 *   armazenar os parâmetros, e não o membro params da 
 *   gsl_monte_function, este deve ser reservado para a wparams
 * 
 * 3)Um vetor de domínios dom, o qual será inicializado para ser um 
 *   domínio ou de cubos ou de triangulos, no presente momento
 * 
 * 4)Um valor de cutoff que estabelecará o maior peso (i.e. integral 
 *   da função dada) que um domínio pode ter. Se um dado domínio 
 *   passar deste valor de peso, ele é particionado pelo código.
 * 
 * Todas as funções a serem executadas retornam um int, que é 0 se 
 * a função foi corretamente executada, e 1 caso contrário. Os erros
 * ainda não foram catalogados, mas 0 é sucesso de execução, caso
 * não seja 0, procure no código da função o que representa o erro.
 * 
 * A inicialização por cubos exigem os intervalos aonde está definido
 * O cubo, o que quer dizer:
 * 
 *       xu[1]  _______________________
 *             |                       |
 *             |                       |
 *             |                       |
 *             |                       |
 *             |                       |
 *       xl[1] |_______________________|
 *            xl[0]                   xu[0]
 * 
 * A inicialização por triangulos exige um array de triangulos, com
 * as posições de cada triangulo.
 * 
 * Um dado triangulo é caracterizado por
 * 
 *                 (x[2],x[3])
 *                      /\
 *                     /  \
 *                    /    \
 *                   /      \
 *                  /        \
 *                 /          \
 *                /            \
 *               /______________\
 *          (x[0],x[1])     (x[4],x[5])
 * 
 *  
 * Quando há vários triângulos, numere cada triangulo, e armazene
 * sequencialmente os vertices de cada um deles separadamente, mesmo
 * que haja vértices em comum, isto não é um problema.
 *                  _______________  ______________
 *                 /\              /\              /\
 *                /  \            /  \            /  \
 *               /    \    1     /    \    3     /    \
 *              /      \        /      \        /      \
 *             /   0    \      /   2    \      /    4   \
 *            /          \    /          \    /          \
 *           /            \  /            \  /            \
 *          /______________\/______________\/______________\
 * 
 * Portanto, o elemento do array que dá, em d dimensões 
 * (i.e., há d+1 vertices), a coordenada l do vertice v do triangulo t
 * é:
 * 
 *  x[t*(d*(d+1)) + d*v + l]
 *  
 * Exemplos estão feitos na função unit_hexagon e unit2_hexagon 
 * 
 * Depois disto, as funções domain_split, clean_domain e print_sph
 * devem ser chamadas nesta ordem.
 * 
 * Se clean_domain for chamado antes de domain_split, ela apenas não
 * terá efeito algum, mas se print_sph for chamado antes de clean_domain
 * , ela imprimirá apenas os domínios de input, enquanto se ela for
 * chamada antes de clean_domain mas depois de domain_split, ela irá
 * imprimir partículas SPH vindas de domínios espúrios.
 * 
 * Numa revisão posterior, clean_domain deverá ser imbutida dentro de 
 * domain_split
 */
 
double tkgauss(double x[],size_t dim,void *par){
  int l,err;
  double X[dim],Xmax[dim],funct,dist=0.;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
    
  Xmax[0]=3.0;
  if(dim>=2){
    Xmax[1]=2.0;
    if(dim>=3){
      Xmax[2]=1.0;
      for(l=3;l<dim;l+=1)
        Xmax[l]=0.0;
    }	
  } 
  
  for(l=0;l<dim;l+=1)
    dist+= (x[l]/Xmax[l])*(x[l]/Xmax[l]);
  funct=exp(-15.0*dist);
  
  return funct;
}

int main(int argc,char **argv){
  const int D=2,Ntri=6,split_type=0;
  int err,l; 
  double cutoff=0.001,xi[D],xf[D],xv[Ntri*(D+1)*D];
  wparams par;
  vector <domain> dom;
  gsl_monte_function F;
  
  F.f = &tkgauss; F.dim=D;par.p=NULL;F.params=&par;
  
  if(split_type==0){
    cout << "init\n";
    for(l=0;l<D;l+=1){xi[l]=-1.;xf[l]=1.;}  
    err=init_cube(xi,xf,dom,D);if(err!=0) return err; 
  }
  else if(split_type==1){
    cout << "unit\n";
    err=unit2_hexagon(Ntri,D,xv);if(err!=0) return err;
    
    cout << "init\n";
    err=init_triangle(D,Ntri,xv,dom);if(err!=0) return err;
  } 
  
  cout << "split\n";
  err=domain_split(D,cutoff,dom,F); if(err!=0){ cout << "out: " << err << endl;return err;}
  
  cout << "clean\n";
  err=clean_domain(dom);if(err!=0) return err;
  
  cout << "print\n";  
  err=print_sph(D,"SPH-particles.dat",dom); if(err!=0) return err;
  
  return 0;
}
