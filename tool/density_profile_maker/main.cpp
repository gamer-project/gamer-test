//This program is an example for producing a density profile
//The scenario is an NFW model with scaling radius r0 = 0.01, and the peak density is rho0 = 1.0
#include<stdlib.h>
#include<iostream>
#include <fstream>
using namespace std;

double nfw(double x){
  return 1/(x*(1+x)*(1+x));
}
int main()
{
  double rho0 = 1.0;
  double r0 = 0.01;
  fstream file;
  file.open("profile.txt",ios::out);
  file<<"radius"<<'\t'<<"density"<<endl;
  double dx = 20./1000;
  for(int k=1;k<1001;k++){
    double x= dx*k;
    file<<x*r0<<'\t'<<rho0*nfw(x)<<endl;
    //In this case, radius is at 0-th column, and density is at 1-th column,
    //and you should input Models_r_col as 0, and Models_rho_col as 1 in your Init_TestProb file
  }
  file.close();
  
  

}