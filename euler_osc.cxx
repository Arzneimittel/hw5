#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

void analytic( const double dt1,const double dt2,  const double N1, const double N2,double* a1,double* a2){
   for(int i=0;i<N1;i++){
     a1[i] = cos(i * dt1);
  }   
   for(int i=0;i<N2;i++){
     a2[i] = cos(i * dt2); 
  } 
}

void write( const double N1,const double N2,double* a1,double* a2,double* f1,double* f2,double* b1,double* b2,const double dt1,const double dt2){
    ofstream out1 ("Euler_OS_dt1.txt");
  for(int i=0; i<N1; i++){
    out1 << dt1 *i << "\t" <<a1[i] <<"\t"<< f1[i] <<"\t"<< b1[i]<<endl;
  }
    ofstream out2 ("Euler_OS_dt2.txt");
  for(int i=0; i<N2; i++){
    out2 << dt2 *i << "\t" <<a2[i] <<"\t"<< f2[i] <<"\t"<< b2[i]<<endl;
  }
}
void forward(const double dt1,const double dt2,const int N1,const int N2,double* f1, double* f2){
   double x[2];
   x[0] = 1;
   x[1] = 0;
   f1[0] = x[0];
   f2[0] = x[0];
   for(int i=1; i<N1; i++){
     double temp = x[0];
     x[0] = x[0] + x[1] * dt1; 
     x[1] = x[1] - temp * dt1;
     f1[i] = x[0];
  }
   x[0] = 1;
   x[1] = 0;
   for(int i=1; i<N2; i++){
     double temp = x[0];
     x[0] = x[0] + x[1] * dt2; 
     x[1] = x[1] - temp * dt2;
     f2[i] = x[0];
  } 
}

void backward(const double dt1,const double dt2,const int N1,const int N2,double* b1, double* b2){
  double x[2];
  x[0] = 1;
  x[1] = 0;
  b1[0] = x[0];
  b2[0] = x[0];
    for(int i=1; i<N1; i++){
     double temp = x[0];
     x[0] = x[0] + x[1] * dt1; 
     x[0] /= (1+ (dt1*dt1));
     x[1] = x[1] - temp * dt1;
     x[1] /= (1+ (dt1*dt1));
     b1[i] = x[0];
}
    x[0] = 1;
    x[1] = 0;
    for(int i=1; i<N2; i++){
     double temp = x[0];
     x[0] = x[0] + x[1] * dt2; 
     x[0] /= (1+ (dt2*dt2));
     x[1] = x[1] - temp * dt2;
     x[1] /= (1+ (dt2*dt2));
     b2[i] = x[0];
}
}

int main(void){
  const double dt1 = M_PI / 10;
  const double dt2 = M_PI / 100;
  const int N1 = (20 * M_PI) / dt1;
  const int N2 = (20 * M_PI) / dt2;
  double a1[N1],f1[N1],b1[N1];
  double a2[N2],f2[N2],b2[N2];
  analytic(dt1,dt2,N1,N2,a1,a2);
  forward(dt1,dt2,N1,N2,f1,f2);
  backward(dt1,dt2,N1,N2,b1,b2);
  write(N1,N2,a1,a2,f1,f2,b1,b2,dt1,dt2);
  
  return 0;
}
