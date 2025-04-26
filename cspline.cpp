#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<vector>
#include<algorithm>
using namespace std;


void finalans(vector<double>& X,vector<double>& F, vector<double>& b,vector<double>& c,vector<double>& d, ofstream &op){
int ind;
double xdif,Sx,lb,step=0.01,p_error;
for(double i=X[0]+step;i<X[X.size()-1];i+=step){ 
    ind=(lower_bound(X.begin(), X.end(), i))-X.begin();
   lb=X[ind-1];
   xdif=i-lb;
   Sx = F[ind-1] + b[ind-1]*xdif + c[ind-1]*xdif*xdif + d[ind-1]*xdif*xdif*xdif;
   p_error=abs(((0.9*i*i*i-3*i*i+9*i+5)-Sx)/(0.9*i*i*i-3*i*i+9*i+5))*100; //sample function used is 0.9x^3-3x^2+9x+5
  op<<fixed<<setprecision(6)<<setw(15)<<i;
  op<<fixed<<setprecision(10)<<setw(20)<<Sx<<setw(20)<<p_error<<endl;

}
}
void matrixB(int n, vector<double>& B, vector<double>& F,vector<double>& h){
    for(int i=1;i<n;i++){
        B[i]=3*((F[i+1]-F[i])/h[i]-(F[i]-F[i-1])/h[i-1]);
      }
}

void solveC(int n, vector<vector<double>>& A, vector<double>& B,vector<double>& b,vector<double>& c,vector<double>& d,vector<double>& F,vector<double>& h){
    double m;
    for(int i=1;i<n;i++){
        m=1/(A[i][i]-(A[i][i-1]*A[i-1][i]));
        A[i][i+1]=A[i][i+1]*m;
        B[i]=m*(B[i]-A[i][i-1]*B[i-1]);
        }        
        c[n]=B[n];
        for(int i=n-1;i>=0;i--){
            c[i]=B[i]-(c[i+1]*A[i][i+1]);
            b[i]=(1.0/h[i])*(F[i+1]-F[i])-(h[i]/3)*(2*c[i]+c[i+1]);
            d[i]=(c[i+1]-c[i])/(3.0*h[i]);
        
        }
}

void resizing(vector<vector<double>>& A, vector<double>& B,vector<double>& b,vector<double>& c,vector<double>& d,vector<double>& h,int n){
A.resize(n+1,vector<double>(n+1,0.0));
B.resize(n+1,0.0);
c.resize(n+1,0.0);
b.resize(n+1,0.0);
d.resize(n+1,0.0);
h.resize(n,0.0);
}

void matrixA(int n,vector<vector<double>>& A, vector<double>& h){
    for(int i=0;i<=n;i++){
        for(int j=0;j<=n;j++){
            if(i==j && i!=0 && i!=n) A[i][j]=2*(h[i-1]+h[i]);
            else if(j==i-1) A[i][j]=h[j];
            else if(j==i+1) A[i][j]=h[i];
            else A[i][j]=0;
        }
    }
}

int main(){
    ifstream ip("cub2.txt");//sample inpu file name
if (!ip) {
    cout << "Input file couldn't be opened!" << endl;
    return 1;
}
ofstream op("cub_cs2.txt");//sample output file name

vector<double> X,F,c,B,b,d,h;//Ac=b
vector<vector<double>> A;
double x,f;

while(ip>>x>>f){
    X.push_back(x);
    F.push_back(f);
}
int n=X.size()-1;
resizing(A,B,b,c,d,h,n);
for(int i=0;i<n;i++){
    h[i]=X[i+1]-X[i];
}
op<<"\nnatural spline"<<endl;
op<<setw(15)<<"x "<<setw(20)<<"Natural CS "<<setw(20)<<"N_Error(%) "<<endl;
matrixA(n,A,h);
for(int j=0;j<=n;j++){
    A[0][j]=0.0;
    A[n][j]=0.0;
}
A[0][0]=1.0;
A[n][n]=1.0;

matrixB(n,B,F,h);
B[0]=0;B[n]=0.0;

solveC(n,A,B,b,c,d,F,h);
finalans(X,F,b,c,d,op);

resizing(A,B,b,c,d,h,n);
op<<"\nClamped spline"<<endl;
op<<setw(15)<<"x "<<setw(20)<<"Clamped S "<<setw(20)<<"Cd_Error(%) "<<endl;

matrixA(n,A,h);
A[0][0]=2*h[0];
A[n][n]=2*h[n-1];

matrixB(n,B,F,h);
double F1x=(F[1]-F[0])/(X[1]-X[0]); 
double Fnx=(F[n]-F[n-1])/(X[n]-X[n-1]); 
B[0]=3*(((F[1]-F[0])/h[0])-F1x);
B[n]=3*((-1*(F[n]-F[n-1])/h[n-1])+Fnx);

A[0][1]=A[0][1]/A[0][0];
B[0]=B[0]/A[0][0];
solveC(n,A,B,b,c,d,F,h);
finalans(X,F,b,c,d,op);

cout<<"Interpolation completeddd";
}