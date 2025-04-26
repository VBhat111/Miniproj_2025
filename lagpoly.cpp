#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<vector>
#include<algorithm>
using namespace std;

double p_error(double a, double i){
return abs((3*i*i+9*i+5-a)/(3*i*i+9*i+5))*100.0;//sample function 3x^2+9x+5 used
}
int main(){

int near1,near2,near3,near4;
float h=0.01,x,f,newX1,newX2a,newX2b,newX3,L1,L2,L3,L4;
ifstream ip("quad.txt");//sample input filename
if (!ip) {
    cout << "Input file couldn't be opened!" << endl;
    return 1;
}
ofstream op("quad_lp.txt");//sample output filename
op<<setw(15)<<"x "<<setw(15)<<"Linear LP "<<setw(15)<<"L_Error(%) "<<setw(15)<<"Quadratic LP "<<setw(15)<<"Q_Error(%) "<<setw(15)<<"Cubic LP "<<setw(15)<<"C_Error(%) "<<endl;
vector<double> X,F;
while(ip>>x>>f){
    X.push_back(x);
    F.push_back(f);
}
for(double i=X[0]+h;i<X[X.size()-1];i+=h){
    near2=(lower_bound(X.begin(), X.end(), i))-X.begin();
    
    near1=near2-1;
    near3=near2-2;
    near4=near2+1;
    
    

    op<<fixed<<setprecision(2)<<setw(15)<<i;
    op<<fixed<<setprecision(6);
    L1=(i-X[near2])/(X[near1]-X[near2]);
    L2=(i-X[near1])/(X[near2]-X[near1]);
    newX1=F[near1]*L1+F[near2]*L2;
    op<<setw(15)<<newX1<<setw(15)<<p_error(newX1,i);
    if(X.size()==2) { op<<endl; continue;}
    bool check=false;
    if(near3>=0){
        L1=((i-X[near1])*(i-X[near2]))/((X[near3]-X[near1])*(X[near3]-X[near2]));
        L2=((i-X[near3])*(i-X[near2]))/((X[near1]-X[near3])*(X[near1]-X[near2]));
        L3=((i-X[near3])*(i-X[near1]))/((X[near2]-X[near3])*(X[near2]-X[near1]));
        newX2a= F[near3]*L1+F[near1]*L2+F[near2]*L3;
        
    }
    else{
        check=true;
    }
    
    if(near4<X.size()){
        L1=((i-X[near2])*(i-X[near4]))/((X[near1]-X[near2])*(X[near1]-X[near4]));
        L2=((i-X[near1])*(i-X[near4]))/((X[near2]-X[near1])*(X[near2]-X[near4]));
        L3=((i-X[near1])*(i-X[near2]))/((X[near4]-X[near1])*(X[near4]-X[near2]));
        newX2b= F[near1]*L1+F[near2]*L2+F[near4]*L3;
       if(check==false) {
    if(p_error(newX2a,i)<p_error(newX2b,i)){
        op<<setw(15)<<newX2a<<setw(15)<<p_error(newX2a,i);
    }
    else{
        op<<setw(15)<<newX2b<<setw(15)<<p_error(newX2b,i);
    }}
    else{op<<setw(15)<<newX2b<<setw(15)<<p_error(newX2b,i);}
}
 else{
    op<<setw(15)<<newX2a<<setw(15)<<p_error(newX2a,i);
}
if(X.size()==3) { op<<endl; continue;}

 if(near4>=X.size()){
    near3--,near2--,near1--,near4--;
}
else if(near3<0){
    near3++;near1++;near2++;near4++;
}

 
        L1=((i-X[near1])*(i-X[near2])*(i-X[near4]))/((X[near3]-X[near1])*(X[near3]-X[near2])*(X[near3]-X[near4]));
        L2=((i-X[near3])*(i-X[near2])*(i-X[near4]))/((X[near1]-X[near3])*(X[near1]-X[near2])*(X[near1]-X[near4]));
        L3=((i-X[near3])*(i-X[near1])*(i-X[near4]))/((X[near2]-X[near3])*(X[near2]-X[near1])*(X[near2]-X[near4]));
        L4=((i-X[near3])*(i-X[near1])*(i-X[near2]))/((X[near4]-X[near3])*(X[near4]-X[near1])*(X[near4]-X[near2]));
        newX3= F[near3]*L1+F[near1]*L2+F[near2]*L3+F[near4]*L4;
         op<<setw(15)<<newX3<<setw(15)<<p_error(newX3,i);

    op<<endl;

}
cout<<"Interpolationn completed";
return 0;
}
