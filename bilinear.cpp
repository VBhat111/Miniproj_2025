#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<map>
#include<algorithm>
#include<fstream>
using namespace std;

double p_error(double a, double i,double j){
return abs((sin(i)+cos(j)-a)/(sin(i)+cos(j)))*100.0;  //sample function sinx+cosy used
}

int main(){
int x1,x2,y1,y2,a,b;
double h=0.01,ans,x,y,f,L1,L2;
vector<double> X, Y;
map<pair<double, double>, double> temp_F;
ifstream ip("sincos.txt");  //sample input filename
if (!ip) {
    cout << "Input file couldn't be opened!" << endl;
    return 1;
}
ofstream op("sincos_bilinear.txt");//sample output filename
op<<setw(15)<<"x "<<setw(15)<<"y "<<setw(15)<<"BiLinear "<<setw(15)<<"Error(%) "<<endl;
   while (ip >> x >> y >> f) {
        temp_F[{x, y}] = f;  
    }

    
    for (const auto& entry : temp_F) {
        if (find(X.begin(), X.end(), entry.first.first) == X.end()) {
            X.push_back(entry.first.first);  // Add X if not already present
        }
        if (find(Y.begin(), Y.end(), entry.first.second) == Y.end()) {
            Y.push_back(entry.first.second); // Add Y if not already present
        }
    }


    vector<vector<double>> F(X.size(), vector<double>(Y.size(), 0));

    
    for (size_t i = 0; i < X.size(); ++i) {
        for (size_t j = 0; j < Y.size(); ++j) {
            pair<double, double> key = {X[i], Y[j]};
            if (temp_F.find(key) != temp_F.end()) {
                F[i][j] = temp_F[key];
            }
        }
    }
for(double i=X[0]+h;i<X[X.size()-1];i+=h){
 for(double j=Y[0]+h;j<Y[Y.size()-1];j+=h){
      x2=(lower_bound(X.begin(), X.end(), i))-X.begin();
      x1=x2-1;
      y2=(lower_bound(Y.begin(), Y.end(), j))-Y.begin();
      y1=y2-1;
      op<<fixed<<setprecision(2)<<setw(15)<<i<<setw(15)<<j;
      op<<fixed<<setprecision(6);
      L1=((i-X[x2])/(X[x1]-X[x2]))*F[x1][y1]+((i-X[x1])/(X[x2]-X[x1]))*F[x2][y1]    ;
      L2=((i-X[x2])/(X[x1]-X[x2]))*F[x1][y2]+((i-X[x1])/(X[x2]-X[x1]))*F[x2][y2]    ;
      ans=L1*((j-Y[y2])/(Y[y1]-Y[y2]))+L2*((j-Y[y1])/(Y[y2]-Y[y1]));
      op<<setw(15)<<ans<<setw(15)<<p_error(ans,i,j);
op<<endl;
}
op<<endl;
}
cout<<"Interpolation completed";

return 0;}
