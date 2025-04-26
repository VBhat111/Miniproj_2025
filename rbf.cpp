#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<vector>
#include<algorithm>
#include<eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;
double func_select(int choice, double xi, double xj, double yi, double yj,double r0){
    double rsq=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj);
   
    if(choice==1){
        return sqrt(rsq+(r0*r0));
    }
    else if(choice==2){
        return 1/sqrt(rsq+(r0*r0));
    }
    else if(choice==3){
        return (rsq==0)?0:rsq*log(sqrt(rsq)/r0);
    }
    else if(choice==4){
        return exp(-1*rsq)/(2*r0*r0);
    }
    else{
        cout<<endl<<"!!!!INVALID CHOICE!!!!"<<endl;
        return 0;
    }
}
int main(){
    int option,count=0;
    double x,y,f,d_avg,d_tot=0.0,r0,ans=0.0,xmin,xmax,ymin,ymax,error;
    vector<double> X,Y,F;
    cout<<"Choose RBF function:"<<endl<<"options:"<<endl<<"1)Multiquadratic\n2)Inverse Multiquadratic\n3)Thin plate spline\n4)Gaussian"<<endl;
    cin>>option;

    ifstream ip("add.txt"); //sample input filename
    if (!ip) {
        cout<<"Input file couldn't be opened!"<<endl;
        return 1;
    }
    ofstream op("add_rbf_mq.txt");//sample output filename when multiquadratic option is chosen

    while(ip>>x>>y>>f){
     X.push_back(x);
     Y.push_back(y);
     F.push_back(f);
    }
    for(int i=0;i<X.size();i++){
        for(int j=0;j<X.size();j++){
            d_tot+=sqrt(pow((X[i]-X[j]),2)+pow((Y[i]-Y[j]),2));
            count++;
        }
    }
    d_avg=d_tot/count;
    r0=d_avg/2;

    int n=X.size();
    MatrixXd A(n,n);
    VectorXd f_new(n);

    
    for(int i=0;i<X.size();i++){
        f_new(i)=F[i];
        for(int j=0;j<X.size();j++){
            A(i,j)= func_select(option,X[i],X[j],Y[i],Y[j],r0);
        }
    }

    VectorXd lambda= A.colPivHouseholderQr().solve(f_new);
    
    xmin=*min_element(X.begin(),X.end());
    xmax=*max_element(X.begin(),X.end());
    ymin=*min_element(Y.begin(),Y.end());
    ymax=*max_element(Y.begin(),Y.end());

    op<<fixed<<setprecision(10);
    op<<setw(20)<<"x"<<setw(20)<<"y"<<setw(20)<<"RBF Interpolated val"<<setw(20)<<"error(%)"<<endl;

    for(double xip=xmin;xip<=xmax;xip+=0.01){
        for(double yip=ymin;yip<=ymax;yip+=0.01){
            ans=0.0;
    for(int i=0;i<n;i++){
        ans+=lambda(i)*func_select(option,xip,X[i],yip,Y[i],r0);
    }
        error=abs(((5*xip+3*yip)-ans)/(5*xip+3*yip))*100;//sample function 5x+3y chosen
        op<<setw(20)<<xip<<setw(20)<<yip<<setw(20)<<ans<<setw(20)<<error<<endl;
        }
    }
    cout<<"Interpolation completed";
  return 0;  
}