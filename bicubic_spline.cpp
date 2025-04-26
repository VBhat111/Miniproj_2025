#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<iomanip>
#include<algorithm>
using namespace std;



struct CubicSpline {
    vector<double> x, f, b, c, d;
    int n;

    void init(vector<double>& X, vector<double>& F) {
        x=X;
        f=F;
        n=x.size()-1;
        vector<double> h(n),alpha(n+1),l(n+1),mu(n+1),z(n+1);
        b.resize(n+1); c.resize(n+1); d.resize(n+1);

        for(int i=0;i<n;++i)
            h[i]=x[i+1]-x[i];

        for(int i=1;i<n;++i)
            alpha[i]=(3/h[i])*(f[i+1]-f[i])-(3/h[i-1])*(f[i]-f[i-1]);

        l[0]=1; 
        mu[0]=z[0]=0;

        for(int i=1;i<n;++i) {
            l[i]=2*(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
            mu[i]=h[i]/l[i];
            z[i]=(alpha[i]-h[i-1]*z[i-1])/l[i];
        }

        l[n]=1; 
        z[n]=c[n]=0;

        for(int j=n-1;j>=0;--j) {
            c[j]=z[j]-mu[j]*c[j+1];
            b[j]=(f[j+1]-f[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
            d[j]=(c[j+1]-c[j])/(3*h[j]);
        }
    }

    double eval(double xi) const {
        int i=lower_bound(x.begin(), x.end(), xi) - x.begin()-1;
        if(i<0)i=0;
        if(i>=n)i=n-1;
        double dx = xi-x[i];
        return f[i]+b[i]*dx+c[i]*dx*dx+d[i]*dx*dx*dx;
    }
};

int main() {
    ifstream ip("sincos.txt");//sample input filename
    if (!ip) {
        cout<<"Input file couldn't be opened!"<<endl;
        return 1;
    }
    ofstream op("sincos_bicubicspline.txt");//sample output filename

    vector<double> x_vals, y_vals;
    vector<vector<double>> Z;

    double x, y, z;
    int x_count=0, y_count=0;
    double prev_x=-1.0;

  
    while(ip>>x>>y>>z) {
        if(x!=prev_x) {
            x_vals.push_back(x);
            Z.push_back(vector<double>());
            prev_x=x;
            x_count++;
        }
        if(x_count==1) y_vals.push_back(y);
        Z.back().push_back(z);
    }
    y_count=y_vals.size();

    op<<setw(10)<<"X"<<setw(10)<<"Y"<<setw(20)<<"Bicubic spline"<<setw(20)<<"error(%)"<<endl;

    double step=0.05;
    for(double xi=x_vals.front();xi<=x_vals.back();xi+=step) {
        
        vector<double> temp_y_interp;
        for(int i=0;i<x_vals.size();++i) {
            CubicSpline spline_y;
            spline_y.init(y_vals, Z[i]);
            temp_y_interp.push_back(spline_y.eval(xi));
        }

        
        CubicSpline spline_x;
        spline_x.init(x_vals,temp_y_interp);
        for (double yi=y_vals.front();yi<=y_vals.back();yi+=step) {
            double interp_val=spline_x.eval(yi);
            double error= abs((sin(xi)+cos(yi)-interp_val)/(sin(xi)+cos(yi)))*100;//sample func sinx+cosy used
            op<<fixed<<setprecision(2)<<setw(10)<<xi<<setw(10)<<yi<<fixed<<setprecision(6)<<setw(20)<<interp_val<<setw(20)<<error<<endl; 
               
        }
    }

    cout<<"Interpolationn completed!"<<endl;
    return 0;
}
