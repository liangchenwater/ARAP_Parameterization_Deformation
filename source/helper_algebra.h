//
//  helper_algebra.h
//  arap
//
//  Created by Liang Chen on 2021/6/14.
//

#ifndef helper_algebra_h

#define helper_algebra_h
#include "global_var.h"

MatrixXd getCovariance3x3(const vector<int>& neighbor,const MatrixXd& verts,const MatrixXd& now_res,const VectorXd& weights,int cur);

MatrixXd getJacobian2x2(const vector<HalfEdge>& half_edges,const MatrixXd& now_res,int cur);

MatrixXd getCovariance2x2(const vector<HalfEdge>& half_edges,const MatrixXd& now_res,const VectorXd& weights,int cur);

void newton_optimizerNto1(VectorXd& x,std::function<VectorXd(const VectorXd&)> compute_J,std::function<MatrixXd(const VectorXd&)> compute_H);

void binary_find_root(double& x, std::function<double(double)>);

class R_Jacobian{
public:
    VectorXd operator()(const VectorXd& x)
    {
        VectorXd J(x.rows());
        J(0)=4*lamda*x(0)*(x(0)*x(0)+x(1)*x(1)-1)+2*C1*x(0)-2*C2;
        J(1)=4*lamda*x(1)*(x(0)*x(0)+x(1)*x(1)-1)+2*C1*x(1)-2*C3;
        return J;
    }
    R_Jacobian(double lamda,double C1,double C2,double C3):lamda(lamda),C1(C1),C2(C2),C3(C3){}
private:
    double C1;
    double C2;
    double C3;
    double lamda;
};

class R_Hessian{
public:
    MatrixXd operator()(const VectorXd& x)
    {
        MatrixXd H(x.rows(),x.rows());
        H(0,0)=4*lamda*(3*x(0)*x(0)+x(1)*x(1)-1)+2*C1;
        H(1,1)=4*lamda*(3*x(1)*x(1)+x(0)*x(0)-1)+2*C1;
        H(0,1)=H(1,0)=8*lamda*x(0)*x(1);
        return H;
    }
    R_Hessian(double lamda,double C1):lamda(lamda),C1(C1){}
private:
    double C1;
    double lamda;
};

#endif /* helper_algebra_h */
