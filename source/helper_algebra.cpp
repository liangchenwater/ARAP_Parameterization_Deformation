//
//  helper_algebra.cpp
//  arap
//
//  Created by Liang Chen on 2021/6/15.
//

#include "helper_algebra.h"
//it's not a general solver, but for the cubic function in the local step of Hybrid method, it surely works
void binary_find_root(double& x, std::function<double(double)> fun)
{
    x=0;
    double y=fun(x);
    double left=0,right=0;
    if(y<-eps){
        left=x;
        right=inf_for_bfr;
    }
    else if(y>eps){
        right=x;
        left=-inf_for_bfr;
    }
    else return;
    while(right-left>eps){
         x=(left+right)/2;
         y=fun(x);
        if(y<-eps) left=x;
        else if(y>eps) right=x;
        else return;
    }
}

void newton_optimizerNto1(VectorXd& x,std::function<VectorXd(const VectorXd&)> compute_J,std::function<MatrixXd(const VectorXd&)> compute_H)
{
    VectorXd prevx(2);
    prevx.setZero();
    int i=0;
    while(i<MAX_NEWTON_ITRS&&(x-prevx).norm()/prevx.norm()>eps){
        VectorXd J=compute_J(x);
        MatrixXd H=compute_H(x);
       // cout<<H.determinant()<<endl;
        prevx=x;
        assert(H.fullPivLu().isInvertible());
        x-=newton_step_size*H.fullPivLu().inverse()*J;
        i++;
    }
}

MatrixXd getCovariance3x3(const vector<int>& neighbor,const MatrixXd& verts,const MatrixXd& now_res,const VectorXd& weights,int cur)
{
    MatrixXd S(3,3);
    S.setZero();
    for(int i=0;i<neighbor.size();i++){
        double w=weights((long long unsigned)cur*verts.cols()+neighbor[i]);
        S+=w*(now_res.row(cur)-now_res.row(neighbor[i])).transpose()* (verts.col(cur)-verts.col(neighbor[i])).transpose();
    }
    return S;
}

MatrixXd getJacobian2x2(const vector<HalfEdge>& half_edges,const MatrixXd& now_res,int cur)
{
    MatrixXd JT(2,2);
    MatrixXd Origin(2,2);
    MatrixXd RHS(2,2);
    int a=half_edges[cur*3].Endpoints(0);
    int b=half_edges[cur*3].Endpoints(1);
    int c=half_edges[cur*3+1].Endpoints(0);
    int d=half_edges[cur*3+1].Endpoints(1);
    Origin.row(0) = half_edges[cur*3].EdgeVec.transpose();
    Origin.row(1)=  half_edges[cur*3+1].EdgeVec.transpose();
    RHS.row(0)=now_res.row(a)-now_res.row(b);
    RHS.row(1)=now_res.row(c)-now_res.row(d);
    JT=Origin.colPivHouseholderQr().solve(RHS);
    return JT.transpose();
}

MatrixXd getCovariance2x2(const vector<HalfEdge>& half_edges,const MatrixXd& now_res,const VectorXd& weights,int cur)
{
    MatrixXd S(2,2);
    S.fill(0);
    for(int i=0;i<3;i++){
        int a=half_edges[cur*3+i].Endpoints(0);
        int b=half_edges[cur*3+i].Endpoints(1);
        Vector2d U=(now_res.row(a)-now_res.row(b)).transpose();
        Vector2d X=half_edges[cur*3+i].EdgeVec;
        S+=weights(cur*3+i)*U*(X.transpose());
    }
    return S;
}
