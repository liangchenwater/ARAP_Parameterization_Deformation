//
//  ARAP.h
//  arap
//
//  Created by Liang Chen on 2021/6/14.
//

#ifndef ARAP_h
#define ARAP_h

#include "helper_algebra.h"
#include "helper_geometry.h"
void getLaplace(const vector<HalfEdge>& half_edges,const VectorXd& weights,const vector<int>& fix, SparseMatrix<double>& Laplace,int func);

void getRHS(MatrixXd* R,const vector<HalfEdge>& half_edges,const VectorXd& weights,const vector<int>& fix, const vector<VectorXd>& fix_vec,MatrixXd& RHS,int func,bool is_first);

void local_phase_deform(MatrixXd* R,vector<int>* neighbors,const MatrixXd& verts,const MatrixXd& res, const VectorXd& weights);

void local_phase_param(MatrixXd* R,const vector<HalfEdge>& half_edges,const MatrixXd& res,const VectorXd& weights,const VectorXd& area,int method,double lamda,double * distortion_per_unit,double* aread_per_unit,double *angled_per_unit,double& distortion,double& aread,double& angled);

void global_phase (MatrixXd* R,const vector<HalfEdge>& half_edges, const VectorXd& weights, const vector<int>& fix,const vector<VectorXd>& fix_vec,MatrixXd& RHS,int func,bool first, MatrixXd& new_res,SimplicialLDLT<SparseMatrix<double> >& dir_solver);

class ARAP_energy{
public:
    double operator()(MatrixXd& res){
        double E=0;
        MatrixXd L(2,2);
        for(int i=0;i<half_edges.size()/3;i++){
            //get cross-covariance matrix or jacobi matrix from (u,v) to (x,y)
            MatrixXd J=getJacobian2x2(half_edges,res,i);
    //            assert(J.determinant()<0);
            // J=U*Sigma*V^T
           JacobiSVD<MatrixXd> SVD_solver;
           if(method==ARAP||method==ASAP) SVD_solver.compute(J,ComputeThinU | ComputeThinV);
           else if(method==Hybrid) SVD_solver.compute(J);
            // R=U*V^T
            Vector2d singular_value=SVD_solver.singularValues();
            MatrixXd SI(2,2);
            SI.fill(0);
            SI(0,0)=SI(1,1)=0.5*(singular_value(0)+singular_value(1));
            if(method==ARAP||method==ASAP){
            if(method==ARAP) L=SVD_solver.matrixU()*SVD_solver.matrixV().transpose();
            else if(method==ASAP) L=SVD_solver.matrixU()*SI*SVD_solver.matrixV().transpose();
          if( L.determinant()<0){
              //cout<<R[i]<<endl;
              Matrix2d newV;
              newV <<SVD_solver.matrixV().transpose()(0, 0), SVD_solver.matrixV().transpose()(0, 1),
                 -(SVD_solver.matrixV().transpose()(1, 0)), -(SVD_solver.matrixV().transpose()(1, 1));
              SI(0,0)=SI(1,1)=0.5*(singular_value(0)-singular_value(1));
              if(method==ARAP) L=SVD_solver.matrixU()*newV;
              else if(method==ASAP) L=SVD_solver.matrixU()*SI*newV;
            }
            }
            if(method==Hybrid){
                double C1=0;
                double C2=0;
                double C3=0;
                for(int j=0;j<3;j++){
                    Vector2d v=half_edges[i*3+j].EdgeVec;
                    int a=half_edges[i*3+j].Endpoints(0);
                    int b=half_edges[i*3+j].Endpoints(1);
                    Vector2d u=(res.row(a)-res.row(b)).transpose();
                    C1+=weights(i*3+j)*(v(0)*v(0)+v(1)*v(1));
                    C2+=weights(i*3+j)*(u(0)*v(0)+u(1)*v(1));
                    C3+=weights(i*3+j)*(u(0)*v(1)-u(1)*v(0));
                }
                /*
                VectorXd init(2);
                init(0)=init(1)=1;
                R_Jacobian R_J(lamda,C1,C2,C3);
                R_Hessian R_H(lamda,C1);
                newton_optimizerNto1(init,R_J,R_H);*/
                double entry_1=0;
                double coe3=2*lamda*(1+C3*C3/C2/C2),coe1=C1-2*lamda,coe0=-C2;
                std::function<double(double)> cubic_fun=[&](double x)->double
                {
                    return coe3*x*x*x+coe1*x+coe0;
                };
                binary_find_root(entry_1,cubic_fun);
                L(0,0)=L(1,1)=entry_1;
                L(0,1)=entry_1*C3/C2;
                L(1,0)=-entry_1*C3/C2;
            }
            E+=area(i)* ( ( (J-L).transpose() )*(J-L) ).trace() ;
        }
        return E;
    }
    ARAP_energy(vector<HalfEdge>& half_edges,VectorXd& weights,VectorXd& area,int method,double lamda):half_edges(half_edges),weights(weights),area(area),method(method),lamda(lamda){}
private:
    VectorXd weights;
    VectorXd area;
    vector<HalfEdge> half_edges;
    double lamda;
    int method=ARAP;
};

#endif /* ARAP_h */
