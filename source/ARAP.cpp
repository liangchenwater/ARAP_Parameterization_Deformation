//
//  ARAP.cpp
//  arap
//
//  Created by Liang Chen on 2021/6/14.
//

#include "ARAP.h"

void getLaplace(const vector<HalfEdge>& half_edges,const VectorXd& weights,const vector<int>& fix, SparseMatrix<double>& Laplace,int func)
{
    vector< Triplet<double> > tripletList;
    VectorXd row_sum;
    row_sum.resize(Laplace.rows());
    row_sum.fill(0);
    for(int i=0;i<fix.size();i++)
    row_sum(fix[i])=2*beta;
    for(int i=0;i<half_edges.size();i++){
        long long unsigned int a=half_edges[i].Endpoints(0);
        long long unsigned int b=half_edges[i].Endpoints(1);
        int inv_idx=half_edges[i].InverseIdx;
        double w=0;
        if(func==PARAM) w=weights(i);
        else if(func==DEFORM) w=weights(a*Laplace.rows()+b);
        if(inv_idx!=-1){
             if(func==PARAM) w+=weights(inv_idx);
              row_sum(a)+=w;
         tripletList.push_back(Triplet<double>(a,b,-w));
        }
        else{
            row_sum(a)+=w;
            row_sum(b)+=w;
            tripletList.push_back(Triplet<double>(a,b,-w));
            tripletList.push_back(Triplet<double>(b,a,-w));
        }
    }
    for(int i=0;i<Laplace.rows();i++) tripletList.push_back(Triplet<double>(i,i,row_sum(i)));
    Laplace.setFromTriplets(tripletList.begin(),tripletList.end());
}

void getRHS(MatrixXd* R,const vector<HalfEdge>& half_edges,const VectorXd& weights,const vector<int>& fix, const vector<VectorXd>& fix_vec,MatrixXd& RHS,int func,bool is_first)
{
    RHS.fill(0);
   for(int i=0;i<fix.size();i++)
    RHS.row(fix[i])=2*beta*fix_vec[i].transpose();
    
    if((func==PARAM&&!is_first)||func==DEFORM){
    for(int i=0;i<half_edges.size();i++){
       long long unsigned int a=half_edges[i].Endpoints(0);
        long long unsigned  int b=half_edges[i].Endpoints(1);
        int inv_idx=half_edges[i].InverseIdx;
        MatrixXd coeR1;
        MatrixXd coeR2;
        if(inv_idx!=-1){
            if(func==PARAM){
            coeR1=weights(i)*R[i/3];
            coeR2=weights(inv_idx)*R[inv_idx/3];
            RHS.row(a)+=(coeR1*half_edges[i].EdgeVec-coeR2*half_edges[inv_idx].EdgeVec).transpose();
            }
            else if(func==DEFORM){
                RHS.row(a)+=(0.5*weights(a*RHS.rows()+b)*(R[a]+R[b])*half_edges[i].EdgeVec).transpose();
            }
            }
        else{
            if(func==PARAM) {
            coeR1=weights(i)*R[i/3];
                RHS.row(a)+=(coeR1*half_edges[i].EdgeVec).transpose();
                RHS.row(b)-=(coeR1*half_edges[i].EdgeVec).transpose();
            }
            else if(func==DEFORM) {
                RHS.row(a)+=(0.5*weights(a*RHS.rows()+b)*(R[a]+R[b])*half_edges[i].EdgeVec).transpose();
                RHS.row(b)-=(0.5*weights(a*RHS.rows()+b)*(R[a]+R[b])*half_edges[i].EdgeVec).transpose();
            }
            }
    }
    }
}

void local_phase_deform(MatrixXd* R,vector<int>* neighbors,const MatrixXd& verts,const MatrixXd& res, const VectorXd& weights)
{
    int vert_num=verts.cols();
    for(int i=0;i<vert_num;i++){
        MatrixXd S=getCovariance3x3(neighbors[i],verts,res,weights,i);
        JacobiSVD<MatrixXd> SVD_solver;
        SVD_solver.compute(S,ComputeThinU | ComputeThinV);
        R[i]=SVD_solver.matrixU()*SVD_solver.matrixV().transpose();
        VectorXd singular_value=SVD_solver.singularValues();
        if(R[i].determinant()<0){
            double smallest_sv=1e16;
            double smallest_sv_idx=-1;
            for(int j=0;j<3;j++)
            if(singular_value(j)<smallest_sv){
                smallest_sv=singular_value(j);
                smallest_sv_idx=j;
            }
            MatrixXd newV=SVD_solver.matrixV().transpose();
            newV.row(smallest_sv_idx)=-1.0*newV.row(smallest_sv_idx);
            R[i]=SVD_solver.matrixU()*newV;
        }
        assert(R[i].determinant()>0);
    }
}

void local_phase_param(MatrixXd* R,const vector<HalfEdge>& half_edges,const MatrixXd& res,const VectorXd& weights,const VectorXd& area,int method,double lamda,double * distortion_per_unit,double* aread_per_unit,double *angled_per_unit,double& distortion,double& aread,double& angled)
{
    distortion=0;
    aread=0;
    angled=0;
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
        if(method==ARAP) R[i]=SVD_solver.matrixU()*SVD_solver.matrixV().transpose();
        else if(method==ASAP) R[i]=SVD_solver.matrixU()*SI*SVD_solver.matrixV().transpose();
      if( R[i].determinant()<0){
          //cout<<R[i]<<endl;
          Matrix2d newV;
          newV <<SVD_solver.matrixV().transpose()(0, 0), SVD_solver.matrixV().transpose()(0, 1),
             -(SVD_solver.matrixV().transpose()(1, 0)), -(SVD_solver.matrixV().transpose()(1, 1));
          SI(0,0)=SI(1,1)=0.5*(singular_value(0)-singular_value(1));
          if(method==ARAP) R[i]=SVD_solver.matrixU()*newV;
          else if(method==ASAP) R[i]=SVD_solver.matrixU()*SI*newV;
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
            R[i](0,0)=R[i](1,1)=entry_1;
            R[i](0,1)=entry_1*C3/C2;
            R[i](1,0)=-entry_1*C3/C2;
        }
        distortion_per_unit[i]= area(i)* ( ( (J-R[i]).transpose() )*(J-R[i]) ).trace() ;
        aread_per_unit[i]=area(i)*(singular_value(0)*singular_value(1)+1.0/singular_value(0)/singular_value(1));
        angled_per_unit[i]=area(i)*(singular_value(0)/singular_value(1)+singular_value(1)/singular_value(0));
        distortion+=distortion_per_unit[i];
        aread+=aread_per_unit[i];
        angled+=angled_per_unit[i];
    }
}

void global_phase(MatrixXd* R,const vector<HalfEdge>& half_edges, const VectorXd& weights, const vector<int>& fix,const vector<VectorXd>& fix_vec,MatrixXd& RHS,int func,bool first, MatrixXd& res,SimplicialLDLT<SparseMatrix<double> >& dir_solver)
{
    getRHS(R,half_edges,weights,fix,fix_vec,RHS,func,first);
    res=dir_solver.solve(RHS);
}
