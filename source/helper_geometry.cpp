//
//  helper_geometry.cpp
//  arap
//
//  Created by Liang Chen on 2021/6/15.
//

#include "helper_geometry.h"

void getNeighbors(const vector<HalfEdge>& half_edges, vector<int>* neighbors)
{
    for(int i=0;i<half_edges.size();i++){
        int a=half_edges[i].Endpoints(0);
        int b=half_edges[i].Endpoints(1);
        neighbors[a].push_back(b);
        if(half_edges[i].InverseIdx==-1) neighbors[b].push_back(a);
    }
}


void getWeights(const vector<HalfEdge>& half_edges,const MatrixXd& verts,VectorXd& weights,int Type)
{
    if(Type==Uniform){
        weights.fill(1);
    }
    else if(Type==Cotangent_1){
        for(int i=0;i<half_edges.size();i++){
            int c=half_edges[i].OppositePoint;
            int a=half_edges[i].Endpoints(0);
            int b=half_edges[i].Endpoints(1);
            VectorXd edge_for_a=verts.col(c)-verts.col(a);
            VectorXd edge_for_b=verts.col(c)-verts.col(b);
            edge_for_a.normalize();
            edge_for_b.normalize();
            double cos_theta=edge_for_a.dot(edge_for_b);
           // if(fabs(cos_theta-1)<1e-3||fabs(cos_theta+1)<1e-3) cos_theta=0;
            weights(i)=cos_theta/sqrt(1-cos_theta*cos_theta);
        }
    }
    else if(Type==Cotangent_2){
        for(int i=0;i<half_edges.size();i++){
            int c=half_edges[i].OppositePoint;
            int inv_idx=half_edges[i].InverseIdx;
           long long unsigned int a=half_edges[i].Endpoints(0);
            long long unsigned int b=half_edges[i].Endpoints(1);
            VectorXd edge_for_a=verts.col(c)-verts.col(a);
            VectorXd edge_for_b=verts.col(c)-verts.col(b);
            edge_for_a.normalize();
            edge_for_b.normalize();
            double cos_theta=fabs(edge_for_a.dot(edge_for_b));
            weights(a*verts.cols()+b)=0.5*cos_theta/sqrt(1-cos_theta*cos_theta);
            if(inv_idx!=-1) {
            int c1=half_edges[inv_idx].OppositePoint;
            VectorXd edge_for_a1=verts.col(c1)-verts.col(a);
            VectorXd edge_for_b1=verts.col(c1)-verts.col(b);
            edge_for_a1.normalize();
            edge_for_b1.normalize();
            double cos_theta1=fabs(edge_for_a1.dot(edge_for_b1));
            weights(a*verts.cols()+b)+=0.5*cos_theta1/sqrt(1-cos_theta1*cos_theta1);
            }
            if(inv_idx==-1) {
                weights(b*verts.cols()+a)=weights(a*verts.cols()+b);
            }
        }
    }
    else if(Type==MeanValue){
        for(int i=0;i<half_edges.size();i++){
        int c=half_edges[i].OppositePoint;
        int a=half_edges[i].Endpoints(0);
        int b=half_edges[i].Endpoints(1);
        double len=half_edges[i].EdgeVec.norm();
            Vector3d v0=verts.col(a)-verts.col(b);
            Vector3d v1=verts.col(a)-verts.col(c);
           v0.normalize();
            v1.normalize();
            double cos=v0.dot(v1);
            double alpha=acos(cos);
            weights(i)=tan(alpha*0.5)/len;
            int inv_idx=half_edges[i].InverseIdx;
            if(inv_idx!=-1){
                c=half_edges[inv_idx].OppositePoint;
                v1=verts.col(a)-verts.col(c);
                v1.normalize();
                cos=v0.dot(v1);
               alpha=acos(cos);
                weights(i)+=tan(alpha*0.5)/len;
            }
        
        }
    }
    //else if...
}



void mapTo2DBoundary(const MatrixXd& verts, const vector<int>& boundary_points, vector<VectorXd>& fix_vec,int alpha)
{
    double L=0;
    int bsize=boundary_points.size();
    //L is the circumference of the boundary polygon
    VectorXd Lk(bsize);
    for(int i=0;i<boundary_points.size();i++){
        double distance=( verts.col(boundary_points[i])-verts.col(boundary_points[(i+1)%bsize]) ).norm();
        L+=pow(distance,alpha);
        Lk(i)=L;
    }
    for(int i=0;i<boundary_points.size();i++) fix_vec.push_back(Vector2d(1.0/sqrt(M_PI)*cos(2.0*M_PI*Lk(i)/L),1.0/sqrt(M_PI)*sin(2.0*M_PI*Lk(i)/L)));
}


void isometricProj(vector<HalfEdge>& half_edges)
{
    for(int i=0;i<half_edges.size();i+=3){
        int a=half_edges[i].Endpoints(0);
        int b=half_edges[i].Endpoints(1);
        int c=half_edges[i+1].Endpoints(0);
        int c1=half_edges[i+1].Endpoints(1);
        bool flag=false;
        if(c==a||c==b) {
            int tmp=c;
            c=c1;
            c1=tmp;
            flag=true;
        }
        double dist_ab=half_edges[i].EdgeVec.norm();
        double Sin=-2,Cos=-2;
        if(c1==a){
            double dist_ac=half_edges[i+1].EdgeVec.norm();
            Cos=half_edges[i].EdgeVec.dot(half_edges[i+1].EdgeVec)/dist_ab/dist_ac;
            if(!flag) Cos*=-1.0;
            Sin=sqrt(1-Cos*Cos);
            half_edges[i].EdgeVec=Vector2d(-dist_ab,0);
            half_edges[i+1].EdgeVec=Vector2d(-dist_ac*Cos,-dist_ac*Sin);
            if(!flag) half_edges[i+1].EdgeVec*=-1.0;
            int q=half_edges[i+2].Endpoints(1);
            half_edges[i+2].EdgeVec=Vector2d(dist_ac*Cos-dist_ab,dist_ac*Sin);
            if(q==c) half_edges[i+2].EdgeVec*=-1.0;
        }
        else if(c1==b){
            double dist_ac=half_edges[i+2].EdgeVec.norm();
            int p=half_edges[i+2].Endpoints(0);
            int q=half_edges[i+2].Endpoints(1);
            Cos=half_edges[i].EdgeVec.dot(half_edges[i+2].EdgeVec)/dist_ab/dist_ac;
            if(p==c) Cos*=-1.0;
            Sin=sqrt(1-Cos*Cos);
            half_edges[i].EdgeVec=Vector2d(-dist_ab,0);
            half_edges[i+2].EdgeVec=Vector2d(dist_ac*Cos,dist_ac*Sin);
            if(q==c) half_edges[i+2].EdgeVec*=-1.0;
            half_edges[i+1].EdgeVec=Vector2d(dist_ac*Cos-dist_ab,dist_ac*Sin);
            if(flag) half_edges[i+1].EdgeVec*=-1.0;
        }
    }
}

void findBoundary(const vector<HalfEdge>& half_edges,const MatrixXi& edges,vector<int>& fix)
{
    int begin=-1,cur=-1,prev=-1;
    for(int i=0;i<half_edges.size();i++){
        if(half_edges[i].InverseIdx==-1){
            int a=half_edges[i].Endpoints(0);
            int b=half_edges[i].Endpoints(1);
            fix.push_back(a);
            fix.push_back(b);
            begin=a;
            prev=a;
            cur=b;
            break;
        }
    }
    while(cur!=begin){
        for(int i=0;i<edges.cols();i++)
        if( i!=prev&&((edges(cur,i)==0&&edges(i,cur)!=0) || (edges(cur,i)!=0&&edges(i,cur)==0) )){
            prev=cur;
            cur=i;
            fix.push_back(cur);
            if(fix.size()>=edges.cols()){
                cout<<"Model is not disk-like!"<<endl;
                assert(0);
            }
            break;
        }
    }
    if(fix.size()<3){
    cout<<"Model is not disk-like!"<<endl;
    assert(0);
    }
   if(!fix.empty()) fix.erase(fix.end()-1);
}

MatrixXd normalize_to_one2D(const MatrixXd& res)
{
    MatrixXd normalized_res(res.rows(),2);
    double xmax=-1e16,ymax=-1e16,xmin=1e16,ymin=1e16;
    //normalize the coordinates to (0,1) x (0,1)
    for(int i=0;i<res.rows();i++){
        if(res(i,0)<xmin) xmin=res(i,0);
        if(res(i,0)>xmax) xmax=res(i,0);
        if(res(i,1)<ymin) ymin=res(i,1);
        if(res(i,1)>ymax) ymax=res(i,1);
    }
    RowVector2d Min(xmin,ymin);
    double scale=1.0/std::max((xmax-xmin),(ymax-ymin));
    for(int i=0;i<res.rows();i++){
        normalized_res.row(i)=(res.row(i)-Min)*scale;
        if(std::max((xmax-xmin),(ymax-ymin))==xmax-xmin) normalized_res(i,1)+=0.5*(1.0-scale*(ymax-ymin));
        else normalized_res(i,0)+=0.5*(1.0-scale*(xmax-xmin));
    }
    return normalized_res;
}
