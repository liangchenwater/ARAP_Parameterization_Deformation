//
//  helper_print.h
//  arap
//
//  Created by Liang Chen on 2021/6/14.
//

#ifndef helper_print_h
#define helper_print_h

#include "global_var.h"
#include "image.h"
#include "helper_geometry.h"

inline void printVTK(const vector<HalfEdge>& half_edges, const MatrixXi& F,double* distortion,
                                             double* aread,double * angled,const MatrixXd& res ,string output_name)
{
    //double sum1=0,sum2=0,sum3=0;
    output_name+=".vtk";
    ofstream vtkfile(output_name);
    vtkfile<<"# vtk DataFile Version 3.0"<<endl;
    vtkfile<<output_name<<endl;
    vtkfile<<"ASCII"<<endl;
    vtkfile<<"DATASET POLYDATA"<<endl;
   // vtkfile<<"DIMENSIONS 1 1"<<endl<<endl;
    vtkfile<<"POINTS "<<res.rows()<<" "<<"double"<<endl;
    for(int i=0;i<res.rows();i++){
        vtkfile<<res.row(i)<<" ";
        vtkfile<<0<<endl;
    }
    vtkfile<<endl;
    vtkfile<<"POLYGONS "<<F.rows()<<" "<<4*F.rows()<<endl;
    for(int i=0;i<F.rows();i++)
    vtkfile<<"3 "<<F.row(i)<<endl;
    vtkfile<<endl;
    vtkfile<<"CELL_DATA "<<F.rows()<<endl;
    vtkfile<<"SCALARS distortion double 1"<<endl;
    vtkfile<<"LOOKUP_TABLE distortion_table"<<endl;
    for(int i=0;i<F.rows();i++){
        vtkfile<<distortion[i]<<endl;
       // sum1+=distortion[i];
    }
    vtkfile<<endl;
    vtkfile<<"SCALARS aread double 1"<<endl;
    vtkfile<<"LOOKUP_TABLE aread_table"<<endl;
    for(int i=0;i<F.rows();i++) {
        vtkfile<<aread[i]<<endl;
       // sum2+=aread[i];
    }
    vtkfile<<endl;
    vtkfile<<"SCALARS angled double 1"<<endl;
    vtkfile<<"LOOKUP_TABLE angled_table"<<endl;
    for(int i=0;i<F.rows();i++){
        vtkfile<<angled[i]<<endl;
        //sum3+=angled[i];
    }
    vtkfile<<endl;
    vtkfile.close();
}

inline void twoDDDA(Vector2d a, Vector2d b, Image& pic,int r)
{
    Vector2d distance=b-a;
    int steps=std::max(abs(distance(0)),abs(distance(1)));
    double x=a(0),y=a(1);
    double xinc=distance(0)/steps;
    double yinc=distance(1)/steps;
    for(int i=0;i<=steps;i++){
        pic.SetPixel(x,y,Vector3d(r,0,0));
        x+=xinc;
        y+=yinc;
    }
}

inline void printImage(const vector<HalfEdge>& half_edges,const MatrixXd& res,string output_name)
{
    int width=4000,height=4000;
    Image pic(width,height);
    MatrixXd normalized_res(res.rows(),2);
    normalized_res=normalize_to_one2D(res);
    Vector3d white(1.0,1.0,1.0);
    pic.SetAllPixels(white);
    //cout<<xmax<<endl<<xmin<<endl<<ymax<<endl<<ymin<<endl;
    double screen_size=std::min(width,height);
    double scale=1.0*(screen_size-1);
    for(int i=0;i<res.rows();i++){
        normalized_res.row(i)=normalized_res.row(i)*scale;
        //cout<<normalized_res.row(i)<<endl;
        if(screen_size==height) normalized_res(i,0)+=0.5*(width-screen_size);
        else normalized_res(i,1)+=0.5*(height-screen_size);
    }
    VectorXi vis(half_edges.size());
    vis.fill(0);
    for(int i=0;i<half_edges.size();i++)
    if(!vis(i)){
        int a=half_edges[i].Endpoints(0);
        int b=half_edges[i].Endpoints(1);
        if(half_edges[i].InverseIdx==-1)  twoDDDA(normalized_res.row(a).transpose(),normalized_res.row(b).transpose(),pic,1);
        else twoDDDA(normalized_res.row(a).transpose(),normalized_res.row(b).transpose(),pic,0);
        vis(i)=1;
        if(half_edges[i].InverseIdx!=-1){
            vis(half_edges[i].InverseIdx)=1;
        }
    }
    output_name+=".tga";
    pic.SaveTGA(output_name.c_str());
}

inline void printFile(bool print_pic,bool print_vtkfile,bool print_txtfile,bool print_each_frame,int now_itr,double distortion,double aread,double angled,string output_name, double* distortion_per_unit, double* aread_per_unit, double * angled_per_unit,const vector<HalfEdge>& half_edges, const MatrixXi& F,const MatrixXd& res,ofstream& distortion_file)
{
    string itr;
    int cnt=now_itr;
    if(cnt==0) itr.push_back('0');
    while(cnt){
        itr.push_back(cnt%10+'0');
        cnt/=10;
    }
    for(int i=0;i<itr.size()/2;i++){
        char tmp=itr[i];
        itr[i]=itr[itr.size()-1-i];
        itr[itr.size()-1-i]=tmp;
    }
    
    if(print_pic&&print_each_frame)  printImage(half_edges,res,output_name+"_"+itr);
    if(print_vtkfile&&print_each_frame) printVTK(half_edges,F,distortion_per_unit,
                                                 aread_per_unit,angled_per_unit,res,output_name+"_"+itr);
      //write file
      if(print_txtfile){
      distortion_file<<"After the "<<now_itr<<"th iterations"<<endl;
      distortion_file<<"Total Energy: "<<distortion<<endl;
      distortion_file<<"Total Area Distortion: "<<aread<<endl;
      distortion_file<<"Total Angle Distortion: "<<angled<<endl<<endl;
      }
}

inline void printObjModel(const MatrixXd& res,const MatrixXi& F,string output_name,int num)
{
    string itr;
    int cnt=num;
    if(cnt==0) itr.push_back('0');
    while(cnt){
        itr.push_back(cnt%10+'0');
        cnt/=10;
    }
    for(int i=0;i<itr.size()/2;i++){
        char tmp=itr[i];
        itr[i]=itr[itr.size()-1-i];
        itr[itr.size()-1-i]=tmp;
    }
    
    ofstream Out(output_name+"_"+itr+".obj");
    for(int i=0;i<res.rows();i++) Out<<"v "<<res.row(i)<<endl;
    for(int i=0;i<F.rows();i++) Out<<"f "<<F.row(i)+RowVector3i(1,1,1)<<endl;
    Out.close();
}
#endif /* helper_print_h */
