//
//  helper_init.h
//  arap
//
//  Created by Liang Chen on 2021/6/14.
//

#ifndef helper_init_h
#define helper_init_h

#include "global_var.h"
#include <sstream>

inline void processArgv(int argc, const char* argv[], string& input_name,int& itrs,int& func,int& method,bool& flip_avoid ,bool& print_txtfile,bool& print_vtkfile,bool& print_pic,bool& print_each_frame,bool& pause,bool& inf_itr,bool& show_texture,string& texture_name,double& lamda,string& slamda)
{
    assert(argc>=1);
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-input")){
            i++;
            assert(i<argc);
            input_name=argv[i];
        }
        else if(!strcmp(argv[i],"-iterations")){
            i++;
            assert(i<argc);
            itrs=atoi(argv[i]);
        }
        else if(!strcmp(argv[i],"-method")){
            i++;
            assert(i<argc);
            if(!strcmp(argv[i],"ARAP")) method=ARAP;
            else if(!strcmp(argv[i],"ASAP")) method=ASAP;
            else if(!strcmp(argv[i],"Hybrid")){
                method=Hybrid;
                i++;
                assert(i<argc);
                slamda=argv[i];
                lamda=atof(argv[i]);
                if(lamda<0){
                    cout<<"Lambda of the Hybrid method must be positive!"<<endl;
                    assert(0);
                }
            }
            else {
                cout<<"Invalid method specified!"<<endl;
                assert(0);
            }
        }
        else if(!strcmp(argv[i],"-function")){
            i++;
            assert(i<argc);
            if(!strcmp(argv[i],"PARAM")) func=PARAM;
            else if(!strcmp(argv[i],"DEFORM")) func=DEFORM;
            else {
                cout<<"Invalid function specified!"<<endl;
                assert(0);
            }
        }
        else if(!strcmp(argv[i],"-flip_avoid")){
            flip_avoid=true;
        }
        else if(!strcmp(argv[i],"-print_pic")){
            print_pic=true;
        }
        else if(!strcmp(argv[i],"-no_txtfile")){
            print_txtfile=false;
        }
        else if(!strcmp(argv[i],"-print_vtkfile")){
            print_vtkfile=true;
        }
        else if(!strcmp(argv[i],"-print_each_frame")){
            print_each_frame=true;
        }
        else if(!strcmp(argv[i],"-pause")){
            pause=true;
        }
        else if(!strcmp(argv[i],"-inf_itr")){
            inf_itr=true;
        }
        else if(!strcmp(argv[i],"-show_texture")){
            show_texture=true;
            i++;
            assert(i<argc);
            texture_name=argv[i];
        }
        else{
            cout<<"Invalid input arguments specified!"<<endl;
        }
    }
    if(func==DEFORM) method=ARAP;
}

inline int readObj(const string& input_name,MatrixXi& F,MatrixXi& edges,MatrixXd& verts,vector<HalfEdge>& half_edges,VectorXd& area)
{
    ifstream input_file(input_name);
    string cur;
    int vert_num=0, face_num=0;
    input_file >> cur;
    /*two pass read*/
    //the first pass, count verts and facets

    while(!input_file.eof()){
        if(cur=="v") {
            vert_num++;
            getline(input_file,cur);
        }
        else if(cur=="f") {
            int vert_num_per_face=0;
            getline(input_file,cur);
            while(cur[0]==' '&&cur.length()>1) cur=cur.substr(1,cur.length());
            while(cur.find(" ")!=string::npos&&cur.find(" ")!=cur.length()-1){
                int i=cur.find(" ");
                cur=cur.substr(i+1,cur.length());
                while(cur[0]==' '&&cur.length()>1) cur=cur.substr(1,cur.length());
                vert_num_per_face++;
            }
            vert_num_per_face++;
            if(vert_num_per_face==3) face_num++;
            else{
                cout<<"Only triangle mesh are supported!"<<endl;
                assert(0);
            }
        }
       else getline(input_file,cur);
        input_file >> cur;
    }
    //initialize matrice
    verts.resize(3,vert_num);
    edges.resize(vert_num,vert_num);
    edges.fill(0);
    F.resize(face_num,3);
    //initialize vectors and scalars
    area.resize(face_num);
    double total_area=0;
    //the second pass, read verts and facets
    input_file.close();
    input_file.open(input_name);
    input_file >> cur;
    int cur_vert_num=0,cur_face_num=0;
    while(!input_file.eof()){
        if(cur=="v") {
        input_file >>verts(0,cur_vert_num)>>verts(1,cur_vert_num)>>verts(2,cur_vert_num);
        cur_vert_num++;
    }
        else if(cur=="f"){
            string as,bs,cs;
            input_file>>as>>bs>>cs;
            int enda=as.length(),endb=bs.length(),endc=cs.length();
            enda=as.find("/");
            endb=bs.find("/");
            endc=cs.find("/");
            int a=stoi(as.substr(0,enda))-1,b=stoi(bs.substr(0,endb))-1,c=stoi(cs.substr(0,endc))-1;
            F(cur_face_num,0)=a;
            F(cur_face_num,1)=b;
            F(cur_face_num,2)=c;
            //no repeatative triangles
            //there is no check for repeatative triangles
            Vector3d ab=verts.col(a)-verts.col(b);
            Vector3d ac=verts.col(a)-verts.col(c);
            area(cur_face_num)=fabs(0.5*ab.cross(ac).norm());
            total_area+=area(cur_face_num);
            HalfEdge tmpE1,tmpE2,tmpE3;
            if(!edges(a,b)){
                tmpE1.Endpoints=Vector2i(a,b);
                tmpE1.OppositePoint=c;
                tmpE1.BelongFacet=cur_face_num;
                tmpE1.EdgeVec=verts.col(a)-verts.col(b);
                half_edges.push_back(tmpE1);
                edges(a,b)=half_edges.size();
                if(edges(b,a)){
                    half_edges[edges(a,b)-1].InverseIdx=edges(b,a)-1;
                    half_edges[edges(b,a)-1].InverseIdx=edges(a,b)-1;
                }
            }
            else if(!edges(b,a)){
                tmpE1.Endpoints=Vector2i(b,a);
                tmpE1.OppositePoint=c;
                tmpE1.BelongFacet=cur_face_num;
                tmpE1.EdgeVec=verts.col(b)-verts.col(a);
                tmpE1.InverseIdx=edges(a,b)-1;
                half_edges.push_back(tmpE1);
                edges(b,a)=half_edges.size();
                half_edges[edges(a,b)-1].InverseIdx=edges(b,a)-1;
            }
            else{
                std::cout<<"Input is not a manifold!"<<std::endl;
                assert(0);
            }
            if(!edges(c,a)){
                tmpE2.Endpoints=Vector2i(c,a);
                tmpE2.OppositePoint=b;
                tmpE2.BelongFacet=cur_face_num;
                tmpE2.EdgeVec=verts.col(c)-verts.col(a);
                half_edges.push_back(tmpE2);
                edges(c,a)=half_edges.size();
                if(edges(a,c)){
                    half_edges[edges(c,a)-1].InverseIdx=edges(a,c)-1;
                    half_edges[edges(a,c)-1].InverseIdx=edges(c,a)-1;
                }
            }
            else if(!edges(a,c)){
                tmpE2.Endpoints=Vector2i(a,c);
                tmpE2.OppositePoint=b;
                tmpE2.BelongFacet=cur_face_num;
                tmpE2.EdgeVec=verts.col(a)-verts.col(c);
                tmpE2.InverseIdx=edges(c,a)-1;
                half_edges.push_back(tmpE2);
                edges(a,c)=half_edges.size();
                half_edges[edges(c,a)-1].InverseIdx=edges(a,c)-1;
            }
            else{
                std::cout<<"Input is not a manifold!"<<std::endl;
                assert(0);
            }
            if(!edges(b,c)){
                tmpE3.Endpoints=Vector2i(b,c);
                tmpE3.OppositePoint=a;
                tmpE3.BelongFacet=cur_face_num;
                tmpE3.EdgeVec=verts.col(b)-verts.col(c);
                half_edges.push_back(tmpE3);
                edges(b,c)=half_edges.size();
                if(edges(c,b)){
                    half_edges[edges(b,c)-1].InverseIdx=edges(c,b)-1;
                    half_edges[edges(c,b)-1].InverseIdx=edges(b,c)-1;
                }
            }
            else if(!edges(c,b)){
                tmpE3.Endpoints=Vector2i(c,b);
                tmpE3.OppositePoint=a;
                tmpE3.BelongFacet=cur_face_num;
                tmpE3.EdgeVec=verts.col(c)-verts.col(b);
                tmpE3.InverseIdx=edges(b,c)-1;
                half_edges.push_back(tmpE3);
                edges(c,b)=half_edges.size();
                half_edges[edges(b,c)-1].InverseIdx=edges(c,b)-1;
            }
            else{
                std::cout<<"Input is not a manifold!"<<std::endl;
                assert(0);
            }
            cur_face_num++;
        }
        else getline(input_file,cur);
        input_file >> cur;
    }
    input_file.close();
    assert(cur_vert_num==vert_num&&cur_face_num==face_num);
    //normalize the area
    double reg_param=sqrt(total_area);
    for(int i=0;i<face_num;i++) area(i)/=total_area;
    for(int i=0;i<half_edges.size();i++){
        half_edges[i].EdgeVec/=reg_param;
    }
    for(int i=0;i<vert_num;i++) verts.col(i)/=reg_param;
    return vert_num;
}

inline string genOutputName(string input_name,int method,string slamda,bool flip_avoid)
{
    string output_name=input_name.substr(0,input_name.find(".obj"));
    if(method==ARAP) output_name+="_ARAP";
    else if(method==ASAP) output_name+="_ASAP";
    else if(method==Hybrid){
    output_name+="_Hybrid_";
    output_name+=slamda;
}
if(flip_avoid) output_name+="_noflip";
    return output_name;
}

#endif /* helper_init_h */
