//
//  main.cpp
//  arap
//
//  Created by Liang Chen on 2021/5/20.
//

#include "global_var.h"
#include "helper_init.h"
#include "helper_print.h"
#include "helper_geometry.h"
#include "ARAP.h"
#include "call_back.h"

int main(int argc, const char * argv[]) {
    /*initialization*/
    string input_name;
    int itrs=4;
    int method=ARAP;
    int func=PARAM;
    bool flip_avoid=false;
    bool print_pic=false;
    bool print_vtkfile=false;
    bool print_txtfile=true;
    bool print_each_frame=false;
    bool pause=false;
    bool inf_itr=false;
    bool show_texture=false;
    string texture_name;
    double lamda=0;
    string slamda="0";
    processArgv(argc,argv,input_name,itrs,func,method,flip_avoid,print_txtfile,print_vtkfile,print_pic,print_each_frame,pause,inf_itr,show_texture,texture_name,lamda,slamda);
    MatrixXd verts;
    MatrixXi edges;
    vector<HalfEdge> half_edges;
    VectorXd area;
    MatrixXi F;
    int vert_num=readObj(input_name,F,edges,verts,half_edges,area);
    
    MatrixXd res;
    MatrixXd RHS;
    VectorXd weights;
    SparseMatrix<double> Laplace;
    Laplace.resize(vert_num,vert_num);
    SimplicialLDLT<SparseMatrix<double> > dir_solver;
    MatrixXd* R=NULL;
    vector<int> fix;
    vector<VectorXd> fix_vec;
    string output_name=genOutputName(input_name,method,slamda,flip_avoid);
    igl::opengl::glfw::Viewer viewer;
    
    /*PARAMETERIZATION*/
    if(func==PARAM){
    /*Get initial solution(Barycentric Mapping)*/
    findBoundary(half_edges,edges,fix);
    mapTo2DBoundary(verts,fix,fix_vec,1);
    edges.resize(0,0);
    weights.resize(half_edges.size());
    getWeights(half_edges,verts,weights,Uniform);
    getLaplace(half_edges,weights,fix,Laplace,func);
    RHS.resize(vert_num,2);
    getRHS(R,half_edges,weights,fix,fix_vec,RHS,func,true);
    dir_solver.analyzePattern(Laplace);
    dir_solver.factorize(Laplace);
    assert(dir_solver.info()==Success);
    res.resize(vert_num,2);
    res=dir_solver.solve(RHS);
        
    /*Preparation for ARAP*/
    isometricProj(half_edges);
    getWeights(half_edges,verts,weights,Cotangent_1);
    R=new MatrixXd[half_edges.size()/3];
    for(int i=0;i<half_edges.size()/3;i++) R[i].resize(2,2);
    double distortion=0;
    double aread=0;
    double angled=0;
    double *distortion_per_unit=new double[half_edges.size()/3];
    double *aread_per_unit= new double[half_edges.size()/3];
    double *angled_per_unit= new double[half_edges.size()/3];
    ofstream distortion_file;
    if(!show_texture&&print_txtfile) distortion_file.open(output_name+".txt");
    fix.clear();
    fix_vec.clear();
    srand((unsigned int)time(NULL));
    fix.push_back(rand()%vert_num);
    for(int i=0;i<fix.size();i++) fix_vec.push_back(res.row(fix[i]).transpose());
    
    /*ARAP*/
    getLaplace(half_edges,weights,fix,Laplace,func);
    dir_solver.factorize(Laplace);
    assert( dir_solver.info()==Success);
    ARAP_energy E(half_edges,weights,area,method,lamda);
    if(show_texture){
        MatrixXd texture(res.rows(),2);
        texture=normalize_to_one2D(res);
        //read img using openCV
        cv::Mat img = imread(texture_name, cv::IMREAD_COLOR);
        if(img.empty())
        {
            std::cout << "Could not read the image: " << texture_name << std::endl;
            assert(0);
        }
        int width=img.cols,height=img.rows;
        Matrix<unsigned char,Dynamic,Dynamic> Red(width,height),Green(width,height),Blue(width,height);
        for(int x=0;x<width;x++)
        for(int y=0;y<height;y++){
            cv::Vec3b color=img.at<cv::Vec3b>(cv::Point(x,y));
            Red(x,y)=color[2];
            Green(x,y)=color[1];
            Blue(x,y)=color[0];
        }
        viewer.callback_key_pressed= [&](igl::opengl::glfw::Viewer &, unsigned int key, int)->bool
        {
            if(key==' '){
            local_phase_param(R,half_edges,res,weights,area,method,lamda,distortion_per_unit,aread_per_unit,angled_per_unit,distortion,aread,angled);
            global_phase(R,half_edges,weights,fix,fix_vec,RHS,func,false,res,dir_solver);
            texture=normalize_to_one2D(res);
            viewer.data().set_uv(texture,F);
            for(int i=0;i<fix.size();i++)
            fix_vec[i]=res.row(fix[i]).transpose();
            return true;
            }
            return false;
        };
        viewer.data().set_mesh(verts.transpose(),F);
        viewer.data().set_vertices(verts.transpose());
        viewer.data().set_uv(texture,F);
        viewer.data().show_lines = false;
        viewer.data().show_texture=true;
        viewer.data().set_texture(Red,Green,Blue);
        viewer.data().set_colors(RowVector3d(0.825,0.825,0.825));
        viewer.core().is_animating = true;
        viewer.launch(true,false,"show_texture",256,256);
    }
    else{
    MatrixXd new_res(vert_num,2);
    for(int now_itr=0;now_itr<itrs;now_itr++){
        //local phase: given p', solve R
        local_phase_param(R,half_edges,res,weights,area,method,lamda,distortion_per_unit,aread_per_unit,angled_per_unit,distortion,aread,angled);
        printFile(print_pic,print_vtkfile,print_txtfile,print_each_frame,now_itr,distortion,aread,angled,output_name,distortion_per_unit,aread_per_unit,angled_per_unit,half_edges,F,res,distortion_file);
        //global phase: given R, solve p'
        global_phase(R,half_edges,weights,fix,fix_vec,RHS,func,false,new_res,dir_solver);
        if(flip_avoid) igl::flip_avoiding_line_search(F,res,new_res,E);
        else res=new_res;
        for(int i=0;i<fix.size();i++)
        fix_vec[i]=res.row(fix[i]).transpose();
    }
    //The final "local optimization" in order to get the final energy
   local_phase_param(R,half_edges,res,weights,area,method,lamda,distortion_per_unit,aread_per_unit,angled_per_unit,distortion,aread,angled);
        print_each_frame=true;
        printFile(print_pic,print_vtkfile,print_txtfile,print_each_frame,itrs,distortion,aread,angled,output_name,distortion_per_unit,aread_per_unit,angled_per_unit,half_edges,F,res,distortion_file);
    }
        
    if(R) delete[] R;
    if(distortion_per_unit)    delete[] distortion_per_unit;
    if(aread_per_unit)    delete[] aread_per_unit;
    if(angled_per_unit)  delete[] angled_per_unit;
   return 0;
}
    
    /*DEFORMATION*/
    if(func==DEFORM) {
    /*Preparation for deformation*/
    bool placing_handles=true;
    bool first=true;
    long sel=-1;
    int now_itr=0;
    int print_num=0;
    bool moved=false;
    bool changed=false;
    RowVector3f last_mouse(0,0,0);
    edges.resize(0,0);
    R=new MatrixXd[vert_num];
    for(int i=0;i<vert_num;i++) {
            R[i].resize(3,3);
            R[i].setZero();
            R[i](0,0)= R[i](1,1)= R[i](2,2)=1;
    }
    weights.resize(vert_num*vert_num);
    getWeights(half_edges,verts,weights,Cotangent_2);
    sort(half_edges.begin(),half_edges.end(),compare);
    vector<int>* neighbors=new vector<int>[vert_num];
    getNeighbors(half_edges,neighbors);
    res.resize(vert_num,3);
    RHS.resize(vert_num,3);
    for(int i=0;i<vert_num;i++) res.row(i)=verts.col(i).transpose();
    //set viewer. The code of initial guess and ARAP optimizations for deformation is in the file call_back.h
    callbackMouseDown mouseDown(&placing_handles,&now_itr,&last_mouse,&res,&fix,&fix_vec,R,neighbors,&verts,&half_edges,&weights,&RHS,&F,&first,&moved,&changed,&sel,&dir_solver);
    callbackMouseMove mouseMove(&fix_vec,&last_mouse,&moved,&sel,&now_itr);
    callbackKeyPressed keyPressed(&placing_handles,&now_itr,&print_num,&res,&fix,&fix_vec,R,neighbors,&verts,&half_edges,&weights,&RHS,&F,&first,&moved,&changed,&Laplace,&dir_solver,output_name);
    callbackPreDraw preDraw(&placing_handles,&now_itr,&res,&fix,&fix_vec,R,neighbors,&verts,&half_edges,&weights,&RHS,&F,&first,&moved,&pause,&inf_itr,&dir_solver);
    viewer.callback_mouse_down=mouseDown;
    viewer.callback_mouse_move=mouseMove;
    viewer.callback_key_pressed=keyPressed;
    viewer.callback_pre_draw=preDraw;
    viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer&, int, int)->bool{ sel = -1; return false;};
    viewer.data().set_mesh(res,F);
    viewer.data().show_lines = false;
    viewer.core().is_animating = true;
    viewer.data().face_based = true;
    viewer.data().set_vertices(res);
    viewer.data().set_colors(RowVector3d(1.0,0.9,0.2));
    viewer.launch(true,false,"deformation",256,256);
    if(R) delete[] R;
    if(neighbors) delete[] neighbors;
    return EXIT_SUCCESS;
    }
}

 
