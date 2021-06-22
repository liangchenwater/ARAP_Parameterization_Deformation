//
//  globar_var.h
//  arap
//
//  Created by Liang Chen on 2021/6/14.
//

#ifndef global_var_h
#define global_var_h
#define STB_IMAGE_IMPLEMENTATION
#include <iostream>
#include <vector>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <string>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <igl/flip_avoiding_line_search.h>
#include <igl/project.h>
#include <igl/unproject.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/snap_points.h>
#include <igl/unproject_onto_mesh.h>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include "image.h"
#define beta 50000
#define eps 1e-6
#define newton_step_size 0.5
#define MAX_NEWTON_ITRS 200
#define inf_for_bfr 1e20

using namespace std;
using namespace Eigen;

struct HalfEdge{
    Vector2i Endpoints;
    VectorXd EdgeVec;
    int OppositePoint;
    int BelongFacet;
    int InverseIdx=-1;
    //friend bool operator < (const HalfEdge& a,const HalfEdge& b);
};

inline bool compare (HalfEdge a,HalfEdge b){
    if(a.Endpoints(0)!=b.Endpoints(0)) return a.Endpoints(0)<b.Endpoints(0);
    else return a.Endpoints(1)<b.Endpoints(1);
}

enum {Uniform,Wachspress,DH,MeanValue,Cotangent_1,Cotangent_2};
enum {ARAP,ASAP,Hybrid};
enum {PARAM,DEFORM};


#endif /* global_var_h */
