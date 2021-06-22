//
//  helpler_geometry.h
//  arap
//
//  Created by Liang Chen on 2021/6/14.
//

#ifndef helper_geometry_h
#define helper_geometry_h

#include "global_var.h"

void getNeighbors(const vector<HalfEdge>& half_edges, vector<int>* neighbors);

void getWeights(const vector<HalfEdge>& half_edges,const MatrixXd& verts,VectorXd& weights,int Type);

void mapTo2DBoundary(const MatrixXd& verts, const vector<int>& boundary_points, vector<VectorXd>& fix_vec,int alpha);

void isometricProj(vector<HalfEdge>& half_edges);

void findBoundary(const vector<HalfEdge>& half_edges,const MatrixXi& edges,vector<int>& fix);

MatrixXd normalize_to_one2D(const MatrixXd& res);
#endif /* helper_geometry_h */
