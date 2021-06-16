# ARAP_Parameterization_Deformation
- Commited to 2020-2021 Spring/Summer Zhejiang University Advanced Computer Graphics course.
- Please do not cheat. For good of both of us.
- ARAP Deformation and Parameterization.
- Deformation is interactive.
- Xcode project. MacOS executable file. Compiled with -ofast optimization.
- Dependencies: Eigen/igl(No core algorithms about ARAP are used. Only used for user interaction and flip avoid line searching.)/OpenGL/GLFW
- Some code for user interaction (in call_back.h) adapted from https://github.com/alecjacobson/geometry-processing-deformation
- image.h and image.cpp adapted from http://10.76.1.181/courses/training/mitF04/assignments/
- Arguments for parameterization are like: '-function PARAM -input bunny.obj -iterations 4 -method Hybrid 1e-6 -print_pic -print_vtkfile -print_each_frame'
- Arguments for deformation are like: '-function DEFORM -input bunny.obj -inf_itr'
- A strict disk-like manifold check is not included. Please be sure that the input is disk-like so that results are reasonable for PARAM.
- User interactions when deformation: 
	- mouse left: select anchor points
	- mouse drag: deform interactively
	- press ' ': change between anchor selection mode and deformation mode
	- press 'u': single update
	- press 'g': print deformed .obj file
- For more data, including results and models, please refer to BAIDU Cloud Disk: https://pan.baidu.com/s/1_ObT-NQUKhNzIlSJ6tfFCw  password: aaaa
