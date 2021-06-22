#ifndef _IMAGE_H_
#define _IMAGE_H_

#include <assert.h>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <Eigen/Core>

// ====================================================================
// ====================================================================
// Simple image class

class Image {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Image(int w, int h) {
    width = w;
    height = h;
    data = new Eigen::Vector3d[width*height]; }
  ~Image() {
    delete [] data; }

  // =========
  // ACCESSORS
  int Width() const { return width; }
  int Height() const { return height; }
  const Eigen::Vector3d& GetPixel(int x, int y) const {
    assert(x >= 0 && x < width);
    assert(y >= 0 && y < height);
    return data[y*width + x]; }
    Eigen::Vector3d GetPixel(int x,int y){
        assert(x >= 0 && x < width);
        assert(y >= 0 && y < height);
        return data[y*width + x];
    }
  // =========
  // MODIFIERS
  void SetAllPixels(const Eigen::Vector3d &color) {
    for (int i = 0; i < width*height; i++) {
      data[i] = color; } }
  void SetPixel(int x, int y, const Eigen::Vector3d &color) {
    assert(x >= 0 && x < width);
    assert(y >= 0 && y < height);
    data[y*width + x] = color; }

  // ===========
  // LOAD & SAVE
  static Image* LoadPPM(const char *filename);
  void SavePPM(const char *filename) const; 
  static Image* LoadTGA(const char *filename);
  void SaveTGA(const char *filename) const; 
    static Image* LoadImage(const char* filename);
  // extension for image comparison
  static Image* Compare(Image* img1, Image* img2);
  
private:

  // ==============
  // REPRESENTATION
  int width;
  int height;
  Eigen::Vector3d *data;

};

// ====================================================================
// ====================================================================

#endif
