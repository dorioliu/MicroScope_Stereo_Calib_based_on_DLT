#include "stereocalib.h"
#include <opencv.hpp>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace cv;
using namespace std;
using namespace Eigen;

void calib()
{
    Mat imL = imread("left.png");
    Mat imR = imread("right.png");

    stereoCalib calib;

    calib.solveCameraParaByDLT(imL, imR);
   
}

int main(int* argc, int** argv)
{
    calib();
	return 0;
}

