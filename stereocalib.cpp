#include "stereocalib.h"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <iostream>
#include <opencv.hpp>

using namespace std;
using namespace Eigen;
using namespace cv;

stereoCalib::stereoCalib()
{
	;
}

stereoCalib::~stereoCalib()
{
	;
}

void stereoCalib::initialWorldPoint()
{
	// 2 initial 3d world control points 
	for (int i = 8; i > 0; i--)       // 
	{
		m_world_pts.push_back(Point3f(24.3f * 6, 0, 24.3f * i));
		m_world_pts.push_back(Point3f(24.3f * 5, 0, 24.3f * i));
		m_world_pts.push_back(Point3f(24.3f * 4, 0, 24.3f * i));
		m_world_pts.push_back(Point3f(24.3f * 3, 0, 24.3f * i));
		m_world_pts.push_back(Point3f(24.3f * 2, 0, 24.3f * i));
		m_world_pts.push_back(Point3f(24.3f * 1, 0, 24.3f * i));

		m_world_pts.push_back(Point3f(0, 24.3f * 1, 24.3f * i));
		m_world_pts.push_back(Point3f(0, 24.3f * 2, 24.3f * i));
		m_world_pts.push_back(Point3f(0, 24.3f * 3, 24.3f * i));
		m_world_pts.push_back(Point3f(0, 24.3f * 4, 24.3f * i));
		m_world_pts.push_back(Point3f(0, 24.3f * 5, 24.3f * i));
		m_world_pts.push_back(Point3f(0, 24.3f * 6, 24.3f * i));
	}
}

void stereoCalib::imagePointDetection(cv::Mat& imgL, cv::Mat& imgR)
{
	m_imSize = imgL.size();

	Mat grayL,grayR;
	cvtColor(imgL, grayL, cv::COLOR_BGR2GRAY);
	cvtColor(imgR, grayR, cv::COLOR_BGR2GRAY);

	

	/*bool retL = cv::findChessboardCorners(grayL,
		cv::Size(12, 8),
		m_img_ptsL,
		cv::CALIB_CB_ADAPTIVE_THRESH);

	bool retR = cv::findChessboardCorners(grayR,
		cv::Size(12, 8),
		m_img_ptsR,
		cv::CALIB_CB_ADAPTIVE_THRESH);*/

	vector<Point2f> mpL, mpR;
	bool retL = cv::findChessboardCorners(grayL,
		cv::Size(12, 8),
		mpL,
		cv::CALIB_CB_ADAPTIVE_THRESH);

	bool retR = cv::findChessboardCorners(grayR,
		cv::Size(12, 8),
		mpR,
		cv::CALIB_CB_ADAPTIVE_THRESH);

	for (int i = 0; i < 12*8; i++)       //
	{
		m_img_ptsL.push_back(mpL[i]);
		m_img_ptsR.push_back(mpR[i]);
	}



	////指定亚像素计算迭代标注
	//cv::TermCriteria criteria = cv::TermCriteria(
	//    cv::TermCriteria::MAX_ITER + cv::TermCriteria::EPS,
	//    40,
	//    0.1);

	////亚像素检测
	//cv::cornerSubPix(grayL, corners, cv::Size(5, 5), cv::Size(-1, -1), criteria);

	//角点绘制
	//cv::drawChessboardCorners(imL, cv::Size(12, 8), corners, ret);
}

void stereoCalib::solveProjectionMatrix(std::vector<cv::Point2f>& img_pts, std::vector<cv::Point3f>& word_pts, Eigen::MatrixXd& solve)
{
	int pts_num = word_pts.size();

	if (word_pts.size() < 7 || img_pts.size() < 7)
	{
		std::cerr << "wrong points num" << std::endl;
		return;
	}

	if (word_pts.size() != img_pts.size())
	{
		std::cerr << "num of img_pts are not equal to the one world pts have" << std::endl;
		return;
	}
		
	// Ax = 0 
	// inital Matrix A
	Eigen::MatrixXd A_h = Eigen::MatrixXd::Zero(2*word_pts.size(), 12);
		
		//Constant(word_pts.size(), 12, 0);

	cout << "rows: " << A_h.rows() << endl;
	cout << "cols: " << A_h.cols() << endl;

	for (int i = 0; i < pts_num; i++)
	{
		MatrixXd m_world_pts(1,4);
		m_world_pts << word_pts[i].x, word_pts[i].y, word_pts[i].z, 1;
		MatrixXd m_img_pts(1,2);
		m_img_pts << img_pts[i].x, img_pts[i].y;

		A_h.block(i * 2, 0, 1, 4) = -1* m_world_pts;
		A_h.block(i * 2, 8, 1, 4) = m_img_pts(0) * m_world_pts;
		A_h.block(i * 2 + 1, 4, 1, 4) = -1 * m_world_pts;
		A_h.block(i * 2 + 1, 8, 1, 4) = m_img_pts(1) * m_world_pts;
	}

	JacobiSVD<Eigen::MatrixXd> svd(A_h, ComputeThinU | ComputeThinV);

	 // A_h = USV^T
	MatrixXd V(12,12), U(2 * word_pts.size(), 12);

	V = svd.matrixV();
	U = svd.matrixU();

	auto S_singular = svd.singularValues();

	MatrixXd S_M = MatrixXd::Zero(12,12);
	S_M.diagonal() = S_singular.transpose();
	
	// 
	/*
	std::cout << "A_h :\n" << A_h << std::endl;
	std::cout << "USVT :\n" << U* S_M*V.transpose() << std::endl;*/
	
	// please don't combine this two process as solve = V.col(11).reshaped<RowMajor>(3, 4).eval();
	solve = V.col(11);
	solve = solve.reshaped<RowMajor>(3, 4).eval();   
	
	solve = solve / solve(2,3);

	//std::cout << "solve :\n" << solve << std::endl;

	/*solveQR(solve);*/
}

void stereoCalib::solveQR(Eigen::MatrixXd& Ho, Eigen::MatrixXd& intriK)
{ 
	// Ho 3by4 with coefficent of last row last col is 1
	// H = KR[I, -C] = [KR, -KRC] = [H1 h] = K[R, -RC] = K[R, t]
	//  

	/*cout << "H\n" << Ho << endl;
	cout << "Hx\n" << Ho.block(0, 0, 3, 3) << endl;*/

	// H1 = KR
	// H1^-1 = R^-1 (K^-1)
	Matrix3d H1 = Ho.block(0, 0, 3, 3).inverse();

	//cout << "H1\n" << H1 << endl;

	HouseholderQR<MatrixXd> qr;
	qr.compute(H1);

	MatrixXd R_qr = qr.matrixQR().triangularView<Eigen::Upper>();
	MatrixXd Q_qr = qr.householderQ();

	MatrixXd K = R_qr.inverse();
	MatrixXd R = Q_qr.inverse();

	//cout << "K\n" << K << endl;

	K = K / K(2, 2);

	/*cout << "K\n" << K << endl;
	cout << "R\n" << R << endl;*/

	intriK = K;

	return;

}

/*********
	* This is an implementation of DLT based camera calibration with Eigen
	*
	* General,we describe a imaging system of pinhole model as:
	* sp = K[R|t]P, s is depth scale coefficient, lower case p
	* represents 2d image point (u,v,1), upper P is a space 3d
	* point (X,Y,Z,1)  in world coordinate.
	*
	* usually, in rigid transformation, we descibe a rigid motion
	* from P1(X1,Y1,Z1) to P2(X2,Y2,Z2) with rotation transformation
	* R and the new camera position is C as P2 = R(P1 - C), also
	* sometimes, we hide description of C and introduces a tanslation
	* vecton t = -RC ,in this case, rigid transformation can also
	* be written down as P2 = RP1 + t.
	*
	* So:
	*     sp = K[R|t]P = [KR, -KRC]P = HP  (1)
	*
	* which C reperents the camera optical center position in world
	* coordinate, H is a preojection matrix with a shape of 3 by 4.
	*
	* The DLT method takes a set of at least 6 pairs of 3d world
	* points and those corrsponding projected 2d iamge points as
	* inputs and outputs the K,R parameters what we want to get
	* by DLT based calibration. please note that all 3d world points
	* must not space on a space planar.
	*
	* So, how to evaluate H and how to do decompsition about H
	* to get K and R ?
	* in general, we calculate H by build over-determined equation
	* based on the relation between p and P.
	*
	* [u]    [ h11, h12, h13, h14 ]
	* [v] =  [ h21, h22, h23, h24 ] * [X, Y,Z,1]^T   (2)
	* [1]    [ h31, h32, h33, h34 ]
	*
	* thus, we get two quations:
	*
	* u = (hllX + h12Y + h13Z + 1)/(h31X + h31Y + h33Z + 1)      (3)
	* v = (h2lX + h22Y + h23Z + 1)/(h31X + h31Y + h33Z + 1)
	*
	* now we have I pairs of that point (ui,vi,1) and (Xi,Yi,Zi,1)
	* each point pair provide two equation.
	*
	* based on (3) we get a linear algebraic equation:
	*
	* [ -X1 -Y1 -Z1 -1 0 0 0 0 x1X1 x1Y1 x1Z1 x1 ]
	* [ 0 0 0 0 -X1 -Y1 -Z1 -1 y1X1 y1Y1 y1Z1 y1 ]
	*                      .....
	* [ -Xi -Yi -Zi -1 0 0 0 0 xiXi xiY1 xiZi xi ]    * [h11, h12, h13, h14,h21, h22, h23, h24,h31, h32, h33, h34 ]^T = 0 (4)
	* [ 0 0 0 0 -Xi -Yi -Zi -1 yiXi yiYi yiZi yi ]
	*                     .....
	* [ -XI -YI -ZI -1 0 0 0 0 xIXI xIYI xIZI xI ]
	* [ 0 0 0 0 -XI -YI -ZI -1 yIXI yIYI yIZI yI ]
	*
	* this linear equation is a traditional Ax = 0 problem.
	* the sulution of x can be evaluated by SVD decomposition
	* (sigular value decomposition)
	* at that case:
	* A = USV^T, U: 2Ix12, V: 12x12, S: 12x12 is diagonal matrix with descendent value
	* the solution is a sigular vector within V corresponding to the minimal
	* sigular value among S diagonal line coefficent， tipically x = V.col(11) as the repository.
	*
	* now we have got H, and we use QR decomposition to calculate K , R ,C
	*
	* H = [KR, -KRC] = [Ho, ho]
	* Ho = KR
	* ho = -KRC
	* so C = -Ho^-1 * ho
	* to our konwledge, in QR decomposition, R is a upper
	* triangular matrix, Q is a orthogonal matrix, and is our Ho =KR
	* matrix K is a upper triangular matrix, R is a orthogonal matrix
	*
	* so  we can make a samll trick H^-1 = R^-1 * K^-1 = QR
	* than we get rotation matrix R = Q^-1, and get K = R^-1
	* finnally we devide K by K(2,3) to get a normalized K.
	*
**********/
void stereoCalib::solveCameraParaByDLT(cv::Mat& img, cv::Mat& imgR)
{
	//inital 3d world points
	initialWorldPoint();

	// deterct cheeseboard
	imagePointDetection(img, imgR);

	// evaluate projection matrix
	Eigen::MatrixXd solveProjectMatrixL, solveProjectMatrixR;
	solveProjectionMatrix(m_img_ptsL, m_world_pts, solveProjectMatrixL);
	solveProjectionMatrix(m_img_ptsR, m_world_pts, solveProjectMatrixR);

	// evaluate camera matirx
	MatrixXd intriKL, intriKR;
	solveQR(solveProjectMatrixL, intriKL);
	solveQR(solveProjectMatrixR, intriKR);


	// evalute the final k,distortion,R,and T based on opencv 
	// with a guassed intrinsic matrix got by DLT
	vector<vector<cv::Point3f>> world_ptss;
	vector<vector<cv::Point2f>> img_ptssL, img_ptssR;

	world_ptss.push_back(m_world_pts);
	img_ptssL.push_back(m_img_ptsL);
	img_ptssR.push_back(m_img_ptsR);

	
	intrinsic_MatrixL = Mat::zeros(3, 3, CV_64F);
	intrinsic_MatrixL.at<double>(0, 0) = abs(intriKL(0,0));
	intrinsic_MatrixL.at<double>(0, 2) = abs(intriKL(0, 2));
	intrinsic_MatrixL.at<double>(1, 1) = abs(intriKL(1, 1));
	intrinsic_MatrixL.at<double>(2, 2) = 1.0;

	double rmsL = calibrateCamera(world_ptss, img_ptssL, m_imSize, intrinsic_MatrixL, distortion_coeffsL, rvecsL, tvecsL, CALIB_USE_INTRINSIC_GUESS);

	cout << "RMSL error: " << rmsL << endl;
	cout << "KL:\n " << intrinsic_MatrixL << endl;
	cout << "distortionL: \n" << distortion_coeffsL << endl;
	cout << "rvecsL: \n" << rvecsL << endl;
	cout << "tvecsL: \n" << tvecsL << endl;

	
	intrinsic_MatrixR = Mat::zeros(3, 3, CV_64F);
	intrinsic_MatrixR.at<double>(0, 0) = abs(intriKR(0, 0));
	intrinsic_MatrixR.at<double>(0, 2) = abs(intriKR(0, 2));
	intrinsic_MatrixR.at<double>(1, 1) = abs(intriKR(1, 1));
	intrinsic_MatrixR.at<double>(1, 2) = abs(intriKR(1, 2));
	intrinsic_MatrixR.at<double>(2, 2) = 1.0;

	double rmsR = calibrateCamera(world_ptss, img_ptssR, m_imSize, intrinsic_MatrixR, distortion_coeffsR, rvecsR, tvecsR, CALIB_USE_INTRINSIC_GUESS);

	cout << "RMSR error: " << rmsR << endl;
	cout << "KR:\n " << intrinsic_MatrixR << endl;
	cout << "distortionR: \n" << distortion_coeffsR << endl;
	cout << "rvecsR: \n" << rvecsR << endl;
	cout << "tvecsR: \n" << tvecsR << endl;


	double rmsStereo = stereoCalibrate(world_ptss, img_ptssL, img_ptssR, intrinsic_MatrixL, distortion_coeffsL,
		intrinsic_MatrixR, distortion_coeffsR, m_imSize, R,
		T, E, F, CALIB_USE_INTRINSIC_GUESS, TermCriteria(TermCriteria::COUNT + TermCriteria::EPS, 30, 1e-6));
		
	cout << "final KL:\n " << intrinsic_MatrixL << endl;
	cout << "final distortionL: \n" << distortion_coeffsL << endl;
	cout << "final KR:\n " << intrinsic_MatrixR << endl;
	cout << "final distortionR: \n" << distortion_coeffsR << endl;

	cout << "R: \n" << R << endl;
	cout << "T: \n" << T << endl;

	
	return;
}

