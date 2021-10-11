#pragma once
#ifndef STEREOCALIB_H
#define STEREOCALIB_H

#include <opencv2/core.hpp>
#include <vector>
#include <Eigen/Dense>


class stereoCalib
{
public:
	stereoCalib();
	~stereoCalib();

	void initialWorldPoint();
	void imagePointDetection(cv::Mat& img, cv::Mat& imgR);
	void solveProjectionMatrix(std::vector<cv::Point2f>& img_pts, std::vector<cv::Point3f>& word_pts, Eigen::MatrixXd& solve);
	void solveQR(Eigen::MatrixXd& Ho, Eigen::MatrixXd& intriK);
	void solveCameraParaByDLT(cv::Mat& img, cv::Mat& imgR);

private:
	std::vector<cv::Point3f> m_world_pts;
	std::vector<cv::Point2f> m_img_ptsL, m_img_ptsR;
	cv::Size m_imSize;
	cv::Mat intrinsic_MatrixL, distortion_coeffsL, rvecsL, tvecsL;
	cv::Mat intrinsic_MatrixR, distortion_coeffsR, rvecsR, tvecsR;

	cv::Mat R, T,E,F;
};

#endif // !STEREOCALIB_H

