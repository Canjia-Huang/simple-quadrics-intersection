#include "simple_quadrics_intersection.h"

namespace QuadricsIntersection 
{
	int get_intersections(
		Plane& P1, Sphere& S1,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between a Plane and a Sphere");

		// init 
		std::vector<Point>().swap(points);
		std::vector<ParameterizationCircle>().swap(circles);

		Eigen::Vector3d S1_cor_to_P1_cor = P1.cor() - S1.cor();
		double P1_S1_cor_dot_P1_nor = S1_cor_to_P1_cor.dot(P1.nor());
		double P1_S1_cor_dis = std::abs(P1_S1_cor_dot_P1_nor);

		if (P1_S1_cor_dis > S1.r() + SQI_EPS) { // not intersect
			SQI_VERBOSE_ONLY_COUT("plane -- not intersect -- sphere");
		}
		else {
			if (P1_S1_cor_dis > S1.r() - SQI_EPS) { // tangent
				SQI_VERBOSE_ONLY_COUT("plane -- tangent -- sphere");

				points.push_back(Point(S1.cor() + P1_S1_cor_dot_P1_nor * P1.nor()));
			}
			else { // intersect
				SQI_VERBOSE_ONLY_COUT("plane -- intersect -- sphere");

				Eigen::Vector3d circle_center = S1.cor() + P1_S1_cor_dot_P1_nor * P1.nor();
				Eigen::Vector3d circle_nor = P1.nor();
				if (P1_S1_cor_dot_P1_nor < 0) circle_nor = -circle_nor;
				double circle_r = std::sqrt(S1.r() * S1.r() - P1_S1_cor_dis * P1_S1_cor_dis);

				circles.push_back(ParameterizationCircle(circle_center, -P1.nor(), circle_r));
			}
		}

		return points.size() + circles.size();
	}

	int get_intersections(
		Cylinder& C1, Sphere& S1,
		std::vector<Point>& points,
		std::vector<ParameterizationCylindricCurve>& curves
	) {
		return get_intersections(S1, C1, points, curves);
	}

	int get_intersections(
		Sphere& S1, Sphere& S2,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between a Sphere and a Sphere");

		// init
		std::vector<Point>().swap(points);
		std::vector<ParameterizationCircle>().swap(circles);

		Eigen::Vector3d S2_cor_to_S1_cor = S1.cor() - S2.cor();
		double centers_dis = S2_cor_to_S1_cor.norm();

		if (centers_dis < SQI_EPS) { // overlap: result nothing
			SQI_VERBOSE_ONLY_COUT("sphere -- overlap -- sphere");
		}
		else {
			Eigen::Vector3d S2_cor_to_S1_cor_unit = S2_cor_to_S1_cor / centers_dis;

			if (centers_dis > S1.r() + S2.r() + SQI_EPS) { // not intersect
				SQI_VERBOSE_ONLY_COUT("sphere -- not intersect -- sphere");
			}
			else if (centers_dis > S1.r() + S2.r() - SQI_EPS) { // tangent: result a point
				SQI_VERBOSE_ONLY_COUT("sphere -- tangent -- sphere");

				points.push_back(Point(S2.cor() + S2.r() * S2_cor_to_S1_cor_unit));
			}
			else { // intersect: result a circle
				SQI_VERBOSE_ONLY_COUT("sphere -- intersect -- sphere");

				Eigen::Vector3d circle_center = S2.cor() + 0.5 * S2_cor_to_S1_cor;
				double circle_r = std::sqrt(S2.r() * S2.r() - 0.25 * centers_dis * centers_dis);

				circles.push_back(
					ParameterizationCircle(circle_center, S2_cor_to_S1_cor_unit, circle_r));
			}
		}

		return points.size() + circles.size();
	}

	int get_intersections(
		Sphere& S1,
		std::vector<Plane>& planes, std::vector<Cylinder>& cylinders, std::vector<Sphere>& spheres,
		std::vector<Point>& res_points, std::vector<ParameterizationCircle>& res_circles, std::vector<ParameterizationCylindricCurve>& res_curves
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between a Sphere and other primitives");

		// init
		std::vector<Point>().swap(res_points);
		std::vector<ParameterizationCircle>().swap(res_circles);
		std::vector<ParameterizationCylindricCurve>().swap(res_curves);

		// get the intersections between sphere S1 and other primitives
		std::vector<std::vector<ParameterizationCircle>> res_c_circles;
		std::vector<std::vector<ParameterizationCylindricCurve>> res_c_curves;
		
		// with planes
		for (int i = 0, i_end = planes.size(); i < i_end; ++i) {
			std::vector<Point> points;
			std::vector<ParameterizationCircle> circles;

			get_intersections(planes[i], S1, points, circles);

			for (int ii = 0, ii_end = points.size(); ii < ii_end; ++ii) res_points.push_back(points[ii]);
			res_c_circles.push_back(circles);
		}

		// with cylinders
		for (int i = 0, i_end = cylinders.size(); i < i_end; ++i) {
			std::vector<Point> points;
			std::vector<ParameterizationCylindricCurve> curves;

			get_intersections(cylinders[i], S1, points, curves);

			for (int ii = 0, ii_end = points.size(); ii < ii_end; ++ii) res_points.push_back(points[ii]);
			res_c_curves.push_back(curves);
		}

		// with spheres
		for (int i = 0, i_end = spheres.size(); i < i_end; ++i) {
			std::vector<Point> points;
			std::vector<ParameterizationCircle> circles;

			get_intersections(spheres[i], S1, points, circles);

			for (int ii = 0, ii_end = points.size(); ii < ii_end; ++ii) res_points.push_back(points[ii]);
			res_c_circles.push_back(circles);
		}

		// process the intersections' intersections

		// curves and curves
		for (int i = 0, i_end = res_c_curves.size(); i < i_end; ++i) {
			for (int j = 0, j_end = res_c_curves.size(); j < j_end; ++j) {
				if (i == j) continue;
				get_intersections(res_c_curves[i], res_c_curves[j]);
			}
		}

		// curves and circles
		for (int i = 0, i_end = res_c_curves.size(); i < i_end; ++i) {
			for (int j = 0, j_end = res_c_circles.size(); j < j_end; ++j) {
				get_intersections(res_c_circles[j], res_c_curves[i]);
			}
		}

		// circles and circles
		for (int i = 0, i_end = res_c_circles.size(); i < i_end; ++i) {
			for (int j = 0, j_end = res_c_circles.size(); j < j_end; ++j) {
				if (i == j) continue;
				get_intersections(res_c_circles[i], res_c_circles[j]);
			}
		}
		
		// points and curves
		for (int i = 0, i_end = res_c_curves.size(); i < i_end; ++i) {
			get_intersections(res_points, res_c_curves[i]);
		}

		// points and circles
		for (int i = 0, i_end = res_c_circles.size(); i < i_end; ++i) {
			get_intersections(res_points, res_c_circles[i]);
		}

		// return
		for (int i = 0, i_end = res_c_circles.size(); i < i_end; ++i) {
			for (int j = 0, j_end = res_c_circles[i].size(); j < j_end; ++j) {
				res_circles.push_back(res_c_circles[i][j]);
			}
		}
		for (int i = 0, i_end = res_c_curves.size(); i < i_end; ++i) {
			for (int j = 0, j_end = res_c_curves[i].size(); j < j_end; ++j) {
				res_curves.push_back(res_c_curves[i][j]);
			}
		}

		return res_points.size() + res_circles.size() + res_curves.size();
	}
}