#include "simple_quadrics_intersection.h"

namespace QuadricsIntersection 
{
	int get_intersections(
		Plane& P1, Sphere& S1,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles
	) {
		SQI_VERBOSE_ONLY_COUT("");

		// init 
		std::vector<Point>().swap(points);
		std::vector<ParameterizationCircle>().swap(circles);

		Eigen::Vector3d S1_cor_to_P1_cor = P1.cor() - S1.cor();
		double P1_S1_cor_dot_P1_nor = S1_cor_to_P1_cor.dot(P1.nor());
		double P1_S1_cor_dis = std::abs(P1_S1_cor_dot_P1_nor);

		if (P1_S1_cor_dis > S1.r() + SQI_EPS) { // not intersect
			SQI_VERBOSE_ONLY_COUT("not intersect");
		}
		else {
			if (P1_S1_cor_dis > S1.r() - SQI_EPS) { // tangent
				SQI_VERBOSE_ONLY_COUT("tangent");

				points.push_back(Point(S1.cor() + P1_S1_cor_dot_P1_nor * P1.nor()));
			}
			else { // intersect
				SQI_VERBOSE_ONLY_COUT("intersect");

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
		std::vector<ParameterizationCylindricPoint>& points,
		std::vector<ParameterizationCylindricCurve>& curves
	) {
		return get_intersections(S1, C1, points, curves);
	}

	int get_intersections(
		Sphere& S1, Sphere& S2,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles
	) {
		SQI_VERBOSE_ONLY_COUT("");

		// init
		std::vector<Point>().swap(points);
		std::vector<ParameterizationCircle>().swap(circles);

		Eigen::Vector3d S2_cor_to_S1_cor = S1.cor() - S2.cor();
		double centers_dis = S2_cor_to_S1_cor.norm();

		if (centers_dis < SQI_EPS) { // overlap: result nothing
			SQI_VERBOSE_ONLY_COUT("overlap");
		}
		else {
			Eigen::Vector3d S2_cor_to_S1_cor_unit = S2_cor_to_S1_cor / centers_dis;

			if (centers_dis > S1.r() + S2.r() + SQI_EPS) { // not intersect
				SQI_VERBOSE_ONLY_COUT("not intersect");
			}
			else if (centers_dis > S1.r() + S2.r() - SQI_EPS) { // tangent: result a point
				SQI_VERBOSE_ONLY_COUT("tangent");

				points.push_back(Point(S2.cor() + S2.r() * S2_cor_to_S1_cor_unit));
			}
			else { // intersect: result a circle
				SQI_VERBOSE_ONLY_COUT("intersect");

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
		std::vector<Plane>& planes,
		std::vector<Cylinder>& cylinders,
		std::vector<Sphere>& spheres,
		std::vector<Point>& res_points,
		std::vector<ParameterizationCircle>& res_circles,
		std::vector<std::vector<ParameterizationCylindricPoint>>& res_c_points,
		std::vector<std::vector<ParameterizationCylindricCurve>>& res_c_curves
	) {
		SQI_VERBOSE_ONLY_COUT("");

		// init
		std::vector<Point>().swap(res_points);
		std::vector<ParameterizationCircle>().swap(res_circles);
		std::vector<std::vector<ParameterizationCylindricPoint>>().swap(res_c_points);
		std::vector<std::vector<ParameterizationCylindricCurve>>().swap(res_c_curves);

		// get the intersections between sphere S1 and other primitives
		std::vector<std::vector<Point>> res_points_array;
		std::vector<std::vector<ParameterizationCircle>> res_circles_array;
		std::vector<std::vector<ParameterizationCylindricPoint>> res_c_points_array;
		std::vector<std::vector<ParameterizationCylindricCurve>> res_c_curves_array;
		for (int i = 0, i_end = planes.size(); i < i_end; ++i) {
			std::vector<Point> points;
			std::vector<ParameterizationCircle> circles;

			get_intersections(planes[i], S1, points, circles);

			res_points_array.push_back(points);
			res_circles_array.push_back(circles);
		}
		for (int i = 0, i_end = spheres.size(); i < i_end; ++i) {
			std::vector<Point> points;
			std::vector<ParameterizationCircle> circles;

			get_intersections(spheres[i], S1, points, circles);

			res_points_array.push_back(points);
			res_circles_array.push_back(circles);
		}
		for (int i = 0, i_end = cylinders.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricPoint> c_points;
			std::vector<ParameterizationCylindricCurve> c_curves;

			get_intersections(cylinders[i], S1, c_points, c_curves);

			res_c_points_array.push_back(c_points);
			res_c_curves_array.push_back(c_curves);
		}

		// process the intersections' intersections
		for (int i = 0, i_end = res_points_array.size(); i < i_end; ++i) {
			res_points.insert(res_points.end(), res_points_array[i].begin(), res_points_array[i].end());
		}
		for (int i = 0, i_end = res_circles_array.size(); i < i_end; ++i) {
			std::vector<ParameterizationCircle>* circles = &res_circles_array[i];

			if (res_circles.size() > 0) {
				get_intersections(*circles, res_circles);
			}

			res_circles.insert(res_circles.end(), circles->begin(), circles->end());
		}
		res_c_points = res_c_points_array;
		for (int i = 0, i_end = res_c_curves_array.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricCurve>* c_curves = &res_c_curves_array[i];

			if (c_curves->size() > 0) {
				get_intersections(res_circles, *c_curves, cylinders[i]);
			}
		}
		res_c_curves = res_c_curves_array;

		// return
		int total_num = res_points.size() + res_circles.size();
		for (auto ps : res_c_points) total_num += ps.size();
		for (auto cs : res_c_curves) total_num += cs.size();

		return total_num;
	}
}