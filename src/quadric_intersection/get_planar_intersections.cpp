#include "simple_quadrics_intersection.h"

namespace QuadricsIntersection
{
	int get_intersections(
		Plane& P1, Plane& P2,
		Line& L
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		if (1 - std::abs(P1.nor().dot(P2.nor())) < SQI_EPS) { // parallel
			return 0;
		}

		Eigen::Vector3d L_nor = (P1.nor().cross(P2.nor())).normalized();

		Eigen::Vector3d P1_vec = (P1.nor().cross(L_nor)).normalized();
		double t = (P2.cor() - P1.cor()).dot(P2.nor()) / (P1_vec.dot(P2.nor()));

		Eigen::Vector3d L_cor = P1.cor() + t * P1_vec;

		L = Line(L_cor, L_nor);

		return 1;
	}

	int get_intersections(
		Plane& P1, Plane& P2,
		std::vector<Line>& Ls
	) {
		// init
		std::vector<Line>().swap(Ls);

		Line L;
		if (get_intersections(P1, P2, L) == 1) Ls.push_back(L);

		return Ls.size();
	}

	int get_intersections(
		Cylinder& C1, Plane& P1,
		std::vector<ParameterizationCylindricLine>& lines,
		std::vector<ParameterizationCylindricCurve>& curves
	) {
		return get_intersections(P1, C1, lines, curves);
	}

	int get_intersections(
		Sphere& S1, Plane& P1,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles
	){
		return get_intersections(P1, S1, points, circles);
	}

	int get_intersections(
		Plane& P1,
		std::vector<Plane>& planes,
		std::vector<Cylinder>& cylinders,
		std::vector<Sphere>& spheres,
		std::vector<Point>& res_points,
		std::vector<Line>& res_lines,
		std::vector<ParameterizationCircle>& res_circles,
		std::vector<std::vector<ParameterizationCylindricLine>>& res_c_lines,
		std::vector<std::vector<ParameterizationCylindricCurve>>& res_c_curves
	) {
		SQI_VERBOSE_ONLY_COUT("");

		// init
		std::vector<Point>().swap(res_points);
		std::vector<Line>().swap(res_lines);
		std::vector<ParameterizationCircle>().swap(res_circles);
		std::vector<std::vector<ParameterizationCylindricLine>>().swap(res_c_lines);
		std::vector<std::vector<ParameterizationCylindricCurve>>().swap(res_c_curves);

		// get the intersections between Plane P1 and other primitives
		std::vector<std::vector<Point>> res_points_array;
		std::vector<std::vector<Line>> res_lines_array;
		std::vector<std::vector<ParameterizationCircle>> res_circles_array;
		std::vector<std::vector<ParameterizationCylindricLine>> res_c_lines_array;
		std::vector<std::vector<ParameterizationCylindricCurve>> res_c_curves_array;
		for (int i = 0, i_end = planes.size(); i < i_end; ++i) {
			std::vector<Line> lines;

			get_intersections(planes[i], P1, lines);

			res_lines_array.push_back(lines);
		}
		for (int i = 0, i_end = spheres.size(); i < i_end; ++i) {
			std::vector<Point> points;
			std::vector<ParameterizationCircle> circles;

			get_intersections(spheres[i], P1, points, circles);

			res_points_array.push_back(points);
			res_circles_array.push_back(circles);
		}
		for (int i = 0, i_end = cylinders.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricLine> c_lines;
			std::vector<ParameterizationCylindricCurve> c_curves;

			get_intersections(cylinders[i], P1, c_lines, c_curves);

			res_c_lines_array.push_back(c_lines);
			res_c_curves_array.push_back(c_curves);
		}

		// process the intersections' intersections
		for (int i = 0, i_end = res_points_array.size(); i < i_end; ++i) {
			res_points.insert(res_points.end(), res_points_array[i].begin(), res_points_array[i].end());
		}
		for (int i = 0, i_end = res_lines_array.size(); i < i_end; ++i) {
			std::vector<Line>* lines = &res_lines_array[i];

			if (res_lines.size() > 0) {
				get_intersections(*lines, res_lines);
			}

			res_lines.insert(res_lines.end(), lines->begin(), lines->end());
		}
		for (int i = 0, i_end = res_c_lines_array.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricLine>* c_lines = &res_c_lines_array[i];

			get_intersections(res_lines, *c_lines, cylinders[i]);
		}
		res_c_lines = res_c_lines_array;
		for (int i = 0, i_end = res_circles_array.size(); i < i_end; ++i) {
			std::vector<ParameterizationCircle>* circles = &res_circles_array[i];

			get_intersections(res_lines, *circles);
			
			res_circles.insert(res_circles.end(), circles->begin(), circles->end());
		}
		for (int i = 0, i_end = res_c_curves_array.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricCurve>* curves = &res_c_curves_array[i];

			get_intersections(res_lines, *curves, cylinders[i]);
		}
		res_c_curves = res_c_curves_array;

		// return

		return 0;
	}
}