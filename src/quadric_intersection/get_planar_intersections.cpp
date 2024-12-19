#include "simple_quadrics_intersection.h"

namespace QuadricsIntersection
{
	int get_intersections(
		Plane& P1, Plane& P2,
		Line& L,
		double limit_angle
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between a Plane and a Plane");

		double P1_nor_dor_P2_nor = P1.nor().dot(P2.nor());
		double P1_nor_P2_nor_angle = safetyAcos(P1_nor_dor_P2_nor);
		
		if (P1_nor_P2_nor_angle < (limit_angle + SQI_EPS) ||
			P1_nor_P2_nor_angle > (180 - limit_angle - SQI_EPS)) { // parallel
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
		std::vector<Line>& Ls,
		double limit_angle
	) {
		// init
		std::vector<Line>().swap(Ls);

		Line L;
		if (get_intersections(P1, P2, L, limit_angle) == 1) Ls.push_back(L);

		return Ls.size();
	}

	int get_intersections(
		Cylinder& C1, Plane& P1,
		std::vector<Line>& lines,
		std::vector<ParameterizationCylindricCurve>& curves,
		double limit_angle
	) {
		return get_intersections(P1, C1, lines, curves, limit_angle);
	}

	int get_intersections(
		Sphere& S1, Plane& P1,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles,
		double limit_angle
	){
		return get_intersections(P1, S1, points, circles, limit_angle);
	}

	int get_intersections(
		Plane& P1,
		std::vector<Plane>& planes, std::vector<Cylinder>& cylinders, std::vector<Sphere>& spheres,
		std::vector<Point>& res_points, std::vector<Line>& res_lines, std::vector<ParameterizationCircle>& res_circles, std::vector<ParameterizationCylindricCurve>& res_curves
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between a Plane and other primitives");

		// init
		std::vector<Point>().swap(res_points);
		std::vector<Line>().swap(res_lines);
		std::vector<ParameterizationCircle>().swap(res_circles);
		std::vector<ParameterizationCylindricCurve>().swap(res_curves);

		// get the intersections between Plane P1 and other primitives
		std::vector<std::vector<Line>> res_c_lines;
		std::vector<std::vector<ParameterizationCircle>> res_c_circles;
		std::vector<std::vector<ParameterizationCylindricCurve>> res_c_curves;

		// with planes
		for (int i = 0, i_end = planes.size(); i < i_end; ++i) {
			std::vector<Line> lines;
			get_intersections(planes[i], P1, lines);
			if (lines.size() > 0) res_c_lines.push_back(lines);
		}

		// with cylinders
		for (int i = 0, i_end = cylinders.size(); i < i_end; ++i) {
			std::vector<Line> lines;
			std::vector<ParameterizationCylindricCurve> curves;

			get_intersections(cylinders[i], P1, lines, curves);

			if (lines.size() > 0) res_c_lines.push_back(lines);
			res_c_curves.push_back(curves);
		}
		
		// with spheres
		for (int i = 0, i_end = spheres.size(); i < i_end; ++i) {
			std::vector<Point> points;
			std::vector<ParameterizationCircle> circles;

			get_intersections(spheres[i], P1, points, circles);

			for (int ii = 0, ii_end = points.size(); ii < ii_end; ++ii) res_points.push_back(points[ii]);
			if (circles.size() > 0) res_c_circles.push_back(circles);
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

		// curves and lines
		for (int i = 0, i_end = res_c_curves.size(); i < i_end; ++i) {
			for (int j = 0, j_end = res_c_lines.size(); j < j_end; ++j) {
				get_intersections(res_c_lines[j], res_c_curves[i]);
			}
		}

		// circles and circles
		for (int i = 0, i_end = res_c_circles.size(); i < i_end; ++i) {
			for (int j = 0, j_end = res_c_circles.size(); j < j_end; ++j) {
				if (i == j) continue;
				get_intersections(res_c_circles[i], res_c_circles[j]);
			}
		}

		// circles and lines
		for (int i = 0, i_end = res_c_lines.size(); i < i_end; ++i) {
			for (int j = 0, j_end = res_c_circles.size(); j < j_end; ++j) {
				get_intersections(res_c_lines[i], res_c_circles[j]);
			}
		}

		// lines and lines
		for (int i = 0, i_end = res_c_lines.size(); i < i_end; ++i) {
			for (int j = 0, j_end = res_c_lines.size(); j < j_end; ++j) {
				if (i == j) continue;
				get_intersections(res_c_lines[i], res_c_lines[j]);
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

		// points and lines
		for (int i = 0, i_end = res_c_lines.size(); i < i_end; ++i) {
			get_intersections(res_points, res_c_lines[i]);
		}

		// return
		for (int i = 0, i_end = res_c_lines.size(); i < i_end; ++i) {
			for (int j = 0, j_end = res_c_lines[i].size(); j < j_end; ++j) {
				res_lines.push_back(res_c_lines[i][j]);
			}
		}
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

		return res_points.size() + res_lines.size() + res_circles.size() + res_curves.size();
	}
}