#include "simple_quadrics_intersection.h"

namespace QuadricsIntersection
{
	int get_intersections(
		ParameterizationCylindricLine& L1, ParameterizationCylindricCurve& PC1,
		std::vector<ParameterizationCylindricLine>& sub_L1s,
		std::vector<ParameterizationCylindricCurve>& sub_PC1s
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		// init
		std::vector<ParameterizationCylindricLine>().swap(sub_L1s);
		std::vector<ParameterizationCylindricCurve>().swap(sub_PC1s);

		if (L1.t() < PC1.t_lb() + SQI_EPS || L1.t() > PC1.t_ub() - SQI_EPS) { // not intersect
			// SQI_VERBOSE_ONLY_COUT("not intersect");
		}
		else { // may intersect
			std::vector<double> ss;
			if (PC1.get_s(L1.t(), ss) == 1) {
				if (ss[0] < L1.s_lb() || ss[0] > L1.s_ub()) { // not intersect
					// SQI_VERBOSE_ONLY_COUT("not intersect");
				}
				else { // intersect
					// SQI_VERBOSE_ONLY_COUT("intersect");

					// cut the curve
					ParameterizationCylindricCurve sub_PC1;
					PC1.separate_t(sub_PC1, L1.t());

					sub_PC1s.push_back(sub_PC1);

					// cut the line
					ParameterizationCylindricLine sub_L1;
					L1.separate_s(sub_L1, ss[0]);

					sub_L1s.push_back(sub_L1);
				}
			}
			else {
				SQI_VERBOSE_ONLY_WARNING("input curves may not monotonicity!");
				return 0;
			}
		}

		return sub_L1s.size() + sub_PC1s.size();
	}

	int get_intersections(
		std::vector<ParameterizationCylindricLine>& L1s,
		std::vector<ParameterizationCylindricCurve>& PC1s
	) {
		SQI_VERBOSE_ONLY_COUT("");

		for (int i = 0; i < L1s.size(); ++i) { // note: the vector L1s's size is dymanic, so cannot use i_end = L1s.size.
			for (int j = 0, j_end = PC1s.size(); j < j_end; ++j) { // note: the vector PC1s's size is dymanic, but the newly added curves do not need to be judged again.
				ParameterizationCylindricLine* L1 = &L1s[i];
				ParameterizationCylindricCurve* PC1 = &PC1s[j];

				std::vector<ParameterizationCylindricLine> sub_L1s;
				std::vector<ParameterizationCylindricCurve> sub_PC1s;
				if (get_intersections(*L1, *PC1, sub_L1s, sub_PC1s) > 0) {

					for (auto l : sub_L1s) L1s.push_back(l);
					for (auto pc : sub_PC1s) PC1s.push_back(pc);
				}
			}
		}

		return L1s.size() + PC1s.size();
	}

	// --------------------------------------------------------------------------------------------

	int get_intersections(
		ParameterizationCylindricCurve& PC1, ParameterizationCylindricCurve& PC2,
		std::vector<ParameterizationCylindricCurve>& sub_PC1s,
		std::vector<ParameterizationCylindricCurve>& sub_PC2s,
		double t_step
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		// init
		std::vector<ParameterizationCylindricCurve>().swap(sub_PC1s);
		std::vector<ParameterizationCylindricCurve>().swap(sub_PC2s);

		// get overlap parameterization region
		double overlap_t_lb = std::max(PC1.t_lb(), PC2.t_lb());
		double overlap_t_ub = std::min(PC1.t_ub(), PC2.t_ub());
		if (overlap_t_lb > overlap_t_ub - SQI_EPS) { // regions dont overlap
			// SQI_VERBOSE_ONLY_COUT("regions dont overlap");
			// do nothing
		}
		else { // regions overlap
			// SQI_VERBOSE_ONLY_COUT("regions overlap");
			int pre_status = 0, status = 0;
			double pre_t;
			double t = overlap_t_lb + SQI_EPS;

			while (t < overlap_t_ub + t_step - SQI_EPS) {
				std::vector<double> C1_ss, C2_ss;
				PC1.get_s(t, C1_ss);
				PC2.get_s(t, C2_ss);
				if (C1_ss.size() == 1 && C2_ss.size() == 1) {
					status = (C1_ss[0] > C2_ss[0]) ? 1 : -1;

					if (pre_status != 0) { // pre_status is meaningful

						if (status == pre_status) { // [pre_t, t] dont have intersection
							// do nothing
						}
						else { // [pre_t, t] have intersection
							// SQI_VERBOSE_ONLY_COUT("get an intersection");

							// do bisection
							double t_lb = pre_t;
							double t_ub = (t > overlap_t_ub) ? (overlap_t_ub) : t;

							while (t_ub - t_lb > SQI_EPS) {
								double mid_t = 0.5 * (t_lb + t_ub);
								// SQI_VERBOSE_ONLY_TEST(pre_t << " " << mid_t << " " << t);

								std::vector<double> ss1, ss2;
								if (PC1.get_s(mid_t, ss1) == 1 && PC2.get_s(mid_t, ss2) == 1) {
									int cur_status = (ss1[0] > ss2[0]) ? 1 : -1;

									if (cur_status == pre_status) t_lb = mid_t;
									else t_ub = mid_t;
								}
								else {
									SQI_VERBOSE_ONLY_WARNING("bisection may error!");
									return 0;
								}
							}

							// cut the curves
							ParameterizationCylindricCurve sub_PC1, sub_PC2;
							if (PC1.separate_t(sub_PC1, t_lb) == 1) sub_PC1s.push_back(sub_PC1);
							if (PC2.separate_t(sub_PC2, t_lb) == 1) sub_PC2s.push_back(sub_PC2);
						}
					}

					pre_status = status;
				}
				else if (C1_ss.size() == 0 || C2_ss.size() == 0) { // may out of region
				}
				else { // may error before
					SQI_VERBOSE_ONLY_WARNING("input curves may not monotonicity!");
					return 0;
				}

				pre_t = t;
				t += t_step;
			}
		}

		return sub_PC1s.size() + sub_PC2s.size();
	}

	int get_intersections(
		std::vector<ParameterizationCylindricCurve>& PC1s,
		std::vector<ParameterizationCylindricCurve>& PC2s,
		double t_step
	) {
		SQI_VERBOSE_ONLY_COUT("");

		for (int i = 0; i < PC1s.size(); ++i) { // note: the vector PC1s's size is dymanic, so cannot use i_end = PC1s.size.
			for (int j = 0, j_end = PC2s.size(); j < j_end; ++j) { // note: the vector PC2s's size is dymanic, but the newly added curves do not need to be judged again.								
				ParameterizationCylindricCurve* PC1 = &PC1s[i];
				ParameterizationCylindricCurve* PC2 = &PC2s[j];

				std::vector<ParameterizationCylindricCurve> sub_PC1s;
				std::vector<ParameterizationCylindricCurve> sub_PC2s;
				if (get_intersections(*PC1, *PC2, sub_PC1s, sub_PC2s, t_step) > 0) {
					for (auto pc : sub_PC1s) PC1s.push_back(pc);
					for (auto pc : sub_PC2s) PC2s.push_back(pc);
				}
			}
		}

		return PC1s.size() + PC2s.size();
	}

	// --------------------------------------------------------------------------------------------

	int get_intersections(
		ParameterizationCircle& c1, ParameterizationCylindricCurve& PC1, Cylinder& C1,
		std::vector<ParameterizationCircle>& sub_c1s,
		std::vector<ParameterizationCylindricCurve>& sub_PC1s
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		// init 
		std::vector<ParameterizationCircle>().swap(sub_c1s);
		std::vector<ParameterizationCylindricCurve>().swap(sub_PC1s);

		Plane P1(c1.cor(), c1.nor()); // proxy plane that is co-faceted with the circle

		std::vector<ParameterizationCylindricLine> lines;
		std::vector<ParameterizationCylindricCurve> curves;
		if (get_intersections(P1, C1, lines, curves) == 0) {
			SQI_VERBOSE_ONLY_WARNING("circle and curve are not on a sphere!");
			return 0;
		}

		if (lines.size() > 0) { // circle is parallel to cylinder's axis
			SQI_VERBOSE_ONLY_COUT("parallel to cylinder's axis");

			std::vector<ParameterizationCylindricCurve> tmp_PC1s = { PC1 };
			get_intersections(lines, tmp_PC1s);
			if (tmp_PC1s.size() > 1) {
				std::vector<Eigen::Vector3d> intersect_points;
				std::vector<double> intersect_cts;
				std::vector<double> intersect_Cts;

				for (int i = 0, i_end = tmp_PC1s.size(); i < i_end; ++i) {
					double t = tmp_PC1s[i].t_ub();

					if (std::abs(t - PC1.t_ub()) < SQI_EPS) continue;

					std::vector<double> ss;
					PC1.get_s(t, ss);
					if (ss.size() != 1) {
						SQI_VERBOSE_ONLY_COUT("may error!");
						return 0;
					}

					Eigen::Vector3d intersect_point = C1.get_point(ss[0], t);

					double ct = c1.get_t(intersect_point);

					if (c1.is_t_valid(ct)) {
						intersect_points.push_back(intersect_point);
						intersect_cts.push_back(ct);
						intersect_Cts.push_back(t);
					}
				}

				// cut the curve
				for (double t : intersect_Cts) {
					ParameterizationCylindricCurve sub_PC1;
					PC1.separate_t(sub_PC1, t);
					sub_PC1s.push_back(sub_PC1);
				}

				// cut the circle
				c1.separate_t(sub_c1s, intersect_cts);
			}
		}
		else if (curves.size() == 1) {
			SQI_VERBOSE_ONLY_COUT("circle is intersect");

			ParameterizationCylindricCurve tmp_PC1 = PC1;
			std::vector<ParameterizationCylindricCurve> sub_tmp_PC1s;
			std::vector<ParameterizationCylindricCurve> sub_tmp_PC2s;
			if (get_intersections(tmp_PC1, curves[0], sub_tmp_PC1s, sub_tmp_PC2s) > 0) { // have intersections, but may not in actual
				std::vector<Eigen::Vector3d> intersect_points;
				std::vector<double> intersect_cts;
				std::vector<double> intersect_Cts;

				for (int i = 0, i_end = sub_tmp_PC1s.size(); i < i_end; ++i) {
					double t = sub_tmp_PC1s[i].t_ub();

					std::vector<double> ss;
					PC1.get_s(t, ss);
					if (ss.size() != 1) {
						SQI_VERBOSE_ONLY_COUT("may error!");
						return 0;
					}

					Eigen::Vector3d intersect_point = C1.get_point(ss[0], t);

					double ct = c1.get_t(intersect_point);

					if (std::abs((intersect_point - c1.cor()).norm() - c1.r())> SQI_EPS) {
						continue;
					}

					if (c1.is_t_valid(ct)) {
						intersect_points.push_back(intersect_point);
						intersect_cts.push_back(ct);
						intersect_Cts.push_back(t);
					}
				}

				// cut the curve
				for (double t : intersect_Cts) {
					ParameterizationCylindricCurve sub_PC1;
					PC1.separate_t(sub_PC1, t);
					sub_PC1s.push_back(sub_PC1);
				}

				// cut the circle
				c1.separate_t(sub_c1s, intersect_cts);
			}
		}
		else {
			SQI_VERBOSE_ONLY_WARNING("may error!");
			return 0;
		}

		return sub_c1s.size() + sub_PC1s.size();
	}

	int get_intersections(
		std::vector<ParameterizationCircle>& c1s,
		std::vector<ParameterizationCylindricCurve>& PC1s, Cylinder& C1
	) {
		SQI_VERBOSE_ONLY_COUT("");

		for (int i = 0; i < c1s.size(); ++i) { // note: the vector c1s's size is dymanic, so cannot use i_end = c1s.size.
			for (int j = 0, j_end = PC1s.size(); j < j_end; ++j) { // note: the vector PC1s's size is dymanic, but the newly added curves do not need to be judged again.				
				ParameterizationCircle* c1 = &c1s[i];
				ParameterizationCylindricCurve* PC1 = &PC1s[j];

				std::vector<ParameterizationCircle> sub_c1s;
				std::vector<ParameterizationCylindricCurve> sub_PC1s;
				if (get_intersections(*c1, *PC1, C1, sub_c1s, sub_PC1s) > 0) {
					for (auto c : sub_c1s) c1s.push_back(c);
					for (auto pc : sub_PC1s) PC1s.push_back(pc);
				}
			}
		}

		return c1s.size() + PC1s.size();
	}

	// --------------------------------------------------------------------------------------------

	int get_intersections(
		ParameterizationCircle& c1, ParameterizationCircle& c2,
		std::vector<ParameterizationCircle>& sub_c1s,
		std::vector<ParameterizationCircle>& sub_c2s
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		// init
		std::vector<ParameterizationCircle>().swap(sub_c1s);
		std::vector<ParameterizationCircle>().swap(sub_c2s);

		if (1 - std::abs(c1.nor().dot(c2.nor())) < SQI_EPS) { // two circle's nors are parallel: result nothing
			SQI_VERBOSE_ONLY_COUT("parallel");
			return 0;
		}

		Plane P1(c1.cor(), c1.nor());
		Plane P2(c2.cor(), c2.nor());
		Line intersect_L;
		get_intersections(P1, P2, intersect_L);

		// find intersections on Circle c1, solve a * t^2 + b * t + c = 0
		Eigen::Vector3d c1_cor_to_L_cor = intersect_L.cor() - c1.cor();
		double a = intersect_L.nor().squaredNorm();
		double b = 2 * c1_cor_to_L_cor.dot(intersect_L.nor());
		double c = c1_cor_to_L_cor.squaredNorm() - c1.r() * c1.r();
		double delta = b * b - 4 * a * c;
		if (delta < -SQI_EPS) { // not intersect
			SQI_VERBOSE_ONLY_COUT("not intersect");
			return 0;
		}
		else if (delta < SQI_EPS) { // tangent: 1 point
			SQI_VERBOSE_ONLY_COUT("tangent");

			double intersect_t = -b / (2 * a);
			Eigen::Vector3d intersect_point = intersect_L.cor() + intersect_t * intersect_L.nor();

			if (std::abs((intersect_point - c2.cor()).squaredNorm() - c2.r() * c2.r()) < SQI_EPS) { // intersect point is on c2
				double c1_t = c1.get_t(intersect_point);
				if (c1.is_t_valid(c1_t) == false) return 0;

				double c2_t = c2.get_t(intersect_point);
				if (c2.is_t_valid(c2_t) == false) return 0;

				// cut the circles
				ParameterizationCircle sub_c1, sub_c2;
				c1.separate_t(sub_c1, c1_t);
				c2.separate_t(sub_c2, c2_t);
				sub_c1s.push_back(sub_c1);
				sub_c2s.push_back(sub_c2);
			}
		}
		else { // intersect: 2 point
			SQI_VERBOSE_ONLY_COUT("intersect");

			double sqrt_delta = std::sqrt(delta);
			double intersect_t1 = (-b + sqrt_delta) / (2 * a);
			double intersect_t2 = (-b - sqrt_delta) / (2 * a);
			Eigen::Vector3d intersect_point1 = intersect_L.cor() + intersect_t1 * intersect_L.nor();
			Eigen::Vector3d intersect_point2 = intersect_L.cor() + intersect_t2 * intersect_L.nor();

			std::vector<double> cutting_c1ts;
			std::vector<double> cutting_c2ts;

			if (std::abs((intersect_point1 - c2.cor()).squaredNorm() - c2.r() * c2.r()) < SQI_EPS) { // intersect point is on c2
				double c1_t = c1.get_t(intersect_point1);
				if (c1.is_t_valid(c1_t) == false) return 0;

				double c2_t = c2.get_t(intersect_point1);
				if (c2.is_t_valid(c2_t) == false) return 0;

				cutting_c1ts.push_back(c1_t);
				cutting_c2ts.push_back(c2_t);
			}

			if (std::abs((intersect_point2 - c2.cor()).squaredNorm() - c2.r() * c2.r()) < SQI_EPS) { // intersect point is on c2
				double c1_t = c1.get_t(intersect_point2);
				if (c1.is_t_valid(c1_t) == false) return 0;

				double c2_t = c2.get_t(intersect_point2);
				if (c2.is_t_valid(c2_t) == false) return 0;

				cutting_c1ts.push_back(c1_t);
				cutting_c2ts.push_back(c2_t);
			}

			// cut!
			c1.separate_t(sub_c1s, cutting_c1ts);
			c2.separate_t(sub_c2s, cutting_c2ts);
		}

		return sub_c1s.size() + sub_c2s.size();
	}

	int get_intersections(
		std::vector<ParameterizationCircle>& c1s,
		std::vector<ParameterizationCircle>& c2s
	) {
		SQI_VERBOSE_ONLY_COUT("");

		for (int i = 0; i < c1s.size(); ++i) { // note: the vector c1s's size is dymanic, so cannot use i_end = c1s.size.
			for (int j = 0, j_end = c2s.size(); j < j_end; ++j) { // note: the vector c2s's size is dymanic, but the newly added circles do not need to be judged again.				
				ParameterizationCircle* c1 = &c1s[i];
				ParameterizationCircle* c2 = &c2s[j];

				std::vector<ParameterizationCircle> sub_c1s;
				std::vector<ParameterizationCircle> sub_c2s;
				if (get_intersections(*c1, *c2, sub_c1s, sub_c2s) > 0) {
					for (auto c : sub_c1s) c1s.push_back(c);
					for (auto c : sub_c2s) c2s.push_back(c);
				}
			}
		}

		return c1s.size() + c2s.size();
	}

	// --------------------------------------------------------------------------------------------





	int get_intersections(
		Cylinder& C1,
		std::vector<Plane>& planes,
		std::vector<Cylinder>& cylinders,
		std::vector<Sphere>& spheres,
		std::vector<ParameterizationCylindricPoint>& res_points,
		std::vector<ParameterizationCylindricLine>& res_lines,
		std::vector<ParameterizationCylindricCurve>& res_curves
	) {
		SQI_VERBOSE_ONLY_COUT("");

		// init
		std::vector<ParameterizationCylindricPoint>().swap(res_points);
		std::vector<ParameterizationCylindricLine>().swap(res_lines);
		std::vector<ParameterizationCylindricCurve>().swap(res_curves);

		// get the intersections between cylinder C1 and other primitives
		std::vector<std::vector<ParameterizationCylindricPoint>> res_points_array;
		std::vector<std::vector<ParameterizationCylindricLine>> res_lines_array;
		std::vector<std::vector<ParameterizationCylindricCurve>> res_curves_array;
		for (int i = 0, i_end = planes.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricLine> lines;
			std::vector<ParameterizationCylindricCurve> curves;

			get_intersections(planes[i], C1, lines, curves);

			res_lines_array.push_back(lines);
			res_curves_array.push_back(curves);
		}
		// SQI_VERBOSE_ONLY_TEST("process cylinders");
		for (int i = 0, i_end = cylinders.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricPoint> points;
			std::vector<ParameterizationCylindricLine> lines;
			std::vector<ParameterizationCylindricCurve> curves;

			get_intersections(cylinders[i], C1, points, lines, curves);

			res_points_array.push_back(points);
			res_lines_array.push_back(lines);
			res_curves_array.push_back(curves);
		}
		// SQI_VERBOSE_ONLY_TEST("process spheres");
		for (int i = 0, i_end = spheres.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricPoint> points;
			std::vector<ParameterizationCylindricCurve> curves;

			get_intersections(spheres[i], C1, points, curves);

			res_points_array.push_back(points);
			res_curves_array.push_back(curves);
		}		

		// process the intersections' intersections
		for (int i = 0, i_end = res_points_array.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricPoint>* points = &(res_points_array[i]);

			res_points.insert(res_points.end(), points->begin(), points->end());
		}
		SQI_VERBOSE_ONLY_TEST("process cylindric curves");
		for (int i = 0, i_end = res_curves_array.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricCurve>* curves = &(res_curves_array[i]);

			if (res_curves.size() > 0) {
				std::cout << get_intersections(*curves, res_curves) << std::endl;
			}

			res_curves.insert(res_curves.end(), curves->begin(), curves->end());
		}
		SQI_VERBOSE_ONLY_TEST("process cylindric line");
		for (int i = 0, i_end = res_lines_array.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricLine>* lines = &(res_lines_array[i]);

			if (res_curves.size() > 0) {
				get_intersections(*lines, res_curves);
			}

			res_lines.insert(res_lines.end(), lines->begin(), lines->end());
		}

		return res_points.size() + res_lines.size() + res_curves.size();
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
		std::vector<std::vector<ParameterizationCylindricCurve>> res_c_circles_array;
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
			res_c_circles_array.push_back(c_curves);
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
		for (int i = 0, i_end = res_c_circles_array.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricCurve>* c_curves = &res_c_circles_array[i];

			if (c_curves->size() > 0) {
				get_intersections(res_circles, *c_curves, cylinders[i]);
			}
		}
		res_c_curves = res_c_circles_array;

		// return
		int total_num = res_points.size() + res_circles.size();
		for (auto ps : res_c_points) total_num += ps.size();
		for (auto cs : res_c_curves) total_num += cs.size();

		return total_num;
	}
}