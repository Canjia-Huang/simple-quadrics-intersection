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
		while (PC1.t_lb() > 360) {
			PC1.t_lb() -= 360; PC1.t_ub() -= 360;
		}
		while (PC2.t_lb() > 360) {
			PC2.t_lb() -= 360; PC2.t_ub() -= 360;
		}
		double overlap_t_lb = std::max(PC1.t_lb(), PC2.t_lb());
		double overlap_t_ub = std::min(PC1.t_ub(), PC2.t_ub());
		if (overlap_t_lb > overlap_t_ub - SQI_EPS) { // regions dont overlap
			// SQI_VERBOSE_ONLY_COUT("regions dont overlap");
		}
		else { // regions overlap
			// SQI_VERBOSE_ONLY_COUT("regions overlap");

			int pre_status = 0, status = 0;
			double t = overlap_t_lb + SQI_EPS;

			while (t < overlap_t_ub + t_step - SQI_EPS) {

				std::vector<double> C1_ss, C2_ss;
				PC1.get_s(t, C1_ss);
				PC2.get_s(t, C2_ss);
				if (C1_ss.size() == 1 && C2_ss.size() == 1) {
					double C1_s = C1_ss[0];
					double C2_s = C2_ss[0];
					status = (C1_s > C2_s) ? 1 : -1;

					if (pre_status != 0) { // pre_status is meaningful
						if (status == pre_status) { // [pre_t, t] dont have intersection
						}
						else { // [pre_t, t] have intersection
							// SQI_VERBOSE_ONLY_COUT("get an intersection");

							// do bisection
							double t_lb = t - t_step;
							double t_ub = t;
							while (t_ub - t_lb < SQI_EPS) {
								double mid_t = 0.5 * (t_lb + t_ub);

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

				t += t_step;
			}
		}

		return sub_PC1s.size() + sub_PC2s.size();
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
		std::vector<std::vector<ParameterizationCylindricPoint>> res_points_array;
		std::vector<std::vector<ParameterizationCylindricLine>> res_lines_array;
		std::vector<std::vector<ParameterizationCylindricCurve>> res_curves_array;

		// get the intersections between cylinder C1 and other primitives
		for (int i = 0, i_end = planes.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricLine> lines;
			std::vector<ParameterizationCylindricCurve> curves;

			get_intersections(planes[i], C1, lines, curves);

			res_lines_array.push_back(lines);
			res_curves_array.push_back(curves);
		}
		for (int i = 0, i_end = cylinders.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricPoint> points;
			std::vector<ParameterizationCylindricLine> lines;
			std::vector<ParameterizationCylindricCurve> curves;

			get_intersections(cylinders[i], C1, points, lines, curves);

			res_points_array.push_back(points);
			res_lines_array.push_back(lines);
			res_curves_array.push_back(curves);
		}
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
		for (int i = 0, i_end = res_curves_array.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricCurve>* curves = &(res_curves_array[i]);

			if (res_curves.size() != 0) {
				get_intersections(*curves, res_curves);
			}

			res_curves.insert(res_curves.end(), curves->begin(), curves->end());
		}
		for (int i = 0, i_end = res_lines_array.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricLine>* lines = &(res_lines_array[i]);

			if (res_curves.size() != 0) {
				get_intersections(*lines, res_curves);
			}

			res_lines.insert(res_lines.end(), lines->begin(), lines->end());
		}

		return res_points.size() + res_lines.size() + res_curves.size();
	}
}