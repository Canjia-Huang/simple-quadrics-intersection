#include "simple_quadrics_intersection.h"
#include <unordered_set>

// dont show the VERBOSE info
#define SQI_VERBOSE_ONLY_TITLE(x)
#define SQI_VERBOSE_ONLY_COUT(x)

namespace QuadricsIntersection
{
	// Point and Point (dont need to)

	// Point and Line

	int get_intersections(
		Point& P1, Line& L1,
		std::vector<Line>& sub_L1s
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between a Point and a Line");

		// init
		std::vector<Line>().swap(sub_L1s);

		Eigen::Vector3d L1_cor_to_P1 = P1.cor() - L1.cor();
		double L1_cor_to_P1_dot_L1_nor = L1_cor_to_P1.dot(L1.nor());
		double P1_L1_sq_dis = L1_cor_to_P1.squaredNorm() - L1_cor_to_P1_dot_L1_nor * L1_cor_to_P1_dot_L1_nor;

		if (P1_L1_sq_dis < SQI_EPS) { // on the line
			SQI_VERBOSE_ONLY_COUT("point -- on -- line");

			double s = L1_cor_to_P1_dot_L1_nor;
			if (L1.is_s_valid(s)) { // cut
				Line sub_L1;
				if (L1.separate_s(sub_L1, s) == 1) sub_L1s.push_back(sub_L1);
			}
		}
		else { // not on the line
			SQI_VERBOSE_ONLY_COUT("point -- not on -- line");
		}

		return sub_L1s.size();
	}

	int get_intersections(
		std::vector<Point>& P1s,
		std::vector<Line>& L1s
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between Points and Lines");

		for (int i = 0; i < P1s.size(); ++i) { // note: this vector's size is dymanic, so cannot use i_end = vector.size.
			for (int j = 0, j_end = L1s.size(); j < j_end; ++j) { // note: this vector's size is dymanic, but the newly added curves do not need to be judged again.
				Point* P1 = &P1s[i];
				Line* L1 = &L1s[j];

				std::vector<Line> sub_L1s;
				if (get_intersections(*P1, *L1, sub_L1s) > 0) {
					for (int ii = 0, ii_end = sub_L1s.size(); ii < ii_end; ++ii) L1s.push_back(sub_L1s[ii]);
				}
			}
		}

		return P1s.size() + L1s.size();
	}

	// Point and Circle

	int get_intersections(
		Point& P1, ParameterizationCircle& c1,
		std::vector<ParameterizationCircle>& sub_c1s
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between a Point and a Circle");

		// init
		std::vector<ParameterizationCircle>().swap(sub_c1s);

		Eigen::Vector3d c1_cor_to_P1 = P1.cor() - c1.cor();

		if (std::abs(c1_cor_to_P1.squaredNorm() - c1.r() * c1.r()) < SQI_EPS &&
			std::abs(c1_cor_to_P1.normalized().dot(c1.nor())) < SQI_EPS) { // on the circle
			SQI_VERBOSE_ONLY_COUT("point -- on -- circle");

			double t = c1.get_t(P1.cor());
			if (c1.is_t_valid(t)) { // cut
				ParameterizationCircle sub_c1;
				if (c1.separate_t(sub_c1, t) == 1) sub_c1s.push_back(sub_c1);
			}
		}

		return sub_c1s.size();
	}

	int get_intersections(
		std::vector<Point>& P1s,
		std::vector<ParameterizationCircle>& c1s
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between Points and Circles");

		for (int i = 0; i < P1s.size(); ++i) { // note: this vector's size is dymanic, so cannot use i_end = vector.size.
			for (int j = 0, j_end = c1s.size(); j < j_end; ++j) { // note: this vector's size is dymanic, but the newly added curves do not need to be judged again.
				Point* P1 = &P1s[i];
				ParameterizationCircle* c1 = &c1s[j];

				std::vector<ParameterizationCircle> sub_c1s;
				if (get_intersections(*P1, *c1, sub_c1s) > 0) {
					for (int ii = 0, ii_end = sub_c1s.size(); ii < ii_end; ++ii) c1s.push_back(sub_c1s[ii]);
				}
			}
		}

		return P1s.size() + c1s.size();
	}

	// Point and CylindricCurve

	int get_intersections(
		Point& P1, ParameterizationCylindricCurve& PC1, Cylinder& C1,
		std::vector<ParameterizationCylindricCurve>& sub_PC1s
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between a Point and a CylindricCurve");

		// init
		std::vector<ParameterizationCylindricCurve>().swap(sub_PC1s);

		Eigen::Vector3d C1_cor_to_P1 = P1.cor() - C1.cor();
		double C1_cor_to_P1_dot_C1_nor = C1_cor_to_P1.dot(C1.nor());
		double P1_C1_sq_dis = C1_cor_to_P1.squaredNorm() - C1_cor_to_P1_dot_C1_nor * C1_cor_to_P1_dot_C1_nor;
		
		if (std::abs(P1_C1_sq_dis - C1.r() * C1.r()) < SQI_EPS) { // point on the cylinder
			SQI_VERBOSE_ONLY_COUT("point -- on -- cylinder");

			double s, t;
			C1.get_s_t(P1.cor(), s, t);

			if (PC1.is_t_valid(t)) {
				std::vector<double> pc_ss;
				PC1.get_s(t, pc_ss);

				if (pc_ss.size() != 1) {
					SQI_VERBOSE_ONLY_WARNING("may error!");
					return 0;
				}

				if (std::abs(pc_ss[0] - s) < SQI_EPS) { // intersect
					ParameterizationCylindricCurve sub_PC1;
					if (PC1.separate_t(sub_PC1, t) == 1) sub_PC1s.push_back(sub_PC1);
				}
			}
		}
		else { // point not on the cylinder
			SQI_VERBOSE_ONLY_COUT("point -- not on -- cylinder");
			// do nothing
		}

		return sub_PC1s.size();
	}

	int get_intersections(
		std::vector<Point>& P1s,
		std::vector<ParameterizationCylindricCurve>& PC1s, Cylinder& C1
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between Points and CylindricCurves");

		for (int i = 0; i < P1s.size(); ++i) { // note: this vector's size is dymanic, so cannot use i_end = vector.size.
			for (int j = 0, j_end = PC1s.size(); j < j_end; ++j) { // note: this vector's size is dymanic, but the newly added curves do not need to be judged again.
				Point* P1 = &P1s[i];
				ParameterizationCylindricCurve* PC1 = &PC1s[j];

				std::vector<ParameterizationCylindricCurve> sub_PC1s;
				if (get_intersections(*P1, *PC1, C1, sub_PC1s) > 0) {
					for (int ii = 0, ii_end = sub_PC1s.size(); ii < ii_end; ++ii) PC1s.push_back(sub_PC1s[ii]);
				}
			}
		}

		return P1s.size() + PC1s.size();
	}

	// Line and Line

	int get_intersections(
		Line& L1, Line& L2,
		std::vector<Line>& sub_L1s,
		std::vector<Line>& sub_L2s
	) {
		// init
		std::vector<Line>().swap(sub_L1s);
		std::vector<Line>().swap(sub_L2s);

		Eigen::Vector3d cross_vec = L1.nor().cross(L2.nor());
		double sin_angle = cross_vec.norm();
		cross_vec /= sin_angle; // normalize
		double L1_L2_dis = (L2.cor() - L1.cor()).dot(cross_vec);

		if (L1_L2_dis > SQI_EPS) { // not intersect
			// do nothing
		}
		else {
			if (1 - std::abs(L2.nor().dot(L1.nor())) < SQI_EPS) { // two lines are parallel
				SQI_VERBOSE_ONLY_COUT("may error");
				// do nothing
			}
			else { // intersect a point
				Eigen::Vector3d L2_cor_to_L1_cor = L1.cor() - L2.cor();
				double L2_cor_to_L1_cor_dot_L2_nor = L2_cor_to_L1_cor.dot(L2.nor());
				Eigen::Vector3d L1_cor_proj_L2 = L2.cor() + L2_cor_to_L1_cor_dot_L2_nor * L2.nor();
				Eigen::Vector3d L1_cor_to_L1_cor_proj_L2 = L1_cor_proj_L2 - L1.cor();
				double L1_cor_L2_dis = L1_cor_to_L1_cor_proj_L2.norm();
				double L1_cor_intersect_point_dis = L1_cor_L2_dis / sin_angle;

				// cut L1
				double s = L1_cor_intersect_point_dis;
				if (L1_cor_to_L1_cor_proj_L2.dot(L1.nor()) < 0) s = -s;

				Line sub_L1;
				if (L1.separate_s(sub_L1, s) == 1)	sub_L1s.push_back(sub_L1);

				// cut L2
				Eigen::Vector3d intersect_point = L1.cor() + s * L1.nor();
				double s2 = (intersect_point - L2.cor()).dot(L2.nor());

				Line sub_L2;
				if (L2.separate_s(sub_L2, s2) == 1)	sub_L2s.push_back(sub_L2);
			}
		}


		return sub_L1s.size() + sub_L2s.size();
	}

	int get_intersections(
		std::vector<Line>& L1s,
		std::vector<Line>& L2s
	) {
		SQI_VERBOSE_ONLY_COUT("");

		for (int i = 0; i < L1s.size(); ++i) { // note: the vector L1s's size is dymanic, so cannot use i_end = L1s.size.
			for (int j = 0, j_end = L2s.size(); j < j_end; ++j) { // note: the vector PC1s's size is dymanic, but the newly added curves do not need to be judged again.
				Line* L1 = &L1s[i];
				Line* L2 = &L2s[j];

				std::vector<Line> sub_L1s;
				std::vector<Line> sub_L2s;
				if (get_intersections(*L1, *L2, sub_L1s, sub_L2s) > 0) {

					for (auto l : sub_L1s) L1s.push_back(l);
					for (auto l : sub_L2s) L2s.push_back(l);
				}
			}
		}

		return L1s.size() + L2s.size();
	}

	// Line and CylindricLine (delete later)

	int get_intersections(
		Line& L1, ParameterizationCylindricLine& PL1, Cylinder& C1,
		std::vector<Line>& sub_L1s,
		std::vector<ParameterizationCylindricLine>& sub_PL1s
	) {
		// init
		std::vector<Line>().swap(sub_L1s);
		std::vector<ParameterizationCylindricLine>().swap(sub_PL1s);

		Line L2(C1.get_point(0, PL1.t()), C1.nor());

		Eigen::Vector3d cross_vec = L1.nor().cross(L2.nor());
		double sin_angle = cross_vec.norm();
		cross_vec /= sin_angle; // normalize
		double L1_L2_dis = (L2.cor() - L1.cor()).dot(cross_vec);

		if (L1_L2_dis > SQI_EPS) { // not intersect
			// do nothing
		}
		else {
			if (1 - std::abs(L2.nor().dot(L1.nor())) < SQI_EPS) { // two lines are parallel
				SQI_VERBOSE_ONLY_COUT("may error");
				// do nothing
			}
			else { // intersect a point
				Eigen::Vector3d L2_cor_to_L1_cor = L1.cor() - L2.cor();
				double L2_cor_to_L1_cor_dot_L2_nor = L2_cor_to_L1_cor.dot(L2.nor());
				Eigen::Vector3d L1_cor_proj_L2 = L2.cor() + L2_cor_to_L1_cor_dot_L2_nor * L2.nor();
				Eigen::Vector3d L1_cor_to_L1_cor_proj_L2 = L1_cor_proj_L2 - L1.cor();
				double L1_cor_L2_dis = L1_cor_to_L1_cor_proj_L2.norm();
				double L1_cor_intersect_point_dis = L1_cor_L2_dis / sin_angle;

				// cut L1
				double s = L1_cor_intersect_point_dis;
				if (L1_cor_to_L1_cor_proj_L2.dot(L1.nor()) < 0) s = -s;

				Line sub_L1;
				if (L1.separate_s(sub_L1, s) == 1)	sub_L1s.push_back(sub_L1);

				// cut PL1
				Eigen::Vector3d intersect_point = L1.cor() + s * L1.nor();
				double Cs = (intersect_point - C1.cor()).dot(C1.nor());

				ParameterizationCylindricLine sub_PL1;
				if (PL1.separate_s(sub_PL1, Cs) == 1) sub_PL1s.push_back(sub_PL1);
			}
		}

		return sub_L1s.size() + sub_PL1s.size();
	}

	int get_intersections(
		std::vector<Line>& L1s,
		std::vector<ParameterizationCylindricLine>& CL1s, Cylinder& C1
	) {
		SQI_VERBOSE_ONLY_COUT("");

		for (int i = 0; i < L1s.size(); ++i) { // note: the vector L1s's size is dymanic, so cannot use i_end = L1s.size.
			for (int j = 0, j_end = CL1s.size(); j < j_end; ++j) { // note: the vector PC1s's size is dymanic, but the newly added curves do not need to be judged again.
				Line* L1 = &L1s[i];
				ParameterizationCylindricLine* CL1 = &CL1s[j];

				std::vector<Line> sub_L1s;
				std::vector<ParameterizationCylindricLine> sub_CL1s;
				if (get_intersections(*L1, *CL1, C1, sub_L1s, sub_CL1s) > 0) {
					for (auto l : sub_L1s) L1s.push_back(l);
					for (auto l : sub_CL1s) CL1s.push_back(l);
				}
			}
		}

		return L1s.size() + CL1s.size();
	}

	// Line and Circle

	int get_intersections(
		Line& L1, ParameterizationCircle& c1,
		std::vector<Line>& sub_L1s,
		std::vector<ParameterizationCircle>& sub_c1s
	) {
		// init
		std::vector<Line>().swap(sub_L1s);
		std::vector<ParameterizationCircle>().swap(sub_c1s);

		Eigen::Vector3d L1_cor_to_c1_cor = c1.cor() - L1.cor();
		double L1_cor_to_c1_cor_dot_L1_nor = L1_cor_to_c1_cor.dot(L1.nor());
		Eigen::Vector3d c1_cor_proj_L1 = L1.cor() + L1_cor_to_c1_cor_dot_L1_nor * L1.nor();
		double c1_cor_L1_sq_dis = (c1.cor() - c1_cor_proj_L1).squaredNorm();
		double c1_sq_r = c1.r() * c1.r();

		if (c1_cor_L1_sq_dis > c1_sq_r + SQI_EPS) { // not intersect
			// SQI_VERBOSE_ONLY_COUT("not intersect");
			// do nothing
		}
		else if (c1_cor_L1_sq_dis > c1_sq_r - SQI_EPS) { // tangent: cut by 1 point
			// SQI_VERBOSE_ONLY_COUT("tangent");

			double s = L1_cor_to_c1_cor_dot_L1_nor;
			double ct = c1.get_t(c1_cor_proj_L1);
			if (L1.is_s_valid(s) && c1.is_t_valid(ct)) {
				// cut L1
				Line sub_L1;
				if (L1.separate_s(sub_L1, s) == 1)	sub_L1s.push_back(sub_L1);

				// cut c1
				ParameterizationCircle sub_c1;
				if (c1.separate_t(sub_c1, ct) == 1) sub_c1s.push_back(sub_c1);
			}
		}
		else { // intersect: cut by 2 points
			// SQI_VERBOSE_ONLY_COUT("intersect");
		
			double move_dis = std::sqrt(c1_sq_r - c1_cor_L1_sq_dis);
			double s1 = L1_cor_to_c1_cor_dot_L1_nor - move_dis;
			double s2 = L1_cor_to_c1_cor_dot_L1_nor + move_dis;
			double ct1 = c1.get_t(c1_cor_proj_L1 - move_dis * L1.nor());
			double ct2 = c1.get_t(c1_cor_proj_L1 + move_dis * L1.nor());
			if (L1.is_s_valid(s1) && c1.is_t_valid(ct1)) {
				// cut L1
				Line sub_L1;
				if (L1.separate_s(sub_L1, s1) == 1)	sub_L1s.push_back(sub_L1);

				// cut c1
				ParameterizationCircle sub_c1;
				if (c1.separate_t(sub_c1, ct1) == 1) sub_c1s.push_back(sub_c1);
			}
			if (L1.is_s_valid(s2)) {
				if (c1.is_t_valid(ct2)) {
					// cut L1
					Line sub_L1;
					if (L1.separate_s(sub_L1, s2) == 1)	sub_L1s.push_back(sub_L1);

					// cut c1
					ParameterizationCircle sub_c1;
					if (c1.separate_t(sub_c1, ct2) == 1) sub_c1s.push_back(sub_c1);
				}
				else if (sub_c1s.size() > 0 && sub_c1s[0].is_t_valid(ct2)) {
					// cut L1
					Line sub_L1;
					if (L1.separate_s(sub_L1, s2) == 1)	sub_L1s.push_back(sub_L1);

					// cut c1
					ParameterizationCircle sub_c1;
				    if (sub_c1s[0].separate_t(sub_c1, ct2) == 1) sub_c1s.push_back(sub_c1);
				}
			}
		}

		return sub_L1s.size() + sub_c1s.size();
	}

	int get_intersections(
		std::vector<Line>& L1s,
		std::vector<ParameterizationCircle>& c1s
	) {
		SQI_VERBOSE_ONLY_COUT("");

		for (int i = 0; i < L1s.size(); ++i) { // note: the vector L1s's size is dymanic, so cannot use i_end = L1s.size.
			for (int j = 0, j_end = c1s.size(); j < j_end; ++j) { // note: the vector PC1s's size is dymanic, but the newly added curves do not need to be judged again.
				Line* L1 = &L1s[i];
				ParameterizationCircle* c1 = &c1s[j];

				std::vector<Line> sub_L1s;
				std::vector<ParameterizationCircle> sub_c1s;
				if (get_intersections(*L1, *c1, sub_L1s, sub_c1s) > 0) {
					for (auto l : sub_L1s) L1s.push_back(l);
					for (auto c : sub_c1s) c1s.push_back(c);
				}
			}
		}

		return L1s.size() + c1s.size();
	}

	// Line and CylindricCurve

	int get_intersections(
		Line& L1, ParameterizationCylindricCurve& PC1, Cylinder& C1,
		std::vector<Line>& sub_L1s,
		std::vector<ParameterizationCylindricCurve>& sub_PC1s
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between a Line and a CylinderCurve");

		// init
		std::vector<Line>().swap(sub_L1s);
		std::vector<ParameterizationCylindricCurve>().swap(sub_PC1s);

		Eigen::Vector3d L1_cor_to_C1_cor = C1.cor() - L1.cor();

		if (1 - std::abs(L1.nor().dot(C1.nor())) < SQI_EPS) { // line parallel to cylinder's axis
			SQI_VERBOSE_ONLY_COUT("line -- parallel -- cylinder's axis");

			double L1_cor_to_C1_cor_dot_C1_nor = L1_cor_to_C1_cor.dot(C1.nor());
			double L1_C1_nor_sq_dis = L1_cor_to_C1_cor.squaredNorm() - L1_cor_to_C1_cor_dot_C1_nor * L1_cor_to_C1_cor_dot_C1_nor;
			double C1_sq_r = C1.r() * C1.r();

			if (std::abs(L1_C1_nor_sq_dis - C1_sq_r) < SQI_EPS) { // on the cylinder: a cylindric line cut the curve
				SQI_VERBOSE_ONLY_COUT("line -- on -- cylinder");

				double s, t;
				C1.get_s_t(L1.cor(), s, t);

				if (PC1.is_t_valid(t)) {
					std::vector<double> pc_ss;
					PC1.get_s(t, pc_ss);

					if (pc_ss.size() != 1) {
						SQI_VERBOSE_ONLY_WARNING("may error!");
						return 0;
					}

					Eigen::Vector3d intersect_point = C1.get_point(pc_ss[0], t);
					double ls = L1.get_s(intersect_point);

					if (L1.is_s_valid(ls)) { // cut!						
						Line sub_L1;
						if (L1.separate_s(sub_L1, ls) == 1) sub_L1s.push_back(sub_L1);

						ParameterizationCylindricCurve sub_PC1;
						if (PC1.separate_t(sub_PC1, t) == 1) sub_PC1s.push_back(sub_PC1);
					}
				}
			}
			else { // not on the cylinder
				SQI_VERBOSE_ONLY_COUT("line -- not on -- cylinder");
				// do nothing
			}
		}
		else { // line not parallel to cylinder's axis
			SQI_VERBOSE_ONLY_COUT("line -- not parallel -- cylinder's axis");

			std::vector<Point> intersect_points;
			get_intersections(L1, C1, intersect_points);

			if (intersect_points.size() == 1) { // tangent: 1 point
				SQI_VERBOSE_ONLY_COUT("tangent");

				double ls = L1.get_s(intersect_points[0].cor());
				if (L1.is_s_valid(ls) && get_intersections(intersect_points[0], PC1, C1, sub_PC1s) > 0) {
					// the curve is cutted, next to cut the line
					Line sub_L1;
					if (L1.separate_s(sub_L1, ls) == 1) sub_L1s.push_back(sub_L1);
				}
			}
			else if (intersect_points.size() == 2) { // intersect: 2 point
				SQI_VERBOSE_ONLY_COUT("intersect");

				double ls1 = L1.get_s(intersect_points[0].cor());
				double ls2 = L1.get_s(intersect_points[1].cor());
				if (ls1 > ls2) { // sort ls1 < ls2
					double tmp = ls1;
					ls1 = ls2; ls2 = tmp;
				}

				std::vector<ParameterizationCylindricCurve> sub_PC1s_1, sub_PC1s_2;
				if (L1.is_s_valid(ls1) && get_intersections(intersect_points[0], PC1, C1, sub_PC1s_1) > 0) {
					for (int i = 0, i_end = sub_PC1s_1.size(); i < i_end; ++i) sub_PC1s.push_back(sub_PC1s_1[i]);

					// the curve is cutted, next to cut the line
					Line sub_L1;
					if (L1.separate_s(sub_L1, ls1) == 1) sub_L1s.push_back(sub_L1);
				}

				if (L1.is_s_valid(ls2)) {
					if (get_intersections(intersect_points[1], PC1, C1, sub_PC1s_2) > 0) {
						for (int i = 0, i_end = sub_PC1s_2.size(); i < i_end; ++i) sub_PC1s.push_back(sub_PC1s_2[i]);

						Line sub_L2;
						if (L1.separate_s(sub_L2, ls2) == 1) sub_L1s.push_back(sub_L2);
					}
					else if (sub_PC1s.size() > 0 &&
						get_intersections(intersect_points[1], sub_PC1s[0], C1, sub_PC1s_2) > 0) {
						for (int i = 0, i_end = sub_PC1s_2.size(); i < i_end; ++i) sub_PC1s.push_back(sub_PC1s_2[i]);

						Line sub_L2;
						if (L1.separate_s(sub_L2, ls2) == 1) sub_L1s.push_back(sub_L2);
					}
				}
			}
			else { // not intersect
				SQI_VERBOSE_ONLY_COUT("not intersect");
			}
		}

		return sub_L1s.size() + sub_PC1s.size();
	}

	int get_intersections(
		std::vector<Line>& L1s,
		std::vector<ParameterizationCylindricCurve>& PC1s, Cylinder& C1
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between Lines and CylindricCurves");

		std::unordered_set<std::pair<int, int>, pair_hash> skip_set;
		for (int i = 0; i < L1s.size(); ++i) { // note: this vector's size is dymanic, so cannot use i_end = vector.size.
			for (int j = 0, j_end = PC1s.size(); j < j_end; ++j) { // note: this vector's size is dymanic, but the newly added curves do not need to be judged again.
				if (skip_set.find(std::make_pair(i, j)) != skip_set.end()) continue;
				
				Line* L1 = &L1s[i];
				ParameterizationCylindricCurve* PC1 = &PC1s[j];

				std::vector<Line> sub_L1s;
				std::vector<ParameterizationCylindricCurve> sub_PC1s;
				if (get_intersections(*L1, *PC1, C1, sub_L1s, sub_PC1s) > 0) {

					int cur_L1s_size = L1s.size();
					int cur_PC1s_size = PC1s.size();
					for (int ii = 0, ii_end = sub_L1s.size(); ii < ii_end; ++ii) {
						skip_set.insert(std::make_pair(cur_L1s_size + ii, j));
						for (int jj = 0, jj_end = sub_PC1s.size(); jj < jj_end; ++jj) {
							skip_set.insert(std::make_pair(cur_L1s_size + ii, cur_PC1s_size + jj));
						}
					}

					for (int ii = 0, ii_end = sub_L1s.size(); ii < ii_end; ++ii) L1s.push_back(sub_L1s[ii]);
					for (int ii = 0, ii_end = sub_PC1s.size(); ii < ii_end; ++ii) PC1s.push_back(sub_PC1s[ii]);
				}
			}
		}

		return L1s.size() + PC1s.size();
	}

	// Circle and Circle

	int get_intersections(
		ParameterizationCircle& c1, ParameterizationCircle& c2,
		std::vector<ParameterizationCircle>& sub_c1s,
		std::vector<ParameterizationCircle>& sub_c2s
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between Circles");

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
				double c2_t = c2.get_t(intersect_point1);
				if (c1.is_t_valid(c1_t) == true && c2.is_t_valid(c2_t) == true) {
					cutting_c1ts.push_back(c1_t);
					cutting_c2ts.push_back(c2_t);
				}
			}

			if (std::abs((intersect_point2 - c2.cor()).squaredNorm() - c2.r() * c2.r()) < SQI_EPS) { // intersect point is on c2
				double c1_t = c1.get_t(intersect_point2);
				double c2_t = c2.get_t(intersect_point2);

				if (c1.is_t_valid(c1_t) == true && c2.is_t_valid(c2_t) == true) {
					cutting_c1ts.push_back(c1_t);
					cutting_c2ts.push_back(c2_t);
				}
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
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between Circles");

		std::unordered_set<std::pair<int, int>, pair_hash> skip_set;
		for (int i = 0; i < c1s.size(); ++i) { // note: this vector's size is dymanic, so cannot use i_end = vector.size.
			for (int j = 0, j_end = c2s.size(); j < j_end; ++j) { // note: this vector's size is dymanic, but the newly added curves do not need to be judged again.
				if (skip_set.find(std::make_pair(i, j)) != skip_set.end()) continue;

				ParameterizationCircle* c1 = &c1s[i];
				ParameterizationCircle* c2 = &c2s[j];

				std::vector<ParameterizationCircle> sub_c1s;
				std::vector<ParameterizationCircle> sub_c2s;
				if (get_intersections(*c1, *c2, sub_c1s, sub_c2s) > 0) {

					int cur_c1s_size = c1s.size();
					int cur_c2s_size = c2s.size();
					for (int ii = 0, ii_end = sub_c1s.size(); ii < ii_end; ++ii) {
						skip_set.insert(std::make_pair(cur_c1s_size + ii, j));
						for (int jj = 0, jj_end = sub_c2s.size(); jj < jj_end; ++jj) {
							skip_set.insert(std::make_pair(cur_c1s_size + ii, cur_c2s_size + jj));
						}
					}

					for (int ii = 0, ii_end = sub_c1s.size(); ii < ii_end; ++ii) c1s.push_back(sub_c1s[ii]);
					for (int ii = 0, ii_end = sub_c2s.size(); ii < ii_end; ++ii) c2s.push_back(sub_c2s[ii]);
				}
			}
		}

		return c1s.size() + c2s.size();
	}

	// Circle and CylindricCurve

	int get_intersections(
		ParameterizationCircle& c1, ParameterizationCylindricCurve& PC1, Cylinder& C1,
		std::vector<ParameterizationCircle>& sub_c1s,
		std::vector<ParameterizationCylindricCurve>& sub_PC1s
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between a Circle and a CylindricCurve");

		// init 
		std::vector<ParameterizationCircle>().swap(sub_c1s);
		std::vector<ParameterizationCylindricCurve>().swap(sub_PC1s);

		Plane P1(c1.cor(), c1.nor()); // proxy plane that is co-faceted with the circle

		std::vector<Line> lines;
		std::vector<ParameterizationCylindricCurve> curves;
		get_intersections(P1, C1, lines, curves);

		if (lines.size() > 0) { // circle is parallel to cylinder's axis
			SQI_VERBOSE_ONLY_COUT("parallel to cylinder's axis");

			std::vector<ParameterizationCylindricCurve> tmp_PC1s = { PC1 };
			get_intersections(lines, tmp_PC1s, C1);
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

					if (std::abs((intersect_point - c1.cor()).norm() - c1.r()) > SQI_EPS) {
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
			SQI_VERBOSE_ONLY_COUT("not intersect");
			return 0;
		}

		return sub_c1s.size() + sub_PC1s.size();
	}

	int get_intersections(
		std::vector<ParameterizationCircle>& c1s,
		std::vector<ParameterizationCylindricCurve>& PC1s, Cylinder& C1
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between Circles and CylindricCurves");

		std::unordered_set<std::pair<int, int>, pair_hash> skip_set;
		for (int i = 0; i < c1s.size(); ++i) { // note: this vector's size is dymanic, so cannot use i_end = vector.size.
			for (int j = 0, j_end = PC1s.size(); j < j_end; ++j) { // note: this vector's size is dymanic, but the newly added curves do not need to be judged again.
				if (skip_set.find(std::make_pair(i, j)) != skip_set.end()) continue;

				ParameterizationCircle* c1 = &c1s[i];
				ParameterizationCylindricCurve* PC1 = &PC1s[j];

				std::vector<ParameterizationCircle> sub_c1s;
				std::vector<ParameterizationCylindricCurve> sub_PC1s;
				if (get_intersections(*c1, *PC1, C1, sub_c1s, sub_PC1s) > 0) {

					int cur_c1s_size = c1s.size();
					int cur_PC1s_size = PC1s.size();
					for (int ii = 0, ii_end = sub_c1s.size(); ii < ii_end; ++ii) {
						skip_set.insert(std::make_pair(cur_c1s_size + ii, j));
						for (int jj = 0, jj_end = sub_PC1s.size(); jj < jj_end; ++jj) {
							skip_set.insert(std::make_pair(cur_c1s_size + ii, cur_PC1s_size + jj));
						}
					}

					for (int ii = 0, ii_end = sub_c1s.size(); ii < ii_end; ++ii) c1s.push_back(sub_c1s[ii]);
					for (int ii = 0, ii_end = sub_PC1s.size(); ii < ii_end; ++ii) PC1s.push_back(sub_PC1s[ii]);
				}
			}
		}

		return c1s.size() + PC1s.size();
	}

	// CylindricLine and CylindricCurve (delete later)

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

	// CylindricCurve and CylindricCurve (on the same cylinder)

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

			// check start endpoint
			/* {
				std::vector<double> ss1, ss2;
				if (PC1.get_s(overlap_t_lb, ss1) == 1 && PC2.get_s(overlap_t_lb, ss2) == 1) {
					if (std::abs(ss1[0] - ss2[0]) < SQI_EPS) { // cut
						ParameterizationCylindricCurve sub_PC1, sub_PC2;
						if (PC1.separate_t(sub_PC1, overlap_t_lb) == 1) sub_PC1s.push_back(sub_PC1);
						if (PC2.separate_t(sub_PC2, overlap_t_lb) == 1) sub_PC2s.push_back(sub_PC2);
					}
				}
			}*/

			int pre_status = 0, status = 0;
			double pre_t;
			double t = overlap_t_lb + SQI_EPS;

			while (t < overlap_t_ub + t_step) {
				double cur_t = (t > overlap_t_ub) ? (overlap_t_ub) : t;

				std::vector<double> C1_ss, C2_ss;
				PC1.get_s(cur_t, C1_ss);
				PC2.get_s(cur_t, C2_ss);
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
							double t_ub = cur_t;

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

			// check end endpoint
			/* {
				std::vector<double> ss1, ss2;
				if (PC1.get_s(overlap_t_ub, ss1) == 1 && PC2.get_s(overlap_t_ub, ss2) == 1) {
					if (std::abs(ss1[0] - ss2[0]) < SQI_EPS) { // cut
						ParameterizationCylindricCurve sub_PC1, sub_PC2;
						if (PC1.separate_t(sub_PC1, overlap_t_ub) == 1) sub_PC1s.push_back(sub_PC1);
						if (PC2.separate_t(sub_PC2, overlap_t_ub) == 1) sub_PC2s.push_back(sub_PC2);
					}
				}
			}*/
		}

		return sub_PC1s.size() + sub_PC2s.size();
	}

	int get_intersections(
		std::vector<ParameterizationCylindricCurve>& PC1s,
		std::vector<ParameterizationCylindricCurve>& PC2s,
		double t_step
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between CylindricCurves");

		std::unordered_set<std::pair<int, int>, pair_hash> skip_set;
		for (int i = 0; i < PC1s.size(); ++i) { // note: this vector's size is dymanic, so cannot use i_end = vector.size.
			for (int j = 0, j_end = PC2s.size(); j < j_end; ++j) { // note: this vector's size is dymanic, but the newly added curves do not need to be judged again.
				if (skip_set.find(std::make_pair(i, j)) != skip_set.end()) continue;

				ParameterizationCylindricCurve* PC1 = &PC1s[i];
				ParameterizationCylindricCurve* PC2 = &PC2s[j];

				std::vector<ParameterizationCylindricCurve> sub_PC1s;
				std::vector<ParameterizationCylindricCurve> sub_PC2s;
				if (get_intersections(*PC1, *PC2, sub_PC1s, sub_PC2s) > 0) {

					int cur_PC1s_size = PC1s.size();
					int cur_PC2s_size = PC2s.size();
					for (int ii = 0, ii_end = sub_PC1s.size(); ii < ii_end; ++ii) {
						skip_set.insert(std::make_pair(cur_PC1s_size + ii, j));
						for (int jj = 0, jj_end = sub_PC2s.size(); jj < jj_end; ++jj) {
							skip_set.insert(std::make_pair(cur_PC1s_size + ii, cur_PC2s_size + jj));
						}
					}

					for (int ii = 0, ii_end = sub_PC1s.size(); ii < ii_end; ++ii) PC1s.push_back(sub_PC1s[ii]);
					for (int ii = 0, ii_end = sub_PC2s.size(); ii < ii_end; ++ii) PC2s.push_back(sub_PC2s[ii]);
				}
			}
		}

		return PC1s.size() + PC2s.size();
	}

	// CylindricCurve and CylindricCurve (on different cylindrs)

	int get_intersections(
		std::vector<ParameterizationCylindricCurve>& PC1s, Cylinder& C1,
		std::vector<ParameterizationCylindricCurve>& PC2s, Cylinder& C2,
		double t_step
	) {
		SQI_VERBOSE_ONLY_TITLE("compute the intersections between CylindricCurves (on different cylinders)");

		// get the intersections between C1 and C2
		std::vector<Point> points;
		std::vector<Line> lines;
		std::vector<ParameterizationCylindricCurve> curves;
		get_intersections(C1, C2, points, lines, curves);

		if (points.size() > 0) { // C1 and C2 are tangent
			if (points.size() != 1) {
				SQI_VERBOSE_ONLY_WARNING("may error!");
				return 0;
			}

			get_intersections(points, PC1s, C1);
			get_intersections(points, PC2s, C2);
		}
		else if (lines.size() > 0) { // C1's axis and C2's axis is parallel
			std::vector<Line> tmp_lines1 = lines;
			get_intersections(tmp_lines1, PC1s, C1);
			
			std::vector<Line> tmp_lines2 = lines;
			get_intersections(tmp_lines2, PC2s, C2);
		}
		else if (curves.size() > 0) {
			// note: curves are on the cylinder C2!
			int original_PC2s_size = PC2s.size();
			get_intersections(curves, PC2s, t_step);

			// cut the PC1
			std::vector<Point> cutting_points;
			for (int i = original_PC2s_size, i_end = PC2s.size(); i < i_end; ++i) {
				ParameterizationCylindricCurve* PC = &PC2s[i];
				double t2 = PC->t_ub();

				std::vector<double> ss2;
				PC->get_s(t2, ss2);
				if (ss2.size() != 1) {
					SQI_VERBOSE_ONLY_WARNING("may error");
					return 0;
				}

				Eigen::Vector3d intersect_point = C2.get_point(ss2[0], t2);
				cutting_points.push_back(Point(intersect_point));
			}

			get_intersections(cutting_points, PC1s, C1);
		}
		else { // not intersect
			// do nothing
		}

		return PC1s.size() + PC2s.size();
	}

	// --------------------------------------------------------------------------------------------













}