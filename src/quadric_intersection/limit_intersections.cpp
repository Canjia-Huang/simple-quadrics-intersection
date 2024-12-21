#include "simple_quadrics_intersection.h"

namespace QuadricsIntersection
{
	bool Point::limited_by(Cylinder& C) {
		Eigen::Vector3d P_cor_to_C_cor = C.cor() - cor_;
		double P_cor_to_C_cor_dot_C_nor = P_cor_to_C_cor.dot(C.nor());
		double P_C_sq_dis = P_cor_to_C_cor.squaredNorm() - P_cor_to_C_cor_dot_C_nor * P_cor_to_C_cor_dot_C_nor;

		if (std::abs(P_C_sq_dis - C.r() * C.r()) > SQI_EPS) { // this Point is not on the Cylinder
			return false;
		}
		return true;
	}

	bool Line::limited_by(Cylinder& C) {
		Eigen::Vector3d L_cor_to_C_cor = C.cor() - cor_;
		double L_cor_to_C_cor_dot_C_nor = L_cor_to_C_cor.dot(C.nor());
		double L_C_sq_dis = L_cor_to_C_cor.squaredNorm() - L_cor_to_C_cor_dot_C_nor * L_cor_to_C_cor_dot_C_nor;

		if (std::abs(L_C_sq_dis - C.r() * C.r()) < SQI_EPS) { // this Line is on Cylinder C
			if (1 - std::abs(nor_.dot(C.nor())) < SQI_EPS) { // this Line is parallel to Cylinder C
				double Cs, Ct;
				C.get_s_t(cor_, Cs, Ct);
				if (Ct < C.t_lb() - SQI_EPS || Ct > C.t_ub() + SQI_EPS) return false;

				// cut this line
				s_lb_ = std::max(s_lb_, C.s_lb());
				s_ub_ = std::min(s_ub_, C.s_ub());

				if (s_ub_ - s_lb_ < SQI_EPS) return false;
			}
			else {
				return false;
			}
		}
		else { // this Line is not on Cylinder C
			return false;
		}

		return true;
	}

	int ParameterizationCylindricCurve::get_t(double s, double& t) {
		double t_lb = t_lb_;
		double t_ub = t_ub_;
		if (t_ub - t_lb > 360 - SQI_EPS) {
			SQI_VERBOSE_ONLY_WARNING("may error!");
			return 0;
		}

		std::vector<double> ss1, ss2;
		if (get_s(t_lb, ss1) != 1 || get_s(t_ub, ss2) != 1) {
			SQI_VERBOSE_ONLY_WARNING("may error!");
			return 0;
		}
		int status_t_lb = (ss1[0] > s) ? 1 : -1;
		int status_t_ub = (ss2[0] > s) ? 1 : -1;

		if (status_t_lb == status_t_ub) { // can not find t
			return 0;
		}

		while (t_ub - t_lb > SQI_EPS) {
			double t_mid = 0.5 * (t_lb + t_ub);
			std::vector<double> ss;
			if (get_s(t_mid, ss) != 1) {
				SQI_VERBOSE_ONLY_WARNING("may error!");
				return 0;
			}
			int status_t_mid = (ss[0] > s) ? 1 : -1;

			if (status_t_mid == status_t_lb) t_lb = t_mid;
			else t_ub = t_mid;
		}

		t = t_lb;

		return 1;
	}

	int ParameterizationCylindricCurve::get_t(Cylinder& C, double s, double& t) {
		double t_lb = t_lb_;
		double t_ub = t_ub_;
		if (t_ub - t_lb > 360 - SQI_EPS) {
			SQI_VERBOSE_ONLY_WARNING("may error!");
			return 0;
		}

		std::vector<double> ss1, ss2;
		if (get_s(t_lb, ss1) != 1 || get_s(t_ub, ss2) != 1) {
			SQI_VERBOSE_ONLY_WARNING("may error!");
			return 0;
		}
		Eigen::Vector3d P1 = C_.get_point(ss1[0], t_lb);
		Eigen::Vector3d P2 = C_.get_point(ss2[0], t_ub);
		double Cs1, Ct1, Cs2, Ct2;
		C.get_s_t(P1, Cs1, Ct1);
		C.get_s_t(P2, Cs2, Ct2);
		int status_t_lb = (Cs1 > s) ? 1 : -1;
		int status_t_ub = (Cs2 > s) ? 1 : -1;

		if (status_t_lb == status_t_ub) { // can not find t
			return 0;
		}

		while (t_ub - t_lb > SQI_EPS) {
			double t_mid = 0.5 * (t_lb + t_ub);
			std::vector<double> ss;
			if (get_s(t_mid, ss) != 1) {
				SQI_VERBOSE_ONLY_WARNING("may error!");
				return 0;
			}

			Eigen::Vector3d P = C_.get_point(ss[0], t_mid);
			double Cs, Ct;
			C.get_s_t(P, Cs, Ct);

			int status_t_mid = (Cs > s) ? 1 : -1;

			if (status_t_mid == status_t_lb) t_lb = t_mid;
			else t_ub = t_mid;
		}

		t = t_lb;

		return 1;
	}

	bool ParameterizationCylindricCurve::limited_by(
	) {
		// limit t
		t_lb_ = std::max(t_lb_, C_.t_lb());
		t_ub_ = std::min(t_ub_, C_.t_ub());
		if (t_ub_ - t_lb_ < SQI_EPS) return false;

		// limit s
		double s_lb_t, s_ub_t;
		int find_status_lb = get_t(C_.s_lb(), s_lb_t);
		int find_status_ub = get_t(C_.s_ub(), s_ub_t);
		// SQI_VERBOSE_ONLY_TEST("status lb:" << find_status_lb << " " << "status ub:" << find_status_ub);
		if (find_status_lb == 1 && find_status_ub == 1) { // get [s_lb_t, s_ub_t]
			// sort to s_lb_t < s_ub_t
			if (s_lb_t > s_ub_t) {
				double tmp_t = s_ub_t;
				s_ub_t = s_lb_t; s_lb_t = tmp_t;
			}

			ParameterizationCylindricCurve sub_PC1, sub_PC2;
			separate_t(sub_PC1, s_lb_t);
			if (separate_t(sub_PC2, s_ub_t) == 1) *this = sub_PC2;
			else return false;
		}
		else {
			std::vector<double> ss;
			get_s(t_lb_ + SQI_EPS, ss);
			if (ss.size() != 1) {
				SQI_VERBOSE_ONLY_WARNING("may error!");
				return false;
			}

			if (find_status_lb == 1) {
				if (ss[0] < C_.s_lb()) { // get [s_lb_t, t_ub]
					ParameterizationCylindricCurve sub_PC1;
					separate_t(sub_PC1, s_lb_t);
				}
				else { // get [t_lb, s_lb_t]
					ParameterizationCylindricCurve sub_PC1;
					if (separate_t(sub_PC1, s_lb_t) == 1) *this = sub_PC1;
					else return false;
				}
			}
			
			if (find_status_ub == 1) {
				if (ss[0] < C_.s_ub()) { // get [t_lb, s_ub_t]
					ParameterizationCylindricCurve sub_PC1;
					if (separate_t(sub_PC1, s_ub_t) == 1) *this = sub_PC1;
					else return false;
				}
				else { // get [s_ub_t, t_ub]
					ParameterizationCylindricCurve sub_PC1;
					separate_t(sub_PC1, s_ub_t);
				}
			}
		}

		return true;
	}

	bool ParameterizationCylindricCurve::limited_by(
		Cylinder& C, 
		std::vector<ParameterizationCylindricCurve>& sub_PCs
	) {
		// init
		std::vector<ParameterizationCylindricCurve>().swap(sub_PCs);
		std::vector<double> cut_ts;

		// cut by C's t limit
		if (C.t_ub() - C.t_lb() > 360 - SQI_EPS) { // do not need to cut by t
			// do nothing
		}
		else {
			Line C_t_lb_line(C.get_point(0, C.t_lb()), C.nor());
			Line C_t_ub_line(C.get_point(0, C.t_ub()), C.nor());
			std::vector<Line> C_t_lines = { C_t_lb_line, C_t_ub_line };
			std::vector<ParameterizationCylindricCurve> tmp_curves = { *this };
			get_intersections(C_t_lines, tmp_curves);

			for (int i = 0, i_end = tmp_curves.size(); i < i_end; ++i) {
				double C_t = tmp_curves[i].t_ub();
				std::vector<double> C_ss;
				get_s(C_t, C_ss);
				if (C_ss.size() != 1) continue;

				Eigen::Vector3d intersected_point = C_.get_point(C_ss[0], C_t);

				double Cs, Ct;
				C.get_s_t(intersected_point, Cs, Ct);
				if (Cs < C.s_lb() || Cs > C.s_ub()) continue;

				cut_ts.push_back(C_t);
			}
		}

		// cut by C's s limit
		if (C.s_lb() < -0.5 * SQI_INFTY || C.s_ub() > 0.5 * SQI_INFTY) { // do not need to cut by s
			// do nothing
		}
		else {
			/* // this way is still not done yet...
			ParameterizationCircle C_s_lb_circle(C.cor() + C.s_lb() * C.nor(), C.nor(), C.r());
			ParameterizationCircle C_s_ub_circle(C.cor() + C.s_ub() * C.nor(), C.nor(), C.r());
			std::vector<ParameterizationCircle> C_s_cirlces = { C_s_lb_circle, C_s_ub_circle };
			std::vector<ParameterizationCylindricCurve> tmp_curves = { *this };
			get_intersections(C_s_cirlces, tmp_curves);
			write_result_points("..//data//intersects.obj", C_s_cirlces, tmp_curves);

			for (int i = 0, i_end = tmp_curves.size(); i < i_end; ++i) {
				double C_t = tmp_curves[i].t_ub();
				std::vector<double> C_ss;
				get_s(C_t, C_ss);
				if (C_ss.size() != 1) continue;

				Eigen::Vector3d intersected_point = C_.get_point(C_ss[0], C_t);

				double Cs, Ct;
				C.get_s_t(intersected_point, Cs, Ct);
				if (Ct < C.t_lb() || Ct > C.t_ub()) continue;

				cut_ts.push_back(C_t);
			}*/

			// another way
			double s_lb_t, s_ub_t;
			if (get_t(C, C.s_lb(), s_lb_t) == 1) {
				std::vector<double> C_ss;
				get_s(s_lb_t, C_ss);
				if (C_ss.size() == 1) {
					Eigen::Vector3d intersected_point = get_point(C_ss[0], s_lb_t);
					double Cs, Ct;
					C.get_s_t(intersected_point, Cs, Ct);
					if (Ct > C.t_lb() && Ct < C.t_ub()) {
						cut_ts.push_back(s_lb_t);
					}
				}
			}
			if (get_t(C, C.s_ub(), s_ub_t) == 1) {
				std::vector<double> C_ss;
				get_s(s_ub_t, C_ss);
				if (C_ss.size() == 1) {
					Eigen::Vector3d intersected_point = get_point(C_ss[0], s_ub_t);
					double Cs, Ct;
					C.get_s_t(intersected_point, Cs, Ct);
					if (Ct > C.t_lb() && Ct < C.t_ub()) {
						cut_ts.push_back(s_ub_t);
					}
				}
			}
		}

		// cut!
		std::sort(cut_ts.begin(), cut_ts.end());
		for (int i = 0, i_end = cut_ts.size(); i < i_end; ++i) {
			ParameterizationCylindricCurve sub_PC;
			if (this->separate_t(sub_PC, cut_ts[i]) == 1) {
				double mid_t = 0.5 * (sub_PC.t_lb() + sub_PC.t_ub());
				std::vector<double> mid_ss;
				sub_PC.get_s(mid_t, mid_ss);
				if (mid_ss.size() != 1) continue;

				Eigen::Vector3d mid_point = sub_PC.get_point(mid_ss[0], mid_t);

				double Cs, Ct;
				C.get_s_t(mid_point, Cs, Ct);
				if (Cs < C.s_lb() || Cs > C.s_ub() || Ct < C.t_lb() || Ct > C.t_ub()) continue;

				sub_PCs.push_back(sub_PC);
			}
		}
		{
			double mid_t = 0.5 * (t_lb_ + t_ub_);
			std::vector<double> mid_ss;
			get_s(mid_t, mid_ss);
			if (mid_ss.size() != 1) return false;

			Eigen::Vector3d mid_point = get_point(mid_ss[0], mid_t);

			double Cs, Ct;
			C.get_s_t(mid_point, Cs, Ct);
			if (Cs < C.s_lb() || Cs > C.s_ub() || Ct < C.t_lb() || Ct > C.t_ub()) return false;
		}

		return true;
	}
}