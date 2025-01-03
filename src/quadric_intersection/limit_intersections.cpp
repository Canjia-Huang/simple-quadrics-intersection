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

	bool Line::limited_by(Plane& P) {
		if (P.constraint_radius > SQI_LOOSE_EPS) { // limit by the ball center at cor_ with radius equal to constraint_radius
			Eigen::Vector3d cor_to_P_cor = P.cor() - cor_;
			double cor_to_P_cor_dot_nor = cor_to_P_cor.dot(nor_); // proj's s
			Eigen::Vector3d P_cor_proj = cor_ + cor_to_P_cor_dot_nor * nor_;
			double P_cor_to_proj_sq_dis = (P.cor() - P_cor_proj).squaredNorm();
			double extend_s = std::sqrt(P.constraint_radius * P.constraint_radius - P_cor_to_proj_sq_dis);

			double limit_s_lb = cor_to_P_cor_dot_nor - extend_s;
			double limit_s_ub = cor_to_P_cor_dot_nor + extend_s;

			s_lb_ = std::max(s_lb_, limit_s_lb);
			s_ub_ = std::min(s_ub_, limit_s_ub);
			if (s_ub_ - s_lb_ < SQI_EPS) return false;
			return true;
		}

		if (P.vertices.size() != 3) return true;

		// limit by the triangle defined by vertices
		Eigen::Vector3d* P1 = &(P.vertices[0]);
		Eigen::Vector3d* P2 = &(P.vertices[1]);
		Eigen::Vector3d* P3 = &(P.vertices[2]);
		
		double P1_cor_dot_nor = (*P1 - cor_).dot(nor_);
		double P2_cor_dot_nor = (*P2 - cor_).dot(nor_);
		double P3_cor_dot_nor = (*P3 - cor_).dot(nor_);

		Eigen::Vector3d P1_proj_L = cor_ + P1_cor_dot_nor * nor_;
		Eigen::Vector3d P2_proj_L = cor_ + P2_cor_dot_nor * nor_;
		Eigen::Vector3d P3_proj_L = cor_ + P3_cor_dot_nor * nor_;
		double P1_proj_L_norm = (*P1 - P1_proj_L).norm();
		double P2_proj_L_norm = (*P2 - P2_proj_L).norm();
		double P3_proj_L_norm = (*P3 - P3_proj_L).norm();

		Eigen::Vector3d Edge12_intersect_point = P2_proj_L_norm * (*P1) + P1_proj_L_norm * (*P2);
		Edge12_intersect_point /= (P1_proj_L_norm + P2_proj_L_norm);
		Eigen::Vector3d Edge23_intersect_point = P3_proj_L_norm * (*P2) + P2_proj_L_norm * (*P3);
		Edge23_intersect_point /= (P2_proj_L_norm + P3_proj_L_norm);
		Eigen::Vector3d Edge31_intersect_point = P1_proj_L_norm * (*P3) + P3_proj_L_norm * (*P1);
		Edge31_intersect_point /= (P3_proj_L_norm + P1_proj_L_norm);

		std::vector<double> intersect_s;
		Eigen::Vector3d cor_to_Edge12_intersect_point = Edge12_intersect_point - cor_;
		Eigen::Vector3d cor_to_Edge23_intersect_point = Edge23_intersect_point - cor_;
		Eigen::Vector3d cor_to_Edge31_intersect_point = Edge31_intersect_point - cor_;
		double s1 = cor_to_Edge12_intersect_point.dot(nor_);
		double s2 = cor_to_Edge23_intersect_point.dot(nor_);
		double s3 = cor_to_Edge31_intersect_point.dot(nor_);
		if (1 - std::abs(s1 / cor_to_Edge12_intersect_point.norm()) < SQI_EPS) intersect_s.push_back(s1);
		if (1 - std::abs(s2 / cor_to_Edge23_intersect_point.norm()) < SQI_EPS) intersect_s.push_back(s2);
		if (1 - std::abs(s3 / cor_to_Edge31_intersect_point.norm()) < SQI_EPS) intersect_s.push_back(s3);

		if (intersect_s.size() == 2) { // intersect
			std::sort(intersect_s.begin(), intersect_s.end());
			s_lb_ = std::max(s_lb_, intersect_s[0]);
			s_ub_ = std::min(s_ub_, intersect_s[1]);

			if (s_ub_ - s_lb_ < SQI_EPS) return false;
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

	bool ParameterizationCircle::limited_by(Sphere& S) {
		double cos_s_ub = std::cos(ang2rad(S.s_ub()));
		Eigen::Vector3d Sc_cor = S.cor() + S.r() * cos_s_ub * S.nor();
		double Sc_r = S.r() * std::sqrt(1 - cos_s_ub * cos_s_ub);

		ParameterizationCircle Sc(Sc_cor, S.nor(), Sc_r);

		std::vector<ParameterizationCircle> sub_cs, sub_Scs;
		get_intersections(*this, Sc, sub_cs, sub_Scs);
		sub_cs.push_back(*this);

		for (int i = 0, i_end = sub_cs.size(); i < i_end; ++i) {
			ParameterizationCircle* cur_s = &(sub_cs[i]);
			Eigen::Vector3d c_p = cur_s->get_point(0.5 * (cur_s->t_lb() + cur_s->t_ub()));

			double Ss, St;
			S.get_s_t(c_p, Ss, St);

			if (Ss > S.s_lb() && Ss < S.s_ub()) {
				*this = *cur_s;
				return true;
			}
		}

		return false;
	}

	bool ParameterizationCylindricCurve::limited_by(
		Plane& P,
		std::vector<ParameterizationCylindricCurve>& sub_PCs
	) {
		if (P.vertices.size() != 3) return true;

		Eigen::Vector3d* P1 = &(P.vertices[0]);
		Eigen::Vector3d* P2 = &(P.vertices[1]);
		Eigen::Vector3d* P3 = &(P.vertices[2]);

		Line L12(*P1, (*P2 - *P1));
		Line L23(*P2, (*P3 - *P2));
		Line L31(*P3, (*P1 - *P3));

		std::vector<Point> L12_intersect_points, L23_intersect_points, L31_intersect_points;
		get_intersections(L12, C_, L12_intersect_points);
		get_intersections(L23, C_, L23_intersect_points);
		get_intersections(L31, C_, L31_intersect_points);

		for (int i = 0, i_end = L12_intersect_points.size(); i < i_end; ++i) {
			Point* P = &(L12_intersect_points[i]);

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
				if (Cs < C.s_lb() || Cs > C.s_ub() || 
					Ct < C.t_lb() || Ct > C.t_ub()) continue;

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