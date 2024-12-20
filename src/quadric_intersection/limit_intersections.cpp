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

	bool ParameterizationCylindricCurve::limited_by(Cylinder& C) {
		/*
		Line L_t_lb(C.get_point(0, C.t_lb()), C.nor());
		Line L_t_ub(C.get_point(0, C.t_ub()), C.nor());

		std::vector<Line> sub_L_t_lb_s;
		std::vector<ParameterizationCylindricCurve> sub_PC;
		get_intersections(L_t_lb, *this, sub_L_t_lb_s, sub_PC);*/

		return true;
	}
}