#include "simple_quadrics_intersection.h"

namespace QuadricsIntersection 
{
	int get_intersections(
		Line& L1, Cylinder& C1,
		std::vector<ParameterizationCylindricPoint>& points
	) {
		// init
		std::vector<ParameterizationCylindricPoint>().swap(points);
		points.reserve(2);

		if (std::abs(std::abs(L1.nor().dot(C1.nor())) - 1) < SQI_EPS) { // line vector is parallel to the cylinder's axis
			SQI_VERBOSE_ONLY_COUT("axis parallel, may error");
		}
		else {
			Eigen::Vector3d L1_cor_C1_cor = L1.cor() - C1.cor();
			double L1_cor_C1_cor_dot_L1_nor = L1_cor_C1_cor.dot(L1.nor());
			double L1_cor_C1_cor_dot_C1_nor = L1_cor_C1_cor.dot(C1.nor());
			double L1_nor_dot_C1_nor = L1.nor().dot(C1.nor());

			// quadratic function
			double a = 1 - L1_nor_dot_C1_nor * L1_nor_dot_C1_nor;
			double b = 2 * L1_cor_C1_cor_dot_L1_nor - 2 * L1_cor_C1_cor_dot_C1_nor * L1_nor_dot_C1_nor;
			double c = L1_cor_C1_cor.squaredNorm() - L1_cor_C1_cor_dot_C1_nor * L1_cor_C1_cor_dot_C1_nor - C1.r() * C1.r();

			double delta = b * b - 4 * a * c;
			if (delta < -SQI_EPS) {
				SQI_VERBOSE_ONLY_COUT("not intersect");
			}
			else if (delta < SQI_EPS) {
				SQI_VERBOSE_ONLY_COUT("tangent");

				double w = -b / (2 * a);

				double s, t;
				C1.get_s_t(L1.cor() + w * L1.nor(), s, t);
				points.push_back(ParameterizationCylindricPoint(s, t));
			}
			else {
				SQI_VERBOSE_ONLY_COUT("two-point intersection");

				double sqrt_delta = std::sqrt(delta);
				double w1 = (-b + sqrt_delta) / (2 * a);
				double w2 = (-b - sqrt_delta) / (2 * a);

				double s1, t1, s2, t2;
				C1.get_s_t(L1.cor() + w1 * L1.nor(), s1, t1);
				C1.get_s_t(L1.cor() + w2 * L1.nor(), s2, t2);
				points.push_back(ParameterizationCylindricPoint(s1, t1));
				points.push_back(ParameterizationCylindricPoint(s2, t2));
			}
		}

		return points.size();
	}

	int get_intersections(
		Plane& P1, Cylinder& C1,
		std::vector<ParameterizationCylindricLine>& lines,
		std::vector<ParameterizationCylindricCurve>& curves
	) {
		SQI_VERBOSE_ONLY_COUT("");
		// SQI_VERBOSE_ONLY_COUT("plane" << " " << "cor:" << P1.cor().transpose() << " " << "nor:" << P1.nor().transpose());
		// SQI_VERBOSE_ONLY_COUT("cylinder" << " " << "cor:" << C1.cor().transpose() << " " << "nor:" << C1.nor().transpose() << " " << "r:" << C1.r());

		// init 
		std::vector<ParameterizationCylindricLine>().swap(lines);
		std::vector<ParameterizationCylindricCurve>().swap(curves);
		Eigen::Vector3d P1_cor_to_C1_cor = C1.cor() - P1.cor();

		if (P1.nor().dot(C1.nor()) < SQI_EPS) { // plane is parallel to cylinder	
			double P1_cor_to_C1_cor_dot_P1_nor = P1_cor_to_C1_cor.dot(P1.nor());
			double P1_cor_to_C1_dis = std::abs(P1_cor_to_C1_cor_dot_P1_nor);

			if (P1_cor_to_C1_dis > C1.r() + SQI_EPS) { // not intersect
				SQI_VERBOSE_ONLY_COUT("not interest");
				return 0;
			}
			else if (P1_cor_to_C1_dis < C1.r() - SQI_EPS) { // intersect: result 2 lines
				SQI_VERBOSE_ONLY_COUT("intersect");

				double rot_angle = safetyAcos(P1_cor_to_C1_dis / C1.r());

				Eigen::Vector3d mid_point = C1.cor() - P1_cor_to_C1_cor_dot_P1_nor * P1.nor();
				double s, t;
				C1.get_s_t(mid_point, s, t);

				lines.push_back(ParameterizationCylindricLine(t + rot_angle));
				lines.push_back(ParameterizationCylindricLine(t - rot_angle));
			}
			else { // tangent: result 1 line
				SQI_VERBOSE_ONLY_COUT("tangent");

				Eigen::Vector3d intersect_point = C1.cor() - P1_cor_to_C1_cor_dot_P1_nor * P1.nor();
				double s, t;
				C1.get_s_t(intersect_point, s, t);
				lines.push_back(ParameterizationCylindricLine(t));
			}
		}
		else { // not parallel
			std::vector<double> a_t = {
					0
			};
			SQI_VERBOSE_ONLY_COUT("a_t:" << " " << a_t[0]);

			std::vector<double> b_t = {
					0,
					0,
					C1.nor().dot(P1.nor())
			};
			SQI_VERBOSE_ONLY_COUT("b_t:" << " " << b_t[0] << " " << b_t[1] << " " << b_t[2]);

			std::vector<double> c_t = {
					0,
					0,
					0,
					C1.r() * C1.u().dot(P1.nor()),
					C1.r() * C1.v().dot(P1.nor()),
					P1_cor_to_C1_cor.dot(P1.nor())
			};
			SQI_VERBOSE_ONLY_COUT("c_t:" << " " << c_t[0] << " " << c_t[1] << " " << c_t[2] << " " << c_t[3] << " " << c_t[4] << " " << c_t[5]);

			ParameterizationCylindricCurve PC(
				a_t, b_t, c_t,
				-SQI_INFTY, SQI_INFTY, 0, 360
			);

			curves.push_back(PC);
		}

		return lines.size() + curves.size();
	}

	int get_intersections(
		Cylinder& C1, Cylinder& C2,
		std::vector<ParameterizationCylindricPoint>& points,
		std::vector<ParameterizationCylindricLine>& lines,
		std::vector<ParameterizationCylindricCurve>& curves
	) {
		SQI_VERBOSE_ONLY_COUT("");
		// SQI_VERBOSE_ONLY_COUT("cylinder1" << " " << "cor:" << C1.cor().transpose() << " " << "nor:" << C1.nor().transpose() << " " << "r:" << C1.r());
		// SQI_VERBOSE_ONLY_COUT("cylinder2" << " " << "cor:" << C2.cor().transpose() << " " << "nor:" << C2.nor().transpose() << " " << "r:" << C2.r());

		// init
		std::vector<ParameterizationCylindricPoint>().swap(points);
		std::vector<ParameterizationCylindricLine>().swap(lines);
		std::vector<ParameterizationCylindricCurve>().swap(curves);
		double C1_C2_sq_r = C1.r() + C2.r();
		C1_C2_sq_r *= C1_C2_sq_r;

		double C1_nor_dot_C2_nor = C1.nor().dot(C2.nor());

		if (std::abs(std::abs(C1_nor_dot_C2_nor) - 1) < SQI_EPS) { // 2 axes are parallel
			SQI_VERBOSE_ONLY_COUT("cylinders axes are parallel");

			// check whether intersect
			Eigen::Vector3d C2_to_C1 = C1.cor() - C2.cor();
			double C1_C2_dot_C1_nor_sq_dis = C2_to_C1.dot(C1.nor());
			C1_C2_dot_C1_nor_sq_dis *= C1_C2_dot_C1_nor_sq_dis;
			double axes_sq_dis = C2_to_C1.squaredNorm() - C1_C2_dot_C1_nor_sq_dis;
			// SQI_VERBOSE_ONLY_COUT("axes distance:" << std::sqrt(axes_sq_dis));

			if (axes_sq_dis > C1_C2_sq_r + SQI_EPS) { // not intersect
				SQI_VERBOSE_ONLY_COUT("not intersect");
				return 0;
			}

			if (axes_sq_dis < SQI_EPS) { // overlap
				SQI_VERBOSE_ONLY_COUT("cylinders are overlap, may error!");
				return 0;
			}
			else if (axes_sq_dis < C1_C2_sq_r - SQI_EPS) { // inside, result 2 lines
				SQI_VERBOSE_ONLY_COUT("intersect");

				Eigen::Vector3d center_p = 0.5 * (C1.cor() + C2.cor());
				Eigen::Vector3d perpendicular_v = (C1.nor().cross(C2_to_C1)).normalized();
				double move_dis = std::sqrt(C1.r() * C1.r() - 0.25 * axes_sq_dis);

				Eigen::Vector3d intersect_point1 = center_p + move_dis * perpendicular_v;
				Eigen::Vector3d intersect_point2 = center_p - move_dis * perpendicular_v;

				double s11, t11, s12, t12, s21, t21, s22, t22;
				C1.get_s_t(intersect_point1, s11, t11);
				C1.get_s_t(intersect_point2, s12, t12);

				lines.push_back(ParameterizationCylindricLine(t11));
				lines.push_back(ParameterizationCylindricLine(t12));
			}
			else { // tangent, result 1 line
				SQI_VERBOSE_ONLY_COUT("tangent");

				Eigen::Vector3d intersect_point = 0.5 * (C1.cor() + C2.cor());

				double s, t;
				C1.get_s_t(intersect_point, s, t);

				lines.push_back(ParameterizationCylindricLine(t));
			}
		}
		else { // 2 axes are not parallel
			SQI_VERBOSE_ONLY_COUT("cylinders axes are not parallel");

			// check whether intersect
			Eigen::Vector3d perpendicular_v = (C1.nor().cross(C2.nor())).normalized();
			double axes_dis = (C2.cor() - C1.cor()).dot(perpendicular_v);
			if (axes_dis > C1.r() + C2.r() + SQI_EPS) {
				SQI_VERBOSE_ONLY_COUT("not intersect");
				return 0;
			}

			// C1 is not being used as the parameterization cylinder, C2 is used as parameterization cylinder
			Line L1(
				C1.cor() + C1.r() * perpendicular_v,
				C1.nor());
			Line L2(
				C1.cor() - C1.r() * perpendicular_v,
				C1.nor());

			// check the intersections between lines and cylinder2, lines must not parallel to cylinder2's axis
			std::vector<ParameterizationCylindricPoint> L1_C2_intersect_points;
			int L1_status = get_intersections(L1, C2, L1_C2_intersect_points); // the number of intersect points

			std::vector<ParameterizationCylindricPoint> L2_C2_intersect_points;
			int L2_status = get_intersections(L2, C2, L2_C2_intersect_points); // the number of intersect points

			if (L1_status == 1 && L2_status == 0) { // tangent + do not intersect: result a point
				points.push_back(L1_C2_intersect_points[0]);
			}
			else if (L1_status == 0 && L2_status == 1) { // do not intersect + tangent: result a point
				points.push_back(L2_C2_intersect_points[0]);
			}
			else { // may result parameterization curves

				// compute a_t
				double C1_nor_dot_C2_nor = C1.nor().dot(C2.nor());
				std::vector<double> a_t = {
					1 - C1_nor_dot_C2_nor * C1_nor_dot_C2_nor
				};
				SQI_VERBOSE_ONLY_COUT("a_t:" << " " << a_t[0]);

				// compute b_t
				Eigen::Vector3d b_tmp = 2 * (C2.nor() - C1_nor_dot_C2_nor * C1.nor());
				Eigen::Vector3d C2_cor_C1_cor = C2.cor() - C1.cor();
				std::vector<double> b_t_C2 = {
					C2.r() * b_tmp.dot(C2.u()),
					C2.r() * b_tmp.dot(C2.v()),
					b_tmp.dot(C2_cor_C1_cor)
				};
				SQI_VERBOSE_ONLY_COUT("b_t:" << " " << b_t_C2[0] << " " << b_t_C2[1] << " " << b_t_C2[2]);

				//compute c_t
				double C2_sq_r = C2.r() * C2.r();
				double C2_u_dot_C1_nor = C2.u().dot(C1.nor());
				double C2_v_dot_C1_nor = C2.v().dot(C1.nor());
				double C2_nor_C1_nor_dot_C2_u = C2_cor_C1_cor.dot(C2.u());
				double C2_nor_C1_nor_dot_C2_v = C2_cor_C1_cor.dot(C2.v());
				double C2_nor_C1_nor_dot_C1_nor = C2_cor_C1_cor.dot(C1.nor());
				std::vector<double> c_t_C2 = {
					C2_sq_r * (1 - C2_u_dot_C1_nor * C2_u_dot_C1_nor),
					C2_sq_r * (1 - C2_v_dot_C1_nor * C2_v_dot_C1_nor),
					-2 * C2_sq_r * C2_u_dot_C1_nor * C2_v_dot_C1_nor,
					2 * C2.r() * (C2_nor_C1_nor_dot_C2_u - C2_nor_C1_nor_dot_C1_nor * C2_u_dot_C1_nor),
					2 * C2.r() * (C2_nor_C1_nor_dot_C2_v - C2_nor_C1_nor_dot_C1_nor * C2_v_dot_C1_nor),
					C2_cor_C1_cor.dot(C2_cor_C1_cor) - C2_nor_C1_nor_dot_C1_nor * C2_nor_C1_nor_dot_C1_nor - C1.r() * C1.r()
				};
				SQI_VERBOSE_ONLY_COUT("c_t:" << " " << c_t_C2[0] << " " << c_t_C2[1] << " " << c_t_C2[2] << " " << c_t_C2[3] << " " << c_t_C2[4] << " " << c_t_C2[5]);

				// build parameterization curves
				ParameterizationCylindricCurve PC_C2(
					a_t, b_t_C2, c_t_C2,
					-SQI_INFTY, SQI_INFTY, 0, 360
				);

				// try to cut this parameterization curve
				if (L1_status == 1 && L2_status == 1) { // two tangent: result two ellipse
					double t1, t2;
					t1 = L1_C2_intersect_points[0].t();
					t2 = L2_C2_intersect_points[0].t();

					if (t1 < t2) { // sort to let t1 > t2
						double tmp_t = t1;
						t1 = t2; t2 = tmp_t;
					}

					// cut curves into 4 (or more) parts
					ParameterizationCylindricCurve PC_C2_1, PC_C2_2;
					if (PC_C2.separate_t(PC_C2_1, t2) == 1) {
						ParameterizationCylindricCurve PC_C2_1s;
						PC_C2_1.separate_s(PC_C2_1s);

						curves.push_back(PC_C2_1); curves.push_back(PC_C2_1s);
					}
					if (PC_C2.separate_t(PC_C2_2, t1) == 1) {
						ParameterizationCylindricCurve PC_C2_2s;
						PC_C2_2.separate_s(PC_C2_2s);

						curves.push_back(PC_C2_2); curves.push_back(PC_C2_2s);
					}
					ParameterizationCylindricCurve PC_C2_s;
					PC_C2.separate_s(PC_C2_s);

					curves.push_back(PC_C2); curves.push_back(PC_C2_s);
				}
				else if (L1_status == 2 || L2_status == 2) { // two-point intersection + do not intersect: result 1 parameterization curve
					std::vector<ParameterizationCylindricPoint>* intersection_points;
					if (L1_status == 2) intersection_points = &L1_C2_intersect_points;
					if (L2_status == 2) intersection_points = &L2_C2_intersect_points;

					double t1 = (*intersection_points)[0].t();
					double t2 = (*intersection_points)[1].t();

					if (t2 < t1) { // sort to let t1 < t2
						double tmp_t = t1;
						t1 = t2; t2 = tmp_t;
					}

					if (axes_dis > C1.r() + SQI_EPS) { // circle center angle < 180
						if (t2 - t1 > 180) { // set to [t2, t1+360]
							PC_C2.t_lb() = t2; PC_C2.t_ub() = t1 + 360;
						}
						else { // set to [t1, t2]
							PC_C2.t_lb() = t1; PC_C2.t_ub() = t2;
						}
					}
					else if (axes_dis < C1.r() - SQI_EPS) { // circle center angle > 180
						if (t2 - t1 > 180) { // set to [t1, t2]
							PC_C2.t_lb() = t1; PC_C2.t_ub() = t2;
						}
						else { // set to [t2, t1+360]
							PC_C2.t_lb() = t2; PC_C2.t_ub() = t1 + 360;
						}
					}

					// cut curves into 4 parts
					ParameterizationCylindricCurve PC_C2_mid; // [lb, center]
					PC_C2.separate_t(PC_C2_mid, 0.5 * (PC_C2.t_lb() + PC_C2.t_ub())); // PC_C2 [center, ub]

					// all set to [0, 360]
					if (PC_C2.t_lb() > 360) {
						PC_C2.t_lb() -= 360; PC_C2.t_ub() -= 360;
					}
					else if (PC_C2.t_ub() > 360) { // cut to [lb, 360] and [0, ub - 360]
						ParameterizationCylindricCurve PC_C2_1 = PC_C2;
						PC_C2.t_ub() = 360;
						PC_C2_1.t_lb() = 0; PC_C2_1.t_ub() -= 360;
						curves.push_back(PC_C2_1);
					}

					if (PC_C2_mid.t_lb() > 360) {
						PC_C2_mid.t_lb() -= 360; PC_C2_mid.t_ub() -= 360;
					}
					else if (PC_C2_mid.t_ub() > 360) { // cut to [lb, 360] and [0, ub - 360]
						ParameterizationCylindricCurve PC_C2_mid_1 = PC_C2_mid;
						PC_C2_mid.t_ub() = 360;
						PC_C2_mid_1.t_lb() = 0; PC_C2_mid_1.t_ub() -= 360;
						curves.push_back(PC_C2_mid_1);
					}
					curves.push_back(PC_C2_mid);
					
					// cut according to s
					for (int i = 0, i_end = curves.size(); i < i_end; ++i) {
						ParameterizationCylindricCurve PC_s;
						curves[i].separate_s(PC_s);
						curves.push_back(PC_s);
					}
					ParameterizationCylindricCurve PC_s;
					PC_C2.separate_s(PC_s);
					curves.push_back(PC_C2); curves.push_back(PC_s);
				}
				else {
					SQI_VERBOSE_ONLY_COUT("may error?");
				}
			}
		}

		return points.size() + lines.size() + curves.size();
	}

	int get_intersections(
		Sphere& S1, Cylinder& C1,
		std::vector<ParameterizationCylindricPoint>& points,
		std::vector<ParameterizationCylindricCurve>& curves
	) {
		SQI_VERBOSE_ONLY_COUT("");
		// SQI_VERBOSE_ONLY_COUT("cylinder1" << " " << "cor:" << C1.cor().transpose() << " " << "nor:" << C1.nor().transpose() << " " << "r:" << C1.r());
		// SQI_VERBOSE_ONLY_COUT("sphere1" << " " << "cor:" << S1.cor().transpose() << " " << "r:" << S1.r());

		// init
		std::vector<ParameterizationCylindricPoint>().swap(points);
		std::vector<ParameterizationCylindricCurve>().swap(curves);

		double r1r2_sq_dis = (C1.r() + S1.r()) * (C1.r() + S1.r());

		Eigen::Vector3d C1_cor_to_S1_cor = S1.cor() - C1.cor();
		double center_to_axis_dot_dis = C1_cor_to_S1_cor.dot(C1.nor());
		double center_to_axis_sq_dis = C1_cor_to_S1_cor.squaredNorm() - center_to_axis_dot_dis * center_to_axis_dot_dis;
		Eigen::Vector3d S1_cor_proj_C1_nor_point = C1.cor() + center_to_axis_dot_dis * C1.nor();

		if (center_to_axis_sq_dis > r1r2_sq_dis + SQI_EPS) { // not intersect
			SQI_VERBOSE_ONLY_COUT("not intersect");
		}
		else if (center_to_axis_sq_dis > r1r2_sq_dis - SQI_EPS) { // tangent: result a point
			SQI_VERBOSE_ONLY_COUT("tangent");

			double s, t;
			C1.get_s_t(S1.cor(), s, t);
			points.push_back(ParameterizationCylindricPoint(s, t));
		}
		else if (center_to_axis_sq_dis < SQI_EPS) { // overlap: result a circle
			SQI_VERBOSE_ONLY_COUT("overlap");

			std::vector<double> a_t = { 0 };
			std::vector<double> b_t = { 0, 0, 1 };
			std::vector<double> c_t = { 0, 0, 0, 0, 0, -center_to_axis_dot_dis };
			ParameterizationCylindricCurve PC(
				a_t, b_t, c_t,
				-SQI_INFTY, SQI_INFTY, 0, 360
			);
			curves.push_back(PC);
		}
		else { // intersect: result a parameterization curve
			SQI_VERBOSE_ONLY_COUT("intersect");

			std::vector<double> a_t = {
					1
			};
			SQI_VERBOSE_ONLY_COUT("a_t:" << " " << a_t[0]);

			std::vector<double> b_t = {
					0,
					0,
					-2 * C1_cor_to_S1_cor.dot(C1.nor())
			};
			SQI_VERBOSE_ONLY_COUT("b_t:" << " " << b_t[0] << " " << b_t[1] << " " << b_t[2]);

			double C1_sq_r = C1.r() * C1.r();
			double S1_sq_r = S1.r() * S1.r();
			std::vector<double> c_t = {
					C1_sq_r,
					C1_sq_r,
					0,
					-2 * C1.r() * C1_cor_to_S1_cor.dot(C1.u()),
					-2 * C1.r() * C1_cor_to_S1_cor.dot(C1.v()),
					C1_cor_to_S1_cor.squaredNorm() - S1_sq_r
			};
			SQI_VERBOSE_ONLY_COUT("c_t:" << " " << c_t[0] << " " << c_t[1] << " " << c_t[2] << " " << c_t[3] << " " << c_t[4] << " " << c_t[5]);

			ParameterizationCylindricCurve PC(
				a_t, b_t, c_t,
				-SQI_INFTY, SQI_INFTY, 0, 360
			);

			// try to cut this parameterization curves
			double s, t;
			Eigen::Vector3d S1_cor_proj_C1_point = S1_cor_proj_C1_nor_point + C1.r() * (S1.cor() - S1_cor_proj_C1_nor_point).normalized();
			C1.get_s_t(S1_cor_proj_C1_point, s, t);
			double rot_angle = safetyAcos(0.5 * std::sqrt(center_to_axis_sq_dis) / C1.r()) + SQI_EPS;
			// SQI_VERBOSE_ONLY_COUT("test" << " " << s << " " << t << " " << rot_angle);
			double t1 = t - rot_angle;
			double t2 = t + rot_angle;

			if (t1 < 0) { // divide to [0, t2] and [t1 + 360, 360]
				ParameterizationCylindricCurve PC_tmp1, PC_tmp2;
				if (PC.separate_t(PC_tmp1, t2) == 1) {
					ParameterizationCylindricCurve PC_tmp1_s;
					PC_tmp1.separate_s(PC_tmp1_s);
					curves.push_back(PC_tmp1); curves.push_back(PC_tmp1_s);
				}
				if (PC.separate_t(PC_tmp2, t1 + 360) == 1) {
					ParameterizationCylindricCurve PC_s;
					PC.separate_s(PC_s);
					curves.push_back(PC); curves.push_back(PC_s);
				}
			}
			else if (t2 > 360) { // divide to [t1, 360] and [0, t2-360]
				ParameterizationCylindricCurve PC_tmp1, PC_tmp2;
				if (PC.separate_t(PC_tmp2, t2 - 360) == 1) {
					ParameterizationCylindricCurve PC_tmp2_s;
					PC_tmp2.separate_s(PC_tmp2_s);
					curves.push_back(PC_tmp2); curves.push_back(PC_tmp2_s);
				}
				if (PC.separate_t(PC_tmp1, t1) == 1) {
					ParameterizationCylindricCurve PC_s;
					PC.separate_s(PC_s);
					curves.push_back(PC); curves.push_back(PC_s);
				}
			}
			else { // divide to [t1, t2]
				ParameterizationCylindricCurve PC_tmp1, PC_tmp2;
				PC.separate_t(PC_tmp1, t1); // if failed, is also ok
				if (PC.separate_t(PC_tmp2, t2) == 1) {
					ParameterizationCylindricCurve PC_tmp2_s;
					PC_tmp2.separate_s(PC_tmp2_s);
					curves.push_back(PC_tmp2); curves.push_back(PC_tmp2_s);
				}
				else {
					ParameterizationCylindricCurve PC_s;
					PC.separate_s(PC_s);
					curves.push_back(PC); curves.push_back(PC_s);
				}
			}
		}

		return points.size() + curves.size();
	}
}