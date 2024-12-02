#ifndef SIMPLE_QUADRICS_INTERSECTION_HPP_
#define SIMPLE_QUADRICS_INTERSECTION_HPP_

#include <iostream>
#include <vector>
#include <math.h>
// Eigen
#include <Eigen/Dense>

#define SIMPLE_QUADRICS_INTERSECTION_VERBOSE_
#ifdef SIMPLE_QUADRICS_INTERSECTION_VERBOSE_
#define SQI_VERBOSE_ONLY_COUT(x) std::cout << "[" << __FUNCTION__ << "]" << " " << x << std::endl
#else
#define SQI_VERBOSE_ONLY_COUT(x)
#endif

#define SQI_EPS	1e-12
#define SQI_INFTY 1e12
#define SQI_PI 3.14159265

namespace QuadricsIntersection {
	// ----------------------------other functions----------------------------
	
	// get an arbitrary unit vector that is orthogonal to the given vector v
	inline Eigen::Vector3d get_perpendicular_normal(
		Eigen::Vector3d& v
	) {
		Eigen::Vector3d n(0, 0, 0);
		if (abs(v.x()) < SQI_EPS) {
			n.x() = 1;
		}
		else {
			if (abs(v.y()) < SQI_EPS) {
				n.y() = 1;
			}
			else {
				if (abs(v.z()) < SQI_EPS) {
					n.z() = 1;
				}
				else {
					n.x() = 1;
					n.z() = -v.x() / v.z();
					n.normalize();
				}
			}
		}

		return n;
	}

	// rotate a vector f by angle ang about the axis r
	inline Eigen::Vector3d slerp(
		Eigen::Vector3d& f,
		Eigen::Vector3d& r,
		double ang
	) {
		if (ang < SQI_EPS) {
			return f;
		}
		double rad_ang = ang / 180 * SQI_PI;
		double cos_ang = std::cos(rad_ang);
		double sin_ang = std::sin(rad_ang);
		return (cos_ang * f + (1 - cos_ang) * f.dot(r) * r + sin_ang * r.cross(f)).normalized();
	}

	// ----------------------------the quadric primitives----------------------------
	class Plane {
	public:
		Plane() {};
		~Plane() {};
		Plane(Eigen::Vector3d cor, Eigen::Vector3d nor) {
			cor_ = cor; nor_ = nor.normalized();
		}

		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }

		void output_model(
			std::vector<Eigen::Vector3d>& points,
			std::vector<Eigen::Vector3i>& faces,
			double w = 1., double h = 1.
		) {
			SQI_VERBOSE_ONLY_COUT("");

			// init
			std::vector<Eigen::Vector3d>().swap(points);
			std::vector<Eigen::Vector3i>().swap(faces);
			points.reserve(4);
			faces.reserve(2);

			Eigen::Vector3d u = get_perpendicular_normal(nor_);
			Eigen::Vector3d v = nor_.cross(u).normalized();
			points.push_back(cor_ + w * u + h * v);
			points.push_back(cor_ + w * u - h * v);
			points.push_back(cor_ - w * u + h * v);
			points.push_back(cor_ - w * u - h * v);
			faces.push_back(Eigen::Vector3i(0, 2, 1));
			faces.push_back(Eigen::Vector3i(1, 2, 3));
		}
	private:
		Eigen::Vector3d cor_; // a point on the plane
		Eigen::Vector3d nor_; // the unit normal of the plane
	};

	class Cylinder {
	public:
		Cylinder() {};
		~Cylinder() {};
		Cylinder(Eigen::Vector3d cor, Eigen::Vector3d nor, double r) {
			cor_ = cor; nor_ = nor.normalized(); r_ = abs(r);
			u_ = get_perpendicular_normal(nor);
			v_ = nor.cross(u_).normalized();
		}
		Eigen::Vector3d get_point(double s, double t) {
			double rad_t = t * SQI_PI / 180;
			return cor_ + r_ * (std::cos(rad_t) * u_ + std::sin(rad_t) * v_) + s * nor_;
		}
		void get_s_t(Eigen::Vector3d p, double& s, double& t) {
			s = (p - cor_).dot(nor_);
			// ...
		}

		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }
		double& r() { return r_; }
		Eigen::Vector3d& u() { return u_; }
		Eigen::Vector3d& v() { return v_; }

		void output_model(
			std::vector<Eigen::Vector3d>& points,
			std::vector<Eigen::Vector3i>& faces,
			int seg = 32, double h = 5.
		) {
			if (seg < 3) {
				SQI_VERBOSE_ONLY_COUT("input seg num is invalid!");
				return;
			}
			SQI_VERBOSE_ONLY_COUT("");

			// init
			std::vector<Eigen::Vector3d>().swap(points);
			std::vector<Eigen::Vector3i>().swap(faces);
			points.reserve(seg * 2);
			faces.reserve(seg * 2);

			Eigen::Vector3d u = get_perpendicular_normal(nor_);
			double rot_angle = 360. / seg;
			for (int i = 0; i < seg; ++i) {
				Eigen::Vector3d ru = slerp(u, nor_, i * rot_angle);
				points.push_back(cor_ - 0.5 * h * nor_ + ru);
				points.push_back(cor_ + 0.5 * h * nor_ + ru);
			}
			for (int i = 0; i < seg; ++i) {
				int ri = 2 * i;
				if (i == seg - 1) {
					faces.push_back(Eigen::Vector3i(ri, 0, ri + 1));
					faces.push_back(Eigen::Vector3i(ri + 1, 0, 1));
				}
				else {
					faces.push_back(Eigen::Vector3i(ri, ri + 2, ri + 1));
					faces.push_back(Eigen::Vector3i(ri + 1, ri + 2, ri + 3));
				}
			}
		}
	private:
		Eigen::Vector3d cor_; // a point at the cylinder's axis
		Eigen::Vector3d nor_; // the unit vector of the cylinder's axis
		double r_; // the radius of the cylinder
		Eigen::Vector3d u_, v_; // the coordinate axes of the plane which perpendicular to the cylinder's axis, used for parameterization representation
	};

	class Sphere {
	public:
		Sphere() {};
		~Sphere() {};
		Sphere(Eigen::Vector3d cor, double r) {
			cor_ = cor; r_ = abs(r);
		}

		Eigen::Vector3d& cor() { return cor_; }
		double& r() { return r_; }
	private:
		Eigen::Vector3d cor_; // the center of the sphere
		double r_; // the radius of the sphere
	};

	// ----------------------------the intersection primitives----------------------------
	class Point {
	public:
		Point() {};
		~Point() {};
		Point(Eigen::Vector3d cor) {
			cor_ = cor;
		}

		Eigen::Vector3d& cor() { return cor_; }
	private:
		Eigen::Vector3d cor_;
	};
	
	class Line {
	public:
		Line() {};
		~Line() {};
		Line(Eigen::Vector3d cor, Eigen::Vector3d nor) {
			cor_ = cor; nor_ = nor.normalized();
		}

		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }

		void output_model(
			std::vector<Eigen::Vector3d>& points,
			std::vector<Eigen::Vector2i>& lines,
			double l = 5.
		) {
			SQI_VERBOSE_ONLY_COUT("");

			// init
			std::vector<Eigen::Vector3d>().swap(points);
			std::vector<Eigen::Vector2i>().swap(lines);
			points.reserve(2);
			lines.reserve(1);

			points.push_back(cor_ + 0.5 * l * nor_);
			points.push_back(cor_ - 0.5 * l * nor_);
			lines.push_back(Eigen::Vector2i(0, 1));
		}
	private:
		Eigen::Vector3d cor_; // a point at the line
		Eigen::Vector3d nor_; // the vector of the line
	};

	class Circle {
	public:
		Circle() {};
		~Circle() {};
		Circle(Eigen::Vector3d cor, Eigen::Vector3d nor, double r) {
			cor_ = cor; nor_ = nor.normalized(); r_ = r;
		}

		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }
		double& r() { return r_; }

	private:
		Eigen::Vector3d cor_; // the center of the circle
		Eigen::Vector3d nor_; // the vector that perpendicular to the plane that the circle lies
		double r_;
	};

	class ParameterizationCurve {
		// a(t) * s^2 + b(t) * s + c(t) = 0
	public:
		ParameterizationCurve() {
			std::vector<double>(1, 0).swap(a_t_);
			std::vector<double>(3, 0).swap(b_t_);
			std::vector<double>(6, 0).swap(c_t_);
			s_lb_ = -1; s_ub_ = 1;
			t_lb_ = -1; t_ub_ = 1;
		}
		~ParameterizationCurve() {
			std::vector<double>().swap(a_t_);
			std::vector<double>().swap(b_t_);
			std::vector<double>().swap(c_t_);
		};
		ParameterizationCurve(
			std::vector<double>& a_t, std::vector<double>& b_t, std::vector<double>& c_t,
			double s_lb, double s_ub, double t_lb, double t_ub
		) {
			a_t_ = a_t; b_t_ = b_t; c_t_ = c_t;
			s_lb_ = s_lb; s_ub_ = s_ub; t_lb_ = t_lb; t_ub_ = t_ub;
		}
		int get_s(
			double t,
			std::vector<double>& s
		) {
			// init
			std::vector<double>().swap(s);

			double rad_t = t * SQI_PI / 180.;
			double cos_t = std::cos(rad_t);
			double sin_t = std::sin(rad_t);

			double a = a_t_[0];
			double b = b_t_[0] * cos_t + b_t_[1] * sin_t + b_t_[2];
			double c = c_t_[0] * cos_t * cos_t + c_t_[1] * sin_t * sin_t + c_t_[2] * sin_t * cos_t + c_t_[3] * cos_t + c_t_[4] * sin_t + c_t_[5];
			double delta = b * b - 4 * a * c;
			if (delta < -SQI_EPS) {
			}
			else if (delta < SQI_EPS) {
				s.push_back(-b / (2 * a));
			}
			else {
				double sqrt_delta = sqrt(delta);
				s.push_back((-b + sqrt_delta) / (2 * a));
				s.push_back((-b - sqrt_delta) / (2 * a));
			}

			return s.size();
		}
		std::vector<double>& a_t() { return a_t_; }
		std::vector<double>& b_t() { return b_t_; }
		std::vector<double>& c_t() { return c_t_; }
		double& s_lb() { return s_lb_; }
		double& s_ub() { return s_ub_; }
		double& t_lb() { return t_lb_; }
		double& t_ub() { return t_ub_; }

		void output_model(
			Cylinder& C,
			std::vector<Eigen::Vector3d>& points,
			double t_step = 10.
		) {
			double t_lb = t_lb_;
			double t_ub = t_ub_;

			for (double cur_t = t_lb; cur_t < t_ub; cur_t += t_step) {
				std::vector<double> ss;
				
				get_s(cur_t, ss);
				for (double s : ss) {
					points.push_back(C.get_point(s, cur_t));
				}
			}
		}
	private:
		std::vector<double> a_t_; // = [0]
		std::vector<double> b_t_; // = [0] * cos(t) + [1] * sin(t) + [2]
		std::vector<double> c_t_; // = [0] * cos(t)^2 + [1] * sin(t)^2 + [2] * sin(t) * cos(t) + [3] * cos(t) + [4] * sin(t) + [5]
		double s_lb_, s_ub_, t_lb_, t_ub_;
	};

	// ----------------------------assessment of the intersections----------------------------
	// get the branches that primitives intersected, return the number of branches
	int get_intersections(
		Line L1, Cylinder C1,
		std::vector<Point>& points
	) {
		SQI_VERBOSE_ONLY_COUT("line" << " " << "cor:" << L1.cor().transpose() << " " << "nor:" << L1.nor().transpose());
		SQI_VERBOSE_ONLY_COUT("cylinder" << " " << "cor:" << C1.cor().transpose() << " " << "nor:" << C1.nor().transpose() << " " << "r:" << C1.r());

		// init
		std::vector<Point>().swap(points);
		points.reserve(2);

		if (abs(abs(L1.nor().dot(C1.nor())) - 1) < SQI_EPS) { // line vector is parallel to the cylinder's axis
			SQI_VERBOSE_ONLY_COUT("axis parallel");
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

				double t = -b / (2 * a);
				points.push_back(
					Point(
						L1.cor() + t * L1.nor()
					)
				);
			}
			else {
				SQI_VERBOSE_ONLY_COUT("two-point intersection");

				double sqrt_delta = sqrt(delta);
				double t1 = (-b + sqrt_delta) / (2 * a);
				double t2 = (-b - sqrt_delta) / (2 * a);
				points.push_back(
					Point(
						L1.cor() + t1 * L1.nor()
					)
				);
				points.push_back(
					Point(
						L1.cor() + t2 * L1.nor()
					)
				);
			}
		}

		return points.size();
	}
	int get_intersections(
		Cylinder C1, Line L1,
		std::vector<Point>& points
	) {
		return get_intersections(L1, C1, points);
	}

	int get_intersections(
		Cylinder C1, Cylinder C2,
		std::vector<Point>& points,
		std::vector<Line>& lines,
		std::vector<Circle>& circles,
		std::vector<ParameterizationCurve>& parameterization_curves
	) {
		SQI_VERBOSE_ONLY_COUT("cylinder1" << " " << "cor:" << C1.cor().transpose() << " " << "nor:" << C1.nor().transpose() << " " << "r:" << C1.r());
		SQI_VERBOSE_ONLY_COUT("cylinder2" << " " << "cor:" << C2.cor().transpose() << " " << "nor:" << C2.nor().transpose() << " " << "r:" << C2.r());

		// init
		std::vector<Point>().swap(points);
		std::vector<Line>().swap(lines);
		std::vector<Circle>().swap(circles);
		std::vector<ParameterizationCurve>().swap(parameterization_curves);
		double C1_C2_sq_r = C1.r() + C2.r();
		C1_C2_sq_r *= C1_C2_sq_r;

		// check whether intersect
		Eigen::Vector3d C2_to_C1 = C1.cor() - C2.cor();
		double C1_C2_dot_sq_dis = C2_to_C1.dot(C1.nor());
		C1_C2_dot_sq_dis *= C1_C2_dot_sq_dis;
		double axes_sq_dis = C2_to_C1.squaredNorm() - C1_C2_dot_sq_dis;
		if (axes_sq_dis > C1_C2_sq_r + SQI_EPS) { // not intersect
			SQI_VERBOSE_ONLY_COUT("not intersect");
			return 0;
		}
		double C1_nor_dot_C2_nor = C1.nor().dot(C2.nor());

		if (abs(abs(C1_nor_dot_C2_nor) - 1) < SQI_EPS) { // 2 axes are parallel
			SQI_VERBOSE_ONLY_COUT("cylinders axes are parallel");

			if (axes_sq_dis < SQI_EPS) { // overlap
				SQI_VERBOSE_ONLY_COUT("cylinders are overlap, may error!");
				return 0;
			}
			else if (axes_sq_dis < C1_C2_sq_r - SQI_EPS) { // inside, result 2 lines
				SQI_VERBOSE_ONLY_COUT("intersect");

				Eigen::Vector3d center_p = 0.5 * (C1.cor() + C2.cor());
				Eigen::Vector3d perpendicular_v = (C1.nor().cross(C2_to_C1)).normalized();
				double move_dis = sqrt(C1.r() * C1.r() - 0.25 * axes_sq_dis);

				lines.push_back(
					Line(
						center_p + move_dis * perpendicular_v,
						C1.nor()
					));
				lines.push_back(
					Line(
						center_p - move_dis * perpendicular_v,
						C1.nor()
					));
			}
			else { // tangent, result 1 line
				SQI_VERBOSE_ONLY_COUT("tangent");

				lines.push_back(
					Line(
						0.5 * (C1.cor() + C2.cor()),
						C1.nor()
					));
			}
		}
		else { // 2 axes are not parallel
			SQI_VERBOSE_ONLY_COUT("cylinders axes are not parallel");

			// C1 is not being used as the parameterization cylinder, C2 is used as parameterization cylinder
			Eigen::Vector3d perpendicular_v = (C1.nor().cross(C2.nor())).normalized();
			Line L1(
				C1.cor() + C1.r() * perpendicular_v, 
				C1.nor());
			Line L2(
				C1.cor() - C1.r() * perpendicular_v,
				C1.nor());

			// check the intersections between lines and cylinder2, lines must not parallel to cylinder2's axis
			std::vector<Point> L1_C2_intersect_points;
			int L1_status = get_intersections(L1, C2, L1_C2_intersect_points); // the number of intersecttttt points

			std::vector<Point> L2_C2_intersect_points;
			int L2_status = get_intersections(L2, C2, L2_C2_intersect_points); // the number of intersecttttt points

			if (L1_status == 1 && L2_status == 0) { // tangent + do not intersect: result a point
				points.push_back(L1_C2_intersect_points[0]);
			}
			else if (L1_status == 0 && L2_status == 1) { // do not intersect + tangent: result a point
				points.push_back(L2_C2_intersect_points[0]);
			}
			else { // may result parameterization curves
				double C1_nor_dot_C2_nor = C1.nor().dot(C2.nor());
				std::vector<double> a_t = { 
					1 - C1_nor_dot_C2_nor * C1_nor_dot_C2_nor 
				};
				SQI_VERBOSE_ONLY_COUT("a_t:" << " " << a_t[0]);

				Eigen::Vector3d b_tmp = 2 * (C2.nor() - C1_nor_dot_C2_nor * C1.nor());
				std::vector<double> b_t = {
					C2.r() * b_tmp.dot(C2.u()),
					C2.r() * b_tmp.dot(C2.v()),
					b_tmp.dot(C2.cor() - C1.cor()) 
				};
				SQI_VERBOSE_ONLY_COUT("b_t:" << " " << b_t[0] << " " << b_t[1] << " " << b_t[2]);

				double C2_sq_r = C2.r() * C2.r();
				double C2_u_dot_C1_nor = C2.u().dot(C1.nor());
				double C2_v_dot_C1_nor = C2.v().dot(C1.nor());
				double C2_u_dot_C1_cor = C2.u().dot(C1.cor());
				double C2_v_dot_C1_cor = C2.v().dot(C1.cor());
				double C1_cor_dot_C1_nor = C1.cor().dot(C1.nor());
				std::vector<double> c_t = {
					C2_sq_r * (1 - C2_u_dot_C1_nor * C2_u_dot_C1_nor),
					C2_sq_r * (1 - C2_v_dot_C1_nor * C2_v_dot_C1_nor),
					-2 * C2_sq_r * C2_u_dot_C1_nor * C2_v_dot_C1_nor,
					2 * C2.r() * (C2_u_dot_C1_nor * C1_cor_dot_C1_nor - C2_u_dot_C1_cor),
					2 * C2.r() * (C2_v_dot_C1_nor * C1_cor_dot_C1_nor - C2_v_dot_C1_cor),
					C1.cor().dot(C1.cor()) - C1_cor_dot_C1_nor * C1_cor_dot_C1_nor
				};
				SQI_VERBOSE_ONLY_COUT("c_t:" << " " << c_t[0] << " " << c_t[1] << " " << c_t[2] << " " << c_t[3] << " " << c_t[4] << " " << c_t[5]);

				ParameterizationCurve PC(
					a_t, b_t, c_t,
					-SQI_INFTY, SQI_INFTY, 0, 360
				);
				parameterization_curves.push_back(PC);

				if (L1_status == 1 && L2_status == 1) { // two tangent: result two ellipse

				}
				else if (L1_status == 2) { // two-point intersection + do not intersect:

				}
				else if (L2_status == 2) { // do not intersect + two_point intersection:

				}
				else {
					SQI_VERBOSE_ONLY_COUT("may error?");
				}
			}
		}

		return points.size() + lines.size() + circles.size() + parameterization_curves.size();
	}
}

#endif
