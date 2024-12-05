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
#define SQI_PI 3.141592653589793
#define SQI_REC_PI 0.318309886183791

namespace QuadricsIntersection {
	// ----------------------------other functions----------------------------

	inline double rad2ang(double rad) {
		return rad * 180 * SQI_REC_PI;
	}

	inline double ang2rad(double ang) {
		return ang * 0.005555555555556 * SQI_PI;
	}

	inline double safetyAsin(double value) {
		if (value < -1.) {
			rad2ang(asin(-1.));
		}
		if (value > 1.) {
			rad2ang(asin(1.));
		}
		return rad2ang(asin(value));
	}

	inline double safetyAcos(double value) {
		if (value < -1.) {
			return rad2ang(acos(-1.));
		}
		if (value > 1.) {
			return rad2ang(acos(1.));
		}
		return rad2ang(acos(value));
	}

	// get an arbitrary unit vector that is orthogonal to the given vector v
	inline Eigen::Vector3d get_perpendicular_normal(
		Eigen::Vector3d& v
	) {
		Eigen::Vector3d n(0, 0, 0);
		if (std::abs(v.x()) < SQI_EPS) {
			n.x() = 1;
		}
		else {
			if (std::abs(v.y()) < SQI_EPS) {
				n.y() = 1;
			}
			else {
				if (std::abs(v.z()) < SQI_EPS) {
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
		double rad_ang = ang2rad(ang);
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
			cor_ = cor; nor_ = nor.normalized(); r_ = std::abs(r);
			u_ = get_perpendicular_normal(nor);
			v_ = nor.cross(u_).normalized();
		}
		Eigen::Vector3d get_point(double s, double t) {
			double rad_t = ang2rad(t);
			return cor_ + r_ * (std::cos(rad_t) * u_ + std::sin(rad_t) * v_) + s * nor_;
		}
		void get_s_t(Eigen::Vector3d p, double& s, double& t) {
			s = (p - cor_).dot(nor_);
			Eigen::Vector3d plane_cor = p - s * nor_;
			Eigen::Vector3d plane_cor_cor = plane_cor - cor_;
			double tu_component = plane_cor_cor.dot(u_) / r_;
			double tv_component = plane_cor_cor.dot(v_) / r_;
			t = safetyAcos(tu_component);
			if (std::abs(std::sin(ang2rad(t)) - tv_component) > SQI_EPS) {
				t = 360 - t;
			}
		}

		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }
		double& r() { return r_; }
		Eigen::Vector3d& u() { return u_; }
		Eigen::Vector3d& v() { return v_; }

		void output_model(
			std::vector<Eigen::Vector3d>& points,
			std::vector<Eigen::Vector3i>& faces,
			int seg = 32, double h = 20.
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

			double rot_angle = 360. / seg;
			for (int i = 0; i < seg; ++i) {
				Eigen::Vector3d ru = slerp(u_, nor_, i * rot_angle);
				points.push_back(cor_ - 0.5 * h * nor_ + r_ * ru);
				points.push_back(cor_ + 0.5 * h * nor_ + r_ * ru);
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
			cor_ = cor; r_ = std::abs(r);
		}

		Eigen::Vector3d& cor() { return cor_; }
		double& r() { return r_; }

		void output_model(
			std::vector<Eigen::Vector3d>& points,
			std::vector<Eigen::Vector3i>& faces,
			int h_seg = 32, int r_seg = 32
		) {
			if (h_seg < 3) {
				SQI_VERBOSE_ONLY_COUT("input h_seg num is invalid!");
				return;
			}
			if (r_seg < 3) {
				SQI_VERBOSE_ONLY_COUT("input r_seg num is invalid!");
				return;
			}
			SQI_VERBOSE_ONLY_COUT("");

			// init
			std::vector<Eigen::Vector3d>().swap(points);
			std::vector<Eigen::Vector3i>().swap(faces);
			points.reserve((h_seg - 2) * r_seg + 2);
			faces.reserve((h_seg - 2) * r_seg * 2);
			double h_angle = 180. / (h_seg - 1);
			double r_angle = 360. / r_seg;

			for (double phi = 0, phi_end = 180 + SQI_EPS; phi < phi_end; phi += h_angle) {
				if (phi < SQI_EPS) {
					points.push_back(cor_ + r_ * Eigen::Vector3d(0, 0, 1));
				}
				else if (std::abs(phi - 180) < SQI_EPS) {
					points.push_back(cor_ + r_ * Eigen::Vector3d(0, 0, -1));
				}
				else {
					for (double theta = 0, theta_end = 360 - SQI_EPS; theta < theta_end; theta += r_angle) {
						points.push_back(
							cor_ + Eigen::Vector3d(
								r_ * std::sin(ang2rad(phi)) * std::cos(ang2rad(theta)),
								r_ * std::sin(ang2rad(phi)) * std::sin(ang2rad(theta)),
								r_ * std::cos(ang2rad(phi))
							));
					}
				}
			}

			for (int h = 1; h < h_seg; ++h) {
				for (int r = 0; r < r_seg; ++r) {
					int cur_p = (h - 1) * r_seg + r + 1;

					if (h == 1) {
						if (r == r_seg - 1) {
							faces.push_back(Eigen::Vector3i(cur_p, (h - 1) * r_seg + 1, 0));
						}
						else {
							faces.push_back(Eigen::Vector3i(cur_p, cur_p + 1, 0));
						}
					}
					else if (h == h_seg - 1) {
						cur_p = (h - 2) * r_seg + r + 1;

						if (r == r_seg - 1) {
							faces.push_back(Eigen::Vector3i((h - 2) * r_seg + 1, cur_p, (h_seg - 2) * r_seg + 1));
						}
						else {
							faces.push_back(Eigen::Vector3i(cur_p + 1, cur_p, (h_seg - 2) * r_seg + 1));
						}
					}
					else {
						if (r == r_seg - 1) {
							faces.push_back(Eigen::Vector3i(cur_p, (h - 1) * r_seg + 1, cur_p - r_seg));
							faces.push_back(Eigen::Vector3i((h - 2) * r_seg + 1, cur_p - r_seg, (h - 1) * r_seg + 1));
						}
						else {
							faces.push_back(Eigen::Vector3i(cur_p, cur_p + 1, cur_p - r_seg));
							faces.push_back(Eigen::Vector3i(cur_p - r_seg + 1, cur_p - r_seg, cur_p + 1));
						}
					}
				}
			}
		}
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
			double l = 30.
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
		void output_model(
			std::vector<Eigen::Vector3d>& points,
			double l = 30., int seg = 100
		) {
			SQI_VERBOSE_ONLY_COUT("");

			// init
			std::vector<Eigen::Vector3d>().swap(points);

			double l_step = l / seg;
			for (double cur_l = -0.5 * l, cur_l_end = 0.5 * l + SQI_EPS; cur_l < cur_l_end; cur_l += l_step) {
				points.push_back(cor_ + cur_l * nor_);
			}
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
			t_lb_ = 0; t_ub_ = SQI_PI;
			s_part_ = 0;
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
			s_part_ = 0;
		}
		ParameterizationCurve(
			std::vector<double>& a_t, std::vector<double>& b_t, std::vector<double>& c_t,
			double s_lb, double s_ub, double t_lb, double t_ub,
			int s_part
		) {
			a_t_ = a_t; b_t_ = b_t; c_t_ = c_t;
			s_lb_ = s_lb; s_ub_ = s_ub; t_lb_ = t_lb; t_ub_ = t_ub;
			s_part_ = s_part;
		}
		ParameterizationCurve& operator =(ParameterizationCurve& PC) {
			if (this != &PC) {
				this->a_t_ = PC.a_t_;
				this->b_t_ = PC.b_t_;
				this->c_t_ = PC.c_t_;
				this->s_lb_ = PC.s_lb_; this->s_ub_ = PC.s_ub_; this->t_lb_ = PC.t_lb_; this->t_ub_ = PC.t_ub_;
				this->s_part_ = PC.s_part_;
			}
			return *this;
		}
		// divide the curve into two parts according to s
		void separate_s(ParameterizationCurve& PC) {
			PC = *this;
			this->s_part_ = 1;
			PC.s_part_ = -1;
		}
		// devide the curve into two parts according to t1-t2, modify self [t1, t2], return [t2, t1+360]
		void separate_t(ParameterizationCurve& PC, 
			double t1, double t2) {
			PC = *this;
			if (t1 < t2) { // sort to let t1 > t2
				double tmp_t = t1;
				t1 = t2; t2 = tmp_t;
			}
			this->t_lb_ = t2; this->t_ub_ = t1;
			PC.t_lb_ = t1; PC.t_ub_ = t2 + 360;
		}
		int get_s(
			double t,
			std::vector<double>& s) {
			// init
			std::vector<double>().swap(s);

			double rad_t = ang2rad(t);
			double cos_t = std::cos(rad_t);
			double sin_t = std::sin(rad_t);

			double a = a_t_[0];
			double b = b_t_[0] * cos_t + b_t_[1] * sin_t + b_t_[2];
			double c = c_t_[0] * cos_t * cos_t + c_t_[1] * sin_t * sin_t + c_t_[2] * sin_t * cos_t + c_t_[3] * cos_t + c_t_[4] * sin_t + c_t_[5];

			if (std::abs(a) < SQI_EPS) { // not a bisecond equation
				if (std::abs(b) < SQI_EPS) { // not a first-degree equation
					SQI_VERBOSE_ONLY_COUT("cannot get the value of s");
				}
				else { // is a firt-degree equation
					s.push_back(-c / b);
				}
			}
			else { // is a bisecond equation
				double delta = b * b - 4 * a * c;
				if (delta < -SQI_EPS) {
				}
				else if (delta < SQI_EPS) {
					s.push_back(-b / (2 * a));
				}
				else {
					double sqrt_delta = std::sqrt(delta);
					if (std::abs(s_part_) < SQI_EPS || s_part_ > 0) s.push_back((-b + sqrt_delta) / (2 * a));
					if (std::abs(s_part_) < SQI_EPS || s_part_ < 0) s.push_back((-b - sqrt_delta) / (2 * a));
				}
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
		int& s_part() { return s_part_; }

		void verbose() {
			SQI_VERBOSE_ONLY_COUT("a_t_:" << " " << a_t_[0]);
			SQI_VERBOSE_ONLY_COUT("b_t_:" << " " << b_t_[0] << " " << b_t_[1] << " " << b_t_[2]);
			SQI_VERBOSE_ONLY_COUT("c_t_:" << " " << c_t_[0] << " " << c_t_[1] << " " << c_t_[2] << " " << c_t_[3] << " " << c_t_[4] << " " << c_t_[5]);
			SQI_VERBOSE_ONLY_COUT("s_bound:" << " " << s_lb_ << "--" << s_ub_ << " " << "s_part:" << " " << s_part_);
			SQI_VERBOSE_ONLY_COUT("t_bound:" << " " << t_lb_ << "--" << t_ub_);
		}
		void output_model(
			Cylinder& C,
			std::vector<Eigen::Vector3d>& points,
			int t_cnt = 300
		) {
			SQI_VERBOSE_ONLY_COUT("");

			std::vector<Eigen::Vector3d>().swap(points);

			double t_step = (t_ub_ - t_lb_) / t_cnt;
			if (t_step < SQI_EPS) return;
			for (double cur_t = t_lb_, cur_t_end = t_ub_ + SQI_EPS; cur_t < cur_t_end; cur_t += t_step) {
				std::vector<double> ss;
				get_s(cur_t, ss);

				for (double s : ss) {
					if (s >= s_lb_ && s <= s_ub_) {
						points.push_back(C.get_point(s, cur_t));
					}
				}
			}
		}
	private:
		std::vector<double> a_t_; // = [0]
		std::vector<double> b_t_; // = [0] * cos(t) + [1] * sin(t) + [2]
		std::vector<double> c_t_; // = [0] * cos(t)^2 + [1] * sin(t)^2 + [2] * sin(t) * cos(t) + [3] * cos(t) + [4] * sin(t) + [5]
		double s_lb_, s_ub_, t_lb_, t_ub_;
		int s_part_; // decide which part of the curve is need, 0: all, 1: upper, -1: lower
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

		if (std::abs(std::abs(L1.nor().dot(C1.nor())) - 1) < SQI_EPS) { // line vector is parallel to the cylinder's axis
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

				double sqrt_delta = std::sqrt(delta);
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
		Plane P1, Plane P2,
		std::vector<Line>& lines
	) {
		SQI_VERBOSE_ONLY_COUT("plane1" << " " << "cor:" << P1.cor().transpose() << " " << "nor:" << P1.nor().transpose());
		SQI_VERBOSE_ONLY_COUT("plane2" << " " << "cor:" << P2.cor().transpose() << " " << "nor:" << P2.nor().transpose());

		// init
		std::vector<Line>().swap(lines);

		if (std::abs(std::abs(P1.nor().dot(P2.nor())) - 1) < SQI_EPS) {
			SQI_VERBOSE_ONLY_COUT("planes parallel");
		}
		else {
			// ...
		}

		return lines.size();
	}

	// using C1 as the parameterization surface
	int get_intersetions(
		Cylinder C1, Plane P1,
		std::vector<Point>& points,
		std::vector<Line>& lines,
		std::vector<Circle>& circles,
		std::vector<ParameterizationCurve>& parameterization_curves
	) {

	}

	// using C2 as the parameterization surface
	int get_intersections(
		Cylinder& C1, Cylinder& C2,
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

		double C1_nor_dot_C2_nor = C1.nor().dot(C2.nor());

		if (std::abs(std::abs(C1_nor_dot_C2_nor) - 1) < SQI_EPS) { // 2 axes are parallel
			SQI_VERBOSE_ONLY_COUT("cylinders axes are parallel");

			// check whether intersect
			Eigen::Vector3d C2_to_C1 = C1.cor() - C2.cor();
			double C1_C2_dot_C1_nor_sq_dis = C2_to_C1.dot(C1.nor());
			C1_C2_dot_C1_nor_sq_dis *= C1_C2_dot_C1_nor_sq_dis;
			double axes_sq_dis = C2_to_C1.squaredNorm() - C1_C2_dot_C1_nor_sq_dis;
			SQI_VERBOSE_ONLY_COUT("axes distance:" << std::sqrt(axes_sq_dis));
			SQI_VERBOSE_ONLY_COUT(C2_to_C1.transpose() << " " << C2_to_C1.norm());

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
			std::vector<Point> L1_C2_intersect_points;
			int L1_status = get_intersections(L1, C2, L1_C2_intersect_points); // the number of intersect points

			std::vector<Point> L2_C2_intersect_points;
			int L2_status = get_intersections(L2, C2, L2_C2_intersect_points); // the number of intersect points

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
				Eigen::Vector3d C2_cor_C1_cor = C2.cor() - C1.cor();
				std::vector<double> b_t = {
					C2.r() * b_tmp.dot(C2.u()),
					C2.r() * b_tmp.dot(C2.v()),
					b_tmp.dot(C2_cor_C1_cor)
				};
				SQI_VERBOSE_ONLY_COUT("b_t:" << " " << b_t[0] << " " << b_t[1] << " " << b_t[2]);

				double C2_sq_r = C2.r() * C2.r();
				double C2_u_dot_C1_nor = C2.u().dot(C1.nor());
				double C2_v_dot_C1_nor = C2.v().dot(C1.nor());
				double C2_nor_C1_nor_dot_C2_u = C2_cor_C1_cor.dot(C2.u());
				double C2_nor_C1_nor_dot_C2_v = C2_cor_C1_cor.dot(C2.v());
				double C2_nor_C1_nor_dot_C1_nor = C2_cor_C1_cor.dot(C1.nor());
				std::vector<double> c_t = {
					C2_sq_r * (1 - C2_u_dot_C1_nor * C2_u_dot_C1_nor),
					C2_sq_r * (1 - C2_v_dot_C1_nor * C2_v_dot_C1_nor),
					-2 * C2_sq_r * C2_u_dot_C1_nor * C2_v_dot_C1_nor,
					2 * C2.r() * (C2_nor_C1_nor_dot_C2_u - C2_nor_C1_nor_dot_C1_nor * C2_u_dot_C1_nor),
					2 * C2.r() * (C2_nor_C1_nor_dot_C2_v - C2_nor_C1_nor_dot_C1_nor * C2_v_dot_C1_nor),
					C2_cor_C1_cor.dot(C2_cor_C1_cor) - C2_nor_C1_nor_dot_C1_nor * C2_nor_C1_nor_dot_C1_nor - C1.r() * C1.r()
				};
				SQI_VERBOSE_ONLY_COUT("c_t:" << " " << c_t[0] << " " << c_t[1] << " " << c_t[2] << " " << c_t[3] << " " << c_t[4] << " " << c_t[5]);

				ParameterizationCurve PC(
					a_t, b_t, c_t,
					-SQI_INFTY, SQI_INFTY, 0, 360
				);

				// try to cut this parameterization curve
				if (L1_status == 1 && L2_status == 1) { // two tangent: result two ellipse
					double s1, t1, s2, t2;
					C2.get_s_t(L1_C2_intersect_points[0].cor(), s1, t1);
					C2.get_s_t(L2_C2_intersect_points[0].cor(), s2, t2);
					SQI_VERBOSE_ONLY_COUT("t1:" << t1 << " " << "t2:" << t2);

					// cut curves into 4 parts
					ParameterizationCurve PC2;
					PC.separate_t(PC2, t1, t2);

					ParameterizationCurve PC11;
					PC.separate_s(PC11);
					ParameterizationCurve PC21;
					PC2.separate_s(PC21);
					
					parameterization_curves.push_back(PC);
					parameterization_curves.push_back(PC2);
					parameterization_curves.push_back(PC11);
					parameterization_curves.push_back(PC21);
				}
				else if (L1_status == 2 || L2_status == 2) { // two-point intersection + do not intersect: result 1 parameterization curve
					std::vector<Point>* intersection_points;
					if (L1_status == 2) intersection_points = &L1_C2_intersect_points;
					if (L2_status == 2) intersection_points = &L2_C2_intersect_points;

					double s1, t1, s2, t2;
					C2.get_s_t((*intersection_points)[0].cor(), s1, t1);
					C2.get_s_t((*intersection_points)[1].cor(), s2, t2);

					if (t1 < t2) { // sort to let t1 > t2
						double tmp_t = t1;
						t1 = t2; t2 = tmp_t;
					}
					SQI_VERBOSE_ONLY_COUT("t1:" << t1 << " " << "t2:" << t2);

					if (axes_dis > C1.r() + SQI_EPS) { // circle center angle > 180
						if (t1 - t2 > 180) {
							PC.t_lb() = t1; PC.t_ub() = t2 + 360;
						}
						else {
							PC.t_lb() = t2; PC.t_ub() = t1;
						}
					}
					else if (axes_dis < C1.r() - SQI_EPS) { // circle center angle < 180
						if (t1 - t2 < 180) {
							PC.t_lb() = t1; PC.t_ub() = t2 + 360;
						}
						else {
							PC.t_lb() = t2; PC.t_ub() = t1;
						}
					}

					// cut curves into 2 parts
					ParameterizationCurve PC2;
					PC.separate_s(PC2);
					parameterization_curves.push_back(PC);
					parameterization_curves.push_back(PC2);
				}
				else {
					SQI_VERBOSE_ONLY_COUT("may error?");
				}
			}
		}

		return points.size() + lines.size() + circles.size() + parameterization_curves.size();
	}

	// using C1 as the parameterization surface
	int get_intersections(
		Cylinder& C1, Sphere& S1,
		std::vector<Point>& points,
		std::vector<Line>& lines,
		std::vector<Circle>& circles,
		std::vector<ParameterizationCurve>& parameterization_curves
	) {
		SQI_VERBOSE_ONLY_COUT("cylinder1" << " " << "cor:" << C1.cor().transpose() << " " << "nor:" << C1.nor().transpose() << " " << "r:" << C1.r());
		SQI_VERBOSE_ONLY_COUT("sphere1" << " " << "cor:" << S1.cor().transpose() << " " << "r:" << S1.r());

		// init
		std::vector<Point>().swap(points);
		std::vector<Line>().swap(lines);
		std::vector<Circle>().swap(circles);
		std::vector<ParameterizationCurve>().swap(parameterization_curves);
		
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
			
			points.push_back(Point(S1_cor_proj_C1_nor_point + C1.r() * (S1.cor() - S1_cor_proj_C1_nor_point).normalized()));
		}
		else if (center_to_axis_sq_dis < SQI_EPS) { // overlap: result a circle
			SQI_VERBOSE_ONLY_COUT("overlap");

			std::vector<double> a_t = { 0 };
			std::vector<double> b_t = { 0, 0, 1 };
			std::vector<double> c_t = { 0, 0, 0, 0, 0, -center_to_axis_dot_dis };
			ParameterizationCurve PC(
				a_t, b_t, c_t,
				-SQI_INFTY, SQI_INFTY, 0, 360
			);
			parameterization_curves.push_back(PC);
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

			ParameterizationCurve PC(
				a_t, b_t, c_t,
				-SQI_INFTY, SQI_INFTY, 0, 360
			);

			// try to cut this parameterization curves
			double s, t;
			Eigen::Vector3d S1_cor_proj_C1_point = S1_cor_proj_C1_nor_point + C1.r() * (S1.cor() - S1_cor_proj_C1_nor_point).normalized();
			C1.get_s_t(S1_cor_proj_C1_point, s, t);
			double rot_angle = safetyAcos(0.5 * std::sqrt(center_to_axis_sq_dis) / C1.r()) + SQI_EPS;
			SQI_VERBOSE_ONLY_COUT("test" << " " << s << " " << t << " " << rot_angle);
			double t1 = t - rot_angle;
			double t2 = t + rot_angle;

			ParameterizationCurve PC_tmp;
			PC.separate_t(PC_tmp, t1, t2);

			ParameterizationCurve PC2;
			PC.separate_s(PC2);

			parameterization_curves.push_back(PC);
			parameterization_curves.push_back(PC2);
		}

		return points.size() + lines.size() + circles.size() + parameterization_curves.size();
	}
	int get_intersections(
		Sphere& S1, Cylinder& C1,
		std::vector<Point>& points,
		std::vector<Line>& lines,
		std::vector<Circle>& circles,
		std::vector<ParameterizationCurve>& parameterization_curves
	) {
		return get_intersections(C1, S1, points, lines, circles, parameterization_curves);
	}

	// ----------------------------output for debug----------------------------

	// output primitives to an OBJ-type file
	void write_result_mesh(
		std::string output_file_path,
		std::vector<Plane>& planes,
		std::vector<Cylinder>& cylinders,
		std::vector<Sphere>& spheres
	) {
		SQI_VERBOSE_ONLY_COUT("");

		std::ofstream out(output_file_path);

		std::vector<Eigen::Vector3d> points;
		std::vector<Eigen::Vector3i> faces;
		int total_points_nb = 0;
		for (Plane P : planes) {
			P.output_model(points, faces);

			for (Eigen::Vector3d p : points) {
				out << "v" << " " << p.transpose() << std::endl;
			}
			for (Eigen::Vector3i f : faces) {
				out << "f" << " " << (f + Eigen::Vector3i(total_points_nb + 1, total_points_nb + 1, total_points_nb + 1)).transpose() << std::endl;
			}
			total_points_nb += points.size();
		}
		for (Cylinder C : cylinders) {
			C.output_model(points, faces);

			for (Eigen::Vector3d p : points) {
				out << "v" << " " << p.transpose() << std::endl;
			}
			for (Eigen::Vector3i f : faces) {
				out << "f" << " " << (f + Eigen::Vector3i(total_points_nb + 1, total_points_nb + 1, total_points_nb + 1)).transpose() << std::endl;
			}
			total_points_nb += points.size();
		}
		for (Sphere S: spheres) {
			S.output_model(points, faces);

			for (Eigen::Vector3d p : points) {
				out << "v" << " " << p.transpose() << std::endl;
			}
			for (Eigen::Vector3i f : faces) {
				out << "f" << " " << (f + Eigen::Vector3i(total_points_nb + 1, total_points_nb + 1, total_points_nb + 1)).transpose() << std::endl;
			}
			total_points_nb += points.size();
		}
		out.close();
	}

	// output points of the curves on the cylinder C
	void write_result_points(
		std::string output_file_path,
		Cylinder& C,
		std::vector<Point>& points,
		std::vector<Line>& lines,
		std::vector<ParameterizationCurve>& parameterization_curves
	) {
		SQI_VERBOSE_ONLY_COUT("");

		std::ofstream out(output_file_path);
		for (Point P : points) {
			Eigen::Vector3d color(
				double(std::rand()) / (RAND_MAX + 1) * 255,
				double(std::rand()) / (RAND_MAX + 1) * 255,
				double(std::rand()) / (RAND_MAX + 1) * 255);

			out << "v" << " " << P.cor().transpose() << " " << color.transpose() << std::endl;
		}
		for (Line L : lines) {
			Eigen::Vector3d color(
				double(std::rand()) / (RAND_MAX + 1) * 255,
				double(std::rand()) / (RAND_MAX + 1) * 255,
				double(std::rand()) / (RAND_MAX + 1) * 255);

			std::vector<Eigen::Vector3d> output_points;
			L.output_model(output_points);

			for (Eigen::Vector3d p : output_points) {
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
		}
		for (ParameterizationCurve PC : parameterization_curves) {
			Eigen::Vector3d color(
				double(std::rand()) / (RAND_MAX + 1) * 255,
				double(std::rand()) / (RAND_MAX + 1) * 255,
				double(std::rand()) / (RAND_MAX + 1) * 255);

			std::vector<Eigen::Vector3d> output_points;
			PC.output_model(C, output_points);

			for (Eigen::Vector3d p : output_points) {
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
		}
		out.close();
	}
}

#endif
