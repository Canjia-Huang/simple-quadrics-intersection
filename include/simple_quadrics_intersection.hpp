#ifndef SIMPLE_QUADRICS_INTERSECTION_HPP_
#define SIMPLE_QUADRICS_INTERSECTION_HPP_

#include <iostream>
#include <vector>
#include <math.h>
#include <functional>
#include <unordered_map>
// Eigen
#include <Eigen/Dense>

#define SIMPLE_QUADRICS_INTERSECTION_VERBOSE_
#ifdef SIMPLE_QUADRICS_INTERSECTION_VERBOSE_
#define SQI_VERBOSE_ONLY_COUT(x) std::cout << "\033[33m" << "[" << __FUNCTION__ << "]" << "\033[0m" << " " << x << std::endl // [yellow] white cout
#define SQI_VERBOSE_ONLY_WARNING(x) std::cout << "\033[33m" << "[" << __FUNCTION__ << "]" << "\033[0m" << " " << "\033[31m" << x << "\033[0m" << std::endl // [yellow] red cout
#else
#define SQI_VERBOSE_ONLY_COUT(x)
#define SQI_VERBOSE_ONLY_WARNING(x)
#endif

#define SQI_EPS	1e-12
#define SQI_INFTY 1e12
#define SQI_PI 3.141592653589793
#define SQI_REC_PI 0.318309886183791

namespace QuadricsIntersection {
	// ----------------------------other functions----------------------------

	struct pair_hash {
		template <class T1, class T2>
		size_t operator () (std::pair<T1, T2> const& pair
			) const {
			size_t h1 = std::hash<T1>()(pair.first);
			size_t h2 = std::hash<T2>()(pair.second);
			return h1 ^ h2;
		}
	};

	std::vector<Eigen::Vector3d> rand_color_bar;
	void create_rand_color_bar(int num) {
		rand_color_bar.reserve(num);
		for (int i = 0; i < num; ++i) {
			rand_color_bar.push_back(
				Eigen::Vector3d(
				double(std::rand()) / (RAND_MAX + 1) * 255,
				double(std::rand()) / (RAND_MAX + 1) * 255,
				double(std::rand()) / (RAND_MAX + 1) * 255));
		}
	}

	inline double rad2ang(double rad) {
		return rad * 180 * SQI_REC_PI;
	}

	inline double ang2rad(double ang) {
		return ang * 0.005555555555556 * SQI_PI;
	}

	// return the ang of std::asin
	inline double safetyAsin(double value) {
		if (value < -1.) {
			rad2ang(asin(-1.));
		}
		if (value > 1.) {
			rad2ang(asin(1.));
		}
		return rad2ang(asin(value));
	}

	// return the ang of std::acos
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
		Eigen::Vector3d v
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
			double w = 10., double h = 10.
		) {
			// SQI_VERBOSE_ONLY_COUT("");

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

	/* Parameterization cylinder
	* C(s, t) = cor + r * (cos(t) * u + sin(t) * v) + s * nor 
	* t = [0, 360], s = [-infty, +infty] */
	class Cylinder {
	public:
		Cylinder() {};
		~Cylinder() {};
		Cylinder(Eigen::Vector3d cor, Eigen::Vector3d nor, double r) {
			nor.normalize();
			cor_ = cor; nor_ = nor; 
			r_ = std::abs(r);
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
			Eigen::Vector3d plane_cor_cor = (plane_cor - cor_).normalized();
			double tu_component = plane_cor_cor.dot(u_);
			double tv_component = plane_cor_cor.dot(v_);
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
				SQI_VERBOSE_ONLY_WARNING("input seg num is invalid!");
				return;
			}
			// SQI_VERBOSE_ONLY_COUT("");

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

	/* Parameterization sphere
	* S(s, t) = cor + r * (cos(t) * u + sin(t) * v) * sin(s) + r * cos(s) * nor 
	* t = [0, 360], s = [0, 180] */
	class Sphere {
	public:
		Sphere() {};
		~Sphere() {};
		Sphere(Eigen::Vector3d cor, double r) {
			cor_ = cor; nor_ = Eigen::Vector3d(0, 0, 1);
			r_ = std::abs(r);
			u_ = get_perpendicular_normal(nor_);
			v_ = nor_.cross(u_).normalized();
		}
		Sphere(Eigen::Vector3d cor, Eigen::Vector3d nor, double r) {
			nor.normalize();
			cor_ = cor; nor_ = nor;
			r_ = std::abs(r);
			u_ = get_perpendicular_normal(nor);
			v_ = nor.cross(u_).normalized();
		}
		Eigen::Vector3d get_point(double s, double t) {
			double rad_s = ang2rad(s);
			double rad_t = ang2rad(t);
			return cor_ + r_ * (std::cos(rad_t) * u_ + std::sin(rad_t) * v_) * std::sin(rad_s) + r_ * std::cos(rad_s) * nor_;
		}
		void get_s_t(Eigen::Vector3d p, double& s, double& t) {
			s = safetyAcos((p - cor_).dot(nor_) / r_);
			Eigen::Vector3d plane_cor = p - r_ * std::cos(ang2rad(s)) * nor_;
			Eigen::Vector3d plane_cor_cor = (plane_cor - cor_).normalized();
			double tu_component = plane_cor_cor.dot(u_);
			double tv_component = plane_cor_cor.dot(v_);
			t = safetyAcos(tu_component);
			if (std::abs(std::sin(ang2rad(t)) - tv_component) > SQI_EPS) {
				t = 360 - t;
			}
		}

		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }
		Eigen::Vector3d& u() { return u_; }
		Eigen::Vector3d& v() { return v_; }
		double& r() { return r_; }

		void output_model(
			std::vector<Eigen::Vector3d>& points,
			std::vector<Eigen::Vector3i>& faces,
			int h_seg = 32, int r_seg = 32
		) {
			if (h_seg < 3) {
				SQI_VERBOSE_ONLY_WARNING("input h_seg num is invalid!");
				return;
			}
			if (r_seg < 3) {
				SQI_VERBOSE_ONLY_WARNING("input r_seg num is invalid!");
				return;
			}
			// SQI_VERBOSE_ONLY_COUT("");

			// init
			std::vector<Eigen::Vector3d>().swap(points);
			std::vector<Eigen::Vector3i>().swap(faces);
			points.reserve((h_seg - 2) * r_seg + 2);
			faces.reserve((h_seg - 2) * r_seg * 2);
			double h_angle = 180. / (h_seg - 1);
			double r_angle = 360. / r_seg;

			for (double phi = 0, phi_end = 180 + SQI_EPS; phi < phi_end; phi += h_angle) {
				if (phi < SQI_EPS) {
					points.push_back(get_point(phi, 0));
					// points.push_back(cor_ + r_ * Eigen::Vector3d(0, 0, 1));
				}
				else if (std::abs(phi - 180) < SQI_EPS) {
					points.push_back(get_point(phi, 0));
					// points.push_back(cor_ + r_ * Eigen::Vector3d(0, 0, -1));
				}
				else {
					for (double theta = 0, theta_end = 360 - SQI_EPS; theta < theta_end; theta += r_angle) {
						points.push_back(get_point(phi, theta));
						/*
						points.push_back(
							cor_ + Eigen::Vector3d(
								r_ * std::sin(ang2rad(phi)) * std::cos(ang2rad(theta)),
								r_ * std::sin(ang2rad(phi)) * std::sin(ang2rad(theta)),
								r_ * std::cos(ang2rad(phi))
							));*/
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
		Eigen::Vector3d nor_; // the unit vector of the sphere's axis
		Eigen::Vector3d u_, v_; // the coordinate axes of the plane which perpendicular to the cylinder's axis, used for parameterization representation
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
			// SQI_VERBOSE_ONLY_COUT("");

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

	class ParameterizationCircle {
		// C(t) = cor_ + r_ * (cos(t) * u_ + sin(t) * v)
	public:
		ParameterizationCircle() {
			t_lb_ = 0; t_ub_ = 360;
		};
		~ParameterizationCircle() {};
		ParameterizationCircle(Eigen::Vector3d cor, Eigen::Vector3d nor, double r) {
			nor.normalize();
			cor_ = cor; nor_ = nor; r_ = r;
			u_ = get_perpendicular_normal(nor_);
			v_ = nor_.cross(u_).normalized();
			t_lb_ = 0; t_ub_ = 360;
		}
		Eigen::Vector3d get_point(double t) {
			double rad_t = ang2rad(t);
			return cor_ + r_ * (std::cos(rad_t) * u_ + std::sin(rad_t) * v_);
		}

		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }
		Eigen::Vector3d& u() { return u_; }
		Eigen::Vector3d& v() { return v_; }
		double& r() { return r_; }
		double& t_lb() { return t_lb_; }
		double& t_ub() { return t_ub_; }

		void output_points(
			std::vector<Eigen::Vector3d>& points,
			double t_step = 0.5) {
			std::vector<Eigen::Vector3d>().swap(points);

			for (double t = t_lb_, t_end = t_ub_ - SQI_EPS; t < t_end; t += t_step) {
				points.push_back(get_point(t));
			}
			points.push_back(get_point(t_ub_));
		}
	private:
		Eigen::Vector3d cor_; // circle's center
		Eigen::Vector3d nor_; // the vector perpendicular to the circle
		Eigen::Vector3d u_, v_; // the coordinate axes of the plane which perpendicular to the circle's nor_, used for parameterization representation
		double r_; // circle's radius
		double t_lb_, t_ub_;
	};

	class ParameterizationCylindricPoint {
	public:
		ParameterizationCylindricPoint() {};
		~ParameterizationCylindricPoint() {};
		ParameterizationCylindricPoint(double s, double t) {
			s_ = s; t_ = t;
		}
		Eigen::Vector3d get_point(Cylinder& C) {
			return C.get_point(s_, t_);
		}
		double& s() { return s_; }
		double& t() { return t_; }

		void output_points(
			Cylinder& C,
			std::vector<Eigen::Vector3d>& points
		) {
			// SQI_VERBOSE_ONLY_COUT("");

			std::vector<Eigen::Vector3d>().swap(points);
			points.push_back(C.get_point(s_, t_));
		}
	private:
		double s_, t_;
	};

	class ParameterizationCylindricLine {
	public:
		ParameterizationCylindricLine() {
			s_lb_ = -SQI_INFTY; s_ub_ = SQI_INFTY;
		};
		~ParameterizationCylindricLine() {};
		ParameterizationCylindricLine(double t) {
			t_ = t;
			s_lb_ = -SQI_INFTY; s_ub_ = SQI_INFTY;
		}
		ParameterizationCylindricLine(double t, double s_lb, double s_ub) {
			t_ = t;
			s_lb_ = s_lb; s_ub_ = s_ub;
		}
		// divide the line into two parts according to s, modify self [s, ub], return [lb, s]
		void separate_s(ParameterizationCylindricLine& L, double s) {
			if (s < s_lb_ || s > s_ub_) {
				SQI_VERBOSE_ONLY_WARNING("input s out of L's bound");
				return;
			}
			L = *this;
			L.s_ub_ = s;
			this->s_lb_ = s;
		}
		double& t() { return t_; }
		double& s_lb() { return s_lb_; }
		double& s_ub() { return s_ub_; }

		void output_points(
			Cylinder& C,
			std::vector<Eigen::Vector3d>& points,
			double h = 20., int h_seg = 100
		) {
			// SQI_VERBOSE_ONLY_COUT("");

			std::vector<Eigen::Vector3d>().swap(points);

			double h_step = h / h_seg;
			if (h_step < SQI_EPS) {
				SQI_VERBOSE_ONLY_WARNING("h or h_seg is invalid!");
				return;
			}
			double center_s = 0;
			Eigen::Vector3d cor = C.get_point(center_s, t_);
			for (double cur_s = -0.5 * h, cur_s_end = 0.5 * h + SQI_EPS; cur_s < cur_s_end; cur_s += h_step) {
				if (center_s + cur_s >= s_lb_ && center_s + cur_s <= s_ub_) {
					points.push_back(cor + (center_s + cur_s) * C.nor());
				}
			}
		}
	private:
		double t_;
		double s_lb_, s_ub_;
	};

	class ParameterizationCylindricCurve {
		// a(t) * s^2 + b(t) * s + c(t) = 0
	public:
		ParameterizationCylindricCurve() {
			std::vector<double>(1, 0).swap(a_t_);
			std::vector<double>(3, 0).swap(b_t_);
			std::vector<double>(6, 0).swap(c_t_);
			s_lb_ = -SQI_INFTY; s_ub_ = SQI_INFTY;
			t_lb_ = 0; t_ub_ = 360;
			s_part_ = 0;
		}
		~ParameterizationCylindricCurve() {
			std::vector<double>().swap(a_t_);
			std::vector<double>().swap(b_t_);
			std::vector<double>().swap(c_t_);
		};
		ParameterizationCylindricCurve(
			std::vector<double>& a_t, std::vector<double>& b_t, std::vector<double>& c_t,
			double s_lb, double s_ub, double t_lb, double t_ub
		) {
			a_t_ = a_t; b_t_ = b_t; c_t_ = c_t;
			s_lb_ = s_lb; s_ub_ = s_ub; t_lb_ = t_lb; t_ub_ = t_ub;
			s_part_ = 0;
		}
		ParameterizationCylindricCurve(
			std::vector<double>& a_t, std::vector<double>& b_t, std::vector<double>& c_t,
			double s_lb, double s_ub, double t_lb, double t_ub,
			int s_part
		) {
			a_t_ = a_t; b_t_ = b_t; c_t_ = c_t;
			s_lb_ = s_lb; s_ub_ = s_ub; t_lb_ = t_lb; t_ub_ = t_ub;
			s_part_ = s_part;
		}
		ParameterizationCylindricCurve(const ParameterizationCylindricCurve& PC) {
			this->a_t_ = PC.a_t_;
			this->b_t_ = PC.b_t_;
			this->c_t_ = PC.c_t_;
			this->s_lb_ = PC.s_lb_; this->s_ub_ = PC.s_ub_; this->t_lb_ = PC.t_lb_; this->t_ub_ = PC.t_ub_;
			this->s_part_ = PC.s_part_;
		}
		ParameterizationCylindricCurve& operator =(const ParameterizationCylindricCurve& PC) {
			if (this != &PC) {
				this->a_t_ = PC.a_t_;
				this->b_t_ = PC.b_t_;
				this->c_t_ = PC.c_t_;
				this->s_lb_ = PC.s_lb_; this->s_ub_ = PC.s_ub_; this->t_lb_ = PC.t_lb_; this->t_ub_ = PC.t_ub_;
				this->s_part_ = PC.s_part_;
			}
			return *this;
		}
		bool operator ==(const ParameterizationCylindricCurve& PC) const {
			if (this == &PC) return true;
			if (std::equal(this->a_t_.begin(), this->a_t_.end(), PC.a_t_.begin(), PC.a_t_.end()) == false) return false;
			if (std::equal(this->b_t_.begin(), this->b_t_.end(), PC.b_t_.begin(), PC.b_t_.end()) == false) return false;
			if (std::equal(this->c_t_.begin(), this->c_t_.end(), PC.c_t_.begin(), PC.c_t_.end()) == false) return false;
			if (this->s_lb_ != PC.s_lb_ || this->s_ub_ != PC.s_ub_ ||
				this->t_lb_ != PC.t_lb_ || this->t_ub_ != PC.t_ub_ ||
				this->s_part_ != PC.s_part_) return false;
			return true;
		}
		/* Divide this curve into 2 parts according to s, modify self 1, output -1. */
		void separate_s(ParameterizationCylindricCurve& PC) {
			PC = *this;
			this->s_part_ = 1;
			PC.s_part_ = -1;
		}
		/* Divide this curve into 2 parts according to t, modify self [t, ub], output [lb, t] */
		int separate_t(ParameterizationCylindricCurve& PC, double t) {
			if (t < t_lb_ || t > t_ub_) {
				SQI_VERBOSE_ONLY_WARNING("input t out of PC's bound");
				return 0;
			}
			PC = *this;
			PC.t_ub_ = t;
			this->t_lb_ = t;
			return 1;
		}
		/* Divide this curve into two parts according to t1-t2 in the total parameterization region, modify self [t1, t2], output [t2, t1+360] */
		void separate_t(ParameterizationCylindricCurve& PC, double t1, double t2) {
			PC = *this;
			if (t1 < t2) { // sort to let t1 > t2
				double tmp_t = t1;
				t1 = t2; t2 = tmp_t;
			}
			this->t_lb_ = t2; this->t_ub_ = t1;
			PC.t_lb_ = t1; PC.t_ub_ = t2 + 360;
		}
		int get_s(double t, std::vector<double>& s) {
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
					// SQI_VERBOSE_ONLY_COUT("cannot get the value of s");
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
		bool is_valid() {
			if (a_t_.size() != 1 || b_t_.size() != 3 || c_t_.size() != 6)	return false;
			return true;
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
		void output_points(
			Cylinder& C,
			std::vector<Eigen::Vector3d>& points,
			double t_step = 0.5
		) {
			// SQI_VERBOSE_ONLY_COUT("");

			std::vector<Eigen::Vector3d>().swap(points);

			for (double cur_t = t_lb_; cur_t < t_ub_ + t_step + SQI_EPS; cur_t += t_step) {
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

	// --------------assessment of the intersection primitives' intersections--------------

	/*	Get the results of the intersection between a parameterization Line L1 and a parameterization Curve PC1.
	*	This function will involve clipping the original input, and outputting other clipped lines and curves through std::vector.
	*	Return The number of newly generated lines and curves. 
	*/
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

	/*	Get the results of the intersection between a parameterization Curve PC1 and a parameterization Curve PC2.
	*	This function will involve clipping the original input, and outputting other clipped curves through std::vector.
	*	Return the number of newly generated curves. 
	*/
	int get_intersections(
		ParameterizationCylindricCurve& PC1, ParameterizationCylindricCurve& PC2,
		std::vector<ParameterizationCylindricCurve>& sub_PC1s,
		std::vector<ParameterizationCylindricCurve>& sub_PC2s,
		double t_step = 1.
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

	/*	Get the results of the intersection between parameterization Lines L1s and parameterization Curves PC1s, and the results are stored in the original input std::vector.
	*	Return the number of total lines and curves. 
	*/
	int get_intersections(
		std::vector<ParameterizationCylindricLine>& L1s,
		std::vector<ParameterizationCylindricCurve>& PC1s
	) {
		SQI_VERBOSE_ONLY_COUT("");

		for (int i = 0; i < L1s.size(); ++i) { // note: the vector L1s's size is dymanic, so cannot use i_end = L1s.size.
			for (int j = 0, j_end = PC1s.size(); j <j_end; ++j) { // note: the vector PC1s's size is dymanic, but the newly added curves do not need to be judged again.
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

	/*	Get the results of the intersection between parameterization Curves PC1s and parameterization Curves PC2s, and the results are stored in the original input std::vector.
	*   Return the number of total curves. 
	*/
	int get_intersections(
		std::vector<ParameterizationCylindricCurve>& PC1s,
		std::vector<ParameterizationCylindricCurve>& PC2s,
		double t_step = 1.
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

	// ----------------------------assessment of the intersections----------------------------

	/* Get the intersections between Line L1 and Cylinder C1.
	/note: using Cylinder C1 as the parameterizatoin surface.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Line L1, Cylinder C1,
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

	// ****** using plane as the parameterization surface ******

	// ****** using cylinder as the parameterization surface ******

	/* Get the intersections between Plane P1 and Cylinder C1.
	/note: using Cylinder C1 as the parameterizatoin surface.
	*  Return the number of intersections.
	*/
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
		double P1_cor_to_C1_cor_dot_P1_nor = P1_cor_to_C1_cor.dot(P1.nor());
		double P1_cor_to_C1_dis = std::abs(P1_cor_to_C1_cor_dot_P1_nor);

		if (P1.nor().dot(C1.nor()) < SQI_EPS) { // plane is parallel to cylinder			
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
					C1.r() * C1.u().dot(P1.nor()),
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

	/* Get the intersections between Cylinder C1 and Cylinder C2.
	/note: using Cylinder C2 as the parameterizatoin surface.
	*  Return the number of intersections.
	*/
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

					// cut curves into 4 parts
					ParameterizationCylindricCurve PC_C2_2;
					PC_C2.separate_t(PC_C2_2, t1, t2);
					ParameterizationCylindricCurve PC_C2_11;
					PC_C2.separate_s(PC_C2_11);
					ParameterizationCylindricCurve PC_C2_21;
					PC_C2_2.separate_s(PC_C2_21);
					
					curves.push_back(PC_C2);
					curves.push_back(PC_C2_2);
					curves.push_back(PC_C2_11);
					curves.push_back(PC_C2_21);
				}
				else if (L1_status == 2 || L2_status == 2) { // two-point intersection + do not intersect: result 1 parameterization curve
					std::vector<ParameterizationCylindricPoint>* intersection_points;
					if (L1_status == 2) intersection_points = &L1_C2_intersect_points;
					if (L2_status == 2) intersection_points = &L2_C2_intersect_points;

					double t1 = (*intersection_points)[0].t();
					double t2 = (*intersection_points)[1].t();

					if (t1 < t2) { // sort to let t1 > t2
						double tmp_t = t1;
						t1 = t2; t2 = tmp_t;
					}

					if (axes_dis > C1.r() + SQI_EPS) { // circle center angle < 180
						if (t1 - t2 > 180) {
							PC_C2.t_lb() = t1; PC_C2.t_ub() = t2 + 360;
						}
						else {
							PC_C2.t_lb() = t2; PC_C2.t_ub() = t1;
						}
					}
					else if (axes_dis < C1.r() - SQI_EPS) { // circle center angle > 180
						if (t1 - t2 > 180) {
							PC_C2.t_lb() = t2; PC_C2.t_ub() = t1;
						}
						else {
							PC_C2.t_lb() = t1; PC_C2.t_ub() = t2 + 360;
						}
					}
					double center_t = 0.5 * (PC_C2.t_lb() + PC_C2.t_ub());

					// cut curves into 4 parts
					ParameterizationCylindricCurve PC_C2_s, PC_C2_st, PC_C2_t;
					PC_C2.separate_s(PC_C2_s);
					PC_C2.separate_t(PC_C2_t, center_t);
					PC_C2_s.separate_t(PC_C2_st, center_t);

					curves.push_back(PC_C2);
					curves.push_back(PC_C2_s);
					curves.push_back(PC_C2_st);
					curves.push_back(PC_C2_t);
				}
				else {
					SQI_VERBOSE_ONLY_COUT("may error?");
				}
			}
		}

		return points.size() + lines.size() + curves.size();
	}

	/* Get the intersections between Sphere S1 and Cylinder C1.
	/note: using Cylinder C1 as the parameterizatoin surface.
	*  Return the number of intersections.
	*/
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
			SQI_VERBOSE_ONLY_COUT("test" << " " << s << " " << t << " " << rot_angle);
			double t1 = t - rot_angle;
			double t2 = t + rot_angle;

			ParameterizationCylindricCurve PC_tmp;
			PC.separate_t(PC_tmp, t1, t2);

			ParameterizationCylindricCurve PC2;
			PC.separate_s(PC2);

			curves.push_back(PC);
			curves.push_back(PC2);
		}

		return points.size() + curves.size();
	}

	// ****** using sphere as the parameterization surface ******

	/* Get the intersections between Plane P1 and Sphere S1.
	/note: using the vector from S2's center to a point on P1 as the parameterizatoin circle's normal.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Plane& P1, Sphere& S1,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles
	) {
		SQI_VERBOSE_ONLY_COUT("");

		// init 
		std::vector<Point>().swap(points);
		std::vector<ParameterizationCircle>().swap(circles);

		Eigen::Vector3d S1_cor_to_P1_cor = P1.cor() - S1.cor();
		double P1_S1_cor_dot_P1_nor = S1_cor_to_P1_cor.dot(P1.nor());
		double P1_S1_cor_dis = std::abs(P1_S1_cor_dot_P1_nor);

		if (P1_S1_cor_dis > S1.r() + SQI_EPS) { // not intersect
			SQI_VERBOSE_ONLY_COUT("not intersect");
		}
		else {
			if (P1_S1_cor_dis > S1.r() - SQI_EPS) { // tangent
				SQI_VERBOSE_ONLY_COUT("tangent");

				points.push_back(Point(S1.cor() + P1_S1_cor_dot_P1_nor * P1.nor()));
			}
			else { // intersect
				SQI_VERBOSE_ONLY_COUT("intersect");

				Eigen::Vector3d circle_center = S1.cor() + P1_S1_cor_dot_P1_nor * P1.nor();
				Eigen::Vector3d circle_nor = P1.nor();
				if (P1_S1_cor_dot_P1_nor < 0) circle_nor = -circle_nor;
				double circle_r = std::sqrt(S1.r() * S1.r() - P1_S1_cor_dis * P1_S1_cor_dis);

				circles.push_back(ParameterizationCircle(circle_center, -P1.nor(), circle_r));
			}
		}

		return points.size() + circles.size();
	}

	/* Get the intersections between Cylinder C1 and Sphere S1.
	/note: using Cylinder C1 as the parameterizatoin surface.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Cylinder& C1, Sphere& S1,
		std::vector<ParameterizationCylindricPoint>& points,
		std::vector<ParameterizationCylindricCurve>& curves
	) {
		return get_intersections(S1, C1, points, curves);
	}

	/* Get the intersections between Sphere S1 and Sphere S2.
	/note: using the vector from S2's center to S1's center as the parameterizatoin circle's normal.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Sphere& S1, Sphere& S2,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles
	) {
		SQI_VERBOSE_ONLY_COUT("");

		// init
		std::vector<Point>().swap(points);
		std::vector<ParameterizationCircle>().swap(circles);

		Eigen::Vector3d S2_cor_to_S1_cor = S1.cor() - S2.cor();
		double centers_dis = S2_cor_to_S1_cor.norm();

		if (centers_dis < SQI_EPS) { // overlap: result nothing
			SQI_VERBOSE_ONLY_COUT("overlap");
		}
		else {
			Eigen::Vector3d S2_cor_to_S1_cor_unit = S2_cor_to_S1_cor / centers_dis;

			if (centers_dis > S1.r() + S2.r() + SQI_EPS) { // not intersect
				SQI_VERBOSE_ONLY_COUT("not intersect");
			}
			else if (centers_dis > S1.r() + S2.r() - SQI_EPS) { // tangent: result a point
				SQI_VERBOSE_ONLY_COUT("tangent");

				points.push_back(Point(S2.cor() + S2.r() * S2_cor_to_S1_cor_unit));
			}
			else { // intersect: result a circle
				SQI_VERBOSE_ONLY_COUT("intersect");

				Eigen::Vector3d circle_center = S2.cor() + 0.5 * S2_cor_to_S1_cor;
				double circle_r = std::sqrt(S2.r() * S2.r() - 0.25 * centers_dis * centers_dis);

				circles.push_back(
					ParameterizationCircle(circle_center, S2_cor_to_S1_cor_unit, circle_r));
			}
		}

		return points.size() + circles.size();
	}

	// ----------------------------use for testing----------------------------

	/* Get the intersections between Cylinder C1 and other primitives.
	/note: using Cylinder C1 as the parameterization surface.
	Return the number of all intersections.
	*/
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

	// ----------------------------output for debug----------------------------

	/* Output the primitives as the meshes and export to an OBJ file. 
	*/
	void write_result_mesh(
		std::string output_file_path,
		std::vector<Plane>& planes,
		std::vector<Cylinder>& cylinders,
		std::vector<Sphere>& spheres
	) {
		// SQI_VERBOSE_ONLY_COUT("");

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

	/* Output the intersection primitives as points and export to an OBJ file. 
	/note: using Cylinder C as the parameterization surface.
	*/
	void write_cylinder_result_points(
		std::string output_file_path,
		Cylinder& C,
		std::vector<ParameterizationCylindricPoint>& points,
		std::vector<ParameterizationCylindricLine>& lines,
		std::vector<ParameterizationCylindricCurve>& curves
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		int primitive_num = points.size() + lines.size() + curves.size();
		if (rand_color_bar.size() < primitive_num) {
			create_rand_color_bar(primitive_num);
		}

		std::ofstream out(output_file_path);
		int color_num = 0;
		
		for (ParameterizationCylindricPoint P : points) {
			Eigen::Vector3d color = rand_color_bar[color_num++];

			out << "v" << " " << P.get_point(C).transpose() << " " << color.transpose() << std::endl;
		}

		for (ParameterizationCylindricLine L : lines) {
			Eigen::Vector3d color = rand_color_bar[color_num++];

			std::vector<Eigen::Vector3d> output_points;
			L.output_points(C, output_points);

			for (Eigen::Vector3d p : output_points) {
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
		}

		for (ParameterizationCylindricCurve PC : curves) {
			Eigen::Vector3d color = rand_color_bar[color_num++];

			std::vector<Eigen::Vector3d> output_points;
			PC.output_points(C, output_points);

			for (Eigen::Vector3d p : output_points) {
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
		}
		out.close();
	}

	/* Output the intersection primitives as points and export to an OBJ file.
	*/
	void write_sphere_result_points(
		std::string output_file_path,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles
	) {
		SQI_VERBOSE_ONLY_COUT("");

		std::ofstream out(output_file_path);
		for (Point p : points) {
			out << "v" << " " << p.cor().transpose() << std::endl;
		}
		for (ParameterizationCircle c : circles) {
			std::vector<Eigen::Vector3d> ops;
			c.output_points(ops);

			for (Eigen::Vector3d p : ops) {
				out << "v" << " " << p.transpose() << std::endl;
			}
		}
		out.close();
	}
}

#endif
