#ifndef SIMPLE_QUADRICS_INTERSECTION_H_
#define SIMPLE_QUADRICS_INTERSECTION_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <functional>
#include <unordered_map>
#include <fstream>
#include <stdio.h>
// Eigen
#include <Eigen/Dense>

#define USE_FOR_OFFSET_MESH_GENERATION

#define SIMPLE_QUADRICS_INTERSECTION_VERBOSE_
#ifdef SIMPLE_QUADRICS_INTERSECTION_VERBOSE_
#define SQI_VERBOSE_ONLY_TITLE(x) std::cout << "\033[32m" << "[" << __FUNCTION__ << "]" << "\033[0m" << " " << x << std::endl // [green] white cout
#define SQI_VERBOSE_ONLY_COUT(x) std::cout << "\033[33m" << "[" << __FUNCTION__ << "]" << "\033[0m" << " " << x << std::endl // [yellow] white cout
#define SQI_VERBOSE_ONLY_WARNING(x) std::cout << "\033[33m" << "[" << __FILE__ << " " << __LINE__ << "]" << "\033[0m" << " " << "\033[31m" << x << "\033[0m" << std::endl // [yellow] red cout
#define SQI_VERBOSE_ONLY_TEST(x) std::cout << "\033[33m" << "[" << __FUNCTION__ << "]" << "\033[0m" << " " << "\033[32m" << x << "\033[0m" << std::endl // [yellow] green cout
#else
#define SQI_VERBOSE_ONLY_TITLE(x)
#define SQI_VERBOSE_ONLY_COUT(x)
#define SQI_VERBOSE_ONLY_WARNING(x)
#define SQI_VERBOSE_ONLY_TEST(x)
#endif

#define SQI_EPS	1e-12
#define SQI_INFTY 1e12
#define SQI_PI 3.141592653589793
#define SQI_REC_PI 0.318309886183791

namespace QuadricsIntersection {
	// -----------------------other functions

	struct pair_hash {
		template <class T1, class T2>
		size_t operator () (std::pair<T1, T2> const& pair
			) const {
			size_t h1 = std::hash<T1>()(pair.first);
			size_t h2 = std::hash<T2>()(pair.second);
			return h1 ^ h2;
		}
	};

	void create_rand_color_bar(int num);

	inline double rad2ang(double rad) {
		return rad * 180 * SQI_REC_PI;
	}

	inline double ang2rad(double ang) {
		return ang * 0.005555555555556 * SQI_PI;
	}

	// return the ang of std::asin
	inline double safetyAsin(double value) {
		if (value < -1.) return -90.;
		if (value > 1.) return 90.;
		return rad2ang(asin(value));
	}

	// return the ang of std::acos
	inline double safetyAcos(double value) {
		if (value < -1.) return 360.;
		if (value > 1.) return 0.;
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

	// -----------------------the quadric primitives

	class Plane {
	public:
		Plane() {};
		~Plane() {};
		Plane(Eigen::Vector3d cor, Eigen::Vector3d nor) {
			nor.normalize();
			cor_ = cor; nor_ = nor;
		}

		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }

		void output_model(
			std::vector<Eigen::Vector3d>& points,
			std::vector<Eigen::Vector3i>& faces,
			double w = 10., double h = 10.);
	private:
		Eigen::Vector3d cor_; // a point on the plane
		Eigen::Vector3d nor_; // the unit normal of the plane
	};

	/* Parameterization cylinder
	* C(s, t) = cor + r * (cos(t) * u + sin(t) * v) + s * nor
	* t = [-180, 180], s = [-infty, +infty] */
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
		Cylinder(Eigen::Vector3d cor, Eigen::Vector3d nor, Eigen::Vector3d u, 
			double r, double t_lb, double t_ub, double s_lb, double s_ub) {
			nor.normalize(); u.normalize();
			cor_ = cor; nor_ = nor;
			r_ = std::abs(r);
			u_ = u;
			v_ = nor.cross(u_).normalized();
			t_lb_ = t_lb; t_ub_ = t_ub;
			s_lb_ = s_lb; s_ub_ = s_ub;
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
			// if (tv_component < 0) t = 360 - t; // t = [0, 360]
			if (tv_component < 0) t = -t; // t = [-180, 180]
		}
		Cylinder& operator =(const Cylinder& C) {
			if (this != &C) {
				this->cor_ = C.cor_;
				this->nor_ = C.nor_;
				this->r_ = C.r_;
				this->u_ = C.u_; this->v_ = C.v_;
				this->t_lb_ = C.t_lb_; this->t_ub_ = C.t_ub_;
				this->s_lb_ = C.s_lb_; this->s_ub_ = C.s_ub_;
			}
			return *this;
		}
		bool operator ==(const Cylinder& C) {
			if (this == &C) return true;
			if ((this->cor_ - C.cor_).squaredNorm() < SQI_EPS &&
				1 - this->nor_.dot(C.nor_) < SQI_EPS) return true;
			return false;
		}
		bool compare(Cylinder& C) {
			if (this == &C) return true;
			if ((this->cor_ - C.cor_).squaredNorm() < SQI_EPS &&
				1 - this->nor_.dot(C.nor_) < SQI_EPS) return true;
			return false;
		}

		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }
		double& r() { return r_; }
		Eigen::Vector3d& u() { return u_; }
		Eigen::Vector3d& v() { return v_; }
		double& t_lb() { return t_lb_; }
		double& t_ub() { return t_ub_; }
		double& s_lb() { return s_lb_; }
		double& s_ub() { return s_ub_; }

		void output_model(
			std::vector<Eigen::Vector3d>& points,
			std::vector<Eigen::Vector3i>& faces,
			int seg = 32, double h = 20.);
	private:
		Eigen::Vector3d cor_; // a point at the cylinder's axis
		Eigen::Vector3d nor_; // the unit vector of the cylinder's axis
		double r_; // the radius of the cylinder
		Eigen::Vector3d u_, v_; // the coordinate axes of the plane which perpendicular to the cylinder's axis, used for parameterization representation
		double t_lb_ = -180, t_ub_ = 180;
		double s_lb_ = -SQI_INFTY, s_ub_ = SQI_INFTY;
	};

	/* Parameterization sphere
	* S(s, t) = cor + r * (cos(t) * u + sin(t) * v) * sin(s) + r * cos(s) * nor
	* t = [-180, 180], s = [0, 180] */
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
		Sphere(Eigen::Vector3d cor, Eigen::Vector3d nor, double r, double s_lb, double s_ub) {
			nor.normalize();
			cor_ = cor; nor_ = nor;
			r_ = std::abs(r);
			u_ = get_perpendicular_normal(nor);
			v_ = nor.cross(u_).normalized();
			s_lb_ = s_lb; s_ub_ = s_ub;
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
			if (tv_component < 0) t = -t;
		}

		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }
		Eigen::Vector3d& u() { return u_; }
		Eigen::Vector3d& v() { return v_; }
		double& r() { return r_; }
		double& s_lb() { return s_lb_; }
		double& s_ub() { return s_ub_; }

		void output_model(
			std::vector<Eigen::Vector3d>& points,
			std::vector<Eigen::Vector3i>& faces,
			int h_seg = 32, int r_seg = 32);
	private:
		Eigen::Vector3d cor_; // the center of the sphere
		Eigen::Vector3d nor_; // the unit vector of the sphere's axis
		Eigen::Vector3d u_, v_; // the coordinate axes of the plane which perpendicular to the cylinder's axis, used for parameterization representation
		double r_; // the radius of the sphere
		double s_lb_ = 0, s_ub_ = 180;
	};

	// -----------------------the intersection primitives

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
		// cor + s * nor
	public:
		Line() {};
		~Line() {};
		Line(Eigen::Vector3d cor, Eigen::Vector3d nor) {
			nor.normalize();
			cor_ = cor; nor_ = nor;
			s_lb_ = -SQI_INFTY; s_ub_ = SQI_INFTY;
		}
		Eigen::Vector3d get_point(double s) {
			return cor_ + s * nor_;
		}
		double get_s(Eigen::Vector3d p) {
			return (p - cor_).dot(nor_);
		}
		bool is_s_valid(double s) {
			if (s > s_lb_ + SQI_EPS && s < s_ub_ - SQI_EPS) return true;
			return false;
		}
		/* separate the line according to s, modify self [s, s_ub], output [s_lb, s].*/
		int separate_s(Line& L, double s) {
			if (s < s_lb_ + SQI_EPS || s > s_ub_ - SQI_EPS) return 0;
			L = *this;
			this->s_lb_ = s; L.s_ub_ = s;
			return 1;
		}

		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }
		double& s_lb() { return s_lb_; }
		double& s_ub() { return s_ub_; }

		void output_model(
			std::vector<Eigen::Vector3d>& points,
			std::vector<Eigen::Vector2i>& lines,
			double l = 20.);
		void output_points(
			std::vector<Eigen::Vector3d>& points,
			double l = 50., int seg = 200);
	private:
		Eigen::Vector3d cor_; // a point at the line
		Eigen::Vector3d nor_; // the vector of the line
		double s_lb_, s_ub_;
	};

	/* Parameterization circle
	* c(t) = cor_ + r_ * (cos(t) * u_ + sin(t) * v)
	* t = [-180, 180] */
	class ParameterizationCircle {
	public:
		ParameterizationCircle() {};
		~ParameterizationCircle() {};
		ParameterizationCircle(Eigen::Vector3d cor, Eigen::Vector3d nor, double r) {
			nor.normalize();
			cor_ = cor; nor_ = nor; r_ = r;
			u_ = get_perpendicular_normal(nor_);
			v_ = nor_.cross(u_).normalized();
		}
		Eigen::Vector3d get_point(double t) {
			double rad_t = ang2rad(t);
			return cor_ + r_ * (std::cos(rad_t) * u_ + std::sin(rad_t) * v_);
		}
		double get_t(Eigen::Vector3d p) {
			p -= (p - cor_).dot(nor_) * nor_; // project to the circle's plane
			Eigen::Vector3d center_to_p = (p - cor_).normalized();
			double tu_component = center_to_p.dot(u_);
			double tv_component = center_to_p.dot(v_);
			double t = safetyAcos(tu_component);
			if (tv_component < 0) t = -t;
			return t;
		}
		bool is_t_valid(double t) {
			if (t > t_lb_ + SQI_EPS && t < t_ub_ - SQI_EPS) return true;
			return false;
		}
		bool is_on(Eigen::Vector3d p) {
			Eigen::Vector3d p_cor = p - cor_;
			if (std::abs(p_cor.squaredNorm() - r_ * r_) < SQI_EPS &&
				std::abs(p_cor.normalized().dot(nor_)) < SQI_EPS) return true;
			return false;
		}
		/* separate the circle according to t, modify self [t, t_ub], output [t_lb, t].*/
		int separate_t(ParameterizationCircle& C, double t) {
			if (t < t_lb_ + SQI_EPS || t > t_ub_ - SQI_EPS) return 0;
			C = *this;
			this->t_lb_ = t; C.t_ub_ = t;
			return 1;
		}
		/* separate the circle according to some t */
		void separate_t(std::vector<ParameterizationCircle>& Cs, std::vector<double>& ts) {
			std::vector<ParameterizationCircle>().swap(Cs);

			std::sort(ts.begin(), ts.end());
			for (double t : ts) {
				ParameterizationCircle sub_C;
				if (this->separate_t(sub_C, t) == 1) Cs.push_back(sub_C);
			}
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
			double t_step = 0.5);
	private:
		Eigen::Vector3d cor_; // circle's center
		Eigen::Vector3d nor_; // the vector perpendicular to the circle
		Eigen::Vector3d u_, v_; // the coordinate axes of the plane which perpendicular to the circle's nor_, used for parameterization representation
		double r_; // circle's radius
		double t_lb_ = -180, t_ub_ = 180;
	};

	/* Parameterization cylinder (of a cylinder)
	* a(t) * s^2 + b(t) * s + c(t) = 0
	* t = [-180, 180], s = [-INFTY, INFTY] (also determine by s_part) */
	class ParameterizationCylindricCurve {
	public:
		ParameterizationCylindricCurve() {
			std::vector<double>(1, 0).swap(a_t_);
			std::vector<double>(3, 0).swap(b_t_);
			std::vector<double>(6, 0).swap(c_t_);
			// t_lb_ = 0; t_ub_ = 360;
			s_part_ = 0;
			C_ = Cylinder();
		}
		~ParameterizationCylindricCurve() {
			std::vector<double>().swap(a_t_);
			std::vector<double>().swap(b_t_);
			std::vector<double>().swap(c_t_);
		};
		ParameterizationCylindricCurve(
			std::vector<double>& a_t, std::vector<double>& b_t, std::vector<double>& c_t,
			Cylinder& C
		) {
			a_t_ = a_t; b_t_ = b_t; c_t_ = c_t;
			C_ = C;
		}
		ParameterizationCylindricCurve(
			std::vector<double>& a_t, std::vector<double>& b_t, std::vector<double>& c_t,
			double t_lb, double t_ub,
			Cylinder& C
		) {
			a_t_ = a_t; b_t_ = b_t; c_t_ = c_t;
			t_lb_ = t_lb; t_ub_ = t_ub;
			C_ = C;
		}
		ParameterizationCylindricCurve(
			std::vector<double>& a_t, std::vector<double>& b_t, std::vector<double>& c_t,
			double s_lb, double s_ub, double t_lb, double t_ub,
			Cylinder& C
		) {
			a_t_ = a_t; b_t_ = b_t; c_t_ = c_t;
			s_lb_ = s_lb; s_ub_ = s_ub; t_lb_ = t_lb; t_ub_ = t_ub;
			s_part_ = 0;
			C_ = C;
		}
		ParameterizationCylindricCurve(
			std::vector<double>& a_t, std::vector<double>& b_t, std::vector<double>& c_t,
			double s_lb, double s_ub, double t_lb, double t_ub,
			int s_part,
			Cylinder& C
		) {
			a_t_ = a_t; b_t_ = b_t; c_t_ = c_t;
			s_lb_ = s_lb; s_ub_ = s_ub; t_lb_ = t_lb; t_ub_ = t_ub;
			s_part_ = s_part;
			C_ = C;
		}
		ParameterizationCylindricCurve(const ParameterizationCylindricCurve& PC) {
			this->a_t_ = PC.a_t_;
			this->b_t_ = PC.b_t_;
			this->c_t_ = PC.c_t_;
			this->s_lb_ = PC.s_lb_; this->s_ub_ = PC.s_ub_; this->t_lb_ = PC.t_lb_; this->t_ub_ = PC.t_ub_;
			this->s_part_ = PC.s_part_;
			this->C_ = PC.C_;
		}
		ParameterizationCylindricCurve& operator =(const ParameterizationCylindricCurve& PC) {
			if (this != &PC) {
				this->a_t_ = PC.a_t_;
				this->b_t_ = PC.b_t_;
				this->c_t_ = PC.c_t_;
				this->s_lb_ = PC.s_lb_; this->s_ub_ = PC.s_ub_; this->t_lb_ = PC.t_lb_; this->t_ub_ = PC.t_ub_;
				this->s_part_ = PC.s_part_;
				this->C_ = PC.C_;
			}
			return *this;
		}
		/* Divide this curve into 2 parts according to s, modify self 1, output -1. */
		void separate_s(ParameterizationCylindricCurve& PC) {
			PC = *this;
			this->s_part_ = 1;
			PC.s_part_ = -1;
		}
		/* Divide this curve into 2 parts according to t, modify self [t, ub], output [lb, t] */
		int separate_t(ParameterizationCylindricCurve& PC, double t) {
			if (t < t_lb_ + SQI_EPS || t > t_ub_ - SQI_EPS) return 0;
			PC = *this;
			PC.t_ub_ = t;
			this->t_lb_ = t;
			return 1;
		}
		Eigen::Vector3d get_point(double s, double t) {
			return C_.get_point(s, t);
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
		bool is_t_valid(double t) {
			if (t > t_lb_ - SQI_EPS && t < t_ub_ + SQI_EPS) return true;
			return false;
		}
		bool is_t_strict_valid(double t) {
			if (t > t_lb_ + SQI_EPS && t < t_ub_ - SQI_EPS) return true;
			return false;
		}
		bool is_s_t_valid(double s, double t) {
			if (is_t_valid(t) == false) return false;
			std::vector<double> ss;
			get_s(t, ss);
			for (int i = 0, i_end = ss.size(); i < i_end; ++i) {
				if (std::abs(ss[i] - s) < SQI_EPS) return true;
			}
			return false;
		}
		bool compare(ParameterizationCylindricCurve& PC) {
			if (this == &PC) return true;

			double overlap_t_lb = std::max(this->t_lb(), PC.t_lb());
			double overlap_t_ub = std::min(this->t_ub(), PC.t_ub());
			if (overlap_t_ub > overlap_t_lb - SQI_EPS) return false;

			for (int i = 0; i < 3; ++i) {
				double r = double(rand()) / RAND_MAX;
				double t = overlap_t_lb + (overlap_t_ub - overlap_t_lb) * r;
				std::vector<double> ss1, ss2;
				this->get_s(t, ss1); PC.get_s(t, ss2);
				if (ss1.size() != 1 || ss2.size() != 1) return false;
				if (std::abs(ss1[0] - ss2[0]) > SQI_EPS) return true;
			}

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
		Cylinder& C() { return C_; }

		void verbose() {
			SQI_VERBOSE_ONLY_COUT("a_t_:" << " " << a_t_[0]);
			SQI_VERBOSE_ONLY_COUT("b_t_:" << " " << b_t_[0] << " " << b_t_[1] << " " << b_t_[2]);
			SQI_VERBOSE_ONLY_COUT("c_t_:" << " " << c_t_[0] << " " << c_t_[1] << " " << c_t_[2] << " " << c_t_[3] << " " << c_t_[4] << " " << c_t_[5]);
			SQI_VERBOSE_ONLY_COUT("s_bound:" << " " << s_lb_ << "--" << s_ub_ << " " << "s_part:" << " " << s_part_);
			SQI_VERBOSE_ONLY_COUT("t_bound:" << " " << t_lb_ << "--" << t_ub_);
		}
		void output_points(
			std::vector<Eigen::Vector3d>& points,
			double t_step = 0.5);
	private:
		std::vector<double> a_t_; // = [0]
		std::vector<double> b_t_; // = [0] * cos(t) + [1] * sin(t) + [2]
		std::vector<double> c_t_; // = [0] * cos(t)^2 + [1] * sin(t)^2 + [2] * sin(t) * cos(t) + [3] * cos(t) + [4] * sin(t) + [5]
		double s_lb_ = -SQI_INFTY, s_ub_ = SQI_INFTY, t_lb_ = -180, t_ub_ = 180;
		int s_part_; // decide which part of the curve is need, 0: all, 1: upper, -1: lower
		Cylinder C_;
	};

	// -----------------------assessment of the intersections between primitives

	// -----------------------about plane

	/* Get the intersecting line between Plane P1 and Plane P2.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Plane& P1, Plane& P2,
		Line& L,
		double limit_angle = 0.);
	int get_intersections(
		Plane& P1, Plane& P2,
		std::vector<Line>& Ls,
		double limit_angle = 0.);

	/* Get the intersecting line between Cylinder C1 and Plane P1.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Cylinder& C1, Plane& P1,
		std::vector<Line>& lines,
		std::vector<ParameterizationCylindricCurve>& curves,
		double limit_angle = 0.);

	/* Get the intersecting line between Sphere S1 and Plane P1.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Sphere& S1, Plane& P1,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles);

	// -----------------------about cylinder

	/* Get the intersections between Line L1 and Cylinder C1.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Line& L1, Cylinder& C1,
		std::vector<Point>& points);

	/* Get the intersections between Plane P1 and Cylinder C1.
	*  Return the number of intersections.
	/note: something not done yet!
	*/
	int get_intersections(
		Plane& P1, Cylinder& C1,
		std::vector<Line>& lines,
		std::vector<ParameterizationCylindricCurve>& curves,
		double limit_angle = 0.);

	/* Get the intersections between Cylinder C1 and Cylinder C2.
	/note: the result curve are on cylinder C2's parameterization surface
	*  Return the number of intersections.
	*/
	int get_intersections(
		Cylinder& C1, Cylinder& C2,
		std::vector<Point>& points,
		std::vector<Line>& lines,
		std::vector<ParameterizationCylindricCurve>& curves);

	/* Get the intersections between Sphere S1 and Cylinder C1.
	*  Return the number of intersections.
	/note: something not done yet!
	*/
	int get_intersections(
		Sphere& S1, Cylinder& C1,
		std::vector<Point>& points,
		std::vector<ParameterizationCylindricCurve>& curves,
		double limit_angle = 0.);

	// -----------------------about sphere

	/* Get the intersections between Plane P1 and Sphere S1.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Plane& P1, Sphere& S1,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles);

	/* Get the intersections between Cylinder C1 and Sphere S1.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Cylinder& C1, Sphere& S1,
		std::vector<Point>& points,
		std::vector<ParameterizationCylindricCurve>& curves);

	/* Get the intersections between Sphere S1 and Sphere S2.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Sphere& S1, Sphere& S2,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles);

	// -----------------------assessment of the intersections between primitives' intersections

	/* Get the intersections between Point P1 and Line L1.
	*  This function will involve clipping the original input, and outputting other clipped primitives through std::vector.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Point& P1, Line& L1,
		std::vector<Line>& sub_L1s);
	int get_intersections(
		std::vector<Point>& P1s,
		std::vector<Line>& L1s);

	/* Get the intersections between Point P1 and Circle c1.
	*  This function will involve clipping the original input, and outputting other clipped primitives through std::vector.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Point& P1, ParameterizationCircle& c1,
		std::vector<ParameterizationCircle>& sub_c1s);
	int get_intersections(
		std::vector<Point>& P1s,
		std::vector<ParameterizationCircle>& c1s);

	/* Get the intersections between Point P1 and CylindricCurve PC1.
	*  This function will involve clipping the original input, and outputting other clipped primitives through std::vector.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Point& P1, ParameterizationCylindricCurve& PC1,
		std::vector<ParameterizationCylindricCurve>& sub_PC1s);
	int get_intersections(
		std::vector<Point>& P1s,
		std::vector<ParameterizationCylindricCurve>& PC1s);

	/* Get the intersections between Line L1 and Line L2.
	*  This function will involve clipping the original input, and outputting other clipped lines through std::vector.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Line& L1, Line& L2,
		std::vector<Line>& sub_L1s,
		std::vector<Line>& sub_L2s);
	int get_intersections(
		std::vector<Line>& L1s,
		std::vector<Line>& L2s);

	/* Get the intersections between Line L1 and Circle c1.
	*  This function will involve clipping the original input, and outputting other clipped lines through std::vector.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Line& L1, ParameterizationCircle& c1,
		std::vector<Line>& sub_L1s,
		std::vector<ParameterizationCircle>& sub_c1s);
	int get_intersections(
		std::vector<Line>& L1s,
		std::vector<ParameterizationCircle>& c1s);

	/* Get the intersections between Line L1 and CylindricCurve PC1.
	*  This function will involve clipping the original input, and outputting other clipped lines through std::vector.
	*  Return the number of intersections.
	*/
	int get_intersections(
		Line& L1, ParameterizationCylindricCurve& PC1,
		std::vector<Line>& sub_L1s,
		std::vector<ParameterizationCylindricCurve>& sub_PC1s);
	int get_intersections(
		std::vector<Line>& L1s,
		std::vector<ParameterizationCylindricCurve>& PC1s);

	/* Get the intersections between Circle c1 and Circle c2.
	*  This function will involve clipping the original input, and outputting other clipped lines through std::vector.
	*  Return the number of intersections.
	*/
	int get_intersections(
		ParameterizationCircle& c1, ParameterizationCircle& c2,
		std::vector<ParameterizationCircle>& sub_c1s,
		std::vector<ParameterizationCircle>& sub_c2s);
	int get_intersections(
		std::vector<ParameterizationCircle>& c1s,
		std::vector<ParameterizationCircle>& c2s);

	/* Get the intersections between Circle c1 and CylindricCurve PC1.
	*  This function will involve clipping the original input, and outputting other clipped lines through std::vector.
	*  Return the number of intersections.
	*/
	int get_intersections(
		ParameterizationCircle& c1, ParameterizationCylindricCurve& PC1,
		std::vector<ParameterizationCircle>& sub_c1s,
		std::vector<ParameterizationCylindricCurve>& sub_PC1s);
	int get_intersections(
		std::vector<ParameterizationCircle>& c1s,
		std::vector<ParameterizationCylindricCurve>& PC1s);

	/* Get the intersections between Cylindric PC1 and CylindricCurve PC2.
	*	This function will involve clipping the original input, and outputting other clipped lines through std::vector.
	*  Return the number of intersections.
	*/
	int get_intersections( // PC1 and PC2 must on the same cylinder
		ParameterizationCylindricCurve& PC1, ParameterizationCylindricCurve& PC2,
		std::vector<ParameterizationCylindricCurve>& sub_PC1s,
		std::vector<ParameterizationCylindricCurve>& sub_PC2s,
		double t_step = 1.);
	int get_intersections(
		std::vector<ParameterizationCylindricCurve>& PC1s,
		std::vector<ParameterizationCylindricCurve>& PC2s,
		double t_step = 1.);

	// -----------------------just use for testing

	/* Get the intersections between Plane P1 and other primitives.
	*  Return the number of all intersections.
	*/
	int get_intersections(
		Plane& P1,
		std::vector<Plane>& planes, std::vector<Cylinder>& cylinders, std::vector<Sphere>& spheres,
		std::vector<Point>& res_points, std::vector<Line>& res_lines, std::vector<ParameterizationCircle>& res_circles, std::vector<ParameterizationCylindricCurve>& res_curves);

	/* Get the intersections between Cylinder C1 and other primitives.
	*  Return the number of all intersections.
	*/
	int get_intersections(
		Cylinder& C1,
		std::vector<Plane>& planes, std::vector<Cylinder>& cylinders, std::vector<Sphere>& spheres,
		std::vector<Point>& res_points, std::vector<Line>& res_lines, std::vector<ParameterizationCylindricCurve>& res_curves);

	/* Get the intersections between Sphere S1 and other primitives.
	*  Return the number of all intersections.
	*/
	int get_intersections(
		Sphere& S1,
		std::vector<Plane>& planes, std::vector<Cylinder>& cylinders, std::vector<Sphere>& spheres,
		std::vector<Point>& res_points, std::vector<ParameterizationCircle>& res_circles, std::vector<ParameterizationCylindricCurve>& res_curves);

	// -----------------------output for debug

	/* Output the primitives as the meshes and export to an OBJ file.
	*/
	void write_result_mesh(
		std::string output_file_path,
		std::vector<Plane>& planes, std::vector<Cylinder>& cylinders, std::vector<Sphere>& spheres,
		double scale = 1.);
	void write_result_mesh(
		std::string output_file_path,
		std::vector<Plane>& planes,
		double scale = 1.);
	void write_result_mesh(
		std::string output_file_path,
		std::vector<Cylinder>& cylinders,
		double scale = 1.);
	void write_result_mesh(
		std::string output_file_path,
		std::vector<Sphere>& spheres,
		double scale = 1.);

	/* Output the primitives as the points and export to an OBJ file.
	*/
	void write_result_points(
		std::string output_file_path,
		std::vector<Point>& points, std::vector<Line>& lines, std::vector<ParameterizationCircle>& circles, std::vector<ParameterizationCylindricCurve>& curves,
		double scale = 1.);
	void write_result_points(
		std::string output_file_path,
		std::vector<Point>& points, std::vector<Line>& lines, std::vector<ParameterizationCylindricCurve>& curves,
		double scale = 1.);

}

#endif
