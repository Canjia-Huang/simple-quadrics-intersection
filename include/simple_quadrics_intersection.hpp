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
		}
		Eigen::Vector3d& cor() { return cor_; }
		Eigen::Vector3d& nor() { return nor_; }
		double& r() { return r_; }

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
	public:
	private:
	};

	// ----------------------------assessment of the intersections----------------------------
	// get the branches that primitives intersected, return the number of branches
	int get_intersection(
		Plane P1, Cylinder C1,
		std::vector<Line>& lines,
		std::vector<Circle>& circles
	) {
		SQI_VERBOSE_ONLY_COUT("plane" << " " << "cor:" << P1.cor().transpose() << " " << "nor:" << P1.nor().transpose());
		SQI_VERBOSE_ONLY_COUT("cylinder" << " " << "cor:" << C1.cor().transpose() << " " << "nor:" << C1.nor().transpose() << " " << "r:" << C1.r());



		return lines.size() + circles.size();
	}

	int get_intersections(
		Cylinder C1, Cylinder C2,
		std::vector<Line>& lines,
		std::vector<Circle>& circles,
		std::vector<ParameterizationCurve>& parameterization_curves
	) {
		SQI_VERBOSE_ONLY_COUT("cylinder1" << " " << "cor:" << C1.cor().transpose() << " " << "nor:" << C1.nor().transpose() << " " << "r:" << C1.r());
		SQI_VERBOSE_ONLY_COUT("cylinder2" << " " << "cor:" << C2.cor().transpose() << " " << "nor:" << C2.nor().transpose() << " " << "r:" << C2.r());

		// init
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
			return 0;
		}

		if (abs(abs(C1.nor().dot(C2.nor())) - 1) < SQI_EPS) { // 2 axes are parallel
			if (axes_sq_dis < SQI_EPS) { // overlap
				SQI_VERBOSE_ONLY_COUT("two cylinders are overlap, may error!");
				return 0;
			}
			else if (axes_sq_dis < C1_C2_sq_r - SQI_EPS) { // inside, result 2 lines
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
				lines.push_back(
					Line(
						0.5 * (C1.cor() + C2.cor()),
						C1.nor()
					));
			}
		}
		else { // 2 axes are not parallel
			Eigen::Vector3d perpendicular_v = (C1.nor().cross(C2.nor())).normalized();
			Line L1(
				C1.cor() + C1.r() * perpendicular_v, 
				C1.nor());
			Line L2(
				C1.cor() - C1.r() * perpendicular_v,
				C1.nor());

			// check the intersections between line1 and cylinder2
			double L1_axis_dis = abs((L1.nor().cross(C2.nor())).normalized().dot(L1.cor() - C2.cor()));
			double L2_axis_dis = abs((L2.nor().cross(C2.nor())).normalized().dot(L2.cor() - C2.cor()));
			int L1_status, L2_status; // 0: do not intersect, 1: tangent, 2: two-point intersection

			if (L1_axis_dis > C2.r() + EPS) L1_status = 0;
			else if (L1_axis_dis < C2.r() - EPS) L1_status = 2;
			else L1_status = 1;

			if (L2_axis_dis > C2.r() + EPS) L2_status = 0;
			else if (L2_axis_dis < C2.r() - EPS) L2_status = 2;
			else L2_status = 1;
		}

		return lines.size() + circles.size() + parameterization_curves.size();
	}
}

#endif
