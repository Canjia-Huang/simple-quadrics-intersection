#include "simple_quadrics_intersection.h"

#include <iostream>
#include <fstream>

namespace QuadricsIntersection 
{
	std::vector<Eigen::Vector3d> rand_color_bar;

	void create_rand_color_bar(int num) {
		srand(int(time(NULL)));
		rand_color_bar.reserve(num);
		for (int i = 0; i < num; ++i) {
			rand_color_bar.push_back(
				Eigen::Vector3d(
					double(std::rand()) / (RAND_MAX + 1) * 255,
					double(std::rand()) / (RAND_MAX + 1) * 255,
					double(std::rand()) / (RAND_MAX + 1) * 255));
		}
	}

	void Plane::output_model(
		std::vector<Eigen::Vector3d>& points,
		std::vector<Eigen::Vector3i>& faces,
		double w, double h
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

	void Cylinder::output_model(
		std::vector<Eigen::Vector3d>& points,
		std::vector<Eigen::Vector3i>& faces,
		int seg, double h
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

	void Sphere::output_model(
		std::vector<Eigen::Vector3d>& points,
		std::vector<Eigen::Vector3i>& faces,
		int h_seg, int r_seg
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

	void Line::output_model(
		std::vector<Eigen::Vector3d>& points,
		std::vector<Eigen::Vector2i>& lines,
		double l
	) {
		SQI_VERBOSE_ONLY_COUT("");

		// init
		std::vector<Eigen::Vector3d>().swap(points);
		std::vector<Eigen::Vector2i>().swap(lines);
		points.reserve(2);
		lines.reserve(1);

		points.push_back(get_point(std::min(0.5 * l, s_ub_)));
		points.push_back(get_point(std::max(-0.5 * l, s_lb_)));
		lines.push_back(Eigen::Vector2i(0, 1));
	}

	void Line::output_points(
		std::vector<Eigen::Vector3d>& points,
		double l, int seg
	) {
		// init
		std::vector<Eigen::Vector3d>().swap(points);

		double l_step = l / seg;
		for (double cur_l = std::max(-0.5 * l, s_lb_), cur_l_end = std::min(0.5 * l, s_ub_); cur_l < cur_l_end; cur_l += l_step) {
			points.push_back(get_point(cur_l));
		}
		points.push_back(get_point(std::min(0.5 * l, s_ub_)));
	}

	void ParameterizationCircle::output_points(
		std::vector<Eigen::Vector3d>& points,
		double t_step
	) {
		std::vector<Eigen::Vector3d>().swap(points);

		for (double t = t_lb_, t_end = t_ub_ - SQI_EPS; t < t_end; t += t_step) {
			points.push_back(get_point(t));
		}
		points.push_back(get_point(t_ub_));
	}

	void ParameterizationCylindricPoint::output_points(
		Cylinder& C,
		std::vector<Eigen::Vector3d>& points
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		std::vector<Eigen::Vector3d>().swap(points);
		points.push_back(C.get_point(s_, t_));
	}

	void ParameterizationCylindricLine::output_points(
		Cylinder& C,
		std::vector<Eigen::Vector3d>& points,
		double h, int h_seg 
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

	void ParameterizationCylindricCurve::output_points(
		Cylinder& C,
		std::vector<Eigen::Vector3d>& points,
		double t_step
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		std::vector<Eigen::Vector3d>().swap(points);

		for (double cur_t = t_lb_ + SQI_EPS; cur_t < t_ub_; cur_t += t_step) {
			std::vector<double> ss;
			get_s(cur_t, ss);

			for (double s : ss) {
				if (s >= s_lb_ && s <= s_ub_) {
					points.push_back(C.get_point(s, cur_t));
				}
			}
		}

		std::vector<double> ss;
		get_s(t_ub_ - SQI_EPS, ss);
		for (double s : ss) {
			if (s >= s_lb_ && s <= s_ub_) {
				points.push_back(C.get_point(s, t_ub_ - SQI_EPS));
			}
		}
	}

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
		for (Sphere S : spheres) {
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

	void write_planar_result_points(
		std::string output_file_path,
		std::vector<Point>& points,
		std::vector<Line>& lines
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		int primitive_num = points.size() + lines.size();
		if (rand_color_bar.size() < primitive_num) {
			create_rand_color_bar(primitive_num);
		}

		std::ofstream out(output_file_path);
		int color_num = 0;

		for (Point P : points) {
			Eigen::Vector3d color = rand_color_bar[color_num++];

			out << "v" << " " << P.cor().transpose() << " " << color.transpose() << std::endl;
		}

		for (Line L : lines) {
			Eigen::Vector3d color = rand_color_bar[color_num++];

			std::vector<Eigen::Vector3d> output_points;
			L.output_points(output_points);

			for (Eigen::Vector3d p : output_points) {
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
		}
	}

	void write_cylinder_result_points(
		std::string output_file_path,
		Cylinder& C,
		std::vector<Point>& points,
		std::vector<Line>& lines,
		std::vector<ParameterizationCylindricCurve>& curves
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		int primitive_num = points.size() + lines.size() + curves.size();
		if (rand_color_bar.size() < primitive_num) {
			create_rand_color_bar(primitive_num);
		}

		std::ofstream out(output_file_path);
		int color_num = 0;

		for (Point P : points) {
			Eigen::Vector3d color = rand_color_bar[color_num++];

			out << "v" << " " << P.cor().transpose() << " " << color.transpose() << std::endl;
		}

		for (Line L : lines) {
			Eigen::Vector3d color = rand_color_bar[color_num++];

			std::vector<Eigen::Vector3d> output_points;
			L.output_points(output_points);

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

	void write_cylinder_result_points(
		std::string output_file_path,
		Cylinder& C,
		std::vector<Line>& lines
	) {
		std::vector<Point> points;
		std::vector<ParameterizationCylindricCurve> curves;
		write_cylinder_result_points(output_file_path, C, points, lines, curves);
	}

	void write_cylinder_result_points(
		std::string output_file_path,
		Cylinder& C,
		std::vector<ParameterizationCylindricCurve>& curves
	) {
		std::vector<Point> points;
		std::vector<Line> lines;
		write_cylinder_result_points(output_file_path, C, points, lines, curves);
	}

	void write_cylinders_result_points(
		std::string output_file_path,
		std::vector<std::vector<ParameterizationCylindricCurve>>& c_curves, std::vector<Cylinder>& cylinders
	) {
		// SQI_VERBOSE_ONLY_TITLE("");

		int primitive_num = 0;
		for (int i = 0, i_end = c_curves.size(); i < i_end; ++i) primitive_num += c_curves[i].size();
		if (rand_color_bar.size() < primitive_num) {
			create_rand_color_bar(primitive_num);
		}

		std::ofstream out(output_file_path);
		int color_num = 0;

		for (int i = 0, i_end = c_curves.size(); i < i_end; ++i) {
			std::vector<ParameterizationCylindricCurve>* curves = &c_curves[i];

			for (ParameterizationCylindricCurve PC : *curves) {
				Eigen::Vector3d color = rand_color_bar[color_num++];

				std::vector<Eigen::Vector3d> output_points;
				PC.output_points(cylinders[i], output_points);

				for (Eigen::Vector3d p : output_points) {
					out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
				}
			}
		}

		out.close();
	}

	void write_sphere_result_points(
		std::string output_file_path,
		std::vector<Point>& points,
		std::vector<ParameterizationCircle>& circles
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		int primitive_num = points.size() + circles.size();
		if (rand_color_bar.size() < primitive_num) {
			create_rand_color_bar(primitive_num);
		}

		std::ofstream out(output_file_path);
		int color_num = 0;

		for (Point p : points) {
			Eigen::Vector3d color = rand_color_bar[color_num++];

			out << "v" << " " << p.cor().transpose() << " " << color.transpose() << std::endl;
		}
		for (ParameterizationCircle c : circles) {
			Eigen::Vector3d color = rand_color_bar[color_num++];

			std::vector<Eigen::Vector3d> ops;
			c.output_points(ops);

			for (Eigen::Vector3d p : ops) {
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
		}

		out.close();
	}
}