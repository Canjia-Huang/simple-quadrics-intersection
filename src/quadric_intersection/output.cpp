#include "simple_quadrics_intersection.h"

#include <iostream>
#include <fstream>
#include <random>

namespace QuadricsIntersection 
{
	std::vector<Eigen::Vector3d> rand_color_bar;

	void create_rand_color_bar(int num) {
		std::default_random_engine random(time(NULL));
		std::uniform_real_distribution<double> rand(0., 1.);
		rand_color_bar.reserve(num);
		for (int i = 0; i < num; ++i) {
			rand_color_bar.push_back(
				Eigen::Vector3d(
					rand(random) * 255,
					rand(random) * 255,
					rand(random) * 255));
		}
	}

	void Plane::output_model(
		std::vector<Eigen::Vector3d>& points,
		std::vector<Eigen::Vector3i>& faces,
		double w, double h
	) {
		// init
		std::vector<Eigen::Vector3d>().swap(points);
		std::vector<Eigen::Vector3i>().swap(faces);

		if (vertices.size() == 3) {
			points = vertices;
			faces.push_back(Eigen::Vector3i(0, 1, 2));
		}
		else {
			Eigen::Vector3d u = get_perpendicular_normal(nor_);
			Eigen::Vector3d v = nor_.cross(u).normalized();
			points.push_back(cor_ + w * u + h * v);
			points.push_back(cor_ + w * u - h * v);
			points.push_back(cor_ - w * u + h * v);
			points.push_back(cor_ - w * u - h * v);
			faces.push_back(Eigen::Vector3i(0, 2, 1));
			faces.push_back(Eigen::Vector3i(1, 2, 3));
		}
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

		// init
		std::vector<Eigen::Vector3d>().swap(points);
		std::vector<Eigen::Vector3i>().swap(faces);
		points.reserve(seg * 2);
		faces.reserve(seg * 2);

		double rot_angle = (t_ub_ - t_lb_) / seg;
		double lower_h = (s_lb_ < -0.5 * SQI_INFTY) ? -0.5 * h : s_lb_;
		double upper_h = (s_ub_ > 0.5 * SQI_INFTY) ? 0.5 * h : s_ub_;
		for (int i = 0; i < seg + 1; ++i) {
			// Eigen::Vector3d ru = slerp(u_, nor_, i * rot_angle);
			double angle = t_lb_ + i * rot_angle;
			Eigen::Vector3d ru = u_ * std::cos(ang2rad(angle)) + v_ * std::sin(ang2rad(angle));
			points.push_back(cor_ + lower_h * nor_ + r_ * ru);
			points.push_back(cor_ + upper_h * nor_ + r_ * ru);
		}
		for (int i = 0; i < seg; ++i) {
			int ri = 2 * i;

			if (std::abs(t_ub_ - t_lb_ - 360) < SQI_EPS && i == seg - 1) {
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

		// init
		std::vector<Eigen::Vector3d>().swap(points);
		std::vector<Eigen::Vector3i>().swap(faces);
		// points.reserve((h_seg - 2) * r_seg + 2);
		// faces.reserve((h_seg - 2) * r_seg * 2);
		double h_angle = (s_ub_ - s_lb_) / (h_seg - 1);
		double r_angle = 360. / r_seg;

		for (double phi = s_lb_, phi_end = s_ub_ + SQI_EPS; phi < phi_end; phi += h_angle) {
			if (phi < SQI_EPS) {
				points.push_back(get_point(phi, 0));
			}
			else if (std::abs(phi - 180) < SQI_EPS) {
				points.push_back(get_point(phi, 0));
			}
			else {
				for (double theta = 0, theta_end = 360 - SQI_EPS; theta < theta_end; theta += r_angle) {
					points.push_back(get_point(phi, theta));
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
				else if (std::abs(180 - (s_ub_ - s_lb_)) < SQI_EPS && h == h_seg - 1) {
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

	void Line::output_connections(
		std::vector<Eigen::Vector3d>& points,
		double sample_dis
	) {
		// init
		std::vector<Eigen::Vector3d>().swap(points);

		double line_length = s_ub_ - s_lb_;
		double sample_num = std::ceil(line_length / sample_dis);
		
		for (int i = 0; i < sample_num; ++i) {
			points.push_back(get_point(s_lb_ + (double(i) / sample_num * line_length)));
		}
		points.push_back(get_point(s_ub_));
	}

	void Line::output_points(
		std::vector<Eigen::Vector3d>& points,
		double l, int seg
	) {
		// init
		std::vector<Eigen::Vector3d>().swap(points);

		double s_lb = (s_lb_ < -0.5 * SQI_INFTY) ? -0.5 * l : s_lb_;
		double s_ub = (s_ub_ > 0.5 * SQI_INFTY) ? 0.5 * l : s_ub_;

		double l_step = (s_ub - s_lb) / seg;
		for (double cur_s = s_lb; cur_s < s_ub; cur_s += l_step) {
			points.push_back(get_point(cur_s));
		}
		points.push_back(get_point(s_ub));
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

	void ParameterizationCylindricCurve::output_connections(
		std::vector<Eigen::Vector3d>& points,
		double sample_dis
	) {
		double t_length = t_ub_ - t_lb_;
		double sample_num = std::ceil(t_length / sample_dis);

		for (int i = 0; i < sample_num; ++i) {
			double t = t_lb_ + (double(i) / sample_num * t_length);
			std::vector<double> ss;
			get_s(t, ss);

			for (double s : ss) {
				if (s >= s_lb_ && s <= s_ub_
					&& s >= C_.s_lb() && s <= C_.s_ub()
					) points.push_back(get_point(s, t));
			}
		}
		{
			std::vector<double> ss;
			get_s(t_ub_, ss);
			if (ss.size() == 0) return;
			if (ss[0] >= s_lb_ && ss[0] <= s_ub_
				&& ss[0] >= C_.s_lb() && ss[0] <= C_.s_ub()
				) points.push_back(get_point(ss[0], t_ub_));
		}
	}

	void ParameterizationCylindricCurve::output_points(
		std::vector<Eigen::Vector3d>& points,
		double t_step
	) {
		std::vector<Eigen::Vector3d>().swap(points);

		double cur_t_lb = t_lb_;// std::max(t_lb_, C_.t_lb());
		double cur_t_ub = t_ub_;// std::min(t_ub_, C_.t_ub());

		for (double cur_t = cur_t_lb + SQI_EPS; cur_t < cur_t_ub; cur_t += t_step) {
			std::vector<double> ss;
			get_s(cur_t, ss);

			for (double s : ss) {
				if (s >= s_lb_ && s <= s_ub_ 
					&& s >= C_.s_lb() && s <= C_.s_ub()
					) {
					points.push_back(get_point(s, cur_t));
				}
			}
		}

		std::vector<double> ss;
		get_s(cur_t_ub - SQI_EPS, ss);
		for (double s : ss) {
			if (s >= s_lb_ && s <= s_ub_
				&& s >= C_.s_lb() && s <= C_.s_ub()
				) {
				points.push_back(get_point(s, cur_t_ub - SQI_EPS));
			}
		}
	}

	void write_result_mesh(
		std::string output_file_path,
		std::vector<Plane>& planes, std::vector<Cylinder>& cylinders, std::vector<Sphere>& spheres,
		double scale
	) {
		SQI_VERBOSE_ONLY_TITLE("outputting result mesh");

		int primitive_num = planes.size() + cylinders.size() + spheres.size();
		if (rand_color_bar.size() < primitive_num) {
			create_rand_color_bar(primitive_num);
		}
		int color_num = 0;

		std::ofstream out(output_file_path);

		std::vector<Eigen::Vector3d> points;
		std::vector<Eigen::Vector3i> faces;
		int total_points_nb = 0;
		for (Plane P : planes) {
			Eigen::Vector3d color = rand_color_bar[color_num++];
			P.output_model(points, faces, scale, scale);

			for (Eigen::Vector3d p : points) {
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
			for (Eigen::Vector3i f : faces) {
				out << "f" << " " << (f + Eigen::Vector3i(total_points_nb + 1, total_points_nb + 1, total_points_nb + 1)).transpose() << std::endl;
			}
			total_points_nb += points.size();
		}
		for (Cylinder C : cylinders) {
			Eigen::Vector3d color = rand_color_bar[color_num++];
			C.output_model(points, faces, 32, 2 * scale);

			for (Eigen::Vector3d p : points) {
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
			for (Eigen::Vector3i f : faces) {
				out << "f" << " " << (f + Eigen::Vector3i(total_points_nb + 1, total_points_nb + 1, total_points_nb + 1)).transpose() << std::endl;
			}
			total_points_nb += points.size();
		}
		for (Sphere S : spheres) {
			Eigen::Vector3d color = rand_color_bar[color_num++];
			S.output_model(points, faces);

			for (Eigen::Vector3d p : points) {
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
			for (Eigen::Vector3i f : faces) {
				out << "f" << " " << (f + Eigen::Vector3i(total_points_nb + 1, total_points_nb + 1, total_points_nb + 1)).transpose() << std::endl;
			}
			total_points_nb += points.size();
		}
		out.close();
	}

	void write_result_mesh(
		std::string output_file_path,
		std::vector<Plane>& planes,
		double scale
	) {
		std::vector<Cylinder> cylinders; std::vector<Sphere> spheres;
		write_result_mesh(output_file_path,
			planes, cylinders, spheres,
			scale);
	}

	void write_result_mesh(
		std::string output_file_path,
		std::vector<Cylinder>& cylinders,
		double scale
	) {
		std::vector<Plane> planes; std::vector<Sphere> spheres;
		write_result_mesh(output_file_path,
			planes, cylinders, spheres,
			scale);
	}

	void write_result_mesh(
		std::string output_file_path,
		std::vector<Sphere>& spheres,
		double scale
	) {
		std::vector<Plane> planes; std::vector<Cylinder> cylinders;
		write_result_mesh(output_file_path,
			planes, cylinders, spheres,
			scale);
	}

	void write_result_points(
		std::string output_file_path,
		std::vector<Point>& points, std::vector<Line>& lines, std::vector<ParameterizationCircle>& circles, std::vector<ParameterizationCylindricCurve>& curves,
		double scale
	) {
		SQI_VERBOSE_ONLY_TITLE("outputting result points");

		int primitive_num = points.size() + lines.size() + circles.size() + curves.size();
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
			L.output_points(output_points, scale, 200);

			for (Eigen::Vector3d p : output_points) {
#ifdef USE_FOR_OFFSET_MESH_GENERATION
				if (std::abs(p.x()) > 1 || std::abs(p.y()) > 1 || std::abs(p.z()) > 1) continue;
#endif
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
		}

		for (ParameterizationCircle c : circles) {
			Eigen::Vector3d color = rand_color_bar[color_num++];

			std::vector<Eigen::Vector3d> ops;
			c.output_points(ops);

			for (Eigen::Vector3d p : ops) {
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
		}

		for (ParameterizationCylindricCurve PC : curves) {
			Eigen::Vector3d color = rand_color_bar[color_num++];

			std::vector<Eigen::Vector3d> output_points;
			PC.output_points(output_points);

			for (Eigen::Vector3d p : output_points) {
				out << "v" << " " << p.transpose() << " " << color.transpose() << std::endl;
			}
		}
	}

	void write_result_points(
		std::string output_file_path,
		std::vector<Point>& points, std::vector<Line>& lines, std::vector<ParameterizationCylindricCurve>& curves,
		double scale) {
		std::vector<ParameterizationCircle> circles;
		write_result_points(output_file_path, points, lines, circles, curves, scale);
	}

	void write_result_points(
		std::string output_file_path,
		std::vector<ParameterizationCircle>& circles, std::vector<ParameterizationCylindricCurve>& curves,
		double scale) {
		std::vector<Point> points;
		std::vector<Line> lines;
		write_result_points(output_file_path, points, lines, circles, curves, scale);
	}

	void write_result_points(
		std::string output_file_path,
		std::vector<Line>& lines, std::vector<ParameterizationCylindricCurve>& curves,
		double scale) {
		std::vector<Point> points;
		std::vector<ParameterizationCircle> circles;
		write_result_points(output_file_path, points, lines, circles, curves, scale);
	}

	void write_result_points(
		std::string output_file_path,
		std::vector<Line>& lines,
		double scale) {
		std::vector<Point> points;
		std::vector<ParameterizationCircle> circles;
		std::vector<ParameterizationCylindricCurve> curves;
		write_result_points(output_file_path, points, lines, circles, curves, scale);
	}

	void write_result_points(
		std::string output_file_path,
		std::vector<ParameterizationCircle>& circles,
		double scale) {
		std::vector<Point> points;
		std::vector<Line> lines;
		std::vector<ParameterizationCylindricCurve> curves;
		write_result_points(output_file_path, points, lines, circles, curves, scale);
	}

	void write_result_points(
		std::string output_file_path,
		std::vector<ParameterizationCylindricCurve>& curves,
		double scale) {
		std::vector<Point> points;
		std::vector<Line> lines;
		std::vector<ParameterizationCircle> circles;
		write_result_points(output_file_path, points, lines, circles, curves, scale);
	}
}