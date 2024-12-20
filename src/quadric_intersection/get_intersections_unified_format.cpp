#include "simple_quadrics_intersection.h"

namespace QuadricsIntersection
{
	int get_intersections(
		Plane& P1, Plane& P2,
		std::vector<Point>& Ps, std::vector<Line>& Ls, std::vector<ParameterizationCircle>& cs, std::vector<ParameterizationCylindricCurve>& Cs,
		double limit_angle
	) {
		std::vector<Point>().swap(Ps);
		std::vector<ParameterizationCircle>().swap(cs);
		std::vector<ParameterizationCylindricCurve>().swap(Cs);

		return get_intersections(P1, P2, Ls, limit_angle);
	}

	int get_intersections(
		Cylinder& C1, Plane& P1,
		std::vector<Point>& Ps, std::vector<Line>& Ls, std::vector<ParameterizationCircle>& cs, std::vector<ParameterizationCylindricCurve>& Cs,
		double limit_angle
	) {
		std::vector<Point>().swap(Ps);
		std::vector<ParameterizationCircle>().swap(cs);

		return get_intersections(C1, P1, Ls, Cs, limit_angle);
	}

	int get_intersections(
		Sphere& S1, Plane& P1,
		std::vector<Point>& Ps, std::vector<Line>& Ls, std::vector<ParameterizationCircle>& cs, std::vector<ParameterizationCylindricCurve>& Cs,
		double limit_angle
	) {
		std::vector<Line>().swap(Ls);
		std::vector<ParameterizationCylindricCurve>().swap(Cs);

		return get_intersections(S1, P1, Ps, cs, limit_angle);
	}

	int get_intersections(
		Plane& P1, Cylinder& C1,
		std::vector<Point>& Ps, std::vector<Line>& Ls, std::vector<ParameterizationCircle>& cs, std::vector<ParameterizationCylindricCurve>& Cs,
		double limit_angle
	) {
		std::vector<Point>().swap(Ps);
		std::vector<ParameterizationCircle>().swap(cs);

		return get_intersections(P1, C1, Ls, Cs, limit_angle);
	}

	int get_intersections(
		Cylinder& C1, Cylinder& C2,
		std::vector<Point>& Ps, std::vector<Line>& Ls, std::vector<ParameterizationCircle>& cs, std::vector<ParameterizationCylindricCurve>& Cs,
		double limit_angle
	) {
		std::vector<ParameterizationCircle>().swap(cs);

		return get_intersections(C1, C2, Ps, Ls, Cs);
	}

	int get_intersections(
		Sphere& S1, Cylinder& C1,
		std::vector<Point>& Ps, std::vector<Line>& Ls, std::vector<ParameterizationCircle>& cs, std::vector<ParameterizationCylindricCurve>& Cs,
		double limit_angle
	) {
		std::vector<Line>().swap(Ls);
		std::vector<ParameterizationCircle>().swap(cs);

		return get_intersections(S1, C1, Ps, Cs, limit_angle);
	}

	int get_intersections(
		Plane& P1, Sphere& S1,
		std::vector<Point>& Ps, std::vector<Line>& Ls, std::vector<ParameterizationCircle>& cs, std::vector<ParameterizationCylindricCurve>& Cs,
		double limit_angle
	) {
		std::vector<Line>().swap(Ls);
		std::vector<ParameterizationCylindricCurve>().swap(Cs);

		return get_intersections(P1, S1, Ps, cs, limit_angle);
	}

	int get_intersections(
		Cylinder& C1, Sphere& S1,
		std::vector<Point>& Ps, std::vector<Line>& Ls, std::vector<ParameterizationCircle>& cs, std::vector<ParameterizationCylindricCurve>& Cs,
		double limit_angle
	) {
		std::vector<Line>().swap(Ls);
		std::vector<ParameterizationCircle>().swap(cs);

		return get_intersections(C1, S1, Ps, Cs);
	}

	int get_intersections(
		Sphere& S1, Sphere& S2,
		std::vector<Point>& Ps, std::vector<Line>& Ls, std::vector<ParameterizationCircle>& cs, std::vector<ParameterizationCylindricCurve>& Cs,
		double limit_angle
	) {
		std::vector<Line>().swap(Ls);
		std::vector<ParameterizationCylindricCurve>().swap(Cs);

		return get_intersections(S1, S2, Ps, cs, limit_angle);
	}
}