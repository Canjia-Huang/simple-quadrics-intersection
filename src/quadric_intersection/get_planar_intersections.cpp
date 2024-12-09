#include "simple_quadrics_intersection.h"

namespace QuadricsIntersection
{
	int get_intersections(
		Plane& P1, Plane& P2,
		Line& L
	) {
		// SQI_VERBOSE_ONLY_COUT("");

		if (1 - std::abs(P1.nor().dot(P2.nor())) < SQI_EPS) { // parallel
			return 0;
		}

		Eigen::Vector3d L_nor = (P1.nor().cross(P2.nor())).normalized();

		Eigen::Vector3d P1_vec = (P1.nor().cross(L_nor)).normalized();
		double t = (P2.cor() - P1.cor()).dot(P2.nor()) / (P1_vec.dot(P2.nor()));

		Eigen::Vector3d L_cor = P1.cor() + t * P1_vec;

		L = Line(L_cor, L_nor);

		return 1;
	}
}