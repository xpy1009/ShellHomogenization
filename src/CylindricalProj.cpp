#include <CylindricalProj.h>


Eigen::MatrixX3d CylindricalProj::proj(const Eigen::MatrixX3d& w, double k, double alpha)
{
	if (std::abs(k) < m_sEps) {
		return w;
    }

	Eigen::Matrix3d R;
	R << cos(alpha),sin(alpha),0,
		-sin(alpha),cos(alpha),0,
		 0,0,1;

	Eigen::MatrixX3d x(w.rows(), w.cols());
	for(int i = 0; i < w.rows(); ++i) {
		const Eigen::RowVector3d wi = w.row(i);
		const Eigen::RowVector3d wir = wi * R.transpose();
		Eigen::RowVector3d xir;
		xir << (k * wir(2) + 1) / k * sin(k * wir(0)), 
            wir(1), 
            ((k * wir(2) + 1) * cos(k * wir(0)) - 1) / k;
		const Eigen::RowVector3d xi = xir * R;
		x.row(i) = xi;
	}

	return x;
}

Eigen::SparseMatrix<double> CylindricalProj::gradient(const Eigen::MatrixX3d& w, double k, double alpha)
{
    const int n = w.rows();
    Eigen::SparseMatrix<double> jacobian(3 * n, 3 * n);
	if (std::abs(k) < m_sEps) {
		jacobian.setIdentity();
		return jacobian;
	}

	Eigen::Matrix3d R;
	R << cos(alpha),sin(alpha),0,
		-sin(alpha),cos(alpha),0,
		 0,0,1;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n * 3 * 3);
    for(int i = 0; i < n; ++i) {
		const Eigen::RowVector3d wi = w.row(i);
		const Eigen::RowVector3d wir = wi * R.transpose();
        Eigen::Matrix3d jr;
		jr << (k * wir(2) + 1) * cos(k * wir(0)), 0, sin(k * wir(0)),
			0, 1, 0,
			-(k * wir(2) + 1) * sin(k * wir(0)), 0, cos(k * wir(0));
		jr = R.transpose() * jr * R;

        for(int d = 0; d < 3; ++d) {
            for(int dd = 0; dd < 3; ++dd) {
                triplets.emplace_back(3 * i + d, 3 * i + dd, jr(d, dd));
            }
        }

        // Note: Jacobian wrt curvature seems unnecessary
		// Eigen::Vector3d dxrdk;
		// dxrdk << ((wir(2) * sin(k * wir(0)) + wir(0) * cos(k * wir(0)) * (k * wir(2) + 1)) * k - (k * wir(2) + 1) * sin(k * wir(0))) / (k * k), 0, 
		// 		((wir(2) * cos(k * wir(0)) - wir(0) * sin(k * wir(0)) * (k * wir(2) + 1)) * k - ((k * wir(2) + 1) * cos(k * wir(0)) - 1)) / (k * k);
		// Eigen::Vector3d dxdk = R.transpose() * dxrdk;
		// for (int d = 0; d < 3; d++)
		// 	triplets.push_back(Eigen::Triplet<double>(i+d, w.size(), dxdk(d)));
    }

    jacobian.setFromTriplets(triplets.begin(), triplets.end());

    return jacobian;

}

Eigen::SparseMatrix<double> CylindricalProj::hessian(const Eigen::MatrixX3d& w, double k, double alpha, const Eigen::VectorXd& gradient)
{
    const int n = w.rows();
    Eigen::SparseMatrix<double> hessian(3 * n, 3 * n);
	if (std::abs(k) < m_sEps) {
		hessian.setZero();
		return hessian;
	}

	Eigen::Matrix3d R;
	R << cos(alpha),sin(alpha),0,
		-sin(alpha),cos(alpha),0,
		 0,0,1;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n * 3 * 3);
	for(int i = 0; i < n; ++i)
	{
		const Eigen::RowVector3d wi = w.row(i);
		const Eigen::RowVector3d wir = wi * R.transpose();

		const Eigen::Matrix3d dwrdw = R;

		Eigen::Matrix3d d2xr1dwr2, d2xr2dwr2, d2xr3dwr2;
		d2xr1dwr2 << -k * (k * wir(2) + 1) * sin(k * wir(0)), 0, k * cos(k * wir(0)),
					0, 0, 0, 
                    k * cos(k * wir(0)), 0, 0;
		d2xr2dwr2.setZero();
		d2xr3dwr2 << -k * (k * wir(2) + 1) * cos(k * wir(0)), 0, -k * sin(k * wir(0)),
					0, 0, 0, 
                    -k * sin(k * wir(0)), 0, 0;
		const Eigen::Matrix3d d2xr1dw2 = dwrdw.transpose() * d2xr1dwr2 * dwrdw;
		const Eigen::Matrix3d d2xr2dw2 = dwrdw.transpose() * d2xr2dwr2 * dwrdw;
		const Eigen::Matrix3d d2xr3dw2 = dwrdw.transpose() * d2xr3dwr2 * dwrdw;

		const Eigen::Matrix3d dxdxr = R.transpose();

		const Eigen::Matrix3d d2x1dw2 = dxdxr(0,0) * d2xr1dw2 + dxdxr(0,1) * d2xr2dw2 + dxdxr(0,2) * d2xr3dw2;
		const Eigen::Matrix3d d2x2dw2 = dxdxr(1,0) * d2xr1dw2 + dxdxr(1,1) * d2xr2dw2 + dxdxr(1,2) * d2xr3dw2;
		const Eigen::Matrix3d d2x3dw2 = dxdxr(2,0) * d2xr1dw2 + dxdxr(2,1) * d2xr2dw2 + dxdxr(2,2) * d2xr3dw2;
		
		const Eigen::Matrix3d h1 = gradient[3*i] * d2x1dw2;
		const Eigen::Matrix3d h2 = gradient[3*i+1] * d2x2dw2;
		const Eigen::Matrix3d h3 = gradient[3*i+2] * d2x3dw2;

		const Eigen::Matrix3d h = h1 + h2 + h3;

		for(int d = 0; d < 3; ++d) {
			for(int dd = 0; dd < 3; ++dd) {
				triplets.emplace_back(3 * i + d, 3 * i + dd, h(d, dd));
            }
        }

        // Note: Hessian wrt curvature seems unnecessary
		// Eigen::Vector3d d2xrdk2;
		// d2xrdk2 << -wr(0)*wr(0)*wr(2)*sin(k*wr(0)) - wr(0)*wr(0)*sin(k*wr(0))/k - 2*wr(0)*cos(k*wr(0))/k/k + 2*sin(k*wr(0))/(k*k*k), 0,
		// 	(-k*k*k*wr(0)*wr(0)*wr(2)*cos(k*wr(0)) - k*k*wr(0)*wr(0)*cos(k*wr(0)) + 2*k*wr(0)*sin(k*wr(0)) + 2*cos(k*wr(0)) - 2)/(k*k*k);
		// Eigen::Vector3d d2xdk2 = R.transpose() * d2xrdk2;
		// double dedx_d2xdk2 = d2xdk2.dot(g.segment(i,3));
		// triplets.push_back(Eigen::Triplet<double>(beta.size(), beta.size(), dedx_d2xdk2));

		// Eigen::Matrix3d d2xrdwdk;
		// d2xrdwdk << wr(2)*cos(k*wr(0))-wr(0)*sin(k*wr(0))*(k*wr(2)+1), 0, wr(0)*cos(k*wr(0)),
		// 			0,0,0,
		// 		-wr(2)*sin(k*wr(0))-wr(0)*cos(k*wr(0))*(k*wr(2)+1), 0, -wr(0)*sin(k*wr(0));
		// Eigen::Matrix3d d2xdwdk = R.transpose() * d2xrdwdk * R;
		// Eigen::Vector3d dedx_d2xdwdk = d2xdwdk.transpose() * g.segment(i,3);
		// for (int d = 0; d < 3; d++)
		// {
		// 	triplets.push_back(Eigen::Triplet<double>(i+d, beta.size(), dedx_d2xdwdk(d)));
		// 	triplets.push_back(Eigen::Triplet<double>(beta.size(), i+d, dedx_d2xdwdk(d)));
		// }
		
	}

	hessian.setFromTriplets(triplets.begin(), triplets.end());

    return hessian;
}
