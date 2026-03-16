#include <GeometryLib.h>

#include <Eigen/Dense>
#include <iostream>

Eigen::MatrixX3d GeometryLib::cylindricalProject(const Eigen::MatrixX3d& w, double k, double alpha)
{
	if (abs(k) < 1e-10) {
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

Eigen::SparseMatrix<double> GeometryLib::cylindricalJacobian(const Eigen::MatrixX3d& w, double k, double alpha)
{
    const int n = w.rows();
    Eigen::SparseMatrix<double> jacobian(3 * n, 3 * n);
	if (abs(k) < 1e-10) {
		jacobian.setIdentity();
		return jacobian;
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

Eigen::SparseMatrix<double> GeometryLib::cylindricalHessian(const Eigen::MatrixX3d& w, double k, double alpha, const Eigen::VectorXd& gradient)
{
    const int n = w.rows();
    Eigen::SparseMatrix<double> hessian(3 * n, 3 * n);
	if (abs(k) < 1e-10) {
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

// double GeometryLib::area(const Eigen::Matrix<double,3,2>& X)
// {
// 	// https://mathworld.wolfram.com/TriangleArea.html (17)
// 	return 0.5 * abs(X(0,0) * (X(1,1) - X(2,1)) + X(1,0) * (X(2,1) - X(0,1)) + X(2,0) * (X(0,1) - X(1,1)));
// }

// Eigen::Matrix<double,3,2> GeometryLib::deformationGradient(const Eigen::Matrix3d& x, 
// 														   const Eigen::Matrix<double,3,2>& X,
// 														   Eigen::Matrix<double,6,9>* dFdx)
// {
// 	// F = dx/duv = dx/dq dq/dX dX/duv
// 	// dX/duv = I for orthogonal materials

// 	const Eigen::Vector3d v1 = (x.row(1) - x.row(0)).transpose();
// 	const Eigen::Vector3d v2 = (x.row(2) - x.row(0)).transpose();
// 	const Eigen::Vector2d V1 = (X.row(1) - X.row(0)).transpose();
// 	const Eigen::Vector2d V2 = (X.row(2) - X.row(0)).transpose();

// 	Eigen::Matrix<double,3,2> dxdq;
//     dxdq << v1, v2;
// 	Eigen::Matrix2d dXdq;
// 	dXdq << V1, V2;
    
// 	const Eigen::Matrix2d dqdX = dXdq.inverse();
//     const Eigen::Matrix<double,3,2> F = dxdq * dqdX;

// 	if (dFdx != nullptr) { // F in colume major
// 		Eigen::Matrix3d I; I.setIdentity();
// 		*dFdx << -dqdX.col(0).sum()*I, dqdX(0,0)*I, dqdX(1,0)*I,
// 				  -dqdX.col(1).sum()*I, dqdX(0,1)*I, dqdX(1,1)*I;
// 	}

//     return F;
// }

// Eigen::Matrix3d GeometryLib::crossMatrix(const Eigen::RowVector3d& v)
// {
// 	Eigen::Matrix3d ret;
// 	ret << 0, -v[2], v[1],
// 		   v[2], 0, -v[0],
// 		   -v[1], v[0], 0;
// 	return ret;
// }

// double GeometryLib::exteriorDihedralAngle(const Eigen::RowVector3d &x0, 
// 										  const Eigen::RowVector3d &x1, 
// 										  const Eigen::RowVector3d &x2, 
// 										  const Eigen::RowVector3d &x3,
// 										  Eigen::Matrix<double,12,1>* dtdx,
// 										  Eigen::Matrix<double,12,12>* d2tdx2)
// {
//     const Eigen::RowVector3d e0 = x1 - x0;

//     const Eigen::RowVector3d n0 = (x0 - x2).cross(x1 - x2);
//     const Eigen::RowVector3d n1 = (x1 - x3).cross(x0 - x3);

// 	Eigen::Matrix<double,9,1> dadv;
// 	Eigen::Matrix<double, 9, 9> d2adv2;
//     const double theta = angle(n0, n1, e0, (dtdx!=nullptr || d2tdx2!=nullptr) ? &dadv : nullptr, d2tdx2!=nullptr ? &d2adv2 : nullptr);

// 	if (dtdx != nullptr) {
// 		dtdx->setZero();
// 		dtdx->segment(0,3) += crossMatrix(x1 - x2) * dadv.segment(0,3);
// 		dtdx->segment(3,3) += crossMatrix(x2 - x0) * dadv.segment(0,3);
// 		dtdx->segment(6,3) += crossMatrix(x0 - x1) * dadv.segment(0,3);

// 		dtdx->segment(0,3) += crossMatrix(x3 - x1) * dadv.segment(3,3);
// 		dtdx->segment(3,3) += crossMatrix(x0 - x3) * dadv.segment(3,3);
// 		dtdx->segment(9,3) += crossMatrix(x1 - x0) * dadv.segment(3,3);
// 	}

// 	if (d2tdx2 != nullptr) {
// 		d2tdx2->setZero();
// 		const std::array<Eigen::Matrix3d, 3> vqm = {crossMatrix(x2 - x1), crossMatrix(x0 - x2), crossMatrix(x1 - x0)}; // dn0dx
// 		const std::array<Eigen::Matrix3d, 3> wqm = {crossMatrix(x1 - x3), crossMatrix(x3 - x0), crossMatrix(x0 - x1)}; // dn1dx

// 		const std::array<int, 3> vindices = { 0, 3, 6 };
// 		const std::array<int, 3> windices = { 0, 3, 9 };

// 		for (int i = 0; i < 3; ++i) {
// 			for (int j = 0; j < 3; ++j) {
// 				d2tdx2->block(vindices[i], vindices[j], 3, 3) += vqm[i].transpose() * d2adv2.block(0, 0, 3, 3) * vqm[j];
// 				d2tdx2->block(vindices[i], windices[j], 3, 3) += vqm[i].transpose() * d2adv2.block(0, 3, 3, 3) * wqm[j];
// 				d2tdx2->block(windices[i], vindices[j], 3, 3) += wqm[i].transpose() * d2adv2.block(3, 0, 3, 3) * vqm[j];
// 				d2tdx2->block(windices[i], windices[j], 3, 3) += wqm[i].transpose() * d2adv2.block(3, 3, 3, 3) * wqm[j];
// 			}

// 			d2tdx2->block(vindices[i], 3, 3, 3) += vqm[i].transpose() * d2adv2.block(0, 6, 3, 3);
// 			d2tdx2->block(3, vindices[i], 3, 3) += d2adv2.block(6, 0, 3, 3) * vqm[i];
// 			d2tdx2->block(vindices[i], 0, 3, 3) += -vqm[i].transpose() * d2adv2.block(0, 6, 3, 3);
// 			d2tdx2->block(0, vindices[i], 3, 3) += -d2adv2.block(6, 0, 3, 3) * vqm[i];

// 			d2tdx2->block(windices[i], 3, 3, 3) += wqm[i].transpose() * d2adv2.block(3, 6, 3, 3);
// 			d2tdx2->block(3, windices[i], 3, 3) += d2adv2.block(6, 3, 3, 3) * wqm[i];
// 			d2tdx2->block(windices[i], 0, 3, 3) += -wqm[i].transpose() * d2adv2.block(3, 6, 3, 3);
// 			d2tdx2->block(0, windices[i], 3, 3) += -d2adv2.block(6, 3, 3, 3) * wqm[i];

// 		}

// 		const Eigen::RowVector3d dang1 = dadv.segment(0,3).transpose();
// 		const Eigen::RowVector3d dang2 = dadv.segment(3,3).transpose();

// 		const Eigen::Matrix3d dang1mat = crossMatrix(dang1);
// 		const Eigen::Matrix3d dang2mat = crossMatrix(dang2);

// 		d2tdx2->block<3, 3>(6, 3) += dang1mat;
// 		d2tdx2->block<3, 3>(0, 3) -= dang1mat;
// 		d2tdx2->block<3, 3>(0, 6) += dang1mat;
// 		d2tdx2->block<3, 3>(3, 0) += dang1mat;
// 		d2tdx2->block<3, 3>(3, 6) -= dang1mat;
// 		d2tdx2->block<3, 3>(6, 0) -= dang1mat;

// 		d2tdx2->block<3, 3>(9, 0) += dang2mat;
// 		d2tdx2->block<3, 3>(3, 0) -= dang2mat;
// 		d2tdx2->block<3, 3>(3, 9) += dang2mat;
// 		d2tdx2->block<3, 3>(0, 3) += dang2mat;
// 		d2tdx2->block<3, 3>(0, 9) -= dang2mat;
// 		d2tdx2->block<3, 3>(9, 3) -= dang2mat;
// 	}

// 	return theta;
// }

// double GeometryLib::angle(const Eigen::RowVector3d &v, 
// 						  const Eigen::RowVector3d &w, 
// 						  const Eigen::RowVector3d &axis,
// 						  Eigen::Matrix<double,9,1>* dadv,
// 						  Eigen::Matrix<double,9,9>* d2adv2)
// {
//     const double theta = 2.0 * atan2((v.cross(w).dot(axis) / axis.norm()), v.dot(w) + v.norm() * w.norm());

// 	if (dadv != nullptr) {
// 		dadv->segment(0,3) = -axis.cross(v) / v.squaredNorm() / axis.norm();
// 		dadv->segment(3,3) = axis.cross(w) / w.squaredNorm() / axis.norm();
// 		dadv->segment(6,3).setZero();
// 	}
// 	if (d2adv2 != nullptr) {
// 		d2adv2->setZero();
// 		d2adv2->block(0, 0, 3, 3) += 2.0 * axis.cross(v).transpose() * v / v.squaredNorm() / v.squaredNorm() / axis.norm();
// 		d2adv2->block(3, 3, 3, 3) += -2.0 * axis.cross(w).transpose() * w / w.squaredNorm() / w.squaredNorm() / axis.norm();
// 		d2adv2->block(0, 0, 3, 3) += -crossMatrix(axis) / v.squaredNorm() / axis.norm();
// 		d2adv2->block(3, 3, 3, 3) += crossMatrix(axis) / w.squaredNorm() / axis.norm();

// 		const Eigen::Matrix3d dahat = (Eigen::Matrix3d::Identity() / axis.norm() - axis.transpose() * axis / axis.norm() / axis.norm() / axis.norm());

// 		d2adv2->block(0, 6, 3, 3) += crossMatrix(v) * dahat / v.squaredNorm();
// 		d2adv2->block(3, 6, 3, 3) += -crossMatrix(w) * dahat / w.squaredNorm();
// 	}
//     return theta;
// }

// // double GeometryLib::area(const Eigen::RowVector3d &p0, const Eigen::RowVector3d &p1, const Eigen::RowVector3d &p2)
// // {
// //     return 0.5 * (p1 - p0).cross(p2 - p0).norm();
// // }

// double GeometryLib::height(const Eigen::RowVector3d &x0, 
// 						   const Eigen::RowVector3d &x1, 
// 						   const Eigen::RowVector3d &x2,
// 						   Eigen::Matrix<double,9,1>* dhdx)
// {
// 	const Eigen::RowVector3d n = (x1 - x0).cross(x2 - x0);
// 	const double AA = n.norm();
// 	const Eigen::RowVector3d ev = x2 - x1;
// 	const double e = ev.norm();

// 	if (dhdx != nullptr) {
// 		dhdx->setZero();
// 		Eigen::Matrix<double,3,9> dndx;
// 		dndx << crossMatrix(x2-x1), crossMatrix(x0-x2), crossMatrix(x1-x0);
// 		for (int i = 0; i < 3; ++i) {
// 			*dhdx += dndx.row(i).transpose() * n[i] / AA / e;
// 		}
// 		dhdx->segment(6, 3) += -AA / e / e / e * ev.transpose();
// 		dhdx->segment(3, 3) += AA / e / e / e * ev.transpose();
// 	}

//     return AA / e;
// }

// Eigen::Vector3d GeometryLib::shapeOperator(const std::array<Eigen::RowVector3d, 6> &ps,  
// 										   const std::array<Eigen::RowVector2d, 3> &ts,
// 										   Eigen::Matrix<double,3,18>* dkdx)
// {
// 	const Eigen::RowVector3d& p1 = ps[0];
// 	const Eigen::RowVector3d& p2 = ps[1];
// 	const Eigen::RowVector3d& p3 = ps[2];
// 	const Eigen::RowVector3d& p4 = ps[3];
// 	const Eigen::RowVector3d& p5 = ps[4];
// 	const Eigen::RowVector3d& p6 = ps[5];

// 	constexpr std::array<std::array<int,4>,3> flaps = {{{1,2,0,3}, {2,0,1,4}, {0,1,2,5}}};
// 	Eigen::Matrix2d S; S.setZero();
// 	if (dkdx != nullptr) {
// 		dkdx->setZero();
// 	}
// 	for (size_t i = 0; i < 3; ++i) {
// 		const std::array<int,4>& fi = flaps[i];
// 		const Eigen::RowVector2d& ti = ts[i];
// 		const Eigen::RowVector3d& x0 = ps[fi[0]];
// 		const Eigen::RowVector3d& x1 = ps[fi[1]];
// 		const Eigen::RowVector3d& x2 = ps[fi[2]];
// 		const Eigen::RowVector3d& x3 = ps[fi[3]];

// 		Eigen::Matrix<double,12,1> dthetadx;
// 		const double thetai = ti.isZero() ? 0.0 : exteriorDihedralAngle(x0,x1,x2,x3, (dkdx!=nullptr) ? &dthetadx : nullptr);

// 		Eigen::Matrix<double,9,1> dhdx;
// 		const double hi = height(x2, x0, x1, (dkdx!=nullptr) ? &dhdx : nullptr);
// 		// double hi = 1;
		
// 		Eigen::Matrix2d tmp = 1.0 / hi * ti.transpose() * ti;
// 		S += thetai * tmp;

// 		if (dkdx != nullptr) {
// 			const Eigen::Vector3d dkdtheta(tmp(0,0), tmp(1,1), 2*tmp(0,1));
// 			for (int j = 0; j < 4; ++j) {
// 				dkdx->block(0,3*fi[j],3,3) += dkdtheta * dthetadx.segment(3*j,3).transpose();
// 			}
// 			tmp *= -thetai / hi;
// 			const Eigen::Vector3d dkdh(tmp(0,0), tmp(1,1), 2*tmp(0,1));
// 			for (int j = 0; j < 3; ++j) {
// 				const int xid = (fi[j] + 2) % 3;
// 				dkdx->block(0,3*xid,3,3) += dkdh * dhdx.segment(3*j,3).transpose();
// 			}
// 		}
// 	}
	
// 	// const Eigen::RowVector2d& t1 = ts[0];
// 	// const Eigen::RowVector2d& t2 = ts[1];
// 	// const Eigen::RowVector2d& t3 = ts[2];

// 	// Eigen::Matrix<double,12,1> dt1dx, dt2dx, dt3dx;
// 	// const double theta1 = t1.isZero() ? 0.0 : exteriorDihedralAngle(p2,p3,p1,p4);
//     // const double theta2 = t2.isZero() ? 0.0 : exteriorDihedralAngle(p3,p1,p2,p5);
//     // const double theta3 = t3.isZero() ? 0.0 : exteriorDihedralAngle(p1,p2,p3,p6);

//     // const double h1 = height(p1, p2, p3);
//     // const double h2 = height(p2, p3, p1);
//     // const double h3 = height(p3, p1, p2);

//     // const Eigen::Matrix2d S = theta1 / h1 * t1.transpose() * t1 + 
// 	// 						  theta2 / h2 * t2.transpose() * t2 + 
// 	// 						  theta3 / h3 * t3.transpose() * t3;

//     const Eigen::Vector3d S_voigt(S(0,0), S(1,1), 2.0 * S(0,1));

    
//     return S_voigt;
// }
