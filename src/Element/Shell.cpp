#include <Element/Shell.h>

#include <igl/triangle_triangle_adjacency.h>
#include <igl/doublearea.h>

Shell::Shell(const Eigen::MatrixX3i& F) : m_F(F)
{
    Eigen::MatrixXi TT, TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);
    m_flaps.resize(F.rows());
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            const int eid = (j + 1) % 3;
            const int fid = TT(i, eid);
            m_flaps[i][j] = fid < 0 ? -1 : F(fid, (TTi(i, eid) + 2) % 3);
        }
    }
}

Eigen::VectorXd Shell::areas(const Eigen::MatrixX2d& X) const
{
    Eigen::VectorXd AA;
    igl::doublearea(X, m_F, AA);
    return 0.5 * AA;
}

Eigen::Matrix<double, 3, 2> Shell::deformationGradient(
    int i,
    const Eigen::MatrixX2d& X,
    const Eigen::MatrixX3d& x, 
    Eigen::Matrix<double, 6, 9>* dFdxi) const
{
    const Eigen::Matrix3d xi = x(m_F.row(i), Eigen::placeholders::all);
    const Eigen::Matrix<double, 3, 2> Xi = X(m_F.row(i), Eigen::placeholders::all);

	// F = dx/duv = dx/dq dq/dX dX/duv
	// dX/duv = I for orthogonal materials

	const Eigen::Vector3d v1 = (xi.row(1) - xi.row(0)).transpose();
	const Eigen::Vector3d v2 = (xi.row(2) - xi.row(0)).transpose();
	const Eigen::Vector2d V1 = (Xi.row(1) - Xi.row(0)).transpose();
	const Eigen::Vector2d V2 = (Xi.row(2) - Xi.row(0)).transpose();

	Eigen::Matrix<double, 3, 2> dxdq;
    dxdq << v1, v2;
	Eigen::Matrix2d dXdq;
	dXdq << V1, V2;
    
	const Eigen::Matrix2d dqdX = dXdq.inverse();
    const Eigen::Matrix<double, 3, 2> F = dxdq * dqdX;

	if (dFdxi != nullptr) { // F in colume major
		Eigen::Matrix3d I; I.setIdentity();
		*dFdxi << -dqdX.col(0).sum() * I, dqdX(0, 0) * I, dqdX(1, 0) * I,
				  -dqdX.col(1).sum() * I, dqdX(0, 1) * I, dqdX(1, 1) * I;
	}

    return F;
}

Eigen::Vector3d Shell::GreenStrainVoigt(
    int i,
    const Eigen::MatrixX2d& X,
    const Eigen::MatrixX3d& x, 
    Eigen::Matrix<double, 3, 9>* dEdxi,
    std::array<Eigen::Matrix<double, 9, 9>, 3>* d2Edxi2) const
{
    Eigen::Matrix<double, 6, 9> dFdxi;
    const Eigen::Matrix<double, 3, 2> F = deformationGradient(i, X, x, (dEdxi!=nullptr || d2Edxi2!=nullptr) ? &dFdxi : nullptr);

    const Eigen::Matrix2d E = 0.5 * (F.transpose() * F - Eigen::Matrix2d::Identity());
    const Eigen::Vector3d E_voigt(E(0, 0), E(1, 1), 2.0 * E(0, 1));

    if (dEdxi != nullptr) {
        Eigen::Matrix<double, 3, 6> dEdF;
        dEdF << F.col(0).transpose(), 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, F.col(1).transpose(),
                F.col(1).transpose(), F.col(0).transpose();
        *dEdxi = dEdF * dFdxi;
    }

    if (d2Edxi2 != nullptr) {
        std::array<Eigen::Matrix<double, 6, 6>, 3> d2EdF2;
        d2EdF2[0].setZero(); d2EdF2[0](0, 0) = d2EdF2[0](1, 1) = d2EdF2[0](2, 2) = 1.0;
        d2EdF2[1].setZero(); d2EdF2[1](3, 3) = d2EdF2[1](4, 4) = d2EdF2[1](5, 5) = 1.0;
        d2EdF2[2].setZero(); d2EdF2[2](0, 3) = d2EdF2[2](3, 0) = d2EdF2[2](1, 4) = d2EdF2[2](4, 1) = d2EdF2[2](2, 5) = d2EdF2[2](5, 2) = 1.0;
        for (size_t j = 0; j < 3; ++j) {
            d2Edxi2->at(j) = dFdxi.transpose() * d2EdF2[j] * dFdxi;
        }
        
    }
    return E_voigt;
}

Eigen::Matrix3d Shell::crossMatrix(const Eigen::RowVector3d& v)
{
	Eigen::Matrix3d ret;
	ret << 0.0, -v[2], v[1],
		   v[2], 0.0, -v[0],
		   -v[1], v[0], 0.0;
	return ret;
}

double Shell::height(
    const Eigen::RowVector3d &x0, 
    const Eigen::RowVector3d &x1, 
    const Eigen::RowVector3d &x2,
    Eigen::Matrix<double, 9, 1>* dhdx,
    Eigen::Matrix<double, 9, 9>* d2hdx2)
{
	const Eigen::RowVector3d n = (x1 - x0).cross(x2 - x0);
	const double AA = n.norm();
	const Eigen::RowVector3d ev = x2 - x1;
	const double e = ev.norm();

	if (dhdx != nullptr) {
		dhdx->setZero();
		Eigen::Matrix<double, 3, 9> dndx;
		dndx << crossMatrix(x2-x1), crossMatrix(x0-x2), crossMatrix(x1-x0);
		for (int i = 0; i < 3; ++i) {
			*dhdx += dndx.row(i).transpose() * n[i] / AA / e;
		}
		dhdx->segment(6, 3) += -AA / e / e / e * ev.transpose();
		dhdx->segment(3, 3) += AA / e / e / e * ev.transpose();
	}
    if (d2hdx2 != nullptr) {
        d2hdx2->setZero();
        Eigen::Matrix<double, 3, 9> dndx;
		dndx << crossMatrix(x2-x1), crossMatrix(x0-x2), crossMatrix(x1-x0);
        for (int i = 0; i < 3; ++i) {
            Eigen::Matrix<double, 9, 9> d2nidx2;
            d2nidx2.setZero();
            Eigen::Vector3d ei(0.0, 0.0, 0.0);
            ei[i] = 1.0;
            const Eigen::Matrix3d eic = crossMatrix(ei);
            d2nidx2.block(0, 3, 3, 3) -= eic;
            d2nidx2.block(0, 6, 3, 3) += eic;
            d2nidx2.block(3, 6, 3, 3) -= eic;
            d2nidx2.block(3, 0, 3, 3) += eic;
            d2nidx2.block(6, 0, 3, 3) -= eic;
            d2nidx2.block(6, 3, 3, 3) += eic;

            *d2hdx2 += d2nidx2 * n(i) / AA / e;
        }
        const Eigen::Matrix3d P = Eigen::Matrix3d::Identity() / AA - n.transpose() * n / AA / AA / AA;
        *d2hdx2 += dndx.transpose() * P * dndx / e;
        d2hdx2->block(6, 0, 3, 9) += -ev.transpose() * n * dndx / AA / e / e / e;
        d2hdx2->block(3, 0, 3, 9) += ev.transpose() * n * dndx / AA / e / e / e;
        d2hdx2->block(0, 6, 9, 3) += -dndx.transpose() * n.transpose() * ev / AA / e / e / e;
        d2hdx2->block(0, 3, 9, 3) += dndx.transpose() * n.transpose() * ev / AA / e / e / e;
        d2hdx2->block(6, 6, 3, 3) += -AA / e / e / e * Eigen::Matrix3d::Identity();
        d2hdx2->block(6, 3, 3, 3) += AA / e / e / e * Eigen::Matrix3d::Identity();
        d2hdx2->block(3, 6, 3, 3) += AA / e / e / e * Eigen::Matrix3d::Identity();
        d2hdx2->block(3, 3, 3, 3) += -AA / e / e / e * Eigen::Matrix3d::Identity();
        const Eigen::Matrix3d outer = ev.transpose() * ev * 3.0 * AA / e / e / e / e / e;
        d2hdx2->block(6, 6, 3, 3) += outer;
        d2hdx2->block(6, 3, 3, 3) += -outer;
        d2hdx2->block(3, 6, 3, 3) += -outer;
        d2hdx2->block(3, 3, 3, 3) += outer;
    }

    return AA / e;
}

double Shell::angle(
    const Eigen::RowVector3d &v, 
    const Eigen::RowVector3d &w, 
    const Eigen::RowVector3d &axis,
    Eigen::Matrix<double,9,1>* dadv,
    Eigen::Matrix<double,9,9>* d2adv2)
{
    const double theta = 2.0 * atan2(v.cross(w).dot(axis) / axis.norm(), v.dot(w) + v.norm() * w.norm());

	if (dadv != nullptr) {
		dadv->segment(0, 3) = -axis.cross(v) / v.squaredNorm() / axis.norm();
		dadv->segment(3, 3) = axis.cross(w) / w.squaredNorm() / axis.norm();
		dadv->segment(6, 3).setZero();
	}
	if (d2adv2 != nullptr) {
		d2adv2->setZero();
		d2adv2->block(0, 0, 3, 3) += 2.0 * axis.cross(v).transpose() * v / v.squaredNorm() / v.squaredNorm() / axis.norm();
		d2adv2->block(3, 3, 3, 3) += -2.0 * axis.cross(w).transpose() * w / w.squaredNorm() / w.squaredNorm() / axis.norm();
		d2adv2->block(0, 0, 3, 3) += -crossMatrix(axis) / v.squaredNorm() / axis.norm();
		d2adv2->block(3, 3, 3, 3) += crossMatrix(axis) / w.squaredNorm() / axis.norm();

		const Eigen::Matrix3d dahat = (Eigen::Matrix3d::Identity() / axis.norm() - axis.transpose() * axis / axis.norm() / axis.norm() / axis.norm());

		d2adv2->block(0, 6, 3, 3) += crossMatrix(v) * dahat / v.squaredNorm();
		d2adv2->block(3, 6, 3, 3) += -crossMatrix(w) * dahat / w.squaredNorm();
	}
    return theta;
}

double Shell::exteriorDihedralAngle(
    const Eigen::RowVector3d &x0, 
    const Eigen::RowVector3d &x1, 
    const Eigen::RowVector3d &x2, 
    const Eigen::RowVector3d &x3,
    Eigen::Matrix<double,12,1>* dtdx,
    Eigen::Matrix<double,12,12>* d2tdx2)
{
    const Eigen::RowVector3d e0 = x1 - x0;

    const Eigen::RowVector3d n0 = (x0 - x2).cross(x1 - x2);
    const Eigen::RowVector3d n1 = (x1 - x3).cross(x0 - x3);

	Eigen::Matrix<double, 9, 1> dadv;
	Eigen::Matrix<double, 9, 9> d2adv2;
    const double theta = angle(n0, n1, e0, (dtdx!=nullptr || d2tdx2!=nullptr) ? &dadv : nullptr, d2tdx2!=nullptr ? &d2adv2 : nullptr);

	if (dtdx != nullptr) {
		dtdx->setZero();
		dtdx->segment(0, 3) += crossMatrix(x1 - x2) * dadv.segment(0, 3);
		dtdx->segment(3, 3) += crossMatrix(x2 - x0) * dadv.segment(0, 3);
		dtdx->segment(6, 3) += crossMatrix(x0 - x1) * dadv.segment(0, 3);

		dtdx->segment(0, 3) += crossMatrix(x3 - x1) * dadv.segment(3, 3);
		dtdx->segment(3, 3) += crossMatrix(x0 - x3) * dadv.segment(3, 3);
		dtdx->segment(9, 3) += crossMatrix(x1 - x0) * dadv.segment(3, 3);
	}

	if (d2tdx2 != nullptr) {
		d2tdx2->setZero();
		const std::array<Eigen::Matrix3d, 3> vqm = {crossMatrix(x2 - x1), crossMatrix(x0 - x2), crossMatrix(x1 - x0)}; // dn0dx
		const std::array<Eigen::Matrix3d, 3> wqm = {crossMatrix(x1 - x3), crossMatrix(x3 - x0), crossMatrix(x0 - x1)}; // dn1dx

		const std::array<int, 3> vindices = { 0, 3, 6 };
		const std::array<int, 3> windices = { 0, 3, 9 };

		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				d2tdx2->block(vindices[i], vindices[j], 3, 3) += vqm[i].transpose() * d2adv2.block(0, 0, 3, 3) * vqm[j];
				d2tdx2->block(vindices[i], windices[j], 3, 3) += vqm[i].transpose() * d2adv2.block(0, 3, 3, 3) * wqm[j];
				d2tdx2->block(windices[i], vindices[j], 3, 3) += wqm[i].transpose() * d2adv2.block(3, 0, 3, 3) * vqm[j];
				d2tdx2->block(windices[i], windices[j], 3, 3) += wqm[i].transpose() * d2adv2.block(3, 3, 3, 3) * wqm[j];
			}

			d2tdx2->block(vindices[i], 3, 3, 3) += vqm[i].transpose() * d2adv2.block(0, 6, 3, 3);
			d2tdx2->block(3, vindices[i], 3, 3) += d2adv2.block(6, 0, 3, 3) * vqm[i];
			d2tdx2->block(vindices[i], 0, 3, 3) += -vqm[i].transpose() * d2adv2.block(0, 6, 3, 3);
			d2tdx2->block(0, vindices[i], 3, 3) += -d2adv2.block(6, 0, 3, 3) * vqm[i];

			d2tdx2->block(windices[i], 3, 3, 3) += wqm[i].transpose() * d2adv2.block(3, 6, 3, 3);
			d2tdx2->block(3, windices[i], 3, 3) += d2adv2.block(6, 3, 3, 3) * wqm[i];
			d2tdx2->block(windices[i], 0, 3, 3) += -wqm[i].transpose() * d2adv2.block(3, 6, 3, 3);
			d2tdx2->block(0, windices[i], 3, 3) += -d2adv2.block(6, 3, 3, 3) * wqm[i];

		}

		const Eigen::RowVector3d dang1 = dadv.segment(0, 3).transpose();
		const Eigen::RowVector3d dang2 = dadv.segment(3, 3).transpose();

		const Eigen::Matrix3d dang1mat = crossMatrix(dang1);
		const Eigen::Matrix3d dang2mat = crossMatrix(dang2);

		d2tdx2->block(6, 3, 3, 3) += dang1mat;
		d2tdx2->block(0, 3, 3, 3) -= dang1mat;
		d2tdx2->block(0, 6, 3, 3) += dang1mat;
		d2tdx2->block(3, 0, 3, 3) += dang1mat;
		d2tdx2->block(3, 6, 3, 3) -= dang1mat;
		d2tdx2->block(6, 0, 3, 3) -= dang1mat;

		d2tdx2->block(9, 0, 3, 3) += dang2mat;
		d2tdx2->block(3, 0, 3, 3) -= dang2mat;
		d2tdx2->block(3, 9, 3, 3) += dang2mat;
		d2tdx2->block(0, 3, 3, 3) += dang2mat;
		d2tdx2->block(0, 9, 3, 3) -= dang2mat;
		d2tdx2->block(9, 3, 3, 3) -= dang2mat;
	}

	return theta;
}

Eigen::Vector3d Shell::shapeOperator(
    int idx,
    const Eigen::MatrixX2d& X,
    const Eigen::MatrixX3d& x, 
    Eigen::Matrix<double, 3, 18>* dkdx,
    std::array<Eigen::Matrix<double, 18, 18>, 3>* d2kdx2) const
{
	Eigen::Vector3d S; S.setZero();
	if (dkdx != nullptr) {
		dkdx->setZero();
	}
    if (d2kdx2 != nullptr) {
        for (size_t i = 0; i < 3; ++i) {
            d2kdx2->at(i).setZero();
        }
    }

    const Eigen::RowVector2d Xi0 = X.row(m_F(idx, 0));
    const Eigen::RowVector2d Xi1 = X.row(m_F(idx, 1));
    const Eigen::RowVector2d Xi2 = X.row(m_F(idx, 2));
    const std::array<Eigen::RowVector2d, 3> vs = {Xi2 - Xi1, Xi0 - Xi2, Xi1 - Xi0};

	for (int i = 0; i < 3; ++i) {
        if (m_flaps[idx][i] < 0) {
            continue;
        }
		const std::array<int, 4> av = {(i + 1) % 3, (i + 2) % 3, i, i + 3};
        const Eigen::Vector2d ei = X.row(m_F(idx, av[1])) - X.row(m_F(idx, av[0]));
        const Eigen::Vector2d ti = Eigen::Vector2d(ei.y(), -ei.x()).normalized();
		const Eigen::RowVector3d& x0 = x.row(m_F(idx, av[0]));
		const Eigen::RowVector3d& x1 = x.row(m_F(idx, av[1]));
		const Eigen::RowVector3d& x2 = x.row(m_F(idx, av[2]));
		const Eigen::RowVector3d& x3 = x.row(m_flaps[idx][i]);

		Eigen::Matrix<double, 12, 1> dthetadx;
		Eigen::Matrix<double, 12, 12> d2thetadx2;
		const double thetai = ti.isZero() ? 0.0 : exteriorDihedralAngle(x0, x1, x2, x3, 
                                                (dkdx!=nullptr || d2kdx2!=nullptr) ? &dthetadx : nullptr, 
                                                d2kdx2!=nullptr ? &d2thetadx2 : nullptr);

        const std::array<int, 3> hv = {av[2], av[0], av[1]};
		Eigen::Matrix<double, 9, 1> dhdx;
        Eigen::Matrix<double, 9, 9> d2hdx2;
		const double hi = height(x2, x0, x1, (dkdx!=nullptr || d2kdx2!=nullptr) ? &dhdx : nullptr, d2kdx2!=nullptr ? &d2hdx2 : nullptr);
		
        // voigt notation
		Eigen::Vector3d tmp(ti(0) * ti(0), ti(1) * ti(1), 2.0 * ti(0) * ti(1));
        tmp /= hi;
		S += thetai * tmp;

		if (dkdx != nullptr) {
			const Eigen::Vector3d& dkdtheta = tmp;
			for (int j = 0; j < 4; ++j) {
				dkdx->block(0, 3 * av[j], 3, 3) += dkdtheta * dthetadx.segment(3 * j, 3).transpose();
			}
			const Eigen::Vector3d dkdh = -thetai / hi * tmp;
			for (int j = 0; j < 3; ++j) {
				dkdx->block(0, 3 * hv[j], 3, 3) += dkdh * dhdx.segment(3 * j, 3).transpose();
			}
		}
        if (d2kdx2 != nullptr) {
			const Eigen::Vector3d& dkdtheta = tmp;
            for (size_t d = 0; d < 3; ++d) {
                for (size_t k = 0; k < 4; ++k) {
                    for (size_t j = 0; j < 4; ++j) {
                        d2kdx2->at(d).block(3 * av[j], 3 * av[k], 3, 3) += dkdtheta(d) * d2thetadx2.block(3 * j, 3 * k, 3, 3);
                    }
                }
            }
            const Eigen::Vector3d dkdh = -thetai / hi * tmp;
            const Eigen::Vector3d d2kdh2 = 2.0 * thetai / hi / hi * tmp;
            for (size_t d = 0; d < 3; ++d) {
                for (size_t k = 0; k < 3; ++k) {
                    for (size_t j = 0; j < 3; ++j) {
                        d2kdx2->at(d).block(3 * hv[j], 3 * hv[k], 3, 3) += dkdh(d) * d2hdx2.block(3 * j, 3 * k, 3, 3);
                        d2kdx2->at(d).block(3 * hv[j], 3 * hv[k], 3, 3) += d2kdh2(d) * dhdx.segment(3 * j, 3) * dhdx.segment(3 * k, 3).transpose();
                    }
                }
            }
            const Eigen::Vector3d d2kdthetadh = -tmp / hi;
            for (size_t d = 0; d < 3; ++d) {
                for (size_t k = 0; k < 3; ++k) {
                    for (size_t j = 0; j < 4; ++j) {
                        d2kdx2->at(d).block(3 * av[j], 3 * hv[k], 3, 3) += d2kdthetadh(d) * dthetadx.segment(3 * j, 3) * dhdx.segment(3 * k, 3).transpose();
                        d2kdx2->at(d).block(3 * hv[k], 3 * av[j], 3, 3) += d2kdthetadh(d) * dhdx.segment(3 * k, 3) * dthetadx.segment(3 * j, 3).transpose();
                    }
                }
            }
        }
	}
	    
    return S;
}

std::array<double, 2> Shell::maxStrain(const Eigen::MatrixX2d& X, const Eigen::MatrixX3d& x) const
{
    std::array<double, 2> strain = {0.0, 0.0};
    for (int i = 0; i < m_F.rows(); ++i) {
        const Eigen::Matrix<double, 3, 2> F = deformationGradient(i, X, x);
        const Eigen::Matrix2d E = 0.5 * (F.transpose() * F - Eigen::Matrix2d::Identity());
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(E);
        const Eigen::Vector2d ev = es.eigenvalues();
        strain[0] = std::min(strain[0], ev(0));
        strain[1] = std::max(strain[1], ev(1));
    }
    return strain;
}

std::array<double, 2> Shell::maxCurvature(const Eigen::MatrixX2d& X, const Eigen::MatrixX3d& x) const
{
    std::array<double, 2> curvature = {0.0, 0.0};
    for (int i = 0; i < m_F.rows(); ++i) {
        const Eigen::Vector3d Sv = shapeOperator(i, X, x);
        Eigen::Matrix2d S;
        S << Sv(0), 0.5 * Sv(2), 0.5 * Sv(2), Sv(1);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(S);
        const Eigen::Vector2d ev = es.eigenvalues();
        curvature[0] = std::min(curvature[0], ev(0));
        curvature[1] = std::max(curvature[1], ev(1));
    }
    return curvature;
}

