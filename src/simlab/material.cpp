#include "simlab/material.hpp"

namespace simlab
{

    namespace material
    {

        Eigen::MatrixXd elastic_isotropic_plane_stress(const Eigen::VectorXd &properties)
        {
            Eigen::MatrixXd D(3, 3);
            D.setZero();

            double E = properties(0);
            double nu = properties(1);

            double Dn = E / (1.0 - nu * nu);
            double Ds = Dn * nu;
            double G = E / (2.0 * (1.0 + nu));

            D << Dn, Ds, 0.0,
                Ds, Dn, 0.0,
                0.0, 0.0, G;

            return D;
        }

    } // namespace material

} // namespace simlab