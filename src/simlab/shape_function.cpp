#include "simlab/shape_function.hpp"

namespace simlab
{

    namespace shape_function
    {

        ShapeFunctionData quad_4_isoparametric(const Eigen::VectorXd &point, const Eigen::MatrixXd &nodes)
        {
            ShapeFunctionData data;

            data.N.resize(4);
            data.dNdx.resize(4);
            data.dNdy.resize(4);

            double xi = point(0);
            double eta = point(1);

            data.N << 0.25 * (1 - xi) * (1 - eta),
                0.25 * (1 + xi) * (1 - eta),
                0.25 * (1 + xi) * (1 + eta),
                0.25 * (1 - xi) * (1 + eta);

            Eigen::MatrixXd dN(2, 4);
            dN.row(0) << -0.25 * (1 - eta),
                0.25 * (1 - eta),
                0.25 * (1 + eta),
                -0.25 * (1 + eta);
            dN.row(1) << -0.25 * (1 - xi),
                -0.25 * (1 + xi),
                0.25 * (1 + xi),
                0.25 * (1 - xi);

            Eigen::MatrixXd J = dN * nodes;
            data.detJ = J.determinant();

            dN = J.inverse() * dN;
            data.dNdx = dN.row(0);
            data.dNdy = dN.row(1);

            return data;
        }

    } // namespace shape_function

} // namespace simlab