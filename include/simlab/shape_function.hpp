#ifndef SIMLAB_SHAPE_FUNCTION_HPP
#define SIMLAB_SHAPE_FUNCTION_HPP

#include <Eigen/Dense>

namespace simlab
{

    struct ShapeFunctionData
    {
        Eigen::VectorXd N;
        Eigen::VectorXd dNdx;
        Eigen::VectorXd dNdy;
        double detJ;
    };

    namespace shape_function
    {

        ShapeFunctionData quad_4_isoparametric(const Eigen::VectorXd &point, const Eigen::MatrixXd &nodes);

    } // namespace shape_function

} // namespace simlab

#endif // !SIMLAB_SHAPE_FUNCTION_HPP