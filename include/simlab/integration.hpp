#ifndef SIMLAB_INTEGRATION_HPP
#define SIMLAB_INTEGRATION_HPP

#include <Eigen/Dense>

namespace simlab
{

    struct IntegrationData
    {
        Eigen::VectorXd weights;
        Eigen::MatrixXd points;
    };

    namespace integration
    {

        IntegrationData gauss_2D_1x1();
        IntegrationData gauss_2D_2x2();

    } // namespace integration

} // namespace simlab

#endif // !SIMLAB_INTEGRATION_HPP