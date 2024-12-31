#ifndef SIMLAB_BOUNDARY_CONDITION_HPP
#define SIMLAB_BOUNDARY_CONDITION_HPP

#include <Eigen/Dense>

namespace simlab
{

    struct SupportBC
    {
        Eigen::Index node;
        bool directions[3];
    };

    struct PointLoadBC
    {
        Eigen::Index node;
        double value = 0.0;
        bool directions[3];
    };

} // namespace simlab

#endif // !SIMLAB_BOUNDARY_CONDITION_HPP