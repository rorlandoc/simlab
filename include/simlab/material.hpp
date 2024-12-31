#ifndef SIMLAB_MATERIAL_HPP
#define SIMLAB_MATERIAL_HPP

#include <Eigen/Dense>

namespace simlab
{

    struct MaterialData
    {
        Eigen::VectorXd properties;
    };

    namespace material
    {

        Eigen::MatrixXd elastic_isotropic_plane_stress(const Eigen::VectorXd &properties);

    } // namespace material

} // namespace simlab

#endif // !SIMLAB_MATERIAL_HPP