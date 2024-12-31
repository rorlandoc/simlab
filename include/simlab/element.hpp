#ifndef SIMLAB_ELEMENT_HPP
#define SIMLAB_ELEMENT_HPP

#include <Eigen/Dense>

#include "simlab/material.hpp"

namespace simlab
{

    struct Element
    {
        enum class Type
        {
            QUAD_4_PS,
        };

        Type type;
        Eigen::MatrixXd nodal_coordinates;
        MaterialData material;

        Element(Element::Type type, const Eigen::MatrixXd &coords, const MaterialData &material);

        Eigen::Index count_nodes() const;
        Eigen::Index count_dofs_per_node() const;
        Eigen::Index count_dofs() const;

        static Eigen::MatrixXd stiffness(const Element &element);
    };

} // namespace simlab

#endif // !SIMLAB_ELEMENT_HPP