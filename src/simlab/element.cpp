#include "simlab/element.hpp"
#include "simlab/integration.hpp"
#include "simlab/shape_function.hpp"

#include <map>

namespace simlab
{

    namespace internal
    {
        std::map<Element::Type, Eigen::Index> count_dofs_per_node = {
            {Element::Type::QUAD_4_PS, 2}};

        std::map<Element::Type, Eigen::Index> count_nodes = {
            {Element::Type::QUAD_4_PS, 4}};

        Eigen::MatrixXd quad_4_ps(const Element &element)
        {
            IntegrationData integration = integration::gauss_2D_2x2();
            Eigen::MatrixXd D = material::elastic_isotropic_plane_stress(element.material.properties);

            int num_points = integration.weights.size();

            Eigen::MatrixXd K = Eigen::MatrixXd::Zero(8, 8);
            for (int ipt = 0; ipt < num_points; ipt++)
            {
                Eigen::VectorXd point = integration.points.row(ipt);
                double weight = integration.weights(ipt);

                ShapeFunctionData shape = shape_function::quad_4_isoparametric(point, element.nodal_coordinates);

                Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8);
                for (int i = 0; i < 4; i++)
                {
                    B(0, 2 * i) = shape.dNdx(i);
                    B(1, 2 * i + 1) = shape.dNdy(i);
                    B(2, 2 * i) = shape.dNdy(i);
                    B(2, 2 * i + 1) = shape.dNdx(i);
                }

                K += weight * B.transpose() * D * B * shape.detJ;
            }

            return K;
        }

    } // namespace internal

    Element::Element(Element::Type type, const Eigen::MatrixXd &coords, const MaterialData &material)
        : type(type), nodal_coordinates(coords), material(material)
    {
        assert(nodal_coordinates.rows() == count_nodes());
        assert(nodal_coordinates.cols() == 2);
        assert(material.properties.size() == 2);
    }

    Eigen::Index Element::count_nodes() const
    {
        return internal::count_nodes[type];
    }

    Eigen::Index Element::count_dofs_per_node() const
    {
        return internal::count_dofs_per_node[type];
    }

    Eigen::Index Element::count_dofs() const
    {
        return count_nodes() * count_dofs_per_node();
    }

    Eigen::MatrixXd Element::stiffness(const Element &element)
    {
        switch (element.type)
        {
        case Element::Type::QUAD_4_PS:
            return internal::quad_4_ps(element);
        default:
            return Eigen::MatrixXd{};
        }
    }

} // namespace simlab