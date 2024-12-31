#ifndef SIMLAB_MESH_HPP
#define SIMLAB_MESH_HPP

#include <map>
#include <set>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "simlab/boundary_condition.hpp"
#include "simlab/element.hpp"

namespace simlab
{

    struct Mesh
    {
        std::vector<Element> elements;
        std::vector<SupportBC> supports;
        std::vector<PointLoadBC> point_loads;

        Eigen::MatrixXi connectivity;

        void initialize();
        void apply_boundary_conditions();

        Eigen::SparseMatrix<double> stiffness() const;
        Eigen::VectorXd force() const;
        Eigen::VectorXd rhs() const;
        Eigen::VectorXd u(const Eigen::VectorXd &solution) const;

    private:
        Eigen::Index m_num_nodes = 0;
        Eigen::Index m_num_dofs = 0;
        Eigen::Index m_num_elements = 0;

        Eigen::VectorXi active_dofs;

        std::set<Eigen::Index> fixed_dofs;

        std::map<Eigen::Index, std::vector<Eigen::Index>> node_to_dof_map;
        std::map<Eigen::Index, std::vector<Eigen::Index>> element_to_dof_map;
        std::map<Eigen::Index, Eigen::Index> active_dof_map;
    };

} // namespace simlab

#endif // !SIMLAB_MESH_HPP