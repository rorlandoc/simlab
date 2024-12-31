#include "simlab/mesh.hpp"

#include <map>
#include <set>

namespace simlab
{

    void Mesh::initialize()
    {
        m_num_elements = elements.size();

        std::set<Eigen::Index> nodes;
        std::map<Eigen::Index, std::set<Eigen::Index>> node_to_element;
        std::map<Eigen::Index, Eigen::Index> node_to_dof_count;

        for (Eigen::Index irow = 0; irow < connectivity.rows(); irow++)
        {
            Eigen::Index dofs_per_node = elements[irow].count_dofs_per_node();
            for (Eigen::Index icol = 0; icol < connectivity.cols(); icol++)
            {
                auto node_id = connectivity(irow, icol);
                nodes.insert(node_id);
                node_to_element[node_id].insert(irow);
                node_to_dof_count[node_id] = std::max(node_to_dof_count[node_id], dofs_per_node);
            }
        }
        m_num_nodes = nodes.size();

        for (auto &node : nodes)
        {
            m_num_dofs += node_to_dof_count[node];
        }

        Eigen::Index dof_counter = 0;
        for (auto &node : nodes)
        {
            node_to_dof_map[node] = std::vector<Eigen::Index>(node_to_dof_count[node]);
            for (Eigen::Index i = 0; i < node_to_dof_count[node]; i++)
            {
                node_to_dof_map[node][i] = dof_counter++;
            }
        }

        for (Eigen::Index irow = 0; irow < connectivity.rows(); irow++)
        {
            auto element_dofs = element_to_dof_map[irow];
            element_dofs = std::vector<Eigen::Index>{};
            for (Eigen::Index icol = 0; icol < connectivity.cols(); icol++)
            {
                auto node_id = connectivity(irow, icol);
                auto dofs = node_to_dof_map[node_id];
                element_dofs.insert(element_dofs.end(), dofs.begin(), dofs.end());
            }
            element_to_dof_map[irow] = element_dofs;
        }

        active_dofs = Eigen::VectorXi::LinSpaced(m_num_dofs, 0, m_num_dofs - 1);
    }

    void Mesh::apply_boundary_conditions()
    {
        Eigen::VectorXi temp(m_num_dofs);
        Eigen::Index counter = 0;
        for (auto &support : supports)
        {
            auto dofs = node_to_dof_map[support.node];
            for (Eigen::Index i = 0; i < dofs.size(); i++)
            {
                if (support.directions[i])
                {
                    fixed_dofs.insert(dofs[i]);
                }
            }
        }
        for (Eigen::Index i = 0; i < m_num_dofs; i++)
        {
            if (!fixed_dofs.contains(i))
            {
                temp(counter++) = i;
            }
        }
        active_dofs = temp.head(counter);

        for (Eigen::Index i = 0; i < active_dofs.size(); i++)
        {
            active_dof_map[active_dofs[i]] = i;
        }
    }

    Eigen::SparseMatrix<double> Mesh::stiffness() const
    {
        std::vector<Eigen::Triplet<double>> triplets;
        std::map<std::pair<Eigen::Index, Eigen::Index>, Eigen::Index> cache;

        for (Eigen::Index i = 0; i < m_num_elements; i++)
        {
            const auto &element = elements[i];
            const auto &dofs = element_to_dof_map.at(i);
            Eigen::MatrixXd Ke = Element::stiffness(element);

            for (Eigen::Index i = 0; i < dofs.size(); i++)
            {
                if (fixed_dofs.contains(dofs[i]))
                {
                    continue;
                }
                for (Eigen::Index j = 0; j < dofs.size(); j++)
                {
                    if (fixed_dofs.contains(dofs[j]))
                    {
                        continue;
                    }
                    Eigen::Index idof = active_dof_map.at(dofs[i]);
                    Eigen::Index jdof = active_dof_map.at(dofs[j]);
                    auto ids = std::make_pair(idof, jdof);
                    if (cache.contains(ids))
                    {
                        auto triplet = triplets[cache.at(ids)];
                        triplets[cache.at(ids)] = Eigen::Triplet<double>(idof, jdof, triplet.value() + Ke(i, j));
                    }
                    else
                    {
                        triplets.push_back(Eigen::Triplet<double>(idof, jdof, Ke(i, j)));
                        cache[ids] = triplets.size() - 1;
                    }
                }
            }
        }
        Eigen::SparseMatrix<double> K(active_dofs.size(), active_dofs.size());
        K.setFromTriplets(triplets.begin(), triplets.end());
        return K;
    }

    Eigen::VectorXd Mesh::force() const
    {
        Eigen::VectorXd f = Eigen::VectorXd::Zero(m_num_dofs);
        for (auto &load : point_loads)
        {
            auto dofs = node_to_dof_map.at(load.node);
            for (Eigen::Index i = 0; i < dofs.size(); i++)
            {
                if (load.directions[i])
                {
                    f(dofs[i]) = load.value;
                }
            }
        }
        return f;
    }

    Eigen::VectorXd Mesh::rhs() const
    {
        return force()(active_dofs);
    }

    Eigen::VectorXd Mesh::u(const Eigen::VectorXd &solution) const
    {
        Eigen::VectorXd u = Eigen::VectorXd::Zero(m_num_dofs);
        u(active_dofs) = solution;
        return u;
    }

} // namespace simlab