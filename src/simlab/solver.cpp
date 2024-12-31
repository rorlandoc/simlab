#include "simlab/solver.hpp"
#include "simlab/mesh.hpp"

namespace simlab
{

    DirectSolver::DirectSolver(std::shared_ptr<Mesh> mesh)
        : m_mesh(mesh)
    {
    }

    std::optional<Eigen::VectorXd> DirectSolver::solve()
    {
        Eigen::SparseMatrix<double> K = m_mesh->stiffness();
        Eigen::VectorXd f = m_mesh->rhs();

        m_solver.compute(K);
        if (m_solver.info() != Eigen::Success)
        {
            return std::nullopt;
        }

        Eigen::VectorXd result = m_solver.solve(f);
        if (m_solver.info() != Eigen::Success)
        {
            return std::nullopt;
        }

        return result;
    }

    IterativeSolver::IterativeSolver(std::shared_ptr<Mesh> mesh)
        : m_mesh(mesh)
    {
    }

    std::optional<Eigen::VectorXd> IterativeSolver::solve()
    {
        Eigen::SparseMatrix<double> K = m_mesh->stiffness();
        Eigen::VectorXd f = m_mesh->force();

        m_solver.compute(K);
        if (m_solver.info() != Eigen::Success)
        {
            return std::nullopt;
        }

        Eigen::VectorXd result = m_solver.solve(f);
        if (m_solver.info() != Eigen::Success)
        {
            return std::nullopt;
        }

        return result;
    }

} // namespace simlab
