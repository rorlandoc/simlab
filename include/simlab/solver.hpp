#ifndef SIMLAB_SOLVER_HPP
#define SIMLAB_SOLVER_HPP

#include <memory>
#include <optional>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

namespace simlab
{
    struct Mesh;

    class DirectSolver
    {
    public:
        DirectSolver(std::shared_ptr<Mesh> mesh);

        std::optional<Eigen::VectorXd> solve();

    private:
        Eigen::SparseLU<Eigen::SparseMatrix<double>> m_solver;

        std::shared_ptr<Mesh> m_mesh;
    };

    class IterativeSolver
    {
    public:
        IterativeSolver(std::shared_ptr<Mesh> mesh);

        std::optional<Eigen::VectorXd> solve();

    private:
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> m_solver;

        std::shared_ptr<Mesh> m_mesh;
    };

} // namespace simlab

#endif // !SIMLAB_SOLVER_HPP