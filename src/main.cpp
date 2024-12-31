#include <iostream>
#include <memory>
#include <string_view>

#include "simlab/boundary_condition.hpp"
#include "simlab/element.hpp"
#include "simlab/mesh.hpp"
#include "simlab/solver.hpp"

template <typename T>
void print(std::string_view name, const T &m);
#define print_var(x) ::print(#x, x)

int main()
{
    auto mesh = std::make_shared<simlab::Mesh>();
    simlab::MaterialData material;

    Eigen::MatrixXi connectivity(2, 4);
    Eigen::MatrixXd nodes(6, 2);

    connectivity.row(0) << 0, 1, 2, 3;
    connectivity.row(1) << 1, 4, 5, 2;

    nodes.row(0) << 0, 0;
    nodes.row(1) << 1, 0;
    nodes.row(2) << 1, 1;
    nodes.row(3) << 0, 1;
    nodes.row(4) << 2, 0;
    nodes.row(5) << 2, 1;

    material.properties = Eigen::VectorXd{{1.0, 0.3}};

    for (Eigen::Index iel = 0; iel < connectivity.rows(); iel++)
    {
        simlab::Element element(simlab::Element::Type::QUAD_4_PS,
                                nodes(connectivity.row(iel), Eigen::all),
                                material);
        mesh->elements.push_back(element);
    }
    mesh->connectivity = connectivity;
    mesh->initialize();

    mesh->supports.push_back(simlab::SupportBC{0, {true, true, false}});  // fixed in x and y
    mesh->supports.push_back(simlab::SupportBC{1, {false, true, false}}); // free in x and fixed in y
    mesh->supports.push_back(simlab::SupportBC{4, {false, true, false}}); // free in x and fixed in y
    mesh->supports.push_back(simlab::SupportBC{3, {true, false, false}}); // fixed in x and free in y

    mesh->point_loads.push_back(simlab::PointLoadBC{4, 0.5, {true, false, false}}); // 0.5 in x
    mesh->point_loads.push_back(simlab::PointLoadBC{5, 0.5, {true, false, false}}); // 0.5 in x

    mesh->apply_boundary_conditions();

    simlab::DirectSolver direct_solver(mesh);

    auto opt_result = direct_solver.solve();
    if (opt_result.has_value())
    {
        print_var(mesh->u(opt_result.value()));
    }
    else
    {
        std::cout << "Solver failed" << std::endl;
    }
}

template <typename T>
void print(std::string_view name, const T &m)
{
    std::cout << name << ":\n"
              << m << std::endl;
}