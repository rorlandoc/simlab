#include "simlab/integration.hpp"

namespace simlab
{

    namespace integration
    {

        IntegrationData gauss_2D_1x1()
        {
            IntegrationData data;
            data.weights.resize(1);
            data.points.resize(1, 2);
            data.weights << 4.0;
            data.points << 0.0, 0.0;
            return data;
        }

        IntegrationData gauss_2D_2x2()
        {
            IntegrationData data;
            data.weights.resize(4);
            data.points.resize(4, 2);
            data.weights << 1.0, 1.0, 1.0, 1.0;
            data.points << -0.5773502691896257, -0.5773502691896257,
                0.5773502691896257, -0.5773502691896257,
                0.5773502691896257, 0.5773502691896257,
                -0.5773502691896257, 0.5773502691896257;
            return data;
        }

    } // namespace integration

} // namespace simlab