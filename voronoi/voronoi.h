// PUBLIC DOMAIN.
#include <vector>
#include <cstdint>

namespace onart{
    using real = float;
    struct vec2{
        real x,y;
        inline bool operator<(const vec2& other){return y < other.y;}
    };

    // base board range: [0,1]^2 (left-top side (0,0))
    struct VoronoiDiagram{
        std::vector<vec2> sites; // input

        std::vector<vec2> areaV; // intersection (Vertices)
        std::vector<std::vector<uint32_t>> areaE; // cell (set of Edges)
        void calculate();
    };
}
