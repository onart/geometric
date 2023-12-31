#ifndef __OA_VORONOI_H__
#define __OA_VORONOI_H__
#include <vector>
#include <cstdint>

namespace geom{
    using real = float;
    struct vec2{
        real x,y;
    };

    // 보드의 범위는 [0,1]^2 입니다. (좌측 상단이 (0,0))
    struct VoronoiDiagram{
        std::vector<vec2> sites; // input

        std::vector<vec2> areaV; // intersection (Vertices)
        std::vector<std::vector<uint32_t>> areaI; // cell (set of Indices per face, ccw)

        /// @brief calculate voronoi diagram for current sites. no duplicate vertices
        void calculate();
        /// @brief duplicate vertices so that each faces have very their vertices
        void makeVertexPerFace();
    };
}
#endif
