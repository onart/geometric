// PUBLIC DOMAIN.
// PROGRESSING...

#include "voronoi.h"
#include <cmath>
#include <list>
#include <algorithm>
#include <queue>
#include <map>

/*
Fortune's algorithm
[https://jacquesheunis.com/post/fortunes-algorithm/]
포춘 알고리즘은 스윕라인 접근을 기반으로 함
각 기준점(site)을 순서대로 고려하고 기준점을 중심으로 한 '칸(cell)'을 확장
칸이 다른 모든 칸에 완전히 둘러싸인 경우 그 이상 확장될 수 없음
스윕라인을 따라 아직 확장 중인 칸(beachline)만 고려하면 되는 것

확장 중인 칸은 포물선(기준점을 초점, 스윕라인을 준선으로)들에 둘러싸인 형태
여기서는 스윕라인이 위에서 아래로 향하는 형태로

스윕라인의 이동에 따라 생각할 때 실제로 각 객체의 관계가 변경되는 때(event)는
- 스윕라인이 기준점을 만나는 때(site event)
- 2개의 칸이 성장하다 사이의 다른 칸을 덮는 때(edge-intersection event)
알고리즘은 site event들을 순회하면서 시작됨

효율을 위해 beachline의 칸은 x좌표 순으로 정렬됨 (초점의 x좌표 순이 아님. 후술)

포물선끼리 만나는 점은 기준점과의 거리가 각각 준선과의 거리로 동일한 것이므로 칸 사이 경계의 궤적이 됨

site event
- 새로운 칸이 beachline에 추가됨
- 먼저 어떤 기존의 포물선이 새로운 포물선 아래에 끼게 되는지 확인.
    초점과 준선이 일치한 상태일 때 포물선은 해당 초점을 지나는 수선.
    해당 수선과 만나는 기존 포물선은 둘로 나뉘고 새로운 포물선은 그 사이에 삽입
    기존 포물선 아래의 칸과 이 기준점에 대한 칸은 인접, 사이의 간선이 존재
    포물선을 나누는, 확장 중인 간선과 완성된 간선이 구분됨
    확장 중인 간선은 2개의 반직선 비슷한 것으로 시작점(방금 수선과의 교점)으로부터 일정한 방향으로 뻗음
    해당 방향은 물론 교점에서의 접선
    이 2개의 반직선을 순서대로 beachline에는
    {기존 포물선} 제거,
    {잘린 왼쪽 / 반직선 왼쪽 / 새 포물선 / 반직선 오른쪽 / 잘린 오른쪽} 추가
    이때 잘린 왼쪽 / 잘린 오른쪽이 이미 edge-intersection이 되었는지 확인할 필요가 있음
    - 일단 바로 왼쪽/오른쪽 반직선의 성장을 검사 (자체가 없으면 이벤트 추가는 없는 것)
    반직선이 방향 그대로 확장하여 만나는 곳이 3개 기준점에 대한 포물선의 교점
    이 이벤트가 발생하는 스윕라인 위치는 해당 교점 y좌표에 한 초점과의 거리를 더한 값
    (+y 방향이 아래, 즉 스윕 진행 방향인 경우.)

edge-intersection event
- 없어지는 포물선과 만난 반직선(즉 이제는 선분)은 beachline에서 제거
    선분이 된 반직선은 끝점이 생기며 이 끝점은 그것이 속한 칸의 꼭짓점에 해당. 그것을 기록
    없어지는 포물선에 인접한 포물선들을 찾아, 사이에 이어지는 간선을 생성
    간선의 방향은 두 포물선의 초점을 이은 선분에 수직하고 아래를 향하는 방향이어야 함
    추가 후 두 포물선에 대하여 새 간선 교점 이벤트를 테스트

모든 이벤트가 끝난 경우 아직 완전히 끝나지 않은 간선을 처리해야 함
영역의 끝을 만날 때까지 확장
*/

namespace onart{

    void VoronoiDiagram::calculate(){
        // sturctures only used in here
        struct event{
            enum { SITE = 0, INTERSECT = 1 } type;
            union{
                int32_t site;
            };
            real t;
            inline bool operator<(const event& r){ return t < r.t; }
            inline event(){};
        };

        struct beachline{
            struct ray{
                vec2 origin;
                vec2 direction;
                inline real x(real y){ // y: current sweepline pos
                    real dy = y - origin.y;
                    return origin.x + dy * direction.x;
                }
            };
            std::vector<uint32_t> parabolaFocuses;
            std::vector<ray> rays;
            // edge 들어갈 자리를 리턴
            inline int bisect(const vec2& xy){
                if(rays.empty()){ return 0; }
                return bisect(xy, 0, rays.size() - 1);
            }
        private:
            inline int bisect(const vec2& xy, int begin, int end){
                if(begin == end) { return begin; }
                int mid = begin + end / 2;
                if(xy.x > rays[mid].x(xy.y)){
                    return bisect(xy, mid + 1, end);
                }
                else{ // TODO: ==인 경우 0으로 나누기가 생김
                    return bisect(xy, begin, mid - 1);
                }
            }
        };

        // for when reused
        areaV.clear();
        areaE.clear();
        
        std::vector<vec2> verts;
        std::vector<std::vector<uint32_t>> areas(sites.size()); // for areaV's locality in GPU
        areaE.resize(sites.size());

        std::priority_queue<event> evs;
        beachline diag;

        for(uint32_t i=0;i<sites.size();i++){
            event ev;
            ev.type = event::SITE;
            ev.site = i;
            ev.t = sites[i].y;
            evs.push(ev);
        }

        // process start
        while(!evs.empty()){
            auto ev = evs.top(); evs.pop();
            if(ev.type == event::SITE){
                // new parabola
                vec2 site = sites[ev.site];
                int pos = diag.bisect(site);
                vec2 originalFocus = sites[diag.parabolaFocuses[pos]];
                // focus (a, b), directrix (y = c) parabola: (x-a)^2 = 4(b-c)(y-c)
                // a = originalFocus.x, b = originalFocus.y, c = site.y
                vec2 intersection;
                intersection.x = site.x;
                intersection.y = 
                    (site.x - originalFocus.x) * (site.x - originalFocus.x) 
                    / (real(4.0) * (originalFocus.y - site.y)) 
                    + site.y;
                // slope = (x - a) / (2(b-c))
                real slope = real(2.0) * (intersection.x - originalFocus.x) / (real(4.0) * (originalFocus.y - site.y));
                beachline::ray ray0, ray1;
                ray0.origin = ray1.origin = intersection;
                if(originalFocus.y != site.y){
                    ray0.direction = vec2{-1.0f, slope};
                    ray1.direction = vec2{1.0f, slope};
                }
                else{ // (almost)equal site y's -> slope vertical

                }
                diag.rays.insert(diag.rays.cbegin() + pos, {});
            }
            else{
                
            }
        }
    }
}
