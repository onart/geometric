// PUBLIC DOMAIN.
#include "voronoi.h"
#include <memory>
#include <cmath>
#include <list>
#include <algorithm>
#include <queue>
#include <map>
#include <set>

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

namespace geom{

        inline bool operator<(const vec2& lhs, const vec2& rhs) { return lhs.y < rhs.y; }

        inline real distance2(const vec2& a, const vec2& b){
            return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y);
        }

        inline real distance(const vec2& a, const vec2& b){
            return std::sqrt(distance2(a, b));
        }

        inline real ccwDirection(const vec2& site, const vec2& vertex){
            return std::atan2(vertex.y - site.y, vertex.x - site.x);
        }

        inline bool outside(const vec2& point) {
            return point.x > real(1.0) || point.x < real(0.0) || point.y > real(1.0) || point.y < real(0.0);
        }

        struct index_t{
            real ccwdir;
            uint32_t index;
            inline bool operator<(const index_t& rhs) { return ccwdir < rhs.ccwdir; }
        };

        struct beachline{
            struct ray{
                vec2 origin;
                vec2 direction;
                inline real x(real y){ // y: current sweepline pos
                    real dy = y - origin.y;
                    return origin.x + dy * direction.x;
                }
                inline real intersectT(const ray& other){
                    if(direction.x == real(0.0)) {
                        if(other.direction.x == real(0.0)) return real(-1.0);
                        real rooty = (origin.x - other.origin.x) * other.direction.x;
                        return (rooty - origin.y) * direction.y; // diretion.y always 1 or -1
                    }
                    // origin (d,e), direction (1,f) ray: y - e = f(x - d) (x >= d)
                    // origin (d,e), direction (-1,-f) ray: y - e = f(x - d) (x <= d)
                    // d = r.origin.x, e = r.origin.y,  f = r.direction.y * r.direction.x
                    real d1 = origin.x, d2 = other.origin.x;
                    real e1 = origin.y, e2 = other.origin.y;
                    real f1 = direction.x * direction.y, f2 = other.direction.x * other.direction.y;
                    if(f1 == f2) return real(-1.0);
                    real rootx = (e2 - e1 + f1 * d1 - f2 * d2) / (f1 - f2);
                    return (rootx - origin.x) * direction.x;
                }
                inline vec2 intersectionWithBoundary(){
                    real t(FLT_MAX);
                    if(direction.x){
                        real tx0 = -origin.x / direction.x;
                        if(tx0 >= real(0.0)) t = tx0;
                        real tx1 = (1 - origin.x) / direction.x;
                        if(tx1 >= real(0.0)) t = std::min(t, tx1);
                    }
                    if(direction.y){
                        real ty0 = -origin.y / direction.y;
                        if(ty0 >= real(0.0)) t = std::min(t, ty0);
                        real ty1 = (1 - origin.y) / direction.y;
                        if(ty1 >= real(0.0)) t = std::min(t, ty1);
                    }
                    return vec2{origin.x + direction.x * t, origin.y + direction.y * t};
                }

                inline void intersectionsWithBoundary(vec2* points) {
                    if (direction.x) {
                        real t = -origin.x / direction.x;
                        if (t >= real(0.0)) {
                            *points = { origin.x + direction.x * t, origin.y + direction.y * t };
                            if (!outside(*points)) { points++; }
                        }
                        t = (1 - origin.x) / direction.x;
                        if (t >= real(0.0)) {
                            *points = { origin.x + direction.x * t, origin.y + direction.y * t };
                            if (!outside(*points)) { points++; }
                        }
                    }
                    if (direction.y) {
                        real t = -origin.y / direction.y;
                        if (t >= real(0.0)) {
                            *points = { origin.x + direction.x * t, origin.y + direction.y * t };
                            if (!outside(*points)) { points++; }
                        }
                        t = (1 - origin.y) / direction.y;
                        if (t >= real(0.0)) {
                            *points = { origin.x + direction.x * t, origin.y + direction.y * t };
                        }
                    }
                }
            };
            struct parabola {
                vec2 focus;
                uint32_t focusIndex;
                mutable std::shared_ptr<ray> eleft, eright;
                mutable real y;
                mutable real lowx, highx;
                inline bool operator<(const parabola& p) const {
                    p.updatey(std::max(y, p.y));
                    updatey(std::max(y, p.y));
                    return (lowx < p.lowx) || (lowx == p.lowx && highx < p.highx);
                }
                private:
                inline void updatey(real s) const {
                    if(y == s) return;
                    lowx = eleft ? beachline::findx(focus, s, *eleft) : real(-FLT_MAX);
                    highx = eright ? beachline::findx(focus, s, *eright) : real(FLT_MAX);
                    y = s;
                }
            };

            const std::vector<vec2>& sites;
            const real& ycontext;

            inline auto insertSite(uint32_t id){
                parabola newp;
                newp.focusIndex = id;
                newp.focus = sites[id];
                newp.y = ycontext;
                newp.lowx = newp.highx = newp.focus.x;
                if(fronts.empty()) { 
                    return fronts.insert(newp).first;
                }
                auto it = fronts.upper_bound(newp);
                it = std::prev(it);
                vec2 intersection;
                // focus (a, b), directrix (y = c) parabola: (x-a)^2 = (b-c)(2y-b-c)
                // a = it->focus.x, b = it->focus.y, c = newp.focus.y
                intersection.x = newp.focus.x;
                intersection.y =
                    ((newp.focus.x - it->focus.x) * (newp.focus.x - it->focus.x)
                    / (it->focus.y - newp.focus.y)
                    + newp.focus.y + it->focus.y) * real(0.5);
                // slope = (x - a) / (b - c)
                real slope = (intersection.x - it->focus.x) / (it->focus.y - newp.focus.y);
                newp.eleft = std::make_shared<ray>();
                newp.eright = std::make_shared<ray>();
                newp.eleft->origin = newp.eright->origin = intersection;
                if(it->focus.y != newp.focus.y){
                    newp.eleft->direction = vec2{real(-1.0), -slope};
                    newp.eright->direction = vec2{real(1.0), slope};
                }
                else{ // (almost)equal site y's -> slope vertical
                    newp.eleft->direction = vec2{0, -1};
                    newp.eright->direction = vec2{0, -1};
                }
                parabola splitp(*it);
                it->highx = newp.lowx;
                std::shared_ptr<ray> oldRight = it->eright;
                it->eright = newp.eleft;
                splitp.lowx = newp.highx;
                splitp.eleft = newp.eright;
                splitp.eright = oldRight;
                fronts.insert(splitp);
                return fronts.insert(newp).first;
            }

            inline void eraseArc(std::set<parabola>::iterator& it){
                fronts.erase(it);
            }

            inline static ray verticalBisector(const vec2& p1, const vec2& p2, const vec2& origin){
                // downwards: y >= 0
                ray ret;
                ret.origin = origin;
                ret.direction.x = p1.y-p2.y;
                ret.direction.y = p2.x-p1.x;
                if(ret.direction.x == real(0.0)){
                    ret.direction.y = real(1.0);
                }
                else{
                    ret.direction.y /= ret.direction.x;
                    ret.direction.x = real(1.0);
                    if(ret.direction.y < 0){
                        ret.direction.x = -ret.direction.x;
                        ret.direction.y = -ret.direction.y;
                    }
                }
                return ret;
            }

            inline beachline(const std::vector<vec2>& sites, const real& y):sites(sites),ycontext(y){}

            std::set<parabola> fronts;
        private:
            inline static real findx(const vec2& focus, real directrixY, const ray& r) {
                if(r.direction.x == real(0.0)) { return r.origin.x; } // vertical
                // else r.direction.x is +1 or -1, meaning x range is greater or less than origin x respectively

                // focus (a, b), directrix (y = c) parabola: (x-a)^2 = (b-c)(2y-b-c)
                // a = focus.x, b = focus.y, c = directrixY
                // origin (d,e), direction (1,f) ray: y - e = f(x - d) (x >= d)
                // origin (d,e), direction (-1,-f) ray: y - e = f(x - d) (x <= d)
                // d = r.origin.x, e = r.origin.y,  f = r.direction.y * r.direction.x

                // intersection eq: ((x-a)^2 / (b-c) + b + c) / 2 = fx - fd + e
                real slope = r.direction.y * r.direction.x;
                real coef2 = real(1.0) / (real(2.0) * (focus.y - directrixY));
                real coef1 = -real(2.0) * focus.x / (real(2.0) * (focus.y - directrixY)) - slope;
                real constant = focus.x * focus.x / (real(2.0) * (focus.y - directrixY)) + (focus.y + directrixY) * real(0.5) + slope * r.origin.x - r.origin.y;
                real rootp1 = -coef1 / (real(2.0) * coef2);
                real rootp2 = std::abs(std::sqrt((coef1 * coef1) - real(4.0) * coef2 * constant) / (real(2.0) * coef2)); // voronoi: there would always be root(s) -> determinant always non negative
                //rootp1 + rootp2 or rootp1 - rootp2
                return rootp1 + rootp2 * r.direction.x;
            }
        };

        struct event{
            enum { SITE = 0, INTERSECT = 1 } type;
            union{
                int32_t site;
                struct{
                    std::set<beachline::parabola>::iterator intersection;
                    vec2 vert;
                    std::weak_ptr<beachline::ray> leftEdge, rightEdge;
                };
            };
            real t;
            inline bool operator<(const event& r) const { return t < r.t; }
            inline bool operator>(const event& r) const { return t > r.t; }
            inline event& operator=(const event& r) { std::memcpy(this, &r, sizeof(event));  return *this; }
            inline event(){};
            inline event(const event& e){
                std::memcpy(this,&e,sizeof(event));
            }
            inline event(decltype(intersection)& it) : intersection(it), leftEdge(it->eleft), rightEdge(it->eright) {}
            inline ~event(){}
        };

    void VoronoiDiagram::calculate(){
        // for when reused
        areaV.clear();
        areaI.clear();
        
        std::vector<vec2> verts;
        std::vector<std::vector<index_t>> areas(sites.size()); // for areaV's locality in GPU
        areaI.resize(sites.size());

        verts.reserve(sites.size());

        real y(0.0);

        std::priority_queue<event, std::vector<event>, std::greater<event>> evs;
        beachline diag(sites, y);

        uint32_t _00CornerOwner = ~0U;
        uint32_t _01CornerOwner = ~0U;
        uint32_t _10CornerOwner = ~0U;
        uint32_t _11CornerOwner = ~0U;
        real _00dist2 = FLT_MAX, _01dist2 = FLT_MAX, _10dist2 = FLT_MAX, _11dist2 = FLT_MAX;

        for(uint32_t i=0;i<sites.size();i++){
            event ev;
            ev.type = event::SITE;
            ev.site = i;
            ev.t = sites[i].y;
            evs.push(ev);
            // insert 4 corners
            real _00 = distance2(sites[i], vec2{0,0});
            real _01 = distance2(sites[i], vec2{0,1});
            real _10 = distance2(sites[i], vec2{1,0});
            real _11 = distance2(sites[i], vec2{1,1});
            if(_00dist2 > _00){ _00dist2 = _00; _00CornerOwner = i; }
            if(_01dist2 > _01){ _01dist2 = _01; _01CornerOwner = i; }
            if(_10dist2 > _10){ _10dist2 = _10; _10CornerOwner = i; }
            if(_11dist2 > _11){ _11dist2 = _11; _11CornerOwner = i; }
        }

        verts.push_back(vec2{0,0});
        verts.push_back(vec2{0,1});
        verts.push_back(vec2{1,0});
        verts.push_back(vec2{1,1});

        areas[_00CornerOwner].push_back({ccwDirection(sites[_00CornerOwner], vec2{0,0}), 0});
        areas[_01CornerOwner].push_back({ccwDirection(sites[_01CornerOwner], vec2{0,1}), 1});
        areas[_10CornerOwner].push_back({ccwDirection(sites[_10CornerOwner], vec2{1,0}), 2});
        areas[_11CornerOwner].push_back({ccwDirection(sites[_11CornerOwner], vec2{1,1}), 3});

        // process start
        while(!evs.empty()){
            auto ev = evs.top(); evs.pop(); // top can be changed when pushing (same y) -> pop immediately
            y = ev.t;
            if(ev.type == event::SITE){
                auto it = diag.insertSite(ev.site);
                if (it != diag.fronts.begin()) { // size == 1 or same y
                    auto prev = std::prev(it);
                    if (prev->eleft) {
                        real intersectT = it->eleft->intersectT(*prev->eleft);
                        if (intersectT >= 0) {
                            vec2 vert = {
                                it->eleft->origin.x + it->eleft->direction.x * intersectT,
                                it->eleft->origin.y + it->eleft->direction.y * intersectT
                            };
                            event iev(prev);
                            iev.type = event::INTERSECT;
                            iev.t = vert.y + distance(vert, it->focus);
                            iev.vert = vert;
                            evs.push(iev);
                        }
                    }
                }
                auto next = std::next(it);
                if (next == diag.fronts.end()) continue;
                if(next->eright){
                    real intersectT = it->eright->intersectT(*next->eright);
                    if(intersectT >= 0) {
                        vec2 vert = {
                            it->eright->origin.x + it->eright->direction.x * intersectT,
                            it->eright->origin.y + it->eright->direction.y * intersectT
                        };
                        event iev(next);
                        iev.type = event::INTERSECT;
                        iev.t = vert.y + distance(vert,it->focus);
                        iev.vert = vert;
                        evs.push(iev);
                    }
                }
            }
            else{
                if (ev.leftEdge.expired() || ev.rightEdge.expired()) {
                    continue;
                }
                auto prev = std::prev(ev.intersection);
                auto next = std::next(ev.intersection);
                if (prev->eright.get() != ev.leftEdge.lock().get() || next->eleft.get() != ev.rightEdge.lock().get()) {
                    continue;
                }
                // 3 arcs met -> the vertex is always inside the [0,1]^2 area
                if (outside(prev->eright->origin)) {
                    beachline::ray outward;
                    outward.origin = ev.vert;
                    outward.direction.x = -prev->eright->direction.x;
                    outward.direction.y = -prev->eright->direction.y;

                    vec2 borderVert = outward.intersectionWithBoundary();
                    uint32_t index = verts.size();
                    areas[prev->focusIndex].push_back({ ccwDirection(prev->focus, ev.vert), index });
                    areas[ev.intersection->focusIndex].push_back({ ccwDirection(ev.intersection->focus, ev.vert), index });
                    verts.push_back(borderVert);
                }

                if (outside(next->eleft->origin)) {
                    beachline::ray outward;
                    outward.origin = ev.vert;
                    outward.direction.x = -next->eleft->direction.x;
                    outward.direction.y = -next->eleft->direction.y;

                    vec2 borderVert = outward.intersectionWithBoundary();
                    uint32_t index = verts.size();
                    areas[next->focusIndex].push_back({ ccwDirection(next->focus, ev.vert), index });
                    areas[ev.intersection->focusIndex].push_back({ ccwDirection(ev.intersection->focus, ev.vert), index });
                    verts.push_back(borderVert);
                }
                beachline::ray newEdge = beachline::verticalBisector(prev->focus, next->focus, ev.vert);
                next->eleft = prev->eright = std::make_shared<beachline::ray>(newEdge);
                uint32_t index = verts.size();
                verts.push_back(ev.vert);
                areas[ev.intersection->focusIndex].push_back({ ccwDirection(ev.intersection->focus, ev.vert), index });
                areas[prev->focusIndex].push_back({ ccwDirection(prev->focus, ev.vert), index });
                areas[next->focusIndex].push_back({ ccwDirection(next->focus, ev.vert), index });
                diag.eraseArc(ev.intersection);
            }
        }
        // process remainings in beachline
        for(auto it = diag.fronts.begin(); it != diag.fronts.end();){
            if(auto& ry = it->eright){
                vec2 candids[4]{ {real(-1)},{real(-1)},{real(-1)},{real(-1)} };
                ry->intersectionsWithBoundary(candids);
                auto next = std::next(it);
                auto owner1 = it->focusIndex;
                auto owner2 = next->focusIndex;
                for (int i = 0; i < 4; i++) {
                    vec2& vertex = candids[i];
                    if (outside(vertex)) break;
                    uint32_t index = verts.size();
                    verts.push_back(vertex);
                    areas[owner1].push_back({ ccwDirection(it->focus, vertex), index });
                    areas[owner2].push_back({ ccwDirection(next->focus, vertex), index });
                }
                it = next;
            }
            else{
                ++it; // only for the last one
            }
        }

        // vertex locality
        areaV.reserve(verts.size());
        std::vector<uint32_t> indexMap(verts.size());
        for(auto& i: indexMap) { i = ~0u; }
        int indexIter = 0;
        int faceIter = 0;
        for(auto& area: areas){
            std::sort(area.begin(), area.end());
            for(const index_t& vpos: area){
                if(indexMap[vpos.index] == ~0u){
                    areaV.push_back(verts[vpos.index]);
                    areaI[faceIter].push_back(indexIter);
                    indexMap[vpos.index] = indexIter++;
                }
                else {
                    areaI[faceIter].push_back(indexMap[vpos.index]);
                }
            }
            faceIter++;
        }
    }

    void VoronoiDiagram::makeVertexPerFace() {
        std::vector<vec2> newV;
        std::vector<std::vector<uint32_t>> newI(areaI.size());
        newV.reserve(areaV.size() * 3);
        int face = 0;
        uint32_t vtx = 0;
        for(auto& area: areaI) {
            for(uint32_t idx: area){
                newV.push_back(areaV[area[idx]]);
                newI[face].push_back(vtx++);
            }
        }
        areaV.swap(newV);
        areaI.swap(newI);
    }
}
