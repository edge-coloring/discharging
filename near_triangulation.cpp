#include <spdlog/spdlog.h>
#include <fmt/ranges.h>
#include "near_triangulation.hpp"

using std::ifstream;

Degree::Degree(int lower_deg, int upper_deg) : 
    lower_deg_(lower_deg), 
    upper_deg_(upper_deg) {};

Degree::Degree(int deg) :
    lower_deg_(deg),
    upper_deg_(deg) {};

// "5", "5+", "8-" などの次数を表す文字列から Degree クラスを返す。
Degree Degree::fromString(const string &str) {
    std::size_t i;
    int deg = std::stoi(str, &i);
    if (i == str.size()) {
        return Degree(deg);
    } 
    if (i == str.size() - 1) {
        if (str[i] == '+') return Degree(deg, MAX_DEGREE);
        if (str[i] == '-') return Degree(MIN_DEGREE, deg);
    }
    spdlog::critical("Failed to parse {} as degree", str);
    throw std::runtime_error("Failed to parse " + str + " as degree");
}

int Degree::lower(void) const {
    return lower_deg_;
}

int Degree::upper(void) const {
    return upper_deg_;
}

string Degree::toString(void) const {
    if (fixed()) return std::to_string(lower_deg_);
    if (upper_deg_ == MAX_DEGREE) return std::to_string(lower_deg_) + "+";
    if (lower_deg_ == MIN_DEGREE) return std::to_string(upper_deg_) + "-";
    assert(false);
}

bool Degree::include(const Degree &degree) const {
    return lower_deg_ <= degree.lower_deg_ && degree.upper_deg_ <= upper_deg_;
}

bool Degree::disjoint(const Degree &degree0, const Degree &degree1) {
    return degree0.upper_deg_ < degree1.lower_deg_ || degree1.upper_deg_ < degree0.lower_deg_;
}

bool Degree::fixed(void) const {
    return lower_deg_ == upper_deg_;
}

Degree Degree::intersection(const Degree &degree0, const Degree &degree1) {
    return Degree(std::max(degree0.lower(), degree1.lower()), std::min(degree0.upper(), degree1.upper()));
}

bool Degree::mergable(const optional<Degree> &deg0, const optional<Degree> &deg1) {
    if (deg0 && deg1) {
        return !Degree::disjoint(deg0.value(), deg1.value());
    } 
    // どちらかが次数がまだ定まっていなければ、 mergable とする。
    return true;
}

bool operator==(const Degree& lhs, const Degree& rhs) {
    return lhs.lower() == rhs.lower() && lhs.upper() == rhs.upper();
}

// degree のとりうる範囲 [l, r] を 5, 6, ..., max_degree+ の次数に分ける。
// 例えば、max_degree = 8 の時
// 5+ -> 5, 6, 7, 8+
// 7+ -> 7, 8+ 
// 6- は 5, 6 のように分ける。
vector<Degree> divideDegree(const Degree &degree, int max_degree) {
    vector<Degree> degrees;
    int deg = degree.lower();
    assert(deg <= max_degree);
    while (deg < degree.upper() && deg < max_degree) {
        degrees.push_back(Degree(deg));
        deg++;
    }
    degrees.push_back(Degree(deg, degree.upper()));
    return degrees;
}

vector<pair<int, int>> enumerateEdges(int vertex_size, const vector<set<int>> &VtoV) {
    vector<pair<int, int>> edges;
    for (int v = 0;v < vertex_size; v++) {
        for (int u : VtoV[v]) {
            edges.emplace_back(v, u);
        }
    }
    return edges;
}

map<pair<int, int>, vector<int>> calcDiagonalVertices(const vector<pair<int, int>> &edges, const vector<set<int>> &VtoV) {
    map<pair<int, int>, vector<int>> diagonal_vertices;
    for (const auto &edge : edges) {
        diagonal_vertices[edge] = {};
        auto [v, u] = edge;
        for (int w : VtoV[v]) {
            if (VtoV[u].count(w)) {
                diagonal_vertices[edge].push_back(w);
            }
        }
    }
    return diagonal_vertices;
}

NearTriangulation::NearTriangulation(int vertex_size, const vector<set<int>> &VtoV, const vector<optional<Degree>> &degrees) : 
    vertex_size_(vertex_size),
    degrees_(degrees),
    VtoV_(VtoV) {

    edges_ = enumerateEdges(vertex_size_, VtoV);
    diagonal_vertices_ = calcDiagonalVertices(edges_, VtoV);
    for (const auto &[e, vs] : diagonal_vertices_) {
        if (vs.size() > 2) {
            spdlog::info("debug(): {}", debug());
        }
        assert(vs.size() <= 2);
    }
}

int NearTriangulation::vertexSize(void) const {
    return vertex_size_;
}

const vector<optional<Degree>> &NearTriangulation::degrees(void) const {
    return degrees_;
}

const vector<pair<int, int>> &NearTriangulation::edges(void) const {
    return edges_;
}

const map<pair<int, int>, vector<int>> &NearTriangulation::diagonalVertices(void) const {
    return diagonal_vertices_;
}

const vector<set<int>> &NearTriangulation::VtoV(void) const {
    return VtoV_;
}

int NearTriangulation::numIncidentEdge(int v) const {
    return (int)VtoV_[v].size();
}

void NearTriangulation::setDegree(int v, const optional<Degree> &degree) {
    degrees_[v] = degree;
    return;
}

// 新しい頂点を付け足す
int NearTriangulation::newVertex(void) {
    int new_vertex = vertex_size_;
    vertex_size_++;
    degrees_.push_back(std::nullopt);
    VtoV_.push_back({});
    return new_vertex;
}

// NearTriangulation に辺 uv を足し、それに合わせて NearTriangulation のデータを変更する。
// 辺 uv を足す操作が valid な操作ではない場合は、false を返す。 valid ではない操作とは以下のような操作のこと。
// 1. u == v
// 2. 元から 辺 uv が存在している。
// 3. 辺 uv を足したときに、その diagonal_vertex のサイズが 2 より大きい。(このとき、cyclical 3-cut ができている。)
// 4. 辺 uv を足したときに、u,v の共通の近傍 w について、辺 uw, vw の diagonal_vertex のサイズが 2 より大きい。(このとき、cyclical 3-cut ができている。)
bool NearTriangulation::addEdge(int u, int v) {
    // 1.
    if (u == v) {
        spdlog::info("near_triangulation: {}", debug());
        assert(false);
    }
    // 2.
    if (VtoV_[u].count(v)) {
        spdlog::info("near_triangulation: {}", debug());
        assert(false);
    }
    // 3.
    vector<int> uv_neighbors;
    for (int v_neighbor : VtoV_[v]) {
        if (VtoV_[u].count(v_neighbor)) uv_neighbors.push_back(v_neighbor);
    }
    if (uv_neighbors.size() > 2) return false;
    // 4.
    for (int w : uv_neighbors) {
        auto eu0 = std::make_pair(u, w);
        auto eu1 = std::make_pair(w, u);
        auto ev0 = std::make_pair(v, w);
        auto ev1 = std::make_pair(w, v);
        if (diagonal_vertices_[eu0].size() >= 2 || diagonal_vertices_[eu1].size() >= 2 || 
            diagonal_vertices_[ev0].size() >= 2 || diagonal_vertices_[ev1].size() >= 2) {
            return false;
        }
    }

    edges_.emplace_back(u, v);
    edges_.emplace_back(v, u);
    VtoV_[u].insert(v);
    VtoV_[v].insert(u);
    auto e0 = std::make_pair(u, v), e1 = std::make_pair(v, u);
    diagonal_vertices_[e0] = uv_neighbors;
    diagonal_vertices_[e1] = uv_neighbors;
    for (int w : uv_neighbors) {
        auto eu0 = std::make_pair(u, w);
        auto eu1 = std::make_pair(w, u);
        diagonal_vertices_[eu0].push_back(v);
        diagonal_vertices_[eu1].push_back(v);
        auto ev0 = std::make_pair(v, w);
        auto ev1 = std::make_pair(w, v);
        diagonal_vertices_[ev0].push_back(u);
        diagonal_vertices_[ev1].push_back(u);
    }
    return true;
}

// NearTriangulation の点 u,v を同一視する。
// 操作: 頂点 v に接続している辺を全て消して、元々頂点 v に接続していた辺を u に繋ぎ直す。
// (注意: 頂点番号を元のままにしておきたいので、頂点 v は孤立点として残す。)
// 点 u, v の同一視をする操作が valid な操作ではない場合は、false を返す。 valid ではない操作とは以下のような操作のこと。
// 1. v == u
// 2. 辺 uv が存在している。
// 3. 頂点 u, v の次数が mergable ではない。(注意: u, v どちらかの次数が決まっていない状態でも mergable になる。)
bool NearTriangulation::identify(int u, int v) {
    // 1.
    if (v == u) {
        assert(false);
    }
    // 2.
    if (VtoV_[v].count(u)) {
        assert(false);
    }
    // 3.
    if (!Degree::mergable(degrees_[u], degrees_[v])) {
        return false;
    }

    Degree any_degree = Degree(MIN_DEGREE, MAX_DEGREE);
    Degree deguv = Degree::intersection(degrees_[v].value_or(any_degree), degrees_[u].value_or(any_degree));
    bool is_deguv_undecided = !degrees_[v].has_value() && !degrees_[u].has_value();

    auto newid = [&](int x) {
        if (x == v) return u;
        else return x;
    };

    vector<set<int>> new_VtoV(vertex_size_);

    for (const auto &e: edges_) {
        new_VtoV[newid(e.first)].insert(newid(e.second));
        new_VtoV[newid(e.second)].insert(newid(e.first));
    }
    vector<pair<int, int>> new_edges = enumerateEdges(vertex_size_, new_VtoV);
    map<pair<int, int>, vector<int>> new_diagonal_vertices = calcDiagonalVertices(new_edges, new_VtoV);

    for (const auto &[e, vs]: new_diagonal_vertices) {
        if (vs.size() > 2) return false;
    }

    // degrees
    degrees_[u] = is_deguv_undecided ? std::nullopt : std::make_optional(deguv);
    // edges
    edges_ = new_edges;
    // diagonal_vertices
    diagonal_vertices_ = new_diagonal_vertices;
    // VtoV
    VtoV_ = new_VtoV;

    return true;
}


// isolated_vertices は孤立点で、それらを除去する。
// vids に対応する新しい頂点番号を返す。
vector<int> NearTriangulation::removeVertices(const vector<int> &isolated_vertices, const vector<int> &vids) {
    // 孤立点であることを確かめる。
    for (int v : isolated_vertices) {
        if (VtoV_[v].size() > 0) {
            spdlog::info("debug: {}", debug());
            spdlog::info("v: {}", v);
            assert(false);
        }
    }

    vector<bool> to_be_removed(vertex_size_, false);
    for (int v: isolated_vertices) {
        to_be_removed[v] = true;
    }

    int new_vertex_size = 0;
    vector<int> new_vertex_id(vertex_size_, -1);
    {
        int c = 0;
        for (int v = 0;v < vertex_size_; v++) {
            if (to_be_removed[v]) {
                continue;
            }
            new_vertex_id[v] = c++;
        }
        new_vertex_size = c;
    }

    vector<optional<Degree>> new_degrees(new_vertex_size, std::nullopt);
    {
        for (int v = 0;v < vertex_size_; v++) {
            if (new_vertex_id[v] == -1) continue;
            new_degrees[new_vertex_id[v]] = degrees_[v];
        }
    }

    vector<set<int>> new_VtoV(new_vertex_size);
    {
        for (const auto &e: edges_) {
            new_VtoV[new_vertex_id[e.first]].insert(new_vertex_id[e.second]);
            new_VtoV[new_vertex_id[e.second]].insert(new_vertex_id[e.first]);
        }
    }

    vector<pair<int, int>> new_edges = enumerateEdges(new_vertex_size, new_VtoV);
    map<pair<int, int>, vector<int>> new_diagonal_vertices = calcDiagonalVertices(new_edges, new_VtoV);

    vertex_size_ = new_vertex_size;
    degrees_ = new_degrees;
    edges_ = new_edges;
    diagonal_vertices_ = new_diagonal_vertices;
    VtoV_ = new_VtoV;

    vector<int> new_vids(vids.size());
    {
        for (int i = 0;i < (int)vids.size(); i++) {
            new_vids[i] = new_vertex_id[vids[i]]; 
        }
    }

    return new_vids;
}


const string CUT3_ERR_MSG = "There may be a forbidden 3-cut in near-triangulation.";

// NearTriangulation の outer_face の頂点を順に並べる。
//「diaognal_vertices のサイズが 1 <-> outer_face の辺である」 が成り立っている必要がある。
// そうでない場合は runtime_error を投げる。
vector<int> NearTriangulation::getOuterFace(void) const {
    if (vertex_size_ == 2 && edges_.size() == 2) {
        // 1 つの辺のみからなる near_triangulation の場合
        vector<int> outer_face = {0, 1};
        return outer_face;
    }
    vector<vector<int>> cycle_neighbors(vertex_size_);
    for (const auto &e : edges_) {
        if (diagonal_vertices_.at(e).size() == 1 && e.first < e.second) {
            cycle_neighbors[e.first].push_back(e.second);
            cycle_neighbors[e.second].push_back(e.first);
        }
        if (diagonal_vertices_.at(e).size() > 2) {
            throw std::runtime_error(CUT3_ERR_MSG);
        }
    }
    int si = -1;
    int count = 0;
    for (int i = 0;i < vertex_size_; i++) {
        if (cycle_neighbors[i].size() != 0 && cycle_neighbors[i].size() != 2) {
            throw std::runtime_error(CUT3_ERR_MSG);
        }
        if (cycle_neighbors[i].size() == 2) {
            si = i;
            count++;
        }
    }
    if (si == -1) {
        throw std::runtime_error(CUT3_ERR_MSG);
    }
    vector<int> outer_face;
    {
        outer_face.push_back(si);
        int i = cycle_neighbors[si][0];
        while (i != si) {
            int new_i = -1;
            if (cycle_neighbors[i][0] == outer_face.back()) {
                new_i = cycle_neighbors[i][1];
            } else if (cycle_neighbors[i][1] == outer_face.back()) {
                new_i = cycle_neighbors[i][0];
            } 
            if (new_i == -1) {
                throw std::runtime_error(CUT3_ERR_MSG);
            }
            outer_face.push_back(i);
            i = new_i;
        }
        if (outer_face.back() != cycle_neighbors[si][1]) {
            throw std::runtime_error(CUT3_ERR_MSG);
        }
    }
    if (count != (int)outer_face.size()) {
        throw std::runtime_error(CUT3_ERR_MSG);
    }
    return outer_face;
}


// 平面グラフであることと nearTriangulation であることを仮定する。
// cyclical 3-cut が存在しないと
// + outer_face に接続している辺の diagonal_vertex のサイズが 1
// + outer_face に接続していない辺の diagonal_vertex のサイズが 2
// が成り立つ。
// これをチェックする。
bool NearTriangulation::exist3cut(void) const {
    try {
        auto outer_face = getOuterFace();
    } catch (std::exception &e) {
        string err_msg = e.what();
        if (err_msg == CUT3_ERR_MSG) return true;
        spdlog::error(err_msg);
        std::exit(1);
    }
    return false;
}

// v とその近傍 neighbor_v から v の近傍の平面埋め込みの時計回り順を得る。
// + v は outer_face に接続している頂点で、 v_neighbor は v の近傍のうち、 outer_face に接続している頂点とする。
// + v_neighbor から順に v の近傍の頂点の平面埋め込みを時計回り(or 反時計回り)順に得る。
vector<int> NearTriangulation::getCyclicOrder(int v, int v_neighbor) const {
    assert(VtoV_[v].count(v_neighbor));
    vector<int> neighbors = {v_neighbor};
    pair<int, int> e = std::make_pair(v, v_neighbor);
    if (diagonal_vertices_.at(e).size() != 1) {
        spdlog::info("debug: {}", debug());
        spdlog::info("e: {}", e);
        assert(false);
    }
    v_neighbor = diagonal_vertices_.at(e)[0];
    while (1) {
        pair<int, int> e = std::make_pair(v, v_neighbor);
        if (diagonal_vertices_.at(e).size() == 1) {
            assert(neighbors.back() == diagonal_vertices_.at(e)[0]);
            neighbors.push_back(v_neighbor);
            break;
        } else {
            assert(diagonal_vertices_.at(e).size() == 2);
            int new_v_neighbor = -1;
            if (neighbors.back() == diagonal_vertices_.at(e)[0]) {
                new_v_neighbor = diagonal_vertices_.at(e)[1];
            } else {
                assert(neighbors.back() == diagonal_vertices_.at(e)[1]);
                new_v_neighbor = diagonal_vertices_.at(e)[0];
            }
            neighbors.push_back(v_neighbor);
            v_neighbor = new_v_neighbor;
        }
    }
    assert(neighbors.size() == VtoV_[v].size());
    return neighbors;
}

// outer_face の頂点が新たに接続する辺の本数を計算する。ただし、
// + まだ次数が決定していない or 次数が範囲指定されており、固定されていない(fixed() = false)ときは INF を返す。
vector<int> NearTriangulation::calcNumNewNeighbors(const vector<int> &outer_face) const {
    vector<int> num_new_neighbors(outer_face.size(), 0);
    for (int i = 0;i < (int)outer_face.size(); i++) {
        int v = outer_face[i];
        if (!degrees_[v].has_value() || !degrees_[v].value().fixed()) {
            num_new_neighbors[i] = INF;
        } else {
            num_new_neighbors[i] = degrees_[v].value().lower() - numIncidentEdge(v);
        }
    }
    return num_new_neighbors;
}


// NearTriangulation::addEdgeToMatchDegree で必要な構造体を定義する。
// NearTriangulation と outer_face の情報を管理する構造体
struct NearTriangulationWithRingInfo {
    NearTriangulation near_triangulation;
    // num_new_neighbors[i] := outer_face の i 番目の頂点が outer_face にのばす辺の数
    vector<int> num_new_neighbors;
    // before[i] := outer_face のサイクル上で i 番目の頂点の "前の" 頂点 (最初は i - 1 に設定されているが、操作をするにつれて変化する)
    vector<int> before;
    // after[i] := outer_face のサイクル上で i 番目の頂点の "後の" 頂点 (最初は i + 1 に設定されているが、操作をするにつれて変化する)
    vector<int> after;
    // is_outer_face[i] := outer_face の i 番目の頂点が 現時点で outer_face に接続しているかどうか (最初は true に設定されているが、操作をするにつれて変化する)
    vector<bool> is_outer_face;
    // tracking_number は頂点番号の変化を track するための変数
    vector<vector<int>> tracking_number;
    // 頂点の identify によって孤立点になり、消す予定の頂点番号
    vector<int> vertices_to_be_removed;
    // 現在の outer_face の長さ、具体的には is_outer_face[i] = true なる i の数
    int len;
    NearTriangulationWithRingInfo(const NearTriangulation &near_triangulation, 
        const vector<int> &num_new_neighbors, const vector<int> &before, const vector<int> &after, const vector<bool> &is_outer_face,
        const vector<vector<int>> &tracking_number, const vector<int> &vertices_to_be_removed, int len): 
            near_triangulation(near_triangulation), num_new_neighbors(num_new_neighbors), before(before), after(after), is_outer_face(is_outer_face),
            tracking_number(tracking_number), vertices_to_be_removed(vertices_to_be_removed), len(len) {}
};


// num_new_neighbors[i] < 0 であるような i について、outer_face の i 番目の頂点の近傍を適切に同一視する。
// 2 頂点以上の同一視が必要になったときに、reducible configurations を含んでいるかどうかを判定する。
bool handleNegativeNewNeighbor(int i, const vector<int> &outer_face, NearTriangulationWithRingInfo &info, const vector<Configuration> &confs) {
    assert(info.num_new_neighbors[i] < 0);
    NearTriangulation &current = info.near_triangulation;
    int v = outer_face[i];
    int num = -info.num_new_neighbors[i];
    vector<int> cyclic_vneighbors = current.getCyclicOrder(v, outer_face[info.before[i]]);
    assert(cyclic_vneighbors[0] == outer_face[info.before[i]]);
    assert(cyclic_vneighbors.back() == outer_face[info.after[i]]);
    int cyclic_vneighbors_size = (int)cyclic_vneighbors.size();

    for (int j = 0;j < num; j++) {
        int u0 = cyclic_vneighbors[j];
        int u1 = cyclic_vneighbors[j + cyclic_vneighbors_size - num];
        if (!Degree::mergable(current.degrees()[u0], current.degrees()[u1])) {
            // 同一視する頂点の次数があっていない。
            return false;
        }
    }

    if (num != 1) {
        // 将来的には辺の次数による高速化をしても良いかもしれない。
        if (containOneofConfs(current, confs)) return false;
        // 2 つ以上の頂点を同一視することはほとんどない。起きた場合はログにはいて、プログラムを終了させる。
        spdlog::info("current.debug: {}", current.debug());
        spdlog::info("outer_face: {}", outer_face);
        spdlog::info("num_new_neighbors: {}", info.num_new_neighbors);
        spdlog::info("v: {}", v);
        spdlog::info("cyclic_vneighbors: {}", cyclic_vneighbors);
        assert(false);
    }

    int w0 = outer_face[info.before[i]]; 
    int w1 = outer_face[info.after[i]]; 
    // 両方とも次数が定まっていないときは any_degree になる。
    Degree any_degree = Degree(MIN_DEGREE, MAX_DEGREE);
    Degree degw = Degree::intersection(current.degrees()[w0].value_or(any_degree), current.degrees()[w1].value_or(any_degree));

    bool res = current.identify(w1, w0); // ここで nearTriangulation が変わる。
    if (!res) return false;
    info.vertices_to_be_removed.push_back(w0);
    if (info.tracking_number[w0].size() > 0) {
        info.tracking_number[w1].insert(info.tracking_number[w1].end(), info.tracking_number[w0].begin(), info.tracking_number[w0].end());
        info.tracking_number[w0].clear();
    }
    // is_outer_face, before, after を変更する。
    info.is_outer_face[i] = false;
    info.is_outer_face[info.before[i]] = false;
    info.before[info.after[i]] = info.before[info.before[i]];
    info.after[info.before[info.before[i]]] = info.after[i];
    info.num_new_neighbors[info.after[i]] = degw.fixed() ? degw.lower() - current.numIncidentEdge(w1) : INF;
    info.len -= 2;
    return res;
}


// num_new_neighbors[i] = 0 であるような i について、outer_face の i 番目の頂点の近傍に辺を結ぶ。
bool handleZeroNewNeighbor(int i, const vector<int> &outer_face, NearTriangulationWithRingInfo &info) {
    assert(info.num_new_neighbors[i] == 0);
    if (info.len <= 3 && info.near_triangulation.vertexSize() > 3) { // 3 角形(以下)の中に 1 頂点以上含まれているため cyclical 3-cut がある。
        return false;
    }
    // 辺を足す。
    bool res = info.near_triangulation.addEdge(outer_face[info.before[i]], outer_face[info.after[i]]); // ここで neartriangulation が変わる。
    if (res) {
        info.before[info.after[i]] = info.before[i];
        info.after[info.before[i]] = info.after[i];
        info.num_new_neighbors[info.before[i]] -= 1;
        info.num_new_neighbors[info.after[i]] -= 1;
        info.is_outer_face[i] = false;
        info.len -= 1;
        assert(info.is_outer_face[info.before[i]] && info.is_outer_face[info.after[i]]);
    }
    return res;
}

bool handleNewNeighbor(int i, const vector<int> &outer_face, NearTriangulationWithRingInfo &info, const vector<Configuration> &confs) {
    assert(info.num_new_neighbors[i] <= 0);
    if (info.len <= 3 && info.near_triangulation.vertexSize() > 3) { // 3 角形(以下)の中に 1 頂点以上含まれているため cyclical 3-cut がある。
        return false;
    }
    if (info.num_new_neighbors[i] < 0) {
        bool valid = handleNegativeNewNeighbor(i, outer_face, info, confs);
        if (!valid) return false;
    } else {
        assert(info.num_new_neighbors[i] == 0);
        bool valid = handleZeroNewNeighbor(i, outer_face, info);
        if (!valid) return false;
    }
    return true;
}

bool handleNewNeighborRecursive(int i, const vector<int> &outer_face, NearTriangulationWithRingInfo &info, const vector<Configuration> &confs) {
    assert(info.is_outer_face[i] && info.num_new_neighbors[i] <= 0);
    bool valid = true;
    valid = valid && handleNewNeighbor(i, outer_face, info, confs);

    if (info.is_outer_face[info.before[i]] && info.num_new_neighbors[info.before[i]] <= 0) {
        valid = valid && handleNewNeighborRecursive(info.before[i], outer_face, info, confs);
    }
    if (info.is_outer_face[info.after[i]] && info.num_new_neighbors[info.after[i]] <= 0) {
        valid = valid && handleNewNeighborRecursive(info.after[i], outer_face, info, confs);
    }

    return valid;
}

bool handleNewNeighbors(int l, const vector<int> &outer_face, NearTriangulationWithRingInfo &info, const vector<Configuration> &confs) {
    bool valid = true;
    for (int i = 0;i < l; i++) {
        if (info.is_outer_face[i] && info.num_new_neighbors[i] <= 0) {
            valid = valid && handleNewNeighborRecursive(i, outer_face, info, confs);
            if (!valid) return false;
        }
    }
    return true;
}

// neartriangulation の outer_face に接続している頂点の 次数 と 接続している辺の数 に応じて、
// 適切に 頂点を同一視する、辺を付け足す などの操作を行う。そのような操作の結果 
// + valid である (e.g. 最小反例にない 2,3-vertex-cut などができたりしない) 
// ような neartriangulation を返す。
bool NearTriangulation::addEdgeToMatchDegree(const vector<Configuration> &confs, vector<int> &vids) {
    vector<int> outer_face = getOuterFace();
    vector<int> num_new_neighbors = calcNumNewNeighbors(outer_face);

    int l = (int)outer_face.size();
    vector<int> before(l, 0), after(l, 0);
    vector<bool> is_outer_face(l, true);
    for (int i = 0;i < l; i++) {
        before[i] = (i == 0 ? l - 1 : i - 1);
        after[i] = (i == l - 1 ? 0 : i + 1);
    }

    // vids の頂点を track する
    vector<vector<int>> tracking_number(vertex_size_);
    for (int i = 0;i < (int)vids.size(); i++) {
        tracking_number[vids[i]].push_back(i);
    }

    // 頂点の同一視をしたことで、孤立点になっている頂点の番号を入れておく。
    vector<int> vertices_to_be_removed;

    NearTriangulationWithRingInfo info = NearTriangulationWithRingInfo(*this, num_new_neighbors, before, after, is_outer_face, tracking_number, vertices_to_be_removed, l);
    bool valid = handleNewNeighbors(l, outer_face, info, confs);
    if (!valid) return false;

    // vids の処理
    vector<int> vid_changed;
    vector<vector<int>> vid_index;
    for (int v = 0;v < vertex_size_; v++) {
        if (info.tracking_number[v].size() > 0) {
            vid_changed.push_back(v);
            vid_index.push_back(info.tracking_number[v]);
        }
    }
    if (info.vertices_to_be_removed.size() > 0) {
        spdlog::debug("vertices_to_be_removed.size() : {}", info.vertices_to_be_removed.size());
        vid_changed = info.near_triangulation.removeVertices(info.vertices_to_be_removed, vid_changed);
    }
    for (int i = 0;i < (int)vid_index.size(); i++) {
        for (int j: vid_index[i]) {
            vids[j] = vid_changed[i];
        }
    }
    for (int v: vids) {
        assert(v != -1);
    }
    if (info.near_triangulation.exist3cut()) return false;
    *this = info.near_triangulation;
    
    return true;
}

// NearTriangulation の outer_face を拡大する。
// 元からある頂点の番号は変わらない。
bool NearTriangulation::extendNeighbor(void) {
    vector<int> outer_face = this->getOuterFace();
    vector<int> num_new_neighbors = this->calcNumNewNeighbors(outer_face);
    for (auto num : num_new_neighbors) {
        if (num <= 0) {
            spdlog::info("debug: \n{}", debug());
            spdlog::info("outer_face: {}", fmt::join(outer_face, ", "));
            spdlog::info("num_new_neighbors: {}", fmt::join(num_new_neighbors, ", "));
            assert(false);
            // extendNeighbor 関数を呼ぶ前には addEdgeToMatchDegree を読んで num > 0 になっていなければならない。
        }
    }
    
    // outer_face の頂点に新しく接続される辺であって、その端点を共有しているものをグループ化する。
    // num_new_neighbors[i] = 1 のときは、
    // left_group[i] = right_group[i] で 1 本の辺が所属しているグループを表し、
    // num_new_neighbors[i] > 1 のときは、
    // left_group[i], right_group[i] で左、右に出ている辺のグループを表す。
    int l = (int)outer_face.size();
    int n_group = 0;
    vector<int> left_group(l, -1), right_group(l, -1);
    {
        for (int i = 0;i < l; i++) {
            if (num_new_neighbors[i] == 1) {
                left_group[i] = n_group;
                right_group[i] = n_group;
            } else {
                left_group[i] = n_group++;
                right_group[i] = n_group;
            }
        }
        assert(n_group > 0);
        int i = l - 1;
        for (;num_new_neighbors[i] == 1;i--) {
            right_group[i] = left_group[0];
            left_group[i] = left_group[0];
        }
        right_group[i] = left_group[0];
    }
    // num_new_neighbors[i] = num_new_neighbors[(i+1)] = INF となっているときは、
    // i, i+1 番目の頂点に隣接する group の neighbor は作る必要がない。
    vector<bool> valid_group(n_group, true);
    {
        for (int i = 0;i < l; i++) {
            int i_next = (i == l - 1 ? 0 : i + 1);
            if (num_new_neighbors[i] == INF && num_new_neighbors[i_next] == INF) {
                assert(right_group[i] == left_group[i_next]);
                valid_group[right_group[i]] = false;
            }
        }
    }

    // それぞれの group に所属する辺の端点となる neighbor を付け足す。
    vector<int> group_neighbors(n_group, -1);
    for (int i = 0;i < n_group; i++) {
        if (valid_group[i]) group_neighbors[i] = newVertex();
    }

    bool valid = true;
    for (int i = 0;i < l; i++) {
        if (num_new_neighbors[i] == 1) {
            assert(left_group[i] == right_group[i]);
            valid = valid && addEdge(outer_face[i], group_neighbors[left_group[i]]);
        } else {
            if (left_group[i] == right_group[i]) {
                return false; // 意図しない 2-cut がある。
            }
            int first = group_neighbors[left_group[i]];
            int last = group_neighbors[right_group[i]];
            if (first != -1) valid = valid && addEdge(outer_face[i], first);
            if (last != -1) valid = valid && addEdge(outer_face[i], last);
            if (num_new_neighbors[i] == INF) {
                // 9+ みたいな次数のとき
                continue;
            }
            for (int c = 0;c < num_new_neighbors[i] - 2; c++) {
                int v = newVertex();
                valid = valid && addEdge(v, first);
                valid = valid && addEdge(v, outer_face[i]);
                first = v;
            }
            valid = valid && addEdge(first, last);
        }
    }

    return valid;
}

// 辺をその端点の次数で分類する。
// min_degree <= deg0, deg1 < max_degree であるような deg0, deg1 について、
// edgeids[deg0 - min_degree][deg1 - min_degree] := {辺 e = uv であって u の次数が deg0, v の次数が deg1 であるようなもの}
// を計算する。
// + 次数が決まっていない頂点 (degree が std::nullopt)
// + fixed ではない頂点
// のときは実装していない。
// そのため、Configuration が含まれる可能性がある (Posssible) の判定のときはうまくいかない。
// Configuration が含まれているかどうかの判定のときはこの実装で良い。
vector<vector<vector<int>>> NearTriangulation::edgeIdsByclassifyingDegree(int min_degree, int max_degree) const {
    vector<vector<vector<int>>> edgeids(max_degree - min_degree, vector<vector<int>>(max_degree - min_degree));
    for (size_t ei = 0;ei < edges_.size();ei++) {
        auto [u, v] = edges_[ei];
        if (!degrees_[u].has_value() || !degrees_[u].value().fixed()) {
            continue;
        }
        if (!degrees_[v].has_value() || !degrees_[v].value().fixed()) {
            continue;
        }
        int degu = degrees_[u].value().lower();
        int degv = degrees_[v].value().lower();
        assert(min_degree <= degu && degu < max_degree && min_degree <= degv && degv < max_degree);
        edgeids[degu - min_degree][degv - min_degree].push_back(ei);
    }
    return edgeids;
}

string NearTriangulation::debug(void) const {
    string buf = "";
    vector<set<int>> VtoV(vertex_size_);
    for (const auto &e : edges_) {
        VtoV[e.first].insert(e.second);
        VtoV[e.second].insert(e.first);
    }
    for (int v = 0;v < vertex_size_; v++) {
        buf += fmt::format("{} {} {}\n", v, degrees_[v] ? degrees_[v].value().toString() : "?", fmt::join(VtoV[v], ", "));
    }
    return buf;
}

