#pragma once
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <fstream>
#include <set>
#include <optional>

using std::vector;
using std::pair;
using std::map;
using std::tuple;
using std::string;
using std::set;
using std::optional;

class Configuration;

// 次数は高々12くらいまでしか考えないので inf として1000を設定
const int MAX_DEGREE = 1000;
const int MIN_DEGREE = 5;
// 無限大
const int INF = 1000000;

class Degree {
private:
    int lower_deg_, upper_deg_;

public:
    Degree(int lower_deg, int upper_deg);
    Degree(int deg);
    static Degree fromString(const string &str);

    int lower(void) const;
    int upper(void) const;

    string toString(void) const;
    bool include(const Degree &degree) const;
    static bool disjoint(const Degree &degree0, const Degree &degree1);
    bool fixed(void) const;
    static Degree intersection(const Degree &degree0, const Degree &degree1);
    static bool mergable(const optional<Degree> &deg0, const optional<Degree> &deg1);
};

bool operator==(const Degree& lhs, const Degree& rhs);

vector<Degree> divideDegree(const Degree &degree, int max_degree);

class NearTriangulation {
private:
    int vertex_size_;
    // 頂点の次数
    // std::nullopt はまだ次数が定まっていない状態を表す。
    vector<optional<Degree>> degrees_;
    // 辺集合
    vector<pair<int, int>> edges_;
    // diagonal_vertices_[e] := 辺 e を含む三角形の頂点であって、 e の端点でない頂点の集合
    map<pair<int, int>, vector<int>> diagonal_vertices_;
    // VtoV_[v] := v に隣接する頂点の集合
    vector<set<int>> VtoV_;

public:
    NearTriangulation(int vertex_size, const vector<set<int>> &VtoV, const vector<optional<Degree>> &degrees);

    int vertexSize(void) const;
    const vector<optional<Degree>> &degrees(void) const;
    const vector<pair<int, int>> &edges(void) const;
    const map<pair<int, int>, vector<int>> &diagonalVertices(void) const;
    const vector<set<int>> &VtoV(void) const;

    // v に接続している辺の数
    int numIncidentEdge(int v) const;
    // v の次数を設定する。
    void setDegree(int v, const optional<Degree> &degree);
    // 新しい頂点を付け足す
    int newVertex(void);
    // 辺 uv を足す。
    bool addEdge(int u, int v);
    // 頂点 u, v を同一視する。
    bool identify(int u, int v);
    // 孤立点である頂点 isolated_vertices を削除する。
    vector<int> removeVertices(const vector<int> &isolated_vertices, const vector<int> &vids);
    // near-triangulation の outer face (3 角形ではない面) の頂点を順に並べる。
    vector<int> getOuterFace(void) const;
    // 3-cut が存在するかどうか
    bool exist3cut(void) const;
    // v とその近傍 neighbor_v から v の近傍の平面埋め込みの時計回り順を得る。
    vector<int> getCyclicOrder(int v, int v_neighbor) const;
    // outer_face の頂点が追加で後何本の辺に接続するかを計算する。
    vector<int> calcNumNewNeighbors(const vector<int> &outer_face) const;
    // outer_face の頂点から辺が少なくとも 1 本は出るように, 辺を足したり頂点を同一視したりしてグラフを修正する。
    bool addEdgeToMatchDegree(const vector<Configuration> &confs, vector<int> &vids);
    // 近傍を拡大する。
    bool extendNeighbor(void);

    // 辺を端点の次数で分類する。
    vector<vector<vector<int>>> edgeIdsByclassifyingDegree(int min_degree, int max_degree) const;
    string debug(void) const;
};

bool containOneofConfs(const NearTriangulation &wheelgraph, const vector<Configuration> &confs);

