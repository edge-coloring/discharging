#include <fstream>
#include <filesystem>
#include <spdlog/spdlog.h>
#include "configuration.hpp"

namespace fs = std::filesystem;

bool confHasCutVertex(int vertex_size, int ring_size, const vector<set<int>> &VtoV) {
    bool has_cutvertex = false;
    int ord = 0;
    vector<int> num(vertex_size, -1);
    vector<int> low(vertex_size, -1);
    auto dfs = [&](auto &&dfs, int v, int par) -> void {
        num[v] = ord++;
        low[v] = num[v];
        int n_child = 0;
        for (int u : VtoV[v]) {
            if (u == par) continue;
            if (u < ring_size) continue;
            if (num[u] != -1) {
                low[v] = std::min(low[v], num[u]);
                continue;
            }
            n_child++;
            dfs(dfs, u, v);
            low[v] = std::min(low[v], low[u]);
            if (par != -1 && num[v] <= low[u]) has_cutvertex = true;
        }
        if (par == -1 && n_child >= 2) has_cutvertex = true;
        return;
    };
    dfs(dfs, ring_size, -1);
    spdlog::trace("num : {}", fmt::join(num, ", "));
    spdlog::trace("low :{}", fmt::join(low, ", "));
    return has_cutvertex;
}


Configuration::Configuration(int ring_size, bool ring_remains, const string &filename, const NearTriangulation &conf) :
    conf_(conf), ring_size_(ring_size), ring_remains_(ring_remains), filename_(filename) {
    // 次数が大きくなるほど出現頻度が下がると仮定して、
    // 高速化のために、端点が辞書順最大になるような辺を inside_edge_id_ として選ぶ。
    int offset = ring_remains ? ring_size_ : 0;
    int degu_max = 0, degv_max = 0;
    inside_edge_id_ = -1;
    for (size_t ei = 0;ei < conf_.edges().size(); ei++) {
        auto [u, v] = conf_.edges()[ei];
        if (u < offset || v < offset) continue;
        int degu = conf_.degrees()[u].value().lower();
        int degv = conf_.degrees()[v].value().lower();
        if (degu > degu_max || (degu == degu_max && degv >= degv_max)) {
            degu_max = degu;
            degv_max = degv;
            inside_edge_id_ = (int)ei;
        }
    }
    assert(inside_edge_id_ != -1);
};

// ring_remains = True のときは必ず ring を残すようにする。
Configuration Configuration::readConfFile(const string &filename, bool ring_remains) {
    std::ifstream ifs(filename);
    if (!ifs) {
        spdlog::critical("Failed to open {} ", filename);
        throw std::runtime_error("Failed to open" + filename);
    }
    string dummy;
    std::getline(ifs, dummy);
    int vertex_size, ring_size;
    ifs >> vertex_size >> ring_size;

    vector<set<int>> VtoV(vertex_size);
    vector<optional<Degree>> degrees(vertex_size, std::nullopt);

    for (int vi = 0;vi < ring_size; vi++) {
        int vip = (vi + 1) % ring_size;
        VtoV[vi].insert(vip);
        VtoV[vip].insert(vi);
    }
    for (int vi = ring_size;vi < vertex_size; vi++) {
        int v, degv;
        ifs >> v >> degv;
        --v;
        assert(v == vi);
        degrees[v] = Degree(degv);
        for (int i = 0;i < degv; i++) {
            int nv;
            ifs >> nv;
            --nv;
            assert(0 <= nv && nv < vertex_size);
            VtoV[v].insert(nv);
            VtoV[nv].insert(v);
        }
    }

    if (ring_remains) {
        return Configuration(ring_size, true, filename, NearTriangulation(vertex_size, VtoV, degrees));
    }

    if (confHasCutVertex(vertex_size, ring_size, VtoV)) {
        spdlog::trace("has cut vertex");
        return Configuration(ring_size, true, filename, NearTriangulation(vertex_size, VtoV, degrees));
    }
    spdlog::trace("has no cut vertex");

    // delete ring 
    vector<set<int>> VtoV2(vertex_size - ring_size);
    for (int v = ring_size;v < vertex_size; v++) {
        for (auto u : VtoV[v]) {
            if (u < ring_size) {
                continue;
            }
            VtoV2[v - ring_size].insert(u - ring_size);
        }
    }
    degrees.erase(degrees.begin(), degrees.begin() + ring_size);
    vertex_size -= ring_size;
    
    return Configuration(ring_size, false, filename, NearTriangulation(vertex_size, VtoV2, degrees));
}

void Configuration::writeConfFile(const string &filename) const {
    assert(ring_remains_);
    vector<vector<int>> VtoV(conf_.vertexSize());
    for (const auto &e : conf_.edges()) {
        if (e.first < ring_size_) continue;
        VtoV[e.first].push_back(e.second + 1);
    }
    string txt = "\n";
    txt += fmt::format("{} {}\n", conf_.vertexSize(), ring_size_);
    for (int v = ring_size_;v < conf_.vertexSize(); v++) {
        assert(conf_.degrees()[v].has_value());
        txt += fmt::format("{} {} {}\n", v+1, conf_.degrees()[v].value().lower(), fmt::join(VtoV[v], " "));
    }
    std::ofstream ofs(filename);
    ofs << txt;
    return;
}

const NearTriangulation &Configuration::nearTriangulation(void) const {
    return conf_;
}

int Configuration::ringSize(void) const {
    return ring_size_;
}

// 両端点がどちらも ring にない辺を一つ返す。
int Configuration::getInsideEdgeId(void) const {
    return inside_edge_id_;
}

bool Configuration::ringRemains(void) const {
    return ring_remains_;
}

const string &Configuration::fileName(void) const {
    return filename_;
}


// configuration の直径を計算する。
// 注意: ring の頂点を通るパスは考えない。
int Configuration::diameter(void) const {
    int vertex_size = conf_.vertexSize();
    int offset = ring_remains_ ? ring_size_ : 0;
    vector<vector<int>> dist(vertex_size - offset, vector<int>(vertex_size - offset, 10000));
    for (int v = offset;v < vertex_size; v++) {
        dist[v - offset][v - offset] = 0;
    }
    for (const auto &e: conf_.edges()) {
        if (e.first - offset < 0 || e.second - offset < 0) continue;
        dist[e.first - offset][e.second - offset] = 1;
    }
    for (int k = 0;k < vertex_size - offset; k++) {
        for (int i = 0;i < vertex_size - offset; i++) {
            for (int j = 0;j < vertex_size - offset; j++) {
                dist[i][j] = std::min(dist[i][j], dist[i][k] + dist[k][j]);
            }
        }
    }
    int diam = 0;
    for (int i = 0;i < vertex_size - offset; i++) {
        for (int j = 0;j < vertex_size - offset; j++) {
            diam = std::max(diam, dist[i][j]);
        }
    }
    return diam;
}

// near_triangulation に ring をつけて configuration にする。
// configuration にするにあたって、次数が std::nullopt の頂点や、次数が fixed ではない頂点は生成する configuration には含めない。
// 注意: 引数に与えた near_triangulation は変化するかもしれない。
// 将来的には、頂点番号の変化をトラックするようにしても良いかもしれない。
optional<Configuration> convertNearTriangulation2Conf(NearTriangulation &near_triangulation) {
    // ring の分の頂点を生成するために extendNeighbor 関数を呼ぶ。
    vector<int> vids;
    if (!near_triangulation.addEdgeToMatchDegree({}, vids)) return std::nullopt;
    if (!near_triangulation.extendNeighbor()) return std::nullopt;

    // outer_face の頂点は nullopt, fixed でないのどちらかであることを確認する。
    vector<int> outer_face = near_triangulation.getOuterFace();
    for (int v : outer_face) {
        assert(!near_triangulation.degrees()[v].has_value() || !near_triangulation.degrees()[v].value().fixed());
    }

    // outer_face の頂点が ring になるように頂点番号を振り直す。
    int ringsize = (int)outer_face.size();
    vector<int> new_index(near_triangulation.vertexSize(), -1);
    {
        for (int i = 0;i < ringsize; i++) {
            new_index[outer_face[i]] = i;
        }
        int c = ringsize;
        for (int v = 0;v < near_triangulation.vertexSize(); v++) {
            if (new_index[v] == -1) new_index[v] = c++;
        }
        assert(c == near_triangulation.vertexSize());
    }

    vector<set<int>> new_VtoV(near_triangulation.vertexSize());
    {
        for (int v = 0;v < near_triangulation.vertexSize(); v++) {
            for (int u : near_triangulation.VtoV()[v]) {
                new_VtoV[new_index[v]].insert(new_index[u]);
            }
        }
    }

    vector<optional<Degree>> new_degrees(near_triangulation.vertexSize(), std::nullopt);
    {
        for (int v = 0;v < near_triangulation.vertexSize(); v++) {
            new_degrees[new_index[v]] = near_triangulation.degrees()[v];
        }
    }

    return Configuration(ringsize, true, "", NearTriangulation(near_triangulation.vertexSize(), new_VtoV, new_degrees));
}

// ディレクトリに含まれる　conf ファイルの configuration を返す。
vector<Configuration> getConfs(const std::string &dirname, bool ring_remains) {
    vector<Configuration> confs;
    spdlog::info("reading confs from {} ...", dirname);
    for (const fs::directory_entry &file : fs::directory_iterator(dirname)) {
        if (file.is_regular_file() && file.path().extension().string<char>() == ".conf") {
            spdlog::trace("reading {}", file.path().string<char>());
            confs.push_back(Configuration::readConfFile(file.path().string<char>(), ring_remains));
        } 
    }
    spdlog::info("finish reading confs from {} ...", dirname);
    return confs;
}