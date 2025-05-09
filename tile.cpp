
#include "basewheel.hpp"
#include "configuration.hpp"
#include "thread_pool.hpp"
#include <spdlog/spdlog.h>
#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <filesystem>
namespace fs = std::filesystem;
const int DEGREE_UPPER_BOUND = 13;

// near_triangulation に conf が含まれているかどうかをチェックする。
// near_triangulation の辺の端点の次数の情報 edgeids を使って高速化している。
bool containConfs(const NearTriangulation &near_triangulation, const vector<vector<vector<int>>> &edgeids, const vector<Configuration> &confs) {
    for (const auto &conf : confs) {
        int edgeid_conf = conf.getInsideEdgeId();
        auto [u, v] = conf.nearTriangulation().edges()[edgeid_conf];
        Degree degu = conf.nearTriangulation().degrees()[u].value();
        Degree degv = conf.nearTriangulation().degrees()[v].value();
        set<int> ring_vertices;
        if (conf.ringRemains()) {
            for (int r = 0;r < conf.ringSize(); r++) ring_vertices.insert(r);
        }
        for (int du = degu.lower();du <= degu.upper(); du++) {
            for (int dv = degv.lower();dv <= degv.upper(); dv++) {
                assert(0 <= du - MIN_DEGREE && du - MIN_DEGREE < (int)edgeids.size());
                assert(0 <= dv - MIN_DEGREE && dv - MIN_DEGREE < (int)edgeids[du - MIN_DEGREE].size());
                for (auto ei : edgeids[du - MIN_DEGREE][dv - MIN_DEGREE]) {
                    if (BaseWheel::numOfSubgraphWithCorrespondingEdge(near_triangulation, conf.nearTriangulation(), ei, edgeid_conf, ring_vertices) > 0) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

// near_triangulation1 と near_triangulation2 から、occupied1 のような頂点の対応関係を考えて、
// 組み合わせることで新しい near_triangulation を生成する。
// occupied1[v1] := (near_triangulation1 の頂点 v1 に対応する near_triangulation2 の頂点番号)
// vids1 に含まれている頂点番号の near_triangulation1 の頂点の番号の変化をトラックする。 
// vids2 に含まれている頂点番号の near_triangulation2 の頂点の番号の変化をトラックする。 
NearTriangulation genNearTriangulationByTiling(const NearTriangulation &near_triangulation1, const NearTriangulation &near_triangulation2,
    const vector<int> &occupied1, vector<int> &vids1, vector<int> &vids2) {
    
    vector<int> new_index1(near_triangulation1.vertexSize(), -1);
    vector<int> new_index2(near_triangulation2.vertexSize(), -1);
    int new_vertex_size = -1;
    {
        int c = 0;
        for (int v1 = 0;v1 < near_triangulation1.vertexSize(); v1++) {
            new_index1[v1] = c;
            if (occupied1[v1] != -1) new_index2[occupied1[v1]] = c;
            c++;
        }
        for (int v2 = 0;v2 < near_triangulation2.vertexSize(); v2++) {
            if (new_index2[v2] == -1) {
                new_index2[v2] = c++;
            }
        }
        new_vertex_size = c;
    }

    vector<optional<Degree>> degrees(new_vertex_size, std::nullopt);
    {
        for (int v1 = 0;v1 < near_triangulation1.vertexSize(); v1++) {
            degrees[new_index1[v1]] = near_triangulation1.degrees()[v1];
        }
        for (int v2 = 0;v2 < near_triangulation2.vertexSize(); v2++) {
            if (!degrees[new_index2[v2]].has_value()) {
                degrees[new_index2[v2]] = near_triangulation2.degrees()[v2];
            } else {
                assert(!near_triangulation2.degrees()[v2].has_value()
                || degrees[new_index2[v2]] == near_triangulation2.degrees()[v2]);
            }
        }
    }

    vector<set<int>> VtoV(new_vertex_size);
    {
        for (const auto &e : near_triangulation1.edges()) {
            int v = new_index1[e.first];
            int u = new_index1[e.second];
            VtoV[v].insert(u);
        }
        for (const auto &e : near_triangulation2.edges()) {
            int v = new_index2[e.first];
            int u = new_index2[e.second];
            VtoV[v].insert(u);
        }
        for (int v = 0;v < new_vertex_size; v++) {
            if (degrees[v].has_value() && degrees[v].value().fixed()) {
                if ((int)VtoV[v].size() != degrees[v].value().lower()) {
                    throw std::runtime_error(fmt::format("degrees[{}] = {}, VtoV[{}] = {}, which is different", v, degrees[v].value().lower(), v, fmt::join(VtoV[v], " ")));
                }
            }
        }
    }

    for (size_t i = 0;i < vids1.size(); i++) {
        vids1[i] = new_index1[vids1[i]];
    }
    for (size_t i = 0;i < vids2.size(); i++) {
        vids2[i] = new_index2[vids2[i]];
    }

    return NearTriangulation(new_vertex_size, VtoV, degrees);
}

// near_triangulation に引数 degree で与えられる次数を設定してできる新しい near_triangulation を返す。
// + near_triangulation と vids は参照で渡して、変更される。
// + 返り値は near_triangulation が次数を設定して妥当なもの (3-cutが存在しない,など) のままだったかどうかを表すフラグ。
bool updateNearTriangulationBySettingDegree(NearTriangulation &near_triangulation, const vector<optional<Degree>> &degree, vector<int> &vids, const vector<Configuration> &reducibles) {
    assert(near_triangulation.vertexSize() == (int)degree.size());
    for (int v = 0;v < near_triangulation.vertexSize(); v++) {
        if (!degree[v].has_value()) continue;
        if (near_triangulation.degrees()[v].has_value() && near_triangulation.degrees()[v].value() == degree[v].value()) continue;
        assert(!near_triangulation.degrees()[v].has_value());
        near_triangulation.setDegree(v, degree[v]);
    }
    bool res = near_triangulation.addEdgeToMatchDegree(reducibles, vids);
    res = res && near_triangulation.extendNeighbor();
    return res;
}

// occupied1[v1] := (near_triangulation1 の頂点 v1 に対応する near_triangulation2 の頂点)
// が与えられるので、
// occupied2[v2] := (near_triangulation2 の頂点 v2 に対応する near_triangulation1 の頂点) 
// を計算する。ただし、対応するものがない場合は -1 とする。
vector<int> calcOccupied2from1(int vertex_size1, int vertex_size2, const vector<int> &occupied1) {
    assert(vertex_size1 == (int)occupied1.size());
    vector<int> occupied2(vertex_size2, -1);
    for (int v1 = 0;v1 < vertex_size1; v1++) {
        if (occupied1[v1] == -1) continue;
        occupied2[occupied1[v1]] = v1;
    }
    return occupied2;
}

// occupied1[v1] := (near_triangulation1 の頂点 v1 に対応する near_triangulation2 の頂点) と
// near_triangulation2 の次数 degree2 が与えられるから、
// near_triangulation1 の頂点に新しく設定する次数 degree1 を計算する。
vector<optional<Degree>> calcDegreeByOccupied(int vertex_size1, const vector<int> &occupied1, const vector<optional<Degree>> &original_degrees1, const vector<optional<Degree>> &degrees2) {
    vector<optional<Degree>> new_degree1(vertex_size1, std::nullopt);
    for (int v1 = 0;v1 < vertex_size1; v1++) {
        if (occupied1[v1] != -1 && degrees2[occupied1[v1]].has_value()) {
            if (original_degrees1[v1].has_value() && original_degrees1[v1].value() == degrees2[occupied1[v1]].value()) continue;
            assert(!original_degrees1[v1].has_value()); 
            new_degree1[v1] = degrees2[occupied1[v1]].value();
        }
    }
    return new_degree1;
}

// near_triangulation1 の頂点 v1, u1
// near_triangulation2 の頂点 v2, u2
// について、 v1 を v2 に、 u1 を u2 にそれぞれ対応させてできる near_triangulation を生成する。 
// + NearTriangulation, (v1,v2 に対応する頂点), (u1,u2 に対応する頂点) を 3 つ組にして返す。
vector<tuple<NearTriangulation, int, int>> tileByCorrespondingEdge(const NearTriangulation &near_triangulation1, const NearTriangulation &near_triangulation2,
    int v1, int u1, int v2, int u2, const vector<Configuration> &reducibles) {
    vector<tuple<NearTriangulation, int, int>> combined_near_triangulations;

    pair<int, int> e1 = std::make_pair(v1, u1);
    pair<int, int> e2 = std::make_pair(v2, u2);
    assert(near_triangulation1.diagonalVertices().at(e1).size() >= 1);
    assert(near_triangulation2.diagonalVertices().at(e2).size() >= 1);
    int w1 = near_triangulation1.diagonalVertices().at(e1)[0];
    int w2 = near_triangulation2.diagonalVertices().at(e2)[0];
    for (int t = 0;t < 2; t++) {
        // t == 0 のときは w1 と w2 が対応しない
        // t == 1 のときは w1 と w2 が対応する。
        // ときをそれぞれ表している。
        NearTriangulation new_near_triangulation1 = near_triangulation1;
        NearTriangulation new_near_triangulation2 = near_triangulation2;
        vector<int> vids1 = {v1, u1, w1};
        vector<int> vids2 = {v2, u2, w2};
        bool to_check = true;
        vector<int> occupied1;
        while (1) {
            auto contains1_2 = BaseWheel::containSubgraphWithCorrespondingEdge(new_near_triangulation1, new_near_triangulation2, 
                std::make_pair(vids1[0], vids1[1]), std::make_pair(vids2[0], vids2[1]), {}, true);
            assert(contains1_2.size() <= 2);
            // 1 つの辺を一致させるようにするとき、2 通りありうるが、
            // t = 0, 1 の条件が満たされるのはそのうち 1 通り(以下)しかない。
            bool seen = false;
            for (const auto &contain1_2 : contains1_2) {
                if (contain1_2.contain == Contain::No) continue;
                if ((t == 0) == (contain1_2.occupied[vids1[2]] == vids2[2])) continue;
                assert(!seen);
                seen = true;
            }
            if (!seen) {
                to_check = false;
                break;
            }
            bool nochange = false;
            for (const auto &contain1_2 : contains1_2) {
                if (contain1_2.contain == Contain::No) continue;
                if ((t == 0) == (contain1_2.occupied[vids1[2]] == vids2[2])) continue;
                // contain1_2.occupied[v1] := (new_near_triangulation1 の頂点 v1 に対応する new_near_triangulation2 の頂点番号)
                // occupied2[v2] := (new_near_triangulation2 の頂点 v2 に対応する new_near_triangulation1 の頂点番号)
                vector<int> occupied2 = calcOccupied2from1(new_near_triangulation1.vertexSize(), new_near_triangulation2.vertexSize(), contain1_2.occupied);
                // degree1[v1] := (new_near_triangulation1 の頂点 v1 に新しく設定する次数)
                vector<optional<Degree>> degree1 = calcDegreeByOccupied(new_near_triangulation1.vertexSize(), contain1_2.occupied, new_near_triangulation1.degrees(), new_near_triangulation2.degrees());
                // degree2[v2] := (new_near_triangulation2 の頂点 v2 に新しく設定する次数)
                vector<optional<Degree>> degree2 = calcDegreeByOccupied(new_near_triangulation2.vertexSize(), occupied2, new_near_triangulation2.degrees(), new_near_triangulation1.degrees());
                // もし次数の変化がなければ break
                auto no_change = [](const optional<Degree> &d) {return !d.has_value();};
                if (std::all_of(degree1.begin(), degree1.end(), no_change) && std::all_of(degree2.begin(), degree2.end(), no_change)) {
                    occupied1 = contain1_2.occupied;
                    nochange = true;
                    break;
                }
                // ring の次数を決定して、それによるグラフの変化を反映した near_triangualtion を作る。
                // new_near_triangulation1
                if (!updateNearTriangulationBySettingDegree(new_near_triangulation1, degree1, vids1, reducibles)) {
                    to_check = false;
                    break;
                }
                // new_near_triangulation2
                if (!updateNearTriangulationBySettingDegree(new_near_triangulation2, degree2, vids2, reducibles)) {
                    to_check = false;
                }
                break;
            }
            if (!to_check) {
                break;
            }
            if (nochange) {
                break;
            }
        }
        if (to_check) {
            // v1, v2 の頂点変化を追う
            vector<int> new_vids1 = {vids1[0], vids1[1]}, new_vids2 = {};
            // 対応関係から新しい configuration を作る。
            NearTriangulation combined_near_triangulation = genNearTriangulationByTiling(new_near_triangulation1, new_near_triangulation2, occupied1, new_vids1, new_vids2);
            // edgeids[i][j] := 頂点の次数が (i+5), (j+5) の辺の集合 を計算する。
            auto edgeids = combined_near_triangulation.edgeIdsByclassifyingDegree(MIN_DEGREE, DEGREE_UPPER_BOUND);
            // reducible configuration を含むかどうか判定する
            bool contain = containConfs(combined_near_triangulation, edgeids, reducibles);
            if (!contain) {
                combined_near_triangulations.emplace_back(combined_near_triangulation, new_vids1[0], new_vids1[1]);
            }
        }
    }
    
    return combined_near_triangulations;
}

// near_triangulation1 の頂点番号 v1 の頂点と near_triangulation2 の頂点番号 v2 の 
// 2 つの頂点が "隣接" するような状況を列挙して、
// + reducibles に含まれる configuration を部分グラフとして含まない。
// ような near_triangulation とその near_triangulation における v1 に対応する頂点番号と v2 に対応する頂点番号からなる 3 つ組を返す。
vector<tuple<NearTriangulation, int, int>> adjacentCenterInternal(const NearTriangulation &near_triangulation1, const NearTriangulation &near_triangulation2, int v1, int v2, const vector<Configuration> &reducibles) {
    // near_triangulation1 の頂点番号 v1 の頂点と near_triangulation2 の頂点番号 v2 の頂点が隣り合っているときの他の頂点の対応関係を調べる。
    vector<tuple<NearTriangulation, int, int>> combined_near_triangulations;
    assert(v1 < (int)near_triangulation1.VtoV().size());
    assert(v2 < (int)near_triangulation2.VtoV().size());
    for (int v_adjacent_v1 : near_triangulation1.VtoV()[v1]) {
        for (int v_adjacent_v2 : near_triangulation2.VtoV()[v2]) {
            auto near_triangulations = tileByCorrespondingEdge(near_triangulation1, near_triangulation2, v1, v_adjacent_v1, v_adjacent_v2, v2, reducibles);
            combined_near_triangulations.insert(combined_near_triangulations.end(), near_triangulations.begin(), near_triangulations.end());
        }
    }

    return combined_near_triangulations;
}

// near_triangulation1 の頂点番号 v1 の頂点と near_triangulation2 の頂点番号 v2 の 
// 2 つの頂点が "一致" するような状況を列挙して、
// + reducibles に含まれる configuration を部分グラフとして含まない。
// ような near_triangulation とその near_triangulation における v1, v2 に対応する頂点番号からなる 2 つ組を返す。
vector<tuple<NearTriangulation, int>> correspondCenterInternal(const NearTriangulation &near_triangulation1, const NearTriangulation &near_triangulation2, int v1, int v2, const vector<Configuration> &reducibles) {
    vector<tuple<NearTriangulation, int>> combined_near_triangulations;
    int v_adjacent_v1 = *near_triangulation1.VtoV()[v1].begin();
    for (int v_adjacent_v2 : near_triangulation2.VtoV()[v2]) {
        auto near_triangulations = tileByCorrespondingEdge(near_triangulation1, near_triangulation2, v1, v_adjacent_v1, v2, v_adjacent_v2, reducibles);
        for (const auto &[near_triangulation, v, v_adjacent] : near_triangulations) {
            combined_near_triangulations.emplace_back(near_triangulation, v);
        }
    }
    return combined_near_triangulations;
}


void adjacentCenter(const NearTriangulation &near_triangulation1, const NearTriangulation &near_triangulation2, 
    const pair<int, int> &vs, const pair<string, string> &basenames, const vector<Configuration> &reducibles,
    const string &output_dirname) {
    vector<Configuration> combined_confs;
    auto [v1, v2] = vs;
    const auto &[basename1, basename2] = basenames;
    // 各スレッドでエラーが出たときに他のスレッドを止めないようにするためにエラーをキャッチする。
    try {
        // near_triangulation の頂点番号 v1 の頂点と near_triangulation の頂点番号 v2 の 
        // 2 つの頂点が隣接するような状況を列挙する。
        vector<tuple<NearTriangulation, int, int>> combined_near_triangulations = 
            adjacentCenterInternal(near_triangulation1, near_triangulation2, v1, v2, reducibles);
        for (auto &[combined_near_triangulation, v1, v2] : combined_near_triangulations) {
            auto combined_conf = convertNearTriangulation2Conf(combined_near_triangulation);
            if (combined_conf.has_value()) {
                combined_confs.push_back(combined_conf.value());
            }
        }
    } catch (std::exception &e) {
        spdlog::critical("An error ({}) occured when combining {} and {}", e.what(), basename1, basename2);
        std::ofstream ofs(fmt::format("{}/{}_{}.err", output_dirname, basename1, basename2));
        ofs << e.what();
        return;
    }
    int count = 0;
    for (const auto &conf : combined_confs) {
        conf.writeConfFile(fmt::format("{}/{}_{}_{}.conf", output_dirname, count, basename1, basename2));
        count++;
    }
    return;
}

void correspondCenter(const NearTriangulation &near_triangulation1, const NearTriangulation &near_triangulation2, 
    const pair<int, int> &vs, const pair<string, string> &basenames, const vector<Configuration> &reducibles,
    const string &output_dirname) {
    vector<Configuration> combined_confs;
    auto [v1, v2] = vs;
    const auto &[basename1, basename2] = basenames;
    // 各スレッドでエラーが出たときに他のスレッドを止めないようにするためにエラーをキャッチする。
    try {
        // near_triangulation の頂点番号 v1 の頂点と near_triangulation の頂点番号 v2 の 
        // 2 つの頂点が一致するような状況を列挙する。
        vector<tuple<NearTriangulation, int>> combined_near_triangulations = 
            correspondCenterInternal(near_triangulation1, near_triangulation2, v1, v2, reducibles);
        for (auto &[combined_near_triangulation, _] : combined_near_triangulations) {
            auto combined_conf = convertNearTriangulation2Conf(combined_near_triangulation);
            if (combined_conf.has_value()) {
                combined_confs.push_back(combined_conf.value());
            }
        }
    } catch (std::exception &e) {
        spdlog::critical("An error ({}) occured when combining {} and {}", e.what(), basename1, basename2);
        std::ofstream ofs(fmt::format("{}/{}_{}.err", output_dirname, basename1, basename2));
        ofs << e.what();
        return;
    }
    int count = 0;
    for (const auto &conf : combined_confs) {
        conf.writeConfFile(fmt::format("{}/{}_{}_{}.conf", output_dirname, count, basename1, basename2));
        count++;
    }
    return;
}

string baseName(const string &path_str) {
    fs::path path = path_str;
    return path.stem();
}

// conf1_filename, conf2_filename の configuration について、
// R1 := (conf1 の ringsize), R2 := (conf2 の ringsize) としたとき、
// conf1 の頂点番号 R1 の頂点と、conf2 の頂点番号 R2 の頂点が隣接するようにタイリングした結果、
// reducible_dirname の configuration を含まないようなものを返す。
void checkAdjacentCenter(const string &conf1_filename, const string &conf2_filename, 
    const string &confs_dirname, const string &output_dirname) {
    Configuration conf1 = Configuration::readConfFile(conf1_filename, true);
    Configuration conf2 = Configuration::readConfFile(conf2_filename, true);
    vector<Configuration> reducibles = getConfs(confs_dirname);
    bool madedir = fs::create_directory(output_dirname);
    if (madedir) spdlog::info("made {} directory", output_dirname);
    string conf1_base = baseName(conf1_filename);
    string conf2_base = baseName(conf2_filename);
    adjacentCenter(conf1.nearTriangulation(), conf2.nearTriangulation(), std::make_pair(conf1.ringSize(), conf2.ringSize()), std::make_pair(conf1_base, conf2_base), reducibles, output_dirname);
    return;
}


// confs_dirname1 に含まれる configuration と
// confs_dirname2 に含まれる configuration を checkAdjacentCenter と同様にタイリングした結果、
// reducible_dirname の configuration を含まないようなものを返す。
// 並列化している。
void checkAdjacentCenterParallel1(const string &confs_dirname1, const string &confs_dirname2, const string &reducible_dirname, const string &output_dirname) {
    vector<Configuration> confs1 = getConfs(confs_dirname1, true);
    vector<Configuration> confs2 = getConfs(confs_dirname2, true);
    vector<Configuration> reducibles = getConfs(reducible_dirname);
    bool madedir = fs::create_directory(output_dirname);
    if (madedir) spdlog::info("made {} directory", output_dirname);
    size_t N1 = confs1.size(), N2 = confs2.size();
    vector<string> basenames1(N1), basenames2(N2);
    for (size_t i = 0;i < N1; i++) {
        basenames1[i] = baseName(confs1[i].fileName());
    }
    for (size_t i = 0;i < N2; i++) {
        basenames2[i] = baseName(confs2[i].fileName());
    }
    spdlog::info("start creating thread pool");
    boost::asio::io_service io_service;
    unsigned int num_concurrent = boost::thread::hardware_concurrency();
    ThreadPool pool(io_service, num_concurrent);
    spdlog::info("finish creating thread pool (the number of thread is {})", num_concurrent);
    spdlog::info("start queueing thread into pool");
    for (size_t i = 0;i < N1; i++) {
        for (size_t j = 0;j < N2; j++) {
            pair<int, int> vs = std::make_pair(confs1[i].ringSize(), confs2[j].ringSize());
            pair<string, string> basename = std::make_pair(basenames1[i], basenames2[j]);
            pool.post(boost::bind(adjacentCenter, boost::cref(confs1[i].nearTriangulation()), boost::cref(confs2[j].nearTriangulation()), vs, basename, boost::cref(reducibles), boost::cref(output_dirname)));
        }
    }
    spdlog::info("finish queueing thread into pool");
    return;
}

// confs_dirname の configuration たちに含まれる 2 つの configuration を checkAdjacentCenter と同様にタイリングした結果、
// reducible_dirname の configuration を含まないようなものを返す。
// 並列化している。
void checkAdjacentCenterParallel2(const string &confs_dirname, const string &reducible_dirname, const string &output_dirname) {
    vector<Configuration> confs = getConfs(confs_dirname, true);
    vector<Configuration> reducibles = getConfs(reducible_dirname);
    bool madedir = fs::create_directory(output_dirname);
    if (madedir) spdlog::info("made {} directory", output_dirname);
    size_t N = confs.size();
    vector<string> basenames(N);
    for (size_t i = 0;i < N; i++) {
        basenames[i] = baseName(confs[i].fileName());
    }
    spdlog::info("start creating thread pool");
    boost::asio::io_service io_service;
    unsigned int num_concurrent = boost::thread::hardware_concurrency();
    ThreadPool pool(io_service, num_concurrent);
    spdlog::info("finish creating thread pool (the number of thread is {})", num_concurrent);
    spdlog::info("start queueing thread into pool");
    for (size_t i = 0;i < N; i++) {
        for (size_t j = 0;j <= i; j++) {
            pair<int, int> vs = std::make_pair(confs[i].ringSize(), confs[j].ringSize());
            pair<string, string> basename = std::make_pair(basenames[i], basenames[j]);
            pool.post(boost::bind(adjacentCenter, boost::cref(confs[i].nearTriangulation()), boost::cref(confs[j].nearTriangulation()), vs, basename, boost::cref(reducibles), boost::cref(output_dirname)));
        }
    }
    spdlog::info("finish queueing thread into pool");
    return;
}

// checkTriangleParallel 関数の中のスレッドで呼ばれる関数。
void checkTriangleParallelSub(const vector<Configuration> &confs, const vector<string> &basenames, int i, int j, 
    const vector<Configuration> &reducibles, const string &output_dirname) {
    try {
        size_t N = confs.size();
        vector<tuple<NearTriangulation, int, int>> adjacent_near_triangulations = 
            adjacentCenterInternal(confs[i].nearTriangulation(), confs[j].nearTriangulation(), confs[i].ringSize(), confs[j].ringSize(), reducibles);
        // 2 個かさねただけで 0 通りになるものは 3 重目のループに入らない。
        if (adjacent_near_triangulations.size() == 0) return;
        boost::asio::io_service io_service;
        // 並列化しすぎると 1 つのプロセスが持てるスレッドの最大値を超えるから 3 にしておく。
        unsigned int num_concurrent = 3;
        ThreadPool pool(io_service, num_concurrent);
        int count = 0;
        for (size_t l = 0;l < adjacent_near_triangulations.size(); l++) {
            const auto &[adjacent_near_triangulation, v1, v2] = adjacent_near_triangulations[l];
            auto e = std::make_pair(v1, v2);
            for (int v_triangle : adjacent_near_triangulation.diagonalVertices().at(e)) {
                for (size_t k = j;k < N; k++) {
                    auto vs = std::make_pair(v_triangle, confs[k].ringSize());
                    auto basename = std::make_pair(std::to_string(count) + "_" + basenames[i] + "_" + basenames[j], basenames[k]);
                    const NearTriangulation &near_triangulation3 = confs[k].nearTriangulation();
                    pool.post(boost::bind(correspondCenter, boost::cref(adjacent_near_triangulation), boost::cref(near_triangulation3), 
                        vs, basename, boost::cref(reducibles), boost::cref(output_dirname)));
                    count++;
                }
            }
        }
    } catch (std::exception &e) {
        spdlog::critical("An error ({}) occured when combining {}, {}", e.what(), basenames[i], basenames[j]);
        std::ofstream ofs(fmt::format("{}/{}_{}.err", output_dirname, basenames[i], basenames[j]));
        ofs << e.what();
        return;
    }
    return;
}

// confs_dirname の configuration たちに含まれる 3 つの configuration を checkTriangle と同様にタイリングした結果、
// reducible_dirname の configuration を含まないようなものを返す。
// 並列化している。
void checkTriangleParallel(const string &confs_dirname, const string &reducibles_dirname, const string &output_dirname) {
    vector<Configuration> confs = getConfs(confs_dirname, true);
    vector<Configuration> reducibles = getConfs(reducibles_dirname);
    bool madedir = fs::create_directory(output_dirname);
    if (madedir) spdlog::info("made {} directory", output_dirname);
    size_t N = confs.size();
    vector<string> basenames(N);
    for (size_t i = 0;i < N; i++) {
        basenames[i] = baseName(confs[i].fileName());
    }
    spdlog::info("start creating thread pool");
    boost::asio::io_service io_service;
    unsigned int num_concurrent = boost::thread::hardware_concurrency() / 2;
    ThreadPool pool(io_service, num_concurrent);
    spdlog::info("finish creating thread pool (the number of thread is {})", num_concurrent);
    spdlog::info("start queueing thread into pool");
    for (size_t i = 0;i < N; i++) {
        for (size_t j = i;j < N; j++) {
            pool.post(boost::bind(checkTriangleParallelSub, boost::cref(confs), boost::cref(basenames), i, j, boost::cref(reducibles), boost::cref(output_dirname)));
        }
    }
    spdlog::info("finish queueing thread into pool");
    return;
}


// checkAngleParallel 関数の中のスレッドで呼ばれる関数。
void checkAngleParallelSub(const vector<Configuration> &confs1, const vector<string> &basenames1,
    const vector<Configuration> &confs2, const vector<string> &basenames2, int i, int j,
    const vector<Configuration> &reducibles, const string &output_dirname, int angle) {
    try {
        size_t N2 = confs2.size();
        vector<tuple<NearTriangulation, int, int>> adjacent_near_triangulations = 
            adjacentCenterInternal(confs1[i].nearTriangulation(), confs2[j].nearTriangulation(), confs1[i].ringSize(), confs2[j].ringSize(), reducibles);
        // 2 個かさねただけで 0 通りになるものは 3 重目のループに入らない。
        if (adjacent_near_triangulations.size() == 0) return;
        boost::asio::io_service io_service;
        // 並列化しすぎると 1 つのプロセスが持てるスレッドの最大値を超えるから 3 にしておく。
        unsigned int num_concurrent = 3;
        ThreadPool pool(io_service, num_concurrent);
        int count = 0;
        for (size_t l = 0;l < adjacent_near_triangulations.size(); l++) {
            const auto &[adjacent_near_triangulation, v1, v2] = adjacent_near_triangulations[l];
            auto e = std::make_pair(v1, v2);
            for (int v_triangle : adjacent_near_triangulation.diagonalVertices().at(e)) {
                // v1 を軸として angle 回回転する
                int v_angle_old = v2;
                int v_angle = v_triangle;
                for (int a = 0;a < angle; a++) {
                    auto ea = std::make_pair(v1, v_angle);
                    const auto &diagonal_ea = adjacent_near_triangulation.diagonalVertices().at(ea);
                    assert(diagonal_ea.size() == 2);
                    int v_angle_new = diagonal_ea[0] == v_angle_old ? diagonal_ea[1] : diagonal_ea[0];
                    assert(v_angle_old != v_angle_new);
                    v_angle_old = v_angle;
                    v_angle = v_angle_new;
                }
                for (size_t k = j;k < N2; k++) {
                    auto vs = std::make_pair(v_angle, confs2[k].ringSize());
                    auto basename = std::make_pair(std::to_string(count) + "_" + basenames1[i] + "_" + basenames2[j], basenames2[k]);
                    const NearTriangulation &near_triangulation3 = confs2[k].nearTriangulation();
                    pool.post(boost::bind(correspondCenter, boost::cref(adjacent_near_triangulation), boost::cref(near_triangulation3), 
                        vs, basename, boost::cref(reducibles), boost::cref(output_dirname)));
                    count++;
                }
            }
        }
    } catch (std::exception &e) {
        spdlog::critical("An error ({}) occured when combining {}, {}", e.what(), basenames1[i], basenames2[j]);
        std::ofstream ofs(fmt::format("{}/{}_{}.err", output_dirname, basenames1[i], basenames2[j]));
        ofs << e.what();
        return;
    }
    return;
}

// confs_dirname1 の conf から 1 つ、 confs_dirname2 の conf から 2 つ の conf を取ってくることで、angle を構成する。具体的には
// ある頂点 v とその近傍の 2 頂点 u1, u2 であって、v の埋め込みの rotation で u1 と u2 が angle 個の頂点を挟んでいるような状態で、
// confs_dirname1 の conf の中心が v, confs_dirname2 の conf の 2 つの中身が u1, u2 であるような場合を試す。
// angle = 0 のときは 3 角形を試していることになる。
void checkAngleParallel(const string &confs_dirname1, const string &confs_dirname2, const string &reducibles_dirname, const string &output_dirname, int angle) {
    vector<Configuration> confs1 = getConfs(confs_dirname1, true);
    vector<Configuration> confs2 = getConfs(confs_dirname2, true);
    vector<Configuration> reducibles = getConfs(reducibles_dirname);
    bool madedir = fs::create_directory(output_dirname);
    if (madedir) spdlog::info("made {} directory", output_dirname);
    size_t N1 = confs1.size();
    vector<string> basenames1(N1);
    for (size_t i = 0;i < N1; i++) {
        basenames1[i] = baseName(confs1[i].fileName());
    }
    size_t N2 = confs2.size();
    vector<string> basenames2(N2);
    for (size_t i = 0;i < N2; i++) {
        basenames2[i] = baseName(confs2[i].fileName());
    }
    spdlog::info("start creating thread pool");
    boost::asio::io_service io_service;
    unsigned int num_concurrent = boost::thread::hardware_concurrency() / 2;
    ThreadPool pool(io_service, num_concurrent);
    spdlog::info("finish creating thread pool (the number of thread is {})", num_concurrent);
    spdlog::info("start queueing thread into pool");
    for (size_t i = 0;i < N1; i++) {
        for (size_t j = 0;j < N2; j++) {
            pool.post(boost::bind(checkAngleParallelSub, boost::cref(confs1), boost::cref(basenames1), boost::cref(confs2), boost::cref(basenames2), i, j, boost::cref(reducibles), boost::cref(output_dirname), angle));
        }
    }
    spdlog::info("finish queueing thread into pool");
    return;
}

int main(const int ac, const char* const* const av) {
    using namespace boost::program_options;
    options_description description("Options");
    description.add_options()
        ("adjacent", "Tile two configurations so that two centers are adjacent")
        ("triangle", "Tile three configurations sot that three centers consist of a triangle")
        ("Angle", "Tile three configurations so that three centers consist of an Angle")
        ("conf1", value<string>(), "A configuration filename")
        ("conf2",  value<string>(), "A configuration filename")
        ("conf3", value<string>(), "A configuration filename")
        ("configurations1", value<string>(), "A directory that contains configuration files")
        ("configurations2", value<string>(), "A directory that contains configuration files")
        ("reducibles,r", value<string>(), "A directory that contains reducible configuration files")
        ("output,o", value<string>(), "Output directory")
        ("angle", value<int>(), "An angle value")
        ("help,H", "Display options");

    variables_map vm;
    store(parse_command_line(ac, av, description), vm);
    notify(vm);

    if (vm.count("help")) {
        description.print(std::cout);
        return 0;
    }

    if (vm.count("adjacent")) {
        if (vm.count("conf1") && vm.count("conf2") && vm.count("reducibles") && vm.count("output")) {
            auto conf1_filename = vm["conf1"].as<string>();
            auto conf2_filename = vm["conf2"].as<string>();
            auto reducibles_dirname = vm["reducibles"].as<string>();
            auto output_dirname = vm["output"].as<string>();
            checkAdjacentCenter(conf1_filename, conf2_filename, reducibles_dirname, output_dirname);
        }
        else if (vm.count("configurations1") && vm.count("configurations2") && vm.count("reducibles") && vm.count("output")) {
            auto confs1_dirname = vm["configurations1"].as<string>();
            auto confs2_dirname = vm["configurations2"].as<string>();
            auto reducibles_dirname = vm["reducibles"].as<string>();
            auto output_dirname = vm["output"].as<string>();
            checkAdjacentCenterParallel1(confs1_dirname, confs2_dirname, reducibles_dirname, output_dirname);
        }
        else if (vm.count("configurations1") && vm.count("reducibles") && vm.count("output")) {
            auto confs_dirname = vm["configurations1"].as<string>();
            auto reducibles_dirname = vm["reducibles"].as<string>();
            auto output_dirname = vm["output"].as<string>();
            checkAdjacentCenterParallel2(confs_dirname, reducibles_dirname, output_dirname);
        }
    } else if (vm.count("triangle")) {
        if (vm.count("configurations1") && vm.count("reducibles") && vm.count("output")) {
            auto confs_dirname = vm["configurations1"].as<string>();
            auto reducibles_dirname = vm["reducibles"].as<string>();
            auto output_dirname = vm["output"].as<string>();
            checkTriangleParallel(confs_dirname, reducibles_dirname, output_dirname);
        }
    } else if (vm.count("Angle")) {
        if (vm.count("configurations1") && vm.count("configurations2") && vm.count("reducibles") && vm.count("output") && vm.count("angle")) {
            // configurations1 も configurations2 も指定されている場合は triangle を作るときに configurations1 から
            // 1 つ、configurations2 から 2 つ configuration をとってきて、angle を作る。
            auto confs_dirname1 = vm["configurations1"].as<string>();
            auto confs_dirname2 = vm["configurations2"].as<string>();
            auto reducibles_dirname = vm["reducibles"].as<string>();
            auto output_dirname = vm["output"].as<string>();
            int angle = vm["angle"].as<int>();
            checkAngleParallel(confs_dirname1, confs_dirname2, reducibles_dirname, output_dirname, angle);
        }
    }
    spdlog::info("finish executing this program!");

    return 0;
}