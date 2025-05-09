#pragma once
#include <string>
#include "near_triangulation.hpp"

using std::string;

class Configuration {
private:
    NearTriangulation conf_;
    int ring_size_;
    int inside_edge_id_;
    bool ring_remains_;
    string filename_;
    
public:
    Configuration(int ring_size, bool ring_remains, const string &filename, const NearTriangulation &conf);
    static Configuration readConfFile(const string &filename, bool ring_remains=false);
    void writeConfFile(const string &filename) const;

    const NearTriangulation &nearTriangulation(void) const;
    int ringSize(void) const;
    bool ringRemains(void) const;
    const string &fileName(void) const;
    int diameter(void) const;
    
    int getInsideEdgeId(void) const;
};

optional<Configuration> convertNearTriangulation2Conf(NearTriangulation &near_triangulation);
vector<Configuration> getConfs(const std::string &dirname, bool ring_remains=false);

