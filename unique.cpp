
#include "configuration.hpp"
#include "basewheel.hpp"
#include <spdlog/spdlog.h>
#include <boost/program_options.hpp>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

bool confContainOneOfConfs(const vector<Configuration> &confs, const Configuration conf) {

    const string &conffile = conf.fileName();
    for (const auto &confi : confs) {
        const string &confifile = confi.fileName();
        if (BaseWheel::containConf(conf.nearTriangulation(), confi)) {
            spdlog::info("{} contains {} as subgraph, you shouldn't include {}", conffile, confifile, conffile);
            return true;
        }
    }

    spdlog::info("You should add {} if {} is reducible", conffile, conffile);
    return false;
}

bool confEqualOneOfConfs(const vector<Configuration> &confs, const Configuration conf) {

    const string &conffile = conf.fileName();
    for (const auto &confi : confs) {
        const string &confifile = confi.fileName();
        if (BaseWheel::containConf(conf.nearTriangulation(), confi)
         && BaseWheel::containConf(confi.nearTriangulation(), conf)) {
            spdlog::info("{} equals {}, you shouldn't include {}", conffile, confifile, conffile);
            return true;
        }
    }

    spdlog::info("You should add {} if {} is reducible", conffile, conffile);
    return false;
}

void uniquConf(const vector<Configuration> &confs) {

    spdlog::info("Fisrt, You have {} confs", confs.size());

    vector<Configuration> confs1;
    for (const auto& conf : confs) {
        if (!confContainOneOfConfs(confs1, conf)) {
            confs1.push_back(conf);
        }
    }

    std::reverse(confs1.begin(), confs1.end());
    vector<Configuration> confs2;
    for (const auto& conf : confs1) {
        if (!confContainOneOfConfs(confs2, conf)) {
            confs2.push_back(conf);
        }
    }

    spdlog::info("Now, You have {} confs (List follows)", confs2.size());
    for (const auto &conf : confs2) {
        spdlog::info("{}", conf.fileName());
    }

    return;
}

void removeConf(const vector<Configuration> &confs1, const vector<Configuration> &confs2) {
    
    spdlog::info("Fisrt, You have {} confs1", confs1.size());

    vector<string> filenames;
    for (const auto& conf1 : confs1) {
        if (!confEqualOneOfConfs(confs2, conf1)) {
            filenames.push_back(conf1.fileName());
        }
    }

    spdlog::info("Now, You have {} confs (List follows)", filenames.size());
    for (const auto &filename: filenames) {
        spdlog::info("{}", filename);
    }

    return;
}

int main(const int ac, const char* const* const av) {
    using namespace boost::program_options;
    options_description description("Options");
    description.add_options()
        ("dirname1,d", value<string>(), "The dirname that contains conf file")
        ("dirname2,D", value<string>(), "The dirname that contains conf file")
        ("conf,c",  value<string>(), "The configuration file")
        ("help,H", "Display options");

    variables_map vm;
    store(parse_command_line(ac, av, description), vm);
    notify(vm);

    if (vm.count("help")) {
        description.print(std::cout);
        return 0;
    }

    if (vm.count("dirname1") && vm.count("conf")) {
        auto dirname = vm["dirname1"].as<string>();
        auto conf_filename = vm["conf"].as<string>();
        vector<Configuration> confs = getConfs(dirname);
        Configuration conf = Configuration::readConfFile(conf_filename);
        confContainOneOfConfs(confs, conf);
        return 0;
    } else if (vm.count("dirname1") && vm.count("dirname2")) {
        auto dirname1 = vm["dirname1"].as<string>();
        auto dirname2 = vm["dirname2"].as<string>();
        vector<Configuration> confs1 = getConfs(dirname1);
        vector<Configuration> confs2 = getConfs(dirname2);
        removeConf(confs1, confs2);
    } else if (vm.count("dirname1")) {
        auto dirname = vm["dirname1"].as<string>();
        vector<Configuration> confs = getConfs(dirname);
        uniquConf(confs);
        return 0;
    } 

    return 0;
}