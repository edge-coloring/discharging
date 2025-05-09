#!/bin/bash

#
# The script is used to enumerate wheels.
# We specify the degree of the hub and the directory that contains configuration files.
# The log file (e.g. torus_deg7.log, ...) is placed in ./log directory. 
# The results (the wheel files) will be placed in ./torus/wheel directory.
#
# Usage)
# bash enum_wheel.sh torus <The degree of the hub> <The directory that contains configurations>
#
# Example)
# bash enum_wheel.sh torus 7 toroidal_configurations/reducible/conf
# 

set -euxo pipefail
cd $(dirname $0)

if [ $# -ne 3 ]; then
    echo -e "\e[31merror:\e[m Please follow the usage 'bash enum_wheel.sh torus <The degree of the hub> <The directory that contains configurations>'"
    exit 1
fi

deg=$2
conf=$3

mkdir -p log

if [ "$1" = "torus" ]; then
    send="./torus/send"
    wheel="./torus/wheel"
    mkdir -p "$wheel"

    ./build/a.out -d "$deg" -c "$conf" -s "$send" -m 9 -o "$wheel" -t -1 -v 1 > ./log/torus_deg$deg.log &
fi