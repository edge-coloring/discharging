#!/bin/bash

# 
# The script is used to tile configurations of final charge 0
# We have to use this script after executing cartwheel2conf.sh and getting configuration files 
# that represent cartwheels whose center vertex's charge is exactly zero.
#
# We have to check the following cases:
# 1) An edge uv such that d(u) = d(v) = 9
# 2) An edge uv such that d(u) = 8, d(v) = 9
# 3) An edge uv such that d(u) = 7, d(v) = 9
# 4) An edge uv such that d(u) = 8, d(v) = 8
# 5) A facial triangle uvw such that d(u) = d(v) = 7, d(w) = 8
# 6) A facial triangle uvw such that d(u) = d(v) = d(w) = 7
# 7) A vertex v such that d(v) = 8, and at least 3 neighbors of v has degree 7
# 8) A vertex v such that d(v) = 7, and clockwise consecutive neighbors u1, u2, u3 of v satisfy d(u1) = d(u3) = 8, d(u2) <= 6
# 9) An edge uv such that d(u) = 7, d(v) = 8
# 10) A vertex v such that d(v) = 7, and at least 3 neighbors of v has degree 7
# 11) An edge uv such that d(u) = d(v) = 7
# 
# CASE = 99|89|79|88|778|777|8_7ge3|7_86m8|78|7_7ge3|77
#
# By combining cartwheels, we find a reducible configuration in case 1,2,3,...,10, so
# no result appears by executing programs in these cases (no result means all combined cartwheels contain a reducilbe configurations.).
# However we have several results in case 11. these combined cartwheels correpsond to graphs represented in Lemma 9.2 
# (we do not make them isomorphic, so the number is different.)
#
# Usage)
# bash tile.sh torus <The directory that contains reducible configuration files> <CASE>
#
# Example)
# bash tile.sh torus toroidal_configurations/reducible/conf 99
#

set -euxo pipefail
cd $(dirname $0)

if [ $# -ne 3 ]; then
    echo -e "\e[31merror:\e[m Please follow the usage 'bash tile.sh torus <The directory that contains reducible configuration files> <CASE>'"
    exit 1
fi 

mkdir -p log
reducibles=$2
CASE=$3

if [ "$1" = "torus" ]; then
    if [ $CASE = '99' ]; then
        ./build/tile --adjacent --configurations1 torus_deg9 -r $reducibles -o torus_combined99 > log/torus_combined99.log &
    fi
    if [ $CASE = '89' ]; then
        ./build/tile --adjacent --configurations1 torus_deg8_9 --configurations2 torus_deg9 -r $reducibles -o torus_combined89 > log/torus_combined89.log &
    fi
    if [ $CASE = '79' ]; then
        ./build/tile --adjacent --configurations1 torus_deg7_9 --configurations2 torus_deg9 -r $reducibles -o torus_combined79 > log/torus_combined79.log &
    fi
    if [ $CASE = '88' ]; then
        ./build/tile --adjacent --configurations1 torus_deg8_8_no_9 -r $reducibles -o torus_combined88 > log/torus_combined88.log &
    fi
    if [ $CASE = '778' ]; then
        ./build/tile --Angle --configurations1 torus_deg8_77_no_9_8 --configurations2 torus_deg7_78_no_9_88 --angle 0 -r $reducibles -o torus_combined778 > log/torus_combined778.log &
    fi
    if [ $CASE = '777' ]; then
        ./build/tile --triangle --configurations1 torus_deg7_77_no_9_88_78 -r $reducibles -o torus_combined777 > log/torus_combined777.log &
    fi
    if [ $CASE = '8_7ge3' ]; then
        ./build/tile --Angle --configurations1 torus_deg8_7ge3_no_9_8_77 --configurations2 torus_deg7_8_no_9_88_78_77 --angle 1 -r $reducibles -o torus_combined8_7ge3 > log/torus_combined8_7ge3.log &
    fi
    if [ $CASE = '7_86m8' ]; then
        ./build/tile --Angle --configurations1 torus_deg7_86m8_no_9_88_78_77 --configurations2 torus_deg8_7_no_9_8_77_7le2 --angle 1 -r $reducibles -o torus_combined7_86m8 > log/torus_combined7_86m8.log &
    fi
    if [ $CASE = '78' ]; then
        ./build/tile --adjacent --configurations1 torus_deg8_7_no_9_8_77_7le2 --configurations2 torus_deg7_8_no_9_88_78_77_86m8 -r $reducibles -o torus_combined78 > log/torus_combined78.log &
    fi
    if [ $CASE = '7_7ge3' ]; then
        ./build/tile --Angle --configurations1 torus_deg7_7ge3_no_9_77_8 --configurations2 torus_deg7_7_no_9_77_8 --angle 1 -r $reducibles -o torus_combined7_7ge3 > log/torus_combined7_7ge3.log &
    fi
    if [ $CASE = '77' ]; then
        ./build/tile --adjacent --configurations1 torus_deg7_7_no_9_77_8_7le2 -r $reducibles -o torus_combined77 > log/torus_combined77.log &
    fi
fi

