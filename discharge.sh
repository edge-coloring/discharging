#!/bin/bash

# The script is used to execute the discharging procedure.
# We specify the degree of the hub(=: d), the smaller index of the range (=: l), the larger index of the range(=: r), 
# the directory that contains rule files, the directory that contains configuration files.
# Then, the script executes the discharging procedure to ./torus/wheel/d_l.wheel, ./torus/wheel/d_{l+1}.wheel ... ./torus/wheel/d_r.wheel
# The log files (e.g. 7_0.wheel.log) are placed in ./torus_log directory.
#
# Usage)
# bash discharge.sh torus <The degree of the hub> <The smaller index of the range> <The larger index of the range> <The directory that contains rule files> <The directory that contains configuration files>
# 
# Example)
# bash discharge.sh torus 7 0 1500 toroidal_configurations/rule toroidal_configurations/reducible/conf
#
#

set -euxo pipefail
cd $(dirname $0)

if [ $# -ne 6 ]; then
    echo -e "\e[31merror:\e[m Please follow the usage 'bash discharge.sh torus <The degree of the hub> <The smaller index of the range> <The larger index of the range> <The directory that contains rule files> <The directory that contains configuration files>'"
    exit 1
fi

l=$3
r=$4
rule=$5
conf=$6


if [ "$1" = "torus" ]; then
    #
    # degree 7: we need to execute ./torus/wheel/7_{0..5401}.wheel
    # degree 8: we need to execute ./torus/wheel/8_{0..8334}.wheel
    # degree 9: we need to execute ./torus/wheel/9_{0..8912}.wheel
    # degree 10: we need to execute ./torus/wheel/10_{0..6796}.wheel
    # degree 11: we need to execute ./torus/wheel/11_{0..3814}.wheel
    #
    threshold=-1
    send="./torus/send"
    mkdir -p torus/log
    for i in $(seq $l $r); do
        if [ -e "./torus/log/$2_${i}_t${threshold}.wheel.log" ]; then
            # log file exists
            n_charge_eq0=$(grep -E "charge \(initial, receive, send, result\) : -?[0-9]+, [0-9]+, [0-9]+, 0" "./torus/log/$2_${i}_t${threshold}.wheel.log" | wc -l) && true
            n_charge_ge0=$(grep -E "charge \(initial, receive, send, result\) : -?[0-9]+, [0-9]+, [0-9]+, [0-9]+" "./torus/log/$2_${i}_t${threshold}.wheel.log" | wc -l) && true
            all_neighbor_56=$(grep -E "^$2( (5|6))+$" "./torus/wheel/$2_${i}.wheel" | wc -l) && true
            if [ $2 -le 9 ] && [ $all_neighbor_56 -eq 0 ]; then
                # when degree <= 9 and not all neighbors' degree is 5 or 6.
                # detect positive charge
                if [ $n_charge_eq0 -ne $n_charge_ge0 ]; then
                    ./build/a.out -w "./torus/wheel/$2_$i.wheel" -r "$rule" -c "$conf" -s "$send" -m 9 -t "$threshold" -v 1 > "./torus/log/$2_${i}_t${threshold}.wheel.log" &
                fi
            else
                # when degree <= 9 and all neighbors' degree is 5 or 6,
                #    or degree >= 10,
                # detect positive or zero charge
                if [ $n_charge_ge0 -ne 0 ]; then
                    ./build/a.out -w "./torus/wheel/$2_$i.wheel" -r "$rule" -c "$conf" -s "$send" -m 9 -t "$threshold" -v 1 > "./torus/log/$2_${i}_t${threshold}.wheel.log" &
                fi
            fi
        else
            ./build/a.out -w "./torus/wheel/$2_$i.wheel" -r "$rule" -c "$conf" -s "$send" -m 9 -t "$threshold" -v 1 > "./torus/log/$2_${i}_t${threshold}.wheel.log" &
        fi
    done
fi
