#!/bin/bash

#
# The script is used to make sure the center of each wheel has charge at most 0 or -1.
# When center's degree is at most 9:
#   If the hub of all wheels has charge at most 0, the message "All finished" is displayed.
# When center's degree is at least 10:
#   If the hub of all wheels has charge at most -1, the message "All finished" is displayed.
#
# 
# Usage)
# bash charge_result.sh torus <The degree of the hub> <The output file>
#
# Example)
# bash charge_result.sh torus 7 out_torus7.txt
# 

set -euxo pipefail
cd $(dirname $0)

if [ $# -ne 3 ]; then
    echo -e "\e[31merror:\e[m Please follow the usage 'bash charge_result.sh torus <The degree of the hub> <The output file>'"
    exit 1
fi


if [ "$1" = "torus" ]; then
    num_wheel=(0 0 0 0 0 0 0 5401 8334 8912 6796 3814)
    finish_sum=0
    sum=0
    for i in $(seq 0 ${num_wheel["$2"]}); do
        finished=$(grep 'the ratio of overcharged cartwheel' "./torus/log/$2_${i}_t-1.wheel.log" | wc -l) && true
        n_charge_eq0=$(grep -E "charge \(initial, receive, send, result\) : -?[0-9]+, [0-9]+, [0-9]+, 0" "./torus/log/$2_${i}_t-1.wheel.log" | wc -l) && true
        n_charge_ge0=$(grep -E "charge \(initial, receive, send, result\) : -?[0-9]+, [0-9]+, [0-9]+, [0-9]+" "./torus/log/$2_${i}_t-1.wheel.log" | wc -l) && true
        all_neighbor_56=$(grep -E "^$2( (5|6))+$" "./torus/wheel/$2_${i}.wheel" | wc -l) && true
        if [ $finished -ne 1 ]; then
            echo "$i has not finished" >> "$3"
            finish_sum=$(($finish_sum+1))
        fi
        if [ $2 -le 9 ] && [ $all_neighbor_56 -eq 0 ]; then
            # when degree <= 9 and not all neighbors' degree is 5 or 6.
            # detect positive charge
            if [ $n_charge_eq0 -ne $n_charge_ge0 ]; then
                echo "$i has overcharged cartwheel" >> "$3"
                sum=$(($sum+$n_charge_ge0-$n_charge_eq0))
            fi
        else
            # when degree <= 9 and all neighbors' degree is 5 or 6,
            #    or degree >= 10,
            # detect positive or zero charge
            if [ $n_charge_ge0 -ne 0 ]; then
                echo "$i has overcharged cartwheel" >> "$3"
                sum=$(($sum+$n_charge_ge0))
            fi
        fi
    done
    if [ $finish_sum -eq 0 ] && [ $sum -eq 0 ]; then
        echo "All finished!"
    else
        echo "$finish_sum cartwheel have not finished"
        echo "$sum cartwheel remained"
    fi
fi

