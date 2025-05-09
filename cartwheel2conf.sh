#!/bin/bash

#
# The script is used to get configurations that represents cartwheels of final charge 0
# We specify the directory that contains wheel files and the directory that contains discharging's log files.
# The results are placed in appropriate named directories.
#
# Usage)
# bash cartwheel2conf.sh torus <The directory that contains wheel files> <The directory that contains discharging's log files>
#
# Example)
# bash cartwheel2conf.sh torus torus/wheel torus/log
#

set -euo pipefail
cd $(dirname $0)

if [ $# -ne 3 ]; then
    echo -e "\e[31merror:\e[m Please follow the usage 'bash cartwheel2conf.sh torus <The directory that contains wheel files> <The directory that contains discharging's log files>'"
    exit 1
fi

wheeldir=$2
logdir=$3

if [ "$1" = "torus" ]; then
    num_wheel=(0 0 0 0 0 0 0 5401 8334 8912)

    # degree 9
    mkdir -p torus_deg9 # 99, 89, 79

    for i in $(seq 0 ${num_wheel[9]}); do
        logFile="${logdir%/}/9_${i}_t-1.wheel.log"
        wheelFile="${wheeldir%/}/9_${i}.wheel"
        python3 cartwheel2conf.py $logFile torus_deg9 --prefix "9_${i}_" --remain &
    done


    # degree 8
    mkdir -p torus_deg8_9 # 89
    mkdir -p torus_deg8_8_no_9 # 88
    mkdir -p torus_deg8_77_no_9_8 # 778
    mkdir -p torus_deg8_7ge3_no_9_8_77 # 8_7ge3
    mkdir -p torus_deg8_7_no_9_8_77_7le2 # 7_86m8, 78

    for i in $(seq 0 ${num_wheel[8]}); do
        logFile="${logdir%/}/8_${i}_t-1.wheel.log"
        wheelFile="${wheeldir%/}/8_${i}.wheel"
        cnt9=$(grep -oi '.9' "$wheelFile" | wc -l) && true
        cnt8=$(grep -oi '.8' "$wheelFile" | wc -l) && true
        cnt7=$(grep -oi '.7' "$wheelFile" | wc -l) && true
        cnt77=$(grep -oi -E '(.7 7)|(^8 7 (.*) 7$)' "$wheelFile" | wc -l) && true
        if [ $cnt9 -ne 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg8_9 --prefix "8_${i}_" --remain &
        fi
        if [ $cnt8 -ne 0 ] && [ $cnt9 -eq 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg8_8_no_9 --prefix "8_${i}_" --remain &
        fi
        if [ $cnt77 -ne 0 ] && [ $cnt9 -eq 0 ] && [ $cnt8 -eq 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg8_77_no_9_8 --prefix "8_${i}_" --remain &
        fi
        if [ $cnt7 -ge 3 ] && [ $cnt9 -eq 0 ] && [ $cnt8 -eq 0 ] && [ $cnt77 -eq 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg8_7ge3_no_9_8_77 --prefix "8_${i}_" --remain &
        fi
        if [ $cnt7 -ne 0 ] && [ $cnt9 -eq 0 ] && [ $cnt8 -eq 0 ] && [ $cnt77 -eq 0 ] && [ $cnt7 -le 2 ]; then
            python3 cartwheel2conf.py $logFile torus_deg8_7_no_9_8_77_7le2 --prefix "8_${i}_" --remain &
        fi
    done


    # degree 7
    mkdir -p torus_deg7_9 # 79
    mkdir -p torus_deg7_78_no_9_88 # 778
    mkdir -p torus_deg7_77_no_9_88_78 # 777
    mkdir -p torus_deg7_8_no_9_88_78_77 # 8_7ge3
    mkdir -p torus_deg7_86m8_no_9_88_78_77 # 7_86m8
    mkdir -p torus_deg7_8_no_9_88_78_77_86m8 # 78
    mkdir -p torus_deg7_7ge3_no_9_77_8 # 7_7ge3
    mkdir -p torus_deg7_7_no_9_77_8 # 7_7ge3
    mkdir -p torus_deg7_7_no_9_77_8_7le2 # 77

    for i in $(seq 0 ${num_wheel[7]}); do
        logFile="${logdir%/}/7_${i}_t-1.wheel.log"
        wheelFile="${wheeldir%/}/7_${i}.wheel"
        cnt9=$(grep -oi '.9' "$wheelFile" | wc -l) && true
        cnt8=$(grep -oi '.8' "$wheelFile" | wc -l) && true
        cnt7=$(grep -oi '.7' "$wheelFile" | wc -l) && true
        cnt78=$(grep -oi -E '(.7 8)|(.8 7)|(^7 7 (.*) 8$)|(^7 8 (.*) 7$)' "$wheelFile" | wc -l) && true
        cnt77=$(grep -oi -E '(.7 7)|(^7 7 (.*) 7$)' "$wheelFile" | wc -l) && true
        cnt88=$(grep -oi -E '(.8 8)|(^7 8 (.*) 8$)' "$wheelFile" | wc -l) && true
        cnt86m8=$(grep -oi -E '(.8 (5|6) 8)|(^7 8 (.*) 8 (5|6)$)|(^7 (5|6) 8 (.*) 8$)' "$wheelFile" | wc -l) && true
        if [ $cnt9 -ne 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg7_9 --prefix "7_${i}_" --remain &
        fi
        if [ $cnt78 -ne 0 ] && [ $cnt9 -eq 0 ] && [ $cnt88 -eq 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg7_78_no_9_88 --prefix "7_${i}_" --remain &
        fi
        if [ $cnt77 -ne 0 ] && [ $cnt9 -eq 0 ] && [ $cnt88 -eq 0 ] && [ $cnt78 -eq 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg7_77_no_9_88_78 --prefix "7_${i}_" --remain &
        fi
        if [ $cnt8 -ne 0 ] && [ $cnt9 -eq 0 ] && [ $cnt88 -eq 0 ] && [ $cnt78 -eq 0 ] && [ $cnt77 -eq 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg7_8_no_9_88_78_77 --prefix "7_${i}_" --remain &
        fi
        if [ $cnt86m8 -ne 0 ] && [ $cnt9 -eq 0 ] && [ $cnt88 -eq 0 ] && [ $cnt78 -eq 0 ] && [ $cnt77 -eq 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg7_86m8_no_9_88_78_77 --prefix "7_${i}_" --remain &
        fi
        if [ $cnt8 -ne 0 ] && [ $cnt9 -eq 0 ] && [ $cnt88 -eq 0 ] && [ $cnt78 -eq 0 ] && [ $cnt77 -eq 0 ] && [ $cnt86m8 -eq 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg7_8_no_9_88_78_77_86m8 --prefix "7_${i}_" --remain &
        fi
        if [ $cnt7 -ge 3 ] && [ $cnt9 -eq 0 ] && [ $cnt77 -eq 0 ] && [ $cnt8 -eq 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg7_7ge3_no_9_77_8 --prefix "7_${i}_" --remain &
        fi
        if [ $cnt7 -ne 0 ] && [ $cnt9 -eq 0 ] && [ $cnt77 -eq 0 ] && [ $cnt8 -eq 0 ]; then
            python3 cartwheel2conf.py $logFile torus_deg7_7_no_9_77_8 --prefix "7_${i}_" --remain &
        fi
        if [ $cnt7 -ne 0 ] && [ $cnt9 -eq 0 ] && [ $cnt77 -eq 0 ] && [ $cnt8 -eq 0 ] && [ $cnt7 -le 2 ]; then
            python3 cartwheel2conf.py $logFile torus_deg7_7_no_9_77_8_7le2 --prefix "7_${i}_" --remain &
        fi
    done


fi
