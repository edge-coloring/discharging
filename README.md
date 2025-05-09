# Discharging

This program is used to prove Lemma regarding the discharging procedure. Specifically, it proves Lemma 4.5, 4.7, 4.11 in Section 4 and Lemma 9.1, 9.2, 9.4 in Section 9.

## Requirements
+ g++, CMake, Boost, spdlog

## Build
```bash
cmake -S . -B build
cmake --build build
```

## Execution
### Lemma 4.7, 4.11
We enumerate the cases where a vertex sends a charge.
We prepared a shell script (```enum_send.sh```), so we only execute the command below. 
```bash
bash enum_send.sh torus toroidal_configurations/rule toroidal_configurations/reducible/conf
```
More detailed information is in ```enum_send.sh```.

### Lemma 4.5, 9.4 (also a part of Lemma 9.1)
We use the results of ```enum_send.sh``` to prove Lemma 4.5, 9.4, and a part of Lemma 9.1 (these Lemmas are prove at the same time), so first we have to execute the above command. After that, we need three steps to get the results.

1. enumerate wheels.\
We prepared a shell script (```enum_wheel.sh```), so we only execute the commands below. It takes about 3 hours to generate wheel whose center vertex's degree is 11.
```bash
bash enum_wheel.sh torus 7 toroidal_configurations/reducible/conf
bash enum_wheel.sh torus 8 toroidal_configurations/reducible/conf
bash enum_wheel.sh torus 9 toroidal_configurations/reducible/conf
bash enum_wheel.sh torus 10 toroidal_configurations/reducible/conf
bash enum_wheel.sh torus 11 toroidal_configurations/reducible/conf
```
More detailed information is in ```enum_wheel.sh```.

2. exectute the discharging procedure to each wheel.\
We prepared a shell script (```discharge.sh```).
There are many wheel files in ```torus/wheel```, which are enumerated in step 1, so we cannot execute the discharging procedure to all wheels at once. We run the shell script in several batches. An example is the following.
```bash
bash discharge.sh torus 7 0 1500 toroidal_configurations/rule toroidal_configurations/reducible/conf
bash discharge.sh torus 7 1501 3000 toroidal_configurations/rule toroidal_configurations/reducible/conf
bash discharge.sh torus 7 3001 4500 toroidal_configurations/rule toroidal_configurations/reducible/conf
bash discharge.sh torus 7 4501 5401 toroidal_configurations/rule toroidal_configurations/reducible/conf
```
We run the shell script when the degree is 8,9,10,11 in the similar way. More detailed information is in ```discharge.sh```

3. make sure the charge of each wheel is at most 0 or less than 0.\
We check the followings for each wheel.
+ When all neighbors of a vertex of degree at most 9 have degree of 5 or 6, we check the final charge of it is less than 0 (Lemma 9.4).
+ When some neighbor of a vertex of degree at most 9 has degree of at least 7, we check the final charge of it is less than or equal 0 (Lemma 4.5).
+ When the degree of a vertex is at least 10, we check the final charge of it is less than 0. (case 1 in Lemma 9.1).\
\
We prepared the shell script (```charge_result.sh```) to check the above claims, so we have only to execute the commands below (After executing step 2 for all wheels).
```bash
bash charge_result.sh torus 7 out_torus7.txt
bash charge_result.sh torus 8 out_torus8.txt
bash charge_result.sh torus 9 out_torus9.txt
bash charge_result.sh torus 10 out_torus10.txt
bash charge_result.sh torus 11 out_torus11.txt
```
More detailed information is in ```charge_result.sh```.

### Lemma 9.1, 9.2
We need 2 steps to check Lemma 9.1, 9.2.

1. Get cartwheels whose center vertex's final charge is 0.\
We need to get cartwheels whose center vertex's final charge is 0. To do that, we need to execute the above commands (step 1, 2 for Lemma 4.5, 9.4 (also a part of Lemma 9.1)) in advance.
We prepared sheel script (```cartwheel2conf.sh```). We only have to run the command below. The results are placed in appropriate named directories (```torus_deg(7|8|9)_*```). More detaied information is in ```cartwheel2conf.sh```.
```bash
bash cartwheel2conf.sh torus torus/wheel torus/log
```

2. Tile several cartwheels and find reducible configuration.\
We tile cartwheels that are obtained in step 1.
We prepared sheel script (```tile.sh```). 
We check eleven cases. This is done by running the following commands. More detailed information is in ```tile.sh```.

```bash
bash tile.sh torus toroidal_configurations/reducible/conf 99
bash tile.sh torus toroidal_configurations/reducible/conf 89
bash tile.sh torus toroidal_configurations/reducible/conf 79
bash tile.sh torus toroidal_configurations/reducible/conf 88
bash tile.sh torus toroidal_configurations/reducible/conf 778
bash tile.sh torus toroidal_configurations/reducible/conf 777
bash tile.sh torus toroidal_configurations/reducible/conf 8_7ge3
bash tile.sh torus toroidal_configurations/reducible/conf 7_86m8
bash tile.sh torus toroidal_configurations/reducible/conf 78
bash tile.sh torus toroidal_configurations/reducible/conf 7_7ge3
bash tile.sh torus toroidal_configurations/reducible/conf 77
```
Note that this program is executed in parallel, so do not run the above command at the same time. We executed this program in 256 multicore environments. Multi core environment is recommended to execute this program. The program finishes in a few hours in our environment.

The results are placed in appropriate named directories (```torus_combined*```).
If there are no files in directories ```torus_combined*``` other than ```torus_combined77``` (which corresponds to Lemma 9.2), this program verifies Lemma 9.1, 9.2. More detaied information is in ```tile.sh```.

## Results
The directory ```torus/send```,```torus/wheel``` already contains the results. All cases that a vertex sends a charge are enumearted in ```torus/send```. The drawings of them are also placed in ```torus/send_pdf``` (the cases that a vertex sends charge $N$ are drawen in ```torus/send_pdf/sendN.pdf```).  All wheels are enumerated in ```torus/wheel```. The directory ```torus/log``` contains one of the results (the result of ```./torus/wheel/7_0.wheel```). The result shows the successful log of discharging check. 
If the final charge of center vertex of cartwheel is at least 0, log (```overcharged cartwheel (for machine) : ...```) appears. 

## File Format
### Configuration File Format
A file whose extension is ```.conf``` represents a configuration. The symbol ```N, R``` denotes the size of vertices of the free completion with the ring, the size of the ring respectively. ```v_{R+1}, v_{R+2}, ..., v_N``` denote the vertices of a configuration.

```
dummy-line
N R
R+1 (The degree of v_{R+1}) (The list of neighbors of v_{R+1}) 
R+2 (The degree of v_{R+2}) (The list of neighbors of v_{R+2})
...
N (The degree of v_N) (The list of neighbors of v_N)
```
An example.
```
001
10 6
7 5 2 8 9 10 1
8 5 2 3 4 9 7
9 5 8 4 5 10 7
10 5 9 5 6 1 7
```

### Wheel File Format
A file whose extension is ```.wheel``` represents a wheel.
```v_1, v_2, ...``` denotes the neighbors of the hub in clockwise order.

```
(The degree of the hub) (The degree of v_1) (The degree of v_2) ... 
```
An example
```
7 5 6 6 5 6 8+ 9+
```

### Rule File Format
A file whose extension is ```.rule``` represents a rule. The symbol ```N, s, t, r``` represents the size of vertices, the vertex that sends charge, the vertex that receives charge, the amount of charge sent respectively. ```v_1, v_2,..., v_N``` denote the vertices of a rule

```
dummy-line
N s t r
1 (The degree of v_1) (The list of neighbors of v_1)
2 (The degree of v_2) (The list of neighbors of v_2)
...
N (The degree of v_N) (The list of neighbors of v_N)
```
An example (rule1).
```
rule1
2 1 2 2
1 5  2
2 5+ 1
```