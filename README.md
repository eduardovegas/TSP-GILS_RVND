# TSP-GILS_RVND
This project implements the [GILS-RVND](https://www.sciencedirect.com/science/article/abs/pii/S037722171200269X) metaheuristic for solving the Traveling Salesman Problem (TSP).

Compiling and Running
----------------------
A serial executable is obtained by using your favourite C++ compiler to compile and link the downloaded files.

On Windows:

**g++ -O3 src/main.cpp src/readData.cpp -o tsp.exe**

After done compiling, execute the 'tsp' program passing the instance name as an argument, like this:

**./tsp instances/name.tsp**
