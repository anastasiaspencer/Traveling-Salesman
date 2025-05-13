# Traveling-Salesman
This includes a Brute Force, Nearest Neighbor and Original solution coined Trivecta along with a detailed report on the Trivecta algorithm.

TRIVECTA

TRIVECTA is a project demonstrating the TRIVECTA algorithm, which provides an optimal solution to the Traveling Salesman Problem.

Installation

Download TriVecta.cpp.
Generate an input graph in the form of a diagonal matrix. The matrix should represent the distances between nodes in the graph. Here's an example format:

0
3032 0
25282 10805 0
99623 67938 24596 0
209741 162683 90162 20927 0

Update the following lines in TriVecta.cpp:
Change the graph size in line 19 to match the number of nodes in the generated graph.
Update the node size in line 18 accordingly.
Change the name of the input file in line 40 to match the generated graph file.

Usage

Compile the program using the following command:

g++ TriVecta.cpp 

Run the program:

./a.out

Example Output

After running the program with a graph of size 30, you'll get an output similar to the following:

Tour: 2 1 0 10 20 21 11 12 22 13 23 24 25 26 27 17 28 18 29 19 9 8 7 16 6 5 15 4 14 3 
Length of the tour: 346645
Time taken by program: 0.041424 sec

This output provides the optimal tour, the length of the tour, and the time taken by the program to execute.