/*
 * Anastasia Spencer
 * TriVecta algorithm for the Traveling Salesman problem
 * This algorithm performs the nearest neighbor on each different city as a starting point and evaluates the top candidates. It then performs simulated annealing on the top candidates and decides the best path. Last, it performs two-opt swap to further locally refine the final path.
 * Reference: “Simualted Annealing – Solving the Travelling Salesman Problem (TSP),” https://www.codeproject.com/Articles/26758/Simulated-Annealing-Solving-the-Travelling-Salesma, [Online]. Available: URL [Accessed: 20-02-2024].
 * Reference: “C++ Implementation of 2-opt to the “Att48” Travelling Salesman Problem,” https://www.technical-recipes.com/2012/applying-c-implementations-of-2-opt-to-travelling-salesman-problems/. [Online]. Available: URL [Accessed: 20-02-2024].
*/
#include <iostream>
#include <vector>
#include <fstream>
#include <utility>
#include <algorithm>
#include <random>

using namespace std;

// Change as needed for size of input graph
int N = 100; 
int graph[100][100];
int candidateSize;

int initialTemp;
double coolingRate = 0.99999; // Remains as close to 1 as possible for a slow cooldown to allow for more exploration
int maxIterations;

/*
 * This function sets the candidate size, initial temp and the number of iterations based on the size of the input graph
*/
void setVariables(){
    candidateSize = N/10;
    initialTemp = N * 100;
    maxIterations = initialTemp;
    return;
}

/*
 * Function to read in distance matrix
*/
void readDistanceMatrix(){
    ifstream file("Size100.graph"); //change file to read in as needed
    if(!file){
        cerr << "File cannot be opened" << endl;
        return;
    }

    for(int i = 0; i < N; i++){
        for(int j = 0; j<= i; j++){
            file >> graph[i][j];
            if(i != j) graph[j][i] = graph[i][j];
        }
    }
    file.close();
}

/*
 * Function to iterate through all unvisited to cities to determine the closest neighbor and returns the nearest neighbor
 * @param hasVisited represents the list of cities that have been visited thus far
 * @param currentCity is the current city we are looking for a neighbor for
*/
int findNeighbor(vector<bool> hasVisited, int currentCity){
    int currentMin = numeric_limits<int>::max();
    int nearestNeighbor = -1;
    //If the city has not been visited and it is not the current city
    for(int i = 0; i < N; i++){
        if(!hasVisited[i] && i != currentCity){
            int dist = graph[currentCity][i];
            //Check if this distance is smaller than the current minimum distance
            if(dist < currentMin){
                currentMin = dist;
                nearestNeighbor = i;
            }
        }
    }

    return nearestNeighbor;
}

/*
 * This function creates the path as the nearest neighbor of each city is determined. It iteratively calls the findNeighbor() function and marks that the new neighbor has been visited
 * This function is run on each city as a new starting point and stores each path in a vector that is returned at the termination of the function
*/
vector< pair<int, vector<int> > > nearestNeighbor(){
    int numCities = N;
    int minDistance = numeric_limits<int>::max();
    vector< pair<int, vector<int> > > paths;

    for(int startCity = 0; startCity < numCities; startCity++){

    pair<int, vector<int> > currentPath;
    int distance = 0;
    vector<bool> hasVisited(N, false);
    vector<int> path;

    int currentCity = startCity;
    path.push_back(currentCity);
    hasVisited[currentCity] = true;

    while(path.size() < N){
        int nearestNeighbor = findNeighbor(hasVisited, currentCity);
        // As long as the nearest neighbor was found (we have not visited each city yet) add to the path and calculate the new distance
        if(nearestNeighbor != -1){
            path.push_back(nearestNeighbor);
            hasVisited[nearestNeighbor] = true;
            distance += graph[currentCity][nearestNeighbor];
            currentCity = nearestNeighbor;
        }else{
            break;
        }
         //path.first
    }
    //adds the path from the last city of the path back to the starting city.
    distance += graph[startCity][currentCity];
    //add to vector storing all the paths 
     currentPath.first = distance;
     currentPath.second = path;
    paths.push_back(currentPath);
    }

    return paths;
}

/*
 * Comparison function as a basis for the parameter needed in the nth_element function (executed in line 134)
 * @param a the beginning of the path
 * @param b the path plus the candidate size
*/
bool comparePairs(const pair<int, vector<int> >& a, const pair<int, vector<int> >& b){
    return a.first < b.first;
}

/*
 * Takes in all of the paths and partially sorts using the nth_element library function
 * Chooses the top N candidates and returns a new vector that only contains the candidate paths
 * @param paths this is a vector that contains pairs of paths and lengths to be partially sorted
*/
vector< pair<int, vector<int> > > selectCandidates( vector < pair< int, vector<int> > > paths ){
    //sort elements based on critera defined in lambda function
    nth_element(paths.begin(), paths.begin() + candidateSize, paths.end(), comparePairs); 


    //Copy the top candidates to a new vector
    vector< pair<int, vector<int> > > candidates(paths.begin(), paths.begin() + candidateSize);

    return candidates;

}
/*
 * This function calculates the tour length
 * @param newTour is the tour to be calculated
*/
int calculateTourLength(vector<int> newTour){
    int length = 0;
    for(int i = 0; i < newTour.size() - 1; i++){
        int from = newTour[i];
        int to = newTour[i + 1];
        length += graph[from][to];
    }

    length += graph[newTour.back()][newTour.front()];
    return length;

}

/*
 * This function reverses a given tour between given indicies
 * @param tour is the path to be swapped
 * i and k are the given indicies
*/
pair< int, vector<int > > swap(vector<int> tour, int i, int k){
    vector<int> newTour = tour;

    //Reverse the subtour between indicies i and k
    reverse(newTour.begin() + i, newTour.begin() + k + 1);

    int length = calculateTourLength(newTour);
    //return the swapped tour
    return make_pair(length, newTour);
}

/*
 * This function facilitates the two opt swap phase of the algorithm
 * It systematically iterates through each local combination of indidcies and evaluates if the new path is an improvement
 * Returns the best performing tour
 * @param tour this is the path to optimize
*/
pair< int, vector<int > > twoOpt(pair<int, vector<int > > tour){
    vector<int> bestTour = tour.second;
    int bestLength = tour.first;

    for(int i = 0; i< N - 1; i++){
        for(int k = i + 1; k < N; k++){
            //This performs 2-opt swap between edges (i, i+1) and (k, k+1)
            pair< int, vector<int> > currentTour = swap(tour.second, i, k);
            //evaluate if the new path is an improvement
            if(currentTour.first < bestLength){
                bestLength = currentTour.first;
                bestTour = currentTour.second;
            }
        }
    }
    //return the best generated tour
    pair< int, vector<int > > topTour;
    topTour.first = bestLength;
    topTour.second = bestTour;
    return topTour;
}

  /*
  * This function performs simulated annealing on the passed in path, it randomly chooses two indicies to swap and will occasionally accept a worse solution based on the metropolis criteron
  * @param initialTour the path to be annealed
  */
  pair<int, vector<int> > simulatedAnnealing(pair<int, vector<int> > initialTour){
        vector<int> currentTour = initialTour.second;
        int currentLength = initialTour.first;
        vector<int> bestTour = currentTour;
        int bestLength = currentLength;

        //for generating random values within a uniform probability density function to aid in achieving the metropolis criterion
        random_device rd;
        default_random_engine rng(rd());
        uniform_real_distribution<double> dist(0.0, 1.0);

        for(int iteration = 0; iteration < maxIterations; iteration++){

            // randomly generate two indicies to make random swap
            int i = rand() % (N - 1);
            int k = i + 1 + rand() % (N - i - 1);

            //swap our randomly generated edges (reverse the tour between them)
            pair< int, vector<int > > swappedTour = swap(currentTour, i, k);

            //calculate the change in cost
            int delta = swappedTour.first - currentLength;

            
            //If the new tour is an improvement from the previous tour, accept or if the new tour meets the metropolis criteron (this is so worse solutions have a probability of being accepted)
            //worse tours will occasionally get accepted to allow for exploration- the higher the temp the more likely a worse solution is to be accepted
            if(delta < 0 || dist(rng) < exp(-delta / initialTemp)){
                currentTour = swappedTour.second;
                currentLength = swappedTour.first;

                if(currentLength < bestLength){
                    bestTour = currentTour;
                    bestLength = currentLength;
                }
            }

            //cool down the temp by a factor of the cooling rate
            initialTemp *= coolingRate;
        }

        return make_pair(bestLength, bestTour);

    }

  /*
  * This function runs simulated annealing on the candidate paths
  * Once the best path has been decided after simulated annealing, it undergoes two-opt
  * @param candidates this is the vector of pairs of the candidate paths with their respective lengths
  */
pair<int, vector<int> > annealAndTwoOpt(vector< pair<int, vector<int> > > candidates){
    int bestLength = numeric_limits<int>::max();
    vector<int> bestPath;
    pair<int, vector<int> > bestTwoOptTour;

    // Loop through each candidate tour and perform a 2-opt swap
    for(const auto& tour : candidates){
        pair<int, vector<int> > annealedTour = simulatedAnnealing(tour);
       
            
            if(annealedTour.first < bestLength){
                bestLength = annealedTour.first;
                bestPath = annealedTour.second;
            }
        }
        //send best path to two-opt swap
        bestTwoOptTour = twoOpt(make_pair(bestLength, bestPath));

        return bestTwoOptTour;
    }

  

int main(){
    clock_t start, end;
    start = clock();
    //Read in input matrix
    readDistanceMatrix();
    //Set any variables that are determined by input graph
    setVariables();
    //Phase I run nearest neighbor on each city as a starting point
    vector < pair< int, vector<int> > > paths = nearestNeighbor();
    //Selects the top candidates to undergo simulated annealing
    vector< pair<int, vector<int> > > candidates = selectCandidates(paths);
    //Within this funciton, simualted annealing and then two opt is performned
    pair<int, vector<int> > finalTour = annealAndTwoOpt(candidates);


   //Print the tour
    std::cout << "Tour: ";
    for (int node : finalTour.second) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    std::cout << "Length of the tour: " << finalTour.first << std::endl;

    end = clock();

    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by program: " << fixed
        << time_taken << setprecision(5);
    cout << " sec " << endl;

    
    return 0;
}