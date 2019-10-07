#include <iostream>  // cout
#include <fstream>   // ifstream
#include <string.h>  // strncpy
#include <stdlib.h>  // rand
#include <math.h>    // sqrt, pow
#include <omp.h>     // OpenMP
#include "Timer.h"
#include "Trip.h"
#include <bits/stdc++.h>
#include <algorithm>
#include <stdlib.h>     /* srand, rand */
#include <time.h>

using namespace std;

// Already implemented. see the actual implementations below
void initialize(Trip trip[CHROMOSOMES], int coordinates[CITIES][2]);
void select(Trip trip[CHROMOSOMES], Trip parents[TOP_X]);
void populate(Trip trip[CHROMOSOMES], Trip offsprings[TOP_X]);
bool compareTwoTrips(Trip firstTrip, Trip secondTrip);
int getCharPosition(const char *array, size_t size, char charToFind);
void insertCity(char charToBeInsertedAttempt1, char charToBeInsertedAttempt2, int indexInOffspring,
                char offspringIternary[], Trip parent);
void compliment(char offspringIternary[], char offspringComplimentItenary[]);

// need to implement for your program 1
extern void evaluate(Trip trip[CHROMOSOMES], int coordinates[CITIES][2]);
extern void crossover(Trip parents[TOP_X], Trip offsprings[TOP_X], int coordinates[CITIES][2]);
extern void mutate(Trip offsprings[TOP_X]);

/*
 * MAIN: usage: Tsp #threads
 */
int main(int argc, char *argv[]) {
  Trip trip[CHROMOSOMES];       // all 50000 different trips (or chromosomes)
  Trip shortest;                // the shortest path so far
  int coordinates[CITIES][2];   // (x, y) coordinates of all 36 cities:
  int nThreads = 1;

  // verify the arguments
  if (argc == 2)
    nThreads = atoi(argv[1]);
  else {
    cout << "usage: Tsp #threads" << endl;
    if (argc != 1)
      return -1;  // wrong arguments
  }
  cout << "# threads = " << nThreads << endl;

  // shortest path not yet initialized
  shortest.itinerary[CITIES] = 0;  // null path
  shortest.fitness = -1.0;         // invalid distance

  // initialize 50000 trips and 36 cities' coordinates
  initialize(trip, coordinates);

  // start a timer 
  Timer timer;
  timer.start();

  // change # of threads
  omp_set_num_threads(nThreads);

  // define TOP_X parents and offsprings.
  Trip parents[TOP_X], offsprings[TOP_X];

  // find the shortest path in each generation
  for (int generation = 0; generation < 150; generation++) {

    // evaluate the distance of all 50000 trips
    evaluate(trip, coordinates);

    // just print out the progress
    if (generation % 20 == 0)
      cout << "generation: " << generation << endl;

    // whenever a shorter path was found, update the shortest path
    if (shortest.fitness < 0 || shortest.fitness > trip[0].fitness) {

      strncpy(shortest.itinerary, trip[0].itinerary, CITIES);
      shortest.fitness = trip[0].fitness;

      cout << "generation: " << generation << " shortest distance = " << shortest.fitness << "\t itinerary = "
           << shortest.itinerary << endl;
    }

    // cleaning up the offsprings in every generation
    memset(offsprings, '\0', sizeof(Trip) * TOP_X);

    // choose TOP_X parents from trip
    select(trip, parents);

    // generates TOP_X offsprings from TOP_X parents
    crossover(parents, offsprings, coordinates);

    // mutate offsprings
    mutate(offsprings);

    // populate the next generation.
    populate(trip, offsprings);
  }

  // stop a timer
  cout << "elapsed time = " << timer.lap() << endl;
  return 0;
}

/*
 * Initializes trip[CHROMOSOMES] with chromosome.txt and coordiantes[CITIES][2] with cities.txt
 *
 * @param trip[CHROMOSOMES]:      50000 different trips
 * @param coordinates[CITIES][2]: (x, y) coordinates of 36 different cities: ABCDEFGHIJKLMNOPQRSTUVWXYZ
 */
void initialize(Trip trip[CHROMOSOMES], int coordinates[CITIES][2]) {
  // open two files to read chromosomes (i.e., trips)  and cities
  ifstream chromosome_file("chromosome.txt");
  ifstream cities_file("cities.txt");

  // read data from the files
  // chromosome.txt:                                                                                           
  //   T8JHFKM7BO5XWYSQ29IP04DL6NU3ERVA1CZG                                                                    
  //   FWLXU2DRSAQEVYOBCPNI608194ZHJM73GK5T                                                                    
  //   HU93YL0MWAQFIZGNJCRV12TO75BPE84S6KXD
  for (int i = 0; i < CHROMOSOMES; i++) {
    chromosome_file >> trip[i].itinerary;
    trip[i].fitness = 0.0;
  }

  // cities.txt:                                                                                               
  // name    x       y                                                                                         
  // A       83      99                                                                                        
  // B       77      35                                                                                        
  // C       14      64                                                                                        
  for (int i = 0; i < CITIES; i++) {
    char city;
    cities_file >> city;
    int index = (city >= 'A') ? city - 'A' : city - '0' + 26;
    cities_file >> coordinates[index][0] >> coordinates[index][1];
  }

  // close the files.
  chromosome_file.close();
  cities_file.close();

  // just for debugging
  if ( DEBUG) {
    for (int i = 0; i < CHROMOSOMES; i++)
      cout << trip[i].itinerary << endl;
    for (int i = 0; i < CITIES; i++)
      cout << coordinates[i][0] << "\t" << coordinates[i][1] << endl;
  }
}

/**
 *  Comparator function for sorting the itinerary
 */
bool compareTwoTrips(Trip firstTrip, Trip secondTrip) {
  return (firstTrip.fitness < secondTrip.fitness);
}

/**
 *  Sort the trips in ascending order based on their fitness ( distance of the trip ).
 *  @param trip[CHROMOSOMES]:      50000 different trips
 *  @param coordinates[CITIES][2]:  (x, y) coordinates of 36 different cities: ABCDEFGHIJKLMNOPQRSTUVWXYZ
 */
void evaluate(Trip trip[CHROMOSOMES], int coordinates[CITIES][2]) {

#pragma omp parallel for

  for (int j = 0; j < CHROMOSOMES; j++) {

    // Initializes every trip fitness to zero.
    trip[j].fitness = 0;
    for (int i = 0; i < CITIES - 1; i++) {

      char city1 = trip[j].itinerary[i];
      char city2 = trip[j].itinerary[i + 1];

      int index_city_1 = (city1 >= 'A') ? city1 - 'A' : city1 - '0' + 26;
      int city_x1 = coordinates[index_city_1][0];
      int city_y1 = coordinates[index_city_1][1];

      int index_city_2 = (city2 >= 'A') ? city2 - 'A' : city2 - '0' + 26;
      int city_x2 = coordinates[index_city_2][0];
      int city_y2 = coordinates[index_city_2][1];

      float distbtwcities = sqrt(pow((city_x2 - city_x1), 2) + pow((city_y2 - city_y1), 2));
      trip[j].fitness += distbtwcities;

    }
  }

  sort(trip, trip + CHROMOSOMES, compareTwoTrips);

}

/*
 * Select the first TOP_X parents from trip[CHROMOSOMES]
 *
 * @param trip[CHROMOSOMES]: all trips
 * @param parents[TOP_X]:    the firt TOP_X parents
 */
void select(Trip trip[CHROMOSOMES], Trip parents[TOP_X]) {
  // just copy TOP_X trips to parents
  for (int i = 0; i < TOP_X; i++)
    strncpy(parents[i].itinerary, trip[i].itinerary, CITIES + 1);
}

/**
 *  Generates TOP_X offsprings from TOP_X parents.
 *  @param parents[TOP_X]:    the firt TOP_X parents
 *  @param coordinates[CITIES][2]:  (x, y) coordinates of 36 different cities: ABCDEFGHIJKLMNOPQRSTUVWXYZ
 */
void crossover(Trip parents[], Trip offsprings[], int coordinates[CITIES][2]) {

#pragma omp parallel for
  for (int i = 0; i < TOP_X; i += 2) {

    // Starting from (0,0) for every iteration
    int city_x0 = 0, city_y0 = 0;
    int j = 0;

    /**
     *  Choosing the first city in child[i] between parent[i] and parent[i+1]
     */

    char city_parent_1 = parents[i].itinerary[j];
    int index_city_1_parent_1 = (city_parent_1 >= 'A') ? city_parent_1 - 'A' : city_parent_1 - '0' + 26;
    int city_x_parent_1 = coordinates[index_city_1_parent_1][0];
    int city_y_parent_1 = coordinates[index_city_1_parent_1][1];

    float distbtwcities_in_parent_1 = sqrt(pow((city_x_parent_1 - city_x0), 2) + pow((city_y_parent_1 - city_y0), 2));

    char city_parent_2 = parents[i + 1].itinerary[j];

    int index_city_1_parent_2 = (city_parent_2 >= 'A') ? city_parent_2 - 'A' : city_parent_2 - '0' + 26;

    int city_x_parent_2 = coordinates[index_city_1_parent_2][0];
    int city_y_parent_2 = coordinates[index_city_1_parent_2][1];

    float distbtwcities_in_parent_2 = sqrt(pow((city_x_parent_2 - city_x0), 2) + pow((city_y_parent_2 - city_y0), 2));

    if (distbtwcities_in_parent_1 <= distbtwcities_in_parent_2) {
      insertCity(city_parent_1, city_parent_2, j, offsprings[i].itinerary, parents[0]);
    } else {
      insertCity(city_parent_2, city_parent_1, j, offsprings[i].itinerary, parents[0]);
    }

    for (j = 1; j < CITIES; j++) {

      float distbtwcities_child_parent1, distbtwcities_child_parent2;
      char city_child = offsprings[i].itinerary[j - 1];
      int index_city_1 = (city_child >= 'A') ? city_child - 'A' : city_child - '0' + 26;
      int city_x_child_1 = coordinates[index_city_1][0];
      int city_y_child_1 = coordinates[index_city_1][1];

      // index of childs city in parent[i]
      int city_parent_1_pos = getCharPosition(parents[i].itinerary, CITIES, city_child);

      // In case selected city is the last city in the chromosome, then
      // we will choose the first city in it.
      if (city_parent_1_pos + 1 >= CITIES) {
        city_parent_1 = parents[i].itinerary[0];
      } else {
        // index of child city's next city in the parent[i]
        city_parent_1 = parents[i].itinerary[city_parent_1_pos + 1];
      }

      int index_city_parent_1 = (city_parent_1 >= 'A') ? city_parent_1 - 'A' : city_parent_1 - '0' + 26;
      int city_x2 = coordinates[index_city_parent_1][0];
      int city_y2 = coordinates[index_city_parent_1][1];
      distbtwcities_child_parent1 = sqrt(pow((city_x2 - city_x_child_1), 2) + pow((city_y2 - city_y_child_1), 2));

      // index of childs city in parent[i+1]
      int city_parent_2_pos = getCharPosition(parents[i + 1].itinerary, CITIES, city_child);

      // In case selected city is the last city in the chromosome, then
      // we will choose the first city in it.
      if (city_parent_2_pos + 1 >= CITIES) {
        city_parent_2 = parents[i + 1].itinerary[0];
      } else {
        // index of child city's next city in the parent[i+1]
        city_parent_2 = parents[i + 1].itinerary[city_parent_2_pos + 1];
      }

      int index_city_parent_2 = (city_parent_2 >= 'A') ? city_parent_2 - 'A' : city_parent_2 - '0' + 26;
      int city_x4 = coordinates[index_city_parent_2][0];
      int city_y4 = coordinates[index_city_parent_2][1];
      distbtwcities_child_parent2 = sqrt(pow((city_x4 - city_x_child_1), 2) + pow((city_y4 - city_y_child_1), 2));

      // We will select the city which makes the child itinerary distance smaller.
      if (distbtwcities_child_parent1 <= distbtwcities_child_parent2) {
        insertCity(city_parent_1, city_parent_2, j, offsprings[i].itinerary, parents[0]);
      } else {
        insertCity(city_parent_2, city_parent_1, j, offsprings[i].itinerary, parents[0]);
      }

    }

    // Taking the compliment of child[i] to get child[i+1]
    compliment(offsprings[i].itinerary, offsprings[i + 1].itinerary);
    offsprings[i].itinerary[CITIES] = '\0';
    offsprings[i + 1].itinerary[CITIES] = '\0';

  }
}


/**
 * Utility function to insert a city in the offspring.
 * @param charToBeInsertedAttempt1: First preference city to be inserted in child.
 * @param charToBeInsertedAttempt2: Second preference city to be inserted in child.
 * @param indexInOffspring: Position in the offspring to insert the choosen city.
 * @param offspringIternary: Itinerary of the offspring.
 * @param Parent: Shortest distance parent in the generation. Used to choose the random city.
 */
void insertCity(char charToBeInsertedAttempt1, char charToBeInsertedAttempt2, int indexInOffspring,
                char offspringItinerary[], Trip parent) {

  if (getCharPosition(offspringItinerary, CITIES, charToBeInsertedAttempt1) < 0) {
    offspringItinerary[indexInOffspring] = charToBeInsertedAttempt1;
    return;
  } else if (getCharPosition(offspringItinerary, CITIES, charToBeInsertedAttempt2) < 0) {
    offspringItinerary[indexInOffspring] = charToBeInsertedAttempt2;
    return;
  } else {

    for (int i = 0; i < CITIES; i++) {

      if ((getCharPosition(offspringItinerary, CITIES, parent.itinerary[i])) < 0) {
        offspringItinerary[indexInOffspring] = parent.itinerary[i];
        return;
      }
    }
  }
}

/**
 * Utilty function to get the position of the city in the offspring.
 */
int getCharPosition(const char array[], const size_t size, const char charToFind) {

  for (int i = 0; i < size; i++) {
    if (array[i] == charToFind) {
      return i;
    }
  }
  return -1;
}

/**
 *  Utility function to get the compliment of a trip.
 *  @param offspringItinerary: Itineraty of child[i]
 *  @param offspringComplimentItinerary: Itinerary of child[i+1]
 */
void compliment(char offspringItinerary[], char offspringComplimentItinerary[]) {

  for (int i = 0; i < CITIES; i++) {
    int index = (offspringItinerary[i] >= 'A') ? offspringItinerary[i] - 'A' : offspringItinerary[i] - '0' + 26;
    offspringComplimentItinerary[i] = city_names[CITIES - 1 - (index)];
  }
}

/**
 *  Mutates the offsprings based on the mutation rate
 *
 *  @param offsprings[TOP_X]:  TOP_X offsprings.
 */
void mutate(Trip offsprings[TOP_X]) {

  srand(time(NULL));

  int random_number = rand() % 100;

  if (random_number < MUTATE_RATE) {

    for (int i = 0; i < TOP_X; i++) {

      int random_number_1 = rand() % 36;
      int random_number_2 = rand() % 36;

      char temp = offsprings[i].itinerary[random_number_1];
      offsprings[i].itinerary[random_number_1] = offsprings[i].itinerary[random_number_2];
      offsprings[i].itinerary[random_number_2] = temp;
    }
  }
}

/*
 * Replace the bottom TOP_X trips with the TOP_X offsprings
 */
void populate(Trip trip[CHROMOSOMES], Trip offsprings[TOP_X]) {
  // just copy TOP_X offsprings to the bottom TOP_X trips.

#pragma omp parallel for

  for (int i = 0; i < TOP_X; i++)
    strncpy(trip[ CHROMOSOMES - TOP_X + i].itinerary, offsprings[i].itinerary,
    CITIES + 1);

  // for debugging
  if (0) {
    for (int chrom = 0; chrom < CHROMOSOMES; chrom++)
      cout << "chrom[" << chrom << "] = " << trip[chrom].itinerary << ", trip distance = " << trip[chrom].fitness
           << endl;
  }
}
