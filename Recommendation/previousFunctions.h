#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MAX_HAMMING_SIZE 70
#define MAX_DIST 99999999
#define W 4
#define MAXDISTANCE 999999999

double hammingDistance(double *item1, double *item2, int dimension);
double eucledianDistance(double *item1, double *item2, double avg1, double avg2, int dimension);
double cosineDistance(double *item1, double *item2, double avg1, double avg2, int dimension);
int* concentrate(double** array,double ** all_distances, char metriki, char type, int dimensions, int plithos_items, double* avgRatingsPerUser, int K);
int* PAM_assignment(double** array, char distance_metric, int dimensions, int plithos_items, int *medoids, double* avgRatingsPerUser, int K);
double* Silhouette(double** array, double** all_distances, int plithos_items, int dimensions, char metriki, int* medoids, int* ass, double* avgRatingsPerUser, int K);
double findDistance(double** array, int item1, int item2, int dimensions, double* avgRatingsPerUser, char metriki);
int* loydsNEW(double** array, int* assignments, int plithos_items, int dimensions, char distance_metric, int K, int* best_meds, double* avgRatingsPerUser);
int* k_medoids(double** array,double ** all_distances, char distance_metric, int dimensions, int plithos_items, int K,double * avgRatingsPerUser);

