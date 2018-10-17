#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ClusterFunctions.h"

int* selectKBits(int N, int k)
{
    int i, j, sameExists;
    int *indexes = malloc(k*sizeof(int));
 
    indexes[0] = rand() % N;            /* Initializes first bit outside of loop so it can compare two of the for equality inside loop */
    for(i = 1; i < k;)
    {
        indexes[i] = rand() % N;
        
        sameExists = 0;
        for(j = 0; j < i; j++)         /* Checking if array indexes selected up to i have been selected before and if so randomly selects again */
        {
            if(indexes[i] == indexes[j])
                sameExists = 1;
        }
        if(sameExists == 0)           /* if same array index doesn't exist in already randomly selected k-bts, continue with next random bit */
            i++;
    }
    return indexes;                     /* returns the array of k-randomly selected bits (in array indexes, not actual hamming bits) to hash point in bucket */
}

int sameMedoids(int *oldMedoids, int *newMedoids, int k)
{
    int i, j, notExists;
    
    notExists = 1;
    /* For each old medoid */
    for(i = 0; i < k; i++)
    {
        notExists = 1;
        /* Check if it has been selected before in any position of newMedoids */
        for(j = 0; j < k; j++)
        {
            if(oldMedoids[i] == newMedoids[j])
            {
                /* If found, continue checking for other old medoids */
                notExists = 0;
                break;
            }
        }
        /* If not found, arrays don't have the same elements, return 0 */
        if(notExists == 1)
            return 0;
    }
    return 1;
}

int* initializationConcentrate(float **distanceMatrix, int itemCount, int k)
{
    int i, j;
    int *initMedoids, *vIndexes;
    double *v, *distSums;
    
    v = malloc(sizeof(double)*itemCount);
    vIndexes = malloc(sizeof(int)*itemCount);
    distSums = malloc(sizeof(double)*itemCount);
    initMedoids = malloc(sizeof(int)*k);
    if(v == NULL || vIndexes == NULL || distSums == NULL || initMedoids == NULL)
    {	perror("Malloc failed"); exit(-1); }
    
    for(i = 0; i < itemCount; i++)
    {
        /* Calculating denominator for vi values (for each item) */
        distSums[i] = 0;
        for (j = 0; j < itemCount; j++)
       		distSums[i] += distanceMatrix[i][j];
    }
    for(i = 0; i < itemCount; i++)
    {
        v[i] = 0;
        vIndexes[i] = i;    /* Keep item index for sorting */
        /* Calculating vi values for each item */
        for (j = 0; j < itemCount; j++)
            v[i] += (double)(distanceMatrix[i][j]/distSums[j]);
    }
    /* Sort both arrays to insure valid vIndexes values */
    SelectionSortTwoArrays(v, vIndexes, itemCount);
    /* Select k items with smallest vi values */
    for (i = 0; i < k; i++)
        initMedoids[i] = vIndexes[i];
    
    free(distSums); free(v); free(vIndexes);
    
    return initMedoids;
}

void SelectionSortTwoArrays(double *array, int *array2, int arrayLength)
{
	double temp;
    int i, j;
    
    for(i = 0; i < arrayLength - 1; i++)
    {
        for(j = i+1; j < arrayLength; j++)
        {
            if (array[i] > array[j])
            {
                temp = array[i];
                array[i] = array[j];
                array[j] = temp;
                
                temp = array2[i];
                array2[i] = array2[j];
                array2[j] = temp;
            }
        }
    }
}

int* assignmentPAM(float **distanceMatrix, int *medoids, int itemCount, int k, double *totalCost)
{
    int i, j;
    int *clustersArray;
    double minDistance;
    
    clustersArray = malloc(sizeof(int)*itemCount);
    
    *totalCost = 0;
    for(i = 0; i < itemCount; i++) /* For each item */
    {
        /* Choose a large number as min distance */
        minDistance = (double)RAND_MAX;
        for(j = 0; j < k; j++) /* For each medoid */
        {
            /* Find medoid with min distance from item */
            if(minDistance > distanceMatrix[i][medoids[j]])
            {
                minDistance = distanceMatrix[i][medoids[j]];
                clustersArray[i] = j;   /* Item will go to medoids j cluster */
            }
        }
        *totalCost += minDistance; /* Calculate total cost */
    }
    return clustersArray;
}

int* updateClarans(float **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k, double *totalCost, int numLocal, int maxNeighbor)
{
    int i, j, m, *currentBestMedoids, *nonMedoids, temp, *assignm, whichChanges, whichGets, isMedoid;
    double optimalCost, localCost, currentCost;
    
    nonMedoids = malloc(sizeof(int)*(itemCount - k));
    currentBestMedoids = malloc(sizeof(int)*k);
    if(nonMedoids == NULL || currentBestMedoids == NULL )
    	{ perror("ClusterFunctions.c: updateClarans: malloc failed"); exit(-1); }
    
    optimalCost = *totalCost;
    
    for(i = 0, j = 0; i < itemCount; i++) /* Creating an array that contains all non medoid items */
    {
        isMedoid = 0;
        for(m = 0; m < k; m++)
        {
            if(medoids[m] == i)
            { isMedoid = 1; break; }
        }
        
        if(isMedoid==0)
        {
            nonMedoids[j]=i;
            j++;
        }
    }
    
    for(i = 0; i < numLocal; i++)
    {
        currentBestMedoids =  selectKBits(itemCount, k); /* Randomly select k new medoids */
        assignm = assignmentPAM(distanceMatrix, currentBestMedoids, itemCount, k, &localCost); /* Make the assignment to update (local) cost */

        for(j = 0; j < maxNeighbor; j++)
        {
            whichChanges = rand()%k; /* Randomly select a medoid to change */
            whichGets = rand()%(itemCount - k); /* Randomly select a non medoid item */
            /* Interchange medoid with non medoid */
            temp = currentBestMedoids[whichChanges];
            currentBestMedoids[whichChanges] = nonMedoids[whichGets];
            nonMedoids[whichGets] = temp;
            
            assignm = assignmentPAM(distanceMatrix, currentBestMedoids, itemCount, k, &currentCost); /* Make the assignment to calculate new cost */
            
            if(currentCost < localCost)
            {
                j = -1; /* If newis lower than previous cost, restart */
                localCost = currentCost;
            }
            else /* Otherwise, revert back */
            {
                temp = currentBestMedoids[whichChanges]; 
                currentBestMedoids[whichChanges] = nonMedoids[whichGets];
                nonMedoids[whichGets] = temp;
            }
        }
        
        if(optimalCost > localCost) /* If optimal cost is greater than local cost */
        {
            optimalCost = localCost; /* Update optimal cost */
            for(j = 0; j < k; j++)     /* Update medoids to new configuration */
                medoids[j] = currentBestMedoids[j];

            assignm = assignmentPAM(distanceMatrix, currentBestMedoids, itemCount, k, &localCost);
        }
    }
    clustersArray = assignmentPAM(distanceMatrix, medoids, itemCount, k, totalCost);
    
    return clustersArray;
}

void medoidDistance(float **distanceMatrix, double *minDistance, double *maxDistance, int *medoids, int k)
{
    int i, j;
    double tempDistance;
    
    *minDistance = (double)RAND_MAX;    /* Set min distance to a large number */
    *maxDistance = 0.0;
    
    for(i = 0; i < k; i++)
    {
        for(j = 0; j < k; j++)
        {
            /* Do not calculate distance from itself */
            if(i == j)
                continue;
            
            tempDistance = distanceMatrix[medoids[i]][medoids[j]];
            if(tempDistance > *maxDistance)
                *maxDistance = tempDistance;
            else if(tempDistance < *minDistance)
                *minDistance = tempDistance;
            else
                continue;
        }
    }
}

double* silhouette(float **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k)
{
    int i, j, m, searchCluster, counter;
    double *a, *b, *s, sum, medoidDistance;
    
    s = malloc(sizeof(double)*(k+1));
    b = malloc(sizeof(double)*itemCount);
    a = malloc(sizeof(double)*itemCount);
    if(a == NULL || b == NULL || s == NULL)
    { printf("Malloc failed in ClusterFunctions.c...Programme will terminate"); exit(-1); }
    
    for(i = 0; i < (k+1); i++)
        s[i] = 0.0;
    for(i = 0; i < itemCount; i++) /* For each item */
    {
        counter = 0;
        sum = 0.0;
        
        for(j = 0; j <itemCount; j++) /* For each item, calculate sum of distances from other items in its cluster */
        {
            if(clustersArray[i] == clustersArray[j])
            {
                counter++;
                sum += distanceMatrix[i][j];
            }
        }
        if(counter==0)
            a[i]=0;
        else
            a[i] = sum/counter; /* Calculate mean distance from items in same cluster */
        
        counter = 0;
        sum = 0.0;
        medoidDistance = (double)RAND_MAX; /* Set min distance to neighbour medoids to large number */
        
        for(m = 0; m < k; m++) /* For each medoid, find medoid with min distance from item */
        {
            if(clustersArray[i] == m)  /* If medoid from same cluster, continue */
                continue;
            
            if((medoidDistance > distanceMatrix[i][medoids[m]]) && (distanceMatrix[i][medoids[m]] != 0))
            {
                medoidDistance = distanceMatrix[i][medoids[m]];
                searchCluster = m;
            }
        }
        
        for(j = 0; j < itemCount; j++) /* Calculate distance from items of searchCluster (cluster with min distance from item) */
        {
            if(clustersArray[j] == searchCluster)
            {
                counter++;
                sum += distanceMatrix[i][j];
            }
        }
        b[i] = sum/counter; /* Calculate mean distance from items of neighbor cluster */
        
        if(b[i] > a[i])
            sum = (b[i] - a[i])/b[i];
        else
            sum = (b[i] - a[i])/a[i];
        s[i] = sum/counter;
        s[k] += sum;
    }
    s[k] /= itemCount;
    
    free(b); free(a);
    
    return s;
}
