#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "previousFunctions.h"

/******************* Distances ******************/
double hammingDistance(double *item1, double *item2, int dimension)
{
	int i, commons = 0;
	double sum = 0;
	
	for(i=0; i<dimension; i++)
	{
		if(item1[i] == 0 || item2[i] == 0) 
			continue;
			
		commons++;
		if(item1[i]>2.5 && item2[i] < 2.5)
			sum++;
		if(item2[i]>2.5 && item1[i] < 2.5)
			sum++;
	}
	
	if(commons == 0)
		return MAX_DIST;
		
	return sum/commons;
}

double eucledianDistance(double *item1, double *item2, double avg1, double avg2, int dimension)
{
	int i, common = 0;
	double sum = 0.0;
	for(i=0; i<dimension; i++)
	{
		if(item1[i] == 0 || item2[i] == 0) 
			continue;
	
		common++;
		sum = sum + (item1[i]-avg1-item2[i]+avg2)*(item1[i]-avg1-item2[i]+avg2);	
	}
	
	if(common==0) return MAX_DIST;
	// printf("sum=%f\n",sum);
// 	printf("sqrt(sum)=%f\n",sqrt(sum));
	if(isnan(sqrt(sum))){
		return MAX_DIST;
	}
	return sqrt(sum);
}

double cosineDistance(double *item1, double *item2, double avg1, double avg2, int dimension)
{
	int i, common = 0;
	double ar = 0.0, par1 = 0.0, par2 = 0.0;
	for(i=0; i<dimension; i++)
	{
		if(item1[i] == 0 || item2[i] == 0) 
			continue;
			
		common++;
		ar = ar + (item1[i]-avg1)*(item2[i]-avg2);
		par1 = par1 + (item1[i]-avg1)*(item1[i]-avg1);
		par2 = par2 + (item2[i]-avg2)*(item2[i]-avg2);
	} 	
	
	if(common == 0){//printf("\nEPISTREFW MAX DIST");
		return MAX_DIST; }
	//printf("\n epistrefw dist %lf\n", (ar / (sqrt(par1)*sqrt(par2)) )	);
	
	if(isnan(ar)){
		printf("\nError: NAN when returning cosine distance 1. Aborting...\n");
		exit(-1);
	}
	
	if(isnan((sqrt(par1)*sqrt(par2)))){
		printf("\nError: NAN when returning cosine distance 2. Aborting...\n");
		exit(-1);
	}
	if(isnan(fabs(ar / (sqrt(par1)*sqrt(par2)) ))){
		//printf("\nError: NAN when returning cosine distance 3.\n Reporting ar=%f and par1=%f and par2=%f and (sqrt(par1)*sqrt(par2))=%f. Aborting...\n",ar,par1,par2,(sqrt(par1)*sqrt(par2)));
		//exit(-1);
		// if(ar<0.00001)
// 			return 0.0;
// 		else
			return MAX_DIST;
	}
	return fabs(ar / (sqrt(par1)*sqrt(par2)) );
}

double findDistance(double** array, int item1, int item2, int dimensions, double* avgRatingsPerUser, char metriki)
{
	double dist;
	if(metriki == 't'){ dist = array[item1][item2];}
	else if(metriki == 'h'){ dist = hammingDistance(array[item1], array[item2], dimensions);}
	else if(metriki == 'e'){ dist = eucledianDistance(array[item1], array[item2], avgRatingsPerUser[item1], avgRatingsPerUser[item2], dimensions);}
	else { dist = cosineDistance(array[item1], array[item2], avgRatingsPerUser[item1], avgRatingsPerUser[item2], dimensions);}
	return dist;
}
/************************* Distances END ***************************/





/************************ Clustering Functions *********************/
int* concentrate(double** array,double ** distances, char metriki, char type, int dimensions, int plithos_items, double* avgRatingsPerUser, int K) 
{
	int i, j, t, *medoids, thesi_min;
	double *vi, *sumDist, sum, sum2, min;
		
	vi = (double*)malloc(plithos_items*sizeof(double));
	sumDist = (double*)malloc(plithos_items*sizeof(double));
	medoids = (int*)malloc(K*sizeof(int));
		
	for(i=0; i<plithos_items; i++)
	{	//printf("\nLOL1  %d", i);
		sumDist[i] = 0.0;
		for(j=0; j<plithos_items; j++){
			//if(distances[i][j]!=MAX_DIST && !isnan(distances[i][j]))
				sumDist[i] = sumDist[i] +  distances[i][j];
		}
	}
	
	
	for(i=0;i<plithos_items;i++)
	{
		vi[i] = 0.0;
		for(j=0;j<plithos_items;j++)
		{
		//printf("cccc %f\n",distances[i][j]);
		//	printf("bbbb %f\n",sumDist[j]);
		//	printf("aaaa  %f\n",distances[i][j]/sumDist[j]);
			//if(distances[i][j]!=MAX_DIST && !isnan(distances[i][j]))
				vi[i] = vi[i] + distances[i][j]/sumDist[j];
			
		//	printf("vi[%d]=%f\n",i,vi[i]);
		}
	}
	
	
	// for(i=0;i<plithos_items;i++)
// 	{
// 		printf("vi[%d]=%f\n",i,vi[i]);
// 	}
// 	int counter=0;
// 	for(j=0; j<plithos_items; j++)
// 		{
// 			if(vi[j]<MAX_DIST && !isnan(vi[j])){
// 				counter++;
// 			}
// 		}
		
	//printf("counter=%d\n",counter);
	//exit(15);	
	for(i=0; i<K; i++)
	{
//		int tyxaio_medoid=rand()%plithos_items;
		min = MAX_DIST;
		thesi_min = -1;
//		min = vi[tyxaio_medoid];
//		thesi_min = tyxaio_medoid;
		for(j=0; j<plithos_items; j++)
		{
			//printf("vi[%d]=%f\n",j,vi[j]);
			if(vi[j] < min)
			{
				min = vi[j];
				thesi_min = j;
			}
		}
		medoids[i] = thesi_min;
// 		if(i==123){
// 			printf("The medoid[123]=user:%d\n",thesi_min);
// 			//exit(-5);
// 		}
		if(thesi_min<0){
				printf("Error: NEGATIVE thesi_min in function concentrate. Aborting...\n");
				exit(-1);
			}
		vi[thesi_min] = MAX_DIST+1;
	}
	return medoids;	
}

int* PAM_assignment(double** array, char distance_metric, int dimensions, int plithos_items, int *medoids, double* avgRatingsPerUser, int K)
{
	int i, j, *assignment;
	double min_dist, dist;
	srand(time(NULL));
	assignment = (int*)malloc(plithos_items*sizeof(int));
	// printf("Dimensions is %d\n",dimensions);
// 	printf("Plithos items is %d\n",plithos_items);
//	exit(1);
	for(i=0; i<plithos_items; i++)
	{
	//	printf("\n to i einai %d", i);
	int	tyxaio_cluster=rand()%K;
		assignment[i] = tyxaio_cluster;	min_dist = MAX_DIST;
		
// 		printf("%c\n",distance_metric);
// 		exit(1);
		for(j=0; j<K; j++)
		{
			if(medoids[j]<0){
				printf("Error: NEGATIVE medoids in function PAM_assignment. Aborting...\n");
				exit(-1);
			}
			if(distance_metric == 'h')
				dist = hammingDistance(array[medoids[j]], array[i], dimensions);
			else if(distance_metric == 'e'){
// 				if(medoids[j]<0)
// 					printf("medoids[j]=%d kai to i einai %d kai to dimensions einai %d\n",medoids[j],i,dimensions);
				dist = eucledianDistance(array[medoids[j]], array[i], avgRatingsPerUser[medoids[j]], avgRatingsPerUser[i],dimensions);
				//printf("Pam dist=%f\n",dist);
				// if(j==123 && (i==546 || i==0)){
// 					printf("dist=%f\n",dist);
// 					//exit(0);
// 				}
			}
			else //if(distance_metric == 'c')
			{//	printf("\nmpika sto teleutaio else");
				//printf("The medoids[%d]=%d\n",j,medoids[j]);
				dist = cosineDistance(array[medoids[j]], array[i], avgRatingsPerUser[medoids[j]], avgRatingsPerUser[i],dimensions);
				
	
			} //printf("\nupologisa to cosine distance = %3lf", dist); }
			// if(i==546){
// 				printf("The i=%d and the dist=%f but the ass=%d\n",i,dist,assignment[i]);
// 			}
			if(dist<min_dist || medoids[j]==i)
			{
				min_dist = dist;
				// if(j==123 || j==0){
// 					printf("dist=%f and the j is %d\n",dist,j);
// 					//exit(0);
// 				}
				assignment[i] = j;
			}
			
		}
	}
	//printf("J is %d\n",K);
	//exit(1);
	for(j=0; j<plithos_items; j++)
	{
		int found=0;
		for(i=0;i<K;i++)
		{
			if(assignment[j]==i){
				found=1;
				break;
			}
		}
		if(!found){
			printf("the assignment of %d was %d\n",j,assignment[j]);
			exit(-2);
		}
		//printf("PAM[%d]=%d\n",j,assignment[j]);
	}
	
	return assignment;
}

double* Silhouette(double** array,double ** all_distances, int plithos_items, int dimensions, char metriki, int* medoids, int* ass, double* avgRatingsPerUser, int K)
{
	int i, j, k, *itemsInCluster, nCluster; //nCluster = neighborCluster
	double min_dist, dist, *a, *b, *s, max; //to max einai gia sugrisi anamesa se a(i) kai b(i)
	a=(double*)malloc(plithos_items*sizeof(double));
	b=(double*)malloc(plithos_items*sizeof(double));
	s=(double*)malloc(plithos_items*sizeof(double));
	itemsInCluster=(int*)malloc(K*sizeof(int));
	
	for(i=0;i<plithos_items;i++){
		b[i]=MAX_DIST;
	}
	
	//upologismos pinaka itemsInCluster
	for(i=0;i<K;i++)
	{
		itemsInCluster[i] = 0;
		for(j=0;j<plithos_items;j++)
		{
			if(ass[j] == i)
				itemsInCluster[i]++;
		}
	}
	
	for(i=0;i<K;i++)
	{
		//printf("the items in cluster: %d\n",itemsInCluster[i]);
	}

	
	//end upologismos pinaka itemsInCluster
				
			
	for(i=0;i<K;i++)
	{
	//printf("to i einai: %d\n",i);
		for(j=0;j<plithos_items;j++) 
		{
			//printf("to j einai: %d\n",j);
			//upologismos a gia to j item
			dist=0.0;
			if(ass[j] == i)
			{
				for(k=0;k<plithos_items;k++)
				{
					if(ass[k] == i /*&& !isnan(all_distances[j][k]) && all_distances[j][k]!=MAX_DIST*/)
						dist+= all_distances[j][k];// findDistance(array, j, k, dimensions, avgRatingsPerUser, metriki);
				}
				a[j]=dist/itemsInCluster[i];
			}//end upologismos a gia to j item
			
			//upologismos nCluster
			min_dist=all_distances[j][medoids[0]]; nCluster=0; dist=0.0;
			for(k=0;k<K;k++)
			{
				if(ass[j] != k )
				{
					//if(!isnan(all_distances[j][medoids[k]]) && all_distances[j][medoids[k]]!=MAX_DIST){
						dist= all_distances[j][medoids[k]];//findDistance(array, j, medoids[k], dimensions, avgRatingsPerUser, metriki);
						if(dist<min_dist)
						{
							min_dist=dist;
							nCluster=k;
						}
					//}
				}
			} //end upologismos nCluster
				
			//upologismos b gia to j item
			dist=0.0;
			if(ass[j] == i)
			{
				for(k=0;k<plithos_items;k++)
				{
					if(ass[k] == nCluster /*&& !isnan(all_distances[j][k]) && all_distances[j][k]!=MAX_DIST*/)
					{
						
						dist+=all_distances[j][k];//findDistance(array, j, k, dimensions, avgRatingsPerUser, metriki);
						//printf("dist==%f\n",dist);
						}
				}
				//printf("ncluster=%d\n",nCluster);
				b[j]=dist/itemsInCluster[nCluster];
			}//end upologismos b gia to j item
			
			max=a[j];
			
			//printf("max is %f\n",max);
			//printf("a[%d] is %f\n",j,a[j]);
			//printf("b[%d] is %f\n",j,b[j]);
			if(isnan(b[j])){
				exit(-1);
			}
			if(a[j] < b[j])
				max=b[j];
				
			if(max != 0.0)
				s[j] = (b[j] - a[j]) / max;	
			else
				s[j] = 0.0;
				
		}
	}
	
	free(a); free(b); free(itemsInCluster);
	return s;
}

int* loydsNEW(double** array, int* assignments, int plithos_items, int dimensions, char distance_metric, int K, int* best_meds, double* avgRatingsPerUser)
{
	int i,j,k, count, different_ass, loops;
	double *centroid, min_dist, sum, avg_centr;
	double dist;
	int *medoids, *new_ass;
	medoids = (int*)malloc(K*sizeof(int));
	centroid = (double*)malloc(dimensions*sizeof(double));
	new_ass = (int*)malloc(plithos_items*sizeof(int));
	printf("\nPerforming loyds Update, please wait");	
	
					int hh,kk;
				for(hh=0;hh<K-1;hh++){
					for(kk=hh+1;kk<K;kk++){
						if(best_meds[hh]==best_meds[kk]){
							printf("idia medoids!!!!!!!!!\n");
							exit(-1);
						}
					}
				}

	
	for(loops=0; loops<3; loops++)
	{
			//printf("loop=%d\n",loops);
			for(i=0;i<K;i++)
			{
				count = 0;
				for(j=0;j<dimensions;j++)
					centroid[j] = 0.0;
					
				for(j=0;j<plithos_items;j++)
				{
					if(assignments[j] == i) //an anikei to item sto cluster
					{
						count++;
						for(k=0;k<dimensions;k++)
							centroid[k] += array[j][k];
					}
				}	
			
				avg_centr = 0.0;
				for(j=0;j<dimensions;j++)
				{
					centroid[j] /= count;
					avg_centr += centroid[j];
				}
				avg_centr = avg_centr/dimensions;
				
				
					
				medoids[i] = -1;  min_dist = MAX_DIST;
				int mpike=0;
				for(j=0; j<plithos_items; j++)
				{
					
					if(assignments[j] == i)
					{
						mpike=1;
						if(distance_metric == 'h') dist = hammingDistance(centroid, array[j], dimensions);
						else if(distance_metric == 'e') dist = eucledianDistance(centroid, array[j], avg_centr, avgRatingsPerUser[j], dimensions);
						else if(distance_metric == 'c') dist = cosineDistance(centroid, array[j], avg_centr, avgRatingsPerUser[j], dimensions);
						//printf("dist is %f\n",dist);
						if(dist<min_dist)
						{
							min_dist = dist;
							medoids[i] = j;
						}
					}
				}
				
				if(!mpike){
					printf("DEN MPIKE\n");
					printf("cluster=%d and plithos_items=%d\n",i,plithos_items);
					for(i=0;i<dimensions;i++){
						printf("centroid[%d]=%f\n",i,centroid[i]);
					}
					exit(-1);
				}
				
			}
			
			int hh,kk;
				for(hh=0;hh<K-1;hh++){
					for(kk=hh+1;kk<K;kk++){
						if(medoids[hh]==medoids[kk]){
							printf("idia medoids medoids[%d]=%d medoids[%d]=%d!!\n",hh,medoids[hh],kk,medoids[kk]);
							exit(-1);
						}
					}
				}

		new_ass = PAM_assignment(array, distance_metric, dimensions, plithos_items, medoids, avgRatingsPerUser, K);
		printf(".");
		different_ass = 0;
		for(i=0;i<plithos_items;i++)
			if(new_ass[i] != assignments[i])
				different_ass = 1;
						
		if(different_ass == 0)
		{
			for(i=0;i<K;i++)
				best_meds[i]=medoids[i];
			free(assignments);
			return new_ass;
		}
		
		for(i=0;i<plithos_items;i++)
			assignments[i] = new_ass[i];
		//printf("META TO PROTO LOOP\n");
	//	for(i=0;i<plithos_items;i++)
		//	printf("ass %d = %d\n", i, assignments[i]);
	}
	return new_ass;
}
int* k_medoids(double** array,double ** all_distances, char distance_metric, int dimensions, int plithos_items, int K,double * avgRatingsPerUser)
{
	int i, *kmeds, thesi_k, j, int_sd;
	double sum_dist, d_r, probability, dist, *min_distances;
	
	kmeds = (int*)malloc(K*sizeof(int));
	min_distances = (double*)malloc(plithos_items*sizeof(double));
	
	kmeds[0] = rand()%plithos_items;
	thesi_k=1;
	
	for(thesi_k=1; thesi_k<K; thesi_k++)
	{
		sum_dist = 0.0;
		for(i=0; i<plithos_items; i++)
		{	 
			min_distances[i] = MAX_DIST;
			for(j=0; j<thesi_k; j++)
			{
				dist=all_distances[kmeds[j]][i];
				// if(distance_metric == 'h')
// 					dist = hammingDistance(array[kmeds[j]], array[i], dimensions);
// 				else if(distance_metric == 't')
// 					dist = array[kmeds[j]][i];
// 				else if(distance_metric == 'e')
// 					dist = eucledianDistance(array[kmeds[j]], array[i], 0.0,avgRatingsPerUser[i], dimensions);
// 					//eucledianDistance(centroid, array[j], avg_centr, avgRatingsPerUser[j], dimensions);
// 				else if(distance_metric == 'c')
// 					dist = cosineDistance(array[kmeds[j]], array[i], 0.0,avgRatingsPerUser[i], dimensions);
				
					
				if(dist<min_distances[i])	
					min_distances[i] = dist;
			}
			sum_dist = sum_dist + min_distances[i]*min_distances[i];
		}
		
		int_sd = (int)(sum_dist);
		d_r = rand()%int_sd;

		probability = 0.0;
		for(i=0; i<plithos_items; i++)
		{
			probability = probability + min_distances[i]*min_distances[i];
			if(probability >= d_r)
			{
				kmeds[thesi_k] = i;
				break;
			}	
		}		
	}
	
	free(min_distances);
	return kmeds;
}
/************************ Clustering Functions END *********************/



