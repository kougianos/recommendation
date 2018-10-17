#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "previousFunctions.h"

void ParseFile(char* inputFile, int* user, int* movie, int* p, int* userGivesP);
void fillMovieRatings(char* inputFile, double** array, int userCount, int movieCount, int userGivesP);
void fillAvgRatingsPerUser(double **movieRatings, double *array, int userCount, int movieCount);
void clustering_recommendation(double **movieRatings, int userCount, int movieCount, double* avgRatingsPerUser, char distance_metric, int P, int argc);

typedef struct{
	int index;
	double rating;
}crossStruct;

int main(int argc, char *argv[])
{
	srand(time(NULL));
	int userCount, movieCount, P, userGivesP, i, j;
	double **movieRatings, *avgRatingsPerUser;
	char inputFile[20];
	strcpy(inputFile, argv[2]);
	ParseFile(inputFile, &userCount, &movieCount, &P, &userGivesP);
	
	movieRatings = (double**)malloc(userCount*sizeof(double*));
	for(i = 0; i<userCount; i++)
		movieRatings[i] = (double*)malloc(movieCount*sizeof(double));
	avgRatingsPerUser = (double*)malloc(userCount*sizeof(double));
		
	printf("Movies = %d    Users = %d     P = %d   givesP = %d", movieCount, userCount, P, userGivesP);
	fillMovieRatings(inputFile, movieRatings, userCount, movieCount, userGivesP);
	fillAvgRatingsPerUser(movieRatings, avgRatingsPerUser, userCount, movieCount);
	clustering_recommendation(movieRatings,userCount,movieCount,avgRatingsPerUser,'e', P, argc);
//	clustering_recommendation(movieRatings,userCount,movieCount,avgRatingsPerUser,'c', P, argc);
//	clustering_recommendation(movieRatings,userCount,movieCount,avgRatingsPerUser,'h', P, argc);

}


void ParseFile(char* inputFile, int* userCount, int* movieCount, int* p, int* userGivesP)
{	
	char c;
	int maxUser, maxMovie, num, user, movie;
	FILE *infile;
	infile = fopen(inputFile, "r");
	fscanf(infile, "%c", &c);
	if(c == 'P') //to input file dinei P
	{
		*userGivesP = 1;
		fscanf(infile, "%c", &c); 
		fscanf(infile, "%d", *(&p));
	}
	else //default P
	{
		*userGivesP = 0;
		fseek(infile, 0, SEEK_SET);
		*p = 20;
	}
		
	maxUser = -1; maxMovie = -1;
	while(!feof(infile))
	{
		fscanf(infile, "%d %d %d", &user, &movie, &num);
		if(user >= maxUser)
			maxUser = user;
		if(movie >= maxMovie)
			maxMovie = movie;
	}
	*userCount = maxUser;
	*movieCount = maxMovie;
	fclose(infile);
	return;
}

void fillMovieRatings(char* inputFile, double** array, int userCount, int movieCount, int userGivesP)
{
	FILE *infile;
	infile = fopen(inputFile, "r");
	char xar;
	int i, j, user, movie;
	double rating;
	for(i=0; i<userCount; i++) //arxikopoiisi olwn twn stoixeiwn tou pinaka sto 0
	{
		for(j=0; j<movieCount; j++)
		{
			array[i][j] = 0.0;
		//	printf("\narray[%d][%d] = %d", i, j, array[i][j]);
		}
	}
	if(userGivesP == 1)
		fscanf(infile, "%c%c	%lf", &xar, &xar, &rating); //diavazw tin proti grammi pou de mou xreiazetai
	
	while(!feof(infile))
	{
		fscanf(infile, "%d	%d	%lf", &user, &movie, &rating);
		array[user-1][movie-1] = rating;
	}
	
	fclose(infile);
	return;
}

void fillAvgRatingsPerUser(double **movieRatings, double *array, int userCount, int movieCount)
{
	int i, j;
	double sum, plithos;
	for(i=0; i<userCount; i++)
	{
		plithos = 0; sum = 0;
		for(j=0; j<movieCount; j++)
		{
			if(movieRatings[i][j] > 0.0) //diladi o xristis exei dosei vathmologia
			{
				plithos++;
				sum += movieRatings[i][j];
			}
		}
		if(plithos>0)
			array[i] = sum / plithos;
		else
			array[i] = 0.0;
	}
	return;
}
	
void clustering_recommendation(double **movieRatings, int userCount, int movieCount, double* avgRatingsPerUser, char distance_metric, int P, int argc)
{
	printf("\n\nInside function clustering_recommendation, calculating distances, please wait");
	fflush(0);
	int K, *con_array, *pam_ass, i, *loyds_ass , j, k, l;
	double *siloueta, oliki_siloueta, new_oliki_siloueta;
			
	con_array = (int*)malloc(userCount*sizeof(int));
	pam_ass = (int*)malloc(userCount*sizeof(int));
	loyds_ass = (int*)malloc(userCount*sizeof(int));
	siloueta = (double*)malloc(userCount*sizeof(double));
	
	double ** all_distances = (double**)malloc(userCount*(sizeof(double*)));
	for(i=0;i<userCount;i++){
		all_distances[i]=(double *)malloc(userCount*sizeof(double));
		if(i%250==0) printf("."); fflush(0);
		for(j=i;j<userCount;j++){
			all_distances[i][j] = findDistance(movieRatings, i, j, movieCount, avgRatingsPerUser, distance_metric);
//				printf("distance[%d][%d]=%.15lf\n",i,j,all_distances[i][j]);
			// if(all_distances[i][j]  > 0.0){
// 					printf("distance[%d][%d]=%.15lf\n",i,j,all_distances[i][j]);					
// 				}
		}
	}
	
	
	
	
	for(i=0;i<userCount;i++){    //logw simmetrikotitas twn pinakwn
		for(j=i;j<userCount;j++){
			all_distances[j][i]=all_distances[i][j];
			
		}
	} 
	printf("\nAll distances between users calculated succesfully!");
	
	printf("\nCalculating optimal K, please wait.\n");
	
	int bestK = 0; double bestSil = -1.0;
	for(K=(int)(userCount/P - 2); K<(int)(userCount/P + 3); K++)
	{	
	
	//printf("Processing K=%d\n",K);
		if(distance_metric == 'e')
			{
			//con_array=k_medoids(movieRatings,all_distances,'e',movieCount, userCount,K,avgRatingsPerUser);
				con_array = concentrate(movieRatings,all_distances, 'e', 'e', movieCount, userCount, avgRatingsPerUser, K);  
			}
		if(distance_metric == 'c'){
			//con_array=k_medoids(movieRatings,all_distances,'c',movieCount, userCount,K,avgRatingsPerUser);
			con_array = concentrate(movieRatings,all_distances, 'c', 'e', movieCount, userCount, avgRatingsPerUser, K);
			}
		if(distance_metric == 'h'){
		//con_array=k_medoids(movieRatings,all_distances,'h',movieCount, userCount,K,avgRatingsPerUser);
			con_array = concentrate(movieRatings,all_distances, 'h', 'h', movieCount, userCount, avgRatingsPerUser, K);
			}
		//printf("Clusters ready!\n");
// 		printf("\nCONCENTRATE COMPLETED");
//	 	printf("for K=%d \n", K);
// 	for(i=0; i<K; i++)
//  		printf("\ncon_array[%d] = %d", i, con_array[i]);
		// printf("The medoid of 0 is %d\n",con_array[0]);
// 		printf("The medoid of 123 is %d\n",con_array[123]);
		//exit(1);
		
		pam_ass = PAM_assignment(movieRatings, distance_metric, movieCount, userCount, con_array, avgRatingsPerUser, K); // printf("\npam ass completed");
		// for(i=0;i<userCount;i++)
// 		printf("pam_ass  %d  =  %d\n", i, pam_ass[i]);
		loyds_ass = loydsNEW(movieRatings, pam_ass, userCount, movieCount, distance_metric, K, con_array, avgRatingsPerUser);//	printf("\nloyds ass completed");
		siloueta = Silhouette(movieRatings,all_distances, userCount, movieCount, distance_metric, con_array, pam_ass, avgRatingsPerUser, K); //printf("\nsiloueta completed");
		
		oliki_siloueta = 0.0;
		for(i=0; i<userCount; i++){
			
			oliki_siloueta += siloueta[i];
			//printf("silouette[%d]=%f\n",i,siloueta[i]);	
		}
		oliki_siloueta /= userCount;
		
		printf("\nSilhouette %.3lf gia K=%d", oliki_siloueta, K);
		
		if(oliki_siloueta >= bestSil)
		{
			bestSil = oliki_siloueta;
			bestK = K;
		}	
	}
	printf("\nPerforming Clustering with optimal K=%d", bestK);
	if(distance_metric == 'e')
		con_array = concentrate(movieRatings,all_distances, 'e', 'e', movieCount, userCount, avgRatingsPerUser, bestK);
	if(distance_metric == 'c')
		con_array = concentrate(movieRatings,all_distances, 'c', 'e', movieCount, userCount, avgRatingsPerUser, bestK);
	if(distance_metric == 'h')
		con_array = concentrate(movieRatings,all_distances, 'h', 'h', movieCount, userCount, avgRatingsPerUser, bestK);
		
	pam_ass = PAM_assignment(movieRatings, distance_metric, movieCount, userCount, con_array, avgRatingsPerUser, bestK);
	loyds_ass = loydsNEW(movieRatings, pam_ass, userCount, movieCount, distance_metric, bestK, con_array, avgRatingsPerUser);
 			
	printf("\nClustering done with %d Clusters, Silhouette = %.3lf", bestK, bestSil);
	
	
	double sumRating, finalRating, **unratedItems;
	int sum;
	unratedItems = (double**)malloc(userCount*sizeof(double*));
	for(i=0;i<userCount;i++)
	{
		unratedItems[i] = (double*)malloc(movieCount*sizeof(double));
		for(j=0;j<movieCount;j++)
			unratedItems[i][j] = -1.0;
	}
	
	printf("\nPredicting ratings of unrated items per User");
	for(i=0; i<K; i++)
	{
		for(j=0; j<userCount; j++)
		{
			if(loyds_ass[j] == i) //an o xristis j anikei sto cluster i
			{
				for(k=0; k<movieCount; k++)
				{
					if(movieRatings[j][k] == 0) //an o xristis j den exei dosei rating stin k tainia
					{
						sumRating = 0.0; sum = 0;
						for(l=0;l<userCount;l++)
						{
							if(loyds_ass[l] == i) //an o xristis l aniksei sto cluster i
							{
								if(movieRatings[l][k] != 0) // an o xristis l exei dosei rating stin k tainia
								{
									sumRating += movieRatings[l][k] - avgRatingsPerUser[l];
									sum++;
								}
							}
						}
						unratedItems[j][k] = (sumRating / sum) + avgRatingsPerUser[j];
						if(unratedItems[j][k] >= 5.0) unratedItems[j][k] = 5.0;
						if(unratedItems[j][k] <= 1.0) unratedItems[j][k] = 1.0;
						//printf("\nUser %d gives estimated rating of %.2f in movie %d", j, movieRatings[j][k], k);
					}
				}
			}
		}
		if(i%20 == 0) printf(".");
	}
	//printing 5 best movies per user to output file
	int maxMovie;
	double maxRating;
	char output[30];
	if(distance_metric == 'e')
		strcpy(output, "Euclidean Clustering.txt");
	if(distance_metric == 'c')
		strcpy(output, "Cosine Clustering.txt");
	if(distance_metric == 'h')
		strcpy(output, "Hamming Clustering.txt");
	FILE *outfile;
	outfile = fopen(output, "w");
	if(distance_metric == 'e')
		fprintf(outfile, "Euclidean Clustering");
	if(distance_metric == 'c')
		fprintf(outfile, "Cosine Clustering");
	if(distance_metric == 'h')
		fprintf(outfile, "Hamming Clustering");
	for(i=0;i<userCount;i++)
	{
		fprintf(outfile,"\nUser %d    ", i+1);
		for(k=0;k<5;k++)
		{
			maxRating = -1.0;  maxMovie = -1;
			for(j=0;j<movieCount;j++)
			{
				if(unratedItems[i][j] >= maxRating)
				{
					maxMovie = j+1;
					maxRating = unratedItems[i][j];
					unratedItems[i][j] = -1.0; //gia na min upologisei pali tin idia tainia os best
				}
			}
			fprintf(outfile,"Movie %d	",maxMovie);
		}
	}
	printf("\nFunction clustering_recommendation completed!	\n");			
	
	if(argc > 3) //cross validation
	{
		printf("\nCross Validation in progress");
		int ratedItemsCount = 0;
		for(i=0;i<movieCount;i++)  //theoroume oti oloi oi xristes exoun vathmologisei ton idio arithmo tainiwn
		{
			if(movieRatings[0][i] > 0)
				ratedItemsCount++;
		}
		crossStruct **myCrossStruct;
		myCrossStruct = (crossStruct**)malloc(userCount*sizeof(crossStruct*));
		for(i=0;i<userCount;i++)
			myCrossStruct[i] = (crossStruct*)malloc(ratedItemsCount*sizeof(crossStruct));
		
		//fill myCrossStruct
		for(i=0;i<userCount;i++)
		{
			for(j=0;j<ratedItemsCount;j++)
			{
				for(k=0;k<movieCount;k++)
				{
					if(movieRatings[i][k] > 0) //an o xristis exei dosei rating sto item
					{
						myCrossStruct[i][j].index = k;
						myCrossStruct[i][j].rating = movieRatings[i][k];
					}
				}
			}
		}
		int itemForTest;
		//midenizoume tixaia to 10% twn ratings gia kathe user
		for(j=0;j<userCount;j++)
		{
			for(k=0;k<ratedItemsCount/10;k++)
			{
				itemForTest = rand()%ratedItemsCount;
				myCrossStruct[j][itemForTest].rating = 0.0;
			}
		}
		double mae = 0.0;
		for(i=0; i<K; i++)
		{
			for(j=0; j<userCount; j++)
			{
				if(loyds_ass[j] == i) //an o xristis j anikei sto cluster i
				{
					for(k=0; k<ratedItemsCount; k++)
					{
						if(myCrossStruct[j][k].rating == 0) //an o xristis j den exei dosei rating stin k tainia
						{
							sumRating = 0.0; sum = 0;
							for(l=0;l<userCount;l++)
							{
								if(loyds_ass[l] == i) //an o xristis l aniksei sto cluster i
								{
									if(movieRatings[l][k] != 0) // an o xristis l exei dosei rating stin k tainia
									{
										sumRating += movieRatings[l][k] - avgRatingsPerUser[l];
										sum++;
									}
								}
							}
							myCrossStruct[j][k].rating = (sumRating / sum) + avgRatingsPerUser[j];
							if(myCrossStruct[j][k].rating >= 5.0) myCrossStruct[j][k].rating = 5.0;
							if(myCrossStruct[j][k].rating <= 1.0) myCrossStruct[j][k].rating = 1.0;
							//printf("\nUser %d gives estimated rating of %.2f in movie %d", j, movieRatings[j][k], k);
						}
					}
				}
			}
		}
		
	}
		
	return;
}
		
	




