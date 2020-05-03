/**
 *
 *      @file usr.c
 *      @brief UltraShape clone module for the molecular library
 *             acortes notes: I have removed the global array for centroids
 *             since we can have multiple calculations
 *
 *      @author Javier Klett <jklett@cbm.uam.es>
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *      @date 10/05/2012 (acortes repackage in the global library)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */

#include "usr.h"

/**
 *      @brief Perform the dot product in a high dimensional space
 *      @author Javier Klett <jklett@cbm.uam.es>
 *
 *      @param v1 Vector 1
 *      @param v2 Vector 2
 *      @return The product
 *
 */
double scalarProduct(double * v1,double * v2){

	double scaProd =0;
	int i;

	for (i=0;i<DIMENSION-1;i++){
		scaProd = scaProd + v1[i]*v2[i];
	}

	scaProd = sqrt(scaProd);
	return scaProd;
}

/**
 *      @brief Perform the cross product in a high dimensional space
 *      @author Javier Klett <jklett@cbm.uam.es>
 *
 *      @param v1 Vector 1
 *      @param v2 Vector 2
 *      @param crossPro Result
 *      @return The result also!
 *
 */
double * crossProduct(double * v1,double * v2, double * crossPro){

	crossPro[0] = v1[1]*v2[2]-v1[2]*v2[1];
/** 	Megabug of shame! Fixed on 13/06/2012 acortes */
/*	crossPro[1] = v1[0]*v2[2]-v1[2]*v2[0];*/
	crossPro[1] = v1[2]*v2[0]-v1[0]*v2[2];
	crossPro[2] = v1[0]*v2[1]-v1[1]*v2[0];

	return crossPro;

}

/**
 *      @brief Get the first centroid for ElectroShape method
 *      @author Javier Klett <jklett@cbm.uam.es>
 *
 *      @param xs Vector of X coordinates
 *      @param ys Vector of Y coordinates 
 *      @param zs Vector of Z coordinates 
 *      @param qs Vector of partial charges
 *      @param n_atoms Number of atoms of the molecule
 *      @param mycentroids Resulting centroids
 *      @return 0 on success
 *
 */
int getCentroid1(float * xs, float * ys, float * zs, float * qs, int n_atoms, double ***mycentroids){

	int i;
	double sumX = 0 , sumY = 0, sumZ = 0, sumQ = 0;
	double **centroids = NULL;

	centroids = *mycentroids;

	for( i = 0; i < n_atoms; i++ ){

		sumX = sumX + xs[i];
		sumY = sumY + ys[i];
		sumZ = sumZ + zs[i];
		sumQ = sumQ + SCALING_CHARGE*qs[i];
	}

	centroids[0][0] = (double) sumX / (double) n_atoms;
	centroids[0][1] = (double) sumY / (double) n_atoms;
	centroids[0][2] = (double) sumZ / (double) n_atoms;
	centroids[0][3] = (double) sumQ / (double) n_atoms;

	return 0;

}

/**
 *      @brief Get the second centroid for ElectroShape method
 *      @author Javier Klett <jklett@cbm.uam.es>
 *      
 *      @param xs Vector of X coordinates 
 *      @param ys Vector of Y coordinates 
 *      @param zs Vector of Z coordinates
 *      @param qs Vector of partial charges
 *      @param n_atoms Number of atoms of the molecule
 *      @param mycentroids Resulting centroids
 *      @return 0 on success
 *      
 */
int getCentroid2(float * xs, float * ys, float * zs, float * qs, int n_atoms, double ***mycentroids){

	int i;
	int id_biggest=0;

	double dst_tmp = 0, dst_biggest = 0;
        double **centroids = NULL;

        centroids = *mycentroids;

	
	for (i = 0; i< n_atoms; i++){

		dst_tmp = sqrt((xs[i] - centroids[0][0])*(xs[i] - centroids[0][0]) + (ys[i] - centroids[0][1])*(ys[i] - centroids[0][1]) + (zs[i] - centroids[0][2])*(zs[i] - centroids[0][2]) + (SCALING_CHARGE*qs[i] - centroids[0][3])*(SCALING_CHARGE*qs[i] - centroids[0][3]));

		if( dst_tmp > dst_biggest ){

			dst_biggest = dst_tmp;	
			id_biggest = i;

		}
	}

	centroids[1][0] = xs[id_biggest];
	centroids[1][1] = ys[id_biggest];
	centroids[1][2] = zs[id_biggest];
	centroids[1][3] = SCALING_CHARGE*qs[id_biggest];

	return 0;
}

/**
 *      @brief Get the third centroid for ElectroShape method
 *      @author Javier Klett <jklett@cbm.uam.es>
 *      
 *      @param xs Vector of X coordinates 
 *      @param ys Vector of Y coordinates 
 *      @param zs Vector of Z coordinates
 *      @param qs Vector of partial charges
 *      @param n_atoms Number of atoms of the molecule
 *      @param mycentroids Resulting centroids
 *      @return 0 on success
 *      
 */
int getCentroid3(float * xs, float * ys, float * zs, float * qs, int n_atoms, double ***mycentroids){

	int i;
	int id_biggest=0;

	double dst_tmp = 0, dst_biggest = 0;
        double **centroids = NULL;

        centroids = *mycentroids;

	
	for (i = 0; i< n_atoms; i++){

		dst_tmp = sqrt((xs[i] - centroids[1][0])*(xs[i] - centroids[1][0]) + (ys[i] - centroids[1][1])*(ys[i] - centroids[1][1]) + (zs[i] - centroids[1][2])*(zs[i] - centroids[1][2])
										+ (SCALING_CHARGE*qs[i] - centroids[1][3])*(SCALING_CHARGE*qs[i] - centroids[1][3]));

		/* printf("dst: %7.3f ->> [%7.3f,%7.3f,%7.3f,%7.3f] vs [%7.3f,%7.3f,%7.3f,%7.3f]\n",dst_tmp,xs[i],ys[i],zs[i],qs[i],centroids[1][0],centroids[1][1],centroids[1][2],centroids[1][3]); */

		if( dst_tmp > dst_biggest ){

			dst_biggest = dst_tmp;	
			id_biggest = i;

			/*printf("%f %d\n",dst_biggest,id_biggest);*/

		}
	}

	/*printf("%d, %7.3f  %7.3f  %7.3f  %7.3f\n",id_biggest,xs[id_biggest],ys[id_biggest],zs[id_biggest],qs[id_biggest]*SCALING_CHARGE);*/

	centroids[2][0] = xs[id_biggest];
	centroids[2][1] = ys[id_biggest];
	centroids[2][2] = zs[id_biggest];
	centroids[2][3] = qs[id_biggest]*SCALING_CHARGE;

	return 0;
}

/**
 *      @brief Get the forth and fifth centroids for ElectroShape method
 *      @author Javier Klett <jklett@cbm.uam.es>
 *      
 *      @param xs Vector of X coordinates 
 *      @param ys Vector of Y coordinates 
 *      @param zs Vector of Z coordinates
 *      @param qs Vector of partial charges
 *      @param n_atoms Number of atoms of the molecule
 *      @param mycentroids Resulting centroids
 *      @return 0 on success
 *      
 */
int getCentroid4_5(float * xs, float * ys, float * zs, float * qs, int n_atoms, double ***mycentroids){

	double A[DIMENSION-1],B[DIMENSION-1],Cvec[DIMENSION-1];
	double moduleA, moduleCrossAB, crossAB[DIMENSION-1],K;
	int i;
	double tmp[4];
	double qMin=0, qMax=0;
        double **centroids = NULL;

        centroids = *mycentroids;
	

	for (i =0;i < DIMENSION-1;i++){
		A[i] = centroids[1][i]-centroids[0][i];
		B[i] = centroids[2][i]-centroids[0][i];
	}

	moduleA = sqrt((centroids[1][0]-centroids[0][0])*(centroids[1][0]-centroids[0][0])+(centroids[1][1]-centroids[0][1])*(centroids[1][1]-centroids[0][1])+(centroids[1][2]-centroids[0][2])*(centroids[1][2]-centroids[0][2])+(centroids[1][3]-centroids[0][3])*(centroids[1][3]-centroids[0][3]));
	crossProduct(A,B,crossAB);

	moduleCrossAB = scalarProduct(crossAB,crossAB);

  K	= moduleA/(2*moduleCrossAB);

	tmp[0] = K*crossAB[0];
	tmp[1] = K*crossAB[1];
	tmp[2] = K*crossAB[2];

	for (i=0; i < n_atoms; i++){
		if (SCALING_CHARGE*qs[i] > qMax){
			qMax = SCALING_CHARGE*qs[i];
		}else if(SCALING_CHARGE*qs[i] < qMin){
			qMin = SCALING_CHARGE*qs[i];
		}
	}

	centroids[3][0] = centroids[0][0] + tmp[0];
	centroids[3][1] = centroids[0][1] + tmp[1];
	centroids[3][2] = centroids[0][2] + tmp[2];
	centroids[3][3] = qMax;

	centroids[4][0] = centroids[0][0] + tmp[0];
	centroids[4][1] = centroids[0][1] + tmp[1];
	centroids[4][2] = centroids[0][2] + tmp[2];
	centroids[4][3] = qMin;

	return 0;
	
}

/**
 *      @brief Get all centroids for ElectroShape method
 *      @author Javier Klett <jklett@cbm.uam.es>
 *      
 *      @param mol MOL2 with the molecule
 *      @param mycentroids Resulting centroids
 *      @return 0 on success
 *      
 */
int getCentroids( MOL2 *mol, double ***mycentroids){

	double **centroids = NULL;
	centroids = *mycentroids;

        getCentroid1( mol->x, mol->y, mol->z, mol->pcharges, mol->n_atoms, &centroids);
        getCentroid2( mol->x, mol->y, mol->z, mol->pcharges, mol->n_atoms, &centroids);
        getCentroid3( mol->x, mol->y, mol->z, mol->pcharges, mol->n_atoms, &centroids);
        getCentroid4_5( mol->x, mol->y, mol->z, mol->pcharges, mol->n_atoms, &centroids);

/*        printf("Centroid 1 = (%.2f, %.2f, %.2f, %.2f)\n",   centroids[0][0], centroids[0][1],centroids[0][2],centroids[0][3]);
        printf("Centroid 2 = (%.2f, %.2f, %.2f, %.2f)\n",   centroids[1][0], centroids[1][1],centroids[1][2],centroids[1][3]);
        printf("Centroid 3 = (%.2f, %.2f, %.2f, %.2f)\n",   centroids[2][0], centroids[2][1],centroids[2][2],centroids[2][3]);
        printf("Centroid 4 = (%.2f, %.2f, %.2f, %.2f)\n",   centroids[3][0], centroids[3][1],centroids[3][2],centroids[3][3]);
        printf("Centroid 5 = (%.2f, %.2f, %.2f, %.2f)\n\n", centroids[4][0], centroids[4][1],centroids[4][2],centroids[4][3]);*/

        return 0;
}

/**
 *
 *	@brief Print out calculated centroids for debug 
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *	@param mol MOL2 with the molecule
 *	@return 0 on success
 *
 */
int dump_centroids(MOL2 *mol){

        double **centroids = NULL;
        int i = 0;

        centroids = (double **) calloc(NUM_CENTROIDS,sizeof(double **));

        for(i = 0; i < NUM_CENTROIDS; i++)
        {
          centroids[i] = (double *) calloc(DIMENSION, sizeof(double *));
        }

        getCentroid1( mol->x, mol->y, mol->z, mol->pcharges, mol->n_atoms, &centroids);
        getCentroid2( mol->x, mol->y, mol->z, mol->pcharges, mol->n_atoms, &centroids);
        getCentroid3( mol->x, mol->y, mol->z, mol->pcharges, mol->n_atoms, &centroids);
        getCentroid4_5( mol->x, mol->y, mol->z, mol->pcharges, mol->n_atoms, &centroids);

        printf("Centroid 1 = (%.2f, %.2f, %.2f, %.2f)\n",   centroids[0][0], centroids[0][1],centroids[0][2],centroids[0][3]);
        printf("Centroid 2 = (%.2f, %.2f, %.2f, %.2f)\n",   centroids[1][0], centroids[1][1],centroids[1][2],centroids[1][3]);
        printf("Centroid 3 = (%.2f, %.2f, %.2f, %.2f)\n",   centroids[2][0], centroids[2][1],centroids[2][2],centroids[2][3]);
        printf("Centroid 4 = (%.2f, %.2f, %.2f, %.2f)\n",   centroids[3][0], centroids[3][1],centroids[3][2],centroids[3][3]);
        printf("Centroid 5 = (%.2f, %.2f, %.2f, %.2f)\n\n", centroids[4][0], centroids[4][1],centroids[4][2],centroids[4][3]);

        for(i = 0; i < NUM_CENTROIDS; i++)
        {
        	free(centroids[i]);
        }
	free(centroids);


	return 0;
}


/**
 *      @brief Calculate Momenta for ElectroShape method
 *      @author Javier Klett <jklett@cbm.uam.es>
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *
 *      @param mol MOL2 with the molecule
 *      @param centroids Centroids
 *      @param myM Resulting momenta
 *      @return 0 on success
 *
 */
int  calculateMoments(MOL2 *mol, double **centroids, float ** myM){

	int i;
	float d0, d1, d2, d3, d4;
	float *M = NULL;

	double sum_d0_1=0,sum_d0_2=0,sum_d0_3=0;
	double sum_d1_1=0,sum_d1_2=0,sum_d1_3=0;
	double sum_d2_1=0,sum_d2_2=0,sum_d2_3=0;
	double sum_d3_1=0,sum_d3_2=0,sum_d3_3=0;
	double sum_d4_1=0,sum_d4_2=0,sum_d4_3=0;
	
	M = *myM;
	
	for(i=0; i<mol->n_atoms; i++){

		d0 = sqrt( (mol->x[i]-centroids[0][0])*(mol->x[i]-centroids[0][0]) +
			(mol->y[i]-centroids[0][1])*(mol->y[i]-centroids[0][1]) + 
			(mol->z[i]-centroids[0][2])*(mol->z[i]-centroids[0][2]) + 
			((mol->pcharges[i]*SCALING_CHARGE)-centroids[0][3])*
			((mol->pcharges[i]*SCALING_CHARGE)-centroids[0][3]));
						
		d1 = sqrt( (mol->x[i]-centroids[1][0])*(mol->x[i]-centroids[1][0]) +
			(mol->y[i]-centroids[1][1])*(mol->y[i]-centroids[1][1]) + 
			(mol->z[i]-centroids[1][2])*(mol->z[i]-centroids[1][2]) + 
			((mol->pcharges[i]*SCALING_CHARGE)-centroids[1][3])*
			((mol->pcharges[i]*SCALING_CHARGE)-centroids[1][3]));

		d2 = sqrt( (mol->x[i]-centroids[2][0])*(mol->x[i]-centroids[2][0]) +
			(mol->y[i]-centroids[2][1])*(mol->y[i]-centroids[2][1]) + 
			(mol->z[i]-centroids[2][2])*(mol->z[i]-centroids[2][2]) + 
			((mol->pcharges[i]*SCALING_CHARGE)-centroids[2][3])*
			((mol->pcharges[i]*SCALING_CHARGE)-centroids[2][3]));

		d3 = sqrt( (mol->x[i]-centroids[3][0])*(mol->x[i]-centroids[3][0]) +
			(mol->y[i]-centroids[3][1])*(mol->y[i]-centroids[3][1]) + 
			(mol->z[i]-centroids[3][2])*(mol->z[i]-centroids[3][2]) + 
			((mol->pcharges[i]*SCALING_CHARGE)-centroids[3][3])*
			((mol->pcharges[i]*SCALING_CHARGE)-centroids[3][3]) );

		d4 = sqrt( (mol->x[i]-centroids[4][0])*(mol->x[i]-centroids[4][0]) +
			(mol->y[i]-centroids[4][1])*(mol->y[i]-centroids[4][1]) + 
			(mol->z[i]-centroids[4][2])*(mol->z[i]-centroids[4][2])  + 
			((mol->pcharges[i]*SCALING_CHARGE)-centroids[4][3])*
			((mol->pcharges[i]*SCALING_CHARGE)-centroids[4][3]));
						

/* Acumulators for Caluculating moments*/
		sum_d0_1 += d0;
		sum_d0_2 += d0*d0; 
		sum_d0_3 += d0*d0*d0;

		sum_d1_1 += d1;
		sum_d1_2 += d1*d1; 
		sum_d1_3 += d1*d1*d1;

		sum_d2_1 += d2;
		sum_d2_2 += d2*d2; 
		sum_d2_3 += d2*d2*d2;

		sum_d3_1 += d3;
		sum_d3_2 += d3*d3; 
		sum_d3_3 += d3*d3*d3;

		sum_d4_1 += d4;
		sum_d4_2 += d4*d4; 
		sum_d4_3 += d4*d4*d4;
	}

/**
 *	On 13/06/2012. acortes. Fixed issues with std and third moment (mostly sqrt and cubic root)
 *	This is 100% compatible with the original method.
 */
	M[0]  = (float) sum_d0_1/mol->n_atoms;
	M[1]  = (float) (sum_d0_2/mol->n_atoms) - M[0]*M[0];
	M[2]  = (float) (sum_d0_3/mol->n_atoms) - 3*M[0]*M[1] - M[0]*M[0]*M[0];
	M[1] = sqrt(M[1]);
	
	M[3]  = (float) sum_d1_1/mol->n_atoms;
	M[4]  = (float) (sum_d1_2/mol->n_atoms) - M[3]*M[3];
	M[5]  = (float) (sum_d1_3/mol->n_atoms) - 3*M[3]*M[4] - M[3]*M[3]*M[3];
	M[4] = sqrt(M[4]);
	
	M[6]  = (float) sum_d2_1/mol->n_atoms;
	M[7]  = (float) (sum_d2_2/mol->n_atoms) - M[6]*M[6];
	M[8]  = (float) (sum_d2_3/mol->n_atoms) - 3*M[6]*M[7] - M[6]*M[6]*M[6];
	M[7] = sqrt(M[7]);

	M[9]  = (float) sum_d3_1/mol->n_atoms;
	M[10] = (float) (sum_d3_2/mol->n_atoms) - M[9]*M[9];
	M[11] = (float) (sum_d3_3/mol->n_atoms) - 3*M[9]*M[10] - M[9]*M[9]*M[9];
	M[10] = sqrt(M[10]);

	M[12] = (float) sum_d4_1/mol->n_atoms;
	M[13] = (float) (sum_d4_2/mol->n_atoms) - M[12]*M[12];
	M[14] = (float) (sum_d4_3/mol->n_atoms) - 3*M[12]*M[13] - M[12]*M[12]*M[12];
	M[13] = sqrt(M[13]);

	if( M[2] < 0)
		M[2] = -1.0*(pow(-1.*M[2],1.0/3.0));
	else
		M[2] = pow(M[2],1./3.);

        if( M[5] < 0)
                M[5] = -1.0*(pow(-1.0*M[5],1.0/3.0));
        else
                M[5] = pow(M[5],1.0/3.0);

        if( M[8] < 0)
                M[8] = -1.0*(pow(-1.0*M[8],1.0/3.0));
        else
                M[8] = pow(M[8],1.0/3.0);

        if( M[11] < 0)
                M[11] = -1.0*(pow(-1.0*M[11],1.0/3.0));
        else
                M[11] = pow(M[11],1.0/3.0);


        if( M[14] < 0)
                M[14] = -1.0*(pow(-1.0*M[14],1.0/3.0));
        else
                M[14] = pow(M[14],1.0/3.0);


	return 0;
}


/**
 *      @brief Calculate the distance between two sets of momenta
 *      @author Javier Klett <jklett@cbm.uam.es>
 *
 *      @param u First momenta
 *      @param v Second momenta
 *      @return The distance
 *
 */
double shapeDistance(float * u, float * v){

	int i;
	float tmp, dist;

	tmp=0;

	for (i=0; i < NUM_MOMENTS ;i++){
		tmp += fabs(u[i]-v[i]);
	}

	dist = (float) 1/((float) 1 + (float)((float)1/(float)NUM_MOMENTS)*tmp);

	return dist;
}

/**
 *      @brief Dump the momenta for staroge
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *
 *      @param mol MOL2 with the molecule
 *      @return null
 *
 */
void dump_moments(MOL2 *mol1)
{
        float *moments_a = NULL;
        double **centroids_a = NULL;
	int i = 0;

        centroids_a = (double **) calloc(NUM_CENTROIDS,sizeof(double **));

        for(i = 0; i < NUM_CENTROIDS; i++)
        {
          centroids_a[i] = (double *) calloc(DIMENSION, sizeof(double *));
	}

        getCentroids(mol1,&centroids_a);
        moments_a = (float *) calloc(NUM_MOMENTS,sizeof(float*));
        calculateMoments(mol1, centroids_a, &moments_a);

	for(i = 0; i < NUM_MOMENTS; i++)
	{
		printf("%.3f;",moments_a[i]);
	}

        for(i = 0; i < NUM_CENTROIDS; i++)
        {
                free(centroids_a[i]);
        }

        free(centroids_a);
        free(moments_a);


}

/**
 *      @brief Compare two molecules using ElectroShape method
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *
 *      @param mol1 MOL2 with the molecule
 *      @param mol1 MOL2 with the molecule
 *      @return Similarity
 *
 */
float USR_compare_molecules(MOL2 *mol1, MOL2 *mol2)
{
        float *moments_a = NULL, *moments_b = NULL;
        double **centroids_a = NULL, **centroids_b = NULL;
	float dist = 0.0f;
        int i = 0;

        centroids_a = (double **) calloc(NUM_CENTROIDS,sizeof(double **));
        centroids_b = (double **) calloc(NUM_CENTROIDS,sizeof(double **));

        for(i = 0; i < NUM_CENTROIDS; i++)
	{
          centroids_a[i] = (double *) calloc(DIMENSION, sizeof(double *));
          centroids_b[i] = (double *) calloc(DIMENSION, sizeof(double *));
	}

        getCentroids(mol1,&centroids_a);
        getCentroids(mol2,&centroids_b);

        moments_a = (float *) calloc(NUM_MOMENTS,sizeof(float*));
        moments_b = (float *) calloc(NUM_MOMENTS,sizeof(float*));

        calculateMoments(mol1, centroids_a, &moments_a);
        calculateMoments(mol2, centroids_b, &moments_b);

	dist = shapeDistance(moments_a,moments_b);

        for(i = 0; i < NUM_CENTROIDS; i++)
        {
		free(centroids_a[i]);
                free(centroids_b[i]);
        }

	free(centroids_a);
	free(centroids_b);
	free(moments_a);
	free(moments_b);

	return dist;
}


