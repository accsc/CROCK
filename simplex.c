/**
 *
 *      @file simplex.c
 *      @brief Rigid-body downhill SIMPLEX or amoeba algorithm for minimization in grid
 *
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *      @date 01/10/2011
 *
 *      This code is based on an old Fortran 77 implementation of numerical recipes (amoeba)
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation version 2 of the License.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 */

#include <stdio.h>

/** Energy convergence cutoff for the SIMPLEX */
#define CONVERGENCE_SIMPLEX 0.05f
/** Maximum number of SIMPLEX iterations */
#define MAXITER_SIMPLEX 250

#define alpha 1.0
#define beta 0.5
#define gamma 2.0

double get_simplex_energy(MOL2 *mol, float *simplex, MOL2 *template, int **type1, int **type2, int pocket);
double get_simplex_energy_ingrid(MOL2 *mol, float *simplex, float ****grid, int *types, int **types2, float min_grid[3], float max_grid[3],int max_x, int max_y, int max_z, float spacing);


float gen_rand_float()
{       
        return rand() / ((float)(RAND_MAX)+1);
}       


/*
 *  Algorithm works that way:
 *
 *
 *  Step 1: Which point is the worst and the next one
 *  Step 2: Algorithm termination check
 *  Step 3: Average vector of simplex, extrapolate and compute energy
 *  Step 4: Is this point the best? --> Go 4a
 *          Is this point worse that the previous to the last one --> Go 4b
 *          Other cases --> Go 4c
 *       4a: Ok we have a good one. Lets try to find a better one. If we fail, use the good one.
 *       4b: We suck. If we totally suck use the better point. If we do not totally suck replace
 *           the sucker and try again with a contraction. If we fail, use the better point.
 *       4c: We are in the middle of something, just accept it and look for something better.
 *  Step 5: Loop
 */


/**
 *
 *	@brief Perform SIMPLEX minimization on already initialized parameters
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *	@param nvar Number of variables for the minimization (must be 6 for rigid block)
 *	@param simplex Pre initialized random rigid block transformations
 *	@param energies_orig Energies of the pre initialized transformations
 *	@param mymol MOL2 with the molecule
 *      @param grids Potential grids
 *      @param types Atom types for force field grids
 *      @param min_grid Minimum grid coordinates
 *      @param max_grid Maximum grid coordinates
 *      @param nx Closest grid point to atoms. X
 *      @param ny Closest grid point to atoms. Y
 *      @param nz Closest grid point to atoms. Z
 *      @param max_x Number of X grid points
 *      @param max_y Number of Y grid points
 *      @param max_z Number of Z grid points
 *      @param spacing Grid spacing
 *      @return 0 on success
 *
 */
int go_simplex_go(int nvar, float simplex[7][7], double **energies_orig,  MOL2 *mymol, MOL2 *template, int **types1, int **types2, int pocket)
{

	int iter = 0, mpts = 0;
	int ilo = 0, ihi = 0, inhi = 0;
	int i = 0, j = 0;
	float *pbar = NULL;
	float *simplexr = NULL;
	float *simplexrr = NULL;
	double energyr = 0.0f, energyrr = 0.0f;
	double deltaE = 0.0f;
	double *energies = NULL;

	energies = *energies_orig;
	mpts = nvar + 1;

	pbar = (float*)calloc(sizeof(float), nvar);
	simplexr = (float*)calloc(sizeof(float), nvar);
	simplexrr = (float*)calloc(sizeof(float), nvar);

	while ( 1 == 1) { /* (C) Crapy-lazzy techniques */

		/* Which point of the simplex is the worst (higher) and the best (lower) */
		if ( energies[0] > energies[1]) {
			ihi = 0;
			inhi = 1;
		}else{
			ihi = 1;
			inhi = 0;
		}

		for ( i = 0; i < mpts; ++i) {

			if ( energies[i] < energies[ilo]) ilo = i;

			if ( energies[i] > energies[ihi]) {
				inhi = ihi;
				ihi = i;
			}else if ( energies[i] > energies[inhi])
				if ( i != ihi) inhi = i;
		}

		deltaE = fabs( energies[ihi] - energies[ilo]);

/*#ifdef DEBUG*/
/*		fprintf(stderr, "GO_SIMPLEX_GO: deltaE: %f Iter: %i\r", deltaE, iter);*/
/*#endif*/



		/* Termination criteria */
		if ( deltaE < CONVERGENCE_SIMPLEX && iter != 0) return 0;
		if ( iter == MAXITER_SIMPLEX - 1) return 0;


		/* NEW ITERATION */

		iter++;

		/* Initialize */
		for (i = 0; i < nvar; ++i)
			pbar[j] = 0;


		/* Sum over vars of simplex to compute the average vector */
		for (i = 0; i < mpts; ++i) {
			if ( i != ihi) /* Exclude highest point */
				for (j = 0; j < nvar; ++j)
					pbar[j] = pbar[j] + simplex[i][j];
		}

		for ( i = 0; i < nvar; ++i) {
			pbar[i] = pbar[i] / (float)nvar;                                /* Average vector finished */
			simplexr[i] = (1 + alpha) * pbar[i] - alpha * simplex[ihi][i];  /* Extrapolate by factor alpha to construct the reflection */
		}

		/* Energy of this new point */
		energyr = get_simplex_energy(mymol, simplexr, template, types1, types2, pocket);

		if ( energyr <= energies[ilo]) { /* Energy better than the best point */

			/* Since we are here, lets try further extrapolation by factor gamma*/
			for (j = 0; j < nvar; ++j)
				simplexrr[j] = gamma * simplexr[j] + (1 - gamma) * pbar[j];

			/* Get this new energy */
			energyrr = get_simplex_energy(mymol, simplexrr, template, types1, types2, pocket);

			/* Oe oe oe oe! We are the best extrapolating. */
			/* This extrapolation is better again that the best point  */
			if ( energyrr < energies[ilo]) {
				for (j = 0; j < nvar; ++j)
					simplex[ihi][j] = simplexrr[j];
				energies[ihi] = energyrr;
			}else{  /* Ok. We just screw up the last move. Lets take the reflected point. */
				for (j = 0; j < nvar; ++j)
					simplex[ihi][j] = simplexr[j];
				energies[ihi] = energyr;
			}

		}else if ( energyr >= energies[inhi] ) { /* The reflected point sux, even more than the second higher */

			if ( energyr < energies[ihi] ) { /* That saved the day. At least we dont suck more than the highest */
				/* Replace the highest with this new fella' */
				for (j = 0; j < nvar; ++j)
					simplex[ihi][j] = simplexr[j];

				energies[ihi] = energyr;
			}

			/* Before giving up, lets try something in the middle with an one-dim contraction */
			for ( j = 0; j < nvar; ++j)
				simplexrr[j] = beta * simplex[ihi][j] + (1 - beta ) * pbar[j];

			energyrr = get_simplex_energy(mymol, simplexrr, template, types1, types2, pocket);

			if ( energyrr < energies[ihi] ) { /* Ok. We did it, that move worked.*/
				for (j = 0; j < nvar; ++j)
					simplex[ihi][j] = simplexrr[j];
				energies[ihi] = energyrr;
			}else{  /* Out of luck today. We suck with the high point. */

				/* Try better with the lowest shit */
				for ( i = 0; i < mpts; ++i) {
					if ( i != ilo ) {
						for ( j = 0; j < nvar; ++j) {
							simplexr[j] = 0.5 * ( simplex[i][j] + simplex[ilo][j] );
							simplex[i][j] = simplexr[j];
						}
						energies[i] = get_simplex_energy(mymol, simplexr, template, types1, types2, pocket); /* Should be better */
					}
				}

			}

		}else{  /* We are in the middle, just eat that and keep going */
			for ( j = 0; j < nvar; ++j)
				simplex[ihi][j] = simplexr[j];
			energies[ihi] = energyr;
		}

	}


	*energies_orig = energies;

	free(pbar);
	free(simplexr);
	free(simplexrr);

	return 0;
}

/**
 *
 *      @brief Apply a rigid block transformation and calculate the energy for the SIMPLEX min
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *      @param mol MOL2 with the molecule
 *      @param simplex Transformation (xyz + 3 euler angles)
 *      @param grids Potential grids
 *      @param types Atom types for force field grids
 *      @param min_grid Minimum grid coordinates
 *      @param max_grid Maximum grid coordinates
 *      @param nx Closest grid point to atoms. X
 *      @param ny Closest grid point to atoms. Y
 *      @param nz Closest grid point to atoms. Z
 *      @param max_x Number of X grid points
 *      @param max_y Number of Y grid points
 *      @param max_z Number of Z grid points
 *      @param spacing Grid spacing
 *      @return Grid energy
 */
double get_simplex_energy(MOL2 *mol, float *simplex, MOL2 *template, int **types1, int **types2, int pocket)
{

	int i = 0, j = 0, real_atoms = 0, l = 0;

	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;

	float *x_b = NULL, *y_b = NULL, *z_b = NULL;

	float ith = 0, iph = 0, ips = 0;
	float theta = 0.0f, phi = 0.0f, psi = 0.0f;
	float sthe = 0.0f, cthe = 0.0f, sphi = 0.0f, cphi = 0.0f, spsi = 0.0f, cpsi = 0.0f;
	float A[3][3], dx = 0.0f, dy = 0.0f, dz = 0.0f;
	const float cdr = 0.0174532777778;

	float comt[3], com[3], Xrot = 0.0f, Yrot = 0.0f, Zrot = 0.0f, VX = 0.0f, VY = 0.0f, VZ = 0.0f;
	int *tor_type = NULL;
	double energy = 0.0;
	int out_flag = 0;

	float aa = 0.0f, ab = 0.0, bb = 0.0;

	for ( i = 0; i < mol->n_atoms; ++i)
		if ( mol->backbone[i] != 1)
			++real_atoms;

	x_b = (float *) calloc(sizeof(float),mol->n_atoms);
        y_b = (float *) calloc(sizeof(float),mol->n_atoms);
        z_b = (float *) calloc(sizeof(float),mol->n_atoms);

        for ( i = 0; i < mol->n_atoms; ++i)
	{
		x_b[i] = mol->x[i];
                y_b[i] = mol->y[i];
                z_b[i] = mol->z[i];
	}

	x = simplex[0];
	y = simplex[1];
	z = simplex[2];

	ith = simplex[3];
	iph = simplex[4];
	ips = simplex[5];


	theta = cdr * ith;
	phi = cdr * iph;
	psi = cdr * ips;
	sthe = sin(theta);
	cthe = cos(theta);
	sphi = sin(phi);
	cphi = cos(phi);
	spsi = sin(psi);
	cpsi = cos(psi);

	A[0][0] = -sphi * spsi + cthe * cphi * cpsi;
	A[0][1] =  cphi * spsi + cthe * sphi * cpsi;
	A[0][2] = -sthe * cpsi;
	A[1][0] = -sphi * cpsi - cthe * cphi * spsi;
	A[1][1] =  cphi * cpsi - cthe * sphi * spsi;
	A[1][2] =  sthe * spsi;
	A[2][0] =  sthe * cphi;
	A[2][1] =  sthe * sphi;
	A[2][2] =  cthe;


	for (i = 0; i < 3; ++i)
		com[i] = 0.0f;

	for (i = 0; i < real_atoms; ++i) {
		com[0] = com[0] + mol->x[i];
		com[1] = com[1] + mol->y[i];
		com[2] = com[2] + mol->z[i];
	}

	for (i = 0; i < 3; ++i)
		com[i] = com[i] / (float)real_atoms;


	comt[0] = com[0] + x;
	comt[1] = com[1] + y;
	comt[2] = com[2] + z;

	dx = com[0] - comt[0];
	dy = com[1] - comt[1];
	dz = com[2] - comt[2];

	for ( i = 0; i < real_atoms; ++i) {
		mol->x[i] = mol->x[i] - dx;
		mol->y[i] = mol->y[i] - dy;
		mol->z[i] = mol->z[i] - dz;
	}

	out_flag = 0;

	for ( i = 0; i < real_atoms; ++i) {

		VX = mol->x[i] - comt[0];
		VY = mol->y[i] - comt[1];
		VZ = mol->z[i] - comt[2];

		Xrot = (A[0][0] * VX + A[1][0] * VY + A[2][0] * VZ) + comt[0];
		Yrot = (A[0][1] * VX + A[1][1] * VY + A[2][1] * VZ) + comt[1];
		Zrot = (A[0][2] * VX + A[1][2] * VY + A[2][2] * VZ) + comt[2];


		mol->x[i] = Xrot;
		mol->y[i] = Yrot;
		mol->z[i] = Zrot;
	}


/*	energy = 0.0f;
	aa = get_volumen_intersection(template, &template);
        bb = get_volumen_intersection(mol, &mol);
	ab = get_volumen_intersection(template, &mol);

        energy = -(ab / (aa+bb-ab));

        aa = get_color_intersection(template, &template);
        bb = get_color_intersection(mol, &mol);
        ab = get_color_intersection(template, &mol);

        energy += -(ab / (aa+bb-ab));*/
	energy = get_volumen_intersection(template, &mol, pocket, 0);
	energy += get_color_intersection(template, &mol, types1, types2, pocket, 0);

/*	fprintf(stderr,"Evaluating energy: %f\n",energy);*/

        for ( i = 0; i < mol->n_atoms; ++i)
        {
                mol->x[i] = x_b[i];
                mol->y[i] = y_b[i];
                mol->z[i] = z_b[i];
        }

	free(x_b); free(y_b); free(z_b);

	return energy;

}


/**
 *
 *	@brief Perform SIMPLEX minimization on already initialized parameters
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *	@param nvar Number of variables for the minimization (must be 6 for rigid block)
 *	@param simplex Pre initialized random rigid block transformations
 *	@param energies_orig Energies of the pre initialized transformations
 *	@param mymol MOL2 with the molecule
 *      @param grids Potential grids
 *      @param types Atom types for force field grids
 *      @param min_grid Minimum grid coordinates
 *      @param max_grid Maximum grid coordinates
 *      @param nx Closest grid point to atoms. X
 *      @param ny Closest grid point to atoms. Y
 *      @param nz Closest grid point to atoms. Z
 *      @param max_x Number of X grid points
 *      @param max_y Number of Y grid points
 *      @param max_z Number of Z grid points
 *      @param spacing Grid spacing
 *      @return 0 on success
 *
 */
int go_simplex_go_ingrid(int nvar, float simplex[7][7], double **energies_orig,  MOL2 *mymol, float ****grids, int *types, int **types2,float min_grid[3], float max_grid[3],int max_x, int max_y, int max_z, float spacing)
{

	int iter = 0, mpts = 0;
	int ilo = 0, ihi = 0, inhi = 0;
	int i = 0, j = 0;
	float *pbar = NULL;
	float *simplexr = NULL;
	float *simplexrr = NULL;
	double energyr = 0.0f, energyrr = 0.0f;
	double deltaE = 0.0f;
	double *energies = NULL;

	energies = *energies_orig;
	mpts = nvar + 1;

	pbar = (float*)calloc(sizeof(float), mpts);
	simplexr = (float*)calloc(sizeof(float), mpts);
	simplexrr = (float*)calloc(sizeof(float), mpts);

	while ( 1 == 1) { /* (C) Crapy-lazzy techniques */

		/* Which point of the simplex is the worst (higher) and the best (lower) */
		if ( energies[0] > energies[1]) {
			ihi = 0;
			inhi = 1;
		}else{
			ihi = 1;
			inhi = 0;
		}

		for ( i = 0; i < mpts; ++i) {

			if ( energies[i] < energies[ilo]) ilo = i;

			if ( energies[i] > energies[ihi]) {
				inhi = ihi;
				ihi = i;
			}else if ( energies[i] > energies[inhi])
				if ( i != ihi) inhi = i;
		}

		deltaE = fabs( energies[ihi] - energies[ilo]);

/*#ifdef DEBUG*/
/*		fprintf(stderr, "GO_SIMPLEX_GO: deltaE: %f Iter: %i\r", deltaE, iter);*/
/*#endif*/



		/* Termination criteria */
		if ( deltaE < CONVERGENCE_SIMPLEX && iter != 0)
		{
			free(simplexr);
			free(simplexrr);
			free(pbar);
			return 0;
		}

		if ( iter == MAXITER_SIMPLEX - 1)
		{
			free(simplexr);
                        free(simplexrr);
                        free(pbar);
			return 0;
		}


		/* NEW ITERATION */

		iter++;

		/* Initialize */
		for (i = 0; i < nvar; ++i)
			pbar[i] = 0;


		/* Sum over vars of simplex to compute the average vector */
		for (i = 0; i < mpts; ++i) {
			if ( i != ihi) /* Exclude highest point */
				for (j = 0; j < nvar; ++j)
					pbar[j] = pbar[j] + simplex[i][j];
		}

		for ( i = 0; i < nvar; ++i) {
			/*pbar[i] = pbar[i] / (float)nvar;*/                              
			pbar[i] = pbar[i] / (float) mpts;                                /* Average vector finished */
			simplexr[i] = (1 + alpha) * pbar[i] - alpha * simplex[ihi][i];  /* Extrapolate by factor alpha to construct the reflection */
		}

		/* Energy of this new point */
		energyr = get_simplex_energy_ingrid(mymol, simplexr, grids, types, types2, min_grid, max_grid, max_x, max_y, max_z, spacing);

		if ( energyr <= energies[ilo]) { /* Energy better than the best point */

			/* Since we are here, lets try further extrapolation by factor gamma*/
			for (j = 0; j < nvar; ++j)
				simplexrr[j] = gamma * simplexr[j] + (1 - gamma) * pbar[j];

			/* Get this new energy */
			energyrr = get_simplex_energy_ingrid(mymol, simplexrr, grids, types, types2, min_grid, max_grid, max_x, max_y, max_z, spacing);

			/* Oe oe oe oe! We are the best extrapolating. */
			/* This extrapolation is better again that the best point  */
			if ( energyrr < energies[ilo]) {
				for (j = 0; j < nvar; ++j)
					simplex[ihi][j] = simplexrr[j];
				energies[ihi] = energyrr;
			}else{  /* Ok. We just screw up the last move. Lets take the reflected point. */
				for (j = 0; j < nvar; ++j)
					simplex[ihi][j] = simplexr[j];
				energies[ihi] = energyr;
			}

		}else if ( energyr >= energies[inhi] ) { /* The reflected point sux, even more than the second higher */

			if ( energyr < energies[ihi] ) { /* That saved the day. At least we dont suck more than the highest */
				/* Replace the highest with this new fella' */
				for (j = 0; j < nvar; ++j)
					simplex[ihi][j] = simplexr[j];

				energies[ihi] = energyr;
			}

			/* Before giving up, lets try something in the middle with an one-dim contraction */
			for ( j = 0; j < nvar; ++j)
				simplexrr[j] = beta * simplex[ihi][j] + (1 - beta ) * pbar[j];

			energyrr = get_simplex_energy_ingrid(mymol, simplexrr, grids, types, types2, min_grid, max_grid, max_x, max_y, max_z, spacing);

			if ( energyrr < energies[ihi] ) { /* Ok. We did it, that move worked.*/
				for (j = 0; j < nvar; ++j)
					simplex[ihi][j] = simplexrr[j];
				energies[ihi] = energyrr;
			}else{  /* Out of luck today. We suck with the high point. */

				/* Try better with the lowest shit */
				for ( i = 0; i < mpts; ++i) {
					if ( i != ilo ) {
						for ( j = 0; j < nvar; ++j) {
							simplexr[j] = 0.5 * ( simplex[i][j] + simplex[ilo][j] );
							simplex[i][j] = simplexr[j];
						}
						energies[i] = get_simplex_energy_ingrid(mymol, simplexr, grids, types, types2, min_grid, max_grid, max_x, max_y, max_z, spacing); /* Should be better */
					}
				}

			}

		}else{  /* We are in the middle, just eat that and keep going */
			for ( j = 0; j < nvar; ++j)
				simplex[ihi][j] = simplexr[j];
			energies[ihi] = energyr;
		}

	}


	*energies_orig = energies;

	free(pbar);
	free(simplexr);
	free(simplexrr);

	return 0;
}



/**
 *
 *      @brief Apply a rigid block transformation and calculate the energy for the SIMPLEX min
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *      @param mol MOL2 with the molecule
 *      @param simplex Transformation (xyz + 3 euler angles)
 *      @param grids Potential grids
 *      @param types Atom types for force field grids
 *      @param min_grid Minimum grid coordinates
 *      @param max_grid Maximum grid coordinates
 *      @param nx Closest grid point to atoms. X
 *      @param ny Closest grid point to atoms. Y
 *      @param nz Closest grid point to atoms. Z
 *      @param max_x Number of X grid points
 *      @param max_y Number of Y grid points
 *      @param max_z Number of Z grid points
 *      @param spacing Grid spacing
 *      @return Grid energy
 */
double get_simplex_energy_ingrid(MOL2 *mol, float *simplex, float ****grids, int *types, int **types2, float min_grid[3], float max_grid[3],int max_x, int max_y, int max_z, float spacing)
{

	int i = 0, j = 0, real_atoms = 0, l = 0;

	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;

	float *x_b = NULL, *y_b = NULL, *z_b = NULL;

	float ith = 0, iph = 0, ips = 0;
	float theta = 0.0f, phi = 0.0f, psi = 0.0f;
	float sthe = 0.0f, cthe = 0.0f, sphi = 0.0f, cphi = 0.0f, spsi = 0.0f, cpsi = 0.0f;
	float A[3][3], dx = 0.0f, dy = 0.0f, dz = 0.0f;
	const float cdr = 0.0174532777778;

	float comt[3], com[3], Xrot = 0.0f, Yrot = 0.0f, Zrot = 0.0f, VX = 0.0f, VY = 0.0f, VZ = 0.0f;
	int *tor_type = NULL;
	double energy = 0.0;
	int out_flag = 0;

	float aa = 0.0f, ab = 0.0, bb = 0.0;

/*	for ( i = 0; i < mol->n_atoms; ++i)
		if ( mol->backbone[i] != 1)
			++real_atoms;*/

	real_atoms = mol->n_atoms;
	x_b = (float *) calloc(sizeof(float),mol->n_atoms);
        y_b = (float *) calloc(sizeof(float),mol->n_atoms);
        z_b = (float *) calloc(sizeof(float),mol->n_atoms);

        for ( i = 0; i < mol->n_atoms; ++i)
	{
		x_b[i] = mol->x[i];
                y_b[i] = mol->y[i];
                z_b[i] = mol->z[i];
	}

	x = simplex[0];
	y = simplex[1];
	z = simplex[2];

	ith = simplex[3];
	iph = simplex[4];
	ips = simplex[5];


	theta = cdr * ith;
	phi = cdr * iph;
	psi = cdr * ips;
	sthe = sin(theta);
	cthe = cos(theta);
	sphi = sin(phi);
	cphi = cos(phi);
	spsi = sin(psi);
	cpsi = cos(psi);

	A[0][0] = -sphi * spsi + cthe * cphi * cpsi;
	A[0][1] =  cphi * spsi + cthe * sphi * cpsi;
	A[0][2] = -sthe * cpsi;
	A[1][0] = -sphi * cpsi - cthe * cphi * spsi;
	A[1][1] =  cphi * cpsi - cthe * sphi * spsi;
	A[1][2] =  sthe * spsi;
	A[2][0] =  sthe * cphi;
	A[2][1] =  sthe * sphi;
	A[2][2] =  cthe;


	for (i = 0; i < 3; ++i)
		com[i] = 0.0f;

	for (i = 0; i < real_atoms; ++i) {
		com[0] = com[0] + mol->x[i];
		com[1] = com[1] + mol->y[i];
		com[2] = com[2] + mol->z[i];
	}

	for (i = 0; i < 3; ++i)
		com[i] = com[i] / (float)real_atoms;


	comt[0] = com[0] + x;
	comt[1] = com[1] + y;
	comt[2] = com[2] + z;

	dx = com[0] - comt[0];
	dy = com[1] - comt[1];
	dz = com[2] - comt[2];

	for ( i = 0; i < real_atoms; ++i) {
		mol->x[i] = mol->x[i] - dx;
		mol->y[i] = mol->y[i] - dy;
		mol->z[i] = mol->z[i] - dz;
	}

	out_flag = 0;

	for ( i = 0; i < real_atoms; ++i) {

		VX = mol->x[i] - comt[0];
		VY = mol->y[i] - comt[1];
		VZ = mol->z[i] - comt[2];

		Xrot = (A[0][0] * VX + A[1][0] * VY + A[2][0] * VZ) + comt[0];
		Yrot = (A[0][1] * VX + A[1][1] * VY + A[2][1] * VZ) + comt[1];
		Zrot = (A[0][2] * VX + A[1][2] * VY + A[2][2] * VZ) + comt[2];


		mol->x[i] = Xrot;
		mol->y[i] = Yrot;
		mol->z[i] = Zrot;
	}

/*	energy = get_volumen_intersection(template, &mol);
	energy += get_color_intersection(template, &mol);*/

	energy = get_volume_intersection_fromgrid(mol, grids, types, min_grid, max_grid, max_x, max_y, max_z, spacing);
	energy += get_color_volume_intersection_fromgrid(mol, grids, types2, min_grid, max_grid, max_x, max_y, max_z, spacing);


/*	fprintf(stderr,"Evaluating energy: %f\n",energy);*/

        for ( i = 0; i < mol->n_atoms; ++i)
        {
                mol->x[i] = x_b[i];
                mol->y[i] = y_b[i];
                mol->z[i] = z_b[i];
        }

	free(x_b); free(y_b); free(z_b);

	return energy;

}


/* END OF MODULE */
