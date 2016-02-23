/**
 *
 *      @file force.c
 *      @brief Calculate first derivative by numericall method
 *	(you know, the crappy stuff of (f(x+h)-f(x-h))/2h
 *
 *      Reason to do this: analytical differentiation with
 *      torsional angles is hell!!
 *
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @date 01/10/2011
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
#include <stdlib.h>

#define DEV_STEP 0.1

/**
 *
 *	@brief Calculate off-grid GAFF numerical derivatives
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mymol MOL2 with the molecule
 *	@return 0 on success
 *
 */
int get_numerical_forces(MOL2 **mymol)
{
	MOL2 *mol = NULL;
	float *xc = NULL, *yc = NULL, *zc = NULL;
	int i = 0;
	double init_tot = 0;
	double move_tot = 0;
	int *tor_type = NULL;


	mol = *mymol;


	if ( (xc = (float*)calloc(sizeof(float), mol->n_atoms)) == NULL)
		return -2;

	if ( (yc = (float*)calloc(sizeof(float), mol->n_atoms)) == NULL) {
		free(xc);
		return -2;
	}

	if ( (zc = (float*)calloc(sizeof(float), mol->n_atoms)) == NULL) {
		free(xc);
		free(yc);
		return -2;
	}

	tor_type = (int*)calloc(sizeof(int), mol->n_bonds * 6);


	for (i = 0; i < mol->n_atoms; ++i) {

		xc[i] = mol->x[i];
		mol->x[i] = mol->x[i] + DEV_STEP;
		get_energy(&mol, &move_tot, &tor_type, 1);
		mol->x[i] = xc[i] - DEV_STEP;
		get_energy(&mol, &init_tot, &tor_type, 1);
		mol->x[i] = xc[i];
		mol->grads_X[i] = ((move_tot - init_tot) / (2 * DEV_STEP));


		yc[i] = mol->y[i];
		mol->y[i] = mol->y[i] + DEV_STEP;
		get_energy(&mol, &move_tot, &tor_type, 1);
		mol->y[i] = yc[i] - DEV_STEP;
		get_energy(&mol, &init_tot, &tor_type, 1);
		mol->y[i] = yc[i];
		mol->grads_Y[i] = ((move_tot - init_tot) / (2 * DEV_STEP));


		zc[i] = mol->z[i];
		mol->z[i] = mol->z[i] + DEV_STEP;
		get_energy(&mol, &move_tot, &tor_type, 1);
		mol->z[i] = zc[i] - DEV_STEP;
		get_energy(&mol, &init_tot, &tor_type, 1);
		mol->z[i] = zc[i];
		mol->grads_Z[i] = ((move_tot - init_tot) / (2 * DEV_STEP));


	}

	free(xc);
	free(yc);
	free(zc);
	free(tor_type);

	*mymol = mol;
	return 0;
}

/**
 *
 *	@brief Numerical derivatives from GAFF grids with finite displacements
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mymol MOL2 with the ligand
 *	@param grids Potential grids
 *	@param types Atom types for the potential grids
 *	@param min_grid Minimum grid coorindates
 *	@param max_grid Maximum grid coordinates
 *	@param nx Closest grid point to atoms. Axis x
 *	@param ny Closest grid point to atoms. Axis y
 *	@param nz Closest grid point to atoms. Axis z
 *	@param max_x Maximum grid points. Axis x
 *	@param max_y Maximum grid points. Axis y
 *	@param max_z Maximum grid points. Axis z
 *	@param spacing Grid spacing
 *	@return 0 on success
 *
 */
int get_numerical_in_grid_forces(MOL2 **mymol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing)
{
	MOL2 *mol = NULL;
	float *xc = NULL, *yc = NULL, *zc = NULL;
	int i = 0;
	double init_tot = 0;
	double move_tot = 0;
	int *tor_type = NULL;


	mol = *mymol;


	if ( (xc = (float*)calloc(sizeof(float), mol->n_atoms)) == NULL)
		return -2;

	if ( (yc = (float*)calloc(sizeof(float), mol->n_atoms)) == NULL) {
		free(xc);
		return -2;
	}

	if ( (zc = (float*)calloc(sizeof(float), mol->n_atoms)) == NULL) {
		free(xc);
		free(yc);
		return -2;
	}

	tor_type = (int*)calloc(sizeof(int), mol->n_bonds * 6);


	for (i = 0; i < mol->n_atoms; ++i) {

		xc[i] = mol->x[i];
		mol->x[i] = mol->x[i] + DEV_STEP;
		get_energy(&mol, &move_tot, &tor_type, 2);
		move_tot += get_ff_grid_energy(mol, grids, types, min_grid, max_grid, nx, ny, nz, max_x, max_y, max_z, spacing);
		mol->x[i] = xc[i] - DEV_STEP;
		get_energy(&mol, &init_tot, &tor_type, 2);
		move_tot += get_ff_grid_energy(mol, grids, types, min_grid, max_grid, nx, ny, nz, max_x, max_y, max_z, spacing);
		mol->x[i] = xc[i];
		mol->grads_X[i] = ((move_tot - init_tot) / (2 * DEV_STEP));


		yc[i] = mol->y[i];
		mol->y[i] = mol->y[i] + DEV_STEP;
		get_energy(&mol, &move_tot, &tor_type, 2);
		move_tot += get_ff_grid_energy(mol, grids, types, min_grid, max_grid, nx, ny, nz, max_x, max_y, max_z, spacing);
		mol->y[i] = yc[i] - DEV_STEP;
		get_energy(&mol, &init_tot, &tor_type, 2);
		move_tot += get_ff_grid_energy(mol, grids, types, min_grid, max_grid, nx, ny, nz, max_x, max_y, max_z, spacing);
		mol->y[i] = yc[i];
		mol->grads_Y[i] = ((move_tot - init_tot) / (2 * DEV_STEP));


		zc[i] = mol->z[i];
		mol->z[i] = mol->z[i] + DEV_STEP;
		get_energy(&mol, &move_tot, &tor_type, 2);
		move_tot += get_ff_grid_energy(mol, grids, types, min_grid, max_grid, nx, ny, nz, max_x, max_y, max_z, spacing);
		mol->z[i] = zc[i] - DEV_STEP;
		get_energy(&mol, &init_tot, &tor_type, 2);
		move_tot += get_ff_grid_energy(mol, grids, types, min_grid, max_grid, nx, ny, nz, max_x, max_y, max_z, spacing);
		mol->z[i] = zc[i];
		mol->grads_Z[i] = ((move_tot - init_tot) / (2 * DEV_STEP));

	}

	free(xc);
	free(yc);
	free(zc);
	free(tor_type);

	*mymol = mol;
	return 0;
}


/**
 *
 *      @brief Interpolate derivatives from GAFF grids. Method 2
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @param mymol MOL2 with the ligand
 *      @param grids Potential grids
 *      @param types Atom types for the potential grids
 *      @param min_grid Minimum grid coorindates
 *      @param max_grid Maximum grid coordinates
 *      @param nx Closest grid point to atoms. Axis x
 *      @param ny Closest grid point to atoms. Axis y
 *      @param nz Closest grid point to atoms. Axis z
 *      @param max_x Maximum grid points. Axis x
 *      @param max_y Maximum grid points. Axis y
 *      @param max_z Maximum grid points. Axis z
 *      @param spacing Grid spacing
 *      @return 0 on success
 *
 */
void add_ff_grid_forces_finite(MOL2 **mymol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing)
{

	int i = 0, real_atoms = 0, out_flag = 0, l = 0;
	float energy = 0.0f;
	float energy2 = 0.0f;
	MOL2 *mol = NULL;

	mol = *mymol;

	for ( i = 0; i < mol->n_atoms; ++i)
		++real_atoms;


	for ( i = 0; i < real_atoms; ++i) {
		if ( mol->x[i] < min_grid[0] || mol->x[i] > max_grid[0]) {
			out_flag = 1;
			break;
		}

		if ( mol->y[i] < min_grid[1] || mol->y[i] > max_grid[1]) {
			out_flag = 1;
			break;
		}

		if ( mol->z[i] < min_grid[2] || mol->z[i] > max_grid[2]) {
			out_flag = 1;
			break;
		}
	}


	if ( out_flag == 0) {


		for ( i = 0; l < real_atoms; ++l) {
			nx[l] = (int)((mol->x[l] - min_grid[0]) / spacing);
			ny[l] = (int)((mol->y[l] - min_grid[1]) / spacing);
			nz[l] = (int)((mol->z[l] - min_grid[2]) / spacing);
		}

		for ( l = 0; l < real_atoms; ++l) {
			energy = 0.0f;
			energy2 = 0.0f;


			if ( types[l] != -1 && nx[l] >= 1 && ny[l] >= 1 &&
			     nz[l] >= 1 && nx[l] < max_x - 1 && ny[l] < max_y - 1 &&
			     nz[l] < max_z - 1) {

				energy =
					grids[types[l]][nx[l] + 1][ny[l]][nz[l]] +
					mol->pcharges[l] *
					grids[8][nx[l] + 1][ny[l]][nz[l]];

				energy2 =
					grids[types[l]][nx[l] - 1][ny[l]][nz[l]] +
					mol->pcharges[l] *
					grids[8][nx[l] - 1][ny[l]][nz[l]];
				mol->grads_X[l] += (energy - energy2) / (2.0f * spacing);

				energy =
					grids[types[l]][nx[l]][ny[l] + 1][nz[l]] +
					mol->pcharges[l] *
					grids[8][nx[l] + 1][ny[l]][nz[l]];

				energy2 =
					grids[types[l]][nx[l]][ny[l] - 1][nz[l]] +
					mol->pcharges[l] *
					grids[8][nx[l]][ny[l] - 1][nz[l]];
				mol->grads_Y[l] += (energy - energy2) / (2.0f * spacing);

				energy =
					grids[types[l]][nx[l]][ny[l]][nz[l] + 1] +
					mol->pcharges[l] *
					grids[8][nx[l]][ny[l]][nz[l] - 1];

				energy2 =
					grids[types[l]][nx[l]][ny[l]][nz[l] - 1] +
					mol->pcharges[l] *
					grids[8][nx[l]][ny[l]][nz[l] + 1];
				mol->grads_Z[l] += (energy - energy2) / (2.0f * spacing);



			}
		}


	}else
		energy = 9999.9;


	return;
}

/**
 *
 *      @brief Interpolate derivatives from GAFF grids. Method 1
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @param mymol MOL2 with the ligand
 *      @param grids Potential grids
 *      @param types Atom types for the potential grids
 *      @param min_grid Minimum grid coorindates
 *      @param max_grid Maximum grid coordinates
 *      @param nx Closest grid point to atoms. Axis x
 *      @param ny Closest grid point to atoms. Axis y
 *      @param nz Closest grid point to atoms. Axis z
 *      @param max_x Maximum grid points. Axis x
 *      @param max_y Maximum grid points. Axis y
 *      @param max_z Maximum grid points. Axis z
 *      @param spacing Grid spacing
 *      @return 0 on success
 *
 */

void add_ff_grid_forces(MOL2 **mymol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing)
{

	int i = 0, real_atoms = 0, out_flag = 0, l = 0, kpot = 0, ix, iy, iz;
	float energy = 0.0f;
	float energy2 = 0.0f, xa, ya, za, xd, yd, zd, q, x, y, z;
	float V000, V100, V010, V110, V001, V101, V011, V111;

	MOL2 *mol = NULL;

	mol = *mymol;

	for ( i = 0; i < mol->n_atoms; ++i)
		++real_atoms;


	for ( i = 0; i < real_atoms; ++i) {
		if ( mol->x[i] < min_grid[0] || mol->x[i] > max_grid[0]) {
			out_flag = 1;
			break;
		}

		if ( mol->y[i] < min_grid[1] || mol->y[i] > max_grid[1]) {
			out_flag = 1;
			break;
		}

		if ( mol->z[i] < min_grid[2] || mol->z[i] > max_grid[2]) {
			out_flag = 1;
			break;
		}
	}


	if ( out_flag == 0) {


		for ( i = 0; i < real_atoms; ++i) {
			energy = 0.0f;
			energy2 = 0.0f;


			kpot = types[i];
			x = mol->x[i];
			y = mol->y[i];
			z = mol->z[i];

			nx[i] = (int)((mol->x[i] - min_grid[0]) / spacing);
			ny[i] = (int)((mol->y[i] - min_grid[1]) / spacing);
			nz[i] = (int)((mol->z[i] - min_grid[2]) / spacing);

			ix = nx[i];
			iy = ny[i];
			iz = nz[i];
			q = mol->pcharges[i];

			if ( types[i] != -1 && ix >= 1 && iy >= 1 &&
			     iz >= 1 && ix < max_x - 1 && iy < max_y - 1 &&
			     iz < max_z - 1 ) {

				V000 = grids[kpot][ix][iy][iz];
				V100 = grids[kpot][ix + 1][iy][iz];
				V010 = grids[kpot][ix][iy + 1][iz];
				V110 = grids[kpot][ix + 1][iy + 1][iz];
				V001 = grids[kpot][ix][iy][iz + 1];
				V101 = grids[kpot][ix + 1][iy][iz + 1];
				V011 = grids[kpot][ix][iy + 1][iz + 1];
				V111 = grids[kpot][ix + 1][iy + 1][iz + 1];

				xa = (ix * spacing) + min_grid[0];
				ya = (iy * spacing) + min_grid[1];
				za = (iz * spacing) + min_grid[2];

				xd = x - xa;
				yd = y - ya;
				zd = z - za;


				/* dVinterpoled/dx */

				mol->grads_X[i] += V000 * -1 * (1 - yd) * (1 - zd) +
						   V100 * (1 - yd) * (1 - zd) +
						   V010 * -1 * yd * (1 - zd) +
						   V110 * yd * (1 - zd) +
						   V001 * -1 * (1 - yd) * zd +
						   V101 * (1 - yd) * zd +
						   V011 * -1 * yd * zd +
						   V111 * yd * zd;

				mol->grads_Y[i] += V000 * (1 - xd) * -1 * (1 - zd) +
						   V100 * xd * -1 * (1 - zd) +
						   V010 * (1 - xd) *  (1 - zd) +
						   V110 * xd * (1 - zd) +
						   V001 * (1 - xd) * -1 * zd +
						   V101 * xd * -1 * zd +
						   V011 * (1 - xd) * zd +
						   V111 * xd * zd;

				mol->grads_Z[i] += V000 * (1 - xd) * (1 - yd) * -1 +
						   V100 * xd * (1 - yd) * -1 +
						   V010 * (1 - xd) * yd * -1 +
						   V110 * xd * yd * -1 +
						   V001 * (1 - xd) * (1 - yd) +
						   V101 * xd * (1 - yd) +
						   V011 * (1 - xd) * yd +
						   V111 * xd * yd;

				kpot = 8;
				V000 = grids[kpot][ix][iy][iz];
				V100 = grids[kpot][ix + 1][iy][iz];
				V010 = grids[kpot][ix][iy + 1][iz];
				V110 = grids[kpot][ix + 1][iy + 1][iz];
				V001 = grids[kpot][ix][iy][iz + 1];
				V101 = grids[kpot][ix + 1][iy][iz + 1];
				V011 = grids[kpot][ix][iy + 1][iz + 1];
				V111 = grids[kpot][ix + 1][iy + 1][iz + 1];

				mol->grads_X[i] += V000 * -1 * (1 - yd) * (1 - zd) +
						   V100 * (1 - yd) * (1 - zd) +
						   V010 * -1 * yd * (1 - zd) +
						   V110 * yd * (1 - zd) +
						   V001 * -1 * (1 - yd) * zd +
						   V101 * (1 - yd) * zd +
						   V011 * -1 * yd * zd +
						   V111 * yd * zd;

				mol->grads_Y[i] += V000 * (1 - xd) * -1 * (1 - zd) +
						   V100 * xd * -1 * (1 - zd) +
						   V010 * (1 - xd) *  (1 - zd) +
						   V110 * xd * (1 - zd) +
						   V001 * (1 - xd) * -1 * zd +
						   V101 * xd * -1 * zd +
						   V011 * (1 - xd) * zd +
						   V111 * xd * zd;

				mol->grads_Z[i] += V000 * (1 - xd) * (1 - yd) * -1 +
						   V100 * xd * (1 - yd) * -1 +
						   V010 * (1 - xd) * yd * -1 +
						   V110 * xd * yd * -1 +
						   V001 * (1 - xd) * (1 - yd) +
						   V101 * xd * (1 - yd) +
						   V011 * (1 - xd) * yd +
						   V111 * xd * yd;



			}

		}

	}else
		energy = 9999.9;


	return;
}

/**
 *
 *      @brief Interpolate derivatives from GAFF grids. Method 3
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @param mymol MOL2 with the ligand
 *      @param grids Potential grids
 *      @param types Atom types for the potential grids
 *      @param min_grid Minimum grid coorindates
 *      @param max_grid Maximum grid coordinates
 *      @param nx Closest grid point to atoms. Axis x
 *      @param ny Closest grid point to atoms. Axis y
 *      @param nz Closest grid point to atoms. Axis z
 *      @param max_x Maximum grid points. Axis x
 *      @param max_y Maximum grid points. Axis y
 *      @param max_z Maximum grid points. Axis z
 *      @param spacing Grid spacing
 *      @return 0 on success
 *
 */

void add_ff_grid_forces_test(MOL2 **mymol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing)
{

	int i = 0, real_atoms = 0, out_flag = 0, l = 0, kpot = 0, ix, iy, iz;
	float energy = 0.0f;
	float energy2 = 0.0f, xa, ya, za, xd, yd, zd, q, x, y, z;
	float V000, V100, V010, V110, V001, V101, V011, V111, e1, e2, e3;
	register double u,   v,   w;
	register double p0u, p0v, p0w;
	register double p1u, p1v, p1w;
	register int u0,  v0,  w0;
	register int u1,  v1,  w1;

	MOL2 *mol = NULL;

	mol = *mymol;

	for ( i = 0; i < mol->n_atoms; ++i)
		++real_atoms;


	for ( i = 0; i < real_atoms; ++i) {
		if ( mol->x[i] < min_grid[0] || mol->x[i] > max_grid[0]) {
			out_flag = 1;
			break;
		}

		if ( mol->y[i] < min_grid[1] || mol->y[i] > max_grid[1]) {
			out_flag = 1;
			break;
		}

		if ( mol->z[i] < min_grid[2] || mol->z[i] > max_grid[2]) {
			out_flag = 1;
			break;
		}
	}


	if ( out_flag == 0) {


		for ( i = 0; i < real_atoms; ++i) {
			energy = 0.0f;
			energy2 = 0.0f;


			kpot = types[i];
			x = mol->x[i];
			y = mol->y[i];
			z = mol->z[i];

			nx[i] = (int)((mol->x[i] - min_grid[0]) / spacing);
			ny[i] = (int)((mol->y[i] - min_grid[1]) / spacing);
			nz[i] = (int)((mol->z[i] - min_grid[2]) / spacing);

			ix = nx[i];
			iy = ny[i];
			iz = nz[i];
			q = mol->pcharges[i];

			if ( types[i] != -1 && ix >= 1 && iy >= 1 &&
			     iz >= 1 && ix < max_x - 1 && iy < max_y - 1 &&
			     iz < max_z - 1 ) {

				xd = min_grid[0];
				yd = min_grid[1];
				zd = min_grid[2];

				w1  = (w0 = (int)(w = (x - xd) / spacing )) + 1;
				p1w = 1.0L - (p0w = w - (float)w0);

				v1  = (v0 = (int)(v = (y - yd) / spacing)) + 1;
				p1v = 1.0L - (p0v = v - (float)v0);

				u1  = (u0 = (int)(u = (z - zd) / spacing)) + 1;
				p1u = 1.0L - (p0u = u - (float)u0);


				e1 = p1u * p1v * grids[kpot][ w0 ][ v0 ][ u0 ];
				e1 += p0u * p1v *  grids[kpot][ w0 ][ v0 ][ u1 ];
				e1 += p1u * p0v *  grids[kpot][ w0 ][ v1 ][ u0 ];
				e1 += p0u * p0v *  grids[kpot][ w0 ][ v1 ][ u1 ];

				mol->grads_X[i] = e1;

				e1 = p1u * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
				e1 += p0u * p1w * grids[kpot][ w0 ][ v0 ][ u1 ];
				e1 += p1u * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
				e1 += p0u * p0w * grids[kpot][ w1 ][ v0 ][ u1 ];

				mol->grads_Y[i] = e1;

				e1 = p1v * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
				e1 += p0v * p1w * grids[kpot][ w0 ][ v1 ][ u0 ];
				e1 += p1v * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
				e1 += p0v * p0w * grids[kpot][ w1 ][ v1 ][ u0 ];

				mol->grads_Z[i] = e1;

				e1 = p1u * p1v * grids[kpot][ w0 ][ v0 ][ u0 ];
				e1 += p0u * p1v *  grids[kpot][ w0 ][ v0 ][ u1 ];
				e1 += p1u * p0v *  grids[kpot][ w0 ][ v1 ][ u0 ];
				e1 += p0u * p0v *  grids[kpot][ w0 ][ v1 ][ u1 ];

				mol->grads_X[i] = e1;

				e1 = p1u * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
				e1 += p0u * p1w * grids[kpot][ w0 ][ v0 ][ u1 ];
				e1 += p1u * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
				e1 += p0u * p0w * grids[kpot][ w1 ][ v0 ][ u1 ];

				mol->grads_Y[i] = e1;

				e1 = p1v * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
				e1 += p0v * p1w * grids[kpot][ w0 ][ v1 ][ u0 ];
				e1 += p1v * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
				e1 += p0v * p0w * grids[kpot][ w1 ][ v1 ][ u0 ];

				mol->grads_Z[i] = e1;




			}

		}

	}else
		energy = 9999.9;


	return;
}
