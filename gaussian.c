/**
 *
 *	@file gaussian.c
 *	@brief Gaussian volume overlapping module with rigid body Broyden-Fletcher-Goldfarb-Shanno (BFGS) minimizer
 *
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@date 05/03/2013
 *
 * 	This program is free software; you can redistribute it and/or modify
 * 	it under the terms of the GNU General Public License as published by
 * 	the Free Software Foundation version 3 of the License.
 *
 * 	This program is distributed in the hope that it will be useful,
 * 	but WITHOUT ANY WARRANTY; without even the implied warranty of
 * 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * 	GNU General Public License for more details.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

float get_volume_intersection_fromgrid(MOL2 *mol, float ****grids, int *types, float min_grid[3], float max_grid[3], int max_x, int max_y, int max_z, float spacing);

float get_color_volume_intersection_fromgrid(MOL2 *mol, float ****grids, int **types, float min_grid[3], float max_grid[3], int max_x, int max_y, int max_z, float spacing);




void gaussian_color_add_grid_forces(MOL2 **mymol, float ****grids, int **types, float min_grid[3], float max_grid[3], int max_x, int max_y, int max_z, float spacing)
{

        int i = 0, real_atoms = 0, out_flag = 0, l = 0, kpot = 0, ix, iy, iz;
        float energy = 0.0f;
        float energy2 = 0.0f, xa, ya, za, xd, yd, zd, q, x, y, z;
        float V000, V100, V010, V110, V001, V101, V011, V111;
	int nx = -1, ny = -1, nz = -1;
	int bestpot = -1;
	float e1 = 0.0, tmpe = 0.0f;
        register double u,   v,   w;
        register double p0u, p0v, p0w;
        register double p1u, p1v, p1w;
        register int u0,  v0,  w0;
        register int u1,  v1,  w1;
	float com[3];

        MOL2 *mol = NULL;

        mol = *mymol;
	real_atoms = mol->n_atoms;

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

	if( out_flag == 1)
		return;

        for ( i = 0; i < real_atoms; ++i) 
	    {
                e1 = 0.0f;
                energy = 0.0f;
                energy2 = 0.0f;

                x = mol->x[i];
                y = mol->y[i];
                z = mol->z[i];

                ix = (int)((mol->x[i] - min_grid[0]) / spacing);
                iy = (int)((mol->y[i] - min_grid[1]) / spacing);
                iz = (int)((mol->z[i] - min_grid[2]) / spacing);

                if ( ix >= 1 && iy >= 1 &&
                     iz >= 1 && ix < (max_x - 1) && iy < (max_y - 1) &&
                     iz < (max_z - 1) ) {

			bestpot = -1;
                        xd = min_grid[0];
                        yd = min_grid[1];
                        zd = min_grid[2];
                        w1  = (w0 = (int)(w = (x - xd) / spacing )) + 1;
                        p1w = 1.0L - (p0w = w - (float)w0);
                        v1  = (v0 = (int)(v = (y - yd) / spacing)) + 1;
                        p1v = 1.0L - (p0v = v - (float)v0);
                        u1  = (u0 = (int)(u = (z - zd) / spacing)) + 1;
                        p1u = 1.0L - (p0u = u - (float)u0);

                        for( l = 0; l < 5; l++)
                        {
                                tmpe = 0.0f;
                                if ( types[i][l] != 0)
                                {
                                        kpot = 6+l;
                                        tmpe = p1u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
                                        tmpe += p0u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u1 ];
                                        tmpe += p1u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u0 ];
                                        tmpe += p0u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u1 ];
                                        tmpe += p1u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
                                        tmpe += p0u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u1 ];
                                        tmpe += p1u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u0 ];
                                        tmpe += p0u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u1 ];
                                        if( tmpe < e1)
					                    {
                     					 bestpot = kpot;
                                         e1 = tmpe;
					                    }
                                }
                        }

			kpot = bestpot;
			if (kpot < 0)
			 return;

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
                }

        }

        for( i = 0; i < mol->n_rings; i++)
        {
                com[0] = com[1] = com[2] = 0.0f;
                for( l = 0; l < mol->rings[i][0]; l++)
                {
                        com[0] += mol->x[mol->rings[i][1+l]];
                        com[1] += mol->y[mol->rings[i][1+l]];
                        com[2] += mol->z[mol->rings[i][1+l]];
                }
                com[0] /= (float) mol->rings[i][0];
                com[1] /= (float) mol->rings[i][0];
                com[2] /= (float) mol->rings[i][0];

                ix = (int)((com[0] - min_grid[0]) / spacing);
                iy = (int)((com[1] - min_grid[1]) / spacing);
                iz = (int)((com[2] - min_grid[2]) / spacing);

                if ( ix >= 1 && iy >= 1 &&
                     iz >= 1 && ix < (max_x - 1) && iy < (max_y - 1) &&
                     iz < (max_z - 1) ) {

                        xd = min_grid[0];
                        yd = min_grid[1];
                        zd = min_grid[2];
                        w1  = (w0 = (int)(w = (x - xd) / spacing )) + 1;
                        p1w = 1.0L - (p0w = w - (float)w0);
                        v1  = (v0 = (int)(v = (y - yd) / spacing)) + 1;
                        p1v = 1.0L - (p0v = v - (float)v0);
                        u1  = (u0 = (int)(u = (z - zd) / spacing)) + 1;
                        p1u = 1.0L - (p0u = u - (float)u0);

                        kpot = 11;
                        tmpe = p1u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
                        tmpe += p0u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u1 ];
                        tmpe += p1u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u0 ];
                        tmpe += p0u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u1 ];
                        tmpe += p1u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
                        tmpe += p0u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u1 ];
                        tmpe += p1u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u0 ];
                        tmpe += p0u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u1 ];

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


	                for( l = 0; l < mol->rings[i][0]; l++)
       	        	{
				V000 /= (float) mol->rings[i][0];
				V100 /= (float) mol->rings[i][0];
				V010 /= (float) mol->rings[i][0];
				V110 /= (float) mol->rings[i][0];
				V001 /= (float) mol->rings[i][0];
				V101 /= (float) mol->rings[i][0];
				V011 /= (float) mol->rings[i][0];
				V111 /= (float) mol->rings[i][0];

	                        x = mol->x[mol->rings[i][1+l]];
	                        y = mol->y[mol->rings[i][1+l]];
	                        z = mol->z[mol->rings[i][1+l]];
		
	                        xd = x - xa;
	                        yd = y - ya;
	                        zd = z - za;

	                        mol->grads_X[l] += V000 * -1 * (1 - yd) * (1 - zd) +
                                           V100 * (1 - yd) * (1 - zd) +
                                           V010 * -1 * yd * (1 - zd) +
                                           V110 * yd * (1 - zd) +
                                           V001 * -1 * (1 - yd) * zd +
                                           V101 * (1 - yd) * zd +
                                           V011 * -1 * yd * zd +
                                           V111 * yd * zd;

	                        mol->grads_Y[l] += V000 * (1 - xd) * -1 * (1 - zd) +
                                           V100 * xd * -1 * (1 - zd) +
                                           V010 * (1 - xd) *  (1 - zd) +
                                           V110 * xd * (1 - zd) +
                                           V001 * (1 - xd) * -1 * zd +
                                           V101 * xd * -1 * zd +
                                           V011 * (1 - xd) * zd +
                                           V111 * xd * zd;

	                        mol->grads_Z[l] += V000 * (1 - xd) * (1 - yd) * -1 +
                                           V100 * xd * (1 - yd) * -1 +
                                           V010 * (1 - xd) * yd * -1 +
                                           V110 * xd * yd * -1 +
                                           V001 * (1 - xd) * (1 - yd) +
                                           V101 * xd * (1 - yd) +
                                           V011 * (1 - xd) * yd +
                                           V111 * xd * yd;

			}
                }
	}

        return;
}


void gaussian_add_grid_forces(MOL2 **mymol, float ****grids, int *types, float min_grid[3], float max_grid[3], int max_x, int max_y, int max_z, float spacing)
{

        int i = 0, real_atoms = 0, out_flag = 0, l = 0, kpot = 0, ix, iy, iz;
        float energy = 0.0f;
        float energy2 = 0.0f, xa, ya, za, xd, yd, zd, q, x, y, z;
        float V000, V100, V010, V110, V001, V101, V011, V111;
	int nx = -1, ny = -1, nz = -1;
        MOL2 *mol = NULL;

        mol = *mymol;
	real_atoms = mol->n_atoms;

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

	if( out_flag == 1)
		return;

        for ( i = 0; i < real_atoms; ++i) 
	{
                energy = 0.0f;
                energy2 = 0.0f;

                kpot = types[i];
                x = mol->x[i];
                y = mol->y[i];
                z = mol->z[i];

                ix = (int)((mol->x[i] - min_grid[0]) / spacing);
                iy = (int)((mol->y[i] - min_grid[1]) / spacing);
                iz = (int)((mol->z[i] - min_grid[2]) / spacing);

                if ( types[i] != -1 && ix >= 1 && iy >= 1 &&
                     iz >= 1 && ix < (max_x - 1) && iy < (max_y - 1) &&
                     iz < (max_z - 1) ) {

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
                }

        }

        return;
}




float get_color_intersection_at_point(MOL2 *a, int **types_a, float x, float y, float z, int type, int pocket )
{
        float vol = 0.0f, tmpvol = 0.0f;
        int i = 0, j = 0, l = 0;
        float alfa_a = 0.0f, alfa_b = 0.0f;
        const float kappa_i = 1.6258;
        const float pi_i = 2.7f;
        float kij = 0.0f;
        float dx = 0.0f, dy = 0.0f, dz = 0.0f, dist = 0.0f;
	int *type_a, type_b[6];
	int polar_a = 0, polar_b = 0;
	float weight = 1.0f;
	int flag = 0;
	float com[3];

	alfa_a = kappa_i / 1.0f;
        alfa_b = kappa_i / 1.0f;


	        for( i = 0; i < a->n_atoms; i++)
	        {

			type_a = types_a[i];
	
	                for( l = 0; l < 6; l++)
	                {
	                        type_b[l] = 0;
	                }

			type_b[type] = 1;


			flag = 0;
			for( l = 0; l < 6; l++)
			{
				if( type_a[l] == type_b[l] && type_a[l] != 0)
				{
				  flag = 1;
				}
			}
	  
	                 if( flag == 1)
			 {
	                        dx = a->x[i] - x;
	                        dy = a->y[i] - y;
	                        dz = a->z[i] - z;
	                        dist = (dx*dx)+(dy*dy)+(dz*dz);
			
	                        kij = exp(-((alfa_a*alfa_b*dist)/(alfa_a+alfa_b)));
				tmpvol = (pi_i*pi_i*kij* cbrt( (3.14159f/(alfa_a+alfa_b))*(3.14159f/(alfa_a+alfa_b))) );
				if( dist < 2.9) /* 1.7 A */
		                        vol += tmpvol;
	
			}
		} /* End of atoms */

		if (type == 5)
		{
	                for( i = 0; i < a->n_rings; i++)
	                {	
				com[0] = com[1] = com[2] = 0.0f;
				for( l = 0; l < a->rings[i][0]; l++)
				{
					com[0] += a->x[a->rings[i][1+l]];
					com[1] += a->y[a->rings[i][1+l]];
					com[2] += a->z[a->rings[i][1+l]];
				}
				com[0] /= (float) a->rings[i][0];
				com[1] /= (float) a->rings[i][0];
				com[2] /= (float) a->rings[i][0];

                                dx = com[0] - x;
                                dy = com[1] - y;
                                dz = com[2] - z;
                                dist = (dx*dx)+(dy*dy)+(dz*dz);

                                kij = exp(-((alfa_a*alfa_b*dist)/(alfa_a+alfa_b)));
                                tmpvol = (pi_i*pi_i*kij* cbrt( (3.14159f/(alfa_a+alfa_b))*(3.14159f/(alfa_a+alfa_b))) );
                                if( dist < 2.9) /* 1.7 A */
                                        vol += tmpvol;

			}
		}


	return -vol;

}



float get_color_intersection(MOL2 *a, MOL2 **b2, int **types_a, int **types_b, int pocket1, int pocket2)
{
        float vol = 0.0f, tmpvol = 0.0f;
        int i = 0, j = 0, l = 0;
        float alfa_a = 0.0f, alfa_b = 0.0f;
/*        const float kappa_i = 0.6604f;*/
        const float kappa_i = 1.6258;
        const float pi_i = 2.7f;
        float kij = 0.0f;
        float dx = 0.0f, dy = 0.0f, dz = 0.0f, dist = 0.0f;
	int *type_a, *type_b;
	float weight = 1.0f;
	int flag = 0;
        MOL2 *b = NULL;
        b = *b2;
	float com_a[3], com_b[3];

	alfa_a = kappa_i / 1.0f;
        alfa_b = kappa_i / 1.0f;

        for( i = 0; i < a->n_atoms; i++)
        {
	      type_a = types_a[i];
	      for( j = 0; j < b->n_atoms; j++)
	      {
		 type_b = types_b[j];

		 flag = 0;
		 for( l = 0; l < 6; l++)
		 {
			if( type_a[l] == type_b[l] && type_a[l] != 0)
			{
			  flag = 1;
			}
		 }
	  
                 if( flag == 1)
		 {
                        dx = a->x[i] - b->x[j];
                        dy = a->y[i] - b->y[j];
                        dz = a->z[i] - b->z[j];
                        dist = (dx*dx)+(dy*dy)+(dz*dz);
				
                        kij = exp(-((alfa_a*alfa_b*dist)/(alfa_a+alfa_b)));
			tmpvol = (pi_i*pi_i*kij* cbrt( (3.14159f/(alfa_a+alfa_b))*(3.14159f/(alfa_a+alfa_b))) );
			if( dist < 2.7)
			{
	                        vol += tmpvol;
	                        b->grads_X[j] += (2.0f*alfa_a*alfa_b*dx*tmpvol)/(alfa_a+alfa_b);
	                        b->grads_Y[j] += (2.0f*alfa_a*alfa_b*dy*tmpvol)/(alfa_a+alfa_b);
	                        b->grads_Z[j] += (2.0f*alfa_a*alfa_b*dz*tmpvol)/(alfa_a+alfa_b);
			}
		}else if( flag == 2){
                        dx = a->x[i] - b->x[j];
                        dy = a->y[i] - b->y[j];
                        dz = a->z[i] - b->z[j];
                        dist = (dx*dx)+(dy*dy)+(dz*dz);
                        kij = exp(-((alfa_a*alfa_b*dist)/(alfa_a+alfa_b)));
                        tmpvol = (pi_i*pi_i*kij* cbrt( (3.14159f/(alfa_a+alfa_b))*(3.14159f/(alfa_a+alfa_b))) );

                        vol -= tmpvol;

                        b->grads_X[j] -= (2.0f*alfa_a*alfa_b*dx*tmpvol)/(alfa_a+alfa_b);
                        b->grads_Y[j] -= (2.0f*alfa_a*alfa_b*dy*tmpvol)/(alfa_a+alfa_b);
                        b->grads_Z[j] -= (2.0f*alfa_a*alfa_b*dz*tmpvol)/(alfa_a+alfa_b);
	
		}
	      }
	}


        for( i = 0; i < a->n_rings; i++)
        {
                com_a[0] = com_a[1] = com_a[2] = 0.0f;
                for( l = 0; l < a->rings[i][0]; l++)
                {
                        com_a[0] += a->x[a->rings[i][1+l]];
                        com_a[1] += a->y[a->rings[i][1+l]];
                        com_a[2] += a->z[a->rings[i][1+l]];
                }
                com_a[0] /= (float) a->rings[i][0];
                com_a[1] /= (float) a->rings[i][0];
                com_a[2] /= (float) a->rings[i][0];


		for( j = 0; j < b->n_rings; j++)
		{

	                com_b[0] = com_b[1] = com_b[2] = 0.0f;
	                for( l = 0; l < b->rings[j][0]; l++)
	                {
	                        com_b[0] += b->x[b->rings[j][1+l]];
	                        com_b[1] += b->y[b->rings[j][1+l]];
	                        com_b[2] += b->z[b->rings[j][1+l]];
	                }
	                com_b[0] /= (float) b->rings[j][0];
	                com_b[1] /= (float) b->rings[j][0];
	                com_b[2] /= (float) b->rings[j][0];


	                dx = com_a[0] - com_b[0];
	                dy = com_a[1] - com_b[1];
	                dz = com_a[2] - com_b[2];
	                dist = (dx*dx)+(dy*dy)+(dz*dz);

		        kij = exp(-((alfa_a*alfa_b*dist)/(alfa_a+alfa_b)));
		        tmpvol = (pi_i*pi_i*kij* cbrt( (3.14159f/(alfa_a+alfa_b))*(3.14159f/(alfa_a+alfa_b))) );
		        if( dist < 2.9)
		        {
		                vol += tmpvol;
	                        for( l = 0; l < b->rings[j][0]; l++)
				{
		                	b->grads_X[b->rings[j][1+l]]+=((2.0f*alfa_a*alfa_b*dx*tmpvol)/(alfa_a+alfa_b))/b->rings[j][0]; 
			                b->grads_Y[b->rings[j][1+l]]+=((2.0f*alfa_a*alfa_b*dy*tmpvol)/(alfa_a+alfa_b))/b->rings[j][0];
			                b->grads_Z[b->rings[j][1+l]]+=((2.0f*alfa_a*alfa_b*dz*tmpvol)/(alfa_a+alfa_b))/b->rings[j][0];
				}
		        }


		}
        }


	*b2 = b;

	return -vol;

}

float get_volumen_intersection_at_point(MOL2 *a, float x, float y, float z, int type, int pocket)
{

        float vol = 0.0f, tmpvol = 0.0f;
        int i = 0, j = 0;
        float *alfa_a = NULL, alfa_b = 0.0f;
        const float kappa_i = 1.6258;
        const float pi_i = 2.7f;
        float kij = 0.0f;
        float dx = 0.0f, dy = 0.0f, dz = 0.0f, dist = 0.0f;

	if( type == 4)
	 return 0.0f;

        alfa_a = (float *) calloc(sizeof(float),a->n_atoms);

	if( pocket == 1)
	{
		
                for( i = 0; i < a->n_atoms; i++)
                {
                if( a->atoms[i] == 4)
                 alfa_a[i] = kappa_i / (1.70*1.70);
                else if( a->atoms[i] == 6)
                 alfa_a[i] = kappa_i / (1.60*1.60);
                else if( a->atoms[i] == 3)
                 alfa_a[i] = kappa_i / (1.65*1.65);
                else if( a->atoms[i] == 2)
                 alfa_a[i] = kappa_i / (1.60*1.60);
                else
                 alfa_a[i] = kappa_i / (1.70*1.70);
		}

	}else{
	        for( i = 0; i < a->n_atoms; i++)
	        {
                if( a->atoms[i] == 1)
                 alfa_a[i] = kappa_i / (1.70*1.70);
                else if( a->atoms[i] == 2)
                 alfa_a[i] = kappa_i / (1.60*1.60);
                else if( a->atoms[i] == 3)
                 alfa_a[i] = kappa_i / (1.65*1.65);
                else if( a->atoms[i] == 4)
                 alfa_a[i] = kappa_i / (1.00*1.00);
                else if( a->atoms[i] == 5)
                 alfa_a[i] = kappa_i / (1.90*1.90);
                else if( a->atoms[i] == 6)
                 alfa_a[i] = kappa_i / (1.90*1.90);
                else if( a->atoms[i] == 10)
                 alfa_a[i] = kappa_i / (1.30*1.30);
                else
                 alfa_a[i] = kappa_i / (1.70*1.70);
        	}
	}

        if( type  == 1)
         alfa_b = kappa_i / (1.70*1.70);
        else if( type == 2)
         alfa_b = kappa_i / (1.60*1.60);
        else if( type == 3)
         alfa_b = kappa_i / (1.65*1.65);
        else if( type == 4)
         alfa_b = kappa_i / (1.00*1.00);
        else if( type == 5)
         alfa_b = kappa_i / (1.90*1.90);
        else if( type == 6)
         alfa_b = kappa_i / (1.90*1.90);
        else if( type == 10)
         alfa_b = kappa_i / (1.30*1.30);
        else
         alfa_b = kappa_i / (1.70*1.70);

	vol = 0.0f;
        for( i = 0; i < a->n_atoms; i++)
        {
                if( a->atoms[i] == 4 && pocket != 1)
                 continue;

 	        dx = a->x[i] - x;
	        dy = a->y[i] - y;
	        dz = a->z[i] - z;
	        dist = (dx*dx)+(dy*dy)+(dz*dz);
	        kij = exp(-((alfa_a[i]*alfa_b*dist)/(alfa_a[i]+alfa_b)));
	        tmpvol = (pi_i*pi_i*kij* cbrt( (3.14159f/(alfa_a[i]+alfa_b))*(3.14159f/(alfa_a[i]+alfa_b))) );
	        vol += tmpvol;
        }

	free(alfa_a);
        return -vol;

}

float get_volumen_intersection(MOL2 *a, MOL2 **b2, int pocket1, int pocket2 )
{
	float vol = 0.0f, tmpvol = 0.0f;
	int i = 0, j = 0;
	float *alfa_a = NULL, *alfa_b = NULL;
/*	const float kappa_i = 0.6604f;*/
	const float kappa_i = 1.6258;
        const float pi_i = 2.7f;
	float kij = 0.0f;
	float dx = 0.0f, dy = 0.0f, dz = 0.0f, dist = 0.0f;
	MOL2 *b = NULL;
	b = *b2;

	alfa_a = (float *) calloc(sizeof(float),a->n_atoms);
        alfa_b = (float *) calloc(sizeof(float),b->n_atoms);



        if( pocket1 == 1)
        {

                for( i = 0; i < a->n_atoms; i++)
                {
/*                if( a->atoms[i] == 4)
                 alfa_a[i] = kappa_i / (1.70*1.70);
                else if( a->atoms[i] == 6)
                 alfa_a[i] = kappa_i / (1.60*1.60);
                else if( a->atoms[i] == 3)
                 alfa_a[i] = kappa_i / (1.65*1.65);
                else if( a->atoms[i] == 2)
                 alfa_a[i] = kappa_i / (1.60*1.60);
                else*/
                 alfa_a[i] = kappa_i / (1.70*1.70);
                }

        }else{

	        for( i = 0; i < a->n_atoms; i++)
	        {
			if( a->atoms[i] == 1)
			 alfa_a[i] = kappa_i / (1.70*1.70);
			else if( a->atoms[i] == 2)
	                 alfa_a[i] = kappa_i / (1.60*1.60);
	                else if( a->atoms[i] == 3)
	                 alfa_a[i] = kappa_i / (1.65*1.65);
	                else if( a->atoms[i] == 4)
	                 alfa_a[i] = kappa_i / (1.00*1.00);
	                else if( a->atoms[i] == 5)
	                 alfa_a[i] = kappa_i / (1.90*1.90);
	                else if( a->atoms[i] == 6)
	                 alfa_a[i] = kappa_i / (1.90*1.90);
	                else if( a->atoms[i] == 10)
	                 alfa_a[i] = kappa_i / (1.30*1.30);
			else
	                 alfa_a[i] = kappa_i / (1.70*1.70);

	        }
	}


        if( pocket2 == 1)
        {

                for( i = 0; i < b->n_atoms; i++)
                {
/*                if( b->atoms[i] == 4)
                 alfa_b[i] = kappa_i / (1.70*1.70);
                else if( b->atoms[i] == 6)
                 alfa_b[i] = kappa_i / (1.60*1.60);
                else if( b->atoms[i] == 3)
                 alfa_b[i] = kappa_i / (1.65*1.65);
                else if( b->atoms[i] == 2)
                 alfa_b[i] = kappa_i / (1.60*1.60);
                else*/
                 alfa_b[i] = kappa_i / (1.70*1.70);

                b->grads_X[i] = 0.0f;
                b->grads_Y[i] = 0.0f;
                b->grads_Z[i] = 0.0f;
       
		}
	}else{

	    if( pocket1 == 1) /* As we cannot be sure of the pocket types everything should be the same */
	    {
                for( i = 0; i < b->n_atoms; i++)
                {
                        alfa_b[i] = kappa_i / (1.70*1.70);
                        b->grads_X[i] = 0.0f;
                        b->grads_Y[i] = 0.0f;
                        b->grads_Z[i] = 0.0f;
		}
            }else{

	        for( i = 0; i < b->n_atoms; i++)
	        {
	                if( b->atoms[i] == 1)
	                 alfa_b[i] = kappa_i / (1.70*1.70);
	                else if( b->atoms[i] == 2)
	                 alfa_b[i] = kappa_i / (1.60*1.60);
	                else if( b->atoms[i] == 3)
	                 alfa_b[i] = kappa_i / (1.65*1.65);
	                else if( b->atoms[i] == 4)
	                 alfa_b[i] = kappa_i / (1.00*1.00);
	                else if( b->atoms[i] == 5)
	                 alfa_b[i] = kappa_i / (1.90*1.90);
	                else if( b->atoms[i] == 6)
	                 alfa_b[i] = kappa_i / (1.90*1.90);
	                else if( b->atoms[i] == 10)
	                 alfa_b[i] = kappa_i / (1.30*1.30);
	                else
	                 alfa_b[i] = kappa_i / (1.70*1.70);
		
			b->grads_X[i] = 0.0f;
	                b->grads_Y[i] = 0.0f;
	                b->grads_Z[i] = 0.0f;

	        }
            }
	}

	for( i = 0; i < a->n_atoms; i++)
	{
		if( a->atoms[i] == 4 && pocket1 != 1)
		 continue;

		for( j = 0; j < b->n_atoms; j++)
		{
	
                 if( b->atoms[j] == 4 && pocket2 != 1)
                  continue;

			dx = a->x[i] - b->x[j];
                        dy = a->y[i] - b->y[j];
                        dz = a->z[i] - b->z[j];
			dist = (dx*dx)+(dy*dy)+(dz*dz);
			kij = exp(-((alfa_a[i]*alfa_b[j]*dist)/(alfa_a[i]+alfa_b[j])));
			tmpvol = (pi_i*pi_i*kij* cbrt( (3.14159f/(alfa_a[i]+alfa_b[j]))*(3.14159f/(alfa_a[i]+alfa_b[j]))) );
/*			fprintf(stderr,"%f %f %f\n",pow(3.14159f/(alfa_a[i]+alfa_b[j]),3.0f/2.0f), cbrt( (3.14159f/(alfa_a[i]+alfa_b[j]))*(3.14159f/(alfa_a[i]+alfa_b[j]))),  cbrt(pow(3.14159f/(alfa_a[i]+alfa_b[j]),2.0f)));*/
			vol += tmpvol;
			b->grads_X[j] += (2.0f*alfa_a[i]*alfa_b[j]*dx*tmpvol)/(alfa_a[i]+alfa_b[j]);
                        b->grads_Y[j] += (2.0f*alfa_a[i]*alfa_b[j]*dy*tmpvol)/(alfa_a[i]+alfa_b[j]);
                        b->grads_Z[j] += (2.0f*alfa_a[i]*alfa_b[j]*dz*tmpvol)/(alfa_a[i]+alfa_b[j]);
		}
	}

	free(alfa_a);
	free(alfa_b);

	*b2 = b;

	return -vol;
}


/**
 *
 *	@brief Calculate the dot product (scalar product) of two vector
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param a Vector a
 *	@param b Vector b
 *	@param n Vector dimensions
 *	@return dot product
 *
 */
float dot_product(float *a, float *b, int n)
{
	int i = 0;
	float res = 0;

	for ( i = 0; i < n; ++i)
		res += a[i] * b[i];

	return res;
}

/**
 *
 *	@brief Update the approximate Heassian for the BFGS algorithm
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param hessian Approximate hessian
 *	@param p Gradient vector
 *	@param y Vars vector
 *	@param alfa step size
 *	@return 0 on success. 1 if not updated
 *
 */
int gau_update_hessian(float hessian[6][6], float p[6], float y[6], float alfa)
{
	float yp = 0.0f, yhy = 0.0f, r = 0.0f;
	float *clone_y = NULL, cum = 0.0f, ycloney = 0.0f;
	int n = 0, i = 0, j = 0;

	yp = dot_product(y, p, 6);

	if ( alfa * yp < 0.000001) return 1;

	clone_y = (float*)calloc( sizeof(float), 6);

	for ( i = 0; i < 6; ++i) {
		cum = 0.0f;
		for ( j = 0; j < 6; ++j)
			cum += hessian[i][j] * y[j];
		clone_y[i] = cum;
	}

	r = 1.0f / (alfa * yp);
	ycloney = -dot_product(y, clone_y, 6);

	for ( i = 0; i < 6; ++i) {
		for ( j = 0; j < 6; ++j)
			hessian[i][j] += alfa * r * ( clone_y[i] * p[j] + clone_y[j] * p[i]) + alfa * alfa * (r * r * ycloney + r) * p[i] * p[j];

	}

	free(clone_y);

	return 0;
}

/**
 *
 *	@brief Perform linear search in gradient oppsite direction
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mmol MOL2 with the molecule to minimize
 *	@param forces Vector with forces
 *	@param product Vector with the product
 *	@param last_energy Last call energy
 *	@param new_gradient To be fill with the new forces
 *	@param new_energy To be fill with the new energy
 *	@param grids Potential grids
 *	@param types Atom types for force field grids
 *	@param min_grid Minimum grid coordinates
 *	@param max_grid Maximum grid coordinates
 *	@param nx Closest grid point to atoms. X
 *	@param ny Closest grid point to atoms. Y
 *	@param nz Closest grid point to atoms. Z
 *	@param max_x Number of X grid points
 *	@param max_y Number of Y grid points
 *	@param max_z Number of Z grid points
 *	@param spacing Grid spacing
 *	@return best step
 *
 *
 */
float gau_get_alpha_with_linear_search(MOL2 *template, MOL2 **mmol, float forces[6], float product[6],
				   float last_energy, float new_gradient[6], float *new_energy, 
				   int **types1, int **types2, int pocket)
{

	MOL2 *mol = NULL;
	int max_steps = 10, i = 0, j = 0;
	float alfa = 1.0f;
	float pdotg = 0.0f;
	float *xc = NULL, *yc = NULL, *zc = NULL;
	float ten_forces[6];
	float force[6], torque[6];

	mol = *mmol;

	xc = (float*)calloc(sizeof(float), mol->n_atoms);
	yc = (float*)calloc(sizeof(float), mol->n_atoms);
	zc = (float*)calloc(sizeof(float), mol->n_atoms);


	pdotg = dot_product(product, forces, 6);
	*new_energy = 0.0f;
	for ( j = 0; j < mol->n_atoms; j++) {
		xc[j] = mol->x[j];
		yc[j] = mol->y[j];
		zc[j] = mol->z[j];
	}


	for ( i = 0; i < max_steps; i++) {
		for ( j = 0; j < mol->n_atoms; j++) {
			mol->x[j] = xc[j];
			mol->y[j] = yc[j];
			mol->z[j] = zc[j];
		}

		for ( j = 0; j < 6; j++)
			ten_forces[j] = product[j] * alfa;

		force[0] = ten_forces[0];
		force[1] = ten_forces[1];
		force[2] = ten_forces[2];
		torque[0] = ten_forces[3];
		torque[1] = ten_forces[4];
		torque[2] = ten_forces[5];

		transformate_rigid(&mol, force, torque, alfa);


                for (j = 0; j < mol->n_atoms; ++j) {
                        mol->grads_X[j] = 0.0;
                        mol->grads_Y[j] = 0.0;
                        mol->grads_Z[j] = 0.0;
                }

		*new_energy = get_volumen_intersection(template, &mol, pocket, 0);
		*new_energy += get_color_intersection(template,&mol,types1,types2, pocket, 0);

		grad_to_rigid(mol, force, torque);

		new_gradient[0] = force[0];
		new_gradient[1] = force[1];
		new_gradient[2] = force[2];
		new_gradient[3] = torque[0];
		new_gradient[4] = torque[1];
		new_gradient[5] = torque[2];


		if ( (*new_energy - last_energy) < 0.0001 * alfa * pdotg)
			break;

		alfa = alfa * 0.5f;

	}

	free(xc); free(yc); free(zc);
	return alfa;
}

float gau_get_alpha_with_linear_search_ingrid(MOL2 **mmol, 
					      float forces[6], float product[6], 
					      float last_energy, float new_gradient[6], 
					      float *new_energy,float ****grids, int *types,
                                              int **types2, float min_grid[3], float max_grid[3],
                                   	      int max_x, int max_y, int max_z,float spacing)
{

	MOL2 *mol = NULL;
	int max_steps = 10, i = 0, j = 0;
	float alfa = 1.0f;
	float pdotg = 0.0f;
	float *xc = NULL, *yc = NULL, *zc = NULL;
	float ten_forces[6];
	float force[6], torque[6];

	mol = *mmol;

	xc = (float*)calloc(sizeof(float), mol->n_atoms);
	yc = (float*)calloc(sizeof(float), mol->n_atoms);
	zc = (float*)calloc(sizeof(float), mol->n_atoms);


	pdotg = dot_product(product, forces, 6);
	*new_energy = 0.0f;
	for ( j = 0; j < mol->n_atoms; j++) {
		xc[j] = mol->x[j];
		yc[j] = mol->y[j];
		zc[j] = mol->z[j];
	}


	for ( i = 0; i < max_steps; i++) {
		for ( j = 0; j < mol->n_atoms; j++) {
			mol->x[j] = xc[j];
			mol->y[j] = yc[j];
			mol->z[j] = zc[j];
		}

		for ( j = 0; j < 6; j++)
			ten_forces[j] = product[j] * alfa;

		force[0] = ten_forces[0];
		force[1] = ten_forces[1];
		force[2] = ten_forces[2];
		torque[0] = ten_forces[3];
		torque[1] = ten_forces[4];
		torque[2] = ten_forces[5];

		transformate_rigid(&mol, force, torque, alfa);


                for (j = 0; j < mol->n_atoms; ++j) {
                        mol->grads_X[j] = 0.0;
                        mol->grads_Y[j] = 0.0;
                        mol->grads_Z[j] = 0.0;
                }

	        *new_energy = get_volume_intersection_fromgrid(mol, grids, types, min_grid, max_grid, max_x, max_y, max_z, spacing);
	        *new_energy += get_color_volume_intersection_fromgrid(mol, grids, types2, min_grid, max_grid, max_x, max_y, max_z, spacing);


/*		*new_energy = get_volumen_intersection(template, &mol);
		*new_energy += get_color_intersection(template,&mol);*/

		gaussian_add_grid_forces(&mol, grids, types, min_grid, max_grid, max_x, max_y,
                                   max_z, spacing);

                gaussian_color_add_grid_forces(&mol, grids, types2, min_grid, max_grid, max_x, max_y,max_z, spacing);

		grad_to_rigid(mol, force, torque);

		new_gradient[0] = force[0];
		new_gradient[1] = force[1];
		new_gradient[2] = force[2];
		new_gradient[3] = torque[0];
		new_gradient[4] = torque[1];
		new_gradient[5] = torque[2];


		if ( (*new_energy - last_energy) < 0.0001 * alfa * pdotg)
			break;

		alfa = alfa * 0.5f;

	}

	free(xc); free(yc); free(zc);
	return alfa;
}




/**
 *
 *	@brief Rigid body minimizer using the BFGS algorithm
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param template MOL2 with the template
 *	@param mymol MOL2 with the molecule
 *	@param defsteps Number of steps to perform
 *	@return 0 on success
 */
int gau_minimizer_bfgs(MOL2 *template, MOL2 **mymol, int defsteps, int **types1, int **types2, int pocket)
{
	MOL2 *mol = NULL;
	int n = 0;
	float *x = NULL, *y = NULL, *z = NULL;
	float hessian[6][6];
	int i = 0, j = 0, k = 0, l = 0;
	float energy0 = 0.0f, energy1 = 0.0f, alfa = 0.0f;
	float original_gradient[6], force[3], torque[3];
	float gradient[6], tmp = 0.0f;
	float product[6], prev_energy = 0.0f;
	float new_gradient[6], new_energy = 0.0f;
	float ygradient[6];

	mol = *mymol;

	n = 6;
	x = (float*)calloc(sizeof(float), mol->n_atoms);
	y = (float*)calloc(sizeof(float), mol->n_atoms);
	z = (float*)calloc(sizeof(float), mol->n_atoms);


	/* Obtain energy and rigid-solid gradients */
        for (j = 0; j < mol->n_atoms; ++j) {
                mol->grads_X[j] = 0.0;
                mol->grads_Y[j] = 0.0;
                mol->grads_Z[j] = 0.0;
                x[j] = mol->x[j];
                y[j] = mol->y[j];
                z[j] = mol->z[j];

        }

	energy0 = get_volumen_intersection(template, &mol, pocket, 0);
	energy0 += get_color_intersection(template,&mol, types1, types2, pocket, 0);

	grad_to_rigid(mol, force, torque);
	original_gradient[0] = force[0];
	original_gradient[1] = force[1];
	original_gradient[2] = force[2];
	original_gradient[3] = torque[0];
	original_gradient[4] = torque[1];
	original_gradient[5] = torque[2];

	for ( j = 0; j < 6; ++j)
		gradient[j] = original_gradient[j];

#ifdef DEBUG
	fprintf(stderr, "Iteration 0. ENERGY: %f. Forces: %f %f %f %f %f %f\n", energy0, gradient[0], gradient[1], gradient[2], gradient[3], gradient[4], gradient[5]);
#endif
	/* Initialize Hessian with Identity matrix */
	for ( j = 0; j < n; ++j) {
		for ( i = 0; i < n; ++i)
			hessian[j][i] = 0.0f;
		hessian[j][j] = 1.0f;
	}


	prev_energy = energy0;
	for ( i = 0; i < defsteps; i++) {
		for (j = 0; j < mol->n_atoms; ++j) {
			x[j] = mol->x[j];
			y[j] = mol->y[j];
			z[j] = mol->z[j];
		}
		for ( j = 0; j < n; ++j) {
			tmp = 0.0f;
			for ( k = 0; k < n; ++k)
				tmp += hessian[j][k] * gradient[k];
			product[j] = -tmp;
		}

		energy1 = 0.0f;
		new_energy = 0.0f;
		alfa = gau_get_alpha_with_linear_search(template, &mol, gradient, product,
						    prev_energy, new_gradient, &new_energy, types1, types2, pocket);


		for ( j = 0; j < 6; ++j)
			ygradient[j] = -gradient[j];


		if ( prev_energy < new_energy ) {
			for (j = 0; j < mol->n_atoms; ++j) {
				mol->x[j] = x[j];
				mol->y[j] = y[j];
				mol->z[j] = z[j];
			}
/*			fprintf(stderr, "Energy rises. Now is: %f. Used to be: %f Alfa was: %f\n",new_energy,prev_energy,alfa);*/
			break;
		}

		prev_energy = new_energy;

#ifdef DEBUG
		fprintf(stderr, "Iteration %i. Energy: %e. Step: %e. Forces: %e %e %e %e %e %e\n", i + 1, prev_energy, alfa, new_gradient[0], new_gradient[1], new_gradient[2], new_gradient[3], new_gradient[4], new_gradient[5]);
#endif

		if (!(sqrt(dot_product(gradient, gradient, n)) >= 1e-5)) break;

		for ( j = 0; j < 6; ++j)
			gradient[j] = new_gradient[j];

		if (i == 0) {
			if ( (tmp = dot_product(ygradient, ygradient, 6)) < 0.00000001)
				for ( j = 0; j < 6; ++j)
					hessian[j][j] = alfa * dot_product(ygradient, product, 6) / tmp;
		}

		gau_update_hessian(hessian, product, ygradient, alfa);

	}


	free(x);
	free(y);
	free(z);
	return 0;
}


int gau_minimizer_bfgs_ingrid(MOL2 **mymol, int defsteps,float ****grids, int *types,int **types2, float min_grid[3], float max_grid[3],int max_x, int max_y, int max_z,float spacing)
{
	MOL2 *mol = NULL;
	int n = 0;
	float *x = NULL, *y = NULL, *z = NULL;
	float hessian[6][6];
	int i = 0, j = 0, k = 0, l = 0;
	float energy0 = 0.0f, energy1 = 0.0f, alfa = 0.0f;
	float original_gradient[6], force[3], torque[3];
	float gradient[6], tmp = 0.0f;
	float product[6], prev_energy = 0.0f;
	float new_gradient[6], new_energy = 0.0f;
	float ygradient[6];

	mol = *mymol;

	n = 6;
	x = (float*)calloc(sizeof(float), mol->n_atoms);
	y = (float*)calloc(sizeof(float), mol->n_atoms);
	z = (float*)calloc(sizeof(float), mol->n_atoms);


	/* Obtain energy and rigid-solid gradients */
        for (j = 0; j < mol->n_atoms; ++j) {
                mol->grads_X[j] = 0.0;
                mol->grads_Y[j] = 0.0;
                mol->grads_Z[j] = 0.0;
                x[j] = mol->x[j];
                y[j] = mol->y[j];
                z[j] = mol->z[j];

        }

/*	energy0 = get_volumen_intersection(template, &mol);
	energy0 += get_color_intersection(template,&mol);*/

        energy0 = get_volume_intersection_fromgrid(mol, grids, types, min_grid, max_grid, max_x, max_y, max_z, spacing);
        energy0 += get_color_volume_intersection_fromgrid(mol, grids, types2, min_grid, max_grid, max_x, max_y, max_z, spacing);

        gaussian_add_grid_forces(&mol, grids, types, min_grid, max_grid, max_x, max_y,
                                   max_z, spacing);

        gaussian_color_add_grid_forces(&mol, grids, types2, min_grid, max_grid, max_x, max_y,
                                   max_z, spacing);

	grad_to_rigid(mol, force, torque);
	original_gradient[0] = force[0];
	original_gradient[1] = force[1];
	original_gradient[2] = force[2];
	original_gradient[3] = torque[0];
	original_gradient[4] = torque[1];
	original_gradient[5] = torque[2];

	for ( j = 0; j < 6; ++j)
		gradient[j] = original_gradient[j];

#ifdef DEBUG
	fprintf(stderr, "Iteration 0. ENERGY: %f. Forces: %f %f %f %f %f %f\n", energy0, gradient[0], gradient[1], gradient[2], gradient[3], gradient[4], gradient[5]);
#endif
	/* Initialize Hessian with Identity matrix */
	for ( j = 0; j < n; ++j) {
		for ( i = 0; i < n; ++i)
			hessian[j][i] = 0.0f;
		hessian[j][j] = 1.0f;
	}


	prev_energy = energy0;
	for ( i = 0; i < defsteps; i++) {
		for (j = 0; j < mol->n_atoms; ++j) {
			x[j] = mol->x[j];
			y[j] = mol->y[j];
			z[j] = mol->z[j];
		}
		for ( j = 0; j < n; ++j) {
			tmp = 0.0f;
			for ( k = 0; k < n; ++k)
				tmp += hessian[j][k] * gradient[k];
			product[j] = -tmp;
		}

		energy1 = 0.0f;
		new_energy = 0.0f;
		alfa = gau_get_alpha_with_linear_search_ingrid(&mol, gradient, product,
						    prev_energy, new_gradient, &new_energy, grids,
						    types, types2, min_grid, max_grid, max_x,
						    max_y, max_z, spacing);


		for ( j = 0; j < 6; ++j)
			ygradient[j] = -gradient[j];


		if ( prev_energy < new_energy ) {
			for (j = 0; j < mol->n_atoms; ++j) {
				mol->x[j] = x[j];
				mol->y[j] = y[j];
				mol->z[j] = z[j];
			}
/*			fprintf(stderr, "Energy rises. Now is: %f. Used to be: %f Alfa was: %f\n",new_energy,prev_energy,alfa);*/
			break;
		}

		prev_energy = new_energy;

#ifdef DEBUG
		fprintf(stderr, "Iteration %i. Energy: %e. Step: %e. Forces: %e %e %e %e %e %e\n", i + 1, prev_energy, alfa, new_gradient[0], new_gradient[1], new_gradient[2], new_gradient[3], new_gradient[4], new_gradient[5]);
#endif

		if (!(sqrt(dot_product(gradient, gradient, n)) >= 1e-5)) break;

		for ( j = 0; j < 6; ++j)
			gradient[j] = new_gradient[j];

		if (i == 0) {
			if ( (tmp = dot_product(ygradient, ygradient, 6)) < 0.00000001)
				for ( j = 0; j < 6; ++j)
					hessian[j][j] = alfa * dot_product(ygradient, product, 6) / tmp;
		}

		gau_update_hessian(hessian, product, ygradient, alfa);

	}


	free(x);
	free(y);
	free(z);
	return 0;
}


/**
 *
 *      @brief Perform a translation operation on a molecule
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @param mmol MOL2 with the molecule
 *      @param x New x coordinate for the center of mass
 *      @param y New y coordinate for the center of mass
 *      @param z New z coordinate for the center of mass
 *      @return 0 on success
 *
 */
int gau_translate_molecule(MOL2 **mmol, float x, float y, float z)
{
        int i = 0;
        float center[3];
        float vec[3];
        MOL2 *mol;

        mol = *mmol;

        center[0] = center[1] = center[2] = 0.0f;
        for ( i = 0; i < mol->n_atoms; ++i) {
                center[0] = center[0] + mol->x[i];
                center[1] = center[1] + mol->y[i];
                center[2] = center[2] + mol->z[i];
        }

        center[0] = center[0] / (float)mol->n_atoms;
        center[1] = center[1] / (float)mol->n_atoms;
        center[2] = center[2] / (float)mol->n_atoms;

        vec[0] = x - center[0];
        vec[1] = y - center[1];
        vec[2] = z - center[2];

        for ( i = 0; i < mol->n_atoms; ++i) {
                mol->x[i] = mol->x[i] + vec[0];
                mol->y[i] = mol->y[i] + vec[1];
                mol->z[i] = mol->z[i] + vec[2];
        }

        *mmol = mol;
        return 0;
}




float get_volume_intersection_fromgrid(MOL2 *mol, float ****grids, int *types, float min_grid[3], float max_grid[3], int max_x, int max_y, int max_z, float spacing)
{

        int i = 0, real_atoms = 0, out_flag = 0, l = 0;
        float energy = 0.0f;
        int icube = 0;
        int kpot = 0, ix = 0, iy = 0, iz = 0;
        float xj = 0.0f, yj = 0.0f, zj = 0.0f, dis = 0.0f;
        float q = 0.0f, xa = 0.0f, ya = 0.0f, za = 0.0f;
        int m = 0, n = 0, ia = 0, ib = 0, ic = 0;
        float hdist = 0.0f, x = 0.0f, y = 0.0f, z = 0.0f;
        float a = 0, b = 0, c = 0, wi = 0.0f;
        float xd, yd, zd, i1, i2, j1, j2, w2, e1, e2, hdis, dx, dy, dz;
        int pi = 0, pj = 0, pk = 0;
        float dist = 0.0f, dist2 = 0.0f,  enorm2 = 0.0f;
        float V000, V100, V010, V110, V001, V101, V011, V111;
        register double u,   v,   w;
        register double p0u, p0v, p0w;
        register double p1u, p1v, p1w;
        register int u0,  v0,  w0;
        register int u1,  v1,  w1;

	real_atoms = mol->n_atoms;

        out_flag = 0;
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


	if( out_flag != 0)
         return 9999.9f;

        energy = 0.0f;
        e1 = e2 = 0.0f;
        for ( i = 0; i < real_atoms; ++i) {
                e1 = 0.0f;
                kpot = types[i];
                x = mol->x[i];
                y = mol->y[i];
                z = mol->z[i];

                ix = (int)roundf((mol->x[i] - min_grid[0]) / spacing);
                iy = (int)roundf((mol->y[i] - min_grid[1]) / spacing);
                iz = (int)roundf((mol->z[i] - min_grid[2]) / spacing);

                if ( !(ix >= 1 && iy >= 1 &&
                       iz >= 1 && ix < (max_x - 1) && iy < (max_y - 1) &&
                       iz < (max_z - 1))) {
                        return 9999.9f;
                }


                if ( types[i] != -1 && (ix >= 1 && iy >= 1 &&
                                        iz >= 1 && ix < max_x - 1 && iy < max_y - 1 &&
                                        iz < max_z - 1)) {

                        xd = min_grid[0];
                        yd = min_grid[1];
                        zd = min_grid[2];
                        w1  = (w0 = (int)(w = (x - xd) / spacing )) + 1;
                        p1w = 1.0L - (p0w = w - (float)w0);
                        v1  = (v0 = (int)(v = (y - yd) / spacing)) + 1;
                        p1v = 1.0L - (p0v = v - (float)v0);
                        u1  = (u0 = (int)(u = (z - zd) / spacing)) + 1;
                        p1u = 1.0L - (p0u = u - (float)u0);

			if( w0 >= max_x || w0 < 0)
				return 9999.9f;
                        if( w1 >= max_x || w1 < 0)
				return 9999.9f;
			if( v0 >= max_y || v0 < 0)
                                return 9999.9f;
			if( v1 >= max_y || v1 < 0)
				return 9999.9f;
			if( u0 >= max_z || u0 < 0)
				return 9999.9f;
			if( u1 >= max_z || u1 < 0)
				return 9999.9f;

                        e1 = p1u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
                        e1 += p0u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u1 ];
                        e1 += p1u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u0 ];
                        e1 += p0u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u1 ];
                        e1 += p1u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
                        e1 += p0u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u1 ];
                        e1 += p1u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u0 ];
                        e1 += p0u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u1 ];
                        energy += e1;
                }
        }


        return energy;
}


float get_color_volume_intersection_fromgrid(MOL2 *mol, float ****grids, int **types, float min_grid[3], float max_grid[3], int max_x, int max_y, int max_z, float spacing)
{

        int i = 0, real_atoms = 0, out_flag = 0, l = 0;
        float energy = 0.0f;
        int icube = 0;
        int kpot = 0, ix = 0, iy = 0, iz = 0;
        float xj = 0.0f, yj = 0.0f, zj = 0.0f, dis = 0.0f;
        float q = 0.0f, xa = 0.0f, ya = 0.0f, za = 0.0f;
        int m = 0, n = 0, ia = 0, ib = 0, ic = 0;
        float hdist = 0.0f, x = 0.0f, y = 0.0f, z = 0.0f;
        float a = 0, b = 0, c = 0, wi = 0.0f;
        float xd, yd, zd, i1, i2, j1, j2, w2, e1, e2, hdis, dx, dy, dz;
        int pi = 0, pj = 0, pk = 0;
        float dist = 0.0f, dist2 = 0.0f,  enorm2 = 0.0f;
        float V000, V100, V010, V110, V001, V101, V011, V111;
        register double u,   v,   w;
        register double p0u, p0v, p0w;
        register double p1u, p1v, p1w;
        register int u0,  v0,  w0;
        register int u1,  v1,  w1;
	float tmpe = 0.0f;
	float com[3];

	real_atoms = mol->n_atoms;

        out_flag = 0;
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


	if( out_flag != 0)
         return 9999.9f;

        energy = 0.0f;
        e1 = e2 = 0.0f;
        for ( i = 0; i < real_atoms; ++i) 
	{
                e1 = 0.0f;
		e2 = 0.0f;
                x = mol->x[i];
                y = mol->y[i];
                z = mol->z[i];

                ix = (int)roundf((mol->x[i] - min_grid[0]) / spacing);
                iy = (int)roundf((mol->y[i] - min_grid[1]) / spacing);
                iz = (int)roundf((mol->z[i] - min_grid[2]) / spacing);

                if ( !(ix >= 1 && iy >= 1 &&
                       iz >= 1 && ix < max_x - 1 && iy < max_y - 1 &&
                       iz < max_z - 1)) {
                        return 9999.9f;
                }


                if ( types[i] != -1 && (ix >= 1 && iy >= 1 &&
                                        iz >= 1 && ix < max_x - 1 && iy < max_y - 1 &&
                                        iz < max_z - 1)) 
		{

                        xd = min_grid[0];
                        yd = min_grid[1];
                        zd = min_grid[2];
                        w1  = (w0 = (int)(w = (x - xd) / spacing )) + 1;
                        p1w = 1.0L - (p0w = w - (float)w0);
                        v1  = (v0 = (int)(v = (y - yd) / spacing)) + 1;
                        p1v = 1.0L - (p0v = v - (float)v0);
                        u1  = (u0 = (int)(u = (z - zd) / spacing)) + 1;
                        p1u = 1.0L - (p0u = u - (float)u0);

            e1 = 0;
			for( l = 0; l < 5; l++)
			{
				tmpe = 0.0f;
				if ( types[i][l] != 0)
				{
					kpot = 6+l;
		                        tmpe = p1u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
		                        tmpe += p0u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u1 ];
		                        tmpe += p1u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u0 ];
		                        tmpe += p0u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u1 ];
		                        tmpe += p1u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
		                        tmpe += p0u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u1 ];
		                        tmpe += p1u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u0 ];
		                        tmpe += p0u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u1 ];
					if( tmpe < e1)
					 e1 = tmpe;
				}
			}
                        energy += e1;
                }
        }

	for( i = 0; i < mol->n_rings; i++)
	{
            e1 = 0;
	        com[0] = com[1] = com[2] = 0.0f;
	        for( l = 0; l < mol->rings[i][0]; l++)
	        {
	                com[0] += mol->x[mol->rings[i][1+l]];
	                com[1] += mol->y[mol->rings[i][1+l]];
	                com[2] += mol->z[mol->rings[i][1+l]];
	        }
	        com[0] /= (float) mol->rings[i][0];
	        com[1] /= (float) mol->rings[i][0];
	        com[2] /= (float) mol->rings[i][0];

                ix = (int)roundf((com[0] - min_grid[0]) / spacing);
                iy = (int)roundf((com[1] - min_grid[1]) / spacing);
                iz = (int)roundf((com[2] - min_grid[2]) / spacing);

                if ( !(ix >= 1 && iy >= 1 &&
                       iz >= 1 && ix < max_x - 1 && iy < max_y - 1 &&
                       iz < max_z - 1)) {
                        return 9999.9f;
                }


                if ( (ix >= 1 && iy >= 1 &&
                   iz >= 1 && ix < max_x - 1 && iy < max_y - 1 &&
                   iz < max_z - 1))
                {

                        xd = min_grid[0];
                        yd = min_grid[1];
                        zd = min_grid[2];
                        w1  = (w0 = (int)(w = (x - xd) / spacing )) + 1;
                        p1w = 1.0L - (p0w = w - (float)w0);
                        v1  = (v0 = (int)(v = (y - yd) / spacing)) + 1;
                        p1v = 1.0L - (p0v = v - (float)v0);
                        u1  = (u0 = (int)(u = (z - zd) / spacing)) + 1;
                        p1u = 1.0L - (p0u = u - (float)u0);

                        tmpe = 0.0f;
                        kpot = 11; /* Rings */
                        tmpe = p1u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
			tmpe += p0u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u1 ];
			tmpe += p1u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u0 ];
			tmpe += p0u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u1 ];
			tmpe += p1u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
			tmpe += p0u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u1 ];
			tmpe += p1u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u0 ];
			tmpe += p0u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u1 ];
                           e1 = tmpe;
                }
                energy += e1;


	}

        return energy;
}

