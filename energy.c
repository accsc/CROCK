/**
 *
 *      @file energy.c
 *      @brief Energy module for grid and off-grid calculations
 *
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @date 01/10/2011
 *
 * 	This program is free software; you can redistribute it and/or modify
 * 	it under the terms of the GNU General Public License as published by
 * 	the Free Software Foundation version 2 of the License.
 *
 * 	This program is distributed in the hope that it will be useful,
 * 	but WITHOUT ANY WARRANTY; without even the implied warranty of
 * 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * 	GNU General Public License for more details.
 *
 */

float get_ff_grid_energy(MOL2 *mol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing);


/**
 * @brief Get the energy of the molecule in the GAFF force field
 *
 * @param mymol Molecule
 * @param tenergy Returned energy of the molecule
 * @param tortype1 Vector for temporal torsional type calculation
 * @param flg Flag for specific calculations. 0. Only vdw+elec without gradient; 1. All terms with gradient; 2. All terms without gradient
 * @return 0 if success. <0 if error.
 *
 */
int get_energy(MOL2 **mymol, double *tenergy, int **tortype1, int flg)
{
	MOL2 *mols;
	float dx, dy, dz;
	float df;
	float ddf1;
	float df1;
	float e1;
	float dgx, dgy, dgz;
	float dtxi, dtxj, dtyi, dtyj, dtzi, dtzj;
	float dfx, dfy, dfz;
	float gaa, ra2, rb2, gbb;
	float fg, hg, fga, hgb;
	float r2, r;
	float RJR, RIR, cst;
	int *tor_type;
	int retval = 0;

	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;

	float pn, pk, phi0;
	int cur_angle = 0;

	int t1, t2, t3, t4;

	float A, B;
	float sigma6;
	float epsilon;
	float tmpE;
	float vdwE;
	float r6, r12;

	float bondE = 0;
	float angleE = 0;
	float torE = 0;
	double totalE = 0;
	double vec1[3], vec2[3];
	double uvec1[3], uvec2[3];
	double *lvec1, *lvec2;
	double tor1;
	double tor2;
	double tor3;
	float thetha;
	float dist = 0.0f;
	float dist2 = 0.0f;
	float distr = 0.0f;
	int ttor = 0;


	/* FLAGS */
	int bondflag;
	int angleflag;
	int torflag;
	int gradflag;
	int vdwflag;
	int verboseflag;

	mols = *mymol;
	tor_type = *tortype1;

	verboseflag = 0;

	if ( flg == 1) { /* Internal energy terms on */
		bondflag = 1;
		torflag = 1;
		vdwflag = 1;
		gradflag = 1;
		angleflag = 1;
	}else if ( flg == 0) { /* Fast energy evaluation for vdw */
		bondflag = 0;
		torflag = 0;
		vdwflag = 1;
		gradflag = 0;
		angleflag = 0;
		verboseflag = 1;
	}else if ( flg == 2) { /* Internal without gradient  */
		bondflag = 1;
		torflag = 1;
		vdwflag = 1;
		gradflag = 0;
		angleflag = 1;
	}else{  /* Only Energy to grads  */
		bondflag = 1;
		torflag = 1;
		vdwflag = 1;
		gradflag = 0;
		angleflag = 1;
	}


/******************************************************************************
 *
 *                                 BONDS STUFF
 *
 *****************************************************************************/

	if ( bondflag == 1) {
		for (j = 0; j < mols->n_bonds; ++j) {

			if ( mols->gaff_types[ mols->bond_a1[j] - 1 ] < mols->gaff_types[ mols->bond_a2[j] - 1 ]) {
				t1 = mols->gaff_types[ mols->bond_a1[j] - 1 ] - 1;
				t2 = mols->gaff_types[ mols->bond_a2[j] - 1 ] - 1;
			}else{
				t1 = mols->gaff_types[ mols->bond_a2[j] - 1 ] - 1;
				t2 = mols->gaff_types[ mols->bond_a1[j] - 1 ] - 1;
			}

#ifdef DEBUG
			printf("%i-%i | %i-%i %f %f.\n", mols->bond_a1[j] - 1, mols->bond_a2[j] - 1, t1, t2, mols->bond_dist[j], bonds[t1][t2][1]);
			printf("%f * %f = %f.\n", bonds[t1][t2][0], pow( mols->bond_dist[j] - bonds[t1][t2][1], 2), bonds[t1][t2][0] * pow( mols->bond_dist[j] - bonds[t1][t2][1], 2));
#endif

			dx = mols->x[ mols->bond_a1[j] - 1 ] - mols->x[ mols->bond_a2[j] - 1 ];
			dy = mols->y[ mols->bond_a1[j] - 1 ] - mols->y[ mols->bond_a2[j] - 1 ];
			dz = mols->z[ mols->bond_a1[j] - 1 ] - mols->z[ mols->bond_a2[j] - 1 ];
			mols->bond_dist[j] = sqrt(dx * dx + dy * dy + dz * dz);
			df = bonds[t1][t2][0] * ( mols->bond_dist[j] - bonds[t1][t2][1]);

			bondE = bondE +  ( df * ( mols->bond_dist[j] - bonds[t1][t2][1]) );


			if ( gradflag == 1) {
				df *= 2 / mols->bond_dist[j];
				dx *= df;
				dy *= df;
				dz *= df;
				mols->grads_X[ mols->bond_a1[j] - 1] += dx;
				mols->grads_X[ mols->bond_a2[j] - 1] -= dx;

				mols->grads_Y[ mols->bond_a1[j] - 1] += dy;
				mols->grads_Y[ mols->bond_a2[j] - 1] -= dy;

				mols->grads_Z[ mols->bond_a1[j] - 1] += dz;
				mols->grads_Z[ mols->bond_a2[j] - 1] -= dz;

			}

		}

	}

/******************************************************************************
*
*                               ANGLES STUFF
*
******************************************************************************/

	if ( angleflag == 1) {
		for (j = 0; j < mols->n_angles; ++j) {


			t1 = mols->gaff_types[ mols->ia[j] ] - 1;
			t2 = mols->gaff_types[ mols->ja[j] ] - 1;
			t3 = mols->gaff_types[ mols->ka[j] ] - 1;

			vec1[0] = mols->x[ mols->ia[j] ] - mols->x[ mols->ja[j] ];
			vec1[1] = mols->y[ mols->ia[j] ] - mols->y[ mols->ja[j] ];
			vec1[2] = mols->z[ mols->ia[j] ] - mols->z[ mols->ja[j] ];

			vec2[0] = mols->x[ mols->ka[j] ] - mols->x[ mols->ja[j] ];
			vec2[1] = mols->y[ mols->ka[j] ] - mols->y[ mols->ja[j] ];
			vec2[2] = mols->z[ mols->ka[j] ] - mols->z[ mols->ja[j] ];

			RIR = sqrt((vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2]));
			RJR = sqrt((vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2]));
			cst = (vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2])  / ( RJR * RIR);
			thetha = acos( cst );
			dx = ( thetha - (angl[t1][t2][t3][1] * G_PI / 180));
			df =  angl[t1][t2][t3][0] * dx;
			angleE = angleE +  ( df * dx);
			if ( gradflag == 1) {
				df *= 2;

				dtxi = (cst * RJR * vec1[0] + RIR * (0 - vec2[0])) / ( sin(thetha) * RIR * RIR * RJR);
				dtyi = (cst * RJR * vec1[1] + RIR * (0 - vec2[1])) / ( sin(thetha) * RIR * RIR * RJR);
				dtzi = (cst * RJR * vec1[2] + RIR * (0 - vec2[2])) / ( sin(thetha) * RIR * RIR * RJR);
				mols->grads_X[ mols->ia[j] ] += ( df  * dtxi);
				mols->grads_Y[ mols->ia[j] ] += ( df  * dtyi);
				mols->grads_Z[ mols->ia[j] ] += ( df  * dtzi);

				dtxi = ((cst * RIR * vec2[0]) + RJR * (0 - vec1[0])) / (sin(thetha) * RIR * RJR * RJR);
				dtyi = ((cst * RIR * vec2[1]) + RJR * (0 - vec1[1]))  / (sin(thetha) * RIR * RJR * RJR);
				dtzi = ((cst * RIR * vec2[2]) + RJR * (0 - vec1[2]))  / (sin(thetha) * RIR * RJR * RJR);

				mols->grads_X[ mols->ka[j] ] += ( df  * dtxi);
				mols->grads_Y[ mols->ka[j] ] += ( df  * dtyi);
				mols->grads_Z[ mols->ka[j] ] += ( df  * dtzi);

				dtxi = ((cst * ( RIR * RIR * (0 - vec2[0]) + RJR * RJR * (0 - vec1[0]))) / (sin(thetha) * RIR * RIR * RJR * RJR)) + (( vec1[0] + vec2[0] ) / (sin(thetha) * RJR * RIR));
				dtyi = ((cst * ( RIR * RIR * (0 - vec2[1]) + RJR * RJR * (0 - vec1[1]))) / (sin(thetha) * RIR * RIR * RJR * RJR)) + (( vec1[1] + vec2[1] ) / (sin(thetha) * RJR * RIR));
				dtzi = ((cst * ( RIR * RIR * (0 - vec2[2]) + RJR * RJR * (0 - vec1[2]))) / (sin(thetha) * RIR * RIR * RJR * RJR)) + (( vec1[2] + vec2[2] ) / (sin(thetha) * RJR * RIR));

				mols->grads_X[ mols->ja[j] ] += ( df  * dtxi);
				mols->grads_Y[ mols->ja[j] ] += ( df  * dtyi);
				mols->grads_Z[ mols->ja[j] ] += ( df  * dtzi);

			}


		}

	}



/******************************************************************************
*
*                            TORSIONAL STUFF
*
******************************************************************************/

	torE = 0.0;

	if ( torflag == 1) {


		/* Improper terms */

		for ( i = 0; i < mols->n_impropers; ++i) {

			vec1[0] = (mols->x[mols->jp[i]] - mols->x[mols->kp[i]]);
			vec1[1] = (mols->y[mols->jp[i]] - mols->y[mols->kp[i]]);
			vec1[2] = (mols->z[mols->jp[i]] - mols->z[mols->kp[i]]);

			vec2[0] = (mols->x[mols->lp[i]] - mols->x[mols->kp[i]]);
			vec2[1] = (mols->y[mols->lp[i]] - mols->y[mols->kp[i]]);
			vec2[2] = (mols->z[mols->lp[i]] - mols->z[mols->kp[i]]);

			uvec1[0] = vec1[1] * vec2[2] - vec2[1] * vec1[2];
			uvec1[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
			uvec1[2] = vec1[0] * vec2[1] - vec2[0] * vec1[1];


			vec1[0] = (mols->x[mols->ip[i]] - mols->x[mols->kp[i]]);
			vec1[1] = (mols->y[mols->ip[i]] - mols->y[mols->kp[i]]);
			vec1[2] = (mols->z[mols->ip[i]] - mols->z[mols->kp[i]]);

			uvec2[0] = vec1[1] * vec2[2] - vec2[1] * vec1[2];
			uvec2[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
			uvec2[2] = vec1[0] * vec2[1] - vec2[0] * vec1[1];


			tor1 = uvec1[0] * uvec2[0] + uvec1[1] * uvec2[1] + uvec1[2] * uvec2[2];
			tor2 = sqrt(pow(uvec1[0], 2) + pow(uvec1[1], 2) + pow(uvec1[2], 2)) * sqrt( pow(uvec2[0], 2) + pow(uvec2[1], 2) + pow(uvec2[2], 2));
			tor3 = acos(tor1 / tor2);

			torE = torE + (impropers[mols->tp[i]][4] * pow(tor3 * 180 / G_PI - impropers[mols->tp[i]][5], 2));

		}


		for ( i = 0; i < mols->n_torsionals; ++i)
			tor_type[i] = 0;  /* Not assigned */


		for ( i = 0; i < mols->n_torsionals; ++i) {

			vec1[0] = (mols->x[mols->ik[i]] - mols->x[mols->jk[i]]);
			vec1[1] = (mols->y[mols->ik[i]] - mols->y[mols->jk[i]]);
			vec1[2] = (mols->z[mols->ik[i]] - mols->z[mols->jk[i]]);

			vec2[0] = (mols->x[mols->jk[i]] - mols->x[mols->kk[i]]);
			vec2[1] = (mols->y[mols->jk[i]] - mols->y[mols->kk[i]]);
			vec2[2] = (mols->z[mols->jk[i]] - mols->z[mols->kk[i]]);

			uvec1[0] = vec1[1] * vec2[2] - vec2[1] * vec1[2];
			uvec1[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
			uvec1[2] = vec1[0] * vec2[1] - vec2[0] * vec1[1];

			if ( gradflag == 1) {
				ra2 = sqrt(pow(vec2[0], 2) + pow(vec2[1], 2) + pow(vec2[2], 2));
				gaa = 0 - (1 / (pow(uvec1[0], 2) + pow(uvec1[1], 2) + pow(uvec1[2], 2))) * ra2;
				ra2 = 1 / ra2;
				fg = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
			}

			vec1[0] = (mols->x[mols->jk[i]] - mols->x[mols->kk[i]]);
			vec1[1] = (mols->y[mols->jk[i]] - mols->y[mols->kk[i]]);
			vec1[2] = (mols->z[mols->jk[i]] - mols->z[mols->kk[i]]);

			vec2[0] = (mols->x[mols->kk[i]] - mols->x[mols->lk[i]]);
			vec2[1] = (mols->y[mols->kk[i]] - mols->y[mols->lk[i]]);
			vec2[2] = (mols->z[mols->kk[i]] - mols->z[mols->lk[i]]);

			uvec2[0] = vec1[1] * vec2[2] - vec2[1] * vec1[2];
			uvec2[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
			uvec2[2] = vec1[0] * vec2[1] - vec2[0] * vec1[1];

			if ( gradflag == 1) {
				hg = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];

				rb2 = sqrt(pow(vec2[0], 2) + pow(vec2[1], 2) + pow(vec2[2], 2));
				gbb = (1 / (pow(uvec2[0], 2) + pow(uvec2[1], 2) + pow(uvec2[2], 2))) * rb2;
				rb2 = 1 / rb2;
			}
			tor1 = uvec1[0] * uvec2[0] + uvec1[1] * uvec2[1] + uvec1[2] * uvec2[2];
			tor2 = sqrt(pow(uvec1[0], 2) + pow(uvec1[1], 2) + pow(uvec1[2], 2)) * sqrt( pow(uvec2[0], 2) + pow(uvec2[1], 2) + pow(uvec2[2], 2));
			tor3 = acos(tor1 / tor2);
			if (tor1 > 0.0)
				tor3 = 0 - tor3;


			/* Check specific quartet */

			for ( j = 0; j < 90; ++j) {
				if ( ((torsionals2[j][1] == mols->gaff_types[mols->jk[i]] && torsionals2[j][2] == mols->gaff_types[mols->kk[i]]) || (torsionals2[j][2] == mols->gaff_types[mols->jk[i]] && torsionals2[j][1] == mols->gaff_types[mols->kk[i]])) && ( (torsionals2[j][0] == mols->gaff_types[mols->ik[i]] && torsionals2[j][3] == mols->gaff_types[mols->lk[i]]) || (torsionals2[j][3] == mols->gaff_types[mols->ik[i]] && torsionals2[j][0] == mols->gaff_types[mols->lk[i]])) ) {
					tor_type[i] = 1;
					pn = fabs(torsionals2[j][4]);
					pk = torsionals2[j][5];
					phi0 = (torsionals2[j][6] * G_PI / 180);

					torE = torE + ( (torsionals2[j][5] / torsionals2[j][7]) * (1.0 + cos( (fabs(torsionals2[j][4]) * tor3) - ( (torsionals2[j][6] * G_PI / 180)))));

					/* TODO. WARNING GRADIENT STUFF SHOULD BE HERE AND NOT IN THE NORMAL TORSIONAL BLOCK */


				}
			}


			if (tor_type[i] == 0) { /* Normal torsional */

				if ( mols->gaff_types[mols->jk[i]]  > mols->gaff_types[mols->kk[i]]) {
					t1 = mols->gaff_types[mols->kk[i]] - 1;
					t2 = mols->gaff_types[mols->jk[i]] - 1;
				}else{
					t2 = mols->gaff_types[mols->kk[i]] - 1;
					t1 = mols->gaff_types[mols->jk[i]] - 1;
				}


				if (  torsionals[t1][t2][0] != 0.0)
					torE = torE + ( (torsionals[t1][t2][2] / torsionals[t1][t2][0]) * (1.0 + cos( ( fabs(torsionals[t1][t2][1]) * tor3) - ( (torsionals[t1][t2][3] * G_PI / 180)))));
				if ( torsionals[t1][t2][4] != 0.0)
					torE = torE + ( (torsionals[t1][t2][6] / torsionals[t1][t2][4]) * (1.0 + cos( ( fabs(torsionals[t1][t2][5]) * tor3) - ( (torsionals[t1][t2][7] * G_PI / 180)))));
				if ( torsionals[t1][t2][8] != 0.0)
					torE = torE + ( (torsionals[t1][t2][10] / torsionals[t1][t2][8]) * (1.0 + cos( ( fabs(torsionals[t1][t2][9]) * tor3) - ( (torsionals[t1][t2][11] * G_PI / 180)))));
				if ( torsionals[t1][t2][12] != 0.0)
					torE = torE + ( (torsionals[t1][t2][14] / torsionals[t1][t2][12]) * (1.0 + cos( ( fabs(torsionals[t1][t2][13]) * tor3) - ( (torsionals[t1][t2][15] * G_PI / 180)))));


			}

			if ( gradflag == 1) {
				e1 = 1.0;
				df1 = 0.0;
				ddf1 = 0.0;
				if ( tor_type[i] == 0) {
					for (j = 0; j < torsionals[t1][t2][1]; ++j) {
						df1 = sin(tor3) + df1 * cos(tor3);
						ddf1 = e1 * cos(tor3) - df1 *  sin(tor3);
						e1 = ddf1;
					}

					e1 = e1 * cos( (torsionals[t1][t2][3] * G_PI / 180)) + df1 * sin( (torsionals[t1][t2][3] * G_PI / 180));
					df1 = df1 * cos( (torsionals[t1][t2][3] * G_PI / 180)) - ddf1 * sin( (torsionals[t1][t2][3] * G_PI / 180));
					df1 = (0 - torsionals[t1][t2][1]) * df1;
					ddf1 = ( 0 - (torsionals[t1][t2][1] * torsionals[t1][t2][1])) * e1;
					e1 = e1 + 1.0;
					dfx = gaa * uvec1[0] * df1 * torsionals[t1][t2][2];
					dfy = gaa * uvec1[1] * df1 * torsionals[t1][t2][2];
					dfz = gaa * uvec1[2] * df1 * torsionals[t1][t2][2];

					mols->grads_X[mols->ik[i]] += dfx;
					mols->grads_Y[mols->ik[i]] += dfy;
					mols->grads_Z[mols->ik[i]] += dfz;

					mols->grads_X[mols->jk[i]] -= dfx;
					mols->grads_Y[mols->jk[i]] -= dfy;
					mols->grads_Z[mols->jk[i]] -= dfz;


					dfx = gbb * uvec2[0] * df1 * torsionals[t1][t2][2];
					dfy = gbb * uvec2[0] * df1 * torsionals[t1][t2][2];
					dfz = gbb * uvec2[0] * df1 * torsionals[t1][t2][2];

					mols->grads_X[mols->lk[i]] += dfx;
					mols->grads_Y[mols->lk[i]] += dfy;
					mols->grads_Z[mols->lk[i]] += dfz;

					mols->grads_X[mols->kk[i]] -= dfx;
					mols->grads_Y[mols->kk[i]] -= dfy;
					mols->grads_Z[mols->kk[i]] -= dfz;


					fga = fg * (1 / (pow(uvec1[0], 2) + pow(uvec1[1], 2) + pow(uvec1[2], 2))) * sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
					hgb = hg * (1 / (pow(uvec2[0], 2) + pow(uvec2[1], 2) + pow(uvec2[2], 2))) * sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));

					dfx = (fga * uvec1[0] - hgb * uvec2[0]) * df1 * torsionals[t1][t2][2];
					dfy = (fga * uvec1[1] - hgb * uvec2[1]) * df1 * torsionals[t1][t2][2];
					dfz = (fga * uvec1[2] - hgb * uvec2[2]) * df1 * torsionals[t1][t2][2];

					mols->grads_X[mols->jk[i]] += dfx;
					mols->grads_Y[mols->jk[i]] += dfy;
					mols->grads_Z[mols->jk[i]] += dfz;
					mols->grads_X[mols->kk[i]] -= dfx;
					mols->grads_Y[mols->kk[i]] -= dfy;
					mols->grads_Z[mols->kk[i]] -= dfz;
				}else if ( tor_type[i] == 1 ) {


					for (j = 0; j < pn; ++j) {
						df1 = sin(tor3) + df1 * cos(tor3);
						ddf1 = e1 * cos(tor3) - df1 *  sin(tor3);
						e1 = ddf1;
					}


					e1 = e1 * cos( (phi0 * G_PI / 180)) + df1 * sin( (phi0 * G_PI / 180));
					df1 = df1 * cos( (phi0 * G_PI / 180)) - ddf1 * sin( (phi0 * G_PI / 180));
					df1 = (0 - pn) * df1;
					ddf1 = ( 0 - (pn * pn)) * e1;
					e1 = e1 + 1.0;

					dfx = gaa * uvec1[0] * df1 * pk;
					dfy = gaa * uvec1[1] * df1 * pk;
					dfz = gaa * uvec1[2] * df1 * pk;

					mols->grads_X[mols->ik[i]] += dfx;
					mols->grads_Y[mols->ik[i]] += dfy;
					mols->grads_Z[mols->ik[i]] += dfz;

					mols->grads_X[mols->jk[i]] -= dfx;
					mols->grads_Y[mols->jk[i]] -= dfy;
					mols->grads_Z[mols->jk[i]] -= dfz;


					dfx = gbb * uvec2[0] * df1 * pk;
					dfy = gbb * uvec2[0] * df1 * pk;
					dfz = gbb * uvec2[0] * df1 * pk;

					mols->grads_X[mols->lk[i]] += dfx;
					mols->grads_Y[mols->lk[i]] += dfy;
					mols->grads_Z[mols->lk[i]] += dfz;

					mols->grads_X[mols->kk[i]] -= dfx;
					mols->grads_Y[mols->kk[i]] -= dfy;
					mols->grads_Z[mols->kk[i]] -= dfz;


					fga = fg * (1 / (pow(uvec1[0], 2) + pow(uvec1[1], 2) + pow(uvec1[2], 2))) * sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
					hgb = hg * (1 / (pow(uvec2[0], 2) + pow(uvec2[1], 2) + pow(uvec2[2], 2))) * sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));

					dfx = (fga * uvec1[0] - hgb * uvec2[0]) * df1 * pk;
					dfy = (fga * uvec1[1] - hgb * uvec2[1]) * df1 * pk;
					dfz = (fga * uvec1[2] - hgb * uvec2[2]) * df1 * pk;


					mols->grads_X[mols->jk[i]] += dfx;
					mols->grads_Y[mols->jk[i]] += dfy;
					mols->grads_Z[mols->jk[i]] += dfz;
					mols->grads_X[mols->kk[i]] -= dfx;
					mols->grads_Y[mols->kk[i]] -= dfy;
					mols->grads_Z[mols->kk[i]] -= dfz;



				}

			}

		}
	}



/*****************************************************************************
 *
 *
 *                             VAN DER WAALS STUFF
 *
 *
 ****************************************************************************/

	vdwE = 0.0;
	if ( vdwflag == 1) {
		for ( j = 0; j < mols->n_pairs; ++j) {

			if ( flg != 0 || ( mols->backbone[mols->vdwpairs_a[j]] == 1 || mols->backbone[mols->vdwpairs_b[j]])) {
				dx = (mols->x[ mols->vdwpairs_a[j]] - mols->x[ mols->vdwpairs_b[j]]);
				dy = (mols->y[ mols->vdwpairs_a[j]] - mols->y[ mols->vdwpairs_b[j]]);
				dz = (mols->z[ mols->vdwpairs_a[j]] - mols->z[ mols->vdwpairs_b[j]]);

				dist = dx * dx + dy * dy + dz * dz;
				dist2 = dist;
				dist = sqrt(dist);
				distr = 1 / dist;

				tmpE = (332.0 * mols->pcharges[ mols->vdwpairs_a[j] ] * mols->pcharges[ mols->vdwpairs_b[j] ]) * distr;
				vdwE += tmpE;


				epsilon = prevdw[ mols->gaff_types[ mols->vdwpairs_a[j]] - 1][mols->gaff_types[ mols->vdwpairs_b[j]] - 1 ][1];
				sigma6 = prevdw[ mols->gaff_types[ mols->vdwpairs_a[j]] - 1][mols->gaff_types[ mols->vdwpairs_b[j]] - 1 ][2];

				r6 = sigma6 / (dist2 * dist2 * dist2);
				r12 = r6 * r6;

				epsilon /= mols->vdw_type[j];

				tmpE = epsilon * (r12 - 2.0 * r6);

				vdwE = vdwE + tmpE;

				if ( gradflag == 1 ) {
					r6 = r6 * (distr);      /* 7*/
					r12 = r12 * (distr);    /* 13*/
					df = (12.0 * epsilon) * ( r6 - r12);
					df += 0 - ((332.0 * mols->pcharges[ mols->vdwpairs_a[j] ] * mols->pcharges[ mols->vdwpairs_b[j] ]) / (dist2));

					dx *= df;
					dy *= df;
					dz *= df;
					mols->grads_X[ mols->vdwpairs_a[j] ] += dx;
					mols->grads_X[ mols->vdwpairs_b[j] ] -= dx;

					mols->grads_Y[ mols->vdwpairs_a[j] ] += dy;
					mols->grads_Y[ mols->vdwpairs_b[j] ] -= dy;

					mols->grads_Z[ mols->vdwpairs_a[j] ] += dz;
					mols->grads_Z[ mols->vdwpairs_b[j] ] -= dz;


				}

			}
		}

	}

	totalE = vdwE + torE + bondE + angleE;
	if ( verboseflag == 1)
		fprintf(stderr, "ENERGY ROUTINE> TERMS. VDW: %f. TOR: %f. BOND: %f. ANGLE: %f.\n", vdwE, torE, bondE, angleE);



	if ( isinf(totalE) != 0 || isnan(totalE) == 0 || totalE > 9999.9) {
		*tenergy = totalE;
		retval = 0;
	}else{
		*tenergy = 9999.9;
		retval = -1;
	}

	if ( *tenergy > 9999.9)
		*tenergy = 9999.9;

	*mymol = mols;
	return retval;
}


/**
 * @brief Interpolate non-bonded force field energy from grids with method 2
 *
 * @param mol Structure with the ligand
 * @param grids Grids with potentials
 * @param types Atom types for the force field grids of the ligand
 * @param min_grid minimum grid coordinates
 * @param max_grid maximum grid coordinates
 * @param nx vector for closest point for atom
 * @param ny vector for closest point for atom
 * @param nz vector for closest point for atom
 * @param max_x number of points of x axis in the grid
 * @param max_y number of points of y axis in the grid
 * @param max_z number of points of z axis in the grid
 * @param spacing grid spacing
 * @return Interpolated energy of the ligand in the grid
 *
 */
float get_ff_grid_energy_ipol2(MOL2 *mol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing)
{

	int i = 0, real_atoms = 0, out_flag = 0, l = 0;
	float energy = 0.0f;
	int ndiag[8][3];
	int icube = 0;
	int kpot = 0, ix = 0, iy = 0, iz = 0;
	float xj = 0.0f, yj = 0.0f, zj = 0.0f, dis = 0.0f;
	float q = 0.0f, xa = 0.0f, ya = 0.0f, za = 0.0f;
	int m = 0, n = 0, ia = 0, ib = 0, ic = 0;
	float hdist = 0.0f, x = 0.0f, y = 0.0f, z = 0.0f;
	float a = 0, b = 0, c = 0, wi = 0.0f;
	float xd, yd, zd, i1, i2, j1, j2, w1, w2, e1, e2;
	int pi = 0, pj = 0, pk = 0;
	float dist = 0.0f, dist2 = 0.0f, enorm = 0.0f, enorm2 = 0.0f;
	float V000, V001, V011, V010, V110, V111, V101, V100;

	for ( i = 0; i < mol->n_atoms; ++i)
		++real_atoms;

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




	if ( out_flag == 0) {
		energy = 0.0f;
		enorm = 0.0f;
		for ( i = 0; i < real_atoms; ++i) {

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


			if ( types[i] != -1 && nx[i] >= 0 && ny[i] >= 0 &&
			     nz[i] >= 0 && nx[i] < max_x && ny[i] < max_y &&
			     nz[i] < max_z) {
				enorm = enorm +
					grids[types[i]][nx[i]][ny[i]][nz[i]] +
					mol->pcharges[i] *
					grids[8][nx[i]][ny[i]][nz[i]];
			}



			if ( types[i] != -1 && ix >= 1 && iy >= 1 &&
			     iz >= 1 && ix < max_x - 1 && iy < max_y - 1 &&
			     iz < max_z - 1 ) {

				pi = pj = pj = 1;
				xa = (ix * spacing) + min_grid[0];
				ya = (iy * spacing) + min_grid[1];
				za = (iz * spacing) + min_grid[2];

				xd = x - xa;
				yd = y - ya;
				zd = z - za;

				xa = ((ix + 1) * spacing) + min_grid[0];
				ya = (iy * spacing) + min_grid[1];
				za = (iz * spacing) + min_grid[2];

				dist = (x - xa) * (x - xa) + (y - ya) * (y - ya) + (z - za) * (z - za);

				xa = ((ix - 1) * spacing) + min_grid[0];
				ya = (iy * spacing) + min_grid[1];
				za = (iz * spacing) + min_grid[2];

				dist2 = (x - xa) * (x - xa) + (y - ya) * (y - ya) + (z - za) * (z - za);

				if ( dist > dist2)
					pi = -1;


				xa = (ix * spacing) + min_grid[0];
				ya = ((iy + 1) * spacing) + min_grid[1];
				za = (iz * spacing) + min_grid[2];

				dist = (x - xa) * (x - xa) + (y - ya) * (y - ya) + (z - za) * (z - za);

				xa = (ix * spacing) + min_grid[0];
				ya = ((iy - 1) * spacing) + min_grid[1];
				za = (iz * spacing) + min_grid[2];

				dist2 = (x - xa) * (x - xa) + (y - ya) * (y - ya) + (z - za) * (z - za);

				if ( dist > dist2)
					pj = -1;

				xa = (ix * spacing) + min_grid[0];
				ya = (iy * spacing) + min_grid[1];
				za = ((iz + 1) * spacing) + min_grid[2];

				dist = (x - xa) * (x - xa) + (y - ya) * (y - ya) + (z - za) * (z - za);

				xa = (ix * spacing) + min_grid[0];
				ya = (iy * spacing) + min_grid[1];
				za = ((iz - 1) * spacing) + min_grid[2];

				dist2 = (x - xa) * (x - xa) + (y - ya) * (y - ya) + (z - za) * (z - za);

				if ( dist > dist2)
					pk = -1;

				i1 = grids[kpot][ix][iy][iz] * (1.0 - zd) + (grids[kpot][ix][iy][iz + pk] * zd);
				i2 = grids[kpot][ix][iy + pj][iz] * (1.0 - zd) + (grids[kpot][ix][iy + pj][iz + pk] * zd);
				j1 = grids[kpot][ix + pi][iy][iz] * (1.0 - zd) + (grids[kpot][ix + pi][iy][iz + pk] * zd);
				j2 = grids[kpot][ix + pi][iy + pj][iz] * (1.0 - zd) + (grids[kpot][ix + pi][iy + pj][iz + pk] * zd);
				w1 = i1 * ( 1 - yd) + i2 * yd;
				w2 = j1 * (1 - yd) + j2 * yd;
				e1 = w1 * (1 - xd) + w2 * xd;
				kpot = 8;

				i1 = grids[kpot][ix][iy][iz] * q * (1.0 - zd) + (grids[kpot][ix][iy][iz + pk] * q * zd);
				i2 = grids[kpot][ix][iy + pj][iz] * q * (1.0 - zd) + (grids[kpot][ix][iy + pj][iz + pk] * q * zd);
				j1 = grids[kpot][ix + pi][iy][iz] * q * (1.0 - zd) + (grids[kpot][ix + pi][iy][iz + pk] * q * zd);
				j2 = grids[kpot][ix + pi][iy + pj][iz] * q * (1.0 - zd) + (grids[kpot][ix + pi][iy + pj][iz + pk] * q * zd);

				w1 = i1 * ( 1 - yd) + i2 * yd;
				w2 = j1 * (1 - yd) + j2 * yd;
				e2 = w1 * (1 - xd) + w2 * xd;

				enorm2 += e1 + e2;

			}


		}


	}else
		energy = 9999.9f;



	return energy;
}


/**
 * @brief Get closest ChemScore energy from grids
 *
 * @param mol Structure with the ligand
 * @param grids Grids with potentials
 * @param types Atom types for ChemScore grids of the ligand
 * @param min_grid minimum grid coordinates
 * @param max_grid maximum grid coordinates
 * @param nx vector for closest point for atom
 * @param ny vector for closest point for atom
 * @param nz vector for closest point for atom
 * @param max_x number of points of x axis in the grid
 * @param max_y number of points of y axis in the grid
 * @param max_z number of points of z axis in the grid
 * @param spacing grid spacing
 * @return Energy of the ligand in the grid
 *
 */

float get_chemscore_grid_energy(MOL2 *mol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing)
{

	int i = 0, real_atoms = 0, out_flag = 0, l = 0;
	float energy = 0.0f;

	for ( i = 0; i < mol->n_atoms; ++i)
		++real_atoms;

	/* Check that the ligand is inside the grid */
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
	/* If the ligand is inside, compute the nearest grid point to each atom */
	if ( out_flag == 0) {
		for ( i = 0; l < real_atoms; ++l) {
			nx[l] = (int)roundf((mol->x[l] - min_grid[0]) / spacing);
			ny[l] = (int)roundf((mol->y[l] - min_grid[1]) / spacing);
			nz[l] = (int)roundf((mol->z[l] - min_grid[2]) / spacing);
		}
		energy = 0.0f;
		/* Sum contribution from all ligand atoms */
		for ( l = 0; l < real_atoms; ++l) {
			if ( types[l] != -1 && nx[l] >= 0 && ny[l] >= 0 &&
			     nz[l] >= 0 && nx[l] < max_x && ny[l] < max_y &&
			     nz[l] < max_z) {
				/* Each contribution + clash term */
				energy = energy +
					 grids[types[l]][nx[l]][ny[l]][nz[l]] +
					 grids[6][nx[l]][ny[l]][nz[l]];

				/* If acceptor or mixed, also metals contribution */
				if ( types[l] == 3 || types[l] == 4)
					energy = energy +
						 grids[5][nx[l]][ny[l]][nz[l]];
			}
		}

	}else
		energy = 9999.9;

	return energy;
}



/**
 * @brief Get closest non-bonded force field energy from grids 
 *
 * @param mol Structure with the ligand
 * @param grids Grids with potentials
 * @param types Atom types for the force field grids of the ligand
 * @param min_grid minimum grid coordinates
 * @param max_grid maximum grid coordinates
 * @param nx vector for closest point for atom
 * @param ny vector for closest point for atom
 * @param nz vector for closest point for atom
 * @param max_x number of points of x axis in the grid
 * @param max_y number of points of y axis in the grid
 * @param max_z number of points of z axis in the grid
 * @param spacing grid spacing
 * @return Energy of the ligand in the grid
 *
 */

float get_ff_grid_energy(MOL2 *mol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing)
{

	int i = 0, real_atoms = 0, out_flag = 0, l = 0;
	float energy = 0.0f;

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
			nx[l] = (int)roundf((mol->x[l] - min_grid[0]) / spacing);
			ny[l] = (int)roundf((mol->y[l] - min_grid[1]) / spacing);
			nz[l] = (int)roundf((mol->z[l] - min_grid[2]) / spacing);
		}

		energy = 0.0f;


		for ( l = 0; l < real_atoms; ++l) {

			if ( types[l] != -1 && nx[l] >= 0 && ny[l] >= 0 &&
			     nz[l] >= 0 && nx[l] < max_x && ny[l] < max_y &&
			     nz[l] < max_z) {
				energy = energy +
					 grids[types[l]][nx[l]][ny[l]][nz[l]] +
					 mol->pcharges[l] *
					 grids[8][nx[l]][ny[l]][nz[l]];
			}
		}


	}else
		energy = 9999.9;

	return energy;
}


/**
 * @brief Interpolate non-bonded force field energy from grids with method 1
 *
 * @param mol Structure with the ligand
 * @param grids Grids with potentials
 * @param types Atom types for the force field grids of the ligand
 * @param min_grid minimum grid coordinates
 * @param max_grid maximum grid coordinates
 * @param nx vector for closest point for atom
 * @param ny vector for closest point for atom
 * @param nz vector for closest point for atom
 * @param max_x number of points of x axis in the grid
 * @param max_y number of points of y axis in the grid
 * @param max_z number of points of z axis in the grid
 * @param spacing grid spacing
 * @return Interpolated energy of the ligand in the grid
 *
 */
float get_ff_grid_energy_ipol(MOL2 *mol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing)
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
	float dist = 0.0f, dist2 = 0.0f, enorm2 = 0.0f;
	float V000, V100, V010, V110, V001, V101, V011, V111;
	register double u,   v,   w;
	register double p0u, p0v, p0w;
	register double p1u, p1v, p1w;
	register int u0,  v0,  w0;
	register int u1,  v1,  w1;

	for ( i = 0; i < mol->n_atoms; ++i)
		++real_atoms;

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




	if ( out_flag == 0) {
		energy = 0.0f;
		e1 = e2 = 0.0f;
		for ( i = 0; i < real_atoms; ++i) {
			e1 = e2 = 0.0f;


			kpot = types[i];
			x = mol->x[i];
			y = mol->y[i];
			z = mol->z[i];

			nx[i] = (int)roundf((mol->x[i] - min_grid[0]) / spacing);
			ny[i] = (int)roundf((mol->y[i] - min_grid[1]) / spacing);
			nz[i] = (int)roundf((mol->z[i] - min_grid[2]) / spacing);

			ix = nx[i];
			iy = ny[i];
			iz = nz[i];
			q = mol->pcharges[i];


			if ( !(ix >= 1 && iy >= 1 &&
			       iz >= 1 && ix < max_x - 1 && iy < max_y - 1 &&
			       iz < max_z - 1)) {
				energy = 9999.9f;
				return energy;
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


				e1 = p1u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
				e2 = p1u * p1v * p1w * grids[8][ w0 ][ v0 ][ u0 ];

				e1 += p0u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u1 ];
				e2 += p0u * p1v * p1w * grids[8][ w0 ][ v0 ][ u1 ];

				e1 += p1u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u0 ];
				e2 += p1u * p0v * p1w * grids[8][ w0 ][ v1 ][ u0 ];

				e1 += p0u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u1 ];
				e2 += p0u * p0v * p1w * grids[8][ w0 ][ v1 ][ u1 ];

				e1 += p1u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
				e2 += p1u * p1v * p0w * grids[8][ w1 ][ v0 ][ u0 ];

				e1 += p0u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u1 ];
				e2 += p0u * p1v * p0w * grids[8][ w1 ][ v0 ][ u1 ];

				e1 += p1u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u0 ];
				e2 += p1u * p0v * p0w * grids[8][ w1 ][ v1 ][ u0 ];

				e1 += p0u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u1 ];
				e2 += p0u * p0v * p0w * grids[8][ w1 ][ v1 ][ u1 ];

				energy = energy + e1; /* + (q * e2);*/

			}


		}


	}else
		energy = 9999.9f;


	return energy;
}

/**
 * @brief Interpolate ChemScore energy from grids with method 1
 *
 * @param mol Structure with the ligand
 * @param grids Grids with potentials
 * @param types Atom types for ChemScore grids of the ligand
 * @param min_grid minimum grid coordinates
 * @param max_grid maximum grid coordinates
 * @param nx vector for closest point for atom
 * @param ny vector for closest point for atom
 * @param nz vector for closest point for atom
 * @param max_x number of points of x axis in the grid
 * @param max_y number of points of y axis in the grid
 * @param max_z number of points of z axis in the grid
 * @param spacing grid spacing
 * @return Interpolated energy of the ligand in the grid
 *
 */

float get_chemscore_grid_energy_ipol(MOL2 *mol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing)
{

	int i = 0, real_atoms = 0, out_flag = 0, l = 0;
	float energy = 0.0f;
	register double u,   v,   w;
	register double p0u, p0v, p0w;
	register double p1u, p1v, p1w;
	register int u0,  v0,  w0;
	register int u1,  v1,  w1;
	float e1, e2, e3, x, y, z, xd, yd, zd;
	int kpot = 0;

	/* Count only ligand atoms. Exclude all possible protein atoms inside mol */
	for ( i = 0; i < mol->n_atoms; ++i)
		++real_atoms;

	/* Check that the ligand is inside the grid */
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
	/* If the ligand is inside, compute the nearest grid point to each atom */
	if ( out_flag == 0) {
		for ( i = 0; l < real_atoms; ++l) {
			nx[l] = (int)roundf((mol->x[l] - min_grid[0]) / spacing);
			ny[l] = (int)roundf((mol->y[l] - min_grid[1]) / spacing);
			nz[l] = (int)roundf((mol->z[l] - min_grid[2]) / spacing);
		}
		energy = 0.0f;
		/* Sum contribution from all ligand atoms */
		for ( l = 0; l < real_atoms; ++l) {
			kpot = types[l];
			if ( types[l] != -1 && nx[l] >= 1 && ny[l] >= 1 &&
			     nz[l] >= 1 && nx[l] < max_x - 1 && ny[l] < max_y - 1 &&
			     nz[l] < max_z - 1) {

				x = mol->x[l];
				y = mol->y[l];
				z = mol->z[l];

				xd = min_grid[0];
				yd = min_grid[1];
				zd = min_grid[2];

				w1  = (w0 = (int)(w = (x - xd) / spacing )) + 1;
				p1w = 1.0L - (p0w = w - (float)w0);

				v1  = (v0 = (int)(v = (y - yd) / spacing)) + 1;
				p1v = 1.0L - (p0v = v - (float)v0);

				u1  = (u0 = (int)(u = (z - zd) / spacing)) + 1;
				p1u = 1.0L - (p0u = u - (float)u0);


				e1 = p1u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
				e2 = p1u * p1v * p1w * grids[6][ w0 ][ v0 ][ u0 ];
				e3 = p1u * p1v * p1w * grids[5][ w0 ][ v0 ][ u0 ];


				e1 += p0u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u1 ];
				e2 += p0u * p1v * p1w * grids[6][ w0 ][ v0 ][ u1 ];
				e3 += p0u * p1v * p1w * grids[5][ w0 ][ v0 ][ u1 ];

				e1 += p1u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u0 ];
				e2 += p1u * p0v * p1w * grids[6][ w0 ][ v1 ][ u0 ];
				e3 += p1u * p0v * p1w * grids[5][ w0 ][ v1 ][ u0 ];


				e1 += p0u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u1 ];
				e2 += p0u * p0v * p1w * grids[6][ w0 ][ v1 ][ u1 ];
				e3 += p0u * p0v * p1w * grids[5][ w0 ][ v1 ][ u1 ];

				e1 += p1u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
				e2 += p1u * p1v * p0w * grids[6][ w1 ][ v0 ][ u0 ];
				e3 += p1u * p1v * p0w * grids[5][ w1 ][ v0 ][ u0 ];

				e1 += p0u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u1 ];
				e2 += p0u * p1v * p0w * grids[6][ w1 ][ v0 ][ u1 ];
				e3 += p0u * p1v * p0w * grids[5][ w1 ][ v0 ][ u1 ];

				e1 += p1u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u0 ];
				e2 += p1u * p0v * p0w * grids[6][ w1 ][ v1 ][ u0 ];
				e3 += p1u * p0v * p0w * grids[5][ w1 ][ v1 ][ u0 ];

				e1 += p0u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u1 ];
				e2 += p0u * p0v * p0w * grids[6][ w1 ][ v1 ][ u1 ];
				e3 += p0u * p0v * p0w * grids[5][ w1 ][ v1 ][ u1 ];


				/* Each contribution + clash term */

				energy = energy + e1;
/*				if( kpot != 3 && kpot != 2 && kpot != 4)
				{*/
					energy = energy + e2;
/*				}*/

				/* If acceptor or mixed, also metals contribution */
				if ( types[l] == 3 || types[l] == 4)
					energy = energy + e3;



			}
		}

	}else
		energy = 9999.9;

	return energy;
}

/**
 *
 *	@brief Calculate Cross product of two vectors
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param a Vector i
 *	@param b Vector j
 *	@param c Resulting vector
 *	@return null
 *
 */
void cross_product(float a[3], float b[3], float c[3])
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

/**
 *
 *	@brief Conver invidual force contributions on atoms to global rigid body forces (translation + rotation)
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 with the molecule and the gradients
 *	@param force Resulting force
 *	@param torque Resulting torque
 *	@return null
 *
 */
void grad_to_rigid(MOL2 *mol, float force[3], float torque[3])
{
	int i = 0;
	float com[3];
	float tmp[3];
	float tmp2[3];
	float tmp3[3];
	float norm = 0;

	torque[0] = 0.0f;
	torque[1] = 0.0f;
	torque[2] = 0.0f;
	force[0] = 0.0f;
	force[1] = 0.0f;
	force[2] = 0.0f;

	com[0] = com[1] = com[2] = 0.0f;
	for ( i = 0; i < mol->n_atoms; ++i) {
		com[0] += mol->x[i];
		com[1] += mol->y[i];
		com[2] += mol->z[i];
	}

	com[0] /= (float)mol->n_atoms;
	com[1] /= (float)mol->n_atoms;
	com[2] /= (float)mol->n_atoms;


	for (i = 0; i < mol->n_atoms; ++i) {
		/* Average forces in the center of mass. Ft = E Fi */
		force[0] -= mol->grads_X[i]; 
		force[1] -= mol->grads_Y[i];
		force[2] -= mol->grads_Z[i];

		tmp2[0] = mol->x[i] - com[0];
		tmp2[1] = mol->y[i] - com[1];
		tmp2[2] = mol->z[i] - com[2];

		tmp3[0] = mol->grads_X[i];
		tmp3[1] = mol->grads_Y[i];
		tmp3[2] = mol->grads_Z[i];

		/* Momentum or torque for rotation. M = rxL */
		cross_product(tmp2, tmp3, tmp);

		torque[0] -= tmp[0];
		torque[1] -= tmp[1];
		torque[2] -= tmp[2];
	}

	return;
}

/**
 *
 *	@brief Apply a rigid body transformation using an alpha parameter to a molecule
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mymol MOL2 with the molecule to transform
 *	@param trans vector of translation
 *	@param rot vector of rotation
 *	@param alfa scalar for multiplying the transformation
 *	@return null
 *
 */
void transformate_rigid(MOL2 **mymol, float *trans, float *rot, float alfa)
{

	int i = 0;
	MOL2 *mol = NULL;
	float com[3];
	float vec[3];
	float tmp[3];
	float tmp2[3];
	float norm = 0.0f;
	float a, b, c, d;
	float aa, ab, ac, ad, bb, bc, bd, cc, cd, dd;
	float rotM[3][3];

	mol = *mymol;

	com[0] = com[1] = com[2] = 0.0f;

	for ( i = 0; i < mol->n_atoms; ++i) {
		com[0] += mol->x[i];
		com[1] += mol->y[i];
		com[2] += mol->z[i];
	}

	com[0] /= (float)mol->n_atoms;
	com[1] /= (float)mol->n_atoms;
	com[2] /= (float)mol->n_atoms;

	rot[0] *= alfa;
	rot[1] *= alfa;
	rot[2] *= alfa;

	norm = sqrt(rot[0] * rot[0] + rot[1] * rot[1] + rot[2] * rot[2]);
	rot[0] /= norm;
	rot[1] /= norm;
	rot[2] /= norm;

/* REF: http://content.gpwiki.org/index.php/OpenGL:Tutorials:Using_Quaternions_to_represent_rotation#Quaternion_from_axis-angle */

	a = cos(norm / 2);
	b = sin(norm / 2) * rot[0];
	c = sin(norm / 2) * rot[1];
	d = sin(norm / 2) * rot[2];



/*          a = 1;
          b = rot[0];
          c = rot[1];
         d = rot[2];*/
/*	a = 1; b = 0; c = 0; d = 0;*/

/* REF: http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Conversion_to_and_from_the_matrix_representation */
	aa = a * a;
	ab = a * b;
	ac = a * c;
	ad = a * d;
	bb = b * b;
	bc = b * c;
	bd = b * d;
	cc = c * c;
	cd = c * d;
	dd = d * d;

	rotM[0][0] = (aa + bb - cc - dd);
	rotM[0][1] = 2 * (-ad + bc);
	rotM[0][2] = 2 * (ac + bd);

	rotM[1][0] = 2 * (ad + bc);
	rotM[1][1] = (aa - bb + cc - dd);
	rotM[1][2] = 2 * (-ab + cd);

	rotM[2][0] = 2 * (-ac + bd);
	rotM[2][1] = 2 * (ab + cd);
	rotM[2][2] = (aa - bb - cc + dd);



	for ( i = 0; i < mol->n_atoms; ++i) {
		mol->x[i] = mol->x[i] - com[0];
		mol->y[i] = mol->y[i] - com[1];
		mol->z[i] = mol->z[i] - com[2];

		tmp[0] = mol->x[i] * rotM[0][0] + mol->y[i] * rotM[1][0]  + mol->z[i] * rotM[2][0];
		tmp[1] = mol->x[i] * rotM[0][1] + mol->y[i] * rotM[1][1]  + mol->z[i] * rotM[2][1];
		tmp[2] = mol->x[i] * rotM[0][2] + mol->y[i] * rotM[1][2]  + mol->z[i] * rotM[2][2];


		mol->x[i] = tmp[0] + com[0] + trans[0] * alfa;
		mol->y[i] = tmp[1] + com[1] + trans[1] * alfa;
		mol->z[i] = tmp[2] + com[2] + trans[2] * alfa;
	}


}


/**
 *	@brief Evaluate a one side block function
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param x actual distance
 *	@param R1 Ideal distance
 *	@param R2 Maximum distance
 *	@return block function coefficient
 *
 */
float eval_block(float x, float R1, float R2)
{
	float retval = 0.0f;

	if ( x <= R1)
		retval = 1.0f;
	else if (  x > R1 && x < R2)
		retval = 1.0f - ((x - R1) / (R2 - R1));

	return retval;
}

/**
 *
 *	@brief Get hydrogen bond interactions energy
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 with the ligand
 *	@param prot MOL2 with the protein
 *	@param min_grid Minimum coordinates of the grid to consider an interaction
 *	@param max_grid Maximum coordinates of the grid to consider an interaction
 *	@return Energy of hydrogen bonds
 *
 */
float get_hbonds_energy(MOL2 *mol, MOL2 *prot, float min_grid[3], float max_grid[3])
{

	float dist = 0.0f, dx = 0.0f, dy = 0.0f, dz = 0.0f;
	int i = 0, j = 0, k = 0;
	int ix1 = 0, ix2 = 0;
	float alfa = 0.0, beta = 0.0, tmp_beta = 0.0f;
	int prot_flag  = 0, lig_flag = 0;
	float h_atom[3], D[3], A[3], X[3];
	float ha[3], XA[3], DH[3], AH[3];
	float norm1 = 0.0f, norm2 = 0.0f, norm3 = 0.0f, norm4 = 0.0f;
	int p_type = 0, l_type = 0;

	float score = 0.0f, tmpbeta = 0.0f;

	for ( j = 0; j < prot->n_atoms; ++j ) {
		/*
		 *
		 *	Just HN, HO, O, OH, N1, N2, NA, NB, NC, ND, NE, NF, NH
		 *
		 *
		 */
		if ( prot->gaff_types[j] == HN ||  prot->gaff_types[j] == HO ||  (prot->gaff_types[j] >= 14 && prot->gaff_types[j] <= 15 ) || (prot->gaff_types[j] >= 37 && prot->gaff_types[j] <= 39 ) || (prot->gaff_types[j] >= 41 && prot->gaff_types[j] <= 47 ) ) {
			/* Inside the grid */
			if ( prot->x[j] >= min_grid[0] && prot->y[j] >= min_grid[1] && prot->z[j] >= min_grid[2] && prot->x[j] <= max_grid[0] && prot->y[j] <= max_grid[1] && prot->z[j] <= max_grid[2]) {
				for ( i = 0; i < mol->n_atoms; i++) {
					if ( mol->gaff_types[i] == HN ||  mol->gaff_types[i] == HO ||  (mol->gaff_types[i] >= 14 && mol->gaff_types[i] <= 15 ) || (mol->gaff_types[i] >= 37 && mol->gaff_types[i] <= 39 ) || (mol->gaff_types[i] >= 41 && mol->gaff_types[i] <= 47 ) ) {


						dx = prot->x[j] - mol->x[i];
						dy = prot->y[j] - mol->y[i];
						dz = prot->z[j] - mol->z[i];
						p_type = prot->gaff_types[j];
						l_type = mol->gaff_types[i];
						dist = dx * dx + dy * dy + dz * dz;
						dist = sqrt(dist);
						prot_flag = -1;
						lig_flag = -1;
						alfa = 0.0f;
						beta = 0.0f;
						h_atom[0] = h_atom[1] = h_atom[2] = 0.0f;
						D[0] = D[1] = D[2] = 0.0f;
						A[0] = A[1] = A[2] = 0.0f;
						X[0] = X[1] = X[2] = 0.0f;

						if ( dist <= 3.0) {

							/* Protein is the donnor */
							if ( prot->gaff_types[j] == HN ||  prot->gaff_types[j] == HO ) {
								prot_flag = 0;
								h_atom[0] = prot->x[j];
								h_atom[1] = prot->y[j];
								h_atom[2] = prot->z[j];
								for ( k = 0; k < prot->n_bonds; k++) {
									if ( (prot->bond_a1[k] - 1) == j) {
										D[0] = prot->x[ (prot->bond_a2[k] - 1)];
										D[1] = prot->y[ (prot->bond_a2[k] - 1)];
										D[2] = prot->z[ (prot->bond_a2[k] - 1)];
									}else if ((prot->bond_a2[k] - 1) == j) {
										D[0] = prot->x[ (prot->bond_a1[k] - 1)];
										D[1] = prot->y[ (prot->bond_a1[k] - 1)];
										D[2] = prot->z[ (prot->bond_a1[k] - 1)];
									}
								}
								/* Ligand is the donnor */
							}else if ( mol->gaff_types[i] == HN ||  mol->gaff_types[i] == HO ) {
								lig_flag = 0;
								h_atom[0] = mol->x[i];
								h_atom[1] = mol->y[i];
								h_atom[2] = mol->z[i];
								for ( k = 0; k < mol->n_bonds; k++) {
									if ( (mol->bond_a1[k] - 1) == i) {
										D[0] = mol->x[ (mol->bond_a2[k] - 1)];
										D[1] = mol->y[ (mol->bond_a2[k] - 1)];
										D[2] = mol->z[ (mol->bond_a2[k] - 1)];
									}else if ((mol->bond_a2[k] - 1) == i) {
										D[0] = mol->x[ (mol->bond_a1[k] - 1)];
										D[1] = mol->y[ (mol->bond_a1[k] - 1)];
										D[2] = mol->z[ (mol->bond_a1[k] - 1)];
									}
								}

							}


							/* Protein is acceptor */
							if ( (prot->gaff_types[j] >= 14 && prot->gaff_types[j] <= 15 ) || (prot->gaff_types[j] >= 37 &&
							   prot->gaff_types[j] <= 39 ) || (prot->gaff_types[j] >= 41 && prot->gaff_types[j] <= 47 ) ) {
								prot_flag = 1;
								A[0] = prot->x[j];
								A[1] = prot->y[j];
								A[2] = prot->z[j];
								for ( k = 0; k < prot->n_bonds; k++) {
									if ( (prot->bond_a1[k] - 1) == j) {
										X[0] = prot->x[ (prot->bond_a2[k] - 1)];
										X[1] = prot->y[ (prot->bond_a2[k] - 1)];
										X[2] = prot->z[ (prot->bond_a2[k] - 1)];
									}else if ((prot->bond_a2[k] - 1) == j) {
										X[0] = prot->x[ (prot->bond_a1[k] - 1)];
										X[1] = prot->y[ (prot->bond_a1[k] - 1)];
										X[2] = prot->z[ (prot->bond_a1[k] - 1)];
									}

									XA[0] = X[0] - A[0];
									XA[1] = X[1] - A[1];
									XA[2] = X[2] - A[2];
									ha[0] = h_atom[0] - A[0];
									ha[1] = h_atom[1] - A[1];
									ha[2] = h_atom[2] - A[2];
									norm1 = XA[0] * XA[0] + XA[1] * XA[1] + XA[2] * XA[2];
									norm1 = sqrt(norm1);
									norm2 = ha[0] * ha[0] + ha[1] * ha[1] + ha[2] * ha[2];
									norm2 = sqrt(norm2);
									norm3 = XA[0] * ha[0] + XA[1] * ha[1] + XA[2] * ha[2];
									tmpbeta = acos( norm3 / (norm1 * norm2));
									if ( (2.35 - fabs(tmpbeta)) < (2.35 - fabs(beta)))
										beta = tmpbeta;

								}



								/* Ligand is the acceptor */
							}else if ( (mol->gaff_types[i] >= 14 && mol->gaff_types[i] <= 15 ) || (mol->gaff_types[i] >= 37 &&
							       mol->gaff_types[i] <= 39 ) || (mol->gaff_types[i] >= 41 && mol->gaff_types[i] <= 47 ) ) {
								lig_flag = 1;
								A[0] = mol->x[i];
								A[1] = mol->y[i];
								A[2] = mol->z[i];
								for ( k = 0; k < mol->n_bonds; k++) {
									if ( (mol->bond_a1[k] - 1) == i) {
										X[0] = mol->x[ (mol->bond_a2[k] - 1)];
										X[1] = mol->y[ (mol->bond_a2[k] - 1)];
										X[2] = mol->z[ (mol->bond_a2[k] - 1)];
									}else if ((mol->bond_a2[k] - 1) == i) {
										X[0] = mol->x[ (mol->bond_a1[k] - 1)];
										X[1] = mol->y[ (mol->bond_a1[k] - 1)];
										X[2] = mol->z[ (mol->bond_a1[k] - 1)];
									}

									XA[0] = X[0] - A[0];
									XA[1] = X[1] - A[1];
									XA[2] = X[2] - A[2];
									ha[0] = h_atom[0] - A[0];
									ha[1] = h_atom[1] - A[1];
									ha[2] = h_atom[2] - A[2];
									norm1 = XA[0] * XA[0] + XA[1] * XA[1] + XA[2] * XA[2];
									norm1 = sqrt(norm1);
									norm2 = ha[0] * ha[0] + ha[1] * ha[1] + ha[2] * ha[2];
									norm2 = sqrt(norm2);
									norm3 = XA[0] * ha[0] + XA[1] * ha[1] + XA[2] * ha[2];
									tmpbeta = acos( norm3 / (norm1 * norm2));
									if ( (3.14159 - tmpbeta) < (3.14159 - beta))
										beta = tmpbeta;

								}

							}

							/* Alfa angle */

							DH[0] = D[0] - h_atom[0];
							DH[1] = D[1] - h_atom[1];
							DH[2] = D[2] - h_atom[2];
							AH[0] = A[0] - h_atom[0];
							AH[1] = A[1] - h_atom[1];
							AH[2] = A[2] - h_atom[2];
							norm1 = DH[0] * DH[0] + DH[1] * DH[1] + DH[2] * DH[2];
							norm1 = sqrt(norm1);
							norm2 = AH[0] * AH[0] + AH[1] * AH[1] + AH[2] * AH[2];
							norm2 = sqrt(norm2);
							norm3 = AH[0] * DH[0] + AH[1] * DH[1] + AH[2] * DH[2];
							alfa = acos( norm3 / (norm1 * norm2));


							/* Acceptor-ligand pair and not NH-NH bond */
							if ( lig_flag != prot_flag && lig_flag != -1 && prot_flag != -1 &&
							     !(p_type == HN && l_type == NH) && !(l_type == HN && p_type == NH))
								score += eval_block(dist, 2.4, 3.0) * eval_block(fabs(3.1416 - fabs(alfa)), 0.0f, 0.7854f) * eval_block(fabs(2.35 - fabs(beta)), 0.0f, 0.7854f);
						}

					}
				}
			}

		}

	}


	return score;
}

/**
 * @brief Interpolate non-bonded force field energy from grids with scaled terms
 *
 * @param mol Structure with the ligand
 * @param grids Grids with potentials
 * @param types Atom types for the force field grids of the ligand
 * @param min_grid minimum grid coordinates
 * @param max_grid maximum grid coordinates
 * @param nx vector for closest point for atom
 * @param ny vector for closest point for atom
 * @param nz vector for closest point for atom
 * @param max_x number of points of x axis in the grid
 * @param max_y number of points of y axis in the grid
 * @param max_z number of points of z axis in the grid
 * @param spacing grid spacing
 * @return Scaled interpolated energy of the ligand in the grid
 *
 */
float get_ff_grid_energy_ipol_scaled(MOL2 *mol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing)
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

	for ( i = 0; i < mol->n_atoms; ++i)
		++real_atoms;

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




	if ( out_flag == 0) {
		energy = 0.0f;
		e1 = e2 = 0.0f;
		for ( i = 0; i < real_atoms; ++i) {
			e1 = e2 = 0.0f;


			kpot = types[i];
			x = mol->x[i];
			y = mol->y[i];
			z = mol->z[i];

			nx[i] = (int)roundf((mol->x[i] - min_grid[0]) / spacing);
			ny[i] = (int)roundf((mol->y[i] - min_grid[1]) / spacing);
			nz[i] = (int)roundf((mol->z[i] - min_grid[2]) / spacing);

			ix = nx[i];
			iy = ny[i];
			iz = nz[i];
			q = mol->pcharges[i];


			if ( !(ix >= 1 && iy >= 1 &&
			       iz >= 1 && ix < max_x - 1 && iy < max_y - 1 &&
			       iz < max_z - 1)) {
				energy = 9999.9f;
				return energy;
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


				e1 = p1u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
				e2 = p1u * p1v * p1w * grids[8][ w0 ][ v0 ][ u0 ];

				e1 += p0u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u1 ];
				e2 += p0u * p1v * p1w * grids[8][ w0 ][ v0 ][ u1 ];

				e1 += p1u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u0 ];
				e2 += p1u * p0v * p1w * grids[8][ w0 ][ v1 ][ u0 ];

				e1 += p0u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u1 ];
				e2 += p0u * p0v * p1w * grids[8][ w0 ][ v1 ][ u1 ];

				e1 += p1u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
				e2 += p1u * p1v * p0w * grids[8][ w1 ][ v0 ][ u0 ];

				e1 += p0u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u1 ];
				e2 += p0u * p1v * p0w * grids[8][ w1 ][ v0 ][ u1 ];

				e1 += p1u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u0 ];
				e2 += p1u * p0v * p0w * grids[8][ w1 ][ v1 ][ u0 ];

				e1 += p0u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u1 ];
				e2 += p0u * p0v * p0w * grids[8][ w1 ][ v1 ][ u1 ];


				energy = energy + (0.065 * e1) + (q * e2 * 0.130);


			}


		}


	}else
		energy = 9999.9f;


	return energy;
}

/**
 * @brief Interpolate non-bonded ChemScore energy from grids with scaled terms
 *
 * @param mol Structure with the ligand
 * @param grids Grids with potentials
 * @param types Atom types for ChemScore grids of the ligand
 * @param min_grid minimum grid coordinates
 * @param max_grid maximum grid coordinates
 * @param nx vector for closest point for atom
 * @param ny vector for closest point for atom
 * @param nz vector for closest point for atom
 * @param max_x number of points of x axis in the grid
 * @param max_y number of points of y axis in the grid
 * @param max_z number of points of z axis in the grid
 * @param spacing grid spacing
 * @return Scaled interpolated energy of the ligand in the grid
 *
 */
float get_chemscore_grid_energy_ipol_scaled(MOL2 *mol, float ****grids, int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing, int metal_flag)
{

	int i = 0, real_atoms = 0, out_flag = 0, l = 0;
	float energy = 0.0f;
	register double u,   v,   w;
	register double p0u, p0v, p0w;
	register double p1u, p1v, p1w;
	register int u0,  v0,  w0;
	register int u1,  v1,  w1;
	float e1, e2, e3, x, y, z, xd, yd, zd;
	int kpot = 0;
	float lp = 0.0f, hb = 0.0f;

	/* Count only ligand atoms. Exclude all possible protein atoms inside mol */
/*	for ( i = 0; i < mol->n_atoms; ++i)*/
	real_atoms = mol->n_atoms;
/*		++real_atoms;*/

	/* Check that the ligand is inside the grid */
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
	/* If the ligand is inside, compute the nearest grid point to each atom */
	if ( out_flag == 0) {
		for ( i = 0; l < real_atoms; ++l) {
			nx[l] = (int)roundf((mol->x[l] - min_grid[0]) / spacing);
			ny[l] = (int)roundf((mol->y[l] - min_grid[1]) / spacing);
			nz[l] = (int)roundf((mol->z[l] - min_grid[2]) / spacing);
		}
		energy = 0.0f;
		/* Sum contribution from all ligand atoms */
		for ( l = 0; l < real_atoms; ++l) {
			kpot = types[l];
			if ( types[l] != -1 && nx[l] >= 1 && ny[l] >= 1 &&
			     nz[l] >= 1 && nx[l] < max_x - 1 && ny[l] < max_y - 1 &&
			     nz[l] < max_z - 1) {

/*				energy = energy + ( grids[1][nx[l]][ny[l]][nz[l]]);*/

				x = mol->x[l];
				y = mol->y[l];
				z = mol->z[l];

				xd = min_grid[0];
				yd = min_grid[1];
				zd = min_grid[2];

				w1  = (w0 = (int)(w = (x - xd) / spacing )) + 1;
				p1w = 1.0L - (p0w = w - (float)w0);

				v1  = (v0 = (int)(v = (y - yd) / spacing)) + 1;
				p1v = 1.0L - (p0v = v - (float)v0);

				u1  = (u0 = (int)(u = (z - zd) / spacing)) + 1;
				p1u = 1.0L - (p0u = u - (float)u0);


				e1 = p1u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u0 ];
				e3 = p1u * p1v * p1w * grids[5][ w0 ][ v0 ][ u0 ];


				e1 += p0u * p1v * p1w * grids[kpot][ w0 ][ v0 ][ u1 ];
				e3 += p0u * p1v * p1w * grids[5][ w0 ][ v0 ][ u1 ];

				e1 += p1u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u0 ];
				e3 += p1u * p0v * p1w * grids[5][ w0 ][ v1 ][ u0 ];


				e1 += p0u * p0v * p1w * grids[kpot][ w0 ][ v1 ][ u1 ];
				e3 += p0u * p0v * p1w * grids[5][ w0 ][ v1 ][ u1 ];

				e1 += p1u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u0 ];
				e3 += p1u * p1v * p0w * grids[5][ w1 ][ v0 ][ u0 ];

				e1 += p0u * p1v * p0w * grids[kpot][ w1 ][ v0 ][ u1 ];
				e3 += p0u * p1v * p0w * grids[5][ w1 ][ v0 ][ u1 ];

				e1 += p1u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u0 ];
				e3 += p1u * p0v * p0w * grids[5][ w1 ][ v1 ][ u0 ];

				e1 += p0u * p0v * p0w * grids[kpot][ w1 ][ v1 ][ u1 ];
				e3 += p0u * p0v * p0w * grids[5][ w1 ][ v1 ][ u1 ];

/*				if( mol->gaff_types[l] == OH)
				  e1 *= 0.5;*/ /* xlourdes advice */
/*				if( kpot == 2 || kpot == 3)
				 hb += e1;
                                else if( kpot == 0 || kpot == 19){
                                 lp += e1;
				 fprintf(stderr,"%i - %f\n",l,e1);
				}*/
				
				if( types[l] == 0)
				energy = energy + e1;


				/* If acceptor or mixed, also metals contribution */

/*				if ( metal_flag == 1)*/
				if ( (types[l] == 3 || types[l] == 4) && mol->gaff_types[l] != N3)
				{
					energy = energy + e3;
                                }/*else if( types[l] == 2){
                                        energy = energy - 0.5f*e3;
				}*/
				
			

	
                        }else if ( mol->gaff_types[l] == SH && nx[l] >= 1 && ny[l] >= 1 &&	
                             nz[l] >= 1 && nx[l] < max_x - 1 && ny[l] < max_y - 1 &&
                             nz[l] < max_z - 1) {


                                x = mol->x[l];
                                y = mol->y[l];
                                z = mol->z[l];

                                xd = min_grid[0];
                                yd = min_grid[1];
                                zd = min_grid[2];

                                w1  = (w0 = (int)(w = (x - xd) / spacing )) + 1;
                                p1w = 1.0L - (p0w = w - (float)w0);

                                v1  = (v0 = (int)(v = (y - yd) / spacing)) + 1;
                                p1v = 1.0L - (p0v = v - (float)v0);

                                u1  = (u0 = (int)(u = (z - zd) / spacing)) + 1;
                                p1u = 1.0L - (p0u = u - (float)u0);


                                e3 = p1u * p1v * p1w * grids[5][ w0 ][ v0 ][ u0 ];
                                e3 += p0u * p1v * p1w * grids[5][ w0 ][ v0 ][ u1 ];
                                e3 += p1u * p0v * p1w * grids[5][ w0 ][ v1 ][ u0 ];
                                e3 += p0u * p0v * p1w * grids[5][ w0 ][ v1 ][ u1 ];
                                e3 += p1u * p1v * p0w * grids[5][ w1 ][ v0 ][ u0 ];
                                e3 += p0u * p1v * p0w * grids[5][ w1 ][ v0 ][ u1 ];
                                e3 += p1u * p0v * p0w * grids[5][ w1 ][ v1 ][ u0 ];
                                e3 += p0u * p0v * p0w * grids[5][ w1 ][ v1 ][ u1 ];


				energy += 8.0*e3;
			}


		}

	}else
		energy = 9999.9;

/*	fprintf(stderr,"Lipo:%f, HB: %f\n",lp,hb);*/
	return energy;
}

/**
 * @brief non-bonded ChemScore energy from grids
 *
 * @param mol MOL2 with the ligand
 * @param prot MOL2 with the protein
 * @param types Atom types for ChemScore grids of the ligand
 *
 */
float get_chemscore_grid_energy_scaled(MOL2 *mol, MOL2 *prot,int *types, float min_grid[3], float max_grid[3], int *nx, int *ny, int *nz, int max_x, int max_y, int max_z, float spacing)
{

        int i = 0, real_atoms = 0, out_flag = 0, l = 0; 
        float energy = 0.0f;
        register double u,   v,   w;
        register double p0u, p0v, p0w;
        register double p1u, p1v, p1w;
        register int u0,  v0,  w0;
        register int u1,  v1,  w1; 
        float e1, e2, e3, x, y, z, xd, yd, zd; 
        int kpot = 0;

	float dx = 0.0f, dy = 0.0f, dz = 0.0f;
	float dist = 0.0f, distr = 0.0f, alfa = 0.0f, beta = 0.0f, gbeta = 0.0f;
	float heavy_nb[3], vector1[3], vector2[3], modulus1 = 0.0f, modulus2 = 0.0f;
	int neig_j = 0, neig_index = 0;

	float lp = 0.0f, hb = 0.0f, glp = 0.0f;



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


        /* If the ligand is inside, compute the nearest grid point to each atom */
        if ( out_flag == 0) {
                for ( i = 0; l < real_atoms; ++l) {
                        nx[l] = (int)roundf((mol->x[l] - min_grid[0]) / spacing);
                        ny[l] = (int)roundf((mol->y[l] - min_grid[1]) / spacing);
                        nz[l] = (int)roundf((mol->z[l] - min_grid[2]) / spacing);
                }

                energy = 0.0f;
                for ( l = 0; l < real_atoms; ++l) {

                        kpot = types[l];
			lp = 0.0f;
                        if ( types[l] != -1 && nx[l] >= 1 && ny[l] >= 1 &&
                             nz[l] >= 1 && nx[l] < max_x - 1 && ny[l] < max_y - 1 &&
                             nz[l] < max_z - 1) {

/*                        if ( types[l] != -1 ){*/
				if( types[l] != 2 && types[l] != 3)
				  continue;
                                x = mol->x[l];
                                y = mol->y[l];
                                z = mol->z[l];
				for( i = 0; i < prot->n_atoms; i++)
				{
					if( prot->ism_selection[i] != 1)
					 continue;

					if( types[l] == 0 || types[l] == 19)  /* Lipophilic */
					{
/*                                      		if( (prot->gaff_types[i] >= 67 && prot->gaff_types[i] <= 70) || prot->gaff_types[i] == 56 ||
			                     prot->gaff_types[i] == 50 || prot->gaff_types[i] == 52 || (prot->gaff_types[i] >= 19 &&
			                     prot->gaff_types[i] <= 23) || ( prot->gaff_types[i] >= 31 && prot->gaff_types[i] <= 35 ) ||
			                     prot->gaff_types[i] == 25 || prot->gaff_types[i] == 27 )
						{
							dx = x - prot->x[i];
                                                        dy = y - prot->y[i];
                                                        dz = z - prot->z[i];
							dist = (dx*dx)+(dy*dy)+(dz*dz);
							dist = sqrt(dist);
							e1 = vdw_r[ mol->gaff_types[l]-1] + vdw_r[ prot->gaff_types[i]-1] + 0.5f;
							energy += (eval_block(dist,e1,e1+3.0f)*-0.117f);
							lp +=  (eval_block(dist,e1,e1+3.0f)*-0.117f);
							glp +=  (eval_block(dist,e1,e1+3.0f)*-0.117f);
						}*/
					}else if( types[l] == 2){

                                              if( (prot->gaff_types[i] == O || prot->gaff_types[i] == OH || (prot->atoms[i] == 3 && has_hydrogens(prot,i) == 0 && prot->gaff_types[i] != N4)))
                                              {

                                                        dx = x - prot->x[i];
                                                        dy = y - prot->y[i];
                                                        dz = z - prot->z[i];
                                                        dist = (dx*dx)+(dy*dy)+(dz*dz);

                                                if (dist > 6.25f)
                                                 continue;

                                                        dist = sqrt(dist);
					

                                                 heavy_nb[0] = mol->x[l];
                                                 heavy_nb[1] = mol->y[l];
                                                 heavy_nb[2] = mol->z[l];

                                                         for(neig_j=0; neig_j < mol->n_bonds; ++neig_j)
                                                         {
                                                           if( mol->bond_a1[neig_j] == (l+1))
                                                           {
                                                                neig_index = mol->bond_a2[neig_j]-1;
                                                      if ( mol->atoms[neig_index] == 3 || mol->atoms[neig_index] == 2)
                                                      {
                                                               heavy_nb[0] = mol->x[neig_index];
                                                               heavy_nb[1] = mol->y[neig_index];
                                                               heavy_nb[2] = mol->z[neig_index];

                                                      }
                                                           }else if( mol->bond_a2[neig_j] == (l+1)){
                                                               neig_index = mol->bond_a1[neig_j]-1;
                                                      if ( mol->atoms[neig_index] == 3 || mol->atoms[neig_index] == 2)
                                                      {
                                                               heavy_nb[0] = mol->x[neig_index];
                                                               heavy_nb[1] = mol->y[neig_index];
                                                               heavy_nb[2] = mol->z[neig_index];


                                                       }

                                                           }
                                                         }


                                                        vector1[0] = heavy_nb[0] - x;
                                                        vector1[1] = heavy_nb[1] - y;
                                                        vector1[2] = heavy_nb[2] - z;

                                                        modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                        vector2[0] = prot->x[i] - x;
                                                        vector2[1] = prot->y[i] - y;
                                                        vector2[2] = prot->z[i] - z;

                                                        modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);


                                                        alfa = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
                                                        alfa = 180.0f*alfa/3.14159f;


	                                                 heavy_nb[0] = prot->x[i];
	                                                 heavy_nb[1] = prot->y[i];
	                                                 heavy_nb[2] = prot->z[i];
							 beta = 0.0f;
							 gbeta = 1.0f;
                                                         for(neig_j=0; neig_j < prot->n_bonds; ++neig_j)
                                                         {
                                                           if( prot->bond_a1[neig_j] == (i+1))
                                                           {
                                                               neig_index = prot->bond_a2[neig_j]-1;
                                                               heavy_nb[0] = prot->x[neig_index];
                                                               heavy_nb[1] = prot->y[neig_index];
                                                               heavy_nb[2] = prot->z[neig_index];

	                                                        vector1[0] = x - prot->x[i];
	                                                        vector1[1] = y - prot->y[i];
	                                                        vector1[2] = z - prot->z[i];
	                                                        modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

	                                                        vector2[0] = heavy_nb[0] - prot->x[i];
	                                                        vector2[1] = heavy_nb[1] - prot->y[i];
	                                                        vector2[2] = heavy_nb[2] - prot->z[i];
	                                                        modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);
	                                                        beta = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
								beta = 180.0f*beta/3.14159f;
								gbeta = gbeta * eval_block( fabs(180.0f-beta), 70.0f, 80.0f);

                                                           }else if( prot->bond_a2[neig_j] == (i+1)){
                                                               neig_index = prot->bond_a1[neig_j]-1;
                                                               heavy_nb[0] = prot->x[neig_index];
                                                               heavy_nb[1] = prot->y[neig_index];
                                                               heavy_nb[2] = prot->z[neig_index];

                                                                vector1[0] = x - prot->x[i];
                                                                vector1[1] = y - prot->y[i];
                                                                vector1[2] = z - prot->z[i];
                                                                modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                                vector2[0] = heavy_nb[0] - prot->x[i];
                                                                vector2[1] = heavy_nb[1] - prot->y[i];
                                                                vector2[2] = heavy_nb[2] - prot->z[i];
                                                                modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);
                                                                beta = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
                                                                beta = 180.0f*beta/3.14159f;
                                                                gbeta = gbeta * eval_block( fabs(180.0f-beta), 70.0f, 80.0f);

                                                           }
                                                         }

                                                        e1 = 0.0f;
                                                        e1 = (eval_block( (dist-1.85f), 0.25f, 0.65f)*-3.34f);
                                                        e1 *= eval_block( fabs(180.0f-alfa), 30.0f, 80.0f);
							e1 *= gbeta;
							hb += e1;

                                                        energy += e1;


                                              }

					}else if( types[l] == 3){

                                              if( (prot->gaff_types[i] == HN || prot->gaff_types[i]  == HO) )
                                              {
                                                        dx = x - prot->x[i];
                                                        dy = y - prot->y[i];
                                                        dz = z - prot->z[i];
                                                        dist = (dx*dx)+(dy*dy)+(dz*dz);
                                                        dist = sqrt(dist);

						if (dist > 6.25f)
						 continue;
						
                                                        dist = sqrt(dist);


                                                 heavy_nb[0] = prot->x[i];
                                                 heavy_nb[1] = prot->y[i];
                                                 heavy_nb[2] = prot->z[i];

                                                         for(neig_j=0; neig_j < prot->n_bonds; ++neig_j)
                                                         {
                                                           if( prot->bond_a1[neig_j] == (i+1))
                                                           {
                                                                neig_index = prot->bond_a2[neig_j]-1;
                                                      if ( prot->atoms[neig_index] == 3 || prot->atoms[neig_index] == 2)
                                                      {
                                                               heavy_nb[0] = prot->x[neig_index];
                                                               heavy_nb[1] = prot->y[neig_index];
                                                               heavy_nb[2] = prot->z[neig_index];

                                                      }
                                                           }else if( prot->bond_a2[neig_j] == (i+1)){
                                                               neig_index = prot->bond_a1[neig_j]-1;
                                                      if ( prot->atoms[neig_index] == 3 || prot->atoms[neig_index] == 2)
                                                      {
                                                               heavy_nb[0] = prot->x[neig_index];
                                                               heavy_nb[1] = prot->y[neig_index];
                                                               heavy_nb[2] = prot->z[neig_index];


                                                       }

                                                           }
                                                         }


                                                        vector1[0] = heavy_nb[0] - prot->x[i];
                                                        vector1[1] = heavy_nb[1] - prot->y[i];
                                                        vector1[2] = heavy_nb[2] - prot->z[i];

                                                        modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                        vector2[0] = dx;
                                                        vector2[1] = dy;
                                                        vector2[2] = dz;

                                                        modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);


                                                        alfa = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
							alfa = 180.0f*alfa/3.14159f;



	                                                 heavy_nb[0] = mol->x[l];
	                                                 heavy_nb[1] = mol->y[l];
	                                                 heavy_nb[2] = mol->z[l];
							 beta = 0.0f;
							 gbeta = 1.0f;
                                                         for(neig_j=0; neig_j < mol->n_bonds; ++neig_j)
                                                         {
                                                           if( mol->bond_a1[neig_j] == (l+1))
                                                           {
                                                               neig_index = mol->bond_a2[neig_j]-1;
                                                               heavy_nb[0] = mol->x[neig_index];
                                                               heavy_nb[1] = mol->y[neig_index];
                                                               heavy_nb[2] = mol->z[neig_index];

                                                                vector1[0] = prot->x[i] - x;
                                                                vector1[1] = prot->y[i] - y;
                                                                vector1[2] = prot->z[i] - z;
                                                                modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                                vector2[0] = heavy_nb[0] - x;
                                                                vector2[1] = heavy_nb[1] - y;
                                                                vector2[2] = heavy_nb[2] - z;
                                                                modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);
                                                                beta = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
                                                                beta = 180.0f*beta/3.14159f;
                                                                gbeta = gbeta * eval_block( fabs(180.0f-beta), 70.0f, 80.0f);
                                                           }else if( mol->bond_a2[neig_j] == (l+1)){
                                                               neig_index = mol->bond_a1[neig_j]-1;
                                                               heavy_nb[0] = mol->x[neig_index];
                                                               heavy_nb[1] = mol->y[neig_index];
                                                               heavy_nb[2] = mol->z[neig_index];


                                                                vector1[0] = prot->x[i] - x;
                                                                vector1[1] = prot->y[i] - y;
                                                                vector1[2] = prot->z[i] - z;
                                                                modulus1 = sqrt(vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

                                                                vector2[0] = heavy_nb[0] - x;
                                                                vector2[1] = heavy_nb[1] - y;
                                                                vector2[2] = heavy_nb[2] - z;
                                                                modulus2 = sqrt(vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);
                                                                beta = acos ( (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]) / (modulus1*modulus2) );
                                                                beta = 180.0f*beta/3.14159f;
                                                                gbeta = gbeta * eval_block( fabs(180.0f-beta), 70.0f, 80.0f);

                                                           }
                                                         }


							e1 = 0.0f;
							e1 = (eval_block( (dist-1.85f), 0.25f, 0.65f)*-3.34f);
							e1 *= eval_block( fabs(180.0f-alfa), 30.0f, 80.0f);
							e1 *= gbeta;
							hb += e1;
							energy += e1;

					       }
					}


				}
/*				if( types[l] == 0 || types[l] == 19)
                                fprintf(stderr,"%i(%i) - %f\n",l,mol->gaff_types[l],lp);*/
                        }


                }

        }else
                energy = 9999.9;

/*        fprintf(stderr,"Lipo:%f, HB: %f\n",glp,hb);*/

        return energy;
}

