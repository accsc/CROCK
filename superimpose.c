/**
 *
 *      @file superimpose.c
 *      @brief Molecule-template superimposition routines
 *
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
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



void qtrfit(double *r, double *f, int size, double u[3][3]);
void jacobi(int n, double *a, double *d, double *v);


/**
 *
 *	@brief Superimpose the ligand interaction points on a template triangle of interaction points
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *	@param mymol MOL2 with the molecule
 *	@param p1 Cartesian coordinates of the first template point
 *	@param p2 Cartesian coordinates of the second template point
 *	@param p3 Cartesian coordinates of the third template point
 *	@param xi Index of the first ligand point
 *	@param yi Index of the second ligand point
 *	@param zi Index of the third ligand point
 *	@return 0 on success
 *
 */
int superimpose_triangle(MOL2 **mymol, float p1[3], float p2[3], float p3[3], int xi, int yi, int zi)
{

	float a[3][3], b[3][3];
	int i = 0, l = 0, j = 0, k = 0;
	double center[3];
	double comt[3];
	double vec[3];
	double center2[3];
	double rotM[3][3];
	double Xrot = 0.0f, Yrot = 0.0f, Zrot = 0.0f;
	double a1, a2, a3;
	double b1, b2, b3, b4, b5, b6, b7, b8, b9;
	double *refcoor = NULL, *mvcoor = NULL;
	double m[9];
	double x, y, z;
	MOL2 *mol = NULL;

	mol = *mymol;

	refcoor = calloc( sizeof(double), 9);
	mvcoor = calloc( sizeof(double), 9);



	/* COM of reference points and translation */
	center[0] = p1[0] + p2[0] + p3[0];
	center[1] = p1[1] + p2[1] + p3[1];
	center[2] = p1[2] + p2[2] + p3[2];

	center[0] /= 3.0f;
	center[1] /= 3.0f;
	center[2] /= 3.0f;

	refcoor[0] = p1[0] - center[0];
	refcoor[1] = p1[1] - center[1];
	refcoor[2] = p1[2] - center[2];

	refcoor[3] = p2[0] - center[0];
	refcoor[4] = p2[1] - center[1];
	refcoor[5] = p2[2] - center[2];

	refcoor[6] = p3[0] - center[0];
	refcoor[7] = p3[1] - center[1];
	refcoor[8] = p3[2] - center[2];


	/* COM of 3 three target points and trasnlate the whole molecule */
	center2[0] = center2[1] = center2[2] = 0;

	center2[0] += mol->x[xi];
	center2[1] += mol->y[xi];
	center2[2] += mol->z[xi];

	center2[0] += mol->x[yi];
	center2[1] += mol->y[yi];
	center2[2] += mol->z[yi];

	center2[0] += mol->x[zi];
	center2[1] += mol->y[zi];
	center2[2] += mol->z[zi];

	center2[0] /= 3.0f;
	center2[1] /= 3.0f;
	center2[2] /= 3.0f;

	for ( i = 0; i < mol->n_atoms; ++i) {
		mol->x[i] -= center2[0];
		mol->y[i] -= center2[1];
		mol->z[i] -= center2[2];
	}

	mvcoor[0] = mol->x[xi];
	mvcoor[1] = mol->y[xi];
	mvcoor[2] = mol->z[xi];

	mvcoor[3] = mol->x[yi];
	mvcoor[4] = mol->y[yi];
	mvcoor[5] = mol->z[yi];

	mvcoor[6] = mol->x[zi];
	mvcoor[7] = mol->y[zi];
	mvcoor[8] = mol->z[zi];

	qtrfit(refcoor, mvcoor, 3, rotM);


	for (k = 0, i = 0; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			m[k++] = rotM[i][j];

	for (i = 0; i < mol->n_atoms; ++i) {
		x = mol->x[i];
		y = mol->y[i];
		z = mol->z[i];
		mol->x[i] = m[0] * x + m[1] * y + m[2] * z;
		mol->y[i] = m[3] * x + m[4] * y + m[5] * z;
		mol->z[i] = m[6] * x + m[7] * y + m[8] * z;
	}

	for ( l = 0; l < mol->n_atoms; ++l) {
		mol->x[l] += center[0];
		mol->y[l] += center[1];
		mol->z[l] += center[2];

	}

	free(refcoor);
	free(mvcoor);

	return 0;

}

/*
 *  This code is extracted from the OpenBabel project and it is under 
 *  the GPLv2 license, fully compatible with this project license.
 *
 *
 *  Open Babel AUTHORS
 *  -------------------
 *
 *  Open Babel would not be what it is without the help of a cast of many.
 *  We are fundamentally a community project and aim to offer open
 *  development, responsive to users and contributors alike.
 *
 *  The Open Babel project is currently loosely managed by:
 *  Michael Banck
 *  Geoff Hutchison
 *  Chris Morley
 *  Noel O'Boyle
 *  Tim Vandermeersch
 *
 *  Full list of contributors can be found on the Open Babel website:
 *  <http://openbabel.org/wiki/THANKS>
 *
 */



/****************/
/* OB2.3.0 code */
/****************/

void qtrfit(double *r, double *f, int size, double u[3][3])
{
	register int i;
	double xxyx, xxyy, xxyz;
	double xyyx, xyyy, xyyz;
	double xzyx, xzyy, xzyz;
	double d[4], q[4];
	double c[16], v[16];
	double rx, ry, rz, fx, fy, fz;

	/* generate the upper triangle of the quadratic form matrix */

	xxyx = 0.0;
	xxyy = 0.0;
	xxyz = 0.0;
	xyyx = 0.0;
	xyyy = 0.0;
	xyyz = 0.0;
	xzyx = 0.0;
	xzyy = 0.0;
	xzyz = 0.0;

	for (i = 0; i < size; ++i) {
		rx = r[i * 3];
		ry = r[i * 3 + 1];
		rz = r[i * 3 + 2];
		fx = f[i * 3];
		fy = f[i * 3 + 1];
		fz = f[i * 3 + 2];

		xxyx += fx * rx;
		xxyy += fx * ry;
		xxyz += fx * rz;
		xyyx += fy * rx;
		xyyy += fy * ry;
		xyyz += fy * rz;
		xzyx += fz * rx;
		xzyy += fz * ry;
		xzyz += fz * rz;
	}
	c[4 * 0 + 0] = xxyx + xyyy + xzyz;

	c[4 * 0 + 1] = xzyy - xyyz;
	c[4 * 1 + 1] = xxyx - xyyy - xzyz;

	c[4 * 0 + 2] = xxyz - xzyx;
	c[4 * 1 + 2] = xxyy + xyyx;
	c[4 * 2 + 2] = xyyy - xzyz - xxyx;

	c[4 * 0 + 3] = xyyx - xxyy;
	c[4 * 1 + 3] = xzyx + xxyz;
	c[4 * 2 + 3] = xyyz + xzyy;
	c[4 * 3 + 3] = xzyz - xxyx - xyyy;

	/* diagonalize c */

	jacobi(4, c, d, v);

	/* extract the desired quaternion */

	q[0] = v[4 * 0 + 3];
	q[1] = v[4 * 1 + 3];
	q[2] = v[4 * 2 + 3];
	q[3] = v[4 * 3 + 3];

	/* generate the rotation matrix */

	u[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
	u[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
	u[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);

	u[0][1] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
	u[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
	u[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);

	u[0][2] = 2.0 * (q[3] * q[1] - q[0] * q[2]);
	u[1][2] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
	u[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
}

void jacobi(int n, double *a, double *d, double *v)
{
	double onorm, dnorm;
	double b, dma, q, t, c, s;
	double atemp, vtemp, dtemp;
	register int i, j, k, l;
	int nrot;

	int MAX_SWEEPS = 50;

	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++)
			v[n * i + j] = 0.0;
		v[n * j + j] = 1.0;
		d[j] = a[n * j + j];
	}

	nrot = MAX_SWEEPS;
	for (l = 1; l <= nrot; l++) {
		dnorm = 0.0;
		onorm = 0.0;
		for (j = 0; j < n; j++) {
			dnorm += (double)fabs(d[j]);
			for (i = 0; i < j; i++)
				onorm += (double)fabs(a[n * i + j]);
		}
		if ((onorm / dnorm) <= 1.0e-12)
			goto Exit_now;

		for (j = 1; j < n; j++) {
			for (i = 0; i <= j - 1; i++) {

				b = a[n * i + j];
				if (fabs(b) > 0.0) {
					dma = d[j] - d[i];
					if ((fabs(dma) + fabs(b)) <=  fabs(dma))
						t = b / dma;
					else{
						q = 0.5 * dma / b;
						t = 1.0 / ((double)fabs(q) + (double)sqrt(1.0 + q * q));
						if (q < 0.0)
							t = -t;
					}

					c = 1.0 / (double)sqrt(t * t + 1.0);
					s = t * c;
					a[n * i + j] = 0.0;

					for (k = 0; k <= i - 1; k++) {
						atemp = c * a[n * k + i] - s * a[n * k + j];
						a[n * k + j] = s * a[n * k + i] + c * a[n * k + j];
						a[n * k + i] = atemp;
					}

					for (k = i + 1; k <= j - 1; k++) {
						atemp = c * a[n * i + k] - s * a[n * k + j];
						a[n * k + j] = s * a[n * i + k] + c * a[n * k + j];
						a[n * i + k] = atemp;
					}

					for (k = j + 1; k < n; k++) {
						atemp = c * a[n * i + k] - s * a[n * j + k];
						a[n * j + k] = s * a[n * i + k] + c * a[n * j + k];
						a[n * i + k] = atemp;
					}

					for (k = 0; k < n; k++) {
						vtemp = c * v[n * k + i] - s * v[n * k + j];
						v[n * k + j] = s * v[n * k + i] + c * v[n * k + j];
						v[n * k + i] = vtemp;
					}

					dtemp = c * c * d[i] + s * s * d[j] - 2.0 * c * s * b;
					d[j] = s * s * d[i] + c * c * d[j] +  2.0 * c * s * b;
					d[i] = dtemp;
				}
			}
		}
	}

Exit_now:
	nrot = l;

	for (j = 0; j < n - 1; j++) {
		k = j;
		dtemp = d[k];
		for (i = j + 1; i < n; i++)
			if (d[i] < dtemp) {
				k = i;
				dtemp = d[k];
			}

		if (k > j) {
			d[k] = d[j];
			d[j] = dtemp;
			for (i = 0; i < n; i++) {
				dtemp = v[n * i + k];
				v[n * i + k] = v[n * i + j];
				v[n * i + j] = dtemp;
			}
		}
	}
}

