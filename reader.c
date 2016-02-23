/**
 *
 *
 *
 *      @file reader.c
 *      @brief Read and allocation routine for PDB and MOL2 molecules
 *
 *        - Automatic parser: atom type, dihedrals, contacts, angles and bonds.
 *        - Import & export topology.
 *        - Modular energy terms.
 *        - Analytical gradients for most of the terms
 *        - Parser based on residue and atom name
 *        - Automatic interactions ligand-protein
 *
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *	@date 01/10/2010
 *
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation version 2 of the License.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

/* Molecular library definitions */
#include <libmol2.h>
/* GAFF force field definitions and constants */
#include <gaff.h>

#include <time.h>

/* General tools for reader */
#include "groups.c"
/* Ring detection routines */
#include "rings.c"

#define G_PI 3.14159265358979323846

/* Energetic calculation routines */
#include "energy.c"

#define MAX_BUFFER 1024


/**
 *
 *	@brief	Read and parse the contents of a PDB file 
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *
 *	@param mymol MOL2 structure to allocate and load
 *	@param finput_name Filename of the PDB
 *	@param import Import topology flag. 0. Do not import topology from files. 1. Import topology from top files
 *
 *      @return 0 if success
 *
 */
int PDB_reader(MOL2 **mymol, char *finput_name, int import)
{

	double elapsed;
	double treadab;
	double tring;
	double tdihe;
	double tiped;

	FILE *input;
	char *line = NULL;
	int molecules = 0;
	int i = 0;
	int i2 = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int step = 0;
	int cpair = 0;
	char tmp_atom[MAX_BUFFER];
	char tmp_bond[MAX_BUFFER];
	char tmp_radius[6];
	char tmp_charge[10];
	int rec_atoms = 0;
	MOL2 *mols = NULL;
	float pcharge = 0.0;
	int current_atom = 0;
	int current_bond = 0;
	int total_ppp = 0;
	char *adj = NULL;
	RING *test = NULL;
	char c = 0;
	char *input_name = NULL;
	int mode = 0;
	float x, y, z;
	int biggrad = 0;
	float lastforce = 0;
	float tmpforce = 0.0;
	float cbiggrad = 0.0f;
	float lastenergy = 999999.9f;
	float bestenergy = 0.0f;
	float tmpgrad = 0.0f;
	float dist = 0.0f;
	char name1[20];
	char name2[20];
	char name3[20];
	char name4[20];
	char name5[20];
	int ppp_count = 0;
	int my_type = 0;
	char myx[12];
	char myy[12];
	char myz[12];
	int bonds1, bonds2, bonds3;
	int flag = 0;
	int flag2 = 0;
	int flag3 = 0;
	int t1, t2, t3, t4;
	double bondE = 0;
	double angleE = 0;
	double torE = 0;
	double totalE = 0;
	double mytot;
	double vec1[3], vec2[3];
	double uvec1[3], uvec2[3];
	double *lvec1, *lvec2;
	double tor1;
	double tor2;
	double tor3;
	float thetha;
	int ttor = 0;
	int ewg = 0;
	int ewg2 = 0;
	int bondis = 0;

	RESIDUE *list_r = NULL;
	RESIDUE tmp_r;
	ROTAMER *tmp_rota = NULL;
	FILE *inres = NULL;
	FILE *rotlib = NULL;
	int tmp_resn_i = 0;
	char tmp_resn[10];
	int total_res = 0;
	int resn[1000];
	float target1[3], target2[3];
	float mx, my, mz;
	float brotE;
	int bestrot = -1;

	int cur_angle = 0;
	int *tor_type;

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

	float *xc, *yc, *zc;
	float *xb, *yb, *zb;

	float idi, pn, pk, phi0;

	double A, B;
	double epsilon;
	double tmpE;
	double vdwE;
	float r6, r12;

	float rx = 0.0f, ry = 0.0f, rz = 0.0f;
	int matchs[3], nmatchs = 0;

	int conformers = 0;
	int current_conformer = 0;
	int no_count = 0;

	/* FLAGS */
	int gradflag;
	int verboseflag;

        char res_num[10];
        char res_type[10];




/* CONTROL FLAGS */
	gradflag = 1;
	verboseflag = 0;



	if ( (input = fopen(finput_name, "r")) == NULL) {
		fprintf(stderr, "Error. Cant open file %s.\n", finput_name);
		fflush(stderr);
		exit(-1);
	}



	if ( (line = (char*)calloc(1, MAX_BUFFER + 2)) == NULL) {
		fprintf(stderr, "Error. Cant al memory.\n");
		fflush(stderr);
		exit(-2);
	}


	while ( fgets(line, MAX_BUFFER, input) ) {
		if ( (strstr(line, "ATOM") != NULL || strstr(line, "HETATM") != NULL) && no_count == 0)
			++molecules;

		if ( (line[0] == 'T' && line[1] == 'E' && line[2] == 'R') || (line[0] == 'E' && line[1] == 'N' && line[2] == 'D' && line[3] != 'M') ) {
			conformers++;
			no_count = 1;
		}

	}

	rewind(input);

	if ( (mols = (MOL2*)calloc(sizeof(MOL2), 1)) == NULL) {
		fprintf(stderr, "Error. Cannot allocate big memory chunk.\n");
		fflush(stderr);
		exit(-2);
	}

	*mymol = mols;

	i = 0;
	mols->n_atoms = molecules;
	mols->x = (float*)calloc(sizeof(float), mols[i].n_atoms);
	mols->y = (float*)calloc(sizeof(float), mols[i].n_atoms);
	mols->z = (float*)calloc(sizeof(float), mols[i].n_atoms);
	mols->pcharges = (float*)calloc(sizeof(float), mols[i].n_atoms);
	mols->radius = (float*)calloc(sizeof(float), mols[i].n_atoms);
	mols->atoms = (int*)calloc(sizeof(int), mols[i].n_atoms);
	mols->ringer = (int*)calloc(sizeof(int), mols[i].n_atoms);
	mols->aromatic = (int*)calloc(sizeof(int), mols[i].n_atoms);
	mols->bond_dist = (float*)calloc(sizeof(float), mols[i].n_atoms * 8);
	mols->grads_X = (float*)calloc(sizeof(float), mols[i].n_atoms);
	mols->grads_Y = (float*)calloc(sizeof(float), mols[i].n_atoms);
	mols->grads_Z = (float*)calloc(sizeof(float), mols[i].n_atoms);
	mols->backbone = (int*)calloc(sizeof(int), mols[i].n_atoms);
	mols->selection = (int*)calloc(sizeof(int), mols->n_atoms + 10);
	mols->fragment_flag = (int*)calloc(sizeof(int), mols[i].n_atoms);
	mols->conformers = (CONFORMER*)calloc(sizeof(CONFORMER), conformers);
	mols->n_fragments = 0;

        mols->res_num = (int*)calloc(sizeof(int*), mols->n_atoms + 10);
        mols->res_names = (char**)calloc(sizeof(char**), mols->n_atoms + 1);
        mols->atom_names = (char**)calloc(sizeof(char**), mols->n_atoms + 1);
        for( i = 0; i < mols->n_atoms; i++)
        {
                mols->res_names[i] = (char *) calloc(sizeof(char*),4);
                mols->atom_names[i] = (char *) calloc(sizeof(char*),4);
        }
        i = 0;
        mols->res_type = (int*)calloc(sizeof(int*), mols->n_atoms + 10);


	for ( i = 0; i < conformers; ++i) {
		mols->conformers[i].x = (float*)calloc(sizeof(float), mols->n_atoms);
		mols->conformers[i].y = (float*)calloc(sizeof(float), mols->n_atoms);
		mols->conformers[i].z = (float*)calloc(sizeof(float), mols->n_atoms);
		mols->conformers[i].pcharges = (float*)calloc(sizeof(float), mols->n_atoms);
		mols->conformers[i].radius = (float*)calloc(sizeof(float), mols->n_atoms);
	}


	i = 0;

	mols->gaff_types = (int*)calloc(sizeof(int), mols[i].n_atoms + 1);
	mols->ism_types = (int *)calloc(sizeof(int), mols[i].n_atoms + 1);
        mols->ism_selection = (int *)calloc(sizeof(int), mols[i].n_atoms + 1);


	if ( verboseflag == 1)
		printf("Cleaning up the gradients.\n");

	if ( gradflag == 1) {
		for (j = 0; j < mols->n_atoms; ++j) {
			mols->grads_X[j] = 0.0;
			mols->grads_Y[j] = 0.0;
			mols->grads_Z[j] = 0.0;
		}
	}
	j = 0;
	current_atom = 0;
	current_conformer = 0;
	mols->nconformers = conformers;
	mols->default_conformer = 0;

	while ( fgets(line, MAX_BUFFER, input) ) {
#ifdef DEBUG
		fprintf(stderr, "Estado: Linea:%s", line);
		fflush(stderr);
#endif

		if ( (line[0] == 'T' && line[1] == 'E' && line[2] == 'R') || (line[0] == 'E' && line[1] == 'N' && line[2] == 'D' && line[3] != 'M') ) {
			current_conformer++;
			current_atom = 0;
		}else if ( strstr(line, "ATOM") != NULL || strstr(line, "HETATM") != NULL) {
			sscanf(line, "%*s %*d %4s", tmp_atom);
			strncpy(tmp_charge, &line[60], 8);
			strncpy(tmp_radius, &line[56], 4);


			if ( current_conformer == 0) {
				mols->pcharges[current_atom] = atof(tmp_charge);
				mols->radius[current_atom] = atof(tmp_radius);

				if ( tmp_atom[0] == 'C' && ( tmp_atom[1] != 'l' && tmp_atom[1] != 'L'))
					mols->atoms[current_atom] = 1;
				else if ( tmp_atom[0] == 'O')
					mols->atoms[current_atom] = 2;
				else if ( tmp_atom[0] == 'N')
					mols->atoms[current_atom] = 3;
				else if ( tmp_atom[0] == 'H')
					mols->atoms[current_atom] = 4;
				else if ( tmp_atom[0] == 'P')
					mols->atoms[current_atom] = 5;
				else if ( tmp_atom[0] == 'S')
					mols->atoms[current_atom] = 6;
				else if ( tmp_atom[0] == 'I')
					mols->atoms[current_atom] = 7;
				else if ( tmp_atom[0] == 'B' && (tmp_atom[1] == 'r' || tmp_atom[1] == 'R'))
					mols->atoms[current_atom] = 8;
				else if ( tmp_atom[0] == 'C' && (tmp_atom[1] == 'l' || tmp_atom[1] == 'L'))
					mols->atoms[current_atom] = 9;
	                        else if( tmp_atom[0] == 'F' &&  tmp_atom[1] != 'E' && tmp_atom[1] != 'e')
	                                mols->atoms[current_atom] = 10;
	                        else if(tmp_atom[0] == 'Z' && ( tmp_atom[1] == 'N' || tmp_atom[1] == 'n' )) 
                                   	mols->atoms[current_atom] = 100;
	                        else if(tmp_atom[0] == 'M' && (tmp_atom[1] == 'N' || tmp_atom[1] == 'n') )
                                   	mols->atoms[current_atom] = 101;
	                        else if(tmp_atom[0] == 'M' && (tmp_atom[1] == 'G' || tmp_atom[1] == 'g' ))
                                   	mols->atoms[current_atom] = 102;
	                        else if(tmp_atom[0] == 'C' && tmp_atom[1] == '0')
                                   	mols->atoms[current_atom] = 105;
	                        else if(tmp_atom[0] == 'F' && (tmp_atom[1] == 'E' || tmp_atom[1] == 'e'))
                                   	mols->atoms[current_atom] = 106;
	                        else if(tmp_atom[0] == 'K')
                                   	mols->atoms[current_atom] = 103;
	                        else if(tmp_atom[0] == 'N' && (tmp_atom[1] == 'A' || tmp_atom[1] == 'a'))
                                   	mols->atoms[current_atom] = 104;
				else{

					if ( tmp_atom[1] == 'H')
						mols->atoms[current_atom] = 4;
					else{
						mols->atoms[current_atom] = 11;
						mols->gaff_types[current_atom] = C3;
						fprintf(stderr, "Reader %i: Error. Do not know the atom type $%s$\n", __LINE__,tmp_atom);
					}
				}

	                        strncpy(res_num,&line[23],6);
                                mols->res_num[current_atom] = atoi(res_num);

                                strncpy(res_type,&line[17],3);
                                mols->res_type[current_atom] = name_to_number_residue(res_type);

                                strncpy(mols->res_names[current_atom],res_type,3);
                                strncpy(mols->atom_names[current_atom],&line[12],4);



			}

			strncpy(myx, &line[29], 10);
			strncpy(myy, &line[38], 10);
			strncpy(myz, &line[46], 10);

			if ( current_conformer == 0) {
				mols->x[current_atom] = atof(myx);
				mols->y[current_atom] = atof(myy);
				mols->z[current_atom] = atof(myz);
			}

			mols->conformers[current_conformer].x[current_atom] = atof(myx);
			mols->conformers[current_conformer].y[current_atom] = atof(myy);
			mols->conformers[current_conformer].z[current_atom] = atof(myz);
			mols->conformers[current_conformer].pcharges[current_atom] = atof(tmp_charge);
			mols->conformers[current_conformer].radius[current_atom] = atof(tmp_radius);

			++current_atom;


		}
	}


	free(line);
	fclose(input);

	mol_percieve(&mols, import);


	return 0;
}

/**
 *
 *	@brief Percieve connectivity of the molecule (bonds, angles, torsionals, rings, aromaticity, etc.)
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *
 *	@param mymol MOL2 structure which contains the molecule
 *	@param import flag to import (1) the topology or not (0)
 *	@return 0 on success
 *
 */
int mol_percieve(MOL2 **mymol, int import)
{

	MOL2 *mols = NULL;
	int i = 0;
	int i2 = 0;
	int j = 0;
	int k = 0;
	int l = 0;

	int t1, t2, t3, t4;
	double bondE = 0;
	double angleE = 0;
	double torE = 0;
	double totalE = 0;
	double mytot;
	double vec1[3], vec2[3];
	double uvec1[3], uvec2[3];
	double *lvec1, *lvec2;
	double tor1;
	double tor2;
	double tor3;
	float thetha;
	int ttor = 0;
	int ewg = 0;
	int ewg2 = 0;
	int bondis = 0;
	int *ringer = NULL;
	int *aro = NULL;

	float rx = 0.0f, ry = 0.0f, rz = 0.0f;
	int step = 0;
	int cpair = 0;
	int *vecinos;
	int *vecinos2;
	int vecinos4[16];
	int bonds1, bonds2, bonds3;
	int flag = 0;
	int flag2 = 0;
	int flag3 = 0;
	float dist = 0.0f;


	int angleflag;
	int torflag;
	int vdwflag;
	int verboseflag;
	int matchs[3], nmatchs = 0;


	int *tor_type;

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

	float *xc, *yc, *zc;
	float *xb, *yb, *zb;

	float idi, pn, pk, phi0;

	double A, B;
	double epsilon;
	double vdwE;
	float r6, r12;

	int cur_angle = 0;

	mols = *mymol;


	torflag = 0;
	vdwflag = 0;
	angleflag = 0;
	verboseflag = 0;



	if ( import != 0) {

		import_topology(&mols);


	}else{
		if ( verboseflag == 1)
			printf("Allocating %i bonds.\n", mols->n_atoms * mols->n_atoms);

		mols->bonds = (int*)realloc(mols->bonds, sizeof(int) * mols->n_atoms * 8);
		mols->bond_a1 = (int*)realloc(mols->bond_a1, sizeof(int) * mols->n_atoms * 8);
		mols->bond_a2 = (int*)realloc(mols->bond_a2, sizeof(int) * mols->n_atoms * 8);

		j = 0;

/******************************************************************************
*                                                                            *
*                             DETECT BONDS                                   *
*                                                                            *
*                                                                            *
******************************************************************************/


/*
 * These numbers are not random. Extracted from CCCBDB. http://cccbdb.nist.gov/
 * Some corrections are applied regarding common sense.
 */
		for (i = 0; i < mols->n_atoms; ++i) {
			for (k = i + 1; k < (mols->n_atoms); k++) {
				/* C - C or C = C or aromatic rings */
				if ( mols->atoms[i] == 1 && mols->atoms[k] == 1) {
					rx = mols->x[i] - mols->x[k];
					ry = mols->y[i] - mols->y[k];
					rz = mols->z[i] - mols->z[k];

					dist = sqrt( rx * rx + ry * ry + rz * rz);

					if ( dist <= 1.26 && dist >= 0.5 ) {
						mols->bonds[j] = 3;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						if ( verboseflag == 1)
							printf("Bond C#C. %i %i\n", i + 1, k + 1);
						j++;
					}else if (dist > 1.26 && dist <= 1.35) {
						mols->bonds[j] = 2;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond C=C. %i %i\n", i + 1, k + 1);
					}else if (dist > 1.35 && dist <= 1.9) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						if ( verboseflag == 1)
							printf("Bond C-C. %i %i\n", i + 1, k + 1);
						j++;
					}

					/* C  - N or C=N*/
				}else if ( (mols->atoms[i] == 1 && mols->atoms[k] == 3) || (mols->atoms[i] == 3 && mols->atoms[k] == 1) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.85 && dist >= 1.33) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						if ( verboseflag == 1)
							printf("Bond C-N. %i %i\n", i + 1, k + 1);
						j++;
					}else if ( dist < 1.33 && dist >= 1.19) {
						mols->bonds[j] = 2;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						if ( verboseflag == 1)
							printf("Bond C=N. %i %i\n", i + 1, k + 1);
						j++;
					}else if ( dist < 1.19) {
						mols->bonds[j] = 3;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						if ( verboseflag == 1)
							printf("Bond C#N. %i %i\n", i + 1, k + 1);
						j++;
					}

					/* C - O or C=O*/
				}else if ( (mols->atoms[i] == 1 && mols->atoms[k] == 2) || (mols->atoms[i] == 2 && mols->atoms[k] == 1) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist > 1.30 && dist <= 1.69 ) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						if ( verboseflag == 1) {
							printf("Bond C-O. %i %i\n", i + 1, k + 1);
							printf("%f.\n", dist);
						}
						j++;

					}else if ( dist <= 1.30) {
						mols->bonds[j] = 2;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1) {
							printf("Bond C=O. %i %i\n", i + 1, k + 1);
							printf("%f.\n", dist);
						}
					}
					/* C - H*/
				}else if ( (mols->atoms[i] == 1 && mols->atoms[k] == 4) || (mols->atoms[i] == 4 && mols->atoms[k] == 1) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.20) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond C-H. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 3 && mols->atoms[k] == 4) || (mols->atoms[i] == 4 && mols->atoms[k] == 3) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.20) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond N-H. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 4 && mols->atoms[k] == 2) || (mols->atoms[i] == 2 && mols->atoms[k] == 4) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.20) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond O-H. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 4 && mols->atoms[k] == 6) || (mols->atoms[i] == 6 && mols->atoms[k] == 4) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.48) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond S-H. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 1 && mols->atoms[k] == 6) || (mols->atoms[i] == 6 && mols->atoms[k] == 1) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.65) {
						mols->bonds[j] = 2;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond C=S. %i %i\n", i + 1, k + 1);
					}else if (dist > 1.65 && dist <= 1.91) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						if ( verboseflag == 1)
							printf("Bond C-S. %i %i\n", i + 1, k + 1);
						j++;
					}
				}else if ( (mols->atoms[i] == 3 && mols->atoms[k] == 6) || (mols->atoms[i] == 6 && mols->atoms[k] == 3) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 2.0) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond N-S. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 2 && mols->atoms[k] == 6) || (mols->atoms[i] == 6 && mols->atoms[k] == 2) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.55) {
						mols->bonds[j] = 2;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond S=O. %i %i\n", i + 1, k + 1);
					}else if ( dist > 1.55 && dist <= 1.8) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond S-O. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 2 && mols->atoms[k] == 3) || (mols->atoms[i] == 3 && mols->atoms[k] == 2) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.26) {
						mols->bonds[j] = 2;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond O=N. %i %i\n", i + 1, k + 1);
					}else if (dist <= 1.51 && dist > 1.26) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond O-N. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 1 && mols->atoms[k] == 10) || (mols->atoms[i] == 10 && mols->atoms[k] == 1) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.6) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond C-F. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 1 && mols->atoms[k] == 8) || (mols->atoms[i] == 8 && mols->atoms[k] == 1) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 2.0) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond C-Br. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 1 && mols->atoms[k] == 9) || (mols->atoms[i] == 9 && mols->atoms[k] == 1) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.9) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond C-Cl. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 1 && mols->atoms[k] == 7) || (mols->atoms[i] == 7 && mols->atoms[k] == 1) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 2.2) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond C-I. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 2 && mols->atoms[k] == 2)) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist >= 1.28 && dist <= 2.1) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond O-O. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 3 && mols->atoms[k] == 3) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist >= 1.4 && dist <= 1.8) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond N-N. %i %i\n", i + 1, k + 1);
					}else if ( dist < 1.4 && dist >= 1.20) {
						mols->bonds[j] = 2;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond N=N. %i %i\n", i + 1, k + 1);
					}else if ( dist < 1.20 && dist >= 0.9) {
						mols->bonds[j] = 3;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond N#N. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 5 && mols->atoms[k] == 2) || (mols->atoms[i] == 2 && mols->atoms[k] == 5) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.9) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond P-O. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 5 && mols->atoms[k] == 3) || (mols->atoms[i] == 3 && mols->atoms[k] == 5) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.99) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond N-P. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 5 && mols->atoms[k] == 1) || (mols->atoms[i] == 1 && mols->atoms[k] == 5) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 2.0) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond C-P. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 5 && mols->atoms[k] == 6) || (mols->atoms[i] == 6 && mols->atoms[k] == 5) ) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 2.0) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond S-P. %i %i\n", i + 1, k + 1);
					}
				}else if ( (mols->atoms[i] == 6 && mols->atoms[k] == 6)) {
					dist = sqrt(pow((mols->x[i] - mols->x[k]), 2) + pow((mols->y[i] - mols->y[k]), 2) + pow((mols->z[i] - mols->z[k]), 2));
					if ( dist <= 1.88 ) {
						mols->bonds[j] = 2;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond S=S. %i %i\n", i + 1, k + 1);
					}else if ( dist > 1.88 && dist < 2.2) {
						mols->bonds[j] = 1;
						mols->bond_a1[j] = i + 1;
						mols->bond_a2[j] = k + 1;
						mols->bond_dist[j] = dist;
						j++;
						if ( verboseflag == 1)
							printf("Bond S-S. %i %i\n", i + 1, k + 1);
					}
				}


			}
		}

		mols->n_bonds = j;

		if ( verboseflag == 1)
			printf("Total Bonds: %i.\n", mols->n_bonds);


		ringer = mols->ringer;
		aro =  mols->aromatic;


/******************************************************************************
 *
 *
 *                              DETECT RINGS
 *
 *
 *****************************************************************************/

		if ( mols->n_atoms <= 300)
			get_number_of_rings2(mols[0], &ringer, &aro);

/*****************************************************************************
 *
 *
 *                                ASSIGN TYPES
 *
 *
 ****************************************************************************/

		for (i = 0; i < mols->n_atoms; ++i)
			mols->gaff_types[i] = -1;


        vecinos = (int *) calloc(sizeof(int *), mols->n_atoms+1);
        vecinos2 = (int *) calloc(sizeof(int *), mols->n_atoms+1);

		for (i = 0; i < mols->n_atoms; ++i) {
/* Get vecinos */
			for (j = 0; j < mols->n_atoms; ++j) {
				vecinos[j] = 0;
				vecinos2[j] = 0;
			}
			k = 0;
			for (j = 0; j < mols->n_bonds; ++j) {
				if ( mols->bond_a1[j] == (i + 1)) {
					vecinos[k] = mols->bond_a2[j] - 1;
					vecinos2[k] = mols->bonds[j];
					k++;
				}else if ( mols->bond_a2[j] == (i + 1)) {
					vecinos[k] = mols->bond_a1[j] - 1;
					vecinos2[k] = mols->bonds[j];
					k++;
				}
			}

/* Assign types */

			if ( mols->atoms[i] == 1 ) {

				if ( mols->aromatic[i] == 1) {
					mols->gaff_types[i] = CA;
					if ( verboseflag == 1)
						printf("%i tipo CA.\n", i);
				}else{
					if ( k == 4) {

						if ( mols->ringer[i] == 3) {
							mols->gaff_types[i] = CX;
							if ( verboseflag == 1)
								printf("%i tipo CX.\n", i);
						}else if ( mols->ringer[i] == 4) {
							mols->gaff_types[i] = CY;
							if ( verboseflag == 1)
								printf("%i tipo CY.\n", i);
						}else{
							mols->gaff_types[i] = C3;
							if ( verboseflag == 1)
								printf("%i tipo C3.\n", i);
						}
					}else if ( k == 3) {
						if ( mols->ringer[i] == 3) {
							mols->gaff_types[i] = CU;
							if ( verboseflag == 1)
								printf("%i tipo CU.\n", i);
						}else if ( mols->ringer[i] == 4) {
							mols->gaff_types[i] = CV;
							if ( verboseflag == 1)
								printf("%i tipo CV.\n", i);
						}else{
							ewg = 0;
							ewg += get_number_bond_by_atom(mols[0], i + 1, 2, 2);
							ewg += get_number_bond_by_atom(mols[0], i + 1, 2, 6);
							ewg2 = get_number_bond_by_atom(mols[0], i + 1, 2, 3);
							if ( ewg != 0) {
								mols->gaff_types[i] = C;
								if ( verboseflag == 1)
									printf("%i tipo C.\n", i);
							}else if ( ewg2 == 3) {
								mols->gaff_types[i] = CZ;
								if ( verboseflag == 1)
									printf("%i tipo CZ.\n", i);
							}else{
								ewg = get_bonds(mols[0], i + 1, 2);
								ewg2 = get_bonds(mols[0], i + 1, 1);
								if ( ewg > 0 && ewg2 > 0) {
									ewg = 0;
									for (k = 0; k < 4; ++k)
										if ( vecinos2[k] == 1)
											ewg += get_bonds(mols[0], vecinos[k] + 1, 2);
									if ( ewg > 0) {

										if ( mols->ringer[i] != 0) {
											mols->gaff_types[i] = CC;
											if ( verboseflag == 1)
												printf("%i tipo CC.\n", i);
										}else{
											mols->gaff_types[i] = CE;
											if ( verboseflag == 1)
												printf("%i tipo CE.\n", i);
										}
									}else{
										mols->gaff_types[i] = C2;
										if ( verboseflag == 1)
											printf("%i tipo C2.\n", i);
									}
								}else{
									mols->gaff_types[i] = C2;
									if ( verboseflag == 1)
										printf("%i tipo C2.\n", i);
								}

							}
						}
					}else if (k == 2) {
						ewg = get_bonds(mols[0], i + 1, 2);
						ewg2 = get_bonds(mols[0], i + 1, 1);
						if ( ewg > 0 && ewg2 > 0) {
							ewg = 0;
							for (k = 0; k < 4; ++k)
								if ( vecinos2[k] == 1)
									ewg += get_bonds(mols[0], vecinos[k] + 1, 2);

							if ( ewg > 0) {
								mols->gaff_types[i] = CG;
								if ( verboseflag == 1)
									printf("%i tipo CG.\n", i);
							}else{
								mols->gaff_types[i] = C1;
								if ( verboseflag == 1)
									printf("%i tipo C1.\n", i);
							}
						}else{
							mols->gaff_types[i] = C1;
							if ( verboseflag == 1)
								printf("%i tipo C1.\n", i);
						}
					}else if (k == 1) {
						mols->gaff_types[i] = C1;
						if ( verboseflag == 1)
							printf("%i tipo C1.\n", i);
					}else                 if ( verboseflag == 1)
						printf("Carbon Unknow type!!!\n");
				}

			}else if ( mols->atoms[i] == 2) {
				if ( k == 1) {
					mols->gaff_types[i] = O;
					if ( verboseflag == 1)
						printf("%i tipo O.\n", i);
				}else{
					flag = 0;
					for (k = 0; k < 4; ++k)
						if ( mols->atoms[vecinos[k]] == 4 )
							flag = 1;

					if ( flag == 1) {
						mols->gaff_types[i] = OH;
						if ( verboseflag == 1)
							printf("%i tipo OH.\n", i);
					}else{
						mols->gaff_types[i] = OS;
						if ( verboseflag == 1)
							printf("%i tipo OS.\n", i);
					}


				}

			}else if ( mols->atoms[i] == 4) {
				if ( k > 1)
					fprintf(stderr,"Warning: H bonded to many atoms.\n");

				for (k = 0; k < 1; ++k) {
					ewg = 0;
					bondis = get_number_any_bond(mols[0], vecinos[k] + 1, 0);
					if ( mols->atoms[vecinos[k]] == 1 ) {
						ewg += get_number_bond_by_atom(mols[0], vecinos[k] + 1, 0, 3);
						ewg += get_number_bond_by_atom(mols[0], vecinos[k] + 1, 0, 2);
						ewg += get_number_bond_by_atom(mols[0], vecinos[k] + 1, 0, 6);
						ewg += get_number_bond_by_atom(mols[0], vecinos[k] + 1, 0, 9);
						ewg += get_number_bond_by_atom(mols[0], vecinos[k] + 1, 0, 8);
						ewg += get_number_bond_by_atom(mols[0], vecinos[k] + 1, 0, 7);
						if ( bondis == 4 && ewg == 0) {
							mols->gaff_types[i] = HC;
							if ( verboseflag == 1)
								printf("%i tipo HC.\n", i);
						}else if ( bondis == 4 && ewg == 1) {
							mols->gaff_types[i] = H1;
							if ( verboseflag == 1)
								printf("%i tipo H1.\n", i);
						}else if ( bondis == 4 && ewg == 2) {
							mols->gaff_types[i] = H2;
							if ( verboseflag == 1)
								printf("%i tipo H2.\n", i);
						}else if ( bondis == 4 && ewg == 3) {
							mols->gaff_types[i] = H3;
							if ( verboseflag == 1)
								printf("%i tipo H3.\n", i);
						}else if ( bondis == 3 && ewg == 1) {
							mols->gaff_types[i] = H4;
							if ( verboseflag == 1)
								printf("%i tipo H4.\n", i);
						}else if ( bondis == 3 && ewg == 2) {
							mols->gaff_types[i] = H5;
							if ( verboseflag == 1)
								printf("%i tipo H5.\n", i);
						}else if ( bondis == 5 && ewg == 4) {
							mols->gaff_types[i] = HX;
							if ( verboseflag == 1)
								printf("%i tipo HX.\n", i);
						}else{
							mols->gaff_types[i] = HA;
							if ( verboseflag == 1)
								printf("%i tipo HA.\n", i);
						}
					}else if ( mols->atoms[vecinos[k]] == 2) {
						mols->gaff_types[i] = HO;
						if ( verboseflag == 1)
							printf("%i tipo HO.\n", i);
					}else if ( mols->atoms[vecinos[k]] == 3) {
						mols->gaff_types[i] = HN;
						if ( verboseflag == 1)
							printf("%i tipo HN.\n", i);
					}else if ( mols->atoms[vecinos[k]] == 5) {
						mols->gaff_types[i] = HP;
						if ( verboseflag == 1)
							printf("%i tipo HP.\n", i);
					}else if ( mols->atoms[vecinos[k]] == 6) {
						mols->gaff_types[i] = HS;
						if ( verboseflag == 1)
							printf("%i tipo HS.\n", i);
					}
				}

			}else if ( mols->atoms[i] == 6) {

				if ( k == 1) {
					mols->gaff_types[i] = S;
					if ( verboseflag == 1)
						printf("%i tipo S.\n", i);

				}else if ( k == 5 || k == 6) {
					mols->gaff_types[i] = S6;
					if ( verboseflag == 1)
						printf("%i tipo S6.\n", i);
				}else if ( k == 2) {
					flag = 0;
					for (j = 0; j < 4; ++j)
						if ( mols->atoms[vecinos[j]] == 4 )
							flag = 1;
					if ( flag == 1) {
						mols->gaff_types[i] = SH;
						if ( verboseflag == 1)
							printf("%i tipo SH.\n", i);
					}else if ( get_bonds(mols[0], i + 1, 2) > 0 || get_bonds(mols[0], i + 1, 3) > 0) {
						mols->gaff_types[i] = S2;
						if ( verboseflag == 1)
							printf("%i tipo S2.\n", i);
					}else{
						mols->gaff_types[i] = SS;
						if ( verboseflag == 1)
							printf("%i tipo SS.\n", i);
					}
				}else if ( k == 3) {
					ewg = get_bonds(mols[0], i + 1, 2);
					ewg2 = get_bonds(mols[0], i + 1, 1);
					if ( ewg > 0 && ewg2 > 0) {
						ewg = 0;
						for (k = 0; k < 4; ++k)
							if ( vecinos2[k] == 1)
								ewg += get_bonds(mols[0], vecinos[k] + 1, 2);
						if ( ewg > 0) {
							mols->gaff_types[i] = SX;
							if ( verboseflag == 1)
								printf("%i tipo SX.\n", i);
						}else{
							mols->gaff_types[i] = S4;
							if ( verboseflag == 1)
								printf("%i tipo S4.\n", i);
						}
					}else{
						mols->gaff_types[i] = S4;
						if ( verboseflag == 1)
							printf("%i tipo S4.\n", i);
					}
				}else if ( k == 4) {
					ewg = get_bonds(mols[0], i + 1, 2);
					ewg2 = get_bonds(mols[0], i + 1, 1);
					if ( ewg > 0 && ewg2 > 0) {
						ewg = 0;
						for (k = 0; k < 4; ++k)
							if ( vecinos2[k] == 1)
								ewg += get_bonds(mols[0], vecinos[k] + 1, 2);
						if ( ewg > 0) {
							mols->gaff_types[i] = SX;
							if ( verboseflag == 1)
								printf("%i tipo SY.\n", i);
						}else{
							mols->gaff_types[i] = S4;
							if ( verboseflag == 1)
								printf("%i tipo S6.\n", i);
						}
					}else{
						mols->gaff_types[i] = S4;
						if ( verboseflag == 1)
							printf("%i tipo S6.\n", i);
					}
				}else{

					if ( verboseflag == 1)
						printf("Unkown type of S.\n");
				}


			}else if ( mols->atoms[i] == 5) {
				if ( k == 1) {
					mols->gaff_types[i] = P2;
					if ( verboseflag == 1)
						printf("%i tipo P2.\n", i);
				}else if ( k == 2) {
					ewg = get_bonds(mols[0], i + 1, 2);
					ewg2 = get_bonds(mols[0], i + 1, 1);
					if ( ewg > 0 && ewg2 > 0) {
						ewg = 0;
						for (k = 0; k < 4; ++k)
							if ( vecinos2[k] == 1)
								ewg += get_bonds(mols[0], vecinos[k] + 1, 2);
						if ( ewg > 0) {

							if ( mols->ringer[i] != 0) {
								mols->gaff_types[i] = PC;
								if ( verboseflag == 1)
									printf("%i tipo PC.\n", i);
							}else{
								mols->gaff_types[i] = PE;
								if ( verboseflag == 1)
									printf("%i tipo PE.\n", i);
							}
						}else{
							mols->gaff_types[i] = P2;
							if ( verboseflag == 1)
								printf("%i tipo P2.\n", i);
						}
					}else{
						mols->gaff_types[i] = P2;
						if ( verboseflag == 1)
							printf("%i tipo P2.\n", i);
					}
				}else if ( k == 3) {
					ewg = 0;
					ewg += get_number_bond_by_atom(mols[0], i + 1, 2, 2);
					ewg += get_number_bond_by_atom(mols[0], i + 1, 2, 6);
					if ( ewg != 0) {
						mols->gaff_types[i] = P4;
						if ( verboseflag == 1)
							printf("%i tipo P4.\n", i);
					}else{
						ewg = get_bonds(mols[0], i + 1, 2);
						ewg2 = get_bonds(mols[0], i + 1, 1);
						if ( ewg > 0 && ewg2 > 0) {
							ewg = 0;
							for (k = 0; k < 4; ++k)
								if ( vecinos2[k] == 1)
									ewg += get_bonds(mols[0], vecinos[k] + 1, 2);
							if ( ewg > 0) {
								mols->gaff_types[i] = PX;
								if ( verboseflag == 1)
									printf("%i tipo PX.\n", i);
							}else{
								mols->gaff_types[i] = P3;
								if ( verboseflag == 1)
									printf("%i tipo P3.\n", i);
							}
						}else{
							mols->gaff_types[i] = P3;
							if ( verboseflag == 1)
								printf("%i tipo P3.\n", i);
						}
					}
				}else if ( k == 4) {
					ewg = get_bonds(mols[0], i + 1, 2);
					ewg2 = get_bonds(mols[0], i + 1, 1);
					if ( ewg > 0 && ewg2 > 0) {
						ewg = 0;
						for (k = 0; k < 4; ++k)
							if ( vecinos2[k] == 1)
								ewg += get_bonds(mols[0], vecinos[k] + 1, 2);
						if ( ewg > 0) {
							mols->gaff_types[i] = PY;
							if ( verboseflag == 1)
								printf("%i tipo PY.\n", i);
						}else{
							mols->gaff_types[i] = P5;
							if ( verboseflag == 1)
								printf("%i tipo P5.\n", i);
						}
					}else{
						mols->gaff_types[i] = P5;
						if ( verboseflag == 1)
							printf("%i tipo P5.\n", i);
					}
				}else if ( k == 5 || k == 6) {
					mols->gaff_types[i] = P5;
					if ( verboseflag == 1)
						printf("%i tipo P5.\n", i);
				}else if ( verboseflag == 1)
					printf("P Unknow type.\n");
			}else if ( mols->atoms[i] == 10 ) {
				mols->gaff_types[i] = F;
				if ( verboseflag == 1)
					printf("%i tipo F.\n", i);
			}else if ( mols->atoms[i] == 9 ) {
				mols->gaff_types[i] = CL;
				if ( verboseflag == 1)
					printf("%i tipo CL.\n", i);
			}else if ( mols->atoms[i] == 8 ) {
				mols->gaff_types[i] = BR;
				if ( verboseflag == 1)
					printf("%i tipo BR.\n", i);
			}else if ( mols->atoms[i] == 7 ) {
				mols->gaff_types[i] = I;
				if ( verboseflag == 1)
					printf("%i tipo I.\n", i);
			}
		}


/* Another loop for N. I need carbonyl checks before typing amides*/
/* Lazzy-crap(tm) technique */

		for (i = 0; i < mols->n_atoms; ++i) {
/* Get vecinos */
			for (j = 0; j < mols->n_atoms; ++j) {
				vecinos[j] = 0;
				vecinos2[j] = 0;
			}
			k = 0;
			for (j = 0; j < mols->n_bonds; ++j) {
				if ( mols->bond_a1[j] == (i + 1)) {
					vecinos[k] = mols->bond_a2[j] - 1;
					vecinos2[k] = mols->bonds[j];
					k++;
				}else if ( mols->bond_a2[j] == (i + 1)) {
					vecinos[k] = mols->bond_a1[j] - 1;
					vecinos2[k] = mols->bonds[j];
					k++;
				}
			}

			if ( mols->atoms[i] == 3 ) {
				if ( k == 1) {
					mols->gaff_types[i] = N1;
					if ( verboseflag == 1)
						printf("%i tipo N1.\n", i);
				}else if ( k == 2) {
					ewg = get_bonds(mols[0], i + 1, 2);
					if ( mols->aromatic[i] == 1) {
						mols->gaff_types[i] = NB;
						if ( verboseflag == 1)
							printf("%i tipo NB.\n", i);
					}else if ( get_bonds(mols[0], i + 1, 3)  > 0 || ewg == 2) {
						mols->gaff_types[i] = N1;
						if ( verboseflag == 1)
							printf("%i tipo N1.\n", i);
					}else{
						ewg2 = get_bonds(mols[0], i + 1, 1);
						if ( ewg > 0 && ewg2 > 0) {
							ewg = 0;
							for (k = 0; k < 4; ++k)
								if ( vecinos2[k] == 1)
									ewg += get_bonds(mols[0], vecinos[k] + 1, 2);
							if ( ewg > 0) {
								if ( mols->ringer[i] != 0) {
									mols->gaff_types[i] = NC;
									if ( verboseflag == 1)
										printf("%i tipo NC.\n", i);
								}else{
									mols->gaff_types[i] = NE;
									if ( verboseflag == 1)
										printf("%i tipo NE.\n", i);
								}
							}else{
								mols->gaff_types[i] = N2;
								if ( verboseflag == 1)
									printf("%i tipo N2.\n", i);
							}
						}else{
							mols->gaff_types[i] = N2;
							if ( verboseflag == 1)
								printf("%i tipo N2.\n", i);
						}
					}
				}else if ( k == 3) {
					if ( mols->aromatic[i] == 1) {
						mols->gaff_types[i] = NA;
						if ( verboseflag == 1)
							printf("%i tipo NA.\n", i);
					}else{

						flag2 = 0;
						flag = 0;
						flag3 = 0;
						for ( j = 0; j < k; ++j) {
							if ( mols->atoms[vecinos[j]] == 2 )
								flag2++;
							else if ( mols->gaff_types[vecinos[j]] == C)
								flag = 1;
							else if ( vecinos2[j] == 2 && (mols->gaff_types[vecinos[j]] == C || mols->atoms[vecinos[j]] == 3 || mols->atoms[vecinos[j]] == 1 ))
								flag3++;
						}

						if ( flag2 == 2 && k == 3) {
							mols->gaff_types[i] = NO;
							if ( verboseflag == 1)
								printf("%i tipo NO.\n", i);
						}else{
							if ( flag == 1) {
								mols->gaff_types[i] = N;
								if ( verboseflag == 1)
									printf("%i tipo N.\n", i);
							}else if ( flag3 == 1) {

								mols->gaff_types[i] = NH;
								if ( verboseflag == 1)
									printf("%i tipo NH.\n", i);

							}else{
								mols->gaff_types[i] = N3;
								if ( verboseflag == 1)
									printf("%i tipo N3.\n", i);
							}
						}
					}

				}else if ( k == 4) {
					mols->gaff_types[i] = N4;
					if ( verboseflag == 1)
						printf("%i tipo N4.\n", i);
				}else if ( verboseflag == 1)
					printf("N  Unknow type.\n");
			}
		}


#ifdef DEBUG
		for (j = 0; j < mols->n_atoms; ++j)
			printf("%i - %i.\n", j, mols->gaff_types[j]);

#endif



		if ( angleflag == 1) {
			cur_angle = 0;
			mols->ia = (int*)realloc(mols->ia, sizeof(int) * mols->n_bonds * mols->n_bonds);
			mols->ja = (int*)realloc(mols->ja, sizeof(int) * mols->n_bonds * mols->n_bonds);
			mols->ka = (int*)realloc(mols->ka, sizeof(int) * mols->n_bonds * mols->n_bonds);
			bonds1 = bonds2 = 0;

			for ( i = 0; i < mols->n_atoms; ++i) {
				bonds1 = 0;
				for (step = 0; step < 4; step++)
					vecinos[step] = 0;

				for (step = 0; step < mols->n_bonds; ++step) {
					if ( mols->bond_a1[step] == (i + 1)) {
						vecinos[bonds1] = mols->bond_a2[step] - 1;
						bonds1++;
					}else if ( mols->bond_a2[step] == (i + 1)) {
						vecinos[bonds1] = mols->bond_a1[step] - 1;
						bonds1++;
					}
				}
				for ( j = 0; j < bonds1; ++j) {
					bonds2 = 0;
					for (step = 0; step < 4; step++)
						vecinos2[step] = 0;

					for (step = 0; step < mols->n_bonds; ++step) {
						if ( mols->bond_a1[step] == (vecinos[j] + 1)) {
							vecinos2[bonds2] = mols->bond_a2[step] - 1;
							bonds2++;
						}else if ( mols->bond_a2[step] == (vecinos[j] + 1)) {
							vecinos2[bonds2] = mols->bond_a1[step] - 1;
							bonds2++;
						}
					}

					for ( l = 0; l < bonds2; l++) {


						if ( i != vecinos[j] && i != vecinos2[l] && vecinos[j] != vecinos2[l]  && check_angle(&mols, cur_angle, i, vecinos[j], vecinos2[l]) == 0  ) {
							mols->ia[cur_angle] = i;
							mols->ja[cur_angle] = vecinos[j];
							mols->ka[cur_angle] = vecinos2[l];
							cur_angle++;
						}
					}


				}


			}
			if ( verboseflag == 1)
				printf("%i angles found.\n", cur_angle);

		}

		mols->n_angles = cur_angle;

		if ( torflag == 1) {

			if ( verboseflag == 1)
				printf("Allocating & finding impropers...\n");

			mols->n_impropers = 0;
			mols->ip = (int*)realloc(mols->ip, sizeof(int) * mols->n_bonds * 6);
			mols->jp = (int*)realloc(mols->jp, sizeof(int) * mols->n_bonds * 6);
			mols->kp = (int*)realloc(mols->kp, sizeof(int) * mols->n_bonds * 6);
			mols->lp = (int*)realloc(mols->lp, sizeof(int) * mols->n_bonds * 6);
			mols->tp = (int*)realloc(mols->tp, sizeof(int) * mols->n_bonds * 6);


			for ( i = 0; i < mols->n_atoms; ++i) {
				if ( verboseflag == 1)
					printf("%i/%i\r", i + 1, mols->n_atoms);

				for (j = 0; j < 40; ++j) {
					if ( mols->gaff_types[i] == impropers[j][2]) { /* Look for central atom */
						for (l = 0; l < 8; ++l) {
							vecinos[l] = 0; vecinos2[l] = 0;
						}

						k = 0;
						for (l = 0; l < mols->n_bonds; ++l) {
							if ( mols->bond_a1[l] == (i + 1)) {
								vecinos[k] = mols->gaff_types[mols->bond_a2[l] - 1];
								vecinos2[k] = mols->bond_a2[l] - 1;
								k++;
							}else if ( mols->bond_a2[l] == (i + 1)) {
								vecinos[k] = mols->gaff_types[mols->bond_a1[l] - 1];
								vecinos2[k] = mols->bond_a1[l] - 1;
								k++;
							}
						}

						nmatchs = 0;
						matchs[0] = 0;
						matchs[1] = 0;
						matchs[2] = 0;

						for ( l = 0; l < k; ++l) {
							/* 71 = wildcard */
							if (  vecinos[l] == impropers[j][0]) {
								matchs[0] = vecinos2[l];
								nmatchs++;
								break;
							}
						}

						for ( l = 0; l < k; ++l) {
							/* 71 = wildcard */
							if (  (vecinos[l] == impropers[j][1] ) && vecinos2[l] != matchs[0] ) {
								matchs[1] = vecinos2[l];
								nmatchs++;
								break;
							}
						}

						for ( l = 0; l < k; ++l) {
							/* 71 = wildcard */
							if (  (vecinos[l] == impropers[j][3] ) && vecinos2[l] != matchs[1] && vecinos2[l] != matchs[1] ) {
								matchs[2] = vecinos2[l];
								nmatchs++;
								break;
							}
						}


						if ( impropers[j][0] == 71) {
							if ( vecinos2[0] != matchs[1] && vecinos2[0] != matchs[2]) {
								matchs[0] = vecinos2[0];
								nmatchs++;
							}else if ( vecinos2[1] != matchs[1] && vecinos2[1] != matchs[2]) {
								matchs[0] = vecinos2[1];
								nmatchs++;
							}else if ( vecinos2[2] != matchs[1] && vecinos2[2] != matchs[2]) {
								matchs[0] = vecinos2[2];
								nmatchs++;
							}else if ( vecinos2[3] != matchs[1] && vecinos2[3] != matchs[2]) {
								matchs[0] = vecinos2[3];
								nmatchs++;
							}

						}


						if ( impropers[j][1] == 71) {
							if ( vecinos2[0] != matchs[0] && vecinos2[0] != matchs[2]) {
								matchs[1] = vecinos2[0];
								nmatchs++;
							}else if ( vecinos2[1] != matchs[0] && vecinos2[1] != matchs[2]) {
								matchs[1] = vecinos2[1];
								nmatchs++;
							}else if ( vecinos2[2] != matchs[0] && vecinos2[2] != matchs[2]) {
								matchs[1] = vecinos2[2];
								nmatchs++;
							}else if ( vecinos2[3] != matchs[0] && vecinos2[3] != matchs[2]) {
								matchs[1] = vecinos2[3];
								nmatchs++;
							}

						}



						if ( impropers[j][3] == 71) {
							if ( vecinos2[0] != matchs[0] && vecinos2[0] != matchs[1]) {
								matchs[2] = vecinos2[0];
								nmatchs++;
							}else if ( vecinos2[1] != matchs[0] && vecinos2[1] != matchs[1]) {
								matchs[2] = vecinos2[1];
								nmatchs++;
							}else if ( vecinos2[2] != matchs[0] && vecinos2[2] != matchs[1]) {
								matchs[2] = vecinos2[2];
								nmatchs++;
							}else if ( vecinos2[3] != matchs[0] && vecinos2[3] != matchs[1]) {
								matchs[2] = vecinos2[3];
								nmatchs++;
							}

						}



						if ( nmatchs == 3) { /* Bingo! We have an improper */

							mols->ip[mols->n_impropers] = matchs[0];
							mols->jp[mols->n_impropers] = matchs[1];
							mols->kp[mols->n_impropers] = i;
							mols->lp[mols->n_impropers] = matchs[2];
							mols->tp[mols->n_impropers] = j;
							mols->n_impropers = mols->n_impropers + 1;

							if ( verboseflag == 1)
								printf("\nImproper: %i-%i-%i-%i tipo %i.\n", matchs[0], matchs[1], i, matchs[2], j);
						}

					}
				}
			}

			if ( verboseflag == 1)
				printf("%i impropers found.\n", mols->n_impropers);

			mols->ik = (int*)realloc(mols->ik, sizeof(int) * mols->n_bonds * 6);
			mols->jk = (int*)realloc(mols->jk, sizeof(int) * mols->n_bonds * 6);
			mols->kk = (int*)realloc(mols->kk, sizeof(int) * mols->n_bonds * 6);
			mols->lk = (int*)realloc(mols->lk, sizeof(int) * mols->n_bonds * 6);

			ttor = 0;
			bonds1 = bonds2 = bonds3 = 0;
			for ( i = 0; i < mols->n_atoms; ++i) {
				bonds1 = 0;
				for (step = 0; step < 4; step++)
					vecinos[step] = 0;

				for (step = 0; step < mols->n_bonds; ++step) {
					if ( mols->bond_a1[step] == (i + 1)) {
						vecinos[bonds1] = mols->bond_a2[step] - 1;
						bonds1++;
					}else if ( mols->bond_a2[step] == (i + 1)) {
						vecinos[bonds1] = mols->bond_a1[step] - 1;
						bonds1++;
					}
				}
				for ( j = 0; j < bonds1; ++j) {
					bonds2 = 0;
					for (step = 0; step < 4; step++)
						vecinos2[step] = 0;

					for (step = 0; step < mols->n_bonds; ++step) {
						if ( mols->bond_a1[step] == (vecinos[j] + 1)) {
							vecinos2[bonds2] = mols->bond_a2[step] - 1;
							bonds2++;
						}else if ( mols->bond_a2[step] == (vecinos[j] + 1)) {
							vecinos2[bonds2] = mols->bond_a1[step] - 1;
							bonds2++;
						}
					}

					for ( k = 0; k < bonds2; ++k) {

						for (step = 0; step < 4; step++)
							vecinos4[step] = 0;

						bonds3 = 0;
						for (step = 0; step < mols->n_bonds; ++step) {
							if ( mols->bond_a1[step] == (vecinos2[k] + 1)) {
								vecinos4[bonds3] = mols->bond_a2[step] - 1;
								bonds3++;
							}else if ( mols->bond_a2[step] == (vecinos2[k] + 1)) {
								vecinos4[bonds3] = mols->bond_a1[step] - 1;
								bonds3++;
							}
						}

						for ( l = 0; l < bonds3; l++) {


							if ( i != vecinos[j] && i != vecinos2[k] && i != vecinos4[l] && vecinos[j] != vecinos2[k] && vecinos[j] != vecinos4[l] && vecinos2[k] != vecinos4[l] && check_dihedral(&mols, ttor, i, vecinos[j], vecinos2[k], vecinos4[l]) == 0  ) {

								mols->ik[ttor] = i;
								mols->jk[ttor] = vecinos[j];
								mols->kk[ttor] = vecinos2[k];
								mols->lk[ttor] = vecinos4[l];
								ttor++;
							}
						}
					}
				}
			}

			if ( verboseflag == 1)
				printf("%i dihedrals found.\n", ttor);

			mols->n_torsionals = ttor;
		}


/* VDW Pairs */


		if ( vdwflag == 1) {
			cpair = 0;
			mols->vdwpairs_a = (int*)realloc(mols->vdwpairs_a, sizeof(int) * mols->n_atoms * mols->n_atoms);
			mols->vdwpairs_b = (int*)realloc(mols->vdwpairs_b, sizeof(int) * mols->n_atoms * mols->n_atoms);
			mols->vdw_type = (int*)realloc(mols->vdw_type, sizeof(int) * mols->n_atoms * mols->n_atoms);

			if ( verboseflag == 1) {
				printf("Allocated %i for pairs.\n", mols->n_atoms * mols->n_atoms);
				printf("Setting up van der Waals pairs.\n\n");
			}
			for ( j = 0; j < mols->n_atoms; ++j) {

				if ( verboseflag == 1)
					printf("\r%i/%i", j + 1, mols->n_atoms);

				for ( k = j + 1; k < mols->n_atoms; ++k) {
					dx = mols->x[j] - mols->x[k];
					dy = mols->y[j] - mols->y[k];
					dz = mols->z[j] - mols->z[k];
					dist = dx * dx + dy * dy + dz * dz;
					if ( dist <= 64) {
						if ( get_number_any_bond(mols[0], j + 1, k + 1) == 0) {
							if ( q13_bond(mols, j, k) == 0) {
								mols->vdwpairs_a[cpair] = j;
								mols->vdwpairs_b[cpair] = k;

								if ( q14_bond(mols, j, k) != 0)
									mols->vdw_type[cpair] = 2;
								else
									mols->vdw_type[cpair] = 1;

								cpair++;
							}
						}
					}
				}
			}

			mols->n_pairs = cpair;
			if ( verboseflag == 1)
				printf("Total of VDW pairs: %i\n", cpair);
		}


	}

	free(vecinos);
	free(vecinos2);

	*mymol = mols;

}

/**
 *
 *	@brief Set the default conformer of the MOL2 structure to the given one
 *	@param mymol MOL2 structure of the molecule
 *	@param conformer number of the conformer starting at 0
 *	@return 0 on success
 */
int set_default_conformer(MOL2 **mymol, int conformer)
{

	MOL2 *mols = NULL;
	int i = 0;

	mols = *mymol;

	if ( mols->nconformers <= conformer)
		return -1;


	for ( i = 0; i < mols->n_atoms; ++i) {
		mols->x[i] = mols->conformers[conformer].x[i];
		mols->y[i] = mols->conformers[conformer].y[i];
		mols->z[i] = mols->conformers[conformer].z[i];
		mols->pcharges[i] = mols->conformers[conformer].pcharges[i];
		mols->radius[i] = mols->conformers[conformer].radius[i];

	}

	mols->default_conformer = conformer;


	*mymol = mols;

	return 0;
}

/**
 *
 *	@brief Syncronizes the coordinates of the current default conformer with the PDB set
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *	@param mymol MOL2 structure of the molecule
 *
 */
int update_default_conformer(MOL2 **mymol)
{

	MOL2 *mols = NULL;
	int i = 0;

	mols = *mymol;

	for ( i = 0; i < mols->n_atoms; ++i) {
		mols->conformers[mols->default_conformer].x[i] = mols->x[i];
		mols->conformers[mols->default_conformer].y[i] = mols->y[i];
		mols->conformers[mols->default_conformer].z[i] = mols->z[i];
	}

	*mymol = mols;

	return 0;
}

int MultiMOL2_reader(MOL2 **mymol, char *finput_name)
{

	MOL2 *mol = NULL;
	char *line = NULL;
	FILE *input = NULL;


        if ( (input = fopen(finput_name, "r")) == NULL) {
                fprintf(stderr, "Error. Cant open file %s.\n", finput_name);
                fflush(stderr);
                return(-1);
        }


        if ( (line = (char*)calloc(1, MAX_BUFFER + 2)) == NULL) {
                fprintf(stderr, "Error. Cant al memory.\n");
                fflush(stderr);
		fclose(input);
                return(-2);
        }


	fclose(input);
	free(line);
}
