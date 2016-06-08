/**
 *
 *      @file groups.c
 *      @brief Group of internal function for the library
 *
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @date 01/10/2010
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_RING_MEMBERS 800    /* Puff! but simple linked lists = hell! */
#define MAX_RINGS 30            /* Little hack to prevent infinite recursive problem */

/**
 *  String to print for each atom type
 */
char *elems[10] = { "C ", "O ", "N ", "H ", "P ", "S ", "I ", "Br", "Cl", "F " };


/**
 *	@brief Transformate the molecule according to a rigid body movement contained in the simplex (traslation + rotation)
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *
 *	@param mymol MOL2 structure which contains the ligand
 *	@param simplex vector of 6 values with the traslation (first three values are x,y,z) and rotation (euler angles in degrees are last three components)
 *	@return 0 on success
 *
 */
int transformate_mol(MOL2 **mymol, float *simplex)
{

	int i = 0, j = 0, real_atoms = 0;

	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;
	MOL2 *mol = NULL;

	float ith = 0, iph = 0, ips = 0;
	float theta = 0.0f, phi = 0.0f, psi = 0.0f;
	float sthe = 0.0f, cthe = 0.0f, sphi = 0.0f, cphi = 0.0f, spsi = 0.0f, cpsi = 0.0f;
	float A[3][3], dx = 0.0f, dy = 0.0f, dz = 0.0f;
	const float cdr = 0.0174532777778;

	float comt[3], com[3], Xrot = 0.0f, Yrot = 0.0f, Zrot = 0.0f, VX = 0.0f, VY = 0.0f, VZ = 0.0f;
	int *tor_type = NULL;
	double energy = 0.0;

	mol = *mymol;

	for ( i = 0; i < mol->n_atoms; ++i)
		if ( mol->backbone[i] != 1)
			++real_atoms;

	x = simplex[0];
	y = simplex[1];
	z = simplex[2];

	ith = simplex[3];
	iph = simplex[4];
	ips = simplex[5];


/**
 *
 *	Transform to radians and calculate sines and cosines for the three angles
 *
 */
	theta = cdr * ith;
	phi = cdr * iph;
	psi = cdr * ips;
	sthe = sin(theta);
	cthe = cos(theta);
	sphi = sin(phi);
	cphi = cos(phi);
	spsi = sin(psi);
	cpsi = cos(psi);

/**
 *
 *	Build rotation matrix from euler angles
 *
 */
	A[0][0] = -sphi * spsi + cthe * cphi * cpsi;
	A[0][1] =  cphi * spsi + cthe * sphi * cpsi;
	A[0][2] = -sthe * cpsi;
	A[1][0] = -sphi * cpsi - cthe * cphi * spsi;
	A[1][1] =  cphi * cpsi - cthe * sphi * spsi;
	A[1][2] =  sthe * spsi;
	A[2][0] =  sthe * cphi;
	A[2][1] =  sthe * sphi;
	A[2][2] =  cthe;


/**
 *
 *	Center of mass calculation
 *
 */

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


/**
 *
 *	Rotation and translation
 *
 */
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


	*mymol = mol;
	return 0;
}


/**
 *
 *	@brief Load rotamers from the Richardson library
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param type Residue type
 *	@param max maximum number of rotamers
 *	@param tmp_rota Structure to load rotamers in
 *	@return 0 on success
 *
 */
int get_romaters(int type, int max, ROTAMER **tmp_rota)
{
	int i = 0;
	ROTAMER *rot = NULL;
	FILE *lib = NULL;
	char *lib_name;
	char *line = NULL;
	char myx[12], myy[12], myz[12];
	rot = *tmp_rota;
	rot->n_atoms = aasize[ type - 1];
	rot->n_rotamers = max;

	lib_name = (char*)calloc(sizeof(char), 100);

	switch (type) {
	case 1: strcpy(lib_name,  "lib/thr.pdb"); break;
	case 2: strcpy(lib_name,  "lib/val.pdb"); break;
	case 3: strcpy(lib_name,  "lib/ser.pdb"); break;
	case 4: strcpy(lib_name,  "lib/his.pdb"); break;
	case 5: strcpy(lib_name,  "lib/lys.pdb"); break;
	case 6: strcpy(lib_name,  "lib/arg.pdb"); break;
	case 7: strcpy(lib_name,  "lib/met.pdb"); break;
	case 8: strcpy(lib_name,  "lib/cys.pdb"); break;
	case 9: strcpy(lib_name,  "lib/glu.pdb"); break;
	case 10: strcpy(lib_name, "lib/gln.pdb"); break;
	case 11: strcpy(lib_name,  "lib/asp.pdb"); break;
	case 12: strcpy(lib_name,  "lib/asn.pdb"); break;
	case 13: strcpy(lib_name,  "lib/phe.pdb"); break;
	case 14: strcpy(lib_name,  "lib/tyr.pdb"); break;
	case 15: strcpy(lib_name,  "lib/trp.pdb"); break;
	case 16: strcpy(lib_name,  "lib/pro.pdb"); break;
	case 17: strcpy(lib_name,  "lib/gly.pdb"); break;
	case 18: strcpy(lib_name,  "lib/ala.pdb"); break;
	case 19: strcpy(lib_name,  "lib/leu.pdb"); break;
	case 20: strcpy(lib_name,  "lib/ile.pdb"); break;
	}

	printf("Library is: %s.\n", lib_name);

	lib = fopen(lib_name, "rb");
	if ( lib == NULL) {
		printf("Error. Cant open rotamer library.\n");
		return -1;
	}

	rot->x = (float*)calloc(sizeof(float), (max * rot->n_atoms) + 10);
	rot->y = (float*)calloc(sizeof(float), (max * rot->n_atoms) + 10);
	rot->z = (float*)calloc(sizeof(float), (max * rot->n_atoms) + 10);
	line = (char*)calloc(sizeof(char), 1024);

	printf("Max is: %i\n", max * rot->n_atoms);
	for ( i = 0; i < max * rot->n_atoms; ++i) {
		fgets(line, 1024, lib);
		strncpy(myx, &line[29], 10);
		strncpy(myy, &line[38], 10);
		strncpy(myz, &line[46], 10);

		rot->x[i] = atof(myx);
		rot->y[i] = atof(myy);
		rot->z[i] = atof(myz);

	}
	fclose(lib);
	free(line);
	free(lib_name);


}


/**
 *
 *	@brief Superimpose library romaters on current aminoacid
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 structure with the protein
 *	@param tmp_r residue superimposed
 *	@param target1 Coordinates of CA atom
 *	@param target2 Coordinates of CB atom
 *	@param ca2 Index of the CA atom
 *	@param cb Index of the CB atom
 *	@return 0 on success
 *
 */
int super_transform(MOL2 *mol, RESIDUE *tmp_r, float target1[3], float target2[3], int ca2, int cb)
{
	RESIDUE *res;
	float ca;
	float ref[2][3];
	float from[3], to[3];
	float v[3], vs[3], vt[3];
	float norm = 0.0;
	float dx, dy, dz;

	res = tmp_r;

	/* CA ? Maybe. Depends on nitrogen protonation state. Dont mind at all.*/
	ref[0][0] = mol->x[ca2];
	ref[0][1] = mol->y[ca2];
	ref[0][2] = mol->z[ca2];
	/* CB  Same shity sutff */
	ref[1][0] = mol->x[cb];
	ref[1][1] = mol->y[cb];
	ref[1][2] = mol->z[cb];

	dx = ref[0][0];
	dy = ref[0][1];
	dz = ref[0][2];

	to[0] = ref[0][0] - ref[1][0];
	to[1] = ref[0][1] - ref[1][1];
	to[2] = ref[0][2] - ref[1][2];

	norm = sqrt(to[0] * to[0] + to[1] * to[1] + to[2] * to[2]);
	to[0] /= norm;
	to[1] /= norm;
	to[2] /= norm;


	from[0] = target1[0] - target2[0];
	from[1] = target1[1] - target2[1];
	from[2] = target1[2] - target2[2];

	norm = sqrt(from[0] * from[0] + from[1] * from[1] + from[2] * from[2]);

	from[0] /= norm;
	from[1] /= norm;
	from[2] /= norm;


	vs[0] = from[1] * to[2] - from[2] * to[1];
	vs[1] = (-from[0] * to[2] + from[2] * to[0]);
	vs[2] = from[0] * to[1] - from[1] * to[0];

	v[0] = vs[0];
	v[1] = vs[1];
	v[2] = vs[2];

	norm = sqrt( v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;

	ca = from[0] * to[0] + from[1] * to[1] + from[2] * to[2];

	vt[0] = v[0] * (1.0f - ca);
	vt[1] = v[1] * (1.0f - ca);
	vt[2] = v[2] * (1.0f - ca);

	res->rotM[0][0] = vt[0] * v[0] + ca;
	res->rotM[1][1] = vt[1] * v[1] + ca;
	res->rotM[2][2] = vt[2] * v[2] + ca;

	vt[0] *= v[1];
	vt[2] *= v[0];
	vt[1] *= v[2];

	res->rotM[0][1] = vt[0] - vs[2];
	res->rotM[0][2] = vt[2] + vs[1];
	res->rotM[1][0] = vt[0] + vs[2];
	res->rotM[1][2] = vt[1] - vs[0];
	res->rotM[2][0] = vt[2] - vs[1];
	res->rotM[2][1] = vt[1] + vs[0];

	res->transvector[0] = dx;
	res->transvector[1] = dy;
	res->transvector[2] = dz;


	return 0;
}

/**
 *
 *	@brief Get the index of a residue name
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param tmp_resn Resiude name
 *	@return index of the residue type
 *
 */
int name_to_number_residue(char tmp_resn[10])
{

	int i;

	for (i = 0; i < 20; ++i)
		if ( residues[i][0]  == tmp_resn[0] && residues[i][1]  == tmp_resn[1] && residues[i][2]  == tmp_resn[2] )
			return i + 1;

	return -1;
}

/**
 *
 *	@brief Check if the residue is defined in the set
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param resn Set of residue indices
 *	@param total_res total number of residues defined
 *	@param tmp_resn_i index of the residue to check
 *	@return index on success. -1 on failure
 *
 */
int check_residue(int *resn, int total_res, int tmp_resn_i)
{
	int i;

	for ( i = 0; i < total_res; ++i)
		if ( resn[i] == tmp_resn_i)
			return i;

	return -1;
}

/**
 *
 *	@brief Check if three atoms form an angle
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mmol MOL2 with the molecule
 *	@param tor number of angles in the molecule
 *	@param i atom i
 *	@param j atom j
 *	@param k atom k
 *	@return 1 on success or 0 on failure.
 *
 */
int check_angle(MOL2 **mmol, int tor, int i, int j, int k)
{
	MOL2 *mol;
	int step = 0;

	mol = *mmol;

	for ( step = 0; step < tor; ++step) {
		if ( (mol->ia[step] == i && mol->ja[step] == j && mol->ka[step] == k) || (mol->ia[step] == k && mol->ja[step] == j && mol->ka[step] == i) )
			return 1;

	}

	return 0;
}

/**
 *
 *	@brief Check if four atoms form a dihedral
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mmol MOL2 with the molecule
 *	@param tor number of torsionals in the molecule
 *	@param i atom i
 *	@param j atom j
 *	@param k atom k
 *	@param l atom l
 *	@return 1 on success. 0 on failure
 *
 */
int check_dihedral(MOL2 **mmol, int tor, int i, int j, int k, int l)
{
	MOL2 *mol;

	int step = 0;
	mol = *mmol;

	for ( step = 0; step < tor; ++step) {
		if ( (mol->ik[step] == i && mol->jk[step] == j && mol->kk[step] == k && mol->lk[step] == l) || (mol->ik[step] == l && mol->jk[step] == k && mol->kk[step] == j && mol->lk[step] == i) )
			return 1;

	}
	return 0;
}

/**
 *
 *	@brief Import the topology (types, bonds, angles, torsionals and impropers) from top files
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mymol MOL2 to add the topology
 *	@return 0 on success.
 *
 */
int import_topology(MOL2 **mymol)
{
	MOL2 *mols;
	FILE *in1 = NULL;
	char *line = NULL;
	int n_atoms = 0;
	int n_angles = 0;
	int n_dihedrals = 0;
	int n_bonds = 0;
	int n_pairs = 0;
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int m = 0;

	mols = *mymol;
	line = (char*)calloc(sizeof(char), 1024);

	in1 = fopen("types.top", "rb");
	if ( in1 == NULL) {
		free(line);
		return -1;
	}
	while ( fgets(line, 1024, in1))
		++n_atoms;
	rewind(in1);

	printf("Old atoms: %i. New atoms: %i.\n", mols->n_atoms, n_atoms);
	n_atoms = 0;
	while ( fgets(line, 1024, in1)) {
		sscanf(line, "%*i %i", &i);
		mols->gaff_types[n_atoms] = i;

		n_atoms++;
	}
	fclose(in1);

	n_atoms = 0;
	in1 = fopen("bonds.top", "rb");
	if ( in1 == NULL) {
		free(line);
		return -1;
	}
	while ( fgets(line, 1024, in1))
		++n_atoms;
	rewind(in1);

	mols->n_bonds = n_atoms;
	mols->bonds = (int*)calloc(sizeof(int), mols->n_bonds);
	mols->bond_a1 = (int*)calloc(sizeof(int), mols->n_bonds);
	mols->bond_a2 = (int*)calloc(sizeof(int), mols->n_bonds);

	n_atoms = 0;
	while ( fgets(line, 1024, in1)) {
		sscanf(line, "%i %i %i", &i, &j, &k);
		mols->bond_a1[n_atoms] = i;
		mols->bond_a2[n_atoms] = j;
		mols->bonds[n_atoms] = k;

		n_atoms++;
	}
	fclose(in1);
	printf("Bonds in.\n");
	n_atoms = 0;
	in1 = fopen("angles.top", "rb");
	if ( in1 == NULL) {
		free(line);
		return -1;
	}
	while ( fgets(line, 1024, in1))
		++n_atoms;
	rewind(in1);

	mols->n_angles = n_atoms;
	mols->ia = (int*)calloc( sizeof(int), mols->n_angles);
	mols->ja = (int*)calloc( sizeof(int), mols->n_angles);
	mols->ka = (int*)calloc( sizeof(int), mols->n_angles);

	n_atoms = 0;
	while ( fgets(line, 1024, in1)) {
		sscanf(line, "%i %i %i", &i, &j, &k);
		mols->ia[n_atoms] = i;
		mols->ja[n_atoms] = j;
		mols->ka[n_atoms] = k;

		n_atoms++;
	}
	fclose(in1);
	printf("Angles in.\n");
	n_atoms = 0;
	in1 = fopen("dihedrals.top", "rb");
	if ( in1 == NULL) {
		free(line);
		return -1;
	}
	while ( fgets(line, 1024, in1))
		++n_atoms;
	rewind(in1);

	mols->n_torsionals = n_atoms;
	mols->ik = (int*)calloc( sizeof(int), mols->n_torsionals);
	mols->jk = (int*)calloc( sizeof(int), mols->n_torsionals);
	mols->kk = (int*)calloc( sizeof(int), mols->n_torsionals);
	mols->lk = (int*)calloc( sizeof(int), mols->n_torsionals);

	n_atoms = 0;
	while ( fgets(line, 1024, in1)) {
		sscanf(line, "%i %i %i %i", &i, &j, &k, &l);
		mols->ik[n_atoms] = i;
		mols->jk[n_atoms] = j;
		mols->kk[n_atoms] = k;
		mols->lk[n_atoms] = l;


		n_atoms++;
	}
	fclose(in1);
	printf("Torsionals in.\n");

	n_atoms = 0;

	in1 = fopen("impropers.top", "rb");
	if ( in1 == NULL) {
		free(line);
		return -1;
	}
	while ( fgets(line, 1024, in1))
		++n_atoms;
	rewind(in1);

	mols->n_impropers = n_atoms;
	mols->ip = (int*)calloc( sizeof(int), mols->n_impropers);
	mols->jp = (int*)calloc( sizeof(int), mols->n_impropers);
	mols->kp = (int*)calloc( sizeof(int), mols->n_impropers);
	mols->lp = (int*)calloc( sizeof(int), mols->n_impropers);
	mols->tp = (int*)calloc( sizeof(int), mols->n_impropers);


	n_atoms = 0;
	while ( fgets(line, 1024, in1)) {
		sscanf(line, "%i %i %i %i %i", &i, &j, &k, &l, &m);
		mols->ip[n_atoms] = i;
		mols->jp[n_atoms] = j;
		mols->kp[n_atoms] = k;
		mols->lp[n_atoms] = l;
		mols->tp[n_atoms] = m;
		n_atoms++;
	}
	fclose(in1);
	printf("Impropers in.\n");



	n_atoms = 0;
	in1 = fopen("vdwpairs.top", "rb");
	if ( in1 == NULL) {
		free(line);
		return -1;
	}
	while ( fgets(line, 1024, in1))
		++n_atoms;
	rewind(in1);

	mols->n_pairs = n_atoms;
	mols->vdwpairs_a = (int*)calloc( sizeof(int), mols->n_pairs);
	mols->vdwpairs_b = (int*)calloc( sizeof(int), mols->n_pairs);
	mols->vdw_type = (int*)calloc( sizeof(int), mols->n_pairs);

	n_atoms = 0;
	while ( fgets(line, 1024, in1)) {
		sscanf(line, "%i %i %i", &i, &j, &k);
		mols->vdwpairs_a[n_atoms] = i;
		mols->vdwpairs_b[n_atoms] = j;
		mols->vdw_type[n_atoms] = k;
		n_atoms++;
	}
	fclose(in1);
	printf("Pairs in.\n");

	printf("Import: \n");
	printf("Atoms: %i.\n", mols->n_atoms);
	printf("Bonds: %i.\n", mols->n_bonds);
	printf("Angles: %i.\n", mols->n_angles);
	printf("Tortsionals: %i.\n", mols->n_torsionals);
	printf("VDW Pairs: %i.\n", mols->n_pairs);

	free(line);
	return 0;
}

/**
 *
 *	@brief Save the topology in a MOL2 to top files
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 structure with the topology
 *	@return 0 on success
 *
 */
int dump_topology(MOL2 mol)
{
	FILE *types = NULL;
	FILE *bonds = NULL;
	FILE *angles = NULL;
	FILE *dihedrals = NULL;
	FILE *pairs = NULL;
	FILE *impropers = NULL;
	int i = 0;

	if ( (bonds = fopen("bonds.top", "w")) == NULL)
		printf("Error writing bonds.\n");
	else{
		for (i = 0; i < mol.n_bonds; ++i)
			fprintf(bonds, "%i %i %i\n", mol.bond_a1[i], mol.bond_a2[i], mol.bonds[i]);
		fclose(bonds);
	}

	if ( (types = fopen("types.top", "w")) == NULL)
		printf("Error writing types.\n");
	else{
		for (i = 0; i < mol.n_atoms; ++i)
			fprintf(types, "%i %i\n", i, mol.gaff_types[i]);
		fclose(types);
	}

	if ( (angles = fopen("angles.top", "w")) == NULL)
		printf("Error writing angles.\n");
	else{
		for (i = 0; i < mol.n_angles; ++i)
			fprintf(angles, "%i %i %i\n", mol.ia[i], mol.ja[i], mol.ka[i]);
		fclose(angles);
	}

	if ( (dihedrals = fopen("dihedrals.top", "w")) == NULL)
		printf("Error writing dihedrals.\n");
	else{
		for (i = 0; i < mol.n_torsionals; ++i)
			fprintf(dihedrals, "%i %i %i %i\n", mol.ik[i], mol.jk[i], mol.kk[i], mol.lk[i]);
		fclose(dihedrals);
	}

	if ( (impropers = fopen("impropers.top", "w")) == NULL)
		printf("Error writing impropers.\n");
	else{
		for (i = 0; i < mol.n_impropers; ++i)
			fprintf(dihedrals, "%i %i %i %i %i\n", mol.ip[i], mol.jp[i], mol.kp[i], mol.lp[i], mol.tp[i]);
		fclose(impropers);
	}



	if ( (pairs = fopen("vdwpairs.top", "w")) == NULL)
		printf("Error writing VDW pairs.\n");
	else{
		for (i = 0; i < mol.n_pairs; ++i)
			fprintf(pairs, "%i %i %i\n", mol.vdwpairs_a[i], mol.vdwpairs_b[i], mol.vdw_type[i]);
		fclose(pairs);
	}

	return 0;
}


/**
 *
 *	@brief Apply a rigid tranformation and print the molecule in PDB format in the default output
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 with the molecule
 *	@param simplex 6 variables for the rigid body tranformation (x,y,z,thetha,phy,psi)
 *	@return 0 on success
 *
 */
int dump_simplex_pdb(MOL2 *mol, float *simplex)
{

	int i = 0, j = 0, real_atoms = 0;

	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;

	float ith = 0, iph = 0, ips = 0;
	float theta = 0.0f, phi = 0.0f, psi = 0.0f;
	float sthe = 0.0f, cthe = 0.0f, sphi = 0.0f, cphi = 0.0f, spsi = 0.0f, cpsi = 0.0f;
	float A[3][3], dx = 0.0f, dy = 0.0f, dz = 0.0f;
	const float cdr = 0.0174532777778;

	float comt[3], com[3], Xrot = 0.0f, Yrot = 0.0f, Zrot = 0.0f, VX = 0.0f, VY = 0.0f, VZ = 0.0f;
	int *tor_type = NULL;
	double energy = 0.0;

	for ( i = 0; i < mol->n_atoms; ++i)
		if ( mol->backbone[i] != 1)
			++real_atoms;

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

	for (i = 0; i < mol->n_atoms; ++i) {

		if ( mol->backbone[i] != 1)
			printf("ATOM %6i  %2s  UNK     0    %8.3f%8.3f%8.3f\n", i + 1, elems[ mol->atoms[i] - 1], mol->x[i], mol->y[i], mol->z[i]);

	}

	printf("TER\n");



}

/**
 *
 *	@brief Dump the whole molecule in PDB format to the default output
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 with the molecule
 *	@return 0 on success
 *
 */
int dump_all_pdb(MOL2 *mol)
{
	int i;

	for (i = 0; i < mol->n_atoms; ++i)
		printf("ATOM %6i  %2s  UNK     0    %8.3f%8.3f%8.3f\n", i + 1, elems[ mol->atoms[i] - 1], mol->x[i], mol->y[i], mol->z[i]);
	printf("TER\n");

	return 0;
}

/**
 *
 *	@brief Dump only the ligand atoms of a MOL2 in PDB format to the default output
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 with the molecule
 *	@return 0 on success
 *
 */
int dump_pdb(MOL2 *mol)
{
	int i;

	for (i = 0; i < mol->n_atoms; ++i) {

		if ( mol->backbone[i] != 1)
			printf("ATOM %6i  %2s  UNK     0    %8.3f%8.3f%8.3f\n", i + 1, elems[ mol->atoms[i] - 1], mol->x[i], mol->y[i], mol->z[i]);

	}

	printf("TER\n");

	return 0;
}

/**
 
    @brief Dump only the ligand atoms of a MOL2 in PDB format to the default output
    @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
    @param mol MOL2 with the molecule
    @param out FILE handler output file
    @return 0 on success
 
 */

int dump_pdb_to_file(MOL2 *mol, FILE *out)
{
        int i;

        for (i = 0; i < mol->n_atoms; ++i) {

                if ( mol->backbone[i] != 1)
                        fprintf(out,"ATOM %6i  %2s  UNK     0    %8.3f%8.3f%8.3f\n", i + 1, elems[ mol->atoms[i] - 1], mol->x[i], mol->y[i], mol->z[i]);

        }

        fprintf(out,"TER\n");

        return 0;
}



/**
 * 
 *     @brief Dump only the ligand atoms of a MOL2 in PDB format to the default output with
 *            charges and radii and keeping the original names
 *     @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *     @param mol MOL2 with the molecule
 *     @return 0 on success
 *
 */
int dump_pdb_conservative(MOL2 *mol)
{
        int i;

        for (i = 0; i < mol->n_atoms; ++i) {

/*                if ( mol->backbone[i] != 1)*/
                        printf("ATOM %6i%3s  %3s     0    %8.3f%8.3f%8.3f  %1.2f%7.4f\n", i + 1, mol->atom_names[i], mol->res_names[i], mol->x[i], mol->y[i], mol->z[i], mol->radius[i],mol->pcharges[i]);

        }

        printf("TER\n");

        return 0;
}

/**
 *
 *  @brief Dump only the ligand atoms of a MOL2 in PDB format to the default output with
 *        charges and radii and keeping the original names
 *  @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *  @param mol MOL2 with the molecule
 *  @param out FILE handler for output file
 *  @return 0 on success
 *
 **/
int dump_pdb_conservative_to_file(MOL2 *mol, FILE *out)
{
        int i;

        for (i = 0; i < mol->n_atoms; ++i) {

/*                if ( mol->backbone[i] != 1)*/
                        fprintf(out,"ATOM %6i%3s  %3s     0    %8.3f%8.3f%8.3f  %1.2f%7.4f\n", i + 1, mol->atom_names[i], mol->res_names[i], mol->x[i], mol->y[i], mol->z[i], mol->radius[i],mol->pcharges[i]);

        }

        fprintf(out,"TER\n");

        return 0;
}



/**
 *
 *      @brief Determine if an atom has an hydrogen attached
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @param mol MOL2 with the molecule
 *      @param atom atom number
 *      @return 1 if it has or 0 if has not
 *
 */

int has_hydrogens( MOL2 *mol, int atom )
{
	int res = 0, j = 0;
        for (j = 0; j < mol->n_bonds; ++j) 
	{
          	if ( mol->bond_a1[j] == (atom + 1)){
			if( mol->atoms[mol->bond_a2[j]-1] == 4)
				res = 1;
		}else if( mol->bond_a2[j] == (atom + 1)){
                        if( mol->atoms[mol->bond_a1[j]-1] == 4)
                                res = 1;
		}
	}	
	return res;
}


/**
 *
 *	@brief Determine if two atoms are 1-4 pairs (i-X-X-j)
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 with the molecule
 *	@param atom First atom
 *	@param second_atom Second atom
 *	@return 1 if they are 1-4 pair or 0 if not
 *
 */
int q14_bond( MOL2 *mol, int atom, int second_atom)
{
	int j, k;
	int flag = 0;
	int vecinos[4];

	for (j = 0; j < 4; ++j)
		vecinos[j] = 0;
	k = 0;
	for (j = 0; j < mol->n_bonds; ++j) {
		if ( mol->bond_a1[j] == (atom + 1)) {
			vecinos[k] = mol->bond_a2[j];
			k++;
		}else if ( mol->bond_a2[j] == (atom + 1)) {
			vecinos[k] = mol->bond_a1[j];
			k++;
		}
	}

	for (j = 0; j < k; ++j)
		if ( q13_bond(mol, vecinos[j] - 1, second_atom) > 0)
			flag = 1;

	return flag;


}

/**
 *
 *	@brief Determine if two atoms are bonded
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 with the molecule
 *	@param atom First atom
 *	@param second_atom Second atom
 *	@return 1 they are bonded, 0 if not
 *
 */
int bonded( MOL2 mol, int atom, int second_atom)
{
	int j;

	for (j = 0; j < mol.n_bonds; ++j)
		if ( (mol.bond_a1[j] == (atom + 1) && mol.bond_a2[j] == (second_atom + 1)) || (mol.bond_a2[j] == (atom + 1) && mol.bond_a1[j] == (second_atom + 1)))
			return 1;
	return 0;

}

/**
 *
 *	@brief Determine if two atoms are 1-3 pair
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 structure with the molecule
 *	@param atom First atom
 *	@param second_atom Second atom
 *	@return 1 if they are a 1-3 pair, 0 if not
 *
 */
int q13_bond( MOL2 *mol, int atom, int second_atom)
{
	int j, k;
	int vecinos[4];
	int flag = 0;

	for (j = 0; j < 4; ++j)
		vecinos[j] = 0;
	k = 0;
	for (j = 0; j < mol->n_bonds; ++j) {
		if ( mol->bond_a1[j] == (atom + 1)) {
			vecinos[k] = mol->bond_a2[j];
			k++;
		}else if ( mol->bond_a2[j] == (atom + 1)) {
			vecinos[k] = mol->bond_a1[j];
			k++;
		}
	}

	for (j = 0; j < k; ++j)
		if ( get_number_any_bond_p(mol, vecinos[j], second_atom + 1) > 0)
			flag = 1;

	return flag;

}

/**
 *
 *	@brief Get number of bonds of an atom
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 structure
 *	@param atom Atom number
 *	@param bond_type Bond type (1 = single, 2 = double, etc .. or 0 for any)
 *	@return number of bonds
 *
 */
int get_bonds( MOL2 mol, int atom, int bond_type)
{
	int res = 0;
	int i = 0;

	for ( i = 0; i < mol.n_bonds; ++i)
		if ( (mol.bond_a1[i] == atom || mol.bond_a2[i] == atom ) && (bond_type == 0 || mol.bonds[i] == bond_type))
			++res;
	return res;
}


/**
 *
 *	@brief Get the number of an atom
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol Pointer to MOL2 structure
 *	@param atom Atom number
 *	@param bond_type Bond type (1 = single, 2 = double, etc .. or 0 for any)
 *	@return number of bonds
 *
 */
int get_bonds_p( MOL2 *mol, int atom, int bond_type)
{
	int res = 0;
	int i = 0;

	for ( i = 0; i < mol->n_bonds; ++i)
		if ( (mol->bond_a1[i] == atom || mol->bond_a2[i] == atom ) && (bond_type == 0 || mol->bonds[i] == bond_type))
			++res;
	return res;

}

/**
 *
 *	@brief Check if two atoms are bonded
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol Pointer to MOL2 with the molecule
 *	@param atom First atom
 *	@param second_atom Second atom (0 is any)
 *	@return number of bonds or 0 if they are not bonded
 *
 */
int get_number_any_bond_p( MOL2 *mol, int atom, int second_atom) /* Second atom = 0 -> Any */
{
	int res = 0;
	int i = 0;

	for ( i = 0; i < mol->n_bonds; ++i)
		if ( (mol->bond_a1[i] == atom &&  (mol->bond_a2[i] == second_atom || second_atom == 0) ) || (mol->bond_a2[i] == atom && (mol->bond_a1[i] == second_atom || second_atom == 0)) )
			++res;
	return res;
}

/**
 *
 *      @brief Check if two atoms are bonded
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @param mol MOL2 with the molecule
 *      @param atom First atom
 *      @param second_atom Second atom (0 is any)
 *      @return number of bonds or 0 if they are not bonded
 *
 */
int get_number_any_bond( MOL2 mol, int atom, int second_atom) /* Second atom = 0 -> Any */
{
	int res = 0;
	int i = 0;

	for ( i = 0; i < mol.n_bonds; ++i)
		if ( (mol.bond_a1[i] == atom &&  (mol.bond_a2[i] == second_atom || second_atom == 0) ) || (mol.bond_a2[i] == atom && (mol.bond_a1[i] == second_atom || second_atom == 0)) )
			++res;
	return res;
}

/**
 *
 *	@brief Get number of bonds of an atom using bond and atom types
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 with the molecule
 *	@param atom Atom
 *	@param bond_type Type of bond
 *	@param atom_type Type of atom bonded
 *	@return number of bonds with the criteria
 *
 */
int get_number_bond_by_atom(MOL2 mol, int atom, int bond_type, int atom_type)
{
	int res = 0;
	int i = 0;

	for ( i = 0; i < mol.n_bonds; ++i) {
		if ( mol.bonds[i] == bond_type || bond_type == 0) {
			if ( mol.bond_a1[i] == atom) {
				if ( mol.atoms[(mol.bond_a2[i] - 1)] == atom_type)
					++res;
			}else if ( mol.bond_a2[i] == atom)
				if ( mol.atoms[(mol.bond_a1[i] - 1)] == atom_type)
					++res;
		}
	}
	return res;
}

/**
 *
 *	@brief Check if two atoms are bonded with an specific type of bond
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 structure with the molecule
 *	@param atom First atom
 *	@param type Bond type
 *	@param second_atom Second atom or 0 if any
 *	@return number of bonds
 *
 */
int get_number_type_bonds( MOL2 mol, int atom, int type, int second_atom) /* Second atom = 0 -> Any */
{
	int res = 0;
	int i = 0;

	for ( i = 0; i < mol.n_bonds; ++i)
		if ( (mol.bond_a1[i] == atom && ( mol.atoms[(mol.bond_a2[i] - 1)] == second_atom || second_atom == 0) && ( mol.bonds[i] == type || type == 0 )) || (mol.bond_a2[i] == atom && ( mol.atoms[(mol.bond_a1[i] - 1)] == second_atom || second_atom == 0) && ( mol.bonds[i] == type || type == 0) ))
			++res;
	return res;
}

/**
 *
 *	@brief Return the index of the last element of the ring graph
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param path vector with the walking path
 *	@return index of the last element
 *
 */
int last_link(int *path)
{
	int i = 0;

	for ( i = 0; i < MAX_RING_MEMBERS; ++i)
		if ( path[i] < 0)
			return i;

	return 0;
}


/* Only for debugging */
/**
 *
 *	@brief Show a walking path for ring graph detector. (For debuggin propouses)
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param path Path
 *	@return 0 on success
 *
 */
int show_path(int *path)
{
	int i = 0;

	for (i = 0; i < MAX_RING_MEMBERS; ++i)
		printf("%i,", path[i]);
	printf("\n");
	return 0;
}


/**
 *
 *
 *	@brief Recursive function to find rings in a graph and dihedral information
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 with the molecule
 *	@param path Current path
 *	@param n_atom Current atom to explore
 *	@param last_n_atom Previous atom in the path
 *	@param num_rings Number of rings found
 *	@param myik Dihedral information i
 *	@param myjk Dihedral information j
 *	@param mykk Dihedral information k
 *	@param mylk Dihedral information l
 *	@param elem Current index of elements
 *
 */
int find_rings(MOL2 mol, int *path, int n_atom, int last_n_atom, int num_rings, int **myik, int **myjk, int **mykk, int **mylk,
	       int *elem)
{
	int *new_path = NULL;
	int children = 0;
	int atom_number = 0;
	int i = 0;
	int j = 0;
	int k = 0;
	int i2 = 0;
	int j2 = 0;
	int k2 = 0;
	int tmp_num_rings = 0;
	float com_x = 0.0f;
	float com_y = 0.0f;
	float com_z = 0.0f;
	int *ik = NULL;
	int *jk = NULL;
	int *kk = NULL;
	int *lk = NULL;
	int myelem;

	int flag = 0;
	int flag2 = 0;
	int dc = 0;

	tmp_num_rings = num_rings;
	ik = *myik;
	jk = *myjk;
	kk = *mykk;
	lk = *mylk;

	/* Append new member */
	k = last_link(path);
	path[k] = n_atom;
	myelem = *elem;
	myelem++;

	if ( tmp_num_rings > MAX_RINGS)
		return 0;
	for ( i2 = 0; i2 < (MAX_RING_MEMBERS ); ++i2) {
		for ( j2 = 0; j2 < MAX_RING_MEMBERS; ++j2) {
			if ( path[i2] == path[j2] && path[i2] != -1 && i2 != j2) { /* Ring closure */

				tmp_num_rings = tmp_num_rings + 1;
				return 1;
			}
		}

	}

	if ( (children = get_number_type_bonds(mol, n_atom, 0, 0)) >= 0 ) { /* Prune */
		for (i = 0; i < mol.n_bonds; ++i) {
			if ( mol.bond_a1[i] == n_atom || mol.bond_a2[i] == n_atom ) {
				if ( mol.bond_a1[i] != last_n_atom && mol.bond_a2[i] != last_n_atom) {
					if ( (new_path = calloc(MAX_RING_MEMBERS, sizeof(int))) == NULL) {
						fprintf(stderr, "Error. Cant allocate memory.\n");
						fflush(stderr);
						return -2;
					}
					for (k = 0; k < MAX_RING_MEMBERS; ++k)
						new_path[k] = path[k];

					if ( mol.bond_a1[i] == n_atom)
						atom_number = mol.bond_a2[i];
					else
						atom_number = mol.bond_a1[i];

					tmp_num_rings++;
					if ( *elem >= 3) {
						flag = flag2 = -1;
						for ( dc = 0; dc < 1000; ++dc) {
							if (ik[dc] == -1) {
								flag2 = dc;
								break;
							}
							if ( jk[dc] == (new_path[myelem - 2] - 1) && kk[dc] == (new_path[myelem - 1] - 1))
								flag = 1;
						}
						if ( flag != 1) {
							ik[flag2] = new_path[myelem - 3] - 1;
							jk[flag2] = new_path[myelem - 2] - 1;
							kk[flag2] = new_path[myelem - 1] - 1;
							lk[flag2] = new_path[myelem] - 1;

						}

					}
					if ( find_rings(mol, new_path, atom_number, n_atom, tmp_num_rings, &ik, &jk, &kk, &lk, &myelem) < 0) {
						fprintf(stderr, "Error. Rescursive fail.\n");
						fflush(stderr);
						return -5;
					}
					atom_number = 0;
					free(new_path);

				}else{}

			}
		}
	}else{
		printf("Atomo %i tiene %i enlaces.\n", n_atom, children);
		return 0;
	}

	return 0;
}

/**
 *
 *	@brief Entry point for the recursive function to find rings in a graph and dihedral information
 *	@author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *	@param mol MOL2 with the molecule
 *	@param myik Dihedral information i
 *	@param myjk Dihedral information j
 *	@param mykk Dihedral information k
 *	@param mylk Dihedral information l
 *	@return number of rings found
 *
 */
int get_number_of_rings(MOL2 mol, int **myik, int **myjk, int **mykk, int **mylk)
{
	int *new_path = NULL;
	int i = 0;
	int *ik, *jk, *kk, *lk;
	int elem = 0;

	ik = *myik;
	jk = *myjk;
	kk = *mykk;
	lk = *mylk;

	for ( i = 0; i < 1000; i++) {
		ik[i] = -1;
		jk[i] = -1;
		kk[i] = -1;
		lk[i] = -1;
	}
	if ( (new_path = calloc(MAX_RING_MEMBERS, sizeof(int))) == NULL) {
		fprintf(stderr, "Error. Cant allocate memory.\n");
		fflush(stderr);
		return -2;
	}

	for (i = 0; i < MAX_RING_MEMBERS; ++i)
		new_path[i] = -1;
	new_path[0] = 0;
	if ( find_rings(mol, new_path, 1, 0, 0, &ik, &jk, &kk, &lk, &elem) < 0) {
		fprintf(stderr, "Error. Recursive fail.\n");
		fflush(stderr);
		return -3;
	}
	free(new_path);

	for (i = 0; i < 1000; ++i)
		if ( ik[i] == -1)
			break;

	return i;
}

int is_planar(MOL2 *mol, int i)
{
        int j = 0, k = 0;
        int vecinos[8], vecinos2[8];
        float vecs[6][3];
        float res = 0, mod1 = 0.0f, mod2 =0.0f;

        for (j = 0; j < 8; ++j) {
              vecinos[j] = 0;
              vecinos2[j] = 0;
        }
        k = 0;
        for (j = 0; j < mol->n_bonds; ++j)
        {
                if ( mol->bond_a1[j] == (i + 1)) {
                        vecinos[k] = mol->bond_a2[j] - 1;
                        vecinos2[k] = mol->bonds[j];
                        k++;
                }else if ( mol->bond_a2[j] == (i + 1)) {
                        vecinos[k] = mol->bond_a1[j] - 1;
                        vecinos2[k] = mol->bonds[j];
                        k++;
                }
        }


        if( k <= 2)
          return 1;

        if( k >= 3)
        {
                vecs[0][0] = mol->x[i] - mol->x[vecinos[0]];
                vecs[0][1] = mol->y[i] - mol->y[vecinos[0]];
                vecs[0][2] = mol->z[i] - mol->z[vecinos[0]];

                vecs[1][0] = mol->x[i] - mol->x[vecinos[1]];
                vecs[1][1] = mol->y[i] - mol->y[vecinos[1]];
                vecs[1][2] = mol->z[i] - mol->z[vecinos[1]];

                vecs[2][0] = mol->x[i] - mol->x[vecinos[2]];
                vecs[2][1] = mol->y[i] - mol->y[vecinos[2]];
                vecs[2][2] = mol->z[i] - mol->z[vecinos[2]];

                vecs[3][0] = vecs[0][1]*vecs[1][2] - vecs[1][1]*vecs[0][2];
                vecs[3][1] = vecs[0][2]*vecs[1][0] - vecs[0][0]*vecs[1][2];
                vecs[3][2] = vecs[0][0]*vecs[1][1] - vecs[0][1]*vecs[1][0];

                res = vecs[2][0]*vecs[3][0] + vecs[2][1]*vecs[3][1] + vecs[2][2]*vecs[3][2];

                mod1 = sqrt(vecs[2][0]*vecs[2][0])+(vecs[2][1]*vecs[2][1])+(vecs[2][2]*vecs[2][2]);
                mod2 = sqrt(vecs[3][0]*vecs[3][0])+(vecs[3][1]*vecs[3][1])+(vecs[3][2]*vecs[3][2]);

                res = res / (mod1*mod2);

/*                printf("\n%f\n\n",res);*/

                if( fabs(res) < 0.1)
                 return 0;
	}
	return -1;
}



/*
 *	Return:
 *
 *	5 if nitro
 *	0 carboxylate
 *	1 sulfate S(O)(O)(O)(O)
 *	2 phosphate P(O)(O)(O)(O)
 *
 */
int gaus_is_carboxylate(MOL2 *mol, int l)
{
        int i = 0, nh = -1, no = 0;
        int val = 4;
        if( mol->gaff_types[l] != O )
         return val;

        for(i = 0; i < mol->n_bonds; i++)
        {
                if( mol->bond_a1[i] == (l+1))
                   nh =  mol->bond_a2[i]-1;
                else if( mol->bond_a2[i] == (l+1))
                   nh =  mol->bond_a1[i]-1;
        }

        if ( nh != -1)
        {
                if( mol->atoms[nh] == 3) /* Nitro */
                {
                  val = 5;
                }else{

	                for(i = 0; i < mol->n_bonds; i++)
	                {
	                        if( mol->bond_a1[i] == (nh+1))
	                        {
	                           if( mol->gaff_types[ mol->bond_a2[i]-1] == O)
	                                ++no;
	                        }else if( mol->bond_a2[i] == (nh+1)){
	                           if( mol->gaff_types[ mol->bond_a1[i]-1] == O)
	                                ++no;
	       	                }
       		        }
		}


        if( mol->atoms[nh] == 1 && no == 2)
          val = 0;
	else if( mol->atoms[nh] == 6 && no == 3)
          val = 1;
	else if (mol->atoms[nh] == 5 && no == 3)
	  val = 2;


	}
        return val;
}

int dump_mol2_to_file(MOL2 *mol, FILE *out)
{
        int i;
	char *elems[10] = { "C ", "O ", "N ", "H ", "P ", "S ", "I ", "Br", "Cl", "F " };
	fprintf(out,"@<TRIPOS>MOLECULE\n");
	fprintf(out,"%s\n",mol->comment);
	fprintf(out,"%4i %4i %4i %4i %4i\n",mol->n_atoms, mol->n_bonds, 0, 0,0);
	fprintf(out,"\n\n\n\n");
	fprintf(out,"@<TRIPOS>ATOM\n");

        for (i = 0; i < mol->n_atoms; ++i) {
		fprintf(out,"%4i %10s %8.4f %8.4f %8.4f %-4s\n", i+1, mol->atom_names[i],  mol->x[i], mol->y[i], mol->z[i],mol->atom_types_sybyl[i]);
        }
	fprintf(out,"@<TRIPOS>BOND\n");
        for (i = 0; i < mol->n_bonds; ++i) {
                fprintf(out,"%4i %4i %4i %i\n",i+1, mol->bond_a1[i], mol->bond_a2[i], mol->bonds[i]);
        }
	fprintf(out,"\n");

        return 0;
}

