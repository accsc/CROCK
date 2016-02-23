/**
 *
 *      @file util2.c
 *      @brief Root mean squared deviation calculation subroutines
 *
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
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

/**
 *
 *	@brief Calculate RMSD (root mean squared deviation)  with reordered atoms
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *	@param mols MOL2 with the molecule
 *	@param natom_ref Number of atoms of the reference molecule
 *	@param order Order of atoms for the candidate molecule
 *	@param oder1 Order of atoms for the reference molecule
 *	@param x Reference molecule coordinates X
 *	@param y Reference molecule coordinates Y
 *	@param z Reference molecule coordinates Z
 *	@return RMSD value
 *
 */
float dump_rmsd(MOL2 mols, int natom_ref, int *order, int *order1, float *x, float *y, float *z)
{

	int compared = 0;
	int i = 0;
	int j = 0;
	float myrmsd = 0;

	for (i = 0; i < mols.n_atoms; ++i) {
		for (j = 0; j < natom_ref; ++j) {
			if ( order1[j] == order[i] && order[i] != 0) {
				printf("%i - %i.\n", order1[j], order[i]);
				compared = compared + 1;
				myrmsd = myrmsd + (((x[j] - mols.x[i]) * (x[j] - mols.x[i])) + ((y[j] - mols.y[i]) * (y[j] - mols.y[i])) + ((z[j] - mols.z[i]) * (z[j] - mols.z[i])));
			}
		}
	}

	return myrmsd;

}


/**
 *
 *	@brief Get the new index for an atom in an ordered molecule
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *	@param mol MOL2 with the molecule
 *	@param index Vector with new indices
 *	@param query Current index for an atom
 *	@return new index
 *
 */
int get_new_index(MOL2 mol, int *index, int query)
{
	int i = 0;

	for ( i = 0; i < mol.n_atoms; ++i)
		if ( index[i] == query)
			return i;

}

/**
 *
 *	@brief Write to the default output the molecule in MOL2 format
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *	@param mol MOl2 with the molecule
 *	@param index Vector with indeces of the ordered molecule
 *	@return null
 *
 */
void dump_mol2(MOL2 mol, int *index)
{


	char *retro[33] = { "C.3", "C.2", "C.1", "C.ar", "C.cat", "Du.C", "H.spc", "H.t3p", "H",   "O.3",  "O.2",  "O.co2", "O.spc", "O.t3p", "S.3", "S.2", "S.O2", "S.O", "P.3", "F",
			    "Cl",  "Br",  "I  ", "LP",	 "Du",	  "Any",  "N.3",   "N.2  ", "N.1", "N.ar", "N.am", "N.pl3", "N.4" };

	char *retro2[5] = { " 1", " 2", " 3", "", "ar" };
	int i = 0;

	printf("@<TRIPOS>MOLECULE\n\n");
	printf("%5i %5i     0     0     0\n", mol.n_atoms, mol.n_bonds);
	printf("SMALL\n");
	printf("USER_CHARGES\n\n");
	printf("@<TRIPOS>ATOM\n");
	for ( i = 0; i < mol.n_atoms; ++i)
		printf("%7i  %c%i        %8.4f  %8.4f  %8.4f %5s    1  LIG         0.0000\n", i + 1, retro[mol.atoms[index[i]] - 1][0], i + 1, mol.x[index[i]], mol.y[index[i]], mol.z[index[i]], retro[mol.atoms[index[i]] - 1]);

	printf("@<TRIPOS>BOND\n");
	for ( i = 0; i < mol.n_bonds; ++i)
		printf("%6i%6i%6i%5s\n", i + 1, get_new_index(mol, index, mol.bond_a1[i]),  get_new_index(mol, index, mol.bond_a2[i]), retro2[mol.bonds[i] - 1]);


}



/* Oh shame on me! */
/**
 *
 *	@brief Translate MOL2 atom types to numbers
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *	@param atom_type Name of the MOL2 type
 *	@return index for that type
 *
 */
int translate_atom(char *atom_type)
{

	if ( strlen(atom_type) <= 0)
		return -1;


	if ( strstr(atom_type, "C.3") != NULL)
		return 1;
	else if ( strstr(atom_type, "C.2") != NULL)
		return 2;
	else if ( strstr(atom_type, "C.1") != NULL)
		return 3;
	else if ( strstr(atom_type, "C.ar") != NULL)
		return 4;
	else if ( strstr(atom_type, "C.cat") != NULL)
		return 5;
	else if ( strstr(atom_type, "Du.C") != NULL)
		return 6;
	else if ( strstr(atom_type, "H.spc") != NULL)
		return 7;
	else if ( strstr(atom_type, "H.t3p") != NULL)
		return 8;
	else if ( strstr(atom_type, "H") != NULL)
		return 9;
	else if ( strstr(atom_type, "O.3") != NULL)
		return 10;
	else if ( strstr(atom_type, "O.2") != NULL)
		return 11;
	else if ( strstr(atom_type, "O.co2") != NULL)
		return 12;
	else if ( strstr(atom_type, "O.spc") != NULL)
		return 13;
	else if ( strstr(atom_type, "O.t3p") != NULL)
		return 14;
	else if ( strstr(atom_type, "S.3") != NULL)
		return 15;
	else if ( strstr(atom_type, "S.2") != NULL)
		return 16;
	else if ( strstr(atom_type, "S.O2") != NULL)
		return 17;
	else if ( strstr(atom_type, "S.O") != NULL)
		return 18;
	else if ( strstr(atom_type, "P.3") != NULL)
		return 19;
	else if ( strstr(atom_type, "F") != NULL)
		return 20;
	else if ( strstr(atom_type, "Cl") != NULL)
		return 21;
	else if ( strstr(atom_type, "Br") != NULL)
		return 22;
	else if ( strstr(atom_type, "I") != NULL)
		return 23;
	else if ( strstr(atom_type, "LP") != NULL)
		return 24;
	else if ( strstr(atom_type, "Du") != NULL)
		return 25;
	else if ( strstr(atom_type, "Any") != NULL)
		return 26;
	else if ( strstr(atom_type, "N.3") != NULL)
		return 27;
	else if ( strstr(atom_type, "N.2") != NULL)
		return 28;
	else if ( strstr(atom_type, "N.1") != NULL)
		return 29;
	else if ( strstr(atom_type, "N.ar") != NULL)
		return 30;
	else if ( strstr(atom_type, "N.am") != NULL)
		return 31;
	else if ( strstr(atom_type, "N.pl3") != NULL)
		return 32;
	else if ( strstr(atom_type, "N.4") != NULL)
		return 33;
	else
		return 0;
}

