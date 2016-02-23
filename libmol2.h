/**
 *
 *      @file libmol2.h
 *      @brief Molecular library definitions and structures
 *
 *      @author Alvaro Cortes Cabrera <alvaro.cortes@uah.es>
 *      @date 01/10/2011
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


#define SINGLE_BOND 1
#define DOUBLE_BOND 2
#define TRIPLE_BOND 3
#define AROMATIC 5


/**
 * Structure to store residues
 */
typedef struct{
/** Number of atoms of the residue*/
int n_atoms;  
/** Atoms of the resiude */
int *atoms;   
/** Residue type */
int type;     
/** Index of the C alpha atom */
int ca; 
/** Index of the C beta atom */
int cb;     
/** Vector of translation */
float transvector[3]; 
/** Rotation matrix to superimpose the residue */
float rotM[3][3];     
} RESIDUE;


/**
 * Structure for a rotamer
 */
typedef struct{
/** Number of atoms for the rotamer */
int n_atoms;
/** Number of rotamer of the library */
int n_rotamers;
/** X coordinates for the rotamer */
float *x;
/** Y coordinates for the rotamer */
float *y;
/** Z coordinates for the rotamer */
float *z;
} ROTAMER;


/**
 * Structure for a conformer
 */
typedef struct{
/** Conformer number */
int conf_number;
/** Coordinates for the conformer. X */
float *x;
/** Coordinates for the conformer. Y */
float *y;
/** Coordinates for the conformer. Z */
float *z;
/** Charges for the conformer */
float *pcharges;
/** Radii for the conformer */
float *radius;
} CONFORMER;


/**
 *	Structure to hold molecular information (ligand, template, protein, etc.)
 *
 */
typedef struct{

/** Number of atoms */
int n_atoms;
/** Number of bonds */
int n_bonds;
/** Comments. Not used */
char *comment;  
/** Global charge. Not used */
int charges;
/** Small flag. Not used */
int small_flag;
/** Bond types vector */
int *bonds;  
/** Bonded atoms vector. First atom */
int *bond_a1;     
/** Bonded atoms vector. Second atom */
int *bond_a2;
/** Bonded atoms vector. Distances */
float *bond_dist;
/** Atoms vector */
int *atoms;   
/** In ring flag for atoms */
int *ringer;
/** Aromatic flag for atoms */
int *aromatic;
/** Angles in the molecule. Atom i*/
int *ia;
/** Angles in the molecule. Atom j*/
int *ja;
/** Angles in the molecule. Atom k*/
int *ka;
/** Number of angles in the molecule */
int n_angles;

/** Torsionals in the molecule. Atom i*/
int *ik;
/** Torsionals in the molecule. Atom j*/
int *jk;
/** Torsionals in the molecule. Atom k*/
int *kk;
/** Torsionals in the molecule. Atom l*/
int *lk;
/** Number of torsional in the molecule */
int n_torsionals;

/** Impropers in the molecule. Atom i*/
int *ip;
/** Impropers in the molecule. Atom j*/
int *jp;
/** Impropers in the molecule. Atom k*/
int *kp;
/** Impropers in the molecule. Atom l*/
int *lp;
/** Impropers in the molecule. Type */
int *tp;
/** Number of impropers */
int n_impropers;

/** van der Waals pairs. Atom a */
int *vdwpairs_a;
/** van der Waals pairs. Atom b */
int *vdwpairs_b;
/** van der Waals pairs. Type of interaction */
int *vdw_type;
/** Number of internal non-bonded interactions */
int n_pairs;

/** Atom types for GAFF */
int *gaff_types;

/** Implicit solvation model atom types */
int *ism_types;

/** Implicit solvation model atom selection */
int *ism_selection;

/** Atomic partial charges for the current conformer */
float *pcharges;
/** Atomic radii for the current conformer */
float *radius;

/** Atomic coordinates for the current conformer. X */
float *x;
/** Atomic coordinates for the current conformer. Y */
float *y;
/** Atomic coordinates for the current conformer. Z */
float *z;

/** Atomic forces for the current conformer. dx */
float *grads_X;
/** Atomic forces for the current conformer. dy */
float *grads_Y;
/** Atomic forces for the current conformer. dz */
float *grads_Z;

/** Fragment belonging flag */
int *fragment_flag;
/** Number of fragments in the molecule */
int n_fragments;
/** Protein flag. Temporal universal flag for the molecule */
int *backbone;
/** Temporal flag */
int *selection;
/** Protein flag  */
int *protein;

/** Reside number */
int *res_num;
/** Residue type */
int *res_type;
/** Residue names */
char **res_names;
/** Atom names */
char **atom_names;

/** Conformers for the current molecule */
CONFORMER *conformers;
/** Number of conformers for the current molecule */
int nconformers;
/** Index of the current conformer in the default structure */
int default_conformer;

} MOL2;


/**
 *	Groups for mistmachts structure
 */

typedef struct{
/** Reference to the MOL2 of these groups  */
MOL2 *mol;
/** Number of groups in the structure */
int n_groups;
/** Atom indexes for groups */
int **atoms;
/** Number of atoms in groups */
int *n_atoms;
/** X coordiantes for the groups */
float *x;
/** Y coordinates for the groups */
float *y;
/** Z coordinates for the groups */
float *z;
/** Group types */
int *types;
} GROUP;

/**
 *	Possible Pharmacophoric Points (PPP) structure
 */
typedef struct{
/** Index of the point */
int index;
/** X coordinate of the point */
float x;
/** Y coordinate of the point */
float y;
/** Z coordinate of the point */
float z;
/** Next element in the single linked list */
struct PPP *next;
} PPP;


/**
 *	Structure for ring links 
 */
typedef struct _mol_link {
	/** Index of the atom */
        int n_atom;
	/** Next element in the single list */
        struct _mol_link *next;
} MOL_LINK;


/**
 *	RING structure
 */
typedef struct{
	/** First element of the ring */
        MOL_LINK *first_link;
	/** Number of members of the ring */
        int members;
 	/** Aromatic flag */
        int saturated;
        /** Next ring in the list */
        struct RING *next;
} RING;


/** Proteic residues name */
char *residues[20] = {"THR","VAL","SER","HIS","LYS","ARG","MET","CYS","GLU","GLN","ASP","ASN","PHE","TYR","TRP","PRO","GLY","ALA","LEU","ILE"};
/** Number of rotamer for each residue in Richardson library */
int Nrotamers[20] =  {  415,   38,  228,  545,  878,  812,  350,  222,  292,  496,  292,  543,  410, 1230,    0,   47,    0,    0,  256,  293};
/** Number of atoms for each residue */
int aasize[ 20 ] =   {   14,   16,   11,   18,   22,   24,   17,   11,   15,   17,   12,   14,   20,   21,   24,   14,    0,    0,   19,   19};

#define THR 1
#define VAL 2
#define SER 3
#define HIS 4
#define LYS 5
#define ARG 6
#define MET 7
#define CYS 8
#define GLU 9
#define GLN 10
#define ASP 11
#define ASN 12
#define PHE 13
#define TYR 14
#define TRP 15
#define PRO 16
#define GLY 17
#define ALA 18
#define LEU 19
#define ILE 20
