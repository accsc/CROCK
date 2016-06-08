/**
 *
 *      @file rings.c
 *      @brief Function to detect rings in a molecule and aromaticity
 *
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *      @date 01/10/2010
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/**
 *	These routines were improved on 3/06/2016
 *	Now we dont do random walks anymore. 
 *	We remember which vertex have been visited
 */

#define MAX_RING_MEMBERS 800    /* This should be number of atoms of the molecule! */

/**
 *
 *	@brief Breadth-first search routine for finding rings in a graph
 *	@author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *	@param mol MOL2 structure with the molecule
 *	@param path Current walking path
 *	@param n_atom Current atom to explore
 *	@param last_n_atom Last atom explored
 *	@param num_rings Number of rings found so far
 *	@param ringer Ring flag for atoms
 *	@param aromatic Aromatic flag for atoms
 *	@return negative number on error. 0 on reach limit and 1 on success
 *
 */
int find_rings2(MOL2 mol, int *path, int n_atom, int last_n_atom, int *nring, int **ringer, int **aromatic, int **visited, int ***rings)
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
	int ll = 0;
	int tmp_num_rings = 0;
	float com_x = 0.0f;
	float com_y = 0.0f;
	float com_z = 0.0f;
	int *myringer = NULL;
	int *myaromatic = NULL;
	int *myvisited = NULL;
	int **myrings = NULL;
	int elec = 0;
	int aroflag = 0;
	int flag = 0;
	int current_ring = 0;

	myringer = *ringer;
	myrings = *rings;
	myaromatic = *aromatic;
	myvisited = *visited;
	/* Append new member */
	k = last_link(path);
	path[k] = n_atom;
	myvisited[n_atom-1] = 1;

	for (i = 0; i < mol.n_bonds; ++i) 
	{
		if ( mol.bond_a1[i] == n_atom || mol.bond_a2[i] == n_atom ) 
		{
			if ( mol.bond_a1[i] != last_n_atom && mol.bond_a2[i] != last_n_atom) 
			{
                                if ( mol.bond_a1[i] == n_atom)
                                        atom_number = mol.bond_a2[i];
                                else
                                        atom_number = mol.bond_a1[i];

				if( myvisited[atom_number-1] == 0 )
				{
				  if( get_number_type_bonds(mol, atom_number, 0, 0) > 1)
                                  {
				   if((new_path=calloc(mol.n_atoms+2,sizeof(int)))==NULL) 
				   {
					fprintf(stderr, "Error. Cant allocate memory.\n");
					fflush(stderr);
					return -2;
				   }
				   for (k = 0; k < mol.n_atoms+2; ++k)
					new_path[k] = path[k];

				   if(find_rings2(mol,new_path,atom_number,n_atom,
				      nring, &myringer, &myaromatic, &myvisited, &myrings) < 0) 
				   {
					fprintf(stderr, "Error. Rescursive fail.\n");
					fflush(stderr);
					return -5;
				   }
				   atom_number = 0;
				   free(new_path);
				  }else{
					myvisited[atom_number-1] = 2;
				  }
			        }else if( myvisited[atom_number-1] == 1){

			        k = last_link(path);
        		        path[k] = atom_number;

			        for(i2=0; i2<mol.n_atoms+2; ++i2) 
				{
			         for(j2=i2+1; j2<mol.n_atoms+2; ++j2) 
				 {
                       		   if(path[i2] == path[j2] && path[i2] != -1 && i2 != j2) 
				   {
                                	tmp_num_rings = tmp_num_rings + 1;
                                	elec = 0;
                                	aroflag = 0;
					current_ring = (*nring);
					myrings[current_ring] = (int *) calloc(sizeof(int), j2-i2+1);
					(*nring)++;
					#ifdef DEBUG
						fprintf(stderr,"\nRing from %i to %i: %i\n\n",i2,j2,j2-i2-1);
						fprintf(stderr,"Prev: %i, adding %i, current %i\n", last_n_atom-1, 
							atom_number-1, n_atom-1);
						fflush(stderr);
	                                        for(k2 = i2; k2 <= j2; ++k2)
	                                        {
							fprintf(stderr,"%i ",path[k2]-1);
						}
						fprintf(stderr,"\n");
	                                        fflush(stderr);
					#endif
					ll = 1;
					myrings[current_ring][0] = j2-i2; /* Size of the ring */
                                	for(k2 = i2; k2 < j2; ++k2) 
					{
                                          if(j2 - i2 > myringer[path[k2] - 1])
                                                myringer[path[k2] - 1] = j2 - i2;

					  myrings[current_ring][ll] = path[k2] - 1; /* Add element to the ring */
					  ++ll;

                                          if ( get_bonds(mol, path[k2], 0) > 3)
                                                aroflag = 1;

                                          if ( mol.atoms[ path[k2] - 1 ] == 1 )
                                                elec++;
                                          else if ( mol.atoms[ path[k2] - 1 ] == 6  )
                                                elec = elec + 4;
                                          else if (mol.atoms[ path[k2] - 1 ] == 3 && 
						get_bonds(mol, path[k2], 0) == 3)
                                                elec = elec + 2;
                                          else if (mol.atoms[ path[k2] - 1 ] == 3 && 
						get_bonds(mol, path[k2], 0) < 3 )
                                                elec++;
                                          else if (mol.atoms[ path[k2] - 1 ] == 6  )
                                                elec = elec + 4;
                                	}

                                	if ( aroflag == 0 &&  (elec - 2) % 4 == 0 )
                                       	 for ( k2 = i2; k2 < j2; ++k2)
                                                myaromatic[path[k2] - 1] = 1;

                                	/*return 0;*/
                        	   }
                                 }
                                } /* Aro-ring flag*/
				path[k] = -1;

			       }

			}

		}
	} 

/*	myvisited[n_atom-1] = 2;*/

	return 0;
}


/**
 *
 *      @brief Entry point of the Breadth-first search routine for finding rings in a graph
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *      @param mol MOL2 structure with the molecule
 *      @param myringer Ring flag for atoms
 *      @param myaromatic Aromatic flag for atoms
 *      @return negative number on error. 0 on success
 *
 */

int get_number_of_rings2(MOL2 mol, int **myringer, int **myaromatic, int ***myrings)
{
	int *new_path = NULL;
	int i = 0, j = 0, n = 0;
	int *ringer = NULL;
	int *aromatic = NULL;
	int **rings = NULL;
	int *visited = NULL;
	int nring = 0;
	int *nr = NULL;

	nr = &nring;

	ringer = *myringer;
	aromatic = *myaromatic;
	rings = *myrings;

	visited = (int *) calloc(sizeof(int), mol.n_atoms);

	if ( (new_path = calloc(mol.n_atoms+2, sizeof(int))) == NULL) {
		fprintf(stderr, "Error. Cant allocate memory.\n");
		fflush(stderr);
		return -2;
	}

        for( j = 0; j< mol.n_atoms; j++)
        {
		if( visited[j] != 0)
		  continue;

                if( (n = get_number_type_bonds(mol, j+1, 0, 0)) < 2)
		{
			visited[j] = 2;
			continue;
		}

		for (i = 0; i < mol.n_atoms+2; ++i)
			new_path[i] = -1;

		if(find_rings2(mol,new_path,j+1,0,nr,&ringer,&aromatic,&visited,&rings) < 0) 
		{
			fprintf(stderr, "Error. Recursive fail.\n");
			fflush(stderr);
			return -3;
		}
	}
	free(new_path);
	free(visited);
	return nring;
}

