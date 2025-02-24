/**
 *
 *      @file groups_chg.c
 *      @brief Charged groups detection
 *
 *      @author Alvaro Cortes Cabrera <acortes@cbm.uam.es>
 *      @date 06/06/2016
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


int detect_groups_heuristic(MOL2 *mol, int ***mytypes)
{

	int i = 0, j = 0, k = 0, c = 0, m = 0, n = 0, l = 0, nh =0;
	float x = 0.0f, y = 0.0f, z = 0.0f;
	int idx = 0;
	int type = 9;
	int **neig = NULL;
	int **neig_type = NULL;
	int *neig_n = NULL;
	int *flag_hydrogen = NULL;
	int **types = NULL;
	types = *mytypes;
	int *hydros = NULL;


	if( mol->n_atoms <= 0)
	{
		fprintf(stderr,"Error. No atoms in the structure\n");
		fflush(stderr);
		return -1;
	}

	/* Convert list of atoms to neigbours */

	neig = (int **) calloc(sizeof(int *), mol->n_atoms);
	neig_type = (int **) calloc(sizeof(int *), mol->n_atoms);
	neig_n = (int *) calloc(sizeof(int), mol->n_atoms);
	flag_hydrogen = (int *) calloc(sizeof(int), mol->n_atoms);

        for( i = 0; i < mol->n_atoms; i++)
	{
		neig[i] = (int *)calloc(sizeof(int), 12);
		neig_type[i] = (int *)calloc(sizeof(int), 12);
	}
	
        for( i = 0; i < mol->n_atoms; i++)
        {
                k = 0;
                for (j = 0; j < mol->n_bonds; ++j) {
                        if ( mol->bond_a1[j] == (i + 1)) {
                                neig[i][k] = mol->bond_a2[j] - 1;
                                neig_type[i][k] = mol->bonds[j];
				if( mol->atoms[neig[i][k]] == 4)
				  flag_hydrogen[i] = flag_hydrogen[i] + 1;
                                k++;
                        }else if ( mol->bond_a2[j] == (i + 1)) {
                                neig[i][k] = mol->bond_a1[j] - 1;
                                neig_type[i][k] = mol->bonds[j];
                                if( mol->atoms[neig[i][k]] == 4)
                                  flag_hydrogen[i] = flag_hydrogen[i] + 1;
                                k++;
                        }
                }
		neig_n[i] = k;
	}


	hydros = (int *) calloc( sizeof(int), mol->n_atoms);

        for( i = 0; i < mol->n_atoms; i++)
        {
		nh = flag_hydrogen[i];

		n = 0;
                for (j = 0; j < neig_n[i]; ++j)
                {
			if( neig_type[i][j] == 2)
			{
				++n;
			}
                }


		if( (mol->atoms[i] == 1 || (mol->atoms[i] == 6 && (k-nh) == 2 && n == 0) || mol->atoms[i] == 9 || mol->atoms[i] == 8) && mol->ringer[i] == 0)
		{
			hydros[i] = (neig_n[i]-nh);
		}
	}

        for( i = 0; i < mol->n_atoms; i++)
        {

		n = 0;
                for (j = 0; j < neig_n[i]; ++j)
                {
                        if( mol->atoms[neig[i][j]] != 1 && mol->atoms[neig[i][j]] != 4 && !(mol->atoms[neig[i][j]] == 6 && hydros[neig[i][j]] == 1))
                        {
                                ++n;
                        }
                }

		if( n > 0)
			hydros[i] = 0;
	}


        for( i = 0; i < mol->n_atoms; i++)
        {
		if( hydros[i] > 0)
		{
			types[i][4] = 1;
		}
		
	}
	

	/* Rest of polar groups */
	for( i = 0; i < mol->n_atoms; i++)
	{

		if( mol->atoms[i] == 3) 
		{
			if(mol->aromatic[i] != 0 &&  neig_n[i] > 2 && flag_hydrogen[i] == 0 && mol->ringer[i] > 0)
			        types[i][2] = 1;
                        if( mol->gaff_types[i] == N4)
                                types[i][2] = 1;
			if( flag_hydrogen[i] > 0 )
                                types[i][0] = 1;
			else if( flag_hydrogen[i] == 0  && is_planar(mol,i) != 0)
			     types[i][1] = 1;

                	if( mol->gaff_types[i] == NH){
	                        c = 0;
		                for( j = 0; j < neig_n[i]; j++)
	       	                {
	       	                         if( mol->gaff_types[neig[i][j]] == N)
	       	                         {
	       	                                 c++;
	       	                         }
	       	                }

	                        if( c > 0)
	       	                {
	                                types[i][2] = 1;
	                        }

			}
	                if( mol->gaff_types[i] != N)
			{

	                        l = 0;
	                        n = 0;

	                        for( j = 0; j < neig_n[i]; j++)
	                        {
	                                if( mol->atoms[neig[i][j]] == 6)
	                                {
	                                        c = 0;
	                                        for (m = 0; m < neig_n[neig[i][j]]; ++m) {
	                                                if ( neig_type[neig[i][j]][m] == 2 ) {
	                                                 c++;
	                                                }
	                                        }
	                                        if( c ==2 )
	                                          ++l;
	                                }else if( mol->atoms[neig[i][j]] == 4){
	                                        ++n;
	                                }
	
	                        }
	
	                        if( l >= 1 && n >= 1 && neig_n[i]-n-1 == 0)
	                        {
	
	                                types[i][3] = 1;
	                        }
	
			}
	

		}else if( mol->atoms[i] == 2){
			if( mol->gaff_types[i] == OH || mol->gaff_types[i] == OW)
			{
		                types[i][0] = 1;
                                types[i][1] = 1;
			}else if( mol->gaff_types[i] == OS ){
                                types[i][1] = 1;
			}else if( mol->gaff_types[i] == O)
			{
                                types[i][1] = 1;
			        if( gaus_is_carboxylate(mol,i) == 0)
			                types[i][3] = 1;
			}
		}else if( mol->atoms[i] == 1){
		
			if( mol->gaff_types[i] == C1)
			{
	                        c = 0;
	                        for( j = 0; j < k; j++)
	                        {
	                                if( mol->gaff_types[neig[i][j]] == N1)
	                                {
	                                        c++;
	                                }
	                        }
	
	                        if( c > 0)
				{
					types[i][2] = 1;
				}
			}

				c = 0;
				l = 0;
				n = 0;
	                        for( j = 0; j < neig_n[i]; j++)
	                        {
	                                if( mol->atoms[neig[i][j]] == 3) 
	                                {
						++l;
	                                        c = 0;
	                                        for (m = 0; m < mol->n_bonds; ++m) {
	                                                if ( (mol->bond_a1[m] == (neig[i][j] + 1) && mol->atoms[mol->bond_a2[m]-1] == 4) || (mol->bond_a2[m] == (neig[i][j] + 1) && mol->atoms[mol->bond_a1[m]-1] == 4 )) {
	                                                 c++;
	                                                }
	                                        }
	
	                                        if( c == 2) 
	                                        {
	                                                n++;
	                                        }
	

	                                }
	                        }


				if( l >= 2 && n >= 2 && k == 3)
				{
	                                types[i][2] = 1;
				}
	

		}

	}

	free(hydros);

        for( i = 0; i < mol->n_atoms; i++)
        {
		free(neig[i]);
		free(neig_type[i]);
        }
	free(neig);
	free(neig_type);
	free(neig_n);

	return 0;
}

