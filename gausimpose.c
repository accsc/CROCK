/**
 *
 *	@file gausimpose.c
 *	@brief Main file of CROCK tool
 *
 *	@author Alvaro Cortes Cabrera <alvarocortesc@gmail.com>
 *	@date 08/02/2016
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


/**
        Rev 5 - December 2015

	  Support to read pockets from cGRILL
          Reads and writes grid calcs for faster screening
          Dumps DX files for visualization
	  Support for Tversky index

        Rev 4 - May 2013

          BFGS grid based potential added to speed up even more calcs
	  ElectroShape prefilter and triple score (Tc+Tc_color+Manhattan ES)

        Rev 3 - April 2013

         Grid based potential evaluation to speed up calcs

	Rev 2 - March 2013
	 
	 Added Simplex optimization together with Quasi-Newton BFGS

	Rev 1 - Initial version

          Quasi-Newton analytical version

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Read/Write compressed grids */
#include <zlib.h>
/* This is for reading molecules and other stuff */
#include <reader.c>
/* Memory clean up*/
#include <waste.c>
/* Gaussian volumne overlap and optimization routines */
#include <gaussian.c>
/* General rigid-body superimposition and diagonalization routines */
#include <superimpose.c>
/* Simplex minimization routines for Guassian stuff */
#include <simplex.c>
/* ElectroShape module */
#include <usr.c>

#define DEFAULT_MIN_ES_SHAPE 0.0

void advertise();

void spinner(int progress) 
{
    static char const spin_chars[] = "/-\\|";
    fprintf(stderr,"Pre-calculating grid [%c]",spin_chars[progress % sizeof(spin_chars)]);
    fflush(stderr);
    fprintf(stderr,"\r");
}

void progress_bar(int n, int total)
{
	
	static int const blocks = 60;
	int i = 0;
	n++;
	fprintf(stderr,"Overall progress [");
	for( i = 0; i < (int) ((float) n/ (float)total * (float) blocks); i++)
	{
		fprintf(stderr,"#");
	}
	for( i = i ; i < blocks; i++)
	{
		fprintf(stderr," ");
	}
        fprintf(stderr,"] %3i%% %i/%i",(int) ((float) n/ (float)total * 100.0f),n,total);
	fprintf(stderr,"\r");
}



/**
 *	@brief Main function
 *	@author Alvaro Cortes Cabrera
 *	@date 05/02/2013
 */
int main(argc, argv)
int argc;
char *argv[];
{

	float tanimoto = 0.0f, aa = 0.0f, bb=0.0f, ab=0.0f;
	float tanimoto_color = 0.0f, aa_color = 0.0f, bb_color = 0.0f, ab_color = 0.0f;
	MOL2 *lig_a = NULL, *lig_b = NULL;
	float force[3], torque[3];

	int i = 0, j = 0, bestani = -1, l = 0, k = 0, ll = 0;
	float com[3], rot[3], bestval = 0.0f, com_template[3];
        float bestval1 = 0.0f, bestval2 = 0.0f;

	float tensor[3][3], tmp = 0.0f;
        double d[3];
        double c[9], v[9];
	float rot_mat[3][3], rot_inv[3][3], Xrot = 0.0f, Yrot = 0.0f, Zrot = 0.0f, com2[3];
	float *x = NULL, *y = NULL, *z = NULL;
        float *xd = NULL, *yd = NULL, *zd = NULL;

	int min_flag = 0;

        float simplex[7][7];
        double *energies = NULL;
        double *simplex_energies = NULL;
        int real_lig_atoms = 0;
        double initial_energy = 0.0;
        double Emin = 0.0, Eold = 0.0, Enew = 0.0;
        int simplex_steps = 0;
        float delta_simplex = 0.0f;
        int bestsimplex = 0;
        float randy = 0.0f;
        float translation_step = 0.0f;
        float rotation_step = 0.0f;

	/* Grid-based stuff */

        float ****grids; /* 4D dynamic array for n grids */
	int x_points = 0, y_points = 0, z_points = 0;
	float spacing = 0.5f, min_grid[3], max_grid[3];
	float point[3], min_x = 9999.9f, min_y = 9999.9f, min_z = 9999.9f;
	float max_x = -9999.9f, max_y = -9999.9f, max_z = -9999.9f;
	int *types = NULL, **types2 = NULL;
        static char *output_grids_names[12] = {"g1.grd","g2.grd","g3.grd","g4.grd","g5.grd","g6.grd","g7.grd","g8.grd","g9.grd","g10.grd","g11.grd","g12.grd"};
        static char *output_names[12] = {"C_vol.dx","O_vol.dx","N_vol.dx","P_vol.dx","S_vol.dx","F_vol.dx","OH_color.dx","Nitro_color.dx","Apolar_color.dx","N_color.dx","O_color.dx","ring_color.dx"};

        gzFile *output_grids[12];
        char number[9];
	char *buffer = NULL;
        char keyword[20];
        int counter_points = 0;
        float *all_points = NULL;
	FILE *output_dx = NULL;
	float min_grid_orig[3], max_grid_orig[3];

	/* ElectroShape stuff */
	float current_dist = 0.0f;
	float max_distance = 0.0f;
	float distance_cutoff = 0.0f;


	advertise();

        srand(time(0));

	com[0] = com[1] = com[2] = 0.0f;

        if( argc < 3)
        {
		fprintf(stderr,"Usage: %s template.pdb candidate.pdb <min type>\n",argv[0]);
		fprintf(stderr,"\ttemplate.pdb -> Rigid molecule\n");
		fprintf(stderr,"\tcandidate.pdb -> Molecule to superimpose\n");
		fprintf(stderr,"\tMin type -> *Optional*. 0 for Nelder–Mead (Default). 1 for BFGS\n\n");
		exit(-1);
        }


        PDB_reader(&lig_a, argv[1], 0);
        PDB_reader(&lig_b, argv[2], 0);



        for( i = 0; i < lig_b->nconformers; i++)
        {
                set_default_conformer(&lig_b, i);
                current_dist = USR_compare_molecules(lig_a, lig_b);
                if( current_dist > max_distance)
                {
                        max_distance = current_dist;
                }
        }

	distance_cutoff = DEFAULT_MIN_ES_SHAPE;

	if( max_distance < distance_cutoff)
	{
	        cleanup(&lig_a);
	        cleanup(&lig_b);
		fprintf(stderr,"All solutions are below the ElectroShape distance cutoff\n");
		fprintf(stderr,"Best distance: %f < %f\n",max_distance,distance_cutoff);
		exit(0);
	}

	fprintf(stderr,"ElectroShape filter best distance: %f\n",max_distance);


	if( argc > 3)
	{
	 min_flag = atoi(argv[3]);
	}

	if( min_flag == 0)
	{
		fprintf(stderr,"\nNelder–Mead optimization selected\n");
	}else{
                fprintf(stderr,"\nBroyden–Fletcher–Goldfarb–Shanno optimization selected\n");
	}

	types = (int *) calloc(sizeof(int),lig_b->n_atoms);
        types2 = (int **) calloc(sizeof(int *),lig_b->n_atoms);

        for( i = 0; i < lig_b->n_atoms; i++)
        {
		types2[i] = (int *) calloc(sizeof(int),6);
	}

	for( i = 0; i < lig_b->n_atoms; i++)
	{
		types[i] = -1;
		if( lig_b->atoms[i] == 1)
		 types[i] = 0;
                else if( lig_b->atoms[i] == 2)
                 types[i] = 1;
                else if( lig_b->atoms[i] == 3)
                 types[i] = 2;
                else if( lig_b->atoms[i] == 4)
                 types[i] = -1;
                else if( lig_b->atoms[i] == 5)
                 types[i] = 3;
                else if( lig_b->atoms[i] == 6)
                 types[i] = 4;
                else if( lig_b->atoms[i] == 10)
                 types[i] = 5;
                else 
                 types[i] = 0;

	}

        for( i = 0; i < lig_b->n_atoms; i++)
        {
                        if( lig_b->atoms[i] == 3 && lig_b->aromatic[i] != 0 &&  get_bonds_p(lig_b,i+1,0) > 2 && has_hydrogens(lig_b,i) == 0 && lig_b->ringer[i] > 0)
                        {
/*                              fprintf(stderr,"Got a aromatic cation for a: %i. %i %i %i %i %i\n",i,a->atoms[i],a->aromatic[i],get_bonds_p(a,i+1,0),a->ringer[i],has_hydrogens(a,i));*/
                                types2[i][3] = 1;
                        }

                        if( lig_b->gaff_types[i] == OH || lig_b->gaff_types[i] == OW || ( lig_b->atoms[i] == 3 && has_hydrogens(lig_b,i) > 0))
                        {
                                        types2[i][0] = 1;
                        }

                        if( lig_b->gaff_types[i] == O || lig_b->gaff_types[i] == OH || ( lig_b->atoms[i] == 3 && has_hydrogens(lig_b,i) == 0))
                        {
                                if( gaus_is_carboxylate(lig_b,i) < 4)
                                {
                                        types2[i][4] = 1;
                                }else
                                        types2[i][1] = 1;
                        }

                        if( lig_b->ringer[i] > 0)
                        {
                                types2[i][5] = 1;
                        }

        }



        simplex_energies = (double*)calloc( sizeof(double ), 6 + 1);
        delta_simplex = 0.1;
        rotation_step = 24.0f;
        translation_step = 1.0;

        translation_step *= 0.5;
        rotation_step *= 0.5;


	/* Final-best coordinates */
	x = (float *) calloc(sizeof(float),lig_b->n_atoms);
        y = (float *) calloc(sizeof(float),lig_b->n_atoms);
        z = (float *) calloc(sizeof(float),lig_b->n_atoms);

	/* Backup coordinates for rotation pi around three axes */
        xd = (float *) calloc(sizeof(float),lig_b->n_atoms);
        yd = (float *) calloc(sizeof(float),lig_b->n_atoms);
        zd = (float *) calloc(sizeof(float),lig_b->n_atoms);

	/* Remove COM from template */
	for( i = 0; i < lig_a->n_atoms; i++)
	{
		com[0] += lig_a->x[i];
                com[1] += lig_a->y[i];
                com[2] += lig_a->z[i];
	}

	com[0] /= (float) lig_a->n_atoms;
        com[1] /= (float) lig_a->n_atoms;
        com[2] /= (float) lig_a->n_atoms;


        for( i = 0; i < lig_a->n_atoms; i++)
        {
                lig_a->x[i] -= com[0];
                lig_a->y[i] -= com[1];
                lig_a->z[i] -= com[2];
        }

	com2[0] = com[0];
	com2[1] = com[1];
	com2[2] = com[2];

	/* Calculate of moment of inertia tensor */
        tensor[0][0] = 0.0f;
        tensor[0][1] = 0.0f;
        tensor[0][2] = 0.0f;

        tensor[1][0] = 0.0f;
        tensor[1][1] = 0.0f;
        tensor[1][2] = 0.0f;

        tensor[2][0] = 0.0f;
        tensor[2][1] = 0.0f;
        tensor[2][2] = 0.0f;


        for( i = 0; i < lig_a->n_atoms; i++)
        {
                tensor[0][0] += ( (lig_a->y[i]*lig_a->y[i]) + (lig_a->z[i]*lig_a->z[i])); /*Ixx*/
                tmp = lig_a->x[i]*lig_a->y[i];
                tensor[0][1] -= tmp; /* Ixy */
                tensor[1][0] -= tmp; /* Iyx */
                tmp = lig_a->x[i]*lig_a->z[i];
                tensor[0][2] -= tmp; /* Ixz */
                tensor[2][0] -= tmp; /* Izx */
                tensor[1][1] += ( (lig_a->x[i]*lig_a->x[i]) + (lig_a->z[i]*lig_a->z[i])); /*Iyy*/
                tmp = lig_a->y[i]*lig_a->z[i];
                tensor[1][2] -= tmp; /* Iyz */
                tensor[2][1] -= tmp; /* Izy */
                tensor[2][2] += ( (lig_a->x[i]*lig_a->x[i]) + (lig_a->y[i]*lig_a->y[i])); /*Izz*/
        }

	/* Diagonalization */
        c[0] = tensor[0][0];
        c[1] = tensor[0][1];
        c[2] = tensor[0][2];
        c[3] = tensor[1][0];
        c[4] = tensor[1][1];
        c[5] = tensor[1][2];
        c[6] = tensor[2][0];
        c[7] = tensor[2][1];
        c[8] = tensor[2][2];
        jacobi(3, c, d, v);

	/* Build rotation matrix since I = R * (lambda*I) * Rt */
	/* and R is my rotation matrix */
        rot_mat[0][0] = v[0];
        rot_mat[0][1] = v[1];
        rot_mat[0][2] = v[2];

        rot_mat[1][0] = v[3];
        rot_mat[1][1] = v[4];
        rot_mat[1][2] = v[5];

        rot_mat[2][0] = v[6];
        rot_mat[2][1] = v[7];
        rot_mat[2][2] = v[8];

	rot_inv[0][0] = v[0];
	rot_inv[1][0] = v[1];
        rot_inv[2][0] = v[2];

        rot_inv[0][1] = v[3];
        rot_inv[1][1] = v[4];
        rot_inv[2][1] = v[5];

        rot_inv[0][2] = v[6];
        rot_inv[1][2] = v[7];
        rot_inv[2][2] = v[8];


#ifdef DEBUG
        fprintf(stderr,"Eigenvals: %f %f %f\n",d[0],d[1],d[2]);
#endif
                
	/* Apply rotation to align axes with x,y,z */
        Xrot = Yrot = Zrot = 0.0f;
        for ( i = 0; i < lig_a->n_atoms; ++i)
        {       
/*	fprintf(stderr,"%i %s\n",i,let[lig_a->gaff_types[i]-1]);*/
         Xrot = (rot_mat[0][0] * lig_a->x[i] + rot_mat[1][0] * lig_a->y[i] + rot_mat[2][0] * lig_a->z[i]);        
         Yrot = (rot_mat[0][1] * lig_a->x[i] + rot_mat[1][1] * lig_a->y[i] + rot_mat[2][1] * lig_a->z[i]);
         Zrot = (rot_mat[0][2] * lig_a->x[i] + rot_mat[1][2] * lig_a->y[i] + rot_mat[2][2] * lig_a->z[i]);
         lig_a->x[i] = Xrot;
         lig_a->y[i] = Yrot;
         lig_a->z[i] = Zrot;
        }

	printf("MODEL 1\n");
 	dump_pdb_conservative(lig_a);
	printf("ENDMDL\n");

	/* Get self-overlap for template */
        aa = -get_volumen_intersection(lig_a,&lig_a,0,0);
        /* Get self-overlap-color for template */
	aa_color = -get_color_intersection(lig_a,&lig_a,0,0);

	
	min_x = max_x = lig_a->x[0];
        min_y = max_y = lig_a->y[0];
        min_z = max_z = lig_a->z[0];

        for ( i = 0; i < lig_a->n_atoms; i++)
	{
		if( lig_a->x[i] < min_x)
			min_x = lig_a->x[i];
                if( lig_a->y[i] < min_y)
                        min_y = lig_a->y[i];
                if( lig_a->z[i] < min_z)
                        min_z = lig_a->z[i];

                if( lig_a->x[i] > max_x)
                        max_x = lig_a->x[i];
                if( lig_a->y[i] > max_y)
                        max_y = lig_a->y[i];
                if( lig_a->z[i] > max_z)
                        max_z = lig_a->z[i];

	}

	/* 6.0A around ligand */
	min_x -= 3.0f;
        min_y -= 3.0f;
        min_z -= 3.0f;

        max_x += 3.0f;
        max_y += 3.0f;
        max_z += 3.0f;

	x_points = ceil(fabs(max_x-min_x)/spacing);
        y_points = ceil(fabs(max_y-min_y)/spacing);
        z_points = ceil(fabs(max_z-min_z)/spacing);
	
	min_grid[0] = min_x;
        min_grid[1] = min_y;
        min_grid[2] = min_z;

        max_grid[0] = max_x;
        max_grid[1] = max_y;
        max_grid[2] = max_z;


	/* Pre computation of overlaps in-grid */

	fprintf(stderr,"Allocating grids...");
	fflush(stderr);
        grids = (float ****)calloc(14,sizeof(float ***));  /* Allocate 7+7 probes */
        for( i = 0 ; i < 14; ++i)
        {
                grids[i] = (float ***) calloc(x_points,sizeof(float **));
                for( j = 0; j < x_points; ++j)
                {
                        grids[i][j] = (float **) calloc(y_points,sizeof(float *));
                        for( k = 0; k < y_points; ++k)
                        {
                                grids[i][j][k] = (float *) calloc(z_points,sizeof(float *));
                        }

                }
        }

	fprintf(stderr,"done\n");
        fflush(stderr);


	if( access( "g1.grd", F_OK) == -1 )
	{

        for( i = 0 ; i < x_points; ++i)
        {
		spinner(i);
                for( j = 0; j < y_points; ++j)
                {
                        for( k = 0; k < z_points; ++k)
                        {
                                point[0] = min_grid[0] + (i * spacing);
                                point[1] = min_grid[1] + (j * spacing);
                                point[2] = min_grid[2] + (k * spacing);


				grids[0][i][j][k] = get_volumen_intersection_at_point(lig_a, 
							point[0], point[1], point[2], 1,0);
                                grids[1][i][j][k] = get_volumen_intersection_at_point(lig_a, 
							point[0], point[1], point[2], 2,0);
                                grids[2][i][j][k] = get_volumen_intersection_at_point(lig_a, 
							point[0], point[1], point[2], 3,0);
                                grids[3][i][j][k] = get_volumen_intersection_at_point(lig_a, 
							point[0], point[1], point[2], 5,0);
                                grids[4][i][j][k] = get_volumen_intersection_at_point(lig_a, 
							point[0], point[1], point[2], 6,0);
                                grids[5][i][j][k] = get_volumen_intersection_at_point(lig_a, 
							point[0], point[1], point[2], 10,0);

				grids[6][i][j][k] = get_color_intersection_at_point(lig_a,
							point[0], point[1], point[2],0,0);
                                grids[7][i][j][k] = get_color_intersection_at_point(lig_a, 
                                                        point[0], point[1], point[2],1,0);
                                grids[8][i][j][k] = get_color_intersection_at_point(lig_a, 
                                                        point[0], point[1], point[2],2,0);
                                grids[9][i][j][k] = get_color_intersection_at_point(lig_a, 
                                                        point[0], point[1], point[2],3,0);
                                grids[10][i][j][k] = get_color_intersection_at_point(lig_a, 
                                                        point[0], point[1], point[2],4,0);
                                grids[11][i][j][k] = get_color_intersection_at_point(lig_a, 
                                                        point[0], point[1], point[2],5,0);


			}
		}

	}

	fprintf(stderr,"\nGrid completed!\n");

        fprintf(stderr,"\nWritting grids... ");
	fflush(stderr);

        if( (all_points = (float *) calloc(x_points*y_points*z_points,sizeof(float))) == NULL)
        {
                fprintf(stderr,"ERROR-MAIN%i: Cannot allocate memory.\n\n",__LINE__);
                fprintf(stderr," *** END ERROR *** \n");
                fflush(stderr);
                return -2;
        }

	for( l = 0; l < 12 ; l++)
	{

	        if( (output_grids[l] = gzopen(output_grids_names[l],"w")) == NULL)
	        {
	                fprintf(stderr,"ERROR-MAIN%i: Cannot open grid file to write.\n\n",__LINE__);
	                fprintf(stderr," *** END ERROR *** \n");
	                fflush(stderr);
	                return -3;
	        }



	        gzprintf(output_grids[l],"GAUS %i %i %i %8.3f%8.3f%8.3f%8.3f\n",x_points,y_points,z_points,spacing,min_grid[0],min_grid[1],min_grid[2]);

		counter_points = 0;
	        for( i = 0 ; i < x_points; ++i)
	        {
	                for( j = 0; j < y_points; ++j)
	                {
	                        for( k = 0; k < z_points; ++k)
	                        {
	                                gzprintf(output_grids[l],"%8.3f", grids[l][i][j][k]);
	                                all_points[counter_points] = grids[l][i][j][k];
	                                ++counter_points;
	                        }
	
	                }
	        }

	        gzclose(output_grids[l]);


	        if( (output_dx = fopen(output_names[l],"w")) == NULL)
	        {
	                fprintf(stderr,"ERROR-MAIN%i: Cannot open grid file to write.\n\n",__LINE__);
	                fprintf(stderr," *** END ERROR *** \n");
	                fflush(stderr);
	                return -3;
	        }

	        fprintf(output_dx,"object 1 class gridpositions counts %i %i %i\n", x_points, y_points, z_points);
	        fprintf(output_dx,"origin %8.3f %8.3f %8.3f\n",min_grid[0],min_grid[1],min_grid[2]);
	        fprintf(output_dx,"delta %8.3f 0.0 0.0\n",spacing);
	        fprintf(output_dx,"delta 0.0 %8.3f 0.0\n", spacing);
	        fprintf(output_dx,"delta 0.0 0.0 %8.3f\n", spacing);
	        fprintf(output_dx,"object 2 class gridconnections count %i %i %i\n",x_points, y_points, z_points);
	        fprintf(output_dx,"object 3 class array type double rank 0 items %i data follows\n",x_points*y_points*z_points);

	        for( i = 0; i < counter_points; ++i)
	        {
	                 fprintf(output_dx,"%18.3f",all_points[i]);

	                 if( (i+1)%3  == 0)
	                 {
	                        fprintf(output_dx,"\n");
	                 }
	        }

	        if( counter_points%3 != 0)
	              fprintf(output_dx,"\n");


	        fprintf(output_dx,"attribute \"dep\" string \"positions\"\n");
	        fprintf(output_dx,"object \"regular positions regular connections\" class field\n");
	        fprintf(output_dx,"component \"positions\" value 1\n");
	        fprintf(output_dx,"component \"connections\" value 2\n");
	        fprintf(output_dx,"component \"data\" value 3\n");

	        fclose(output_dx);
	}
        fprintf(stderr,"done\n");
        fflush(stderr);

	free(all_points);

	}else{ /* I found some grids here */

        fprintf(stderr, "Reading grids ...");
        fflush(stderr);

        for ( l = 0; l < 12; ++l) {
                if ( (output_grids[l] = gzopen(output_grids_names[l], "rb")) == NULL) {
                        fprintf(stderr, "Could not open the grid files.\n");
                        fflush(stderr);
                        return -1;
                }

                buffer = (char*)calloc(1024, sizeof(char));
                gzgets(output_grids[l], buffer, MAX_BUFFER);

                sscanf(buffer, "%s %i %i %i %f %f %f %f", keyword, &x_points, &y_points, &z_points, &spacing, &min_grid[0], &min_grid[1], &min_grid[2]);
                free(buffer);

                max_grid[0] = min_grid[0] + (spacing * x_points);
                max_grid[1] = min_grid[1] + (spacing * y_points);
                max_grid[2] = min_grid[2] + (spacing * z_points);

                if ( !(strstr(keyword, "GAUS") > 0)) {
                        printf("Keyword fail. Corruption in file.\n");
                        return -1;
                }
                for ( i = 0; i < x_points; ++i) {
                        for ( j = 0; j < y_points; ++j) {
                                for ( k = 0; k < z_points; ++k) {
                                        grids[l][i][j][k] = 0.0f;
                                        gzread( output_grids[l], number, 8 );
                                        grids[l][i][j][k] = atof(number);

                                }
                        }

                }

                gzclose(output_grids[l]);

	}

        fprintf(stderr, "done\n");
	fflush(stderr);
	}


	/* Flexibility issue is handled before with pre-computation of conformers */
        for ( ll = 0; ll < lig_b->nconformers; ++ll)  
        {
	progress_bar(ll,lig_b->nconformers);
/*	fprintf(stderr," ******* LIGAND %i ********\n",ll);*/
        set_default_conformer(&lig_b, ll);


	/* Intialize gradient for mobile molecule */
        for( i = 0; i < lig_b->n_atoms; i++)
        {
		lig_b->grads_X[i] = 0.0f;
                lig_b->grads_Y[i] = 0.0f;
                lig_b->grads_Z[i] = 0.0f;
        }


        tensor[0][0] = 0.0f;
        tensor[0][1] = 0.0f;
        tensor[0][2] = 0.0f;

        tensor[1][0] = 0.0f;
        tensor[1][1] = 0.0f;
        tensor[1][2] = 0.0f;

        tensor[2][0] = 0.0f;
        tensor[2][1] = 0.0f;
        tensor[2][2] = 0.0f;

	/* COM removal for mobile molecule */
        com[0] = com[1] = com[2] = 0.0f;
        for( i = 0; i < lig_b->n_atoms; i++)
        {       
                com[0] += lig_b->x[i];
                com[1] += lig_b->y[i];
                com[2] += lig_b->z[i];
        }       
        com[0] /= (float) lig_b->n_atoms;
        com[1] /= (float) lig_b->n_atoms;
        com[2] /= (float) lig_b->n_atoms;
        for( i = 0; i < lig_b->n_atoms; i++)
        {
		lig_b->x[i] -= com[0];
                lig_b->y[i] -= com[1];
                lig_b->z[i] -= com[2];
        }


	/* Calculate of moment of inertia tensor for mobile molecule */
        for( i = 0; i < lig_b->n_atoms; i++)
        {
		tensor[0][0] += ( (lig_b->y[i]*lig_b->y[i]) + (lig_b->z[i]*lig_b->z[i])); /*Ixx*/
		tmp = lig_b->x[i]*lig_b->y[i];
		tensor[0][1] -= tmp; /* Ixy */
		tensor[1][0] -= tmp; /* Iyx */
		tmp = lig_b->x[i]*lig_b->z[i];
 		tensor[0][2] -= tmp; /* Ixz */
                tensor[2][0] -= tmp; /* Izx */
		tensor[1][1] += ( (lig_b->x[i]*lig_b->x[i]) + (lig_b->z[i]*lig_b->z[i])); /*Iyy*/
		tmp = lig_b->y[i]*lig_b->z[i];
		tensor[1][2] -= tmp; /* Iyz */
		tensor[2][1] -= tmp; /* Izy */
		tensor[2][2] += ( (lig_b->x[i]*lig_b->x[i]) + (lig_b->y[i]*lig_b->y[i])); /*Izz*/
        }


#ifdef DEBUG
	fprintf(stderr,"%f %f %f\n",tensor[0][0],tensor[0][1],tensor[0][2]);
        fprintf(stderr,"%f %f %f\n",tensor[1][0],tensor[1][1],tensor[1][2]);
        fprintf(stderr,"%f %f %f\n",tensor[2][0],tensor[2][1],tensor[2][2]);
#endif

	/* Diagonalize tensor */
	c[0] = tensor[0][0];
	c[1] = tensor[0][1];
        c[2] = tensor[0][2];
        c[3] = tensor[1][0];
        c[4] = tensor[1][1];
        c[5] = tensor[1][2];
        c[6] = tensor[2][0];
        c[7] = tensor[2][1];
        c[8] = tensor[2][2];
	jacobi(3, c, d, v);

	/* Build rotation matrix */
	rot_mat[0][0] = v[0];
        rot_mat[0][1] = v[1];
        rot_mat[0][2] = v[2];

        rot_mat[1][0] = v[3];
        rot_mat[1][1] = v[4];
        rot_mat[1][2] = v[5];

        rot_mat[2][0] = v[6];
        rot_mat[2][1] = v[7];
        rot_mat[2][2] = v[8];

#ifdef DEBUG
	fprintf(stderr,"Eigenvals: %f %f %f\n",d[0],d[1],d[2]);
#endif

	/* Apply rotation */
        Xrot = Yrot = Zrot = 0.0f;
        for ( i = 0; i < lig_b->n_atoms; ++i) 
	{
         Xrot = (rot_mat[0][0] * lig_b->x[i] + rot_mat[1][0] * lig_b->y[i] + rot_mat[2][0] * lig_b->z[i]);
         Yrot = (rot_mat[0][1] * lig_b->x[i] + rot_mat[1][1] * lig_b->y[i] + rot_mat[2][1] * lig_b->z[i]);
         Zrot = (rot_mat[0][2] * lig_b->x[i] + rot_mat[1][2] * lig_b->y[i] + rot_mat[2][2] * lig_b->z[i]);
         lig_b->x[i] = Xrot;
         lig_b->y[i] = Yrot;
         lig_b->z[i] = Zrot;
        }

	/* Backup coordinates */
        for( j = 0; j < lig_b->n_atoms; j++)
        {
                xd[j] = lig_b->x[j];
                yd[j] = lig_b->y[j];
                zd[j] = lig_b->z[j];
        }

	/* ?! Not useful anymore but to avoid changing var names by zeros since template is at origin  */
        com_template[0] = com_template[1] = com_template[2] = 0.0f;
	/* Get self-score */
        bb = -get_volumen_intersection(lig_b,&lig_b,0,0);
	bb_color = -get_color_intersection(lig_b,&lig_b,0,0);

	for( k = 0; k < 4; k++) /* Standard rotations loop */
	{

	        for( j = 0; j < lig_b->n_atoms; j++)
	        {
	                lig_b->x[j] = xd[j];
	                lig_b->y[j] = yd[j];
	                lig_b->z[j] = zd[j];
	        }

		if( k == 0){ /* First "standard" orientation */
		}else if( k ==1 ){ /* Rotate pi radians x-axis */
			rot[0] = -3.14159f;
			rot[1] = 0.0f;
			rot[2] = 0.0f;
			transformate_rigid(&lig_b, com_template, rot, 1.0f);
                }else if( k ==2 ){ /* Rotate pi radians y-axis */
                        rot[1] = -3.14159f;
                        rot[0] = 0.0f;
                        rot[2] = 0.0f;
                        transformate_rigid(&lig_b, com_template, rot, 1.0f);
                }else if( k ==1 ){ /* Rotate pi radians z-axis */
                        rot[2] = -3.14159f;
                        rot[1] = 0.0f;
                        rot[0] = 0.0f;
                        transformate_rigid(&lig_b, com_template, rot, 1.0f);
		}

		/* Try to find global minimum (since score is negative) in 200 steps */
/*		gau_minimizer_bfgs(lig_a, &lig_b, 200);*/
		if( min_flag == 1 ){
			gau_minimizer_bfgs_ingrid(&lig_b, 200, grids, types, types2, min_grid, max_grid, x_points, y_points, z_points, spacing);
		}else{

                        for ( l = 0; l < 7; ++l)
                                simplex[0][l] = 0.0f;

/*                        simplex_energies[0] =  get_simplex_energy(lig_b, simplex[0], lig_a);*/
                        simplex_energies[0] =  get_simplex_energy_ingrid(lig_b, simplex[0], grids, types, types2, min_grid, max_grid, x_points, y_points, z_points, spacing);
                        for ( l = 1; l < 7; ++l)
                                simplex_energies[l] = 0.0f;

                        while ( 1 == 1) { /* A little kitty died because of this */

                                Eold = simplex_energies[bestsimplex];
                                for ( j = 0; j < 6; ++j)
                                        simplex[0][j] = simplex[bestsimplex][j];
                                simplex_energies[0] = Eold;
                                for ( j = 1; j < 7; ++j) {
                                        for ( l = 0; l < 6; ++l) {
                                                randy = 2.0 * (gen_rand_float() - 0.5);
                                                if ( l <= 2)
                                                        simplex[j][l] = simplex[0][l] + translation_step * randy;
                                                else{
                                                        randy = randy * rotation_step;
                                                        if ( randy < 0) randy = 360.0f + randy;
                                                        simplex[j][l] = simplex[0][l] + randy;
                                                }
                                        }
/*                                        simplex_energies[j] = get_simplex_energy(lig_b, simplex[j], lig_a);*/
                        		  simplex_energies[j] =  get_simplex_energy_ingrid(lig_b, simplex[j], grids, types, types2, min_grid, max_grid, x_points, y_points, z_points, spacing);


                                }

/*                                go_simplex_go(6, simplex, &simplex_energies, lig_b, lig_a);*/
                        	  go_simplex_go_ingrid(6, simplex, &simplex_energies, lig_b, grids, types, types2, min_grid, max_grid, x_points, y_points, z_points, spacing);
                                simplex_steps++;

                                Emin = 99999.9;
                                bestsimplex = 0;
                                for ( l = 0; l < 7; ++l) {
                                        if ( simplex_energies[l] < Emin) {
                                                bestsimplex = l;
                                                Emin = simplex_energies[l];
                                        }
                                }
/*                                Enew = get_simplex_energy(lig_a, simplex[bestsimplex], lig_a);*/
                        	  Enew = get_simplex_energy_ingrid(lig_b, simplex[bestsimplex], grids, types, types2, min_grid, max_grid, x_points, y_points, z_points, spacing);

                                Enew = Emin;

                                if ( simplex_steps > 20 || fabs(Eold - Enew) < 0.01)
				{
                                        break;
				}
                        }
                        for ( l = 0; l < lig_b->n_atoms; ++l) {
                                lig_b->x[l] = xd[l];
                                lig_b->y[l] = yd[l];
                                lig_b->z[l] = zd[l];
                        }
                        transformate_mol(&lig_b, simplex[bestsimplex]);

		} /* Optimization choice */


		/* Get candidate overlap */
	        ab = -get_volumen_intersection(lig_a,&lig_b,0,0);
		ab_color = -get_color_intersection(lig_a,&lig_b,0,0);

		/* Calculate both tanimoto scores */
		tanimoto = ab/(aa+bb-ab);
		tanimoto_color = ab_color/(aa_color+bb_color-ab_color);

/*                tanimoto = ab/( (aa*0.95)+(bb*0.05)+ab);
                tanimoto_color = ab_color/((aa_color*0.95)+(bb_color*0.05)+ab_color);*/

		
		if ( (tanimoto_color+tanimoto) > bestval)
		{
		    bestani = ll;
		    bestval = tanimoto_color+tanimoto; /* Save "combo" value */
                    bestval1 = tanimoto_color;
                    bestval2 = tanimoto;
		    /* Save best coordinates so far */
		    for( j = 0; j < lig_b->n_atoms; j++)
		    {
			x[j] = lig_b->x[j];
	                y[j] = lig_b->y[j];
	                z[j] = lig_b->z[j];
		    }
		}
		/* Clean up the gradients mesh in lig_b structure */
		for( j = 0; j < lig_b->n_atoms; j++)
		{
	                lig_b->grads_X[j] = 0.0f;
	                lig_b->grads_Y[j] = 0.0f;
	                lig_b->grads_Z[j] = 0.0f;
		}

	} /* End standard rotations loop */


	} /* End Conformers loop */

	
	/* We dont need the grids anymore */
        for( i = 0 ; i < 14; ++i)
        {
                for( j = 0; j < x_points; ++j)
                {
                        for( k = 0; k < y_points; ++k)
                        {
				free(grids[i][j][k]);
                        }
			free(grids[i][j]);

                }
		free(grids[i]);
        }
	free(grids);



	/* Restore best coordinates */
        for( j = 0; j < lig_b->n_atoms; j++)
        {
              lig_b->x[j] = x[j];
              lig_b->y[j] = y[j];
              lig_b->z[j] = z[j];
        }

        printf("MODEL 1\n");
        dump_pdb(lig_b);
        printf("ENDMDL\n");

	/* Apply initial transformation of the template to restore the frame */
        Xrot = Yrot = Zrot = 0.0f;
        for ( i = 0; i < lig_b->n_atoms; ++i)
        {
         Xrot = (rot_inv[0][0] * lig_b->x[i] + rot_inv[1][0] * lig_b->y[i] + rot_inv[2][0] * lig_b->z[i]);
         Yrot = (rot_inv[0][1] * lig_b->x[i] + rot_inv[1][1] * lig_b->y[i] + rot_inv[2][1] * lig_b->z[i]);
         Zrot = (rot_inv[0][2] * lig_b->x[i] + rot_inv[1][2] * lig_b->y[i] + rot_inv[2][2] * lig_b->z[i]);
         lig_b->x[i] = Xrot;
         lig_b->y[i] = Yrot;
         lig_b->z[i] = Zrot;
        }

        for( j = 0; j < lig_b->n_atoms; j++)
        {
              lig_b->x[j] += com2[0];
              lig_b->y[j] += com2[1];
              lig_b->z[j] += com2[2];
        }




	/* Print result! */
	printf("MODEL 1\n");
        printf("REMARK CROCK RESULT - rev5/Dec2015\n");
        /** This REMARK is for integration within VSDMIP */
        printf("REMARK numConformer: %8i\n",bestani+1);
	printf("REMARK Tc: %f\n",bestval);
	printf("REMARK Color: %f\n",bestval1);
	printf("REMARK Shape: %f\n",bestval2);
	printf("REMARK TriScore: %f\n",bestval+USR_compare_molecules(lig_a, lig_b));
        dump_pdb(lig_b);
	printf("ENDMDL\n");
	fflush(stdout);
        fprintf(stderr,"\nBest Tc: %f\n",bestval);
	fprintf(stderr,"TriScore: %f\n",bestval+USR_compare_molecules(lig_a, lig_b));

	free(simplex_energies);
	/* Free backup coordinates */
	free(x); free(y); free(z);
        free(xd); free(yd); free(zd);

	free(types); 
	for(i = 0; i < lig_b->n_atoms; i++)
		free(types2[i]);

	free(types2);

	/* Clean structures */
	cleanup(&lig_a);
	cleanup(&lig_b);

}

void advertise()
{
        fprintf(stderr, "CROCK rev5 Another program to superimpose ligands using gaussians\n");
        fprintf(stderr, "Brought to you by Department of Pharmacology UAH:\n\n ");
        fprintf(stderr, "\tAlvaro Cortes Cabrera <alvarocortesc@gmail.com>\n");
        fprintf(stderr, "\tFederico Gago <federico.gago@uah.es>\n\n");
        fprintf(stderr, "Have fun!\n\n");
}

