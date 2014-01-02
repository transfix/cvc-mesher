#include<cellqueue.h>
#include"octree.h"
#include"e_face.h"
#include"cubes.h"
#include<stdio.h>
#include"vtkMarchingCubesCases.h"
#include"pcio.h"
#include"LBIE_geoframe.h"
//#include <qgl.h>
#include<assert.h>

namespace LBIE
{

void Octree::tetrahedralize(geoframe& geofrm) {

	int x, y, z, tx, ty, tz, valid_leaf, cell_size, level, i, j, k;
	int vtx_num, intersect_id, my_vtx[4], con_id[4];
	unsigned int vtx[4], my_vertex, my_vertex_1;
	float val[8];

	for (k = 0; k < octcell_num; k++) vtx_idx_arr[k] = -1;
	for (k = 0; k < dim[0]*dim[1]*dim[2]; k++) grid_idx_arr[k] = -1;

	for (i = 0; i < leaf_num; i++ ) {

		valid_leaf = cut_array[i] ;
		level = get_level(valid_leaf) ;
		cell_size = (dim[0]-1)/(1<<level);
		octcell2xyz(valid_leaf, x, y, z, level);
		getCellValues(valid_leaf, level, val);

		//if(z*cell_size < 62 || z*cell_size > 64)
		//	continue;

		for (j = 0 ; j < 12 ; j++ ) {

			//if(z*cell_size == 62 && j < 4) continue;
			//if(z*cell_size == 64 && (j > 3 && j < 8)) continue;

			if (is_eflag_on(x,y,z,level,j)) continue;
			intersect_id = is_intersect(val, j);

			// intersection edge
			if (intersect_id == 1 || intersect_id == -1) {
				if (is_min_edge(valid_leaf , j , vtx , vtx_num , intersect_id , geofrm)) {
					eflag_on(x , y , z , level , j);
					geofrm.AddBound(vtx[0], 1);	geofrm.AddBound(vtx[1], 1);
					geofrm.AddBound(vtx[2], 1);	geofrm.AddBound(vtx[3], 1);

					get_min_vertex(j, intersect_id, x, y, z, tx, ty, tz);
					my_vertex = grid_idx_arr[(tz*dim[0]*dim[0]+ty*dim[0]+tx)*cell_size];
					if (my_vertex == -1) {
						add_one_vertex(tx, ty, tz, cell_size, my_vertex, geofrm);
						grid_idx_arr[(tz*dim[0]*dim[0]+ty*dim[0]+tx)*cell_size] = my_vertex;
					}
					//geofrm.AddBound(my_vertex, 1);

					//geofrm.AddPyramid(vtx, my_vertex, vtx_num);
					//geofrm.AddTetra(vtx[0], vtx[2], vtx[1], my_vertex);
					//geofrm.AddTetra(vtx[0], vtx[3], vtx[2], my_vertex);
					geofrm.Add_2_Tetra(vtx, my_vertex);
					//geofrm.Extend_Tetra(vtx);
				}
			}	
			// boundary cell edge whose two end points' values < isovalue 
			//else if((intersect_id == 2 || intersect_id == -2) && (! is_skipcell(valid_leaf))) {
			else if((intersect_id == 2 || intersect_id == -2) && minmax[valid_leaf].min <= iso_val) {
				if (is_min_edge_2(valid_leaf, j, my_vtx, vtx_num, con_id, intersect_id, geofrm)) {
					eflag_on(x, y, z, level, j);

					if( (my_vtx[0] == my_vtx[1] || my_vtx[0] == -1 || my_vtx[1] == -1) && 
						(my_vtx[1] == my_vtx[2] || my_vtx[1] == -1 || my_vtx[2] == -1) && 
						(my_vtx[2] == my_vtx[3] || my_vtx[2] == -1 || my_vtx[3] == -1) && 
						(my_vtx[3] == my_vtx[0] || my_vtx[3] == -1 || my_vtx[0] == -1) )	continue;

					get_min_vertex(j, 1, x, y, z, tx, ty, tz);
					my_vertex = grid_idx_arr[(tz*dim[0]*dim[0]+ty*dim[0]+tx)*cell_size];
					if (my_vertex == -1) {
						add_one_vertex(tx, ty, tz, cell_size, my_vertex, geofrm); // beginning point
						grid_idx_arr[(tz*dim[0]*dim[0]+ty*dim[0]+tx)*cell_size] = my_vertex;
					}

					get_min_vertex(j, -1, x, y, z, tx, ty, tz);
					my_vertex_1 = grid_idx_arr[(tz*dim[0]*dim[0]+ty*dim[0]+tx)*cell_size];
					if (my_vertex_1 == -1) {
						add_one_vertex(tx, ty, tz, cell_size, my_vertex_1, geofrm); // ending point
						grid_idx_arr[(tz*dim[0]*dim[0]+ty*dim[0]+tx)*cell_size] = my_vertex_1;
					}

					//geofrm.AddBound(my_vtx[0], 1);	geofrm.AddBound(my_vtx[1], 1);
					//geofrm.AddBound(my_vtx[2], 1);	geofrm.AddBound(my_vtx[3], 1);
					//geofrm.AddBound(my_vertex, 1);
					//geofrm.AddBound(my_vertex_1, 1);

					
					//if(my_vtx[0] != -1 && my_vtx[1] != -1 && con_id[0] == 1 && my_vtx[0] != my_vtx[1])
					if(my_vtx[0] != -1 && my_vtx[1] != -1 && my_vtx[0] != my_vtx[1])
						geofrm.AddTetra(my_vtx[0], my_vtx[1], my_vertex, my_vertex_1);
					//if(my_vtx[1] != -1 && my_vtx[2] != -1 && con_id[1] == 1 && my_vtx[1] != my_vtx[2])
					if(my_vtx[1] != -1 && my_vtx[2] != -1 && my_vtx[1] != my_vtx[2])
						geofrm.AddTetra(my_vtx[1], my_vtx[2], my_vertex, my_vertex_1);
					//if(my_vtx[2] != -1 && my_vtx[3] != -1 && con_id[2] == 1 && my_vtx[2] != my_vtx[3])
					if(my_vtx[2] != -1 && my_vtx[3] != -1 && my_vtx[2] != my_vtx[3])
						geofrm.AddTetra(my_vtx[2], my_vtx[3], my_vertex, my_vertex_1);
					//if(my_vtx[3] != -1 && my_vtx[0] != -1 && con_id[3] == 1 && my_vtx[3] != my_vtx[0])
					if(my_vtx[3] != -1 && my_vtx[0] != -1 && my_vtx[3] != my_vtx[0])
						geofrm.AddTetra(my_vtx[3], my_vtx[0], my_vertex, my_vertex_1);
				}
			}
			else
				continue;
		}

		// boundary cell -- val(four vertice of one face) < isovalue
		if(minmax[valid_leaf].max >= iso_val && minmax[valid_leaf].min <= iso_val) {
			//add_tetra_face(valid_leaf, level, geofrm);
		}

		// interior cell
		if(minmax[valid_leaf].max < iso_val) {
			//add_tetra_cube_adaptive(valid_leaf, level, geofrm);
		}

	}
}

void Octree::tetrahedralize_interval(geoframe& geofrm) {

	int x, y, z, tx, ty, tz, valid_leaf, cell_size, level, i, j;
	int vtx_num, intersect_id, my_vtx[4], con_id[4];
	unsigned int vtx[4], my_vertex, my_vertex_1;
	float val[8];

	int k;
	for (k = 0;k < octcell_num; k++) {vtx_idx_arr[k] = -1; vtx_idx_arr_in[k] = -1;}
	for (k = 0; k < dim[0]*dim[1]*dim[2]; k++) grid_idx_arr[k] = -1;

	for (i = 0; i < leaf_num; i++ ) {

		valid_leaf = cut_array[i] ;
		level = get_level(valid_leaf) ;
		cell_size = (dim[0]-1)/(1<<level);
		octcell2xyz(valid_leaf, x, y, z, level);
		getCellValues(valid_leaf, level, val);

		for (j = 0 ; j < 12 ; j++ ) {

			if (is_eflag_on(x,y,z,level,j)) continue;
			intersect_id = is_intersect_interval(val, j);

			if(is_skipcell(valid_leaf) == 0) in_out = 0;
			else in_out = 1;
			// intersection edge
			if (intersect_id == 1 || intersect_id == -1) {
				if (is_min_edge(valid_leaf, j, vtx, vtx_num, intersect_id, geofrm)) {
					eflag_on(x , y , z , level , j);

					if(! is_skipcell(valid_leaf)) {
						geofrm.AddBound(vtx[0], 1);	geofrm.AddBound(vtx[1], 1);
						geofrm.AddBound(vtx[2], 1);	geofrm.AddBound(vtx[3], 1);

						get_min_vertex(j, intersect_id, x, y, z, tx, ty, tz);
						add_one_vertex(tx, ty, tz, cell_size, my_vertex, geofrm);

						//geofrm.AddPyramid(vtx, my_vertex, vtx_num);
						//geofrm.AddTetra(vtx[0], vtx[2], vtx[1], my_vertex);
						//geofrm.AddTetra(vtx[0], vtx[3], vtx[2], my_vertex);
						geofrm.Add_2_Tetra(vtx, my_vertex);
					}
					else {
						geofrm.AddBound(vtx[0], -1);	geofrm.AddBound(vtx[1], -1);
						geofrm.AddBound(vtx[2], -1);	geofrm.AddBound(vtx[3], -1);

						get_min_vertex(j, -intersect_id, x, y, z, tx, ty, tz);
						add_one_vertex(tx, ty, tz, cell_size, my_vertex, geofrm);

						//geofrm.AddPyramid(vtx, my_vertex, vtx_num);
						//geofrm.AddTetra(vtx[0], vtx[1], vtx[2], my_vertex);
						//geofrm.AddTetra(vtx[0], vtx[2], vtx[3], my_vertex);
						geofrm.Add_2_Tetra(vtx, my_vertex);
					}
				}
			}	
			// boundary cell edge whose two end points' values < isovalue 
			//if((intersect_id == 2 || intersect_id == -2) && (! is_skipcell(valid_leaf))) {
			else if(intersect_id == 2 || intersect_id == -2) {
				if (is_min_edge_2(valid_leaf, j, my_vtx, vtx_num, con_id, intersect_id, geofrm)) {
					eflag_on(x, y, z, level, j);

					if( (my_vtx[0] == my_vtx[1] || my_vtx[0] == -1 || my_vtx[1] == -1) && 
						(my_vtx[1] == my_vtx[2] || my_vtx[1] == -1 || my_vtx[2] == -1) && 
						(my_vtx[2] == my_vtx[3] || my_vtx[2] == -1 || my_vtx[3] == -1) && 
						(my_vtx[3] == my_vtx[0] || my_vtx[3] == -1 || my_vtx[0] == -1) )	continue;

					get_min_vertex(j, 1, x, y, z, tx, ty, tz);
					my_vertex = grid_idx_arr[(tz*dim[0]*dim[0]+ty*dim[0]+tx)*cell_size];
					if (my_vertex == -1) {
						add_one_vertex(tx, ty, tz, cell_size, my_vertex, geofrm); // beginning point
						grid_idx_arr[(tz*dim[0]*dim[0]+ty*dim[0]+tx)*cell_size] = my_vertex;
					}
					//get_min_vertex(j, 1, x, y, z, tx, ty, tz);
					//add_one_vertex(tx, ty, tz, cell_size, my_vertex, geofrm); // beginning point

					get_min_vertex(j, -1, x, y, z, tx, ty, tz);
					my_vertex_1 = grid_idx_arr[(tz*dim[0]*dim[0]+ty*dim[0]+tx)*cell_size];
					if (my_vertex_1 == -1) {
						add_one_vertex(tx, ty, tz, cell_size, my_vertex_1, geofrm); // ending point
						grid_idx_arr[(tz*dim[0]*dim[0]+ty*dim[0]+tx)*cell_size] = my_vertex_1;
					}
					//get_min_vertex(j, -1, x, y, z, tx, ty, tz);
					//add_one_vertex(tx, ty, tz, cell_size, my_vertex_1, geofrm); // ending point

					//if(my_vtx[0] != -1 && my_vtx[1] != -1 && con_id[0] == 1)
					if(my_vtx[0] != -1 && my_vtx[1] != -1 && my_vtx[0] != my_vtx[1])
						geofrm.AddTetra(my_vtx[0], my_vtx[1], my_vertex, my_vertex_1);
					//if(my_vtx[1] != -1 && my_vtx[2] != -1 && con_id[1] == 1)
					if(my_vtx[1] != -1 && my_vtx[2] != -1 && my_vtx[1] != my_vtx[2])
						geofrm.AddTetra(my_vtx[1], my_vtx[2], my_vertex, my_vertex_1);
					//if(my_vtx[2] != -1 && my_vtx[3] != -1 && con_id[2] == 1)
					if(my_vtx[2] != -1 && my_vtx[3] != -1 && my_vtx[2] != my_vtx[3])
						geofrm.AddTetra(my_vtx[2], my_vtx[3], my_vertex, my_vertex_1);
					//if(my_vtx[3] != -1 && my_vtx[0] != -1 && con_id[3] == 1)
					if(my_vtx[3] != -1 && my_vtx[0] != -1 && my_vtx[3] != my_vtx[0])
						geofrm.AddTetra(my_vtx[3], my_vtx[0], my_vertex, my_vertex_1);
				}
			}
			else
				continue;
		}

		// boundary cell -- val(four vertice of one face) < isovalue
		if( (minmax[valid_leaf].max >= iso_val && minmax[valid_leaf].min <= iso_val) ||
			(minmax[valid_leaf].max >= iso_val_in && minmax[valid_leaf].min <= iso_val_in)) {
			//add_tetra_face_interval(valid_leaf, level, geofrm);
		}

		// interior cell
		if(minmax[valid_leaf].max < iso_val && minmax[valid_leaf].min > iso_val_in) {
			//add_tetra_cube_adaptive(valid_leaf, level, geofrm);
		}

	}
}

int Octree::is_min_edge_2(int oc_id, int e_id, int* vtx, int& vtx_num, int* con_id, int intersect_id, geoframe& geofrm)
{
	int x, y, z, level, i;
	int temp_vtx[4], temp_id[8];
	
	level = get_level(oc_id);
	octcell2xyz(oc_id, x, y, z, level);

	vtx_num = 4;
	
	//if(!(minmax[oc_id].max >= iso_val && minmax[oc_id].min <= iso_val)) {
	//	temp_vtx[0] = -1; temp_vtx[1] = -1; temp_vtx[2] = -1; temp_vtx[3] = -1; 
	//	return 0;
	//}
	vtx[0] = -1;		vtx[1] = -1;		vtx[2] = -1;		vtx[3] = -1; 
	temp_vtx[0] = -1;	temp_vtx[1] = -1;	temp_vtx[2] = -1;	temp_vtx[3] = -1; 
	temp_id[0] = 1;		temp_id[1] = 1;		temp_id[2] = 1;		temp_id[3] = 1;

	assert(! is_refined(x, y, z, level));
	temp_vtx[0] = min_vtx_tetra(x, y, z, e_id, e_id, level, temp_id[0], temp_id[4], geofrm);

	switch (e_id) {
		case 0 : 
			if (is_refined(x,y,z-1,level) || is_refined(x,y-1,z-1,level) || is_refined(x,y-1,z,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x,  y,    z-1,  0,  2,  level,  temp_id[1],  temp_id[5],  geofrm);
			temp_vtx[2] = min_vtx_tetra(x,  y-1,  z-1,  0,  6,  level,  temp_id[2],  temp_id[6],  geofrm);
			temp_vtx[3] = min_vtx_tetra(x,  y-1,  z,    0,  4,  level,  temp_id[3],  temp_id[7],  geofrm);
			break;

		case 1 : 
			if (is_refined(x,y-1,z,level) || is_refined(x+1,y-1,z,level) || is_refined(x+1,y,z,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x+1,  y,    z,  1,  3,  level,  temp_id[1],  temp_id[5], geofrm);
			temp_vtx[2] = min_vtx_tetra(x+1,  y-1,  z,  1,  7,  level,  temp_id[2],  temp_id[6], geofrm);
			temp_vtx[3] = min_vtx_tetra(x,    y-1,  z,  1,  5,  level,  temp_id[3],  temp_id[7], geofrm);
			break;

		case 2 : 
			if (is_refined(x,y,z+1,level) || is_refined(x,y-1,z+1,level) || is_refined(x,y-1,z,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x,  y,    z+1,  2,  0,  level,  temp_id[1],  temp_id[5],  geofrm);
			temp_vtx[2] = min_vtx_tetra(x,  y-1,  z+1,  2,  4,  level,  temp_id[2],  temp_id[6],  geofrm);
			temp_vtx[3] = min_vtx_tetra(x,  y-1,  z,    2,  6,  level,  temp_id[3],  temp_id[7],  geofrm);
			break;

		case 3 : 
			if (is_refined(x,y-1,z,level) || is_refined(x-1,y-1,z,level) || is_refined(x-1,y,z,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x,    y-1,  z,  3,  7, level,  temp_id[1],  temp_id[5],  geofrm);
			temp_vtx[2] = min_vtx_tetra(x-1,  y-1,  z,  3,  5, level,  temp_id[2],  temp_id[6],  geofrm);
			temp_vtx[3] = min_vtx_tetra(x-1,  y,    z,  3,  1, level,  temp_id[3],  temp_id[7],  geofrm);
			break;

		case 4 : 
			if (is_refined(x,y,z-1,level) || is_refined(x,y+1,z-1,level) || is_refined(x,y+1,z,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x,  y+1,  z,    4,  0,  level,  temp_id[1],  temp_id[5],  geofrm);
			temp_vtx[2] = min_vtx_tetra(x,  y+1,  z-1,  4,  2,  level,  temp_id[2],  temp_id[6],  geofrm);
			temp_vtx[3] = min_vtx_tetra(x,  y,    z-1,  4,  6,  level,  temp_id[3],  temp_id[7],  geofrm);
			break;

		case 5 :
			if (is_refined(x,y+1,z,level) || is_refined(x+1,y,z,level) || is_refined(x+1,y+1,z,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x,    y+1,  z,  5,  1, level,  temp_id[1],  temp_id[5],  geofrm);
			temp_vtx[2] = min_vtx_tetra(x+1,  y+1,  z,  5,  3, level,  temp_id[2],  temp_id[6],  geofrm);
			temp_vtx[3] = min_vtx_tetra(x+1,  y,    z,  5,  7, level,  temp_id[3],  temp_id[7],  geofrm);
			break;

		case 6 :
			if (is_refined(x,y+1,z,level) || is_refined(x,y+1,z+1,level) || is_refined(x,y,z+1,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x,  y+1,  z,    6,  2,  level,  temp_id[1],  temp_id[5],  geofrm);
			temp_vtx[2] = min_vtx_tetra(x,  y+1,  z+1,  6,  0,  level,  temp_id[2],  temp_id[6],  geofrm);
			temp_vtx[3] = min_vtx_tetra(x,  y,    z+1,  6,  4,  level,  temp_id[3],  temp_id[7],  geofrm);
			break;

		case 7 :
			if (is_refined(x-1,y,z,level) || is_refined(x-1,y+1,z,level) || is_refined(x,y+1,z,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x-1,  y,    z,  7,  5,  level,  temp_id[1],  temp_id[5],  geofrm);
			temp_vtx[2] = min_vtx_tetra(x-1,  y+1,  z,  7,  1,  level,  temp_id[2],  temp_id[6],  geofrm);
			temp_vtx[3] = min_vtx_tetra(x,    y+1,  z,  7,  3,  level,  temp_id[3],  temp_id[7],  geofrm);
			break;

		case 8 :
			if (is_refined(x,y,z-1,level) || is_refined(x-1,y,z-1,level) || is_refined(x-1,y,z,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x-1,  y,  z,    8,  9,   level,  temp_id[1],  temp_id[5],  geofrm);
			temp_vtx[2] = min_vtx_tetra(x-1,  y,  z-1,  8,  11,  level,  temp_id[2],  temp_id[6],  geofrm);
			temp_vtx[3] = min_vtx_tetra(x,    y,  z-1,  8,  10,  level,  temp_id[3],  temp_id[7],  geofrm);
			break;

		case 9 :
			if (is_refined(x,y,z-1,level) || is_refined(x+1,y,z-1,level) || is_refined(x+1,y,z,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x,    y,  z-1,  9,  11,  level,  temp_id[1],  temp_id[5],  geofrm);
			temp_vtx[2] = min_vtx_tetra(x+1,  y,  z-1,  9,  10,  level,  temp_id[2],  temp_id[6],  geofrm);
			temp_vtx[3] = min_vtx_tetra(x+1,  y,  z,    9,  8,   level,  temp_id[3],  temp_id[7],  geofrm);
			break;

		case 10 :
			if (is_refined(x,y,z+1,level) || is_refined(x-1,y,z+1,level) || is_refined(x-1,y,z,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x,    y,  z+1,  10,  8,   level,  temp_id[1],  temp_id[5],  geofrm);
			temp_vtx[2] = min_vtx_tetra(x-1,  y,  z+1,  10,  9,   level,  temp_id[2],  temp_id[6],  geofrm);
			temp_vtx[3] = min_vtx_tetra(x-1,  y,  z,    10,  11,  level,  temp_id[3],  temp_id[7],  geofrm);
			break;

		case 11 :
			if (is_refined(x,y,z+1,level) || is_refined(x+1,y,z+1,level) || is_refined(x+1,y,z,level)) return 0;
			temp_vtx[1] = min_vtx_tetra(x+1,  y,  z,    11,  10,  level,  temp_id[1],  temp_id[5],  geofrm);
			temp_vtx[2] = min_vtx_tetra(x+1,  y,  z+1,  11,  8,   level,  temp_id[2],  temp_id[6],  geofrm);
			temp_vtx[3] = min_vtx_tetra(x,    y,  z+1,  11,  9,   level,  temp_id[3],  temp_id[7],  geofrm);
			break;
	}

	assert(intersect_id == 2 || intersect_id == -2);

	temp_id[0] = temp_id[0] && temp_id[5];
	temp_id[1] = temp_id[1] && temp_id[6];
	temp_id[2] = temp_id[2] && temp_id[7];
	temp_id[3] = temp_id[3] && temp_id[4];

	if (intersect_id == 2) 
		for (i=0;i<4;i++) {
			vtx[i] = temp_vtx[i];	con_id[i] = temp_id[i];
		}
	else if (intersect_id == -2) {
		for (i=0;i<4;i++) {
			vtx[i] = temp_vtx[3-i];	//con_id[i] = temp_id[3-i];
		}
		con_id[0] = temp_id[2];		con_id[1] = temp_id[1];
		con_id[2] = temp_id[0];		con_id[3] = temp_id[3];
	}

	return 1;

}

void Octree::get_min_vertex(int e_id, int intersect_id, int x, int y, int z, int& tx, int& ty, int& tz)
{
	if(intersect_id == 1) {
		if(e_id == 0) {
			tx = x; ty = y; tz = z;	}
		else if(e_id == 1) {
			tx = x+1; ty = y; tz = z;	}
		else if(e_id == 2) {
			tx = x+1; ty = y; tz = z+1;	}
		else if(e_id == 3) {
			tx = x; ty = y; tz = z;	}
		else if(e_id == 4) {
			tx = x; ty = y+1; tz = z;	}
		else if(e_id == 5) {
			tx = x+1; ty = y+1; tz = z;	}
		else if(e_id == 6) {
			tx = x+1; ty = y+1; tz = z+1;	}
		else if(e_id == 7) {
			tx = x; ty = y+1; tz = z;	}
		else if(e_id == 8) {
			tx = x; ty = y; tz = z;	}
		else if(e_id == 9) {
			tx = x+1; ty = y; tz = z;	}
		else if(e_id == 10) {
			tx = x; ty = y; tz = z+1;	}
		else {	//	if(e_id == 11) 
			tx = x+1; ty = y; tz = z+1;	}
	}
	else {	// intersect_id == -1
		if(e_id == 0) {
			tx = x+1; ty = y; tz = z;	}
		else if(e_id == 1) {
			tx = x+1; ty = y; tz = z+1;	}
		else if(e_id == 2) {
			tx = x; ty = y; tz = z+1;	}
		else if(e_id == 3) {
			tx = x; ty = y; tz = z+1;	}
		else if(e_id == 4) {
			tx = x+1; ty = y+1; tz = z;	}
		else if(e_id == 5) {
			tx = x+1; ty = y+1; tz = z+1;	}
		else if(e_id == 6) {
			tx = x; ty = y+1; tz = z+1;	}
		else if(e_id == 7) {
			tx = x; ty = y+1; tz = z+1;	}
		else if(e_id == 8) {
			tx = x; ty = y+1; tz = z;	}
		else if(e_id == 9) {
			tx = x+1; ty = y+1; tz = z;	}
		else if(e_id == 10) {
			tx = x; ty = y+1; tz = z+1;	}
		else {	//	if(e_id == 11) 
			tx = x+1; ty = y+1; tz = z+1;	}
	}

}

int Octree::min_vtx_tetra(int x, int y, int z, int e_id_0, int e_id, int level, int& con_ind, int& con_ind_1, geoframe& geofrm)
{
	int tx,ty,tz, vert;
	float vtx[3], norm[3], val[8];
	unsigned int center;

	tx = x; ty = y; tz = z;
	assert( tx>=0 && ty>=0 && tz>=0 );
	assert( !is_refined(tx,ty,tz,level) );

	while ( level==0 || !is_refined(tx/2 , ty/2 , tz/2 , level-1 )) {
		tx/=2;
		ty/=2;
		tz/=2;
		level--;
	}

	int oc_id = xyz2octcell(tx, ty, tz, level);
	int cell_size = (dim[0]-1)/(1<<level);
	getCellValues(oc_id, level, val);
/*
	int bool_0 = (val[0] <= iso_val && val[3] <= iso_val && val[4] <= iso_val && val[7] <= iso_val);
	int bool_1 = (val[1] <= iso_val && val[2] <= iso_val && val[5] <= iso_val && val[6] <= iso_val);
	int bool_2 = (val[0] <= iso_val && val[1] <= iso_val && val[2] <= iso_val && val[3] <= iso_val);
	int bool_3 = (val[4] <= iso_val && val[5] <= iso_val && val[6] <= iso_val && val[7] <= iso_val);
	int bool_4 = (val[0] <= iso_val && val[1] <= iso_val && val[4] <= iso_val && val[5] <= iso_val);
	int bool_5 = (val[2] <= iso_val && val[3] <= iso_val && val[6] <= iso_val && val[7] <= iso_val);

	if(flag_type > 3) {
		int bool_00 = (val[0] >= iso_val_in && val[3] >= iso_val_in && val[4] >= iso_val_in && val[7] >= iso_val_in);
		int bool_01 = (val[1] >= iso_val_in && val[2] >= iso_val_in && val[5] >= iso_val_in && val[6] >= iso_val_in);
		int bool_02 = (val[0] >= iso_val_in && val[1] >= iso_val_in && val[2] >= iso_val_in && val[3] >= iso_val_in);
		int bool_03 = (val[4] >= iso_val_in && val[5] >= iso_val_in && val[6] >= iso_val_in && val[7] >= iso_val_in);
		int bool_04 = (val[0] >= iso_val_in && val[1] >= iso_val_in && val[4] >= iso_val_in && val[5] >= iso_val_in);
		int bool_05 = (val[2] >= iso_val_in && val[3] >= iso_val_in && val[6] >= iso_val_in && val[7] >= iso_val_in);

		bool_0 = bool_0 && bool_00;		bool_1 = bool_1 && bool_01;
		bool_2 = bool_2 && bool_02;		bool_3 = bool_3 && bool_03;
		bool_4 = bool_4 && bool_04;		bool_5 = bool_5 && bool_05;
	}
	con_ind = 1;	con_ind_1 = 1;

	switch (e_id) {
	case 0:
		if(bool_2 && bool_4) return -1;
		else {
			//if((e_id_0 == 0 || e_id_0 == 4) && bool_4) con_ind = 0;
			//if((e_id_0 == 2 || e_id_0 == 6) && bool_2) con_ind = 0;
			if(e_id_0 == 0 || e_id_0 == 4) {
				if(bool_4) con_ind = 0;
				if(bool_2) con_ind_1 = 0;
			}
			if(e_id_0 == 2 || e_id_0 == 6) {
				if(bool_2) con_ind = 0;
				if(bool_4) con_ind_1 = 0;
			}
		}
		break;
	case 1:
		if(bool_1 && bool_2) return -1;
		else {
			if(bool_1) con_ind = 0;
			if(bool_2) con_ind_1 = 0;
		}
		break;
	case 2:
		if(bool_2 && bool_5) return -1;
		else {
			//if((e_id_0 == 0 || e_id_0 == 4) && bool_2) con_ind = 0;
			//if((e_id_0 == 2 || e_id_0 == 6) && bool_5) con_ind = 0;
			if(e_id_0 == 0 || e_id_0 == 4) {
				if(bool_2) con_ind = 0;
				if(bool_5) con_ind_1 = 0;
			}
			if(e_id_0 == 2 || e_id_0 == 6) {
				if(bool_5) con_ind = 0;
				if(bool_2) con_ind_1 = 0;
			}
		}
		break;
	case 3:
		if(bool_0 && bool_2) return -1;
		else {
			if(bool_2) con_ind = 0;
			if(bool_0) con_ind_1 = 0;
		}
		break;
	case 4:
		if(bool_3 && bool_4) return -1;
		else {
			//if((e_id_0 == 0 || e_id_0 == 4) && bool_3) con_ind = 0;
			//if((e_id_0 == 2 || e_id_0 == 6) && bool_4) con_ind = 0;
			if(e_id_0 == 0 || e_id_0 == 4) {
				if(bool_3) con_ind = 0;
				if(bool_4) con_ind_1 = 0;
			}
			if(e_id_0 == 2 || e_id_0 == 6) {
				if(bool_4) con_ind = 0;
				if(bool_3) con_ind_1 = 0;
			}
		}
		break;
	case 5:
		if(bool_1 && bool_3) return -1;
		else {
			//if(bool_1 || bool_3) con_ind = 0;
			if(bool_3) con_ind = 0;
			if(bool_1) con_ind_1 = 0;
		}
		break;
	case 6:
		if(bool_3 && bool_5) return -1;
		else {
			//if(bool_5 || bool_3) con_ind = 0;
			//if((e_id_0 == 0 || e_id_0 == 4) && bool_5) con_ind = 0;
			//if((e_id_0 == 2 || e_id_0 == 6) && bool_3) con_ind = 0;
			if(e_id_0 == 0 || e_id_0 == 4) {
				if(bool_5) con_ind = 0;
				if(bool_3) con_ind_1 = 0;
			}
			if(e_id_0 == 2 || e_id_0 == 6) {
				if(bool_3) con_ind = 0;
				if(bool_5) con_ind_1 = 0;
			}
		}
		break;
	case 7:
		if(bool_0 && bool_3) return -1;
		else {
			if(bool_0) con_ind = 0;
			if(bool_3) con_ind_1 = 0;
		}
		break;
	case 8:
		if(bool_0 && bool_4) return -1;
		else {
			if(bool_0) con_ind = 0;
			if(bool_4) con_ind_1 = 0;
		}
		break;
	case 9:
		if(bool_1 && bool_4) return -1;
		else {
			if(bool_4) con_ind = 0;
			if(bool_1) con_ind_1 = 0;
		}
		break;
	case 10:
		if(bool_0 && bool_5) return -1;
		else {
			if(bool_5) con_ind = 0;
			if(bool_0) con_ind_1 = 0;
		}
		break;
	case 11:
		if(bool_1 && bool_5) return -1;
		else {
			if(bool_1) con_ind = 0;
			if(bool_5) con_ind_1 = 0;
		}
		break;
	}
*/
	//if(minmax[oc_id].max <= iso_val) {
	if( (minmax[oc_id].max <= iso_val && flag_type < 4) ||
		(minmax[oc_id].max <= iso_val && minmax[oc_id].min >= iso_val_in && flag_type >= 4) ) {
		if ((center = vtx_idx_arr[xyz2octcell(tx, ty, tz, level)]) == -1) {
			add_middle_vertex(tx, ty, tz, 0.5, 0.5, 0.5, cell_size, center, geofrm);
			vtx_idx_arr[xyz2octcell(tx, ty, tz, level)] = center;
			return center;
		}
		else
			return center;
	}
	else {
		get_vtx(tx, ty, tz, level, vtx);
		getVertGrad(tx*cell_size, ty*cell_size, tz*cell_size, norm);
		if ((vert = vtx_idx_arr[xyz2octcell(tx, ty, tz, level)]) == -1) {
			vert = geofrm.AddVert(vtx, norm);
			vtx_idx_arr[xyz2octcell(tx, ty, tz, level)] = vert;
			return vert;
		}
		else
			return vert;
	}
}

void Octree::add_tetra_face(int oc_id, int level, geoframe& geofrm) {

	int x, y, z, cell_size;
	unsigned int my_vtx[8];
	float val[8];

	cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_id, x, y, z, level);
	getCellValues(oc_id, level, val);

	int	my_vertex = min_vtx(x,y,z,level,geofrm);

	if(val[0] <= iso_val && val[3] <= iso_val && val[4] <= iso_val && val[7] <= iso_val) {
		add_one_vertex(x,   y,   z,   cell_size, my_vtx[0], geofrm);
		add_one_vertex(x,   y,   z+1, cell_size, my_vtx[3], geofrm);
		add_one_vertex(x,   y+1, z+1, cell_size, my_vtx[7], geofrm);
		add_one_vertex(x,   y+1, z,   cell_size, my_vtx[4], geofrm);
		march_one_face(0, oc_id, level, my_vtx, my_vertex, geofrm); }
	if(val[1] <= iso_val && val[2] <= iso_val && val[5] <= iso_val && val[6] <= iso_val) {
		add_one_vertex(x+1, y,   z,   cell_size, my_vtx[1], geofrm);
		add_one_vertex(x+1, y,   z+1, cell_size, my_vtx[2], geofrm);
		add_one_vertex(x+1, y+1, z+1, cell_size, my_vtx[6], geofrm);
		add_one_vertex(x+1, y+1, z,   cell_size, my_vtx[5], geofrm);
		march_one_face(1, oc_id, level, my_vtx, my_vertex, geofrm); }
	if(val[0] <= iso_val && val[1] <= iso_val && val[2] <= iso_val && val[3] <= iso_val) {
		add_one_vertex(x,   y,   z,   cell_size, my_vtx[0], geofrm);
		add_one_vertex(x+1, y,   z,   cell_size, my_vtx[1], geofrm);
		add_one_vertex(x+1, y,   z+1, cell_size, my_vtx[2], geofrm);
		add_one_vertex(x,   y,   z+1, cell_size, my_vtx[3], geofrm);
		march_one_face(2, oc_id, level, my_vtx, my_vertex, geofrm); }
	if(val[4] <= iso_val && val[5] <= iso_val && val[6] <= iso_val && val[7] <= iso_val) {
		add_one_vertex(x,   y+1, z,   cell_size, my_vtx[4], geofrm);
		add_one_vertex(x+1, y+1, z,   cell_size, my_vtx[5], geofrm);
		add_one_vertex(x+1, y+1, z+1, cell_size, my_vtx[6], geofrm);
		add_one_vertex(x,   y+1, z+1, cell_size, my_vtx[7], geofrm);
		march_one_face(3, oc_id, level, my_vtx, my_vertex, geofrm); }
	if(val[0] <= iso_val && val[1] <= iso_val && val[4] <= iso_val && val[5] <= iso_val) {
		add_one_vertex(x,   y,   z,   cell_size, my_vtx[0], geofrm);
		add_one_vertex(x+1, y,   z,   cell_size, my_vtx[1], geofrm);
		add_one_vertex(x,   y+1, z,   cell_size, my_vtx[4], geofrm);
		add_one_vertex(x+1, y+1, z,   cell_size, my_vtx[5], geofrm);
		march_one_face(4, oc_id, level, my_vtx, my_vertex, geofrm); }
	if(val[2] <= iso_val && val[3] <= iso_val && val[6] <= iso_val && val[7] <= iso_val) {
		add_one_vertex(x+1, y,   z+1, cell_size, my_vtx[2], geofrm);
		add_one_vertex(x,   y,   z+1, cell_size, my_vtx[3], geofrm);
		add_one_vertex(x+1, y+1, z+1, cell_size, my_vtx[6], geofrm);
		add_one_vertex(x,   y+1, z+1, cell_size, my_vtx[7], geofrm);
		march_one_face(5, oc_id, level, my_vtx, my_vertex, geofrm); }

}

void Octree::add_tetra_face_interval(int oc_id, int level, geoframe& geofrm) {

	int x, y, z, cell_size;
	unsigned int my_vtx[8];
	float val[8];

	cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_id, x, y, z, level);
	getCellValues(oc_id, level, val);

	int	my_vertex = min_vtx(x,y,z,level,geofrm);

	add_one_vertex(x,   y,   z,   cell_size, my_vtx[0], geofrm);
	add_one_vertex(x,   y,   z+1, cell_size, my_vtx[3], geofrm);
	add_one_vertex(x,   y+1, z+1, cell_size, my_vtx[7], geofrm);
	add_one_vertex(x,   y+1, z,   cell_size, my_vtx[4], geofrm);

	add_one_vertex(x+1, y,   z,   cell_size, my_vtx[1], geofrm);
	add_one_vertex(x+1, y,   z+1, cell_size, my_vtx[2], geofrm);
	add_one_vertex(x+1, y+1, z+1, cell_size, my_vtx[6], geofrm);
	add_one_vertex(x+1, y+1, z,   cell_size, my_vtx[5], geofrm);

	if( val[0] <= iso_val && val[3] <= iso_val && val[4] <= iso_val && val[7] <= iso_val &&
		val[0] >= iso_val_in && val[3] >= iso_val_in && val[4] >= iso_val_in && val[7] >= iso_val_in ) {
		march_one_face(0, oc_id, level, my_vtx, my_vertex, geofrm); }
	if( val[1] <= iso_val && val[2] <= iso_val && val[5] <= iso_val && val[6] <= iso_val &&
		val[1] >= iso_val_in && val[2] >= iso_val_in && val[5] >= iso_val_in && val[6] >= iso_val_in) {
		march_one_face(1, oc_id, level, my_vtx, my_vertex, geofrm); }
	if( val[0] <= iso_val && val[1] <= iso_val && val[2] <= iso_val && val[3] <= iso_val &&
		val[0] >= iso_val_in && val[1] >= iso_val_in && val[2] >= iso_val_in && val[3] >= iso_val_in) {
		march_one_face(2, oc_id, level, my_vtx, my_vertex, geofrm); }
	if( val[4] <= iso_val && val[5] <= iso_val && val[6] <= iso_val && val[7] <= iso_val &&
		val[4] >= iso_val_in && val[5] >= iso_val_in && val[6] >= iso_val_in && val[7] >= iso_val_in) {
		march_one_face(3, oc_id, level, my_vtx, my_vertex, geofrm); }
	if( val[0] <= iso_val && val[1] <= iso_val && val[4] <= iso_val && val[5] <= iso_val &&
		val[0] >= iso_val_in && val[1] >= iso_val_in && val[4] >= iso_val_in && val[5] >= iso_val_in) {
		march_one_face(4, oc_id, level, my_vtx, my_vertex, geofrm); }
	if( val[2] <= iso_val && val[3] <= iso_val && val[6] <= iso_val && val[7] <= iso_val &&
		val[2] >= iso_val_in && val[3] >= iso_val_in && val[6] >= iso_val_in && val[7] >= iso_val_in) {
		march_one_face(5, oc_id, level, my_vtx, my_vertex, geofrm); }

}

void Octree::add_one_vertex(int x, int y, int z, int cell_size, unsigned int& my_vertex, geoframe& geofrm) {

	float pv[3], g[3];
	int tx, ty, tz;

	tx = x*cell_size;		ty = y*cell_size;		tz = z*cell_size;
	pv[0] = tx;				pv[1] = ty;				pv[2] = tz;
	getVertGrad(tx, ty, tz, g);
	my_vertex = geofrm.AddVert(pv, g);

}

void Octree::add_middle_vertex(int x, int y, int z, float dx, float dy, float dz, 
							   int cell_size, unsigned int& my_vertex, geoframe& geofrm) {

	float pv[3], g[3], g1[3], g2[3];
	int tx, ty, tz;

	pv[0] = (x + dx) * cell_size;
	pv[1] = (y + dy) * cell_size;
	pv[2] = (z + dz) * cell_size;
	
	tx = x*cell_size;	ty = y*cell_size;	tz = z*cell_size;
	getVertGrad(tx, ty, tz, g1);
	tx = (x+1)*cell_size;	ty = (y+1)*cell_size;	tz = (z+1)*cell_size;
	getVertGrad(tx, ty, tz, g2);
	g[0] = g1[0] + dx*(g2[0] - g1[0]);
	g[1] = g1[1] + dy*(g2[1] - g1[1]);
	g[2] = g1[2] + dz*(g2[2] - g1[2]);

	my_vertex = geofrm.AddVert(pv, g);

}
/*
void Octree::add_middle_vertex(int x, int y, int z, float dx, float dy, float dz, 
							   int cell_size, unsigned int& my_vertex, geoframe& geofrm) {

	float pv[3], g[3];
	int tx, ty, tz;

	tx = x*cell_size + dx * cell_size;
	ty = y*cell_size + dy * cell_size;
	tz = z*cell_size + dz * cell_size;
	pv[0] = tx;				pv[1] = ty;				pv[2] = tz;
	getVertGrad(tx, ty, tz, g);
	my_vertex = geofrm.AddVert(pv, g);

}
*/
void Octree::add_tetra_cube(int oc_id, int level, geoframe& geofrm) {

	int x, y, z;
	unsigned int my_vtx[8];

	int cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_id, x, y, z, level);

	add_one_vertex(x,   y,   z,   cell_size, my_vtx[0], geofrm);
	add_one_vertex(x+1, y,   z,   cell_size, my_vtx[1], geofrm);
	add_one_vertex(x+1, y,   z+1, cell_size, my_vtx[2], geofrm);
	add_one_vertex(x,   y,   z+1, cell_size, my_vtx[3], geofrm);
	add_one_vertex(x,   y+1, z,   cell_size, my_vtx[4], geofrm);
	add_one_vertex(x+1, y+1, z,   cell_size, my_vtx[5], geofrm);
	add_one_vertex(x+1, y+1, z+1, cell_size, my_vtx[6], geofrm);
	add_one_vertex(x,   y+1, z+1, cell_size, my_vtx[7], geofrm);

	if((x + y + z) % 2 == 0) {
		geofrm.AddTetra(my_vtx[0], my_vtx[1], my_vtx[3], my_vtx[4]);
		geofrm.AddTetra(my_vtx[1], my_vtx[5], my_vtx[6], my_vtx[4]);
		geofrm.AddTetra(my_vtx[3], my_vtx[2], my_vtx[6], my_vtx[1]);
		geofrm.AddTetra(my_vtx[3], my_vtx[6], my_vtx[7], my_vtx[4]);
		geofrm.AddTetra(my_vtx[1], my_vtx[3], my_vtx[4], my_vtx[6]);
	}
	else {
		geofrm.AddTetra(my_vtx[3], my_vtx[2], my_vtx[7], my_vtx[0]);
		geofrm.AddTetra(my_vtx[2], my_vtx[6], my_vtx[7], my_vtx[5]);
		geofrm.AddTetra(my_vtx[0], my_vtx[2], my_vtx[5], my_vtx[1]);
		geofrm.AddTetra(my_vtx[0], my_vtx[5], my_vtx[7], my_vtx[4]);
		geofrm.AddTetra(my_vtx[0], my_vtx[2], my_vtx[7], my_vtx[5]);
	}

}

void Octree::add_tetra_cube_adaptive(int oc_id, int level, geoframe& geofrm) {

	int x, y, z, cell_size;
	unsigned int center;

	int reg_bitmask = get_neighbor_bit(oc_id, level);
	if(reg_bitmask == 0) {
		add_tetra_cube(oc_id, level, geofrm);
	}
	else {
		cell_size = (dim[0]-1)/(1<<level);
		octcell2xyz(oc_id, x, y, z, level);

		// add center point as Steiner point
		add_middle_vertex(x, y, z, 0.5, 0.5, 0.5, cell_size, center, geofrm);

		for(int index = 0; index < 6; index++) {
			march_each_face(oc_id, level, index, center, geofrm);
		}
	}
}


void Octree::march_each_edge(int oc_id, int level, int edge_id, int* my_edge, geoframe& geofrm) 
{

	
	int x, y, z, tx, ty, tz, index, index0;

	int cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_id, x, y, z, level);

	for(index = 0; index < 128; index++)	my_edge[index] = -1;

	int diff_level = oct_depth - level;
	int t = 0;
	int num_edge = 1;
	tx = x;		ty = y;		tz = z;

	switch (abs(edge_id)) {
	case 0:
		for(index = level; index <= oct_depth; index++) {
			for(index0 = 0; index0 < num_edge; index0++) {
				if(is_refined2(tx+index0, ty-1, tz-1, index) 
				|| is_refined2(tx+index0, ty,   tz-1, index)
				|| is_refined2(tx+index0, ty,   tz,   index)
				|| is_refined2(tx+index0, ty-1, tz,   index))	my_edge[t] = 1;
				else	my_edge[t] = 0;
				t++;
			}
			tx *= 2;	ty *= 2;	tz *= 2;	num_edge *= 2;
		}
		break;

	case 100:
		for(index = level; index <= oct_depth; index++) {
			for(index0 = num_edge - 1; index0 >=0; index0--) {
				if(is_refined2(tx+index0, ty-1, tz-1, index) 
				|| is_refined2(tx+index0, ty,   tz-1, index)
				|| is_refined2(tx+index0, ty,   tz,   index)
				|| is_refined2(tx+index0, ty-1, tz,   index))	my_edge[t] = 1;
				else	my_edge[t] = 0;
				t++;
			}
			tx *= 2;	ty *= 2;	tz *= 2;	num_edge *= 2;
		}
		break;

	case 1:
		for(index = level; index <= oct_depth; index++) {
			if(edge_id == 1) {
				for(index0 = 0; index0 < num_edge; index0++) {
					if(is_refined2(tx,   ty-1, tz+index0, index) 
					|| is_refined2(tx+1, ty-1, tz+index0, index)
					|| is_refined2(tx,   ty,   tz+index0,   index)
					|| is_refined2(tx+1, ty,   tz+index0, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			else {
				for(index0 = num_edge - 1; index0 >=0; index0--) {
					if(is_refined2(tx,   ty-1, tz+index0, index) 
					|| is_refined2(tx+1, ty-1, tz+index0, index)
					|| is_refined2(tx,   ty,   tz+index0,   index)
					|| is_refined2(tx+1, ty,   tz+index0, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			tx = 2 * tx + 1;	ty *= 2;	tz *= 2;	num_edge *= 2;
		}
		break;

	case 2:
		for(index = level; index <= oct_depth; index++) {
			if(edge_id == 2) {
				for(index0 = 0; index0 < num_edge; index0++) {
					if(is_refined2(tx+index0, ty-1, tz,   index) 
					|| is_refined2(tx+index0, ty-1, tz+1, index)
					|| is_refined2(tx+index0, ty,   tz,   index)
					|| is_refined2(tx+index0, ty,   tz+1, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			else {
				for(index0 = num_edge - 1; index0 >=0; index0--) {
					if(is_refined2(tx+index0, ty-1, tz,   index) 
					|| is_refined2(tx+index0, ty-1, tz+1, index)
					|| is_refined2(tx+index0, ty,   tz,   index)
					|| is_refined2(tx+index0, ty,   tz+1, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			tx *= 2;	ty *= 2;	tz = 2 * tz + 1;	num_edge *= 2;
		}
		break;

	case 3:
		for(index = level; index <= oct_depth; index++) {
			if(edge_id == 3) {
				for(index0 = 0; index0 < num_edge; index0++) {
					if(is_refined2(tx-1, ty-1, tz+index0, index) 
					|| is_refined2(tx,   ty-1, tz+index0, index)
					|| is_refined2(tx,   ty,   tz+index0,   index)
					|| is_refined2(tx-1, ty,   tz+index0, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			else {
				for(index0 = num_edge - 1; index0 >=0; index0--) {
					if(is_refined2(tx-1, ty-1, tz+index0, index) 
					|| is_refined2(tx,   ty-1, tz+index0, index)
					|| is_refined2(tx,   ty,   tz+index0,   index)
					|| is_refined2(tx-1, ty,   tz+index0, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			tx *= 2;	ty *= 2;	tz *= 2;	num_edge *= 2;
		}
		break;

	case 4:
		for(index = level; index <= oct_depth; index++) {
			if(edge_id == 4) {
				for(index0 = 0; index0 < num_edge; index0++) {
					if(is_refined2(tx+index0, ty,   tz-1, index) 
					|| is_refined2(tx+index0, ty+1, tz-1, index)
					|| is_refined2(tx+index0, ty,   tz,   index)
					|| is_refined2(tx+index0, ty+1, tz,   index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			else {
				for(index0 = num_edge - 1; index0 >=0; index0--) {
					if(is_refined2(tx+index0, ty,   tz-1, index) 
					|| is_refined2(tx+index0, ty+1, tz-1, index)
					|| is_refined2(tx+index0, ty,   tz,   index)
					|| is_refined2(tx+index0, ty+1, tz,   index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			tx *= 2;	ty = 2 * ty + 1;	tz *= 2;	num_edge *= 2;
		}
		break;

	case 5:
		for(index = level; index <= oct_depth; index++) {
			if(edge_id == 5) {
				for(index0 = 0; index0 < num_edge; index0++) {
					if(is_refined2(tx+1, ty,   tz+index0, index) 
					|| is_refined2(tx,   ty+1, tz+index0, index)
					|| is_refined2(tx,   ty,   tz+index0,   index)
					|| is_refined2(tx+1, ty+1, tz+index0, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			else {
				for(index0 = num_edge - 1; index0 >=0; index0--) {
					if(is_refined2(tx+1, ty,   tz+index0, index) 
					|| is_refined2(tx,   ty+1, tz+index0, index)
					|| is_refined2(tx,   ty,   tz+index0,   index)
					|| is_refined2(tx+1, ty+1, tz+index0, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			tx = 2 * tx + 1;	ty = 2 * ty + 1;	tz *= 2;	num_edge *= 2;
		}
		break;

	case 6:
		for(index = level; index <= oct_depth; index++) {
			if(edge_id == 6) {
				for(index0 = 0; index0 < num_edge; index0++) {
					if(is_refined2(tx+index0, ty+1, tz,   index) 
					|| is_refined2(tx+index0, ty,   tz+1, index)
					|| is_refined2(tx+index0, ty,   tz,   index)
					|| is_refined2(tx+index0, ty+1, tz+1, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			else {
				for(index0 = num_edge - 1; index0 >=0; index0--) {
					if(is_refined2(tx+index0, ty+1, tz,   index) 
					|| is_refined2(tx+index0, ty,   tz+1, index)
					|| is_refined2(tx+index0, ty,   tz,   index)
					|| is_refined2(tx+index0, ty+1, tz+1, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			tx *= 2;	ty = 2 * ty + 1;	tz = 2 * tz + 1;	num_edge *= 2;
		}
		break;

	case 7:
		for(index = level; index <= oct_depth; index++) {
			if(edge_id == 7) {
				for(index0 = 0; index0 < num_edge; index0++) {
					if(is_refined2(tx-1, ty,   tz+index0, index) 
					|| is_refined2(tx-1, ty+1, tz+index0, index)
					|| is_refined2(tx,   ty,   tz+index0, index)
					|| is_refined2(tx,   ty+1, tz+index0, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			else {
				for(index0 = num_edge - 1; index0 >=0; index0--) {
					if(is_refined2(tx-1, ty,   tz+index0, index) 
					|| is_refined2(tx-1, ty+1, tz+index0, index)
					|| is_refined2(tx,   ty,   tz+index0, index)
					|| is_refined2(tx,   ty+1, tz+index0, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			tx *= 2;	ty = 2 * ty + 1;	tz *= 2;	num_edge *= 2;
		}
		break;

	case 8:
		for(index = level; index <= oct_depth; index++) {
			if(edge_id == 8) {
				for(index0 = 0; index0 < num_edge; index0++) {
					if(is_refined2(tx-1, ty+index0, tz-1, index) 
					|| is_refined2(tx,   ty+index0, tz-1, index)
					|| is_refined2(tx,   ty+index0, tz,   index)
					|| is_refined2(tx-1, ty+index0, tz,   index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			else {
				for(index0 = num_edge - 1; index0 >=0; index0--) {
					if(is_refined2(tx-1, ty+index0, tz-1, index) 
					|| is_refined2(tx,   ty+index0, tz-1, index)
					|| is_refined2(tx,   ty+index0, tz,   index)
					|| is_refined2(tx-1, ty+index0, tz,   index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			tx *= 2;	ty *= 2;	tz *= 2;	num_edge *= 2;
		}
		break;

	case 9:
		for(index = level; index <= oct_depth; index++) {
			if(edge_id == 9) {
				for(index0 = 0; index0 < num_edge; index0++) {
					if(is_refined2(tx,   ty+index0, tz-1, index) 
					|| is_refined2(tx+1, ty+index0, tz-1, index)
					|| is_refined2(tx,   ty+index0, tz,   index)
					|| is_refined2(tx+1, ty+index0, tz,   index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			else {
				for(index0 = num_edge - 1; index0 >=0; index0--) {
					if(is_refined2(tx,   ty+index0, tz-1, index) 
					|| is_refined2(tx+1, ty+index0, tz-1, index)
					|| is_refined2(tx,   ty+index0, tz,   index)
					|| is_refined2(tx+1, ty+index0, tz,   index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			tx = 2 * tx + 1;	ty *= 2;	tz *= 2;	num_edge *= 2;
		}
		break;

	case 10:
		for(index = level; index <= oct_depth; index++) {
			if(edge_id == 10) {
				for(index0 = 0; index0 < num_edge; index0++) {
					if(is_refined2(tx-1, ty+index0, tz,   index) 
					|| is_refined2(tx-1, ty+index0, tz+1, index)
					|| is_refined2(tx,   ty+index0, tz,   index)
					|| is_refined2(tx,   ty+index0, tz+1, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			else {
				for(index0 = num_edge - 1; index0 >=0; index0--) {
					if(is_refined2(tx-1, ty+index0, tz,   index) 
					|| is_refined2(tx-1, ty+index0, tz+1, index)
					|| is_refined2(tx,   ty+index0, tz,   index)
					|| is_refined2(tx,   ty+index0, tz+1, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			tx *= 2;	ty *= 2;	tz = 2 * tz + 1;	num_edge *= 2;
		}
		break;

	case 11:
		for(index = level; index <= oct_depth; index++) {
			if(edge_id == 11) {
				for(index0 = 0; index0 < num_edge; index0++) {
					if(is_refined2(tx+1, ty+index0, tz,   index) 
					|| is_refined2(tx,   ty+index0, tz+1, index)
					|| is_refined2(tx,   ty+index0, tz,   index)
					|| is_refined2(tx+1, ty+index0, tz+1, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			else {
				for(index0 = num_edge - 1; index0 >=0; index0--) {
					if(is_refined2(tx+1, ty+index0, tz,   index) 
					|| is_refined2(tx,   ty+index0, tz+1, index)
					|| is_refined2(tx,   ty+index0, tz,   index)
					|| is_refined2(tx+1, ty+index0, tz+1, index))	my_edge[t] = 1;
					else	my_edge[t] = 0;
					t++;
				}
			}
			tx = 2 * tx + 1;	ty *= 2;	tz = 2 * tz + 1;	num_edge *= 2;
		}
		break;

	}
}

  

void Octree::face_0(int x, int y, int z, int cell_size, int face_id, int v0, int v1, int v2, int v3, int v4, geoframe& geofrm) {
	int my_bool_0 = ((x + y + z) % 2 == 0) && (face_id == 0 || face_id == 2 || face_id == 4 || face_id == 5);
	int my_bool_1 = ((x + y + z) % 2 == 1) && (face_id == 1 || face_id == 3);

	if(my_bool_0 || my_bool_1) {
	//if((x + y + z) % (2*cell_size) == 0) {
		geofrm.AddTetra(v0, v1, v3, v4);
		geofrm.AddTetra(v1, v2, v3, v4);
	}
	else {
		geofrm.AddTetra(v0, v1, v2, v4);
		geofrm.AddTetra(v0, v2, v3, v4);
	}
}

void Octree::face_1(int v0, int v1, int v2, int v3, int v4, unsigned int* array, int m_id, geoframe& geofrm) {
	int t;
	geofrm.AddTetra(v3, array[m_id], v2, v4);
	geofrm.AddTetra(v0, array[0], v3, v4);
	for(t = 0; t < m_id; t++)
		geofrm.AddTetra(array[t], array[t+1], v3, v4);

	t = m_id;
	while(array[t+1] != 999999) {
		geofrm.AddTetra(array[t], array[t+1], v2, v4);
		t++;
	}
	geofrm.AddTetra(array[t], v1, v2, v4);
}

void Octree::face_2_0(int x, int y, int z, int face_id, int v0, int v1, int v2, int v3, int v4, unsigned int* array0, 
					  unsigned int* array1, int m_id0, int m_id1, geoframe& geofrm) {
	int t;
	int my_bool_0 = ((x + y + z) % 2 == 0) && (face_id == 0 || face_id == 2 || face_id == 4 || face_id == 5);
	int my_bool_1 = ((x + y + z) % 2 == 1) && (face_id == 1 || face_id == 3);

	if((x + y + z) % 2 == 0) {
//	if(my_bool_0 || my_bool_1) {
		geofrm.AddTetra(v0, array0[0], v3, v4);
		for(t = 0; t < m_id0; t++)
			geofrm.AddTetra(array0[t], array0[t+1], v3, v4);

		t = m_id0;
		while(array0[t+1] != 999999) {
			geofrm.AddTetra(array0[t], array0[t+1], array1[m_id1], v4);
			t++;
		}
		geofrm.AddTetra(array0[t], v1, array1[m_id1], v4);

		geofrm.AddTetra(v2, array1[0], v1, v4);
		for(t = 0; t < m_id1; t++)
			geofrm.AddTetra(array1[t], array1[t+1], v1, v4);

		t = m_id1;
		while(array1[t+1] != 999999) {
			geofrm.AddTetra(array1[t], array1[t+1], array0[m_id0], v4);
			t++;
		}
		geofrm.AddTetra(array1[t], v3, array0[m_id0], v4);
	}
	else {
		geofrm.AddTetra(v0, array0[0], array1[m_id1], v4);
		for(t = 0; t < m_id0; t++)
			geofrm.AddTetra(array0[t], array0[t+1], array1[m_id1], v4);

		t = m_id0;
		while(array0[t+1] != 999999) {
			geofrm.AddTetra(array0[t], array0[t+1], v2, v4);
			t++;
		}
		geofrm.AddTetra(array0[t], v1, v2, v4);

		geofrm.AddTetra(v2, array1[0], array0[m_id0], v4);
		for(t = 0; t < m_id1; t++)
			geofrm.AddTetra(array1[t], array1[t+1], array0[m_id0], v4);

		t = m_id1;
		while(array1[t+1] != 999999) {
			geofrm.AddTetra(array1[t], array1[t+1], v0, v4);
			t++;
		}
		geofrm.AddTetra(array1[t], v3, v0, v4);
	}
}

void Octree::face_2_1(int v0, int v1, int v2, int v3, int v4, unsigned int* array0,
					  unsigned int* array1, int m_id0, int m_id1, geoframe& geofrm) {
	int t, tt;

	geofrm.AddTetra(v3, array0[m_id0], array1[m_id1], v4);
	geofrm.AddTetra(v0, array0[0], v3, v4);
	for(t = 0; t < m_id0; t++)
		geofrm.AddTetra(array0[t], array0[t+1], v3, v4);

	t = m_id1;
	while(array1[t+1] != 999999) {
		geofrm.AddTetra(array1[t], array1[t+1], v3, v4);
		t++;
	}
	geofrm.AddTetra(array1[t], v2, v3, v4);

	t = m_id0;
	while(array0[t+1] != 999999) {
		geofrm.AddTetra(array0[t], array0[t+1], array1[m_id1], v4);
		t++;
	}

	tt = t;
	geofrm.AddTetra(v1, array1[0], array0[tt], v4);
	for(t = 0; t < m_id1; t++)
		geofrm.AddTetra(array1[t], array1[t+1], array0[tt], v4);

}

void Octree::face_3(int x, int y, int z, int face_id, int cell_size, int v0, int v1, int v2, int v3, int v4, unsigned int* array0, 
		unsigned int* array1, unsigned int* array2, int m_id0, int m_id1, int m_id2, geoframe& geofrm) {
	int t, tt;

	geofrm.AddTetra(array0[m_id0], array1[m_id1], array2[m_id2], v4);

	t = m_id0;
	while(array0[t+1] != 999999) {
		geofrm.AddTetra(array0[t], array0[t+1], array1[m_id1], v4);
		t++;
	}
	tt = t;
	geofrm.AddTetra(v1, array1[0], array0[tt], v4);
	for(t = 0; t < m_id1; t++)
		geofrm.AddTetra(array1[t], array1[t+1], array0[tt], v4);

	for(t = 0; t < m_id2; t++)
		geofrm.AddTetra(array2[t], array2[t+1], array1[m_id1], v4);
	t = m_id1;
	while(array1[t+1] != 999999) {
		geofrm.AddTetra(array1[t], array1[t+1], array2[0], v4);
		t++;
	}
	geofrm.AddTetra(array1[t], v2, array2[0], v4);

	int my_bool_0 = ((x + y + z) % 2 == 0) && (face_id == 0 || face_id == 2 || face_id == 4 || face_id == 5);
	int my_bool_1 = ((x + y + z) % 2 == 1) && (face_id == 1 || face_id == 3);

	if((x + y + z) % 2 == 0) {
	//if(my_bool_0 || my_bool_1) {
		geofrm.AddTetra(v0, array0[0], v3, v4);
		for(t = 0; t < m_id0; t++)
			geofrm.AddTetra(array0[t], array0[t+1], v3, v4);

		t = m_id2;
		while(array2[t+1] != 999999) {
			geofrm.AddTetra(array2[t], array2[t+1], array0[m_id0], v4);
			t++;
		}
		geofrm.AddTetra(array2[t], v3, array0[m_id0], v4);
	}
	else {
		geofrm.AddTetra(v0, array0[0], array2[m_id2], v4);
		for(t = 0; t < m_id0; t++)
			geofrm.AddTetra(array0[t], array0[t+1], array2[m_id2], v4);

		t = m_id2;
		while(array2[t+1] != 999999) {
			geofrm.AddTetra(array2[t], array2[t+1], v0, v4);
			t++;
		}
		geofrm.AddTetra(array2[t], v3, v0, v4);
	}

}

void Octree::face_4(int x, int y, int z, int face_id, int v0, int v1, int v2, int v3, int v4, unsigned int face_center, unsigned int* array0, 
		unsigned int* array1, unsigned int* array2, unsigned int* array3, int m_id0, int m_id1, int m_id2, int m_id3, geoframe& geofrm) {
	int t, tt;

	geofrm.AddTetra(array0[m_id0], array1[m_id1], face_center, v4);
	geofrm.AddTetra(array1[m_id1], array2[m_id2], face_center, v4);
	geofrm.AddTetra(array2[m_id2], array3[m_id3], face_center, v4);
	geofrm.AddTetra(array0[m_id0], face_center, array3[m_id3], v4);

	int my_bool_0 = ((x + y + z) % 2 == 0) && (face_id == 0 || face_id == 2 || face_id == 4 || face_id == 5);
	int my_bool_1 = ((x + y + z) % 2 == 1) && (face_id == 1 || face_id == 3);

//	if((x + y + z) % 2 == 1) {
	//if(my_bool_0 || my_bool_1) {
		t = m_id0;
		while(array0[t+1] != 999999) {
			geofrm.AddTetra(array0[t], array0[t+1], array1[0], v4);
			t++;
		}
		geofrm.AddTetra(array0[t], v1, array1[0], v4);
		for(t = 0; t < m_id1; t++)
			geofrm.AddTetra(array1[t], array1[t+1], array0[m_id0], v4);

		t = m_id1;
		while(array1[t+1] != 999999) {
			geofrm.AddTetra(array1[t], array1[t+1], array2[m_id2], v4);
			t++;
		}
		tt = t;
		geofrm.AddTetra(v2, array2[0], array1[tt], v4);
		for(t = 0; t < m_id2; t++)
			geofrm.AddTetra(array2[t], array2[t+1], array1[tt], v4);

		t = m_id2;
		while(array2[t+1] != 999999) {
			geofrm.AddTetra(array2[t], array2[t+1], array3[0], v4);
			t++;
		}
		tt = t;
		geofrm.AddTetra(v3, array3[0], array2[tt], v4);
		for(t = 0; t < m_id3; t++)
			geofrm.AddTetra(array3[t], array3[t+1], array2[m_id2], v4);

		t = m_id3;
		while(array3[t+1] != 999999) {
			geofrm.AddTetra(array3[t], array3[t+1], array0[m_id0], v4);
			t++;
		}
		tt = t;
		geofrm.AddTetra(v0, array0[0], array3[tt], v4);
		for(t = 0; t < m_id0; t++)
			geofrm.AddTetra(array0[t], array0[t+1], array3[tt], v4);
//	}
//	else {
//	}
}

void Octree::march_edge(int x, int y, int z, int cell_size, int e_id, int num, 
						int* temp, int* index, int& m_id, unsigned int* my_array, geoframe& geofrm){
	int id, t;

	switch (abs(e_id)) {

	case 0:	
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				add_middle_vertex(x, y, z, (t+1.0)/(num+1.0), 0.0, 0.0, cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 100:	// edge -0
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				add_middle_vertex(x+1, y, z, -(t+1.0)/(num+1.0), 0.0, 0.0, cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 1:
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				if(e_id == 1)
					add_middle_vertex(x+1, y, z, 0.0, 0.0, (t+1.0)/(num+1.0), cell_size, my_array[id], geofrm);
				else
					add_middle_vertex(x+1, y, z+1, 0.0, 0.0, -(t+1.0)/(num+1.0), cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 2:
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				if(e_id == 2)
					add_middle_vertex(x, y, z+1, (t+1.0)/(num+1.0), 0.0, 0.0, cell_size, my_array[id], geofrm);
				else
					add_middle_vertex(x+1, y, z+1, -(t+1.0)/(num+1.0), 0.0, 0.0, cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 3:
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				if(e_id == 3)
					add_middle_vertex(x, y, z, 0.0, 0.0, (t+1.0)/(num+1.0), cell_size, my_array[id], geofrm);
				else
					add_middle_vertex(x, y, z+1, 0.0, 0.0, -(t+1.0)/(num+1.0), cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 4:
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				if(e_id == 4)
					add_middle_vertex(x, y+1, z, (t+1.0)/(num+1.0), 0.0, 0.0, cell_size, my_array[id], geofrm);
				else
					add_middle_vertex(x+1, y+1, z, -(t+1.0)/(num+1.0), 0.0, 0.0, cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 5:
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				if(e_id == 5)
					add_middle_vertex(x+1, y+1, z, 0.0, 0.0, (t+1.0)/(num+1.0), cell_size, my_array[id], geofrm);
				else
					add_middle_vertex(x+1, y+1, z+1, 0.0, 0.0, -(t+1.0)/(num+1.0), cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 6:
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				if(e_id == 6)
					add_middle_vertex(x, y+1, z+1, (t+1.0)/(num+1.0), 0.0, 0.0, cell_size, my_array[id], geofrm);
				else
					add_middle_vertex(x+1, y+1, z+1, -(t+1.0)/(num+1.0), 0.0, 0.0, cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 7:
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				if(e_id == 7)
					add_middle_vertex(x, y+1, z, 0.0, 0.0, (t+1.0)/(num+1.0), cell_size, my_array[id], geofrm);
				else
					add_middle_vertex(x, y+1, z+1, 0.0, 0.0, -(t+1.0)/(num+1.0), cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 8:
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				if(e_id == 8)
					add_middle_vertex(x, y, z, 0.0, (t+1.0)/(num+1.0), 0.0, cell_size, my_array[id], geofrm);
				else
					add_middle_vertex(x, y+1, z, 0.0, -(t+1.0)/(num+1.0), 0.0, cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 9:
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				if(e_id == 9)
					add_middle_vertex(x+1, y, z, 0.0, (t+1.0)/(num+1.0), 0.0, cell_size, my_array[id], geofrm);
				else
					add_middle_vertex(x+1, y+1, z, 0.0, -(t+1.0)/(num+1.0), 0.0, cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 10:
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				if(e_id == 10)
					add_middle_vertex(x, y, z+1, 0.0, (t+1.0)/(num+1.0), 0.0, cell_size, my_array[id], geofrm);
				else
					add_middle_vertex(x, y+1, z+1, 0.0, -(t+1.0)/(num+1.0), 0.0, cell_size, my_array[id], geofrm);
			}
		}
		break;

	case 11:
		id = -1;
		for(t = 0; t < num; t++) {
			if(temp[index[t]] == 1) {
				id++;
				if(index[t] == 0)	m_id = id;
				if(e_id == 11)
					add_middle_vertex(x+1, y, z+1, 0.0, (t+1.0)/(num+1.0), 0.0, cell_size, my_array[id], geofrm);
				else
					add_middle_vertex(x+1, y+1, z+1, 0.0, -(t+1.0)/(num+1.0), 0.0, cell_size, my_array[id], geofrm);
			}
		}
		break;

	}
}

void Octree::get_index_array(int level, int& num, int* index) {
	int t;
	int index1[1] = {0};
	int index2[3] = {1, 0, 2};
	int index3[7] = {3, 1, 4, 0, 5, 2, 6};
	int index4[15] = {7, 3, 8, 1, 9, 4, 10, 0, 11, 5, 12, 2, 13, 6, 14};
	int index5[31] = {15,  7, 16, 3, 17,  8, 18, 1, 19,  9, 20, 4, 21, 10, 22, 0, 
					  23, 11, 24, 5, 25, 12, 26, 2, 27, 13, 28, 6, 29, 14, 30};
	int index6[63] = {31, 15, 32,  7, 33, 16, 34, 3, 35, 17, 36,  8, 37, 18, 38, 1,
					 39, 19, 40,  9, 41, 20, 42, 4, 43, 21, 44, 10, 45, 22, 46, 0,
					 47, 23, 48, 11, 49, 24, 50, 5, 51, 25, 52, 12, 53, 26, 54, 2,
					 55, 27, 56, 13, 57, 28, 58, 6, 59, 29, 60, 14, 61, 30, 62};

	if(oct_depth - level == 1) {num = 1;	for(t = 0; t < num; t++) index[t] = index1[t];}
	if(oct_depth - level == 2) {num = 3;	for(t = 0; t < num; t++) index[t] = index2[t];}
	if(oct_depth - level == 3) {num = 7;	for(t = 0; t < num; t++) index[t] = index3[t];}
	if(oct_depth - level == 4) {num = 15;	for(t = 0; t < num; t++) index[t] = index4[t];}
	if(oct_depth - level == 5) {num = 31;	for(t = 0; t < num; t++) index[t] = index5[t];}
	if(oct_depth - level == 6) {num = 63;	for(t = 0; t < num; t++) index[t] = index6[t];}
}

void Octree::get_middle_array_1(int face_id, int* my_edge0, int* my_edge1, int* my_edge2, int* my_edge3, 
				unsigned int* my_array, int& m_id, int x, int y, int z, int level, geoframe& geofrm) {
	int temp[128], t, e_id, num, index[128];

	int cell_size = (dim[0]-1)/(1<<level);

	get_index_array(level, num, index);

	for(t = 0; t < 128; t++) temp[t] = -1;

	if(face_id == 0) {
		if(my_edge0[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge0[t]; e_id = 3;}
		if(my_edge1[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge1[t]; e_id = 10;}
		if(my_edge2[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge2[t]; e_id = -7;}
		if(my_edge3[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge3[t]; e_id = -8;}
	}
	if(face_id == 1) {
		if(my_edge0[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge0[t]; e_id = 9;}
		if(my_edge1[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge1[t]; e_id = 5;}
		if(my_edge2[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge2[t]; e_id = -11;}
		if(my_edge3[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge3[t]; e_id = -1;}
	}
	if(face_id == 2) {
		if(my_edge0[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge0[t]; e_id = 0;}
		if(my_edge1[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge1[t]; e_id = 1;}
		if(my_edge2[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge2[t]; e_id = -2;}
		if(my_edge3[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge3[t]; e_id = -3;}
	}
	if(face_id == 3) {
		if(my_edge0[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge0[t]; e_id = 7;}
		if(my_edge1[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge1[t]; e_id = 6;}
		if(my_edge2[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge2[t]; e_id = -5;}
		if(my_edge3[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge3[t]; e_id = -4;}
	}
	if(face_id == 4) {
		if(my_edge0[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge0[t]; e_id = 8;}
		if(my_edge1[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge1[t]; e_id = 4;}
		if(my_edge2[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge2[t]; e_id = -9;}
		if(my_edge3[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge3[t]; e_id = -100;}
	}
	if(face_id == 5) {
		if(my_edge0[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge0[t]; e_id = 11;}
		if(my_edge1[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge1[t]; e_id = -6;}
		if(my_edge2[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge2[t]; e_id = -10;}
		if(my_edge3[0] == 1) { for(t = 0; t < num; t++) temp[t] = my_edge3[t]; e_id = 2;}
	}
	
	march_edge(x, y, z, cell_size, e_id, num, temp, index, m_id, my_array, geofrm);
}

void Octree::get_middle_array_2(int face_id, int* my_edge0, int* my_edge1, int* my_edge2, int* my_edge3, 
				unsigned int* my_array0, unsigned int* my_array1, int& m_id0, int& m_id1, 
				int x, int y, int z, int level, geoframe& geofrm) {

	int temp[2][128], my_temp[128], t, e_id[2], num, index[128], s0, s1, s2, s3;

	int cell_size = (dim[0]-1)/(1<<level);

	get_index_array(level, num, index);

	s0 = my_edge0[0];	s1 = my_edge1[0];	s2 = my_edge2[0];	s3 = my_edge3[0];
	for(t = 0; t < 128; t++) {temp[0][t] = -1;	temp[1][t] = -1;}

	if(s0 == 1 && s1 == 1) {
		for(t = 0; t < num; t++) {temp[0][t] = my_edge0[t];	temp[1][t] = my_edge1[t]; }
		if(face_id == 0) {e_id[0] = 3;	e_id[1] = 10;}
		if(face_id == 1) {e_id[0] = 9;	e_id[1] = 5;}
		if(face_id == 2) {e_id[0] = 0;	e_id[1] = 1;}
		if(face_id == 3) {e_id[0] = 7;	e_id[1] = 6;}
		if(face_id == 4) {e_id[0] = 8;	e_id[1] = 4;}
		if(face_id == 5) {e_id[0] = 11;	e_id[1] = -6;}
	} 
	if(s0 == 1 && s2 == 1) {
		for(t = 0; t < num; t++) {temp[0][t] = my_edge0[t];	temp[1][t] = my_edge2[t]; }
		if(face_id == 0) {e_id[0] = 3;	e_id[1] = -7;}
		if(face_id == 1) {e_id[0] = 9;	e_id[1] = -11;}
		if(face_id == 2) {e_id[0] = 0;	e_id[1] = -2;}
		if(face_id == 3) {e_id[0] = 7;	e_id[1] = -5;}
		if(face_id == 4) {e_id[0] = 8;	e_id[1] = -9;}
		if(face_id == 5) {e_id[0] = 11;	e_id[1] = -10;}
	}
	if(s3 == 1 && s0 == 1) {
		for(t = 0; t < num; t++) {temp[0][t] = my_edge3[t];	temp[1][t] = my_edge0[t]; }
		if(face_id == 0) {e_id[0] = -8;	e_id[1] = 3;}
		if(face_id == 1) {e_id[0] = -1;	e_id[1] = 9;}
		if(face_id == 2) {e_id[0] = -3;	e_id[1] = 0;}
		if(face_id == 3) {e_id[0] = -4;	e_id[1] = 7;}
		if(face_id == 4) {e_id[0] = -100;	e_id[1] = 8;}
		if(face_id == 5) {e_id[0] = 2;		e_id[1] = 11;}
	}
	if(s1 == 1 && s2 == 1) {
		for(t = 0; t < num; t++) {temp[0][t] = my_edge1[t];	temp[1][t] = my_edge2[t]; }
		if(face_id == 0) {e_id[0] = 10;	e_id[1] = -7;}
		if(face_id == 1) {e_id[0] = 5;	e_id[1] = -11;}
		if(face_id == 2) {e_id[0] = 1;	e_id[1] = -2;}
		if(face_id == 3) {e_id[0] = 6;	e_id[1] = -5;}
		if(face_id == 4) {e_id[0] = 4;	e_id[1] = -9;}
		if(face_id == 5) {e_id[0] = -6;	e_id[1] = -10;}
	}
	if(s1 == 1 && s3 == 1) {
		for(t = 0; t < num; t++) {temp[0][t] = my_edge1[t];	temp[1][t] = my_edge3[t]; }
		if(face_id == 0) {e_id[0] = 10;	e_id[1] = -8;}
		if(face_id == 1) {e_id[0] = 5;	e_id[1] = -1;}
		if(face_id == 2) {e_id[0] = 1;	e_id[1] = -3;}
		if(face_id == 3) {e_id[0] = 6;	e_id[1] = -4;}
		if(face_id == 4) {e_id[0] = 4;	e_id[1] = -100;}
		if(face_id == 5) {e_id[0] = -6;	e_id[1] = 2;}
	}
	if(s2 == 1 && s3 == 1) {
		for(t = 0; t < num; t++) {temp[0][t] = my_edge2[t];	temp[1][t] = my_edge3[t]; }
		if(face_id == 0) {e_id[0] = -7;		e_id[1] = -8;}
		if(face_id == 1) {e_id[0] = -11;	e_id[1] = -1;}
		if(face_id == 2) {e_id[0] = -2;		e_id[1] = -3;}
		if(face_id == 3) {e_id[0] = -5;		e_id[1] = -4;}
		if(face_id == 4) {e_id[0] = -9;		e_id[1] = -100;}
		if(face_id == 5) {e_id[0] = -10;	e_id[1] = 2;}
	}

	for(t = 0; t < 128; t++)	my_temp[t] = -1;
	for(t = 0; t < num; t++)	my_temp[t] = temp[0][t];
	march_edge(x, y, z, cell_size, e_id[0], num, my_temp, index, m_id0, my_array0, geofrm);
	for(t = 0; t < num; t++)	my_temp[t] = temp[1][t];
	march_edge(x, y, z, cell_size, e_id[1], num, my_temp, index, m_id1, my_array1, geofrm);

}

void Octree::get_middle_array_3(int face_id, int* my_edge0, int* my_edge1, int* my_edge2, int* my_edge3, 
				unsigned int* my_array0, unsigned int* my_array1, unsigned int* my_array2,
				int& m_id0, int& m_id1, int& m_id2, int x, int y, int z, int level, geoframe& geofrm) {

	int temp[3][128], my_temp[128], t, e_id[3], num, index[128], s0, s1, s2, s3;

	int cell_size = (dim[0]-1)/(1<<level);

	get_index_array(level, num, index);

	s0 = my_edge0[0];	s1 = my_edge1[0];	s2 = my_edge2[0];	s3 = my_edge3[0];
	for(t = 0; t < 128; t++) {temp[0][t] = -1;	temp[1][t] = -1;	temp[2][t] = -1;}

	if(s3 == 0) {
		for(t = 0; t < num; t++) {
			temp[0][t] = my_edge0[t];	temp[1][t] = my_edge1[t];	temp[2][t] = my_edge2[t];
		}
		if(face_id == 0) {e_id[0] = 3;		e_id[1] = 10;	e_id[2] = -7;}
		if(face_id == 1) {e_id[0] = 9;		e_id[1] = 5;	e_id[2] = -11;}
		if(face_id == 2) {e_id[0] = 0;		e_id[1] = 1;	e_id[2] = -2;}
		if(face_id == 3) {e_id[0] = 7;		e_id[1] = 6;	e_id[2] = -5;}
		if(face_id == 4) {e_id[0] = 8;		e_id[1] = 4;	e_id[2] = -9;}
		if(face_id == 5) {e_id[0] = 11;		e_id[1] = -6;	e_id[2] = -10;}
	} 
	if(s2 == 0) {
		for(t = 0; t < num; t++) {
			temp[0][t] = my_edge3[t];	temp[1][t] = my_edge0[t];	temp[2][t] = my_edge1[t];
		}
		if(face_id == 0) {e_id[0] = -8;		e_id[1] = 3;	e_id[2] = 10;}
		if(face_id == 1) {e_id[0] = -1;		e_id[1] = 9;	e_id[2] = 5;}
		if(face_id == 2) {e_id[0] = -3;		e_id[1] = 0;	e_id[2] = 1;}
		if(face_id == 3) {e_id[0] = -4;		e_id[1] = 7;	e_id[2] = 6;}
		if(face_id == 4) {e_id[0] = -100;	e_id[1] = 8;	e_id[2] = 4;}
		if(face_id == 5) {e_id[0] = 2;		e_id[1] = 11;	e_id[2] = -6;}
	}
	if(s1 == 0) {
		for(t = 0; t < num; t++) {
			temp[0][t] = my_edge2[t];	temp[1][t] = my_edge3[t];	temp[2][t] = my_edge0[t];
		}
		if(face_id == 0) {e_id[0] = -7;		e_id[1] = -8;	e_id[2] = 3;}
		if(face_id == 1) {e_id[0] = -11;	e_id[1] = -1;	e_id[2] = 9;}
		if(face_id == 2) {e_id[0] = -2;		e_id[1] = -3;	e_id[2] = 0;}
		if(face_id == 3) {e_id[0] = -5;		e_id[1] = -4;	e_id[2] = 7;}
		if(face_id == 4) {e_id[0] = -9;		e_id[1] = -100;	e_id[2] = 8;}
		if(face_id == 5) {e_id[0] = -10;	e_id[1] = 2;	e_id[2] = 11;}
	}
	if(s0 == 0) {
		for(t = 0; t < num; t++) {
			temp[0][t] = my_edge1[t];	temp[1][t] = my_edge2[t];	temp[2][t] = my_edge3[t];
		}
		if(face_id == 0) {e_id[0] = 10;		e_id[1] = -7;	e_id[2] = -8;}
		if(face_id == 1) {e_id[0] = 5;		e_id[1] = -11;	e_id[2] = -1;}
		if(face_id == 2) {e_id[0] = 1;		e_id[1] = -2;	e_id[2] = -3;}
		if(face_id == 3) {e_id[0] = 6;		e_id[1] = -5;	e_id[2] = -4;}
		if(face_id == 4) {e_id[0] = 4;		e_id[1] = -9;	e_id[2] = -100;}
		if(face_id == 5) {e_id[0] = -6;		e_id[1] = -10;	e_id[2] = 2;}
	}

	for(t = 0; t < 128; t++)		my_temp[t] = -1;
	for(t = 0; t < num; t++)	my_temp[t] = temp[0][t];
	march_edge(x, y, z, cell_size, e_id[0], num, my_temp, index, m_id0, my_array0, geofrm);
	for(t = 0; t < num; t++)	my_temp[t] = temp[1][t];
	march_edge(x, y, z, cell_size, e_id[1], num, my_temp, index, m_id1, my_array1, geofrm);
	for(t = 0; t < num; t++)	my_temp[t] = temp[2][t];
	march_edge(x, y, z, cell_size, e_id[2], num, my_temp, index, m_id2, my_array2, geofrm);

}

void Octree::get_middle_array_4(int face_id, int* my_edge0, int* my_edge1, int* my_edge2, int* my_edge3, 
				unsigned int* my_array0, unsigned int* my_array1, unsigned int* my_array2, unsigned int* my_array3,
				int& m_id0, int& m_id1, int& m_id2, int& m_id3, unsigned int& face_center, int x, int y, int z, int level, geoframe& geofrm) {

	int temp[4][128], my_temp[128], t, e_id[4], num, index[128], s0, s1, s2, s3;

	int cell_size = (dim[0]-1)/(1<<level);

	get_index_array(level, num, index);

	s0 = my_edge0[0];	s1 = my_edge1[0];	s2 = my_edge2[0];	s3 = my_edge3[0];
	for(t = 0; t < 128; t++) {temp[0][t] = -1;	temp[1][t] = -1;	temp[2][t] = -1;}

	for(t = 0; t < num; t++) {
		temp[0][t] = my_edge0[t];	temp[1][t] = my_edge1[t];
		temp[2][t] = my_edge2[t];	temp[3][t] = my_edge3[t];
	}
	if(face_id == 0) {
		e_id[0] = 3;		e_id[1] = 10;	e_id[2] = -7;	e_id[3] = -8;
		add_middle_vertex(x, y, z, 0.0, 0.5, 0.5, cell_size, face_center, geofrm);
	}
	if(face_id == 1) {
		e_id[0] = 9;		e_id[1] = 5;	e_id[2] = -11;	e_id[3] = -1;
		add_middle_vertex(x+1, y, z, 0.0, 0.5, 0.5, cell_size, face_center, geofrm);
	}
	if(face_id == 2) {
		e_id[0] = 0;		e_id[1] = 1;	e_id[2] = -2;	e_id[3] = -3;
		add_middle_vertex(x, y, z, 0.5, 0.0, 0.5, cell_size, face_center, geofrm);
	}
	if(face_id == 3) {
		e_id[0] = 7;		e_id[1] = 6;	e_id[2] = -5;	e_id[3] = -4;
		add_middle_vertex(x, y+1, z, 0.5, 0.0, 0.5, cell_size, face_center, geofrm);
	}
	if(face_id == 4) {
		e_id[0] = 8;		e_id[1] = 4;	e_id[2] = -9;	e_id[3] = -100;
		add_middle_vertex(x, y, z, 0.5, 0.5, 0.0, cell_size, face_center, geofrm);
	}
	if(face_id == 5) {
		e_id[0] = 11;		e_id[1] = -6;	e_id[2] = -10;	e_id[3] = 2;
		add_middle_vertex(x, y, z+1, 0.5, 0.5, 0.0, cell_size, face_center, geofrm);
	}

	for(t = 0; t < 128; t++)		my_temp[t] = -1;
	for(t = 0; t < num; t++)	my_temp[t] = temp[0][t];
	march_edge(x, y, z, cell_size, e_id[0], num, my_temp, index, m_id0, my_array0, geofrm);
	for(t = 0; t < num; t++)	my_temp[t] = temp[1][t];
	march_edge(x, y, z, cell_size, e_id[1], num, my_temp, index, m_id1, my_array1, geofrm);
	for(t = 0; t < num; t++)	my_temp[t] = temp[2][t];
	march_edge(x, y, z, cell_size, e_id[2], num, my_temp, index, m_id2, my_array2, geofrm);
	for(t = 0; t < num; t++)	my_temp[t] = temp[3][t];
	march_edge(x, y, z, cell_size, e_id[3], num, my_temp, index, m_id3, my_array3, geofrm);

}

void Octree::permute_1(int& i0, int& i1, int& i2, int& i3, int s0, int s1, int s2, int s3) {
	int t0, t1, t2, t3;

	t0 = i0;	t1 = i1;	t2 = i2;	t3 = i3;

	if(s1 == 1) {
		i0 = t1;	i1 = t2;	i2 = t3;	i3 = t0;
	}
	if(s2 == 1) {
		i0 = t2;	i1 = t3;	i2 = t0;	i3 = t1;
	}
	if(s3 == 1) {
		i0 = t3;	i1 = t0;	i2 = t1;	i3 = t2;
	}
}

void Octree::permute_2(int& i0, int& i1, int& i2, int& i3, int& s0, int& s1, int& s2, int& s3) {
	int t0, t1, t2, t3, ts0, ts1, ts2, ts3;

	t0 = i0;	t1 = i1;	t2 = i2;	t3 = i3;
	ts0 = s0;	ts1 = s1;	ts2 = s2;	ts3 = s3;

	if(ts3 == 1 && ts0 == 1) {
		i0 = t3;	i1 = t0;	i2 = t1;	i3 = t2;
		s0 = ts3;	s1 = ts0;	s2 = ts1;	s3 = ts2;
	}
	if((ts1 == 1 && ts2 == 1) || (ts1 == 1 && ts3 == 1)) {
		i0 = t1;	i1 = t2;	i2 = t3;	i3 = t0;
		s0 = ts1;	s1 = ts2;	s2 = ts3;	s3 = ts0;
	}
	if(ts2 == 1 && ts3 == 1) {
		i0 = t2;	i1 = t3;	i2 = t0;	i3 = t1;
		s0 = ts2;	s1 = ts3;	s2 = ts0;	s3 = ts1;
	}
}

void Octree::permute_3(int& i0, int& i1, int& i2, int& i3, int s0, int s1, int s2, int s3) {
	int t0, t1, t2, t3;

	t0 = i0;	t1 = i1;	t2 = i2;	t3 = i3;

	if(s0 == 0) {
		i0 = t1;	i1 = t2;	i2 = t3;	i3 = t0;
	}
	if(s1 == 0) {
		i0 = t2;	i1 = t3;	i2 = t0;	i3 = t1;
	}
	if(s2 == 0) {
		i0 = t3;	i1 = t0;	i2 = t1;	i3 = t2;
	}
}

void Octree::march_one_face(int face_id, int oc_id, int level, unsigned int* my_vtx, unsigned int center, geoframe& geofrm) {

	int x, y, z, tx, ty, tz, index, my_oc_id, sum_edge, v_id[4], e_id[4], neighbor;
	unsigned int my_array0[128], my_array1[128], my_array2[128], my_array3[128];
	int my_edge0[128], my_edge1[128], my_edge2[128], my_edge3[128];
	int i0, i1, i2, i3, s0, s1, s2, s3, m_id0, m_id1, m_id2;
	//int m_id3, face_center;
	float val[8];

	int cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_id, x, y, z, level);
	getCellValues(oc_id, level, val);

	for(index = 0; index < 128; index++)	{
		my_array0[index] = 999999;	my_array1[index] = 999999;
		my_array2[index] = 999999;	my_array3[index] = 999999;
	}

	if(face_id == 0) {
		v_id[0] = 0;	v_id[1] = 3;	v_id[2] = 7;	v_id[3] = 4;
		e_id[0] = 3;	e_id[1] = 10;	e_id[2] = -7;	e_id[3] = -8;
		tx = x;			ty = y;			tz = z;
		neighbor = 0;
	}
	if(face_id == 1) {
		v_id[0] = 1;	v_id[1] = 5;	v_id[2] = 6;	v_id[3] = 2;
		e_id[0] = 9;	e_id[1] = 5;	e_id[2] = -11;	e_id[3] = -1;
		tx = x + 1;		ty = y;			tz = z;
		neighbor = 0;
	}
	if(face_id == 2) {
		v_id[0] = 0;	v_id[1] = 1;	v_id[2] = 2;	v_id[3] = 3;
		e_id[0] = 0;	e_id[1] = 1;	e_id[2] = -2;	e_id[3] = -3;
		tx = x;			ty = y;			tz = z;
		neighbor = 2;
	}
	if(face_id == 3) {
		v_id[0] = 4;	v_id[1] = 7;	v_id[2] = 6;	v_id[3] = 5;
		e_id[0] = 7;	e_id[1] = 6;	e_id[2] = -5;	e_id[3] = -4;
		tx = x;			ty = y+1;			tz = z;
		neighbor = 2;
	}
	if(face_id == 4) {
		v_id[0] = 0;	v_id[1] = 4;	v_id[2] = 5;	v_id[3] = 1;
		e_id[0] = 8;	e_id[1] = 4;	e_id[2] = -9;	e_id[3] = -100;
		tx = x;			ty = y;			tz = z;
		neighbor = 4;
	}
	if(face_id == 5) {
		v_id[0] = 2;	v_id[1] = 6;	v_id[2] = 7;	v_id[3] = 3;
		e_id[0] = 11;	e_id[1] = -6;	e_id[2] = -10;	e_id[3] = 2;
		tx = x;			ty = y;			tz = z + 1;
		neighbor = 4;
	}

	march_each_edge(oc_id, level, e_id[0],  my_edge0, geofrm);
	march_each_edge(oc_id, level, e_id[1],  my_edge1, geofrm);
	march_each_edge(oc_id, level, e_id[2],  my_edge2, geofrm);
	march_each_edge(oc_id, level, e_id[3],  my_edge3, geofrm);

	i0 = my_vtx[v_id[0]];	i1 = my_vtx[v_id[1]];
	i2 = my_vtx[v_id[2]];	i3 = my_vtx[v_id[3]];
	s0 = my_edge0[0];		s1 = my_edge1[0];	s2 = my_edge2[0];	s3 = my_edge3[0];

	sum_edge = s0 + s1 + s2 + s3;
	assert(sum_edge >= 0 && sum_edge <= 4);

	int my_bool = ( val[v_id[0]] <= iso_val && val[v_id[1]] <= iso_val &&
					val[v_id[2]] <= iso_val && val[v_id[3]] <= iso_val );
	if(sum_edge == 0 && my_bool) {
		face_0(x, y, z, cell_size, face_id, i0, i1, i2, i3, center, geofrm);
	}
	else if(sum_edge == 1 && my_bool) {
		permute_1(i0, i1, i2, i3, s0, s1, s2, s3);
		get_middle_array_1(face_id, my_edge0, my_edge1, my_edge2, my_edge3, my_array0, m_id0, x, y, z, level, geofrm);
		face_1(i0, i1, i2, i3, center, my_array0, m_id0, geofrm);
	}
	else if(sum_edge == 2 && my_bool) {
		permute_2(i0, i1, i2, i3, s0, s1, s2, s3);
		get_middle_array_2(face_id, my_edge0, my_edge1, my_edge2, my_edge3, my_array0, my_array1, 
						m_id0, m_id1, x, y, z, level, geofrm);
		if(s1 == 0)
			face_2_0(x, y, z, face_id, i0, i1, i2, i3, center, my_array0, my_array1, m_id0, m_id1, geofrm);
		else
			face_2_1(i0, i1, i2, i3, center, my_array0, my_array1, m_id0, m_id1, geofrm);
	}
	else if(sum_edge == 3 && my_bool) {
		permute_3(i0, i1, i2, i3, s0, s1, s2, s3);
		get_middle_array_3(face_id, my_edge0, my_edge1, my_edge2, my_edge3, my_array0, my_array1, 
						my_array2, m_id0, m_id1, m_id2, x, y, z, level, geofrm);
		face_3(x, y, z, face_id, cell_size, i0, i1, i2, i3, center, my_array0, my_array1, my_array2, m_id0, m_id1, m_id2, geofrm);
	}
	else if(sum_edge == 4) {			// sum_edge == 4
		//get_middle_array_4(face_id, my_edge0, my_edge1, my_edge2, my_edge3, my_array0, my_array1, 
		//				my_array2, my_array3, m_id0, m_id1, m_id2, m_id3, face_center, x, y, z, level, geofrm);
		//face_4(x, y, z, face_id, i0, i1, i2, i3, center, face_center, my_array0, my_array1, my_array2, my_array3, m_id0, m_id1, m_id2, m_id3, geofrm);
		for(index = 0; index < 4; index++) {
			int new_level = level;
			if(new_level < oct_depth) {
				if(face_id == 0 || face_id == 1) { 
					if(index == 0) my_oc_id = xyz2octcell(2*tx, 2*ty, 2*tz, ++new_level);
					else if(index == 1) my_oc_id = xyz2octcell(2*tx, 2*ty+1, 2*tz, ++new_level);
					else if(index == 2) my_oc_id = xyz2octcell(2*tx, 2*ty, 2*tz+1, ++new_level);
					else my_oc_id = xyz2octcell(2*tx, 2*ty+1, 2*tz+1, ++new_level);
				}
				if(face_id == 2 || face_id == 3) { 
					if(index == 0) my_oc_id = xyz2octcell(2*tx, 2*ty, 2*tz, ++new_level);
					else if(index == 1) my_oc_id = xyz2octcell(2*tx+1, 2*ty, 2*tz, ++new_level);
					else if(index == 2) my_oc_id = xyz2octcell(2*tx, 2*ty, 2*tz+1, ++new_level);
					else my_oc_id = xyz2octcell(2*tx+1, 2*ty, 2*tz+1, ++new_level);
				}
				if(face_id == 4 || face_id == 5) { 
					if(index == 0) my_oc_id = xyz2octcell(2*tx, 2*ty, 2*tz, ++new_level);
					else if(index == 1) my_oc_id = xyz2octcell(2*tx+1, 2*ty, 2*tz, ++new_level);
					else if(index == 2) my_oc_id = xyz2octcell(2*tx, 2*ty+1, 2*tz, ++new_level);
					else my_oc_id = xyz2octcell(2*tx+1, 2*ty+1, 2*tz, ++new_level);
				}
				march_each_face(my_oc_id, new_level, neighbor, center, geofrm);
			}
		}
	}

}

void Octree::march_each_face(int oc_id, int level, int face_id, unsigned int center, geoframe& geofrm) {

	int x, y, z;
	unsigned int my_vtx[8];

	int cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_id, x, y, z, level);

	add_one_vertex(x,   y,   z,   cell_size, my_vtx[0], geofrm);
	add_one_vertex(x,   y,   z+1, cell_size, my_vtx[3], geofrm);
	add_one_vertex(x,   y+1, z+1, cell_size, my_vtx[7], geofrm);
	add_one_vertex(x,   y+1, z,   cell_size, my_vtx[4], geofrm);

	add_one_vertex(x+1, y,   z,   cell_size, my_vtx[1], geofrm);
	add_one_vertex(x+1, y,   z+1, cell_size, my_vtx[2], geofrm);
	add_one_vertex(x+1, y+1, z+1, cell_size, my_vtx[6], geofrm);
	add_one_vertex(x+1, y+1, z,   cell_size, my_vtx[5], geofrm);

	march_one_face(face_id, oc_id, level, my_vtx, center, geofrm);

}


}

