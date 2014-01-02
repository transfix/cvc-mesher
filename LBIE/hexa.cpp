#include<cellqueue.h>
#include"octree.h"
#include"e_face.h"
#include"cubes.h"
#include<stdio.h>
#include"vtkMarchingCubesCases.h"
#include"pcio.h"
#include"LBIE_geoframe.h"
//#include<qgl.h>
#include<assert.h>
#include<boost/array.hpp>
#include<algorithm>

namespace LBIE
{

void Octree::hexahedralize(geoframe& geofrm, float err_tol) {

	int x, y, z, valid_leaf, cell_size, level;
	int oc_id[8], edge_id[6], i, j, k, flag_method;
	unsigned int vtx[8];
	float val[8];

	for(i = 0; i < octcell_num; i++) vtx_idx_arr[i] = -1;

	flag_method = 2;	// 0 - uniform;	1 - Method 1;	2 - Method 2.
	if(flag_method == 2) assign_refine_sign_hexa(geofrm, err_tol);

	for(i = 0; i < leaf_num; i++ ) {

		valid_leaf = cut_array[i] ;
		level = get_level(valid_leaf) ;
		cell_size = (dim[0]-1)/(1<<level);
		octcell2xyz(valid_leaf, x, y, z, level);
		getCellValues(valid_leaf, level, val);

		//interior vertex at a uniform starting level
		for (j = 0; j < 8; j++ ) {
			if (is_vflag_on(x, y, z, level, j)) continue;
			if( minmax[valid_leaf].min <= iso_val && (val[j] < iso_val) ) {
				if(is_min_vertex(valid_leaf, j, vtx, geofrm)) {
					vflag_on(x, y, z, level, j);
					find_oc_id_hexa(x, y, z, level, j, oc_id);
					for(k = 0; k < 6; k++) edge_id[k] = 0;
					find_edge_id_hexa(x, y, z, cell_size, j, edge_id);

					if(flag_method == 0) add_hexa(geofrm, vtx);					// uniform hexa
					else if(flag_method == 1)	
						hexa_adaptive_1(geofrm, oc_id, edge_id, err_tol, vtx);	// method 1
					else 
						hexa_adaptive_2(geofrm, oc_id, edge_id, err_tol, vtx);	// method 2
				}
			}
		}

	}

}

// adaptive hexa mesh generation -- method 1
void Octree::hexa_adaptive_1(geoframe& geofrm, int* oc_id, int* edge_id, float err_tol, unsigned int* vtx) {

	int num_id, num_quad, i, j, k, v[4];
	unsigned int vtx_new[8], vtx_temp[8], vtx_t[8];
	int level, cell_size, xx, yy, zz;
	float x, y, z;

	num_id = 0;
	for(i = 0; i < 8; i++) {
		if(get_err_grad(oc_id[i]) > err_tol) num_id++;
		if(is_skipcell(oc_id[i]) == 0) get_vtx_new(geofrm, oc_id[i], vtx[i]);
	}
	num_quad = 0;
	for(i = 0; i < 6; i++) {
		if(edge_id[i] == 1) num_quad++;
	}
	if(num_id == 0 || num_quad == 0) {
		add_hexa(geofrm, vtx);
	}
	else {

	  geofrm.AddVert_hexa_adaptive_1_center(vtx, vtx_new);
		add_hexa(geofrm, vtx_new);

		for(i = 0; i < 6; i++) {
			if(i == 0) {v[0] = 0;	v[1] = 1;	v[2] = 2;	v[3] = 3;}	// up
			if(i == 1) {v[0] = 4;	v[1] = 7;	v[2] = 6;	v[3] = 5;}	// down
			if(i == 2) {v[0] = 0;	v[1] = 3;	v[2] = 7;	v[3] = 4;}	// left
			if(i == 3) {v[0] = 1;	v[1] = 5;	v[2] = 6;	v[3] = 2;}	// right
			if(i == 4) {v[0] = 1;	v[1] = 0;	v[2] = 4;	v[3] = 5;}	// front
			if(i == 5) {v[0] = 3;	v[1] = 2;	v[2] = 6;	v[3] = 7;}	// back

			for(j = 0; j < 4; j++) {
				vtx_temp[j] = vtx[v[j]];
				vtx_temp[j+4] = vtx_new[v[j]];
			}
			if((get_err_grad(oc_id[v[0]]) > err_tol || get_err_grad(oc_id[v[1]]) > err_tol ||
				get_err_grad(oc_id[v[2]]) > err_tol || get_err_grad(oc_id[v[3]]) > err_tol) &&
				edge_id[i] == 1) {

				geofrm.AddVert_hexa_adaptive_1_top(vtx_temp, vtx_t);
				
				for(k = 0; k < 4; k++) {
					for(j = 0; j < 4; j++) {
						level = get_level(oc_id[v[j]]) ;
						cell_size = (dim[0]-1)/(1<<level);
						octcell2xyz(oc_id[v[j]], xx, yy, zz, level);
						x = geofrm.verts[vtx_t[k]][0]/cell_size - (float)xx;
						y = geofrm.verts[vtx_t[k]][1]/cell_size - (float)yy;
						z = geofrm.verts[vtx_t[k]][2]/cell_size - (float)zz;
						if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
					}
					if(j < 4) get_vtx_new(geofrm, oc_id[v[j]], vtx_t[k]);
				}
				
				hexa_adaptive_1_top(geofrm, vtx_temp, vtx_t);
			}
			else
				add_hexa(geofrm, vtx_temp);
		}
	}

}

// adaptive hexa mesh generation -- method 2
void Octree::hexa_adaptive_2(geoframe& geofrm, int* oc_id, int* edge_id, float err_tol, unsigned int* vtx) {

	int num_id, i, j, my_bool_2, my_bool_4, edge_id_new[6];
	int level, cell_size, xx, yy, zz, vv[4];
	unsigned int vtx_new[64], vtx_temp[8];
	float x, y, z;

	for(i = 0; i < 6; i++) {
		if(i == 0) {vv[0] = 0;	vv[1] = 1;	vv[2] = 2;	vv[3] = 3;}	// up
		if(i == 1) {vv[0] = 4;	vv[1] = 7;	vv[2] = 6;	vv[3] = 5;}	// down
		if(i == 2) {vv[0] = 0;	vv[1] = 3;	vv[2] = 7;	vv[3] = 4;}	// left
		if(i == 3) {vv[0] = 1;	vv[1] = 5;	vv[2] = 6;	vv[3] = 2;}	// right
		if(i == 4) {vv[0] = 1;	vv[1] = 0;	vv[2] = 4;	vv[3] = 5;}	// front
		if(i == 5) {vv[0] = 3;	vv[1] = 2;	vv[2] = 6;	vv[3] = 7;}	// back

		if(edge_id[i] == 1) {
			geofrm.AddBound_edge(vtx[vv[0]], vtx[vv[1]]);
			geofrm.AddBound_edge(vtx[vv[1]], vtx[vv[2]]);
			geofrm.AddBound_edge(vtx[vv[2]], vtx[vv[3]]);
			geofrm.AddBound_edge(vtx[vv[3]], vtx[vv[0]]);
		}
	}

	num_id = 0;
	for(i = 0; i < 8; i++) {
		if(vtx_idx_arr_refine[oc_id[i]] == 1) num_id++;
		if(is_skipcell(oc_id[i]) == 0) get_vtx_new(geofrm, oc_id[i], vtx[i]);
	}

	// detect one edge
	my_bool_2 = (vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[1]] == 1) ||
				(vtx_idx_arr_refine[oc_id[1]] == 1 && vtx_idx_arr_refine[oc_id[2]] == 1) ||
				(vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1) ||
				(vtx_idx_arr_refine[oc_id[3]] == 1 && vtx_idx_arr_refine[oc_id[0]] == 1) ||
				(vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1) ||
				(vtx_idx_arr_refine[oc_id[5]] == 1 && vtx_idx_arr_refine[oc_id[6]] == 1) ||
				(vtx_idx_arr_refine[oc_id[6]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1) ||
				(vtx_idx_arr_refine[oc_id[7]] == 1 && vtx_idx_arr_refine[oc_id[4]] == 1) ||
				(vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[4]] == 1) ||
				(vtx_idx_arr_refine[oc_id[1]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1) ||
				(vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[6]] == 1) ||
				(vtx_idx_arr_refine[oc_id[3]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1);

	// detect one face
	my_bool_4 = (vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[1]] == 1 &&
				 vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1) ||
				(vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1 &&
				 vtx_idx_arr_refine[oc_id[6]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1) ||
				(vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1 &&
				 vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1) ||
				(vtx_idx_arr_refine[oc_id[1]] == 1 && vtx_idx_arr_refine[oc_id[2]] == 1 &&
				 vtx_idx_arr_refine[oc_id[5]] == 1 && vtx_idx_arr_refine[oc_id[6]] == 1) ||
				(vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[1]] == 1 &&
				 vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1) ||
				(vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1 &&
				 vtx_idx_arr_refine[oc_id[6]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1);

	if(num_id == 0) {
		add_hexa(geofrm, vtx);
	}
	else if(num_id == 1) {	// one point
		if(vtx_idx_arr_refine[oc_id[0]] == 1) {
			for(i = 0; i < 8; i++) vtx_temp[i] = vtx[i];
			for(i = 0; i < 6; i++) edge_id_new[i] = edge_id[i];
		}
		else if(vtx_idx_arr_refine[oc_id[1]] == 1) {
			vtx_temp[0] = vtx[1];	vtx_temp[1] = vtx[2];	vtx_temp[2] = vtx[3];	vtx_temp[3] = vtx[0];
			vtx_temp[4] = vtx[5];	vtx_temp[5] = vtx[6];	vtx_temp[6] = vtx[7];	vtx_temp[7] = vtx[4];
			edge_id_new[0] = edge_id[0];	edge_id_new[1] = edge_id[1];	edge_id_new[2] = edge_id[4];
			edge_id_new[3] = edge_id[5];	edge_id_new[4] = edge_id[3];	edge_id_new[5] = edge_id[2];
		}
		else if(vtx_idx_arr_refine[oc_id[2]] == 1) {
			vtx_temp[0] = vtx[2];	vtx_temp[1] = vtx[3];	vtx_temp[2] = vtx[0];	vtx_temp[3] = vtx[1];
			vtx_temp[4] = vtx[6];	vtx_temp[5] = vtx[7];	vtx_temp[6] = vtx[4];	vtx_temp[7] = vtx[5];
			edge_id_new[0] = edge_id[0];	edge_id_new[1] = edge_id[1];	edge_id_new[2] = edge_id[3];
			edge_id_new[3] = edge_id[2];	edge_id_new[4] = edge_id[5];	edge_id_new[5] = edge_id[4];
		}
		else if(vtx_idx_arr_refine[oc_id[3]] == 1) {
			vtx_temp[0] = vtx[3];	vtx_temp[1] = vtx[0];	vtx_temp[2] = vtx[1];	vtx_temp[3] = vtx[2];
			vtx_temp[4] = vtx[7];	vtx_temp[5] = vtx[4];	vtx_temp[6] = vtx[5];	vtx_temp[7] = vtx[6];
			edge_id_new[0] = edge_id[0];	edge_id_new[1] = edge_id[1];	edge_id_new[2] = edge_id[5];
			edge_id_new[3] = edge_id[4];	edge_id_new[4] = edge_id[2];	edge_id_new[5] = edge_id[3];
		}
		else if(vtx_idx_arr_refine[oc_id[4]] == 1) {
			vtx_temp[0] = vtx[4];	vtx_temp[1] = vtx[7];	vtx_temp[2] = vtx[6];	vtx_temp[3] = vtx[5];
			vtx_temp[4] = vtx[0];	vtx_temp[5] = vtx[3];	vtx_temp[6] = vtx[2];	vtx_temp[7] = vtx[1];
			edge_id_new[0] = edge_id[1];	edge_id_new[1] = edge_id[0];	edge_id_new[2] = edge_id[4];
			edge_id_new[3] = edge_id[5];	edge_id_new[4] = edge_id[2];	edge_id_new[5] = edge_id[3];
		}
		else if(vtx_idx_arr_refine[oc_id[5]] == 1) {
			vtx_temp[0] = vtx[5];	vtx_temp[1] = vtx[4];	vtx_temp[2] = vtx[7];	vtx_temp[3] = vtx[6];
			vtx_temp[4] = vtx[1];	vtx_temp[5] = vtx[0];	vtx_temp[6] = vtx[3];	vtx_temp[7] = vtx[2];
			edge_id_new[0] = edge_id[1];	edge_id_new[1] = edge_id[0];	edge_id_new[2] = edge_id[3];
			edge_id_new[3] = edge_id[2];	edge_id_new[4] = edge_id[4];	edge_id_new[5] = edge_id[5];
		}
		else if(vtx_idx_arr_refine[oc_id[6]] == 1) {
			vtx_temp[0] = vtx[6];	vtx_temp[1] = vtx[5];	vtx_temp[2] = vtx[4];	vtx_temp[3] = vtx[7];
			vtx_temp[4] = vtx[2];	vtx_temp[5] = vtx[1];	vtx_temp[6] = vtx[0];	vtx_temp[7] = vtx[3];
			edge_id_new[0] = edge_id[1];	edge_id_new[1] = edge_id[0];	edge_id_new[2] = edge_id[5];
			edge_id_new[3] = edge_id[4];	edge_id_new[4] = edge_id[3];	edge_id_new[5] = edge_id[2];
		}
		else { //if(vtx_idx_arr_refine[oc_id[7]] == 1) 
			vtx_temp[0] = vtx[7];	vtx_temp[1] = vtx[6];	vtx_temp[2] = vtx[5];	vtx_temp[3] = vtx[4];
			vtx_temp[4] = vtx[3];	vtx_temp[5] = vtx[2];	vtx_temp[6] = vtx[1];	vtx_temp[7] = vtx[0];
			edge_id_new[0] = edge_id[1];	edge_id_new[1] = edge_id[0];	edge_id_new[2] = edge_id[2];
			edge_id_new[3] = edge_id[3];	edge_id_new[4] = edge_id[5];	edge_id_new[5] = edge_id[4];
		}
		//add_hexa(geofrm, vtx);
		geofrm.AddVert_hexa_adaptive_2_1(vtx_temp, edge_id_new, vtx_new);
		for(i = 0; i < 7; i++) {
			if(geofrm.bound_sign[vtx_new[i]] != 1) continue;
			for(j = 0; j < 8; j++) {
				if(is_skipcell(oc_id[j])) continue;
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 8) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
		}
		add_hexa_adaptive_2_1(geofrm, vtx_temp, vtx_new);
	}
	else if(num_id == 2 && my_bool_2) {	// one edge
		if(vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[1]] == 1) {
			for(i = 0; i < 8; i++) vtx_temp[i] = vtx[i];
			for(i = 0; i < 6; i++) edge_id_new[i] = edge_id[i];
		}
		else if(vtx_idx_arr_refine[oc_id[1]] == 1 && vtx_idx_arr_refine[oc_id[2]] == 1) {
			vtx_temp[0] = vtx[1];	vtx_temp[1] = vtx[2];	vtx_temp[2] = vtx[3];	vtx_temp[3] = vtx[0];
			vtx_temp[4] = vtx[5];	vtx_temp[5] = vtx[6];	vtx_temp[6] = vtx[7];	vtx_temp[7] = vtx[4];
			edge_id_new[0] = edge_id[0];	edge_id_new[1] = edge_id[1];	edge_id_new[2] = edge_id[4];
			edge_id_new[3] = edge_id[5];	edge_id_new[4] = edge_id[3];	edge_id_new[5] = edge_id[2];
		}
		else if(vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1) {
			vtx_temp[0] = vtx[2];	vtx_temp[1] = vtx[3];	vtx_temp[2] = vtx[0];	vtx_temp[3] = vtx[1];
			vtx_temp[4] = vtx[6];	vtx_temp[5] = vtx[7];	vtx_temp[6] = vtx[4];	vtx_temp[7] = vtx[5];
			edge_id_new[0] = edge_id[0];	edge_id_new[1] = edge_id[1];	edge_id_new[2] = edge_id[3];
			edge_id_new[3] = edge_id[2];	edge_id_new[4] = edge_id[5];	edge_id_new[5] = edge_id[4];
		}
		else if(vtx_idx_arr_refine[oc_id[3]] == 1 && vtx_idx_arr_refine[oc_id[0]] == 1) {
			vtx_temp[0] = vtx[3];	vtx_temp[1] = vtx[0];	vtx_temp[2] = vtx[1];	vtx_temp[3] = vtx[2];
			vtx_temp[4] = vtx[7];	vtx_temp[5] = vtx[4];	vtx_temp[6] = vtx[5];	vtx_temp[7] = vtx[6];
			edge_id_new[0] = edge_id[0];	edge_id_new[1] = edge_id[1];	edge_id_new[2] = edge_id[5];
			edge_id_new[3] = edge_id[4];	edge_id_new[4] = edge_id[2];	edge_id_new[5] = edge_id[3];
		}
		else if(vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1) {
			vtx_temp[0] = vtx[4];	vtx_temp[1] = vtx[7];	vtx_temp[2] = vtx[6];	vtx_temp[3] = vtx[5];
			vtx_temp[4] = vtx[0];	vtx_temp[5] = vtx[3];	vtx_temp[6] = vtx[2];	vtx_temp[7] = vtx[1];
			edge_id_new[0] = edge_id[1];	edge_id_new[1] = edge_id[0];	edge_id_new[2] = edge_id[4];
			edge_id_new[3] = edge_id[5];	edge_id_new[4] = edge_id[2];	edge_id_new[5] = edge_id[3];
		}
		else if(vtx_idx_arr_refine[oc_id[5]] == 1 && vtx_idx_arr_refine[oc_id[4]] == 1) {
			vtx_temp[0] = vtx[5];	vtx_temp[1] = vtx[4];	vtx_temp[2] = vtx[7];	vtx_temp[3] = vtx[6];
			vtx_temp[4] = vtx[1];	vtx_temp[5] = vtx[0];	vtx_temp[6] = vtx[3];	vtx_temp[7] = vtx[2];
			edge_id_new[0] = edge_id[1];	edge_id_new[1] = edge_id[0];	edge_id_new[2] = edge_id[3];
			edge_id_new[3] = edge_id[2];	edge_id_new[4] = edge_id[4];	edge_id_new[5] = edge_id[5];
		}
		else if(vtx_idx_arr_refine[oc_id[6]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1) {
			vtx_temp[0] = vtx[6];	vtx_temp[1] = vtx[5];	vtx_temp[2] = vtx[4];	vtx_temp[3] = vtx[7];
			vtx_temp[4] = vtx[2];	vtx_temp[5] = vtx[1];	vtx_temp[6] = vtx[0];	vtx_temp[7] = vtx[3];
			edge_id_new[0] = edge_id[1];	edge_id_new[1] = edge_id[0];	edge_id_new[2] = edge_id[5];
			edge_id_new[3] = edge_id[4];	edge_id_new[4] = edge_id[3];	edge_id_new[5] = edge_id[2];
		}
		else if(vtx_idx_arr_refine[oc_id[7]] == 1 && vtx_idx_arr_refine[oc_id[6]] == 1) {
			vtx_temp[0] = vtx[7];	vtx_temp[1] = vtx[6];	vtx_temp[2] = vtx[5];	vtx_temp[3] = vtx[4];
			vtx_temp[4] = vtx[3];	vtx_temp[5] = vtx[2];	vtx_temp[6] = vtx[1];	vtx_temp[7] = vtx[0];
			edge_id_new[0] = edge_id[1];	edge_id_new[1] = edge_id[0];	edge_id_new[2] = edge_id[2];
			edge_id_new[3] = edge_id[3];	edge_id_new[4] = edge_id[5];	edge_id_new[5] = edge_id[4];
		}
		else if(vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[4]] == 1) {
			vtx_temp[0] = vtx[0];	vtx_temp[1] = vtx[4];	vtx_temp[2] = vtx[5];	vtx_temp[3] = vtx[1];
			vtx_temp[4] = vtx[3];	vtx_temp[5] = vtx[7];	vtx_temp[6] = vtx[6];	vtx_temp[7] = vtx[2];
			edge_id_new[0] = edge_id[4];	edge_id_new[1] = edge_id[5];	edge_id_new[2] = edge_id[0];
			edge_id_new[3] = edge_id[1];	edge_id_new[4] = edge_id[2];	edge_id_new[5] = edge_id[3];
		}
		else if(vtx_idx_arr_refine[oc_id[1]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1) {
			vtx_temp[0] = vtx[1];	vtx_temp[1] = vtx[5];	vtx_temp[2] = vtx[6];	vtx_temp[3] = vtx[2];
			vtx_temp[4] = vtx[0];	vtx_temp[5] = vtx[4];	vtx_temp[6] = vtx[7];	vtx_temp[7] = vtx[3];
			edge_id_new[0] = edge_id[3];	edge_id_new[1] = edge_id[2];	edge_id_new[2] = edge_id[0];
			edge_id_new[3] = edge_id[1];	edge_id_new[4] = edge_id[4];	edge_id_new[5] = edge_id[5];
		}
		else if(vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[6]] == 1) {
			vtx_temp[0] = vtx[2];	vtx_temp[1] = vtx[6];	vtx_temp[2] = vtx[7];	vtx_temp[3] = vtx[3];
			vtx_temp[4] = vtx[1];	vtx_temp[5] = vtx[5];	vtx_temp[6] = vtx[4];	vtx_temp[7] = vtx[0];
			edge_id_new[0] = edge_id[5];	edge_id_new[1] = edge_id[4];	edge_id_new[2] = edge_id[0];
			edge_id_new[3] = edge_id[1];	edge_id_new[4] = edge_id[3];	edge_id_new[5] = edge_id[2];
		}
		else { //if(vtx_idx_arr_refine[oc_id[3]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1)
			vtx_temp[0] = vtx[3];	vtx_temp[1] = vtx[7];	vtx_temp[2] = vtx[4];	vtx_temp[3] = vtx[0];
			vtx_temp[4] = vtx[2];	vtx_temp[5] = vtx[6];	vtx_temp[6] = vtx[5];	vtx_temp[7] = vtx[1];
			edge_id_new[0] = edge_id[2];	edge_id_new[1] = edge_id[3];	edge_id_new[2] = edge_id[0];
			edge_id_new[3] = edge_id[1];	edge_id_new[4] = edge_id[5];	edge_id_new[5] = edge_id[4];
		}
		//add_hexa(geofrm, vtx);
		geofrm.AddVert_hexa_adaptive_2_2(vtx_temp, edge_id_new, vtx_new);
		for(i = 0; i < 24; i++) {
			if(geofrm.bound_sign[vtx_new[i]] != 1 || i == 0 || i == 3 || i == 10 || i == 11) continue;
			for(j = 0; j < 8; j++) {
				if(is_skipcell(oc_id[j])) continue;
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 8) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
		}
		add_hexa_adaptive_2_2(geofrm, vtx_new);
	}
	else if(num_id == 4 && my_bool_4) {	// one face
		if( vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[1]] == 1 &&
			vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1) {
			for(i = 0; i < 8; i++) vtx_temp[i] = vtx[i];
			for(i = 0; i < 6; i++) edge_id_new[i] = edge_id[i];
		}
		else if(vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1 &&
				vtx_idx_arr_refine[oc_id[6]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1) {
			vtx_temp[0] = vtx[4];	vtx_temp[1] = vtx[7];	vtx_temp[2] = vtx[6];	vtx_temp[3] = vtx[5];
			vtx_temp[4] = vtx[0];	vtx_temp[5] = vtx[3];	vtx_temp[6] = vtx[2];	vtx_temp[7] = vtx[1];
			edge_id_new[0] = edge_id[1];	edge_id_new[1] = edge_id[0];	edge_id_new[2] = edge_id[4];
			edge_id_new[3] = edge_id[5];	edge_id_new[4] = edge_id[2];	edge_id_new[5] = edge_id[3];
		}
		else if(vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1 &&
				vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1) {
			vtx_temp[0] = vtx[3];	vtx_temp[1] = vtx[7];	vtx_temp[2] = vtx[4];	vtx_temp[3] = vtx[0];
			vtx_temp[4] = vtx[2];	vtx_temp[5] = vtx[6];	vtx_temp[6] = vtx[5];	vtx_temp[7] = vtx[1];
			edge_id_new[0] = edge_id[2];	edge_id_new[1] = edge_id[3];	edge_id_new[2] = edge_id[0];
			edge_id_new[3] = edge_id[1];	edge_id_new[4] = edge_id[5];	edge_id_new[5] = edge_id[4];
		}
		else if(vtx_idx_arr_refine[oc_id[1]] == 1 && vtx_idx_arr_refine[oc_id[2]] == 1 &&
				vtx_idx_arr_refine[oc_id[5]] == 1 && vtx_idx_arr_refine[oc_id[6]] == 1) {
			vtx_temp[0] = vtx[1];	vtx_temp[1] = vtx[5];	vtx_temp[2] = vtx[6];	vtx_temp[3] = vtx[2];
			vtx_temp[4] = vtx[0];	vtx_temp[5] = vtx[4];	vtx_temp[6] = vtx[7];	vtx_temp[7] = vtx[3];
			edge_id_new[0] = edge_id[3];	edge_id_new[1] = edge_id[2];	edge_id_new[2] = edge_id[0];
			edge_id_new[3] = edge_id[1];	edge_id_new[4] = edge_id[4];	edge_id_new[5] = edge_id[5];
		}
		else if(vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[1]] == 1 &&
				vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1) {
			vtx_temp[0] = vtx[0];	vtx_temp[1] = vtx[4];	vtx_temp[2] = vtx[5];	vtx_temp[3] = vtx[1];
			vtx_temp[4] = vtx[3];	vtx_temp[5] = vtx[7];	vtx_temp[6] = vtx[6];	vtx_temp[7] = vtx[2];
			edge_id_new[0] = edge_id[4];	edge_id_new[1] = edge_id[5];	edge_id_new[2] = edge_id[0];
			edge_id_new[3] = edge_id[1];	edge_id_new[4] = edge_id[2];	edge_id_new[5] = edge_id[3];
		}
		else {	//if(vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1 &&
				//   vtx_idx_arr_refine[oc_id[6]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1);
			vtx_temp[0] = vtx[2];	vtx_temp[1] = vtx[6];	vtx_temp[2] = vtx[7];	vtx_temp[3] = vtx[3];
			vtx_temp[4] = vtx[1];	vtx_temp[5] = vtx[5];	vtx_temp[6] = vtx[4];	vtx_temp[7] = vtx[0];
			edge_id_new[0] = edge_id[5];	edge_id_new[1] = edge_id[4];	edge_id_new[2] = edge_id[0];
			edge_id_new[3] = edge_id[1];	edge_id_new[4] = edge_id[3];	edge_id_new[5] = edge_id[2];
		}
		//add_hexa(geofrm, vtx);
		geofrm.AddVert_hexa_adaptive_2_4(vtx_temp, edge_id_new, vtx_new);
		for(i = 0; i < 44; i++) {
			if(geofrm.bound_sign[vtx_new[i]] != 1 || i == 0 || i == 3 || i == 12 || i == 15) continue;
			for(j = 0; j < 8; j++) {
				if(is_skipcell(oc_id[j])) continue;
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 8) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
		}
		add_hexa_adaptive_2_4(geofrm, vtx_new);
	} 
	else {	
		assert(num_id == 8);
		//add_hexa(geofrm, vtx);
		geofrm.AddVert_hexa_adaptive_2(vtx, edge_id, vtx_new);
		for(i = 0; i < 64; i++) {
			if(geofrm.bound_sign[vtx_new[i]] != 1 || i%48 == 0 || i%48 == 3 || i%48 == 12 || i%48 == 15) continue;
			for(j = 0; j < 8; j++) {
				if(is_skipcell(oc_id[j])) continue;
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 8) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
		}
		add_hexa_adaptive_2(geofrm, vtx_new);
	}

}

void Octree::hexa_adaptive_1_top(geoframe& geofrm, unsigned int* vtx, unsigned int* vtx_new) {

	unsigned int temp[8];

	add_hexa(geofrm, vtx_new);

	temp[0] = vtx[0]; temp[1] = vtx[1]; temp[2] = vtx_new[1]; temp[3] = vtx_new[0];
	temp[4] = vtx[4]; temp[5] = vtx[5]; temp[6] = vtx_new[5]; temp[7] = vtx_new[4];
	add_hexa(geofrm, temp);

	temp[0] = vtx[1]; temp[1] = vtx[2]; temp[2] = vtx_new[2]; temp[3] = vtx_new[1];
	temp[4] = vtx[5]; temp[5] = vtx[6]; temp[6] = vtx_new[6]; temp[7] = vtx_new[5];
	add_hexa(geofrm, temp);

	temp[0] = vtx[2]; temp[1] = vtx[3]; temp[2] = vtx_new[3]; temp[3] = vtx_new[2];
	temp[4] = vtx[6]; temp[5] = vtx[7]; temp[6] = vtx_new[7]; temp[7] = vtx_new[6];
	add_hexa(geofrm, temp);

	temp[0] = vtx[3]; temp[1] = vtx[0]; temp[2] = vtx_new[0]; temp[3] = vtx_new[3];
	temp[4] = vtx[7]; temp[5] = vtx[4]; temp[6] = vtx_new[4]; temp[7] = vtx_new[7];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[4];	temp[1] = vtx_new[5];	temp[2] = vtx_new[6];	temp[3] = vtx_new[7];
	temp[4] = vtx[4];		temp[5] = vtx[5];		temp[6] = vtx[6];		temp[7] = vtx[7];
	add_hexa(geofrm, temp);

}

void Octree::add_hexa_adaptive_2(geoframe& geofrm, unsigned int* vtx_new) {

	int i, j, k;
	unsigned int temp[8];

	for(k = 0; k < 3; k++) {
		for(j = 0; j < 3; j++) {
			for(i = 0; i < 3; i++) {
				temp[0] = vtx_new[i+4*j+16*k];		temp[1] = vtx_new[i+4*j+16*k+1];
				temp[2] = vtx_new[i+4*j+16*k+5];	temp[3] = vtx_new[i+4*j+16*k+4];
				temp[4] = vtx_new[i+4*j+16*k+16];	temp[5] = vtx_new[i+4*j+16*k+17];
				temp[6] = vtx_new[i+4*j+16*k+21];	temp[7] = vtx_new[i+4*j+16*k+20];
				add_hexa(geofrm, temp);
			}
		}
	}
}

// one vertex needs to be refined
void Octree::add_hexa_adaptive_2_1(geoframe& geofrm, unsigned int* vtx, unsigned int* vtx_temp) {
	
	int i;
	unsigned int temp[8];

	temp[0] = vtx[0];
	for(i = 0; i < 7; i++) temp[i+1] = vtx_temp[i];
	add_hexa(geofrm, temp);

	for(i = 0; i < 4; i++) temp[i] = vtx_temp[i+3];
	for(i = 4; i < 8; i++) temp[i] = vtx[i];
	add_hexa(geofrm, temp);

	for(i = 0; i < 8; i++) temp[i] = vtx[i];
	temp[0] = vtx_temp[0];	temp[3] = vtx_temp[1];
	temp[4] = vtx_temp[4];	temp[7] = vtx_temp[5];
	add_hexa(geofrm, temp);

	for(i = 0; i < 8; i++) temp[i] = vtx[i];
	temp[0] = vtx_temp[2];	temp[1] = vtx_temp[1];
	temp[4] = vtx_temp[6];	temp[5] = vtx_temp[5];
	add_hexa(geofrm, temp);
	
}

// one edge needs to be refined
void Octree::add_hexa_adaptive_2_2(geoframe& geofrm, unsigned int* vtx_new) {
	
	int i;
	unsigned int temp[8];

	for(i = 0; i < 3; i++) {
		temp[0] = vtx_new[i+0];		temp[1] = vtx_new[i+1];
		temp[2] = vtx_new[i+5];		temp[3] = vtx_new[i+4];
		temp[4] = vtx_new[i+12];	temp[5] = vtx_new[i+13];
		temp[6] = vtx_new[i+17];	temp[7] = vtx_new[i+16];
		add_hexa(geofrm, temp);
	}

	temp[0] = vtx_new[13];		temp[1] = vtx_new[14];
	temp[2] = vtx_new[18];		temp[3] = vtx_new[17];
	temp[4] = vtx_new[20];		temp[5] = vtx_new[21];
	temp[6] = vtx_new[23];		temp[7] = vtx_new[22];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[12];		temp[1] = vtx_new[13];
	temp[2] = vtx_new[17];		temp[3] = vtx_new[16];
	temp[4] = vtx_new[24];		temp[5] = vtx_new[20];
	temp[6] = vtx_new[22];		temp[7] = vtx_new[26];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[14];		temp[1] = vtx_new[15];
	temp[2] = vtx_new[19];		temp[3] = vtx_new[18];
	temp[4] = vtx_new[21];		temp[5] = vtx_new[25];
	temp[6] = vtx_new[27];		temp[7] = vtx_new[23];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[20];		temp[1] = vtx_new[21];
	temp[2] = vtx_new[23];		temp[3] = vtx_new[22];
	temp[4] = vtx_new[24];		temp[5] = vtx_new[25];
	temp[6] = vtx_new[27];		temp[7] = vtx_new[26];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[5];		temp[1] = vtx_new[6];
	temp[2] = vtx_new[9];		temp[3] = vtx_new[8];
	temp[4] = vtx_new[17];		temp[5] = vtx_new[18];
	temp[6] = vtx_new[23];		temp[7] = vtx_new[22];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[4];		temp[1] = vtx_new[5];
	temp[2] = vtx_new[8];		temp[3] = vtx_new[10];
	temp[4] = vtx_new[16];		temp[5] = vtx_new[17];
	temp[6] = vtx_new[22];		temp[7] = vtx_new[26];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[6];		temp[1] = vtx_new[7];
	temp[2] = vtx_new[11];		temp[3] = vtx_new[9];
	temp[4] = vtx_new[18];		temp[5] = vtx_new[19];
	temp[6] = vtx_new[27];		temp[7] = vtx_new[23];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[8];		temp[1] = vtx_new[9];
	temp[2] = vtx_new[11];		temp[3] = vtx_new[10];
	temp[4] = vtx_new[22];		temp[5] = vtx_new[23];
	temp[6] = vtx_new[27];		temp[7] = vtx_new[26];
	add_hexa(geofrm, temp);
}

// one face needs to be refined
void Octree::add_hexa_adaptive_2_4(geoframe& geofrm, unsigned int* vtx_new) {
	
	int i, j;
	unsigned int temp[8];

	for(j = 0; j < 3; j++) {
		for(i = 0; i < 3; i++) {
			temp[0] = vtx_new[i+4*j+0];		temp[1] = vtx_new[i+4*j+1];
			temp[2] = vtx_new[i+4*j+5];		temp[3] = vtx_new[i+4*j+4];
			temp[4] = vtx_new[i+4*j+16];	temp[5] = vtx_new[i+4*j+17];
			temp[6] = vtx_new[i+4*j+21];	temp[7] = vtx_new[i+4*j+20];
			add_hexa(geofrm, temp);
		}
	}

	// middle - 4
	temp[0] = vtx_new[21];		temp[1] = vtx_new[22];
	temp[2] = vtx_new[26];		temp[3] = vtx_new[25];
	temp[4] = vtx_new[32];		temp[5] = vtx_new[33];
	temp[6] = vtx_new[35];		temp[7] = vtx_new[34];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[32];		temp[1] = vtx_new[33];
	temp[2] = vtx_new[35];		temp[3] = vtx_new[34];
	temp[4] = vtx_new[36];		temp[5] = vtx_new[37];
	temp[6] = vtx_new[43];		temp[7] = vtx_new[42];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[17];		temp[1] = vtx_new[18];
	temp[2] = vtx_new[22];		temp[3] = vtx_new[21];
	temp[4] = vtx_new[36];		temp[5] = vtx_new[37];
	temp[6] = vtx_new[33];		temp[7] = vtx_new[32];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[25];		temp[1] = vtx_new[26];
	temp[2] = vtx_new[30];		temp[3] = vtx_new[29];
	temp[4] = vtx_new[34];		temp[5] = vtx_new[35];
	temp[6] = vtx_new[43];		temp[7] = vtx_new[42];
	add_hexa(geofrm, temp);

	// bottom - 1
	temp[0] = vtx_new[36];		temp[1] = vtx_new[37];
	temp[2] = vtx_new[43];		temp[3] = vtx_new[42];
	temp[4] = vtx_new[44];		temp[5] = vtx_new[45];
	temp[6] = vtx_new[47];		temp[7] = vtx_new[46];
	add_hexa(geofrm, temp);

	// left - 4
	temp[0] = vtx_new[16];		temp[1] = vtx_new[17];
	temp[2] = vtx_new[21];		temp[3] = vtx_new[20];
	temp[4] = vtx_new[44];		temp[5] = vtx_new[36];
	temp[6] = vtx_new[32];		temp[7] = vtx_new[38];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[20];		temp[1] = vtx_new[21];
	temp[2] = vtx_new[25];		temp[3] = vtx_new[24];
	temp[4] = vtx_new[38];		temp[5] = vtx_new[32];
	temp[6] = vtx_new[34];		temp[7] = vtx_new[40];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[24];		temp[1] = vtx_new[25];
	temp[2] = vtx_new[29];		temp[3] = vtx_new[28];
	temp[4] = vtx_new[40];		temp[5] = vtx_new[34];
	temp[6] = vtx_new[42];		temp[7] = vtx_new[46];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[38];		temp[1] = vtx_new[32];
	temp[2] = vtx_new[34];		temp[3] = vtx_new[40];
	temp[4] = vtx_new[44];		temp[5] = vtx_new[36];
	temp[6] = vtx_new[42];		temp[7] = vtx_new[46];
	add_hexa(geofrm, temp);

	// right - 4
	temp[0] = vtx_new[18];		temp[1] = vtx_new[19];
	temp[2] = vtx_new[23];		temp[3] = vtx_new[22];
	temp[4] = vtx_new[37];		temp[5] = vtx_new[45];
	temp[6] = vtx_new[39];		temp[7] = vtx_new[33];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[22];		temp[1] = vtx_new[23];
	temp[2] = vtx_new[27];		temp[3] = vtx_new[26];
	temp[4] = vtx_new[33];		temp[5] = vtx_new[39];
	temp[6] = vtx_new[41];		temp[7] = vtx_new[35];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[26];		temp[1] = vtx_new[27];
	temp[2] = vtx_new[31];		temp[3] = vtx_new[30];
	temp[4] = vtx_new[35];		temp[5] = vtx_new[41];
	temp[6] = vtx_new[47];		temp[7] = vtx_new[43];
	add_hexa(geofrm, temp);

	temp[0] = vtx_new[33];		temp[1] = vtx_new[39];
	temp[2] = vtx_new[41];		temp[3] = vtx_new[35];
	temp[4] = vtx_new[37];		temp[5] = vtx_new[45];
	temp[6] = vtx_new[47];		temp[7] = vtx_new[43];
	add_hexa(geofrm, temp);
}

void Octree::vflag_clear()
{
  //memset(vbit, 0, octcell_num*4/8);
  std::fill(vbit.begin(),vbit.end(),0);
}

int Octree::is_vflag_on(int x, int y, int z, int level, int v)
{
	int idx;
	switch (v) {
	case 0 : 
		idx = xyz2octcell(x,y,z,level);
		break;
	case 1 :
		idx = xyz2octcell(x+1,y,z,level);
		break;
	case 2 :
		idx = xyz2octcell(x+1,y,z+1,level);
		break;
	case 3 :
		idx = xyz2octcell(x,y,z+1,level);
		break;
	case 4 : 
		idx = xyz2octcell(x,y+1,z,level);
		break;
	case 5 :
		idx = xyz2octcell(x+1,y+1,z,level);
		break;
	case 6 :
		idx = xyz2octcell(x+1,y+1,z+1,level);
		break;
	case 7 :
		idx = xyz2octcell(x,y+1,z+1,level);
		break;
	}
	
	if (vbit[idx/8] & (1 << (idx%8))) return 1;
	else return 0;
	
}

void Octree::vflag_on(int x, int y, int z, int level, int v)
{
	int idx;
	switch (v) {
	case 0 : 
		idx = xyz2octcell(x,y,z,level);
		break;
	case 1 :
		idx = xyz2octcell(x+1,y,z,level);
		break;
	case 2 :
		idx = xyz2octcell(x+1,y,z+1,level);
		break;
	case 3 :
		idx = xyz2octcell(x,y,z+1,level);
		break;
	case 4 : 
		idx = xyz2octcell(x,y+1,z,level);
		break;
	case 5 :
		idx = xyz2octcell(x+1,y+1,z,level);
		break;
	case 6 :
		idx = xyz2octcell(x+1,y+1,z+1,level);
		break;
	case 7 :
		idx = xyz2octcell(x,y+1,z+1,level);
		break;
	}
	
	vbit[idx/8] |= (1 << (idx%8));
	
}


int Octree::is_min_vertex(int oc_id, int v_id, unsigned int* vtx, geoframe& geofrm)
{
	int x, y, z, level, i;

	//assert(! is_skipcell(oc_id));

	level = get_level(oc_id);
	octcell2xyz(oc_id, x, y, z, level);
	
	for(i = 0; i < 8; i++) vtx[i] = -1;

	switch (v_id) {
		case 0 : 
			if (is_refined(x-1,y-1,z-1,level) || is_refined(x,  y-1,z-1,level) ||
				is_refined(x,  y,  z-1,level) || is_refined(x-1,y,  z-1,level) ||
				is_refined(x-1,y-1,z,  level) || is_refined(x,  y-1,z,  level) ||
				is_refined(x,  y,  z,  level) || is_refined(x-1,y,  z,  level)) return 0;
			vtx[0] = min_vtx_hexa(x-1, y-1, z-1, level, geofrm);
			vtx[1] = min_vtx_hexa(x,   y-1, z-1, level, geofrm);
			vtx[2] = min_vtx_hexa(x,   y,   z-1, level, geofrm);
			vtx[3] = min_vtx_hexa(x-1, y,   z-1, level, geofrm);
			vtx[4] = min_vtx_hexa(x-1, y-1, z,   level, geofrm);
			vtx[5] = min_vtx_hexa(x,   y-1, z,   level, geofrm);
			vtx[6] = min_vtx_hexa(x,   y,   z,   level, geofrm);
			vtx[7] = min_vtx_hexa(x-1, y,   z,   level, geofrm);
			break;

		case 1 : 
			if (is_refined(x,  y-1,z-1,level) || is_refined(x+1,y-1,z-1,level) ||
				is_refined(x+1,y,  z-1,level) || is_refined(x,  y,  z-1,level) ||
				is_refined(x,  y-1,z,  level) || is_refined(x+1,y-1,z,  level) ||
				is_refined(x+1,y,  z,  level) || is_refined(x,  y,  z,  level)) return 0;
			vtx[0] = min_vtx_hexa(x,   y-1, z-1, level, geofrm);
			vtx[1] = min_vtx_hexa(x+1, y-1, z-1, level, geofrm);
			vtx[2] = min_vtx_hexa(x+1, y,   z-1, level, geofrm);
			vtx[3] = min_vtx_hexa(x,   y,   z-1, level, geofrm);
			vtx[4] = min_vtx_hexa(x,   y-1, z,   level, geofrm);
			vtx[5] = min_vtx_hexa(x+1, y-1, z,   level, geofrm);
			vtx[6] = min_vtx_hexa(x+1, y,   z,   level, geofrm);
			vtx[7] = min_vtx_hexa(x,   y,   z,   level, geofrm);
			break;

		case 2 : 
			if (is_refined(x,  y-1,z,  level) || is_refined(x+1,y-1,z,  level) ||
				is_refined(x+1,y,  z,  level) || is_refined(x,  y,  z,  level) ||
				is_refined(x,  y-1,z+1,level) || is_refined(x+1,y-1,z+1,level) ||
				is_refined(x+1,y,  z+1,level) || is_refined(x,  y,  z+1,level)) return 0;
			vtx[0] = min_vtx_hexa(x,   y-1, z,   level, geofrm);
			vtx[1] = min_vtx_hexa(x+1, y-1, z,   level, geofrm);
			vtx[2] = min_vtx_hexa(x+1, y,   z,   level, geofrm);
			vtx[3] = min_vtx_hexa(x,   y,   z,   level, geofrm);
			vtx[4] = min_vtx_hexa(x,   y-1, z+1, level, geofrm);
			vtx[5] = min_vtx_hexa(x+1, y-1, z+1, level, geofrm);
			vtx[6] = min_vtx_hexa(x+1, y,   z+1, level, geofrm);
			vtx[7] = min_vtx_hexa(x,   y,   z+1, level, geofrm);
			break;

		case 3 : 
			if (is_refined(x-1,y-1,z,  level) || is_refined(x,  y-1,z,  level) ||
				is_refined(x,  y,  z,  level) || is_refined(x-1,y,  z,  level) ||
				is_refined(x-1,y-1,z+1,level) || is_refined(x,  y-1,z+1,level) ||
				is_refined(x,  y,  z+1,level) || is_refined(x-1,y,  z+1,level)) return 0;
			vtx[0] = min_vtx_hexa(x-1, y-1, z,   level, geofrm);
			vtx[1] = min_vtx_hexa(x,   y-1, z,   level, geofrm);
			vtx[2] = min_vtx_hexa(x,   y,   z,   level, geofrm);
			vtx[3] = min_vtx_hexa(x-1, y,   z,   level, geofrm);
			vtx[4] = min_vtx_hexa(x-1, y-1, z+1, level, geofrm);
			vtx[5] = min_vtx_hexa(x,   y-1, z+1, level, geofrm);
			vtx[6] = min_vtx_hexa(x,   y,   z+1, level, geofrm);
			vtx[7] = min_vtx_hexa(x-1, y,   z+1, level, geofrm);
			break;

		case 4 : 
			if (is_refined(x-1,y,  z-1,level) || is_refined(x,  y,  z-1,level) ||
				is_refined(x,  y+1,z-1,level) || is_refined(x-1,y+1,z-1,level) ||
				is_refined(x-1,y,  z,  level) || is_refined(x,  y,  z,  level) ||
				is_refined(x,  y+1,z,  level) || is_refined(x-1,y+1,z,  level)) return 0;
			vtx[0] = min_vtx_hexa(x-1, y,   z-1, level, geofrm);
			vtx[1] = min_vtx_hexa(x,   y,   z-1, level, geofrm);
			vtx[2] = min_vtx_hexa(x,   y+1, z-1, level, geofrm);
			vtx[3] = min_vtx_hexa(x-1, y+1, z-1, level, geofrm);
			vtx[4] = min_vtx_hexa(x-1, y,   z,   level, geofrm);
			vtx[5] = min_vtx_hexa(x,   y,   z,   level, geofrm);
			vtx[6] = min_vtx_hexa(x,   y+1, z,   level, geofrm);
			vtx[7] = min_vtx_hexa(x-1, y+1, z,   level, geofrm);
			break;

		case 5 :
			if (is_refined(x,  y,  z-1,level) || is_refined(x+1,y,  z-1,level) ||
				is_refined(x+1,y+1,z-1,level) || is_refined(x,  y+1,z-1,level) ||
				is_refined(x,  y,  z,  level) || is_refined(x+1,y,  z,  level) ||
				is_refined(x+1,y+1,z,  level) || is_refined(x,  y+1,z,  level)) return 0;
			vtx[0] = min_vtx_hexa(x,   y,   z-1, level, geofrm);
			vtx[1] = min_vtx_hexa(x+1, y,   z-1, level, geofrm);
			vtx[2] = min_vtx_hexa(x+1, y+1, z-1, level, geofrm);
			vtx[3] = min_vtx_hexa(x,   y+1, z-1, level, geofrm);
			vtx[4] = min_vtx_hexa(x,   y,   z,   level, geofrm);
			vtx[5] = min_vtx_hexa(x+1, y,   z,   level, geofrm);
			vtx[6] = min_vtx_hexa(x+1, y+1, z,   level, geofrm);
			vtx[7] = min_vtx_hexa(x,   y+1, z,   level, geofrm);
			break;

		case 6 :
			if (is_refined(x,  y,  z,  level) || is_refined(x+1,y,  z,  level) ||
				is_refined(x+1,y+1,z,  level) || is_refined(x,  y+1,z,  level) ||
				is_refined(x,  y,  z+1,level) || is_refined(x+1,y,  z+1,level) ||
				is_refined(x+1,y+1,z+1,level) || is_refined(x,  y+1,z+1,level)) return 0;
			vtx[0] = min_vtx_hexa(x,   y,   z,   level, geofrm);
			vtx[1] = min_vtx_hexa(x+1, y,   z,   level, geofrm);
			vtx[2] = min_vtx_hexa(x+1, y+1, z,   level, geofrm);
			vtx[3] = min_vtx_hexa(x,   y+1, z,   level, geofrm);
			vtx[4] = min_vtx_hexa(x,   y,   z+1, level, geofrm);
			vtx[5] = min_vtx_hexa(x+1, y,   z+1, level, geofrm);
			vtx[6] = min_vtx_hexa(x+1, y+1, z+1, level, geofrm);
			vtx[7] = min_vtx_hexa(x,   y+1, z+1, level, geofrm);
			break;

		case 7 :
			if (is_refined(x-1,y,  z,  level) || is_refined(x,  y,  z,  level) ||
				is_refined(x,  y+1,z,  level) || is_refined(x-1,y+1,z,  level) ||
				is_refined(x-1,y,  z+1,level) || is_refined(x,  y,  z+1,level) ||
				is_refined(x,  y+1,z+1,level) || is_refined(x-1,y+1,z+1,level)) return 0;
			vtx[0] = min_vtx_hexa(x-1, y,   z,   level, geofrm);
			vtx[1] = min_vtx_hexa(x,   y,   z,   level, geofrm);
			vtx[2] = min_vtx_hexa(x,   y+1, z,   level, geofrm);
			vtx[3] = min_vtx_hexa(x-1, y+1, z,   level, geofrm);
			vtx[4] = min_vtx_hexa(x-1, y,   z+1, level, geofrm);
			vtx[5] = min_vtx_hexa(x,   y,   z+1, level, geofrm);
			vtx[6] = min_vtx_hexa(x,   y+1, z+1, level, geofrm);
			vtx[7] = min_vtx_hexa(x-1, y+1, z+1, level, geofrm);
			break;
	}

	if( vtx[0] >= 0 && vtx[1] >= 0 && vtx[2] >= 0 && vtx[3] >= 0 &&
		vtx[4] >= 0 && vtx[5] >= 0 && vtx[6] >= 0 && vtx[7] >= 0 ) return 1;
	else	return 0;

}


int Octree::min_vtx_hexa(int x, int y, int z, int level, geoframe& geofrm)
{
	int oc_id, cell_size, tx, ty, tz, vert;
	unsigned int center;
	float vtx[3], norm[3];

	tx = x;	ty = y;	tz = z;

	assert( tx >= 0 && ty >= 0 && tz >= 0 );
	//assert( !is_refined(tx,ty,tz,level) );

	while ( level == 0 || !is_refined(tx/2 , ty/2 , tz/2 , level-1 )) {
		tx /= 2;
		ty /= 2;
		tz /= 2;
		level--;
	}

	x = tx;	y = ty;	z = tz;

	oc_id = xyz2octcell(x, y, z, level);
	cell_size = (dim[0]-1)/(1<<level);

	if(x < 0 || y < 0 || z < 0 || x > dim[0]-1 || y > dim[1]-1 || z > dim[2]-1)
		return -1;
	else if(minmax[oc_id].max <= iso_val) {
		if ((center = vtx_idx_arr[oc_id]) == -1) {
			add_middle_vertex(x, y, z, 0.5, 0.5, 0.5, cell_size, center, geofrm);
			vtx_idx_arr[oc_id] = center;
			return center;
		}
		else
			return center;
	}
	else {
		get_vtx(x, y, z, level, vtx);
		//getVertGrad(x*cell_size, y*cell_size, z*cell_size, norm);
		get_VtxNorm(vtx, norm);

		if ((vert = vtx_idx_arr[oc_id]) == -1) {
			vert = geofrm.AddVert(vtx, norm);
			geofrm.AddBound(vert, 1);
			vtx_idx_arr[oc_id] = vert;
			return vert;
		}
		else
			return vert;
	}

}


void Octree::add_hexa(geoframe& geofrm, unsigned int* vtx) {

	unsigned int my_vtx[4];

	my_vtx[0] = vtx[0];	my_vtx[1] = vtx[3];
	my_vtx[2] = vtx[7];	my_vtx[3] = vtx[4];
	geofrm.AddQuad(my_vtx, 4);

	my_vtx[0] = vtx[2];	my_vtx[1] = vtx[1];
	my_vtx[2] = vtx[5];	my_vtx[3] = vtx[6];
	geofrm.AddQuad(my_vtx, 4);

	my_vtx[0] = vtx[0];	my_vtx[1] = vtx[4];
	my_vtx[2] = vtx[5];	my_vtx[3] = vtx[1];
	geofrm.AddQuad(my_vtx, 4);

	my_vtx[0] = vtx[3];	my_vtx[1] = vtx[2];
	my_vtx[2] = vtx[6];	my_vtx[3] = vtx[7];
	geofrm.AddQuad(my_vtx, 4);

	my_vtx[0] = vtx[0];	my_vtx[1] = vtx[1];
	my_vtx[2] = vtx[2];	my_vtx[3] = vtx[3];
	geofrm.AddQuad(my_vtx, 4);

	my_vtx[0] = vtx[4];	my_vtx[1] = vtx[7];
	my_vtx[2] = vtx[6];	my_vtx[3] = vtx[5];
	geofrm.AddQuad(my_vtx, 4);

	geofrm.numhexas++;
}

void Octree::find_oc_id_hexa(int x, int y, int z, int level, int j, int* oc_id) {

	oc_id[0] = xyz2octcell(x, y, z, level);

	switch (j) {
		case 0 : 
			oc_id[0] = xyz2octcell(x-1, y-1, z-1, level);
			oc_id[1] = xyz2octcell(x,   y-1, z-1, level);
			oc_id[2] = xyz2octcell(x,   y,   z-1, level);
			oc_id[3] = xyz2octcell(x-1, y,   z-1, level);
			oc_id[4] = xyz2octcell(x-1, y-1, z,   level);
			oc_id[5] = xyz2octcell(x,   y-1, z,   level);
			oc_id[6] = xyz2octcell(x,   y,   z,   level);
			oc_id[7] = xyz2octcell(x-1, y,   z,   level);
			break;

		case 1 : 
			oc_id[0] = xyz2octcell(x,   y-1, z-1, level);
			oc_id[1] = xyz2octcell(x+1, y-1, z-1, level);
			oc_id[2] = xyz2octcell(x+1, y,   z-1, level);
			oc_id[3] = xyz2octcell(x,   y,   z-1, level);
			oc_id[4] = xyz2octcell(x,   y-1, z,   level);
			oc_id[5] = xyz2octcell(x+1, y-1, z,   level);
			oc_id[6] = xyz2octcell(x+1, y,   z,   level);
			oc_id[7] = xyz2octcell(x,   y,   z,   level);
			break;

		case 2 : 
			oc_id[0] = xyz2octcell(x,   y-1, z,   level);
			oc_id[1] = xyz2octcell(x+1, y-1, z,   level);
			oc_id[2] = xyz2octcell(x+1, y,   z,   level);
			oc_id[3] = xyz2octcell(x,   y,   z,   level);
			oc_id[4] = xyz2octcell(x,   y-1, z+1, level);
			oc_id[5] = xyz2octcell(x+1, y-1, z+1, level);
			oc_id[6] = xyz2octcell(x+1, y,   z+1, level);
			oc_id[7] = xyz2octcell(x,   y,   z+1, level);
			break;

		case 3 : 
			oc_id[0] = xyz2octcell(x-1, y-1, z,   level);
			oc_id[1] = xyz2octcell(x,   y-1, z,   level);
			oc_id[2] = xyz2octcell(x,   y,   z,   level);
			oc_id[3] = xyz2octcell(x-1, y,   z,   level);
			oc_id[4] = xyz2octcell(x-1, y-1, z+1, level);
			oc_id[5] = xyz2octcell(x,   y-1, z+1, level);
			oc_id[6] = xyz2octcell(x,   y,   z+1, level);
			oc_id[7] = xyz2octcell(x-1, y,   z+1, level);
			break;

		case 4 : 
			oc_id[0] = xyz2octcell(x-1, y,   z-1, level);
			oc_id[1] = xyz2octcell(x,   y,   z-1, level);
			oc_id[2] = xyz2octcell(x,   y+1, z-1, level);
			oc_id[3] = xyz2octcell(x-1, y+1, z-1, level);
			oc_id[4] = xyz2octcell(x-1, y,   z,   level);
			oc_id[5] = xyz2octcell(x,   y,   z,   level);
			oc_id[6] = xyz2octcell(x,   y+1, z,   level);
			oc_id[7] = xyz2octcell(x-1, y+1, z,   level);
			break;

		case 5 :
			oc_id[0] = xyz2octcell(x,   y,   z-1, level);
			oc_id[1] = xyz2octcell(x+1, y,   z-1, level);
			oc_id[2] = xyz2octcell(x+1, y+1, z-1, level);
			oc_id[3] = xyz2octcell(x,   y+1, z-1, level);
			oc_id[4] = xyz2octcell(x,   y,   z,   level);
			oc_id[5] = xyz2octcell(x+1, y,   z,   level);
			oc_id[6] = xyz2octcell(x+1, y+1, z,   level);
			oc_id[7] = xyz2octcell(x,   y+1, z,   level);
			break;

		case 6 :
			oc_id[0] = xyz2octcell(x,   y,   z,   level);
			oc_id[1] = xyz2octcell(x+1, y,   z,   level);
			oc_id[2] = xyz2octcell(x+1, y+1, z,   level);
			oc_id[3] = xyz2octcell(x,   y+1, z,   level);
			oc_id[4] = xyz2octcell(x,   y,   z+1, level);
			oc_id[5] = xyz2octcell(x+1, y,   z+1, level);
			oc_id[6] = xyz2octcell(x+1, y+1, z+1, level);
			oc_id[7] = xyz2octcell(x,   y+1, z+1, level);
			break;

		case 7 :
			oc_id[0] = xyz2octcell(x-1, y,   z,   level);
			oc_id[1] = xyz2octcell(x,   y,   z,   level);
			oc_id[2] = xyz2octcell(x,   y+1, z,   level);
			oc_id[3] = xyz2octcell(x-1, y+1, z,   level);
			oc_id[4] = xyz2octcell(x-1, y,   z+1, level);
			oc_id[5] = xyz2octcell(x,   y,   z+1, level);
			oc_id[6] = xyz2octcell(x,   y+1, z+1, level);
			oc_id[7] = xyz2octcell(x-1, y+1, z+1, level);
			break;
	}

}

// for each interior vertex, check if which one of its 6 edges is sign change edge
// decide which quad is on the boundary surface
void Octree::find_edge_id_hexa(int x, int y, int z, int cell_size, int j, int* edge_id) {

	int tx, ty, tz, i;
	float val[6];

	switch (j) {
		case 0 : 
			tx = x;		ty = y;		tz = z;
			break;
		case 1 : 
			tx = x+1;	ty = y;		tz = z;
			break;
		case 2 : 
			tx = x+1;	ty = y;		tz = z+1;
			break;
		case 3 : 
			tx = x;		ty = y;		tz = z+1;
			break;
		case 4 : 
			tx = x;		ty = y+1;	tz = z;
			break;
		case 5 : 
			tx = x+1;	ty = y+1;	tz = z;
			break;
		case 6 : 
			tx = x+1;	ty = y+1;	tz = z+1;
			break;
		case 7 : 
			tx = x;		ty = y+1;	tz = z+1;
			break;
	}

	val[0] = getValue(tx*cell_size, ty*cell_size, (tz-1)*cell_size);	// up
	val[1] = getValue(tx*cell_size, ty*cell_size, (tz+1)*cell_size);	// down
	val[2] = getValue((tx-1)*cell_size, ty*cell_size, tz*cell_size);	// left
	val[3] = getValue((tx+1)*cell_size, ty*cell_size, tz*cell_size);	// right
	val[4] = getValue(tx*cell_size, (ty-1)*cell_size, tz*cell_size);	// front
	val[5] = getValue(tx*cell_size, (ty+1)*cell_size, tz*cell_size);	// back

	for(i = 0; i < 6; i++) {
		if(val[i] > iso_val) edge_id[i] = 1;
	}
}

void Octree::assign_refine_sign_hexa(geoframe& geofrm, float err_tol) {

	int x, y, z, valid_leaf, cell_size, level, my_bool_2, my_bool_4, num_id_face;
	int oc_id[8], edge_id[6], i, j, k, m, num_id, refine_sign, vv[4];
	unsigned int vtx[8];
	float val[8];

	for(i = 0; i < octcell_num; i++) vtx_idx_arr_refine[i] = -1;

	for(i = 0; i < leaf_num; i++) {

		valid_leaf = cut_array[i] ;
		level = get_level(valid_leaf) ;
		cell_size = (dim[0]-1)/(1<<level);
		octcell2xyz(valid_leaf, x, y, z, level);
		getCellValues(valid_leaf, level, val);

		//interior vertex at a uniform starting level
		for (j = 0; j < 8; j++ ) {
			if (is_vflag_on(x, y, z, level, j)) continue;
			if( minmax[valid_leaf].min <= iso_val && (val[j] < iso_val) ) {
				if(is_min_vertex(valid_leaf, j, vtx, geofrm)) {
					vflag_on(x, y, z, level, j);
					find_oc_id_hexa(x, y, z, level, j, oc_id);
					for(k = 0; k < 6; k++) edge_id[k] = 0;
					find_edge_id_hexa(x, y, z, cell_size, j, edge_id);

					num_id = 0;
					for(k = 0; k < 8; k++) {
						if(get_err_grad(oc_id[k]) > err_tol) {
							num_id++;
							//vtx_idx_arr_refine[oc_id[k]] = 1;
						}
					}

					if(num_id > 0) {// > 3
						//for(k = 0; k < 8; k++) vtx_idx_arr_refine[oc_id[k]] = 1;
						for(k = 0; k < 8; k++) {
							if(get_err_grad(oc_id[k]) > err_tol) vtx_idx_arr_refine[oc_id[k]] = 1;
						}
					}
				}
			}
		}	// end j

	}		// end i

	vflag_clear();

	refine_sign = 1;

	while(refine_sign == 1) {

		refine_sign = 0;

		for (i = 0; i < leaf_num; i++ ) {
			valid_leaf = cut_array[i] ;
			level = get_level(valid_leaf) ;
			cell_size = (dim[0]-1)/(1<<level);
			octcell2xyz(valid_leaf, x, y, z, level);
			getCellValues(valid_leaf, level, val);

			//interior vertex at a uniform starting level
			for (j = 0; j < 8; j++ ) {
				if (is_vflag_on(x, y, z, level, j)) continue;
				if( minmax[valid_leaf].min <= iso_val && (val[j] < iso_val) ) {
					if(is_min_vertex(valid_leaf, j, vtx, geofrm)) {
						vflag_on(x, y, z, level, j);
						find_oc_id_hexa(x, y, z, level, j, oc_id);
						for(k = 0; k < 6; k++) edge_id[k] = 0;
						find_edge_id_hexa(x, y, z, cell_size, j, edge_id);

						num_id = 0;
						for(k = 0; k < 8; k++) {if(vtx_idx_arr_refine[oc_id[k]] == 1) num_id++;}

						my_bool_2 = (vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[1]] == 1) ||
									(vtx_idx_arr_refine[oc_id[1]] == 1 && vtx_idx_arr_refine[oc_id[2]] == 1) ||
									(vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1) ||
									(vtx_idx_arr_refine[oc_id[3]] == 1 && vtx_idx_arr_refine[oc_id[0]] == 1) ||
									(vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1) ||
									(vtx_idx_arr_refine[oc_id[5]] == 1 && vtx_idx_arr_refine[oc_id[6]] == 1) ||
									(vtx_idx_arr_refine[oc_id[6]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1) ||
									(vtx_idx_arr_refine[oc_id[7]] == 1 && vtx_idx_arr_refine[oc_id[4]] == 1) ||
									(vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[4]] == 1) ||
									(vtx_idx_arr_refine[oc_id[1]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1) ||
									(vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[6]] == 1) ||
									(vtx_idx_arr_refine[oc_id[3]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1);

						my_bool_4 = (vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[1]] == 1 &&
									 vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1) ||
									(vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1 &&
									 vtx_idx_arr_refine[oc_id[6]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1) ||
									(vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1 &&
									 vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1) ||
									(vtx_idx_arr_refine[oc_id[1]] == 1 && vtx_idx_arr_refine[oc_id[2]] == 1 &&
									 vtx_idx_arr_refine[oc_id[5]] == 1 && vtx_idx_arr_refine[oc_id[6]] == 1) ||
									(vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[1]] == 1 &&
									 vtx_idx_arr_refine[oc_id[4]] == 1 && vtx_idx_arr_refine[oc_id[5]] == 1) ||
									(vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1 &&
									 vtx_idx_arr_refine[oc_id[6]] == 1 && vtx_idx_arr_refine[oc_id[7]] == 1);

						if(num_id<2 || (num_id==2 && my_bool_2) || (num_id==4 && my_bool_4) || num_id==8) continue;
						else if(num_id == 2) {
							for(k = 0; k < 6; k++) {
								if(k == 0) {vv[0] = 0;	vv[1] = 1;	vv[2] = 2;	vv[3] = 3;}	// up
								if(k == 1) {vv[0] = 4;	vv[1] = 7;	vv[2] = 6;	vv[3] = 5;}	// down
								if(k == 2) {vv[0] = 0;	vv[1] = 3;	vv[2] = 7;	vv[3] = 4;}	// left
								if(k == 3) {vv[0] = 1;	vv[1] = 5;	vv[2] = 6;	vv[3] = 2;}	// right
								if(k == 4) {vv[0] = 1;	vv[1] = 0;	vv[2] = 4;	vv[3] = 5;}	// front
								if(k == 5) {vv[0] = 3;	vv[1] = 2;	vv[2] = 6;	vv[3] = 7;}	// back
								num_id_face = 0;
								for(m = 0; m < 4; m++) {
									if(vtx_idx_arr_refine[oc_id[vv[m]]] == 1) num_id_face++;
								}
								if(num_id_face == 2) break;
							}
							if(num_id_face == 2) {
								for(m = 0; m < 4; m++) {
									if(vtx_idx_arr_refine[oc_id[vv[m]]] != 1) {
										vtx_idx_arr_refine[oc_id[vv[m]]] = 1; refine_sign = 1;}
								}
							}
							else {
								for(k = 0; k < 8; k++) {
									if(vtx_idx_arr_refine[oc_id[k]] != 1) {
										vtx_idx_arr_refine[oc_id[k]] = 1; refine_sign = 1;}
								}
							}
						}
						else if(num_id == 3) {
							for(k = 0; k < 6; k++) {
								if(k == 0) {vv[0] = 0;	vv[1] = 1;	vv[2] = 2;	vv[3] = 3;}	// up
								if(k == 1) {vv[0] = 4;	vv[1] = 7;	vv[2] = 6;	vv[3] = 5;}	// down
								if(k == 2) {vv[0] = 0;	vv[1] = 3;	vv[2] = 7;	vv[3] = 4;}	// left
								if(k == 3) {vv[0] = 1;	vv[1] = 5;	vv[2] = 6;	vv[3] = 2;}	// right
								if(k == 4) {vv[0] = 1;	vv[1] = 0;	vv[2] = 4;	vv[3] = 5;}	// front
								if(k == 5) {vv[0] = 3;	vv[1] = 2;	vv[2] = 6;	vv[3] = 7;}	// back
								num_id_face = 0;
								for(m = 0; m < 4; m++) {
									if(vtx_idx_arr_refine[oc_id[vv[m]]] == 1) num_id_face++;
								}
								if(num_id_face == 3) break;
							}
							if(num_id_face == 3) {
								for(m = 0; m < 4; m++) {
									if(vtx_idx_arr_refine[oc_id[vv[m]]] != 1) {
										vtx_idx_arr_refine[oc_id[vv[m]]] = 1; refine_sign = 1;}
								}
							}
							else {
								for(k = 0; k < 8; k++) {
									if(vtx_idx_arr_refine[oc_id[k]] != 1) {
										vtx_idx_arr_refine[oc_id[k]] = 1; refine_sign = 1;}
								}
							}
						}
						else { 
							for(k = 0; k < 8; k++) {
								if(vtx_idx_arr_refine[oc_id[k]] != 1) {
									vtx_idx_arr_refine[oc_id[k]] = 1; refine_sign = 1;}
							}

						}	// end if(num_id)
					}
				}

			}	// end j

		}		// end i

		vflag_clear();

	}			// end while

}

void Octree::find_oc_id(int x, int y, int z, int level, int j, int intersect_id, int* oc_id) {

	int i, temp_id[4];

	oc_id[0] = xyz2octcell(x, y, z, level);

	switch (j) {
		case 0 : 
			oc_id[1] = xyz2octcell(x,	y,		z-1,	level);
			oc_id[2] = xyz2octcell(x,	y-1,	z-1,	level);
			oc_id[3] = xyz2octcell(x,	y-1,	z,		level);
			break;

		case 1 : 
			oc_id[1] = xyz2octcell(x+1,	y,		z,		level);
			oc_id[2] = xyz2octcell(x+1,	y-1,	z,		level);
			oc_id[3] = xyz2octcell(x,	y-1,	z,		level);
			break;

		case 2 : 
			oc_id[1] = xyz2octcell(x,	y,		z+1,	level);
			oc_id[2] = xyz2octcell(x,	y-1,	z+1,	level);
			oc_id[3] = xyz2octcell(x,	y-1,	z,		level);
			break;

		case 3 : 
			oc_id[1] = xyz2octcell(x,	y-1,	z,		level);
			oc_id[2] = xyz2octcell(x-1,	y-1,	z,		level);
			oc_id[3] = xyz2octcell(x-1,	y,		z,		level);
			break;

		case 4 : 
			oc_id[1] = xyz2octcell(x,	y+1,	z,		level);
			oc_id[2] = xyz2octcell(x,	y+1,	z-1,	level);
			oc_id[3] = xyz2octcell(x,	y,		z-1,	level);
			break;

		case 5 :
			oc_id[1] = xyz2octcell(x,	y+1,	z,		level);
			oc_id[2] = xyz2octcell(x+1,	y+1,	z,		level);
			oc_id[3] = xyz2octcell(x+1,	y,		z,		level);
			break;

		case 6 :
			oc_id[1] = xyz2octcell(x,	y+1,	z,		level);
			oc_id[2] = xyz2octcell(x,	y+1,	z+1,	level);
			oc_id[3] = xyz2octcell(x,	y,		z+1,	level);
			break;

		case 7 :
			oc_id[1] = xyz2octcell(x-1,	y,		z,		level);
			oc_id[2] = xyz2octcell(x-1,	y+1,	z,		level);
			oc_id[3] = xyz2octcell(x,	y+1,	z,		level);
			break;

		case 8 :
			oc_id[1] = xyz2octcell(x-1,	y,		z,		level);
			oc_id[2] = xyz2octcell(x-1,	y,		z-1,	level);
			oc_id[3] = xyz2octcell(x,	y,		z-1,	level);
			break;

		case 9 :
			oc_id[1] = xyz2octcell(x,	y,		z-1,	level);
			oc_id[2] = xyz2octcell(x+1,	y,		z-1,	level);
			oc_id[3] = xyz2octcell(x+1,	y,		z,		level);
			break;

		case 10 :
			oc_id[1] = xyz2octcell(x,	y,		z+1,	level);
			oc_id[2] = xyz2octcell(x-1,	y,		z+1,	level);
			oc_id[3] = xyz2octcell(x-1,	y,		z,		level);
			break;

		case 11 :
			oc_id[1] = xyz2octcell(x+1,	y,		z,		level);
			oc_id[2] = xyz2octcell(x+1,	y,		z+1,	level);
			oc_id[3] = xyz2octcell(x,	y,		z+1,	level);
			break;
		}
					
		for(i = 0; i < 4; i++)	temp_id[i] = oc_id[i];
		if (intersect_id == -1) 
			for (i = 0; i < 4; i++) oc_id[i] = temp_id[3-i];

}

void Octree::find_vtx_new(geoframe& geofrm, int x, int y, int z, int level, int j, int intersect_id, unsigned int* vtx_new) {

	int i, oc_id_new[4], cell_size, vert, tx, ty, tz;
	float vtx[3], norm[3], val[8];

	switch (j) {
		case 0 : 
			oc_id_new[0] = xyz2octcell(2*x, 2*y, 2*z, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x, 2*y, 2*z, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x+1, 2*y, 2*z, level+1, j, intersect_id, oc_id_new);
			break;

		case 1 : 
			oc_id_new[0] = xyz2octcell(2*x+1, 2*y, 2*z, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x+1, 2*y, 2*z, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x+1, 2*y, 2*z+1, level+1, j, intersect_id, oc_id_new);
			break;

		case 2 : 
			oc_id_new[0] = xyz2octcell(2*x, 2*y, 2*z+1, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x, 2*y, 2*z+1, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x+1, 2*y, 2*z+1, level+1, j, intersect_id, oc_id_new);
			break;

		case 3 : 
			oc_id_new[0] = xyz2octcell(2*x, 2*y, 2*z, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x, 2*y, 2*z, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x, 2*y, 2*z+1, level+1, j, intersect_id, oc_id_new);
			break;

		case 4 : 
			oc_id_new[0] = xyz2octcell(2*x, 2*y+1, 2*z, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x, 2*y+1, 2*z, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x+1, 2*y+1, 2*z, level+1, j, intersect_id, oc_id_new);
			break;

		case 5 :
			oc_id_new[0] = xyz2octcell(2*x+1, 2*y+1, 2*z, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x+1, 2*y+1, 2*z, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x+1, 2*y+1, 2*z+1, level+1, j, intersect_id, oc_id_new);
			break;

		case 6 :
			oc_id_new[0] = xyz2octcell(2*x, 2*y+1, 2*z+1, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x, 2*y+1, 2*z+1, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x+1, 2*y+1, 2*z+1, level+1, j, intersect_id, oc_id_new);
			break;

		case 7 :
			oc_id_new[0] = xyz2octcell(2*x, 2*y+1, 2*z, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x, 2*y+1, 2*z, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x, 2*y+1, 2*z+1, level+1, j, intersect_id, oc_id_new);
			break;

		case 8 :
			oc_id_new[0] = xyz2octcell(2*x, 2*y, 2*z, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x, 2*y, 2*z, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x, 2*y+1, 2*z, level+1, j, intersect_id, oc_id_new);
			break;

		case 9 :
			oc_id_new[0] = xyz2octcell(2*x+1, 2*y, 2*z, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x+1, 2*y, 2*z, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x+1, 2*y+1, 2*z, level+1, j, intersect_id, oc_id_new);
			break;

		case 10 :
			oc_id_new[0] = xyz2octcell(2*x, 2*y, 2*z+1, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x, 2*y, 2*z+1, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x, 2*y+1, 2*z+1, level+1, j, intersect_id, oc_id_new);
			break;

		case 11 :
			oc_id_new[0] = xyz2octcell(2*x+1, 2*y, 2*z+1, level+1);
			getCellValues(oc_id_new[0], level+1, val);
			if(is_intersect(val, j) == 1 || is_intersect(val, j) == -1)
				find_oc_id(2*x+1, 2*y, 2*z+1, level+1, j, intersect_id, oc_id_new);
			else
				find_oc_id(2*x+1, 2*y+1, 2*z+1, level+1, j, intersect_id, oc_id_new);
			break;
		}
					
		for(i = 0; i < 4; i++) {
			octcell2xyz(oc_id_new[i], tx, ty, tz, level+1);

			cell_size = (dim[0]-1)/(1<<(level+1));
			//get_vtx(tx, ty, tz, level+1, vtx);
			get_solution(oc_id_new[i], vtx);
			getVertGrad(tx*cell_size, ty*cell_size, tz*cell_size, norm);
			//get_VtxNorm(vtx, norm);
			if(in_out == 0) {
				if ((vert = vtx_idx_arr[xyz2octcell(tx, ty, tz, level+1)]) == -1) {
					vert = geofrm.AddVert(vtx, norm);
					geofrm.AddBound(vert, 1);
					vtx_idx_arr[xyz2octcell(tx, ty, tz, level+1)] = vert;
				}
			}
			else {
				if ((vert = vtx_idx_arr_in[xyz2octcell(tx, ty, tz, level+1)]) == -1) {
					vert = geofrm.AddVert(vtx, norm);
					geofrm.AddBound(vert, -1);
					vtx_idx_arr_in[xyz2octcell(tx, ty, tz, level+1)] = vert;
				}
			}
			vtx_new[i] = (unsigned int) vert;

		}

}

void Octree::get_vtx_new(geoframe& geofrm, int oc_id, unsigned int vtx)
{
/*
	// average all the minimizer points inside this cell
	int xx, yy, zz, tx, ty, tz, level, i, j, my_id, num;
	float pos[3], pos_average[3];

	level = get_level(oc_id) ;
	octcell2xyz(oc_id, xx, yy, zz, level);

	if(i == 0 || i == 3 || i == 4 || i == 7) tx = 2*xx;
	else tx = 2*xx + 1;
	if(i < 4) ty = 2*yy;
	else ty = 2*yy + 1;
	if(i == 0 || i == 1 || i == 4 || i == 5) tz = 2*zz;
	else tz = 2*zz + 1;

	num = 0;
	for(i = 0; i < 3; i++) pos_average[i] = 0.0f;

	for(i = 0; i < 8; i++) {
		my_id = xyz2octcell(tx, ty, tz, level+1);
		if(is_skipcell(my_id) == 0) {
			num++;
			get_solution(my_id, pos);
			for(j = 0; j < 3; j++) pos_average[j] += pos[j];
		}
	}
	if(num > 1) {
		for(j = 0; j < 3; j++) {
			geofrm.verts[vtx][j] = (geofrm.verts[vtx][j] + pos_average[j]/num)/2.0f;
		}
	}

	// find the intersection point with the trilinear function (the current resolution)
	int level, cell_size, xx, yy, zz, i, my_bool;
	float val[8], x, y, z, step, func, func0;

	level = get_level(oc_id) ;
	cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_id, xx, yy, zz, level);
	getCellValues(oc_id, level, val);

	x = geofrm.verts[vtx][0]/cell_size - (float)xx;
	y = geofrm.verts[vtx][1]/cell_size - (float)yy;
	z = geofrm.verts[vtx][2]/cell_size - (float)zz;

	func0 = (1-x)*(1-y)*(1-z)*val[0] + (1-x)*(1-y)*z*val[3] + (1-x)*y*(1-z)*val[4]
				+ x*(1-y)*(1-z)*val[1] + (1-x)*y*z*val[7] + x*(1-y)*z*val[2]
				+ x*y*(1-z)*val[5] + x*y*z*val[6];

	step = 0.001f;
	if(func0 < 0.0f) step = -0.001f;
	for(i = 1; i < 16000; i++) {
		x = x - geofrm.normals[vtx][0]*step;
		y = y - geofrm.normals[vtx][1]*step;
		z = z - geofrm.normals[vtx][2]*step;
		func = (1-x)*(1-y)*(1-z)*val[0] + (1-x)*(1-y)*z*val[3] + (1-x)*y*(1-z)*val[4]
					+ x*(1-y)*(1-z)*val[1] + (1-x)*y*z*val[7] + x*(1-y)*z*val[2]
					+ x*y*(1-z)*val[5] + x*y*z*val[6];
		my_bool = (x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f);
		if(func*func0 <= 0.0f || my_bool == 0) break;
		func0 = func;
	}

	geofrm.verts[vtx][0] = ((float)xx + x)*(float)cell_size;
	geofrm.verts[vtx][1] = ((float)yy + y)*(float)cell_size;
	geofrm.verts[vtx][2] = ((float)zz + z)*(float)cell_size;

	geofrm.normals[vtx][0] = (1-y)*(1-z)*(val[1] - val[0]) + (1-y)*z*(val[2] - val[3])
								+ y*(1-z)*(val[5] - val[4]) + y*z*(val[6] - val[7]);
	geofrm.normals[vtx][1] = (1-x)*(1-z)*(val[4] - val[0]) + (1-x)*z*(val[7] - val[3])
								+ x*(1-z)*(val[5] - val[1]) + x*z*(val[6] - val[2]);
	geofrm.normals[vtx][2] = (1-x)*(1-z)*(val[3] - val[0]) + (1-x)*y*(val[7] - val[4])
								+ x*(1-y)*(val[2] - val[1]) + x*y*(val[6] - val[5]);
*/
	// find the intersection point with the trilinear function (the finest resolution)
	int level, cell_size, xx, yy, zz, tx, ty, tz, i, my_bool, my_id;
	float val[8], x, y, z, mx, my, mz, step, func, func0, n[3], temp;

	level = get_level(oc_id) ;
	cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_id, xx, yy, zz, level);

	mx = geofrm.verts[vtx][0]/cell_size - (float)xx;
	my = geofrm.verts[vtx][1]/cell_size - (float)yy;
	mz = geofrm.verts[vtx][2]/cell_size - (float)zz;

	tx = cell_size*xx + (int)(mx*cell_size);
	ty = cell_size*yy + (int)(my*cell_size);
	tz = cell_size*zz + (int)(mz*cell_size);

	my_id = xyz2octcell(tx, ty, tz, oct_depth);
	getCellValues(my_id, oct_depth, val);

	x = geofrm.verts[vtx][0] - (float)tx;
	y = geofrm.verts[vtx][1] - (float)ty;
	z = geofrm.verts[vtx][2] - (float)tz;

	func0 = (1-x)*(1-y)*(1-z)*val[0] + (1-x)*(1-y)*z*val[3] + (1-x)*y*(1-z)*val[4]
				+ x*(1-y)*(1-z)*val[1] + (1-x)*y*z*val[7] + x*(1-y)*z*val[2]
				+ x*y*(1-z)*val[5] + x*y*z*val[6] - iso_val;

	for(i = 0; i < 3; i++) n[i] = geofrm.normals[vtx][i];
	temp = (float) sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
	if(temp > 0.001f) for(i = 0; i < 3; i++) n[i] /= temp;

	if(fabs(func0) >= 0.001733f) { //0.2
		for(i = 1; i < 1000*cell_size; i++) {

			if(fabs(func0) < 0.001733f) break;

			step = -0.001f;
			if(func0 < 0.0f) step = 0.001f;
			/*		
			n[0] = (1-y)*(1-z)*(val[1] - val[0]) + (1-y)*z*(val[2] - val[3])
										+ y*(1-z)*(val[5] - val[4]) + y*z*(val[6] - val[7]);
			n[1] = (1-x)*(1-z)*(val[4] - val[0]) + (1-x)*z*(val[7] - val[3])
										+ x*(1-z)*(val[5] - val[1]) + x*z*(val[6] - val[2]);
			n[2] = (1-x)*(1-z)*(val[3] - val[0]) + (1-x)*y*(val[7] - val[4])
										+ x*(1-y)*(val[2] - val[1]) + x*y*(val[6] - val[5]);
			temp = (float) sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
			if(temp > 0.001f) for(int j = 0; j < 3; j++) n[j] /= temp;
			*/
			x = x + n[0]*step;
			y = y + n[1]*step;
			z = z + n[2]*step;
			func = (1-x)*(1-y)*(1-z)*val[0] + (1-x)*(1-y)*z*val[3] + (1-x)*y*(1-z)*val[4]
						+ x*(1-y)*(1-z)*val[1] + (1-x)*y*z*val[7] + x*(1-y)*z*val[2]
						+ x*y*(1-z)*val[5] + x*y*z*val[6] - iso_val;
			my_bool = (x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f);
			if(func*func0 <= 0.0f && my_bool) break;
			func0 = func;
			if(my_bool == 0) {
				if(x < 0.0f) {tx--;	x = x + 1.0f;}	else {tx++;	x = x - 1.0f;}
				if(y < 0.0f) {ty--;	y = y + 1.0f;}	else {ty++; y = y - 1.0f;}
				if(z < 0.0f) {tz--;	z = z + 1.0f;}	else {tz++; z = z - 1.0f;}
				my_id = xyz2octcell(tx, ty, tz, oct_depth);
				getCellValues(my_id, oct_depth, val);
			}
		}

		//if(my_bool == 0 && i < 1000*cell_size) {
			geofrm.verts[vtx][0] = (float)tx + x;
			geofrm.verts[vtx][1] = (float)ty + y;
			geofrm.verts[vtx][2] = (float)tz + z;
		//}
	
		/*
		// calculate normal using the trilinear function
		geofrm.normals[vtx][0] = (1-y)*(1-z)*(val[1] - val[0]) + (1-y)*z*(val[2] - val[3])
									+ y*(1-z)*(val[5] - val[4]) + y*z*(val[6] - val[7]);
		geofrm.normals[vtx][1] = (1-x)*(1-z)*(val[4] - val[0]) + (1-x)*z*(val[7] - val[3])
									+ x*(1-z)*(val[5] - val[1]) + x*z*(val[6] - val[2]);
		geofrm.normals[vtx][2] = (1-x)*(1-z)*(val[3] - val[0]) + (1-x)*y*(val[7] - val[4])
									+ x*(1-y)*(val[2] - val[1]) + x*y*(val[6] - val[5]);
		
		temp = sqrt(geofrm.normals[vtx][0]*geofrm.normals[vtx][0] +
				geofrm.normals[vtx][1]*geofrm.normals[vtx][1] + geofrm.normals[vtx][2]*geofrm.normals[vtx][2]);
		for(i = 0; i < 3; i++)	geofrm.normals[vtx][i] /= temp;
		
		// calculate curvature using the trilinear function
		float df2_dxdy = -(1-z)*(val[1] - val[0]) + -z*(val[2] - val[3]) + (1-z)*(val[5] - val[4]) + z*(val[6] - val[7]);
		float df2_dxdz = -(1-y)*(val[1] - val[0]) + (1-y)*(val[2] - val[3]) - y*(val[5] - val[4]) + y*(val[6] - val[7]);
		float df2_dydz = -(1-x)*(val[4] - val[0]) + (1-x)*(val[7] - val[3]) - x*(val[5] - val[1]) + x*(val[6] - val[2]);

		float coef = sqrt((df2_dxdy*df2_dxdy + df2_dxdz*df2_dxdz + df2_dydz*df2_dydz)/3.0f);
		float theta = acos(df2_dxdy*df2_dxdz*df2_dydz/(coef*coef*coef));
		if(fabs(df2_dxdy*df2_dxdz*df2_dydz/(coef*coef*coef)) > 1.0) theta = 0.0f;

		float eigenvalue_0 = coef*2.0f*cos(theta/3.0f);
		float eigenvalue_1 = coef*2.0f*cos((theta + 2.0f*3.1415926f)/3.0f);
		float eigenvalue_2 = coef*2.0f*cos((theta - 2.0f*3.1415926f)/3.0f);

		float kapa_1 = eigenvalue_0;
		float kapa_2 = eigenvalue_0;
		if(kapa_1 < eigenvalue_1) kapa_1 = eigenvalue_1;
		if(kapa_1 < eigenvalue_2) kapa_1 = eigenvalue_2;
		if(kapa_2 > eigenvalue_1) kapa_2 = eigenvalue_1;
		if(kapa_2 > eigenvalue_2) kapa_2 = eigenvalue_2;

		geofrm.curvatures[vtx][0] = kapa_1;
		geofrm.curvatures[vtx][1] = kapa_2;
		*/
	}

}


}
