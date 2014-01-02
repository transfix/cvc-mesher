#include<cellqueue.h>
#include"octree.h"
#include"e_face.h"
#include"cubes.h"
#include<stdio.h>
#include"vtkMarchingCubesCases.h"
#include"pcio.h"
#include"LBIE_geoframe.h"
#include<assert.h>
#include<time.h>
#include<algorithm>

#include"normalspline.h"

namespace LBIE
{

  //int* interpol_bit;

  Octree::Octree()
    : iso_val(0.0), iso_val_in(0.0), leaf_num(0),
      octcell_num(0), cell_num(0), oct_depth(0),
      flag_type(0), flag_normal(0), in_out(0),
      flag_extend(0), prop_flag(0), interior_flag(1),
      nverts(0), ncells(0), vol_min(0.0), vol_max(0.0)
  {
    std::fill(level_res.begin(),level_res.end(),0);
    e_face_initialization();
    std::fill(dim.begin(),dim.end(),0);
    std::fill(olddim.begin(),olddim.end(),0);
    std::fill(orig.begin(),orig.end(),0.0);
    std::fill(span.begin(),span.end(),0.0);
    std::fill(spans.begin(),spans.end(),0.0);
    std::fill(minext.begin(),minext.end(),0.0);
    std::fill(maxext.begin(),maxext.end(),0.0);
  }

  Octree::Octree(const Octree& oc)
  {
    iso_val = oc.iso_val;
    iso_val_in = oc.iso_val_in;
    leaf_num = oc.leaf_num;
    oct_array = oc.oct_array;
    octcell_num = oc.octcell_num;
    cell_num = oc.cell_num;
    oct_depth = oc.oct_depth;
    level_res = oc.level_res;
    cut_array = oc.cut_array;
    flag_type = oc.flag_type;
    flag_normal = oc.flag_normal;
    in_out = oc.in_out;
    flag_extend = oc.flag_extend;
    e_face = oc.e_face;
    qef_array = oc.qef_array;
    qef_array_in = oc.qef_array_in;
    vtx_idx_arr = oc.vtx_idx_arr;
    vtx_idx_arr_in = oc.vtx_idx_arr_in;
    grid_idx_arr = oc.grid_idx_arr;
    vtx_idx_arr_refine = oc.vtx_idx_arr_refine;
    orig_vol = oc.orig_vol;
    ebit = oc.ebit;
    vbit = oc.vbit;
    BSplineCoeff = oc.BSplineCoeff;
    minmax = oc.minmax;
    prop_flag = oc.prop_flag;
    interior_flag = oc.interior_flag;
    prop_fname = oc.prop_fname;
    nverts = oc.nverts;
    ncells = oc.ncells;
    dim = oc.dim;
    olddim = oc.olddim;
    orig = oc.orig;
    span = oc.span;
    vol_min = oc.vol_min;
    vol_max = oc.vol_max;
    spans = oc.spans;
    minext = oc.minext;
    maxext = oc.maxext;
  }

  Octree& Octree::operator=(const Octree& oc)
  {
    iso_val = oc.iso_val;
    iso_val_in = oc.iso_val_in;
    leaf_num = oc.leaf_num;
    oct_array = oc.oct_array;
    octcell_num = oc.octcell_num;
    cell_num = oc.cell_num;
    oct_depth = oc.oct_depth;
    level_res = oc.level_res;
    cut_array = oc.cut_array;
    flag_type = oc.flag_type;
    flag_normal = oc.flag_normal;
    in_out = oc.in_out;
    flag_extend = oc.flag_extend;
    e_face = oc.e_face;
    qef_array = oc.qef_array;
    qef_array_in = oc.qef_array_in;
    vtx_idx_arr = oc.vtx_idx_arr;
    vtx_idx_arr_in = oc.vtx_idx_arr_in;
    grid_idx_arr = oc.grid_idx_arr;
    vtx_idx_arr_refine = oc.vtx_idx_arr_refine;
    orig_vol = oc.orig_vol;
    ebit = oc.ebit;
    vbit = oc.vbit;
    BSplineCoeff = oc.BSplineCoeff;
    minmax = oc.minmax;
    prop_flag = oc.prop_flag;
    interior_flag = oc.interior_flag;
    prop_fname = oc.prop_fname;
    nverts = oc.nverts;
    ncells = oc.ncells;
    dim = oc.dim;
    olddim = oc.olddim;
    orig = oc.orig;
    span = oc.span;
    vol_min = oc.vol_min;
    vol_max = oc.vol_max;
    spans = oc.spans;
    minext = oc.minext;
    maxext = oc.maxext;
    return *this;
  }

Octree::~Octree()
{
  //if(oct_array)
  //free(oct_array);
  //	if(cut_array)
  //	free(cut_array);
  //if(orig_vol)
  //free(orig_vol);
	//if(vtx_idx_arr)
	//free(vtx_idx_arr);
	//if(vtx_idx_arr_in)
	//free(vtx_idx_arr_in);
	//if(grid_idx_arr)
	//free(grid_idx_arr);
  //if(minmax)
  //free(minmax);
  //if(BSplineCoeff)
  //free(BSplineCoeff);
	//if(vtx_idx_arr_refine)
	//free(vtx_idx_arr_refine);
}

//void Octree::setVolume(SimpleVolumeData * volData)
void Octree::setVolume(const VolMagick::Volume& volData)
{
  VolMagick::Volume volumeData;
  volumeData = volData;
  volumeData.voxelType(VolMagick::Float); //I think this library wants a float*
  loadData(volumeData);
  
  int i;
  for (i=0;i<=oct_depth;i++) {
    level_res[i]=(1<<i);
  }
  construct_octree("errfile.rawiv");
  vol_min=minmax[0].min;
  vol_max=minmax[0].max;
}

void Octree::prop_init(const char* rawiv_fname)
{
  //	strcpy(prop_fname,rawiv_fname);
  prop_fname = rawiv_fname;
	
	/*prop_fp = fopen(rawiv_fname,"rb");

	if (prop_fp==NULL) {
		printf("wrong name : %s\n",rawiv_fname);
		return;
	}*/
  prop_flag = 1;
}

void Octree::Octree_loadVolume(const char* dx_fname)
{
  VolMagick::Volume volumeData;
  VolMagick::readVolumeFile(volumeData,dx_fname);
  volumeData.voxelType(VolMagick::Float); //I think this library wants a float*
  loadData(volumeData);
  
  int i;
  for (i=0;i<=oct_depth;i++) {
    level_res[i]=(1<<i);
  }
  construct_octree((char*)dx_fname);
  vol_min=minmax[0].min;
  vol_max=minmax[0].max;
}

//void Octree::loadData(SimpleVolumeData *volumeData, float* vol)
void Octree::loadData(const VolMagick::Volume& volumeData/*, float* vol*/)
{
  leaf_num=0;
  prop_flag = 0;

  minext[0] = volumeData.XMin();
  minext[1] = volumeData.YMin();
  minext[2] = volumeData.ZMin();
  maxext[0] = volumeData.XMax();
  maxext[1] = volumeData.YMax();
  maxext[2] = volumeData.ZMax();

  spans[0] = volumeData.XSpan();
  spans[1] = volumeData.YSpan();
  spans[2] = volumeData.ZSpan();
	
  dim[0] = volumeData.XDim();
  dim[1] = volumeData.YDim();
  dim[2] = volumeData.ZDim();

  //printf("dim:%d %d %d\n",dim[0],dim[1],dim[2]);	

  orig[0] = 0.0f;	orig[1] = 0.0f;	orig[2] = 0.0f;
  span[0] = spans[0]; span[1] = spans[1]; span[2] = spans[2];
  //span[0] = 1.0f;	span[1] = 1.0f;	span[2] = 1.0f;
  //spans[0] = 1.0f;	spans[1] = 1.0f;	spans[2] = 1.0f;


  if(interior_flag)
    {
      nverts = dim[0] * dim[1] * dim[2];
      ncells = (dim[0] -1)*(dim[1] -1)*(dim[2] -1);

      initMem();
    }
  else
    {
      olddim[0] = dim[0];
      olddim[1] = dim[1];
      olddim[2] = dim[2];
      dim[0] = 2*dim[0]-1;
      dim[1] = 2*dim[1]-1;
      dim[2] = 2*dim[2]-1;

      nverts = dim[0] * dim[1] * dim[2];
      ncells = (dim[0] -1)*(dim[1] -1)*(dim[2] -1);

      initMem();
    }

  int i, j, k;
  //interior_flag =0;
  if(interior_flag)
    {	
      for(i = 0; i < dim[0]*dim[1]*dim[2]; i++) {
	//vol[i] = -tmp[i];//volumeData->getValueAt(0,i); //tmp[i];
	orig_vol[i] = -volumeData(i);
      }
    }
  else
    {

      int incr = (dim[0]-olddim[0])/2;
      float center = (dim[1] - 1.0f)/2.0f;

      float r0, r;
      int index,index0;
      r0 = center - 0.001f;
      for(k = 0; k < dim[2]; k++)
	for(j = 0; j < dim[1]; j++)
	  for(i = 0; i < dim[0]; i++) {
	    index0 = (k-incr)*olddim[0]*olddim[1]+(j-incr)*olddim[0]+(i-incr);
	    index = k*dim[0]*dim[1] + j*dim[0] + i;
				

	    if(i > incr && i < olddim[0]+incr &&
	       j > incr && j < olddim[1]+incr &&
	       k > incr && k < olddim[2]+incr)
	      {
		//vol[index] = tmp[index0] - 1.0f;
		orig_vol[index] = volumeData(index0) - 1.0f;
	      }
	    else
	      {
		orig_vol[index] = -1.0f;
	      }
	    r = (float) sqrt((i-center)*(i-center) + (j-center)*(j-center) + (k-center)*(k-center));
	    if(orig_vol[index] < 0.0) {
	      if((r>r0) || ((r0-r) < -orig_vol[index]))	orig_vol[index] = r - r0;
	    }
	  }
    }
}

/*void Octree::Octree_loadObj(const char* obj_fname)
{
	vol_fp=fopen(obj_fname,"rb");
	
	if (vol_fp==NULL) {
		printf("wrong name : %s\n",obj_fname);
		return;
	}
}*/

void Octree::initMem()
{
  /* pre-initialize */
  //oct_array = NULL;
  oct_array.clear();
  //minmax = NULL;
  minmax.clear();
  //  cut_array = NULL;
  cut_array.clear();
  //orig_vol = NULL;
  orig_vol.clear();
  //ebit = NULL;
  ebit.clear();
  //vbit = NULL;
  vbit.clear();
  //vtx_idx_arr = NULL;
  vtx_idx_arr.clear();
  //grid_idx_arr = NULL;
  grid_idx_arr.clear();
  //vtx_idx_arr_in = NULL;
  vtx_idx_arr_in.clear();
  //vtx_idx_arr_refine = NULL;
  vtx_idx_arr_refine.clear();
  //qef_array = NULL;
  //qef_array_in = NULL;
  qef_array.clear();
  qef_array_in.clear();

  /* get sizes */
  oct_depth   = get_depth(dim[0]);
  octcell_num = get_octcell_num(oct_depth);
  cell_num    = (dim[0]-1)*(dim[1]-1)*(dim[2]-1);

  /* initialize arrays */
  if(octcell_num > 0)
    {
      //oct_array   = (octcell*)malloc(sizeof(octcell)*octcell_num);
      //memset(oct_array,0,sizeof(octcell)*octcell_num);
      oct_array.resize(octcell_num);
      //minmax = (MinMax*)malloc(sizeof(MinMax)*octcell_num);
      minmax.resize(octcell_num*2); //why does construct_octree use octcell_num*2 ??
      //memset(minmax,0,sizeof(MinMax)*octcell_num);

      //ebit = (char*)malloc(sizeof(char)*octcell_num*4/8);
      ebit.resize(octcell_num*4/8);
      //vbit = (char*)malloc(sizeof(char)*octcell_num*4/8);
      vbit.resize(octcell_num*4/8);

      //vtx_idx_arr = (int*)malloc(sizeof(int)*octcell_num);
      vtx_idx_arr.resize(octcell_num);
      //vtx_idx_arr_in = (int*)malloc(sizeof(int)*octcell_num);
      vtx_idx_arr_in.resize(octcell_num);
      //vtx_idx_arr_refine = (int*)malloc(sizeof(int)*octcell_num);
      vtx_idx_arr_refine.resize(octcell_num);

      // for(int k=0;k<octcell_num;k++)
// 	{
// 	  vtx_idx_arr[k]=-1;
// 	  vtx_idx_arr_in[k]=-1;
// 	  vtx_idx_arr_refine[k]=-1;
// 	}

      std::fill(vtx_idx_arr.begin(),vtx_idx_arr.end(),-1);
      std::fill(vtx_idx_arr_in.begin(),vtx_idx_arr_in.end(),-1);
      std::fill(vtx_idx_arr_refine.begin(),vtx_idx_arr_refine.end(),-1);

      //qef_array	= (double**) malloc (sizeof(double*) * octcell_num);
      //qef_array_in = (double**) malloc (sizeof(double*) * octcell_num);
      qef_array.resize(octcell_num);
      qef_array_in.resize(octcell_num);
      //assert(qef_array!=NULL);
      //assert(qef_array_in!=NULL);
      //memset(qef_array,0,octcell_num*sizeof(double*));
      //memset(qef_array_in,0,octcell_num*sizeof(double*));

      //memset(ebit,0,octcell_num*4/8*sizeof(char));
      //memset(vbit,0,octcell_num*4/8*sizeof(char));
    }
  if(cell_num > 0)
    {
      //cut_array   = (int*)malloc(sizeof(int)*2*cell_num);
      cut_array.resize(2*cell_num);
    }
  if(dim[0]*dim[1]*dim[2] > 0)
    {
      //orig_vol    = (float*)malloc(sizeof(float)*dim[0]*dim[1]*dim[2]);
      orig_vol.resize(dim[0]*dim[1]*dim[2]);
      //grid_idx_arr = (int*)malloc(sizeof(int)*dim[0]*dim[1]*dim[2]);
      grid_idx_arr.resize(dim[0]*dim[1]*dim[2]);
      //for(int k=0;k<dim[0]*dim[1]*dim[2];k++) grid_idx_arr[k]=-1;
      std::fill(grid_idx_arr.begin(),grid_idx_arr.end(),-1);
    }
}



// given the path of a volume,
// load the volume data and construct octree from that
// allocate associate arrays
/*void Octree::Octree_init(const char* rawiv_fname)
{
	
	vol_fp=fopen(rawiv_fname,"rb");
	
	if (vol_fp==NULL) {
		printf("wrong name : %s\n",rawiv_fname);
		return;
	}

	leaf_num=0;

	read_header();   // read dimension, span, orig and so on
	
	// octree information
	oct_depth   = get_depth(dim[0]);
	octcell_num = get_octcell_num(oct_depth);
	cell_num    = (dim[0]-1)*(dim[1]-1)*(dim[2]-1);

	oct_array   = (octcell*)malloc(sizeof(octcell)*octcell_num);

	memset(oct_array,0,sizeof(octcell)*octcell_num);

	minmax		= (MinMax*)malloc(sizeof(MinMax)*octcell_num);

	memset(minmax,0,sizeof(MinMax)*octcell_num);

	cut_array   = (int*)malloc(sizeof(int)*2*cell_num);


	orig_vol    = (float*)malloc(sizeof(float)*dim[0]*dim[1]*dim[2]);

	ebit		= (char*)malloc(sizeof(char)*octcell_num*4/8);
	vbit		= (char*)malloc(sizeof(char)*octcell_num*4/8);

	vtx_idx_arr = (int*)malloc(sizeof(int)*octcell_num);
	grid_idx_arr = (int*)malloc(sizeof(int)*dim[0]*dim[1]*dim[2]);
	vtx_idx_arr_in = (int*)malloc(sizeof(int)*octcell_num);
	int k;
	for (k=0;k<octcell_num;k++) {
		vtx_idx_arr[k]=-1;	vtx_idx_arr_in[k]=-1;}
	for (k=0;k<dim[0]*dim[1]*dim[2];k++) grid_idx_arr[k]=-1;	

	qef_array	= (double**) malloc (sizeof(double*) * octcell_num);
	qef_array_in = (double**) malloc (sizeof(double*) * octcell_num);

	assert(qef_array!=NULL);
	assert(qef_array_in!=NULL);

	memset(qef_array,0,octcell_num*sizeof(double*));
	memset(qef_array_in,0,octcell_num*sizeof(double*));

	memset(ebit,0,octcell_num*4/8*sizeof(char));
	memset(vbit,0,octcell_num*4/8*sizeof(char));
	
	read_data();

	// Bspline Interpolation
	//BSplineCoeff = (float*)malloc(sizeof(float)*dim[0]*dim[1]*dim[2]);
	//TransImg2Spline(orig_vol, BSplineCoeff, dim[0], dim[1], dim[2]);

	int i;
	for (i=0;i<=oct_depth;i++) {
		level_res[i]=(1<<i);
	}
//	e_face_initialization();
	construct_octree((char*)rawiv_fname);
	vol_min=minmax[0].min;
	vol_max=minmax[0].max;
}
*/

// error computation on each cell of octree
void Octree::construct_octree(char * rawiv_fname)
{
	int oc_idx=0;
	int level;
	FILE* err_fp;
	char err_fname[256];
	float min , max;
	
	strcpy(err_fname,rawiv_fname);
	strcat(err_fname,".err");
	
	err_fp=fopen(err_fname,"rb");
	
	err_fp = 0;
	if (err_fp) {
		
		//getFloat((float*)minmax,octcell_num*2,err_fp);
	  //fread((float*)minmax,sizeof(float),octcell_num*2,err_fp);
	  fread(&(minmax[0]),sizeof(float),octcell_num*2,err_fp);
		fclose(err_fp);
		
	} else {

		while (oc_idx < octcell_num) {

			level=get_level(oc_idx);
			compute_error(oc_idx,level,min,max);
			minmax[oc_idx].min=min;
			minmax[oc_idx].max=max;
			//oct_array[oc_idx].refine_flag=0;
			oc_idx++;

		} 

		/*err_fp=fopen(err_fname, "wb");
		//putFloat((float*)minmax, octcell_num*2, err_fp);
		fwrite((float*)minmax,sizeof(float),octcell_num*2,err_fp);
		fclose(err_fp);*/
	}
	
}

void Octree::collapse()
{	
	CellQueue prev_queue, cur_queue;
	int i, oc_id, level;

	oc_id=0;

	prev_queue.Add(oc_id);
	while (!prev_queue.Empty()) {
		
		//int leaf_idx=leaf_num;
		while (prev_queue.Get(oc_id) >= 0)  {
			level=get_level(oc_id);
			
			//octcell2xyz(oc_id,x,y,z,level);
			
//			if (is_skipcell(oc_id) || level==oct_depth) {
			if (is_skipcell(oc_id) || level==oct_depth || minmax[oc_id].max < iso_val) {

				//cut_array[leaf_num++]=oc_id;
				oct_array[oc_id].refine_flag=0;

			} 
			else {
				oct_array[oc_id].refine_flag=1;
				cur_queue.Add(oc_id);
			}
		}
		
		while (cur_queue.Get(oc_id)>=0) {
			level=get_level(oc_id);
			for (i=0; i<8; i++)
				prev_queue.Add(child(oc_id,level,i));
		}
		
	}
	//return leaf_idx;
}

void Octree::collapse_interval()
{	
	CellQueue prev_queue, cur_queue;
	int i, oc_id, level;

	oc_id=0;

	prev_queue.Add(oc_id);
	while (!prev_queue.Empty()) {
		
		//int leaf_idx=leaf_num;
		while (prev_queue.Get(oc_id) >= 0)  {
			level=get_level(oc_id);
			
			//octcell2xyz(oc_id,x,y,z,level);
			
//			if (is_skipcell(oc_id) || level==oct_depth) {
			if (is_skipcell_interval(oc_id) || level==oct_depth) {

				//cut_array[leaf_num++]=oc_id;
				oct_array[oc_id].refine_flag=0;

			} 
			else {
				oct_array[oc_id].refine_flag=1;
				cur_queue.Add(oc_id);
			}
		}
		
		while (cur_queue.Get(oc_id)>=0) {
			level=get_level(oc_id);
			for (i=0; i<8; i++)
				prev_queue.Add(child(oc_id,level,i));
		}
		
	}
	//return leaf_idx;
}

int Octree::cell_comp(int oc_id, int level, float pt[12][3], float norm[12][3])
{
	float val[8], f1, f2;
	int code, e, i, j, k, edge;

	getCellValues(oc_id,level,val);

	code = 0;
	
	if (val[0] < iso_val) code |= 0x01;
	if (val[1] < iso_val) code |= 0x02;
	if (val[2] < iso_val) code |= 0x04;
	if (val[3] < iso_val) code |= 0x08;
	if (val[4] < iso_val) code |= 0x10;
	if (val[5] < iso_val) code |= 0x20;
	if (val[6] < iso_val) code |= 0x40;
	if (val[7] < iso_val) code |= 0x80;
/*
	if (val[0] > iso_val) code |= 0x01;
	if (val[1] > iso_val) code |= 0x02;
	if (val[2] > iso_val) code |= 0x04;
	if (val[3] > iso_val) code |= 0x08;
	if (val[4] > iso_val) code |= 0x10;
	if (val[5] > iso_val) code |= 0x20;
	if (val[6] > iso_val) code |= 0x40;
	if (val[7] > iso_val) code |= 0x80;
*/
	assert(code!=0 && code!=255);
		
	octcell2xyz(oc_id,i,j,k,level);
		
	for (e=0 ; e < cubeedges[code][0] ; e++) {
		edge = cubeedges[code][1+e];
		EdgeInfo *ei = &edgeinfo[edge];
		
		f1=val[ei->d1];
		f2=val[ei->d2];
	
		switch (ei->dir) {
			case 0:
				interpRect3Dpts_x(i+ei->di, j+ei->dj, k+ei->dk,
					f1, f2, iso_val, pt[e], norm[e], level);
				break;
		
			case 1:
				interpRect3Dpts_y(i+ei->di, j+ei->dj, k+ei->dk,
					f1, f2, iso_val, pt[e], norm[e], level);
				break;
		
			case 2:
				interpRect3Dpts_z(i+ei->di, j+ei->dj, k+ei->dk,
					f1, f2, iso_val, pt[e], norm[e], level);
				break;
		}
	}
	return cubeedges[code][0];
}

int Octree::cell_comp_in(int oc_id, int level, float pt[12][3], float norm[12][3])
{
	float val[8], f1, f2;
	int code, e, i, j, k, edge;

	getCellValues(oc_id,level,val);

	code = 0;
	
	if (val[0] > iso_val_in) code |= 0x01;
	if (val[1] > iso_val_in) code |= 0x02;
	if (val[2] > iso_val_in) code |= 0x04;
	if (val[3] > iso_val_in) code |= 0x08;
	if (val[4] > iso_val_in) code |= 0x10;
	if (val[5] > iso_val_in) code |= 0x20;
	if (val[6] > iso_val_in) code |= 0x40;
	if (val[7] > iso_val_in) code |= 0x80;

	assert(code!=0 && code!=255);
		
	octcell2xyz(oc_id,i,j,k,level);
		
	for (e=0 ; e < cubeedges[code][0] ; e++) {
		edge = cubeedges[code][1+e];
		EdgeInfo *ei = &edgeinfo[edge];
		
		f1=val[ei->d1];
		f2=val[ei->d2];
	
		switch (ei->dir) {
			case 0:
				interpRect3Dpts_x(i+ei->di, j+ei->dj, k+ei->dk,
					f1, f2, iso_val_in, pt[e], norm[e], level);
				break;
		
			case 1:
				interpRect3Dpts_y(i+ei->di, j+ei->dj, k+ei->dk,
					f1, f2, iso_val_in, pt[e], norm[e], level);
				break;
		
			case 2:
				interpRect3Dpts_z(i+ei->di, j+ei->dj, k+ei->dk,
					f1, f2, iso_val_in, pt[e], norm[e], level);
				break;
		}
	}
	return cubeedges[code][0];
}

void Octree::clear (double* a, double* b, double* c)
{
	int i;
	for (i=0;i<3;i++) {
		a[i]=0;
		b[i]=0;
		c[i]=0;
	}

}

void Octree::clear (double* a)
{
	int i;
	for (i=0;i<3;i++) {
		a[i]=0;		
	}

}

void Octree::get_qef(int child_idx,double* temp_sigma_ni_2,double* temp_sigma_ni_2_pi,double* temp_sigma_ni_2_pi_2)
{
	if (qef_array[child_idx]) {
		temp_sigma_ni_2[0]=qef_array[child_idx][0];
		temp_sigma_ni_2[1]=qef_array[child_idx][1];
		temp_sigma_ni_2[2]=qef_array[child_idx][2];
		
		temp_sigma_ni_2_pi[0]=qef_array[child_idx][3];
		temp_sigma_ni_2_pi[1]=qef_array[child_idx][4];
		temp_sigma_ni_2_pi[2]=qef_array[child_idx][5];
		
		temp_sigma_ni_2_pi_2[0]=qef_array[child_idx][6];
		temp_sigma_ni_2_pi_2[1]=qef_array[child_idx][7];
		temp_sigma_ni_2_pi_2[2]=qef_array[child_idx][8];
	} else {
		//assert(qef_array[child_idx]!=NULL);
		temp_sigma_ni_2[0]=0;
		temp_sigma_ni_2[1]=0;
		temp_sigma_ni_2[2]=0;
		
		temp_sigma_ni_2_pi[0]=0;
		temp_sigma_ni_2_pi[1]=0;
		temp_sigma_ni_2_pi[2]=0;
		
		temp_sigma_ni_2_pi_2[0]=0;
		temp_sigma_ni_2_pi_2[1]=0;
		temp_sigma_ni_2_pi_2[2]=0;
	}

}

void Octree::get_qef_in(int child_idx,double* temp_sigma_ni_2,double* temp_sigma_ni_2_pi,double* temp_sigma_ni_2_pi_2)
{
	if (qef_array_in[child_idx]) {
		temp_sigma_ni_2[0]=qef_array_in[child_idx][0];
		temp_sigma_ni_2[1]=qef_array_in[child_idx][1];
		temp_sigma_ni_2[2]=qef_array_in[child_idx][2];
		
		temp_sigma_ni_2_pi[0]=qef_array_in[child_idx][3];
		temp_sigma_ni_2_pi[1]=qef_array_in[child_idx][4];
		temp_sigma_ni_2_pi[2]=qef_array_in[child_idx][5];
		
		temp_sigma_ni_2_pi_2[0]=qef_array_in[child_idx][6];
		temp_sigma_ni_2_pi_2[1]=qef_array_in[child_idx][7];
		temp_sigma_ni_2_pi_2[2]=qef_array_in[child_idx][8];
	} else {
		//assert(qef_array_in[child_idx]!=NULL);
		temp_sigma_ni_2[0]=0;
		temp_sigma_ni_2[1]=0;
		temp_sigma_ni_2[2]=0;
		
		temp_sigma_ni_2_pi[0]=0;
		temp_sigma_ni_2_pi[1]=0;
		temp_sigma_ni_2_pi[2]=0;
		
		temp_sigma_ni_2_pi_2[0]=0;
		temp_sigma_ni_2_pi_2[1]=0;
		temp_sigma_ni_2_pi_2[2]=0;
	}

}


void Octree::compute_qef()
{
	float norm[12][3], p[12][3];
	double sigma_ni_2[3], sigma_ni_2_pi[3], sigma_ni_2_pi_2[3], qef;
	double temp_sigma_ni_2[3], temp_sigma_ni_2_pi[3], temp_sigma_ni_2_pi_2[3],solution[3];
	int child_idx, intersect_num;
	int oc_id, level, i, x,y,z;
	
	for (oc_id=level_id[oct_depth];oc_id<level_id[oct_depth+1];oc_id++) {
		if (is_skipcell(oc_id)) continue;
		level=get_level(oc_id);
		octcell2xyz(oc_id,x,y,z,level);

		clear(sigma_ni_2,sigma_ni_2_pi,sigma_ni_2_pi_2);
		intersect_num = cell_comp(oc_id, level, p, norm);

		int j;
		for (j=0;j<intersect_num;j++) {
			for (i=0;i<3;i++) {
				sigma_ni_2[i]      += norm[j][i]*norm[j][i];
				sigma_ni_2_pi[i]   += norm[j][i]*norm[j][i]*p[j][i];
				sigma_ni_2_pi_2[i] += norm[j][i]*norm[j][i]*p[j][i]*p[j][i];
			}

		}

		for (i=0;i<3;i++)
			solution[i] = sigma_ni_2_pi[i] / sigma_ni_2[i] ;
/*
		// Jessica begin
		int cell_size = (dim[0]-1)/(1<<level);
		if (!(solution[0] > x*cell_size && solution[0] < (x+1)*cell_size) ||
			!(solution[1] > y*cell_size && solution[1] < (y+1)*cell_size) ||
			!(solution[2] > z*cell_size && solution[2] < (z+1)*cell_size)) {
			solution[0] = 0.0;		solution[1] = 0.0;		solution[2] = 0.0;
			for (i = 0; i < intersect_num; i++) {
				solution[0] += p[i][0] / intersect_num;
				solution[1] += p[i][1] / intersect_num;
				solution[2] += p[i][2] / intersect_num;
			}
		}
		// Jessica end
*/
		qef=0;
		for (i=0;i<3;i++)
			qef = sigma_ni_2_pi_2[i] - (sigma_ni_2_pi[i]*sigma_ni_2_pi[i]) / sigma_ni_2[i] ;

		put_qef(oc_id, sigma_ni_2, sigma_ni_2_pi, sigma_ni_2_pi_2, solution, qef);
	}

	int l;
	for (l=oct_depth-1;l>=0;l--) {
		for (oc_id=level_id[l];oc_id<level_id[l+1];oc_id++) {
			level=l;

			if (oct_array[oc_id].refine_flag==0) continue;

			clear (temp_sigma_ni_2,temp_sigma_ni_2_pi,temp_sigma_ni_2_pi_2);
			clear (sigma_ni_2,sigma_ni_2_pi,sigma_ni_2_pi_2);
			clear (solution);

			qef=0;
			
			//child(oc_id,level,cid);
			int k;
			for (k=0;k<8;k++) {
				child_idx=child(oc_id,level,k);
				//if (oct_array[child_idx].refine_flag==0) continue;
				if (is_skipcell(child_idx)) continue;
				get_qef(child_idx,temp_sigma_ni_2,temp_sigma_ni_2_pi,temp_sigma_ni_2_pi_2);
				int i;
				for (i=0;i<3;i++) {
					sigma_ni_2[i]+=temp_sigma_ni_2[i];
					sigma_ni_2_pi[i]+=temp_sigma_ni_2_pi[i];
					sigma_ni_2_pi_2[i]+=temp_sigma_ni_2_pi_2[i];
				}
			}
			for (i=0;i<3;i++)
				solution[i] = sigma_ni_2_pi[i] / sigma_ni_2[i];
			qef=0;
			for (i=0;i<3;i++) 
				qef = sigma_ni_2_pi_2[i] - (sigma_ni_2_pi[i]*sigma_ni_2_pi[i]) / sigma_ni_2[i] ;
			put_qef(oc_id,sigma_ni_2, sigma_ni_2_pi, sigma_ni_2_pi_2, solution, qef);

		}
	}
	
}

void Octree::compute_qef_interval()
{
	float norm[12][3], p[12][3];
	double sigma_ni_2[3], sigma_ni_2_pi[3], sigma_ni_2_pi_2[3], qef;
	double temp_sigma_ni_2[3], temp_sigma_ni_2_pi[3], temp_sigma_ni_2_pi_2[3],solution[3];
	int child_idx, intersect_num, oc_id;
	int level, i, k, x, y, z;
	
	for (oc_id = level_id[oct_depth]; oc_id < level_id[oct_depth+1]; oc_id++) {
		if(is_skipcell_interval(oc_id)) continue;
		level = get_level(oc_id);
		octcell2xyz(oc_id, x, y, z, level);

		clear(sigma_ni_2, sigma_ni_2_pi, sigma_ni_2_pi_2);
		
		if(! is_skipcell(oc_id)) {
			intersect_num = cell_comp(oc_id, level, p, norm);
			int j;
			for (j=0;j<intersect_num;j++) {
				for (i=0;i<3;i++) {
					sigma_ni_2[i]      += norm[j][i]*norm[j][i];
					sigma_ni_2_pi[i]   += norm[j][i]*norm[j][i]*p[j][i];
					sigma_ni_2_pi_2[i] += norm[j][i]*norm[j][i]*p[j][i]*p[j][i];
				}

			}

			for (i=0;i<3;i++)
				solution[i] = sigma_ni_2_pi[i] / sigma_ni_2[i] ;

			qef=0;
			for (i=0;i<3;i++)
				qef = sigma_ni_2_pi_2[i] - (sigma_ni_2_pi[i]*sigma_ni_2_pi[i]) / sigma_ni_2[i] ;

			put_qef(oc_id, sigma_ni_2, sigma_ni_2_pi, sigma_ni_2_pi_2, solution, qef);
		}

		clear(sigma_ni_2, sigma_ni_2_pi, sigma_ni_2_pi_2);

		if(! is_skipcell_in(oc_id)) {
			intersect_num = cell_comp_in(oc_id, level, p, norm);
			int j;
			for (j=0;j<intersect_num;j++) {
				for (i=0;i<3;i++) {
					sigma_ni_2[i]      += norm[j][i]*norm[j][i];
					sigma_ni_2_pi[i]   += norm[j][i]*norm[j][i]*p[j][i];
					sigma_ni_2_pi_2[i] += norm[j][i]*norm[j][i]*p[j][i]*p[j][i];
				}

			}

			for (i=0;i<3;i++)
				solution[i] = sigma_ni_2_pi[i] / sigma_ni_2[i] ;

			qef=0;
			for (i=0;i<3;i++)
				qef = sigma_ni_2_pi_2[i] - (sigma_ni_2_pi[i]*sigma_ni_2_pi[i]) / sigma_ni_2[i] ;

			put_qef_in(oc_id, sigma_ni_2, sigma_ni_2_pi, sigma_ni_2_pi_2, solution, qef);
		}
	}

	int l;
	for (l=oct_depth-1;l>=0;l--) {
		for (oc_id = level_id[l]; oc_id < level_id[l+1]; oc_id++) {
			level=l;

			if (oct_array[oc_id].refine_flag == 0) continue;

			clear (temp_sigma_ni_2,temp_sigma_ni_2_pi,temp_sigma_ni_2_pi_2);
			clear (sigma_ni_2,sigma_ni_2_pi,sigma_ni_2_pi_2);
			clear (solution);

			//child(oc_id,level,cid);
			for (k=0;k<8;k++) {
				child_idx = child(oc_id,level,k);
				//if (oct_array[child_idx].refine_flag==0) continue;
				if (is_skipcell(child_idx)) continue;
				get_qef(child_idx,temp_sigma_ni_2,temp_sigma_ni_2_pi,temp_sigma_ni_2_pi_2);
				for (i=0;i<3;i++) {
					sigma_ni_2[i]+=temp_sigma_ni_2[i];
					sigma_ni_2_pi[i]+=temp_sigma_ni_2_pi[i];
					sigma_ni_2_pi_2[i]+=temp_sigma_ni_2_pi_2[i];
				}
			}
			for (i=0;i<3;i++)
				solution[i] = sigma_ni_2_pi[i] / sigma_ni_2[i];
			qef=0;
			for (i=0;i<3;i++) 
				qef = sigma_ni_2_pi_2[i] - (sigma_ni_2_pi[i]*sigma_ni_2_pi[i]) / sigma_ni_2[i];
			put_qef(oc_id,sigma_ni_2, sigma_ni_2_pi, sigma_ni_2_pi_2, solution, qef);

			clear (temp_sigma_ni_2,temp_sigma_ni_2_pi,temp_sigma_ni_2_pi_2);
			clear (sigma_ni_2,sigma_ni_2_pi,sigma_ni_2_pi_2);
			clear (solution);

			//child(oc_id,level,cid);
			for (k = 0; k < 8; k++) {
				child_idx = child(oc_id,level,k);
				//if (oct_array[child_idx].refine_flag==0) continue;
				if (is_skipcell_in(child_idx)) continue;
				get_qef_in(child_idx,temp_sigma_ni_2,temp_sigma_ni_2_pi,temp_sigma_ni_2_pi_2);
				for (i=0;i<3;i++) {
					sigma_ni_2[i]+=temp_sigma_ni_2[i];
					sigma_ni_2_pi[i]+=temp_sigma_ni_2_pi[i];
					sigma_ni_2_pi_2[i]+=temp_sigma_ni_2_pi_2[i];
				}
			}
			for (i=0;i<3;i++)
				solution[i] = sigma_ni_2_pi[i] / sigma_ni_2[i];
			qef=0;
			for (i=0;i<3;i++) 
				qef = sigma_ni_2_pi_2[i] - (sigma_ni_2_pi[i]*sigma_ni_2_pi[i]) / sigma_ni_2[i];
			put_qef_in(oc_id,sigma_ni_2, sigma_ni_2_pi, sigma_ni_2_pi_2, solution, qef);

		}
	}
	
}
/*
void Octree::compute_qef_interval()
{
	float norm[12][3], p[12][3];
	double sigma_ni_2[3], sigma_ni_2_pi[3], sigma_ni_2_pi_2[3], qef;
	double temp_sigma_ni_2[3], temp_sigma_ni_2_pi[3], temp_sigma_ni_2_pi_2[3],solution[3];
	int child_idx, intersect_num, oc_id;
	int level, i, x, y, z, cell_size;
	
	for (oc_id=level_id[oct_depth];oc_id<level_id[oct_depth+1];oc_id++) {
		if (is_skipcell_interval(oc_id)) continue;
		level=get_level(oc_id);
		octcell2xyz(oc_id,x,y,z,level);

		clear(sigma_ni_2,sigma_ni_2_pi,sigma_ni_2_pi_2);
		if(! is_skipcell(oc_id))
			intersect_num = cell_comp(oc_id, level, p, norm);
		else
			intersect_num = cell_comp_in(oc_id, level, p, norm);
		for (int j=0;j<intersect_num;j++) {
			for (i=0;i<3;i++) {
				sigma_ni_2[i]      += norm[j][i]*norm[j][i];
				sigma_ni_2_pi[i]   += norm[j][i]*norm[j][i]*p[j][i];
				sigma_ni_2_pi_2[i] += norm[j][i]*norm[j][i]*p[j][i]*p[j][i];
			}

		}

		for (i=0;i<3;i++)
			solution[i] = sigma_ni_2_pi[i] / sigma_ni_2[i] ;

		qef=0;
		for (i=0;i<3;i++)
			qef = sigma_ni_2_pi_2[i] - (sigma_ni_2_pi[i]*sigma_ni_2_pi[i]) / sigma_ni_2[i] ;

		put_qef(oc_id, sigma_ni_2, sigma_ni_2_pi, sigma_ni_2_pi_2, solution, qef);
	}

	for (int l=oct_depth-1;l>=0;l--) {
		for (oc_id=level_id[l];oc_id<level_id[l+1];oc_id++) {
			level=l;

			if (oct_array[oc_id].refine_flag==0) continue;

			clear (temp_sigma_ni_2,temp_sigma_ni_2_pi,temp_sigma_ni_2_pi_2);
			clear (sigma_ni_2,sigma_ni_2_pi,sigma_ni_2_pi_2);
			clear (solution);

			qef=0;
			
			//child(oc_id,level,cid);
			for (int k=0;k<8;k++) {
				child_idx=child(oc_id,level,k);
				//if (oct_array[child_idx].refine_flag==0) continue;
				if (is_skipcell_interval(child_idx)) continue;
				get_qef(child_idx,temp_sigma_ni_2,temp_sigma_ni_2_pi,temp_sigma_ni_2_pi_2);
				for (int i=0;i<3;i++) {
					sigma_ni_2[i]+=temp_sigma_ni_2[i];
					sigma_ni_2_pi[i]+=temp_sigma_ni_2_pi[i];
					sigma_ni_2_pi_2[i]+=temp_sigma_ni_2_pi_2[i];
				}
			}
			for (i=0;i<3;i++)
				solution[i] = sigma_ni_2_pi[i] / sigma_ni_2[i];
			qef=0;
			for (i=0;i<3;i++) 
				qef = sigma_ni_2_pi_2[i] - (sigma_ni_2_pi[i]*sigma_ni_2_pi[i]) / sigma_ni_2[i] ;
			put_qef(oc_id,sigma_ni_2, sigma_ni_2_pi, sigma_ni_2_pi_2, solution, qef);

		}
	}
	
}
*/
void Octree::put_qef(int oc_id, double* sigma_ni_2, double* sigma_ni_2_pi,double* sigma_ni_2_pi_2,double* solution,double qef)
{
	//assert(qef_array[oc_id]==NULL);

  if (!qef_array[oc_id]) {
    //qef_array[oc_id]=(double*)malloc(sizeof(double)*(3+3+3+3+1));
    qef_array[oc_id].reset(new double[3+3+3+3+1]);
  }

	qef_array[oc_id][0]=sigma_ni_2[0];
	qef_array[oc_id][1]=sigma_ni_2[1];
	qef_array[oc_id][2]=sigma_ni_2[2];

	qef_array[oc_id][3]=sigma_ni_2_pi[0];
	qef_array[oc_id][4]=sigma_ni_2_pi[1];
	qef_array[oc_id][5]=sigma_ni_2_pi[2];

	qef_array[oc_id][6]=sigma_ni_2_pi_2[0];
	qef_array[oc_id][7]=sigma_ni_2_pi_2[1];
	qef_array[oc_id][8]=sigma_ni_2_pi_2[2];

	qef_array[oc_id][9]=solution[0];
	qef_array[oc_id][10]=solution[1];
	qef_array[oc_id][11]=solution[2];

	qef_array[oc_id][12]=qef;
}

void Octree::put_qef_in(int oc_id, double* sigma_ni_2, double* sigma_ni_2_pi,double* sigma_ni_2_pi_2,double* solution,double qef)
{
	//assert(qef_array_in[oc_id]==NULL);

  if (!qef_array_in[oc_id]) {
    //qef_array_in[oc_id]=(double*)malloc(sizeof(double)*(3+3+3+3+1));
    qef_array_in[oc_id].reset(new double[3+3+3+3+1]);
  }

	qef_array_in[oc_id][0]=sigma_ni_2[0];
	qef_array_in[oc_id][1]=sigma_ni_2[1];
	qef_array_in[oc_id][2]=sigma_ni_2[2];

	qef_array_in[oc_id][3]=sigma_ni_2_pi[0];
	qef_array_in[oc_id][4]=sigma_ni_2_pi[1];
	qef_array_in[oc_id][5]=sigma_ni_2_pi[2];

	qef_array_in[oc_id][6]=sigma_ni_2_pi_2[0];
	qef_array_in[oc_id][7]=sigma_ni_2_pi_2[1];
	qef_array_in[oc_id][8]=sigma_ni_2_pi_2[2];

	qef_array_in[oc_id][9]=solution[0];
	qef_array_in[oc_id][10]=solution[1];
	qef_array_in[oc_id][11]=solution[2];

	qef_array_in[oc_id][12]=qef;
}

float Octree::get_err(int oc_id)
{
	//assert(qef_array[oc_id]!=NULL) ;

	if(is_skipcell(oc_id) == 0) {
	  //if(qef_array[oc_id] != NULL)
	  if(qef_array[oc_id])
	    return (float)qef_array[oc_id][12];
	  else
	    return -1.0;
	}
	else {
	  //if(qef_array_in[oc_id] != NULL)
	  if(qef_array_in[oc_id])
	    return (float)qef_array_in[oc_id][12];
	  else
	    return -1.0;
	}
}
/*
float Octree::get_err_grad(int oc_id)
{
	int level, cell_size, xx, yy, zz, ii, jj, kk;
	float val[8], val_new[13], x, y, z;
	float func_0, func_1, f_x, f_y, f_z, sum, sum_grad, temp;

	level = get_level(oc_id) ;
	cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_id, xx, yy, zz, level);
	getCellValues(oc_id, level, val);

	bool my_bool = (val[0] < iso_val && val[1] < iso_val && val[2] < iso_val && val[3] < iso_val
		&& val[4] < iso_val && val[5] < iso_val && val[6] < iso_val && val[7] < iso_val)
		|| (val[0] > iso_val && val[1] > iso_val && val[2] > iso_val && val[3] > iso_val
		&& val[4] > iso_val && val[5] > iso_val && val[6] > iso_val && val[7] > iso_val);
	if(flag_type > 3) {
		bool my_bool_in = minmax[oc_id].min > iso_val || minmax[oc_id].max < iso_val_in ||
			(minmax[oc_id].min > iso_val_in && minmax[oc_id].max < iso_val);
		my_bool = my_bool_in;
	}

	val_new[0] = getValue(xx*cell_size+cell_size/2, yy*cell_size, zz*cell_size);
	val_new[1] = getValue(xx*cell_size+cell_size,	yy*cell_size, zz*cell_size+cell_size/2);
	val_new[2] = getValue(xx*cell_size+cell_size/2, yy*cell_size, zz*cell_size+cell_size);
	val_new[3] = getValue(xx*cell_size,				yy*cell_size, zz*cell_size+cell_size/2);
	val_new[4] = getValue(xx*cell_size+cell_size/2, yy*cell_size+cell_size, zz*cell_size);
	val_new[5] = getValue(xx*cell_size+cell_size,	yy*cell_size+cell_size, zz*cell_size+cell_size/2);
	val_new[6] = getValue(xx*cell_size+cell_size/2, yy*cell_size+cell_size, zz*cell_size+cell_size);
	val_new[7] = getValue(xx*cell_size,				yy*cell_size+cell_size, zz*cell_size+cell_size/2);
	val_new[8] = getValue(xx*cell_size,				yy*cell_size+cell_size/2, zz*cell_size);
	val_new[9] = getValue(xx*cell_size+cell_size,	yy*cell_size+cell_size/2, zz*cell_size);
	val_new[10]= getValue(xx*cell_size,				yy*cell_size+cell_size/2, zz*cell_size+cell_size);
	val_new[11]= getValue(xx*cell_size+cell_size,	yy*cell_size+cell_size/2, zz*cell_size+cell_size);
	val_new[12]= getValue(xx*cell_size+cell_size/2,	yy*cell_size+cell_size/2, zz*cell_size+cell_size/2);

	sum = 0.0;
	sum_grad = 0.0;
	for(ii = 0; ii < 13; ii++) {
		if(ii == 3 || ii == 7 || ii == 8 || ii == 10)		x = 0.0;
		else if(ii == 1 || ii == 5 || ii == 9 || ii == 11)	x = 1.0;
		else	x = 0.5;

		if(ii < 4)		y = 0.0;
		else if(ii < 8)	y = 1.0;
		else	y = 0.5;

		if(ii == 0 || ii == 4 || ii == 8 || ii == 9)		z = 0.0;
		else if(ii == 2 || ii == 6 || ii == 10 || ii == 11)	z = 1.0;
		else	z = 0.5;

		func_0 = (1-x)*(1-y)*(1-z)*val[0] + (1-x)*(1-y)*z*val[3] + (1-x)*y*(1-z)*val[4]
					+ x*(1-y)*(1-z)*val[1] + (1-x)*y*z*val[7] + x*(1-y)*z*val[2]
					+ x*y*(1-z)*val[5] + x*y*z*val[6];
		func_1 = val_new[ii];

		if(func_1 > func_0)	temp = func_1 - func_0;
		else	temp = func_0 - func_1;

		sum += temp;

		f_x = (1-y)*(1-z)*(val[1] - val[0]) + (1-y)*z*(val[2] - val[3])
					+ y*(1-z)*(val[5] - val[4]) + y*z*(val[6] - val[7]);
		f_y = (1-x)*(1-z)*(val[4] - val[0]) + (1-x)*z*(val[7] - val[3])
					+ x*(1-z)*(val[5] - val[1]) + x*z*(val[6] - val[2]);
		f_z = (1-x)*(1-z)*(val[3] - val[0]) + (1-x)*y*(val[7] - val[4])
					+ x*(1-y)*(val[2] - val[1]) + x*y*(val[6] - val[5]);

		sum_grad += sqrt(f_x*f_x + f_y*f_y + f_z*f_z);
	}

	for(ii = 0; ii < 2; ii++) {
		for(jj = 0; jj < 2; jj++) {
			for(kk = 0; kk < 2; kk ++) {
				x = ii; y = jj; z = kk;

				f_x = (1-y)*(1-z)*(val[1] - val[0]) + (1-y)*z*(val[2] - val[3])
							+ y*(1-z)*(val[5] - val[4]) + y*z*(val[6] - val[7]);
				f_y = (1-x)*(1-z)*(val[4] - val[0]) + (1-x)*z*(val[7] - val[3])
							+x*(1-z)*(val[5] - val[1]) + x*z*(val[6] - val[2]);
				f_z = (1-x)*(1-y)*(val[3] - val[0]) + (1-x)*y*(val[7] - val[4])
							+x*(1-y)*(val[2] - val[1]) + x*y*(val[6] - val[5]);

				sum_grad += sqrt(f_x*f_x + f_y*f_y + f_z*f_z);
			}
		}
	}

	if(!my_bool) {
		return (sum/sum_grad);
	}
	else
		return (-1.0);

}

*/
float Octree::get_err_grad_test(int oc_id)
{
	int level, cell_size, xx, yy, zz;
	float val[8], val_new[19], x, y, z;
	float func_0, func_1, f_x, f_y, f_z, sum_grad, temp;

	level = get_level(oc_id) ;
	cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_id, xx, yy, zz, level);
	getCellValues(oc_id, level, val);

	if(level == oct_depth) 	return (-1.0);

	bool my_bool = (val[0] < iso_val && val[1] < iso_val && val[2] < iso_val && val[3] < iso_val
		&& val[4] < iso_val && val[5] < iso_val && val[6] < iso_val && val[7] < iso_val)
		|| (val[0] > iso_val && val[1] > iso_val && val[2] > iso_val && val[3] > iso_val
		&& val[4] > iso_val && val[5] > iso_val && val[6] > iso_val && val[7] > iso_val);
	if(flag_type > 3) {
		bool my_bool_in = minmax[oc_id].min > iso_val || minmax[oc_id].max < iso_val_in ||
			(minmax[oc_id].min > iso_val_in && minmax[oc_id].max < iso_val);
		my_bool = my_bool_in;
	}

	if(is_skipcell(oc_id) == 0) {
	  //if(qef_array[oc_id] != NULL) {
	  if(qef_array[oc_id]) {
			x = qef_array[oc_id][9] / cell_size - xx;
			y = qef_array[oc_id][10] / cell_size - yy;
			z = qef_array[oc_id][11] / cell_size - zz;
		}
	}

	val_new[0] = getValue(xx*cell_size+cell_size/2, yy*cell_size,				zz*cell_size);
	val_new[1] = getValue(xx*cell_size+cell_size,	yy*cell_size,				zz*cell_size+cell_size/2);
	val_new[2] = getValue(xx*cell_size+cell_size/2, yy*cell_size,				zz*cell_size+cell_size);
	val_new[3] = getValue(xx*cell_size,				yy*cell_size,				zz*cell_size+cell_size/2);
	val_new[4] = getValue(xx*cell_size+cell_size/2, yy*cell_size+cell_size,		zz*cell_size);
	val_new[5] = getValue(xx*cell_size+cell_size,	yy*cell_size+cell_size,		zz*cell_size+cell_size/2);
	val_new[6] = getValue(xx*cell_size+cell_size/2, yy*cell_size+cell_size,		zz*cell_size+cell_size);
	val_new[7] = getValue(xx*cell_size,				yy*cell_size+cell_size,		zz*cell_size+cell_size/2);
	val_new[8] = getValue(xx*cell_size,				yy*cell_size+cell_size/2,	zz*cell_size);
	val_new[9] = getValue(xx*cell_size+cell_size,	yy*cell_size+cell_size/2,	zz*cell_size);
	val_new[10]= getValue(xx*cell_size,				yy*cell_size+cell_size/2,	zz*cell_size+cell_size);
	val_new[11]= getValue(xx*cell_size+cell_size,	yy*cell_size+cell_size/2,	zz*cell_size+cell_size);
	val_new[12]= getValue(xx*cell_size+cell_size/2,	yy*cell_size+cell_size/2,	zz*cell_size+cell_size/2);
	val_new[13]= getValue(xx*cell_size,				yy*cell_size+cell_size/2,	zz*cell_size+cell_size/2);
	val_new[14]= getValue(xx*cell_size+cell_size,	yy*cell_size+cell_size/2,	zz*cell_size+cell_size/2);
	val_new[15]= getValue(xx*cell_size+cell_size/2,	yy*cell_size,				zz*cell_size+cell_size/2);
	val_new[16]= getValue(xx*cell_size+cell_size/2,	yy*cell_size+cell_size,		zz*cell_size+cell_size/2);
	val_new[17]= getValue(xx*cell_size+cell_size/2,	yy*cell_size+cell_size/2,	zz*cell_size);
	val_new[18]= getValue(xx*cell_size+cell_size/2,	yy*cell_size+cell_size/2,	zz*cell_size+cell_size);

	sum_grad = 0.0;

	func_0 = (1-x)*(1-y)*(1-z)*val[0] + (1-x)*(1-y)*z*val[3] + (1-x)*y*(1-z)*val[4]
				+ x*(1-y)*(1-z)*val[1] + (1-x)*y*z*val[7] + x*(1-y)*z*val[2]
				+ x*y*(1-z)*val[5] + x*y*z*val[6];
	//func_1 = val_new[ii];

	int new_oc_id;
	if(x < 0.5 && y < 0.5 && z < 0.5) {
		new_oc_id = xyz2octcell(xx*2, yy*2, zz*2, level+1);
		getCellValues(new_oc_id, level+1, val);
		x = x*2;	y = y*2;	z = z*2;
	}
	else if(x > 0.5 && y < 0.5 && z < 0.5) {
		new_oc_id = xyz2octcell(xx*2+1, yy*2, zz*2, level+1);
		getCellValues(new_oc_id, level+1, val);
		x = x*2-1;	y = y*2;	z = z*2;
	}
	else if(x < 0.5 && y > 0.5 && z < 0.5) {
		new_oc_id = xyz2octcell(xx*2, yy*2+1, zz*2, level+1);
		getCellValues(new_oc_id, level+1, val);
		x = x*2;	y = y*2-1;	z = z*2;
	}
	else if(x > 0.5 && y > 0.5 && z < 0.5) {
		new_oc_id = xyz2octcell(xx*2+1, yy*2+1, zz*2, level+1);
		getCellValues(new_oc_id, level+1, val);
		x = x*2-1;	y = y*2-1;	z = z*2;
	}
	else if(x < 0.5 && y < 0.5 && z > 0.5) {
		new_oc_id = xyz2octcell(xx*2, yy*2, zz*2+1, level+1);
		getCellValues(new_oc_id, level+1, val);
		x = x*2;	y = y*2;	z = z*2-1;
	}
	else if(x > 0.5 && y < 0.5 && z > 0.5) {
		new_oc_id = xyz2octcell(xx*2+1, yy*2, zz*2+1, level+1);
		getCellValues(new_oc_id, level+1, val);
		x = x*2-1;	y = y*2;	z = z*2-1;
	}
	else if(x < 0.5 && y > 0.5 && z > 0.5) {
		new_oc_id = xyz2octcell(xx*2, yy*2+1, zz*2+1, level+1);
		getCellValues(new_oc_id, level+1, val);
		x = x*2;	y = y*2-1;	z = z*2-1;
	}
	else  {	//if(x > 0.5 && y > 0.5 && z > 0.5)
		new_oc_id = xyz2octcell(xx*2+1, yy*2+1, zz*2+1, level+1);
		getCellValues(new_oc_id, level+1, val);
		x = x*2-1;	y = y*2-1;	z = z*2-1;
	}

	func_1 = (1-x)*(1-y)*(1-z)*val[0] + (1-x)*(1-y)*z*val[3] + (1-x)*y*(1-z)*val[4]
				+ x*(1-y)*(1-z)*val[1] + (1-x)*y*z*val[7] + x*(1-y)*z*val[2]
				+ x*y*(1-z)*val[5] + x*y*z*val[6];

	//temp = fabs(func_0);
	if(func_1 > func_0)	temp = func_1 - func_0;
	else	temp = func_0 - func_1;

	f_x = (1-y)*(1-z)*(val[1] - val[0]) + (1-y)*z*(val[2] - val[3])
				+ y*(1-z)*(val[5] - val[4]) + y*z*(val[6] - val[7]);
	f_y = (1-x)*(1-z)*(val[4] - val[0]) + (1-x)*z*(val[7] - val[3])
				+ x*(1-z)*(val[5] - val[1]) + x*z*(val[6] - val[2]);
	f_z = (1-x)*(1-y)*(val[3] - val[0]) + (1-x)*y*(val[7] - val[4])
				+ x*(1-y)*(val[2] - val[1]) + x*y*(val[6] - val[5]);

	sum_grad = temp / sqrt(f_x*f_x + f_y*f_y + f_z*f_z);

	if(!my_bool) {
		return (sum_grad);
	}
	else
		return (-1.0);

}

float Octree::get_err_grad(int oc_id)
{
	int level, cell_size, xx, yy, zz, ii;
	float val[8], val_new[19], x, y, z;
	float func_0, func_1, f_x, f_y, f_z, sum_grad, temp;

	level = get_level(oc_id) ;
	cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_id, xx, yy, zz, level);
	getCellValues(oc_id, level, val);

	bool my_bool = (val[0] < iso_val && val[1] < iso_val && val[2] < iso_val && val[3] < iso_val
		&& val[4] < iso_val && val[5] < iso_val && val[6] < iso_val && val[7] < iso_val)
		|| (val[0] > iso_val && val[1] > iso_val && val[2] > iso_val && val[3] > iso_val
		&& val[4] > iso_val && val[5] > iso_val && val[6] > iso_val && val[7] > iso_val);
	if(flag_type > 3) {
		bool my_bool_in = minmax[oc_id].min > iso_val || minmax[oc_id].max < iso_val_in ||
			(minmax[oc_id].min > iso_val_in && minmax[oc_id].max < iso_val);
		my_bool = my_bool_in;
	}

	val_new[0] = getValue(xx*cell_size+cell_size/2, yy*cell_size,				zz*cell_size);
	val_new[1] = getValue(xx*cell_size+cell_size,	yy*cell_size,				zz*cell_size+cell_size/2);
	val_new[2] = getValue(xx*cell_size+cell_size/2, yy*cell_size,				zz*cell_size+cell_size);
	val_new[3] = getValue(xx*cell_size,				yy*cell_size,				zz*cell_size+cell_size/2);
	val_new[4] = getValue(xx*cell_size+cell_size/2, yy*cell_size+cell_size,		zz*cell_size);
	val_new[5] = getValue(xx*cell_size+cell_size,	yy*cell_size+cell_size,		zz*cell_size+cell_size/2);
	val_new[6] = getValue(xx*cell_size+cell_size/2, yy*cell_size+cell_size,		zz*cell_size+cell_size);
	val_new[7] = getValue(xx*cell_size,				yy*cell_size+cell_size,		zz*cell_size+cell_size/2);
	val_new[8] = getValue(xx*cell_size,				yy*cell_size+cell_size/2,	zz*cell_size);
	val_new[9] = getValue(xx*cell_size+cell_size,	yy*cell_size+cell_size/2,	zz*cell_size);
	val_new[10]= getValue(xx*cell_size,				yy*cell_size+cell_size/2,	zz*cell_size+cell_size);
	val_new[11]= getValue(xx*cell_size+cell_size,	yy*cell_size+cell_size/2,	zz*cell_size+cell_size);
	val_new[12]= getValue(xx*cell_size+cell_size/2,	yy*cell_size+cell_size/2,	zz*cell_size+cell_size/2);
	val_new[13]= getValue(xx*cell_size,				yy*cell_size+cell_size/2,	zz*cell_size+cell_size/2);
	val_new[14]= getValue(xx*cell_size+cell_size,	yy*cell_size+cell_size/2,	zz*cell_size+cell_size/2);
	val_new[15]= getValue(xx*cell_size+cell_size/2,	yy*cell_size,				zz*cell_size+cell_size/2);
	val_new[16]= getValue(xx*cell_size+cell_size/2,	yy*cell_size+cell_size,		zz*cell_size+cell_size/2);
	val_new[17]= getValue(xx*cell_size+cell_size/2,	yy*cell_size+cell_size/2,	zz*cell_size);
	val_new[18]= getValue(xx*cell_size+cell_size/2,	yy*cell_size+cell_size/2,	zz*cell_size+cell_size);

	sum_grad = 0.0;
	for(ii = 0; ii < 19; ii++) {
		if(ii == 3 || ii == 7 || ii == 8 || ii == 10 || ii == 13)		x = 0.0;
		else if(ii == 1 || ii == 5 || ii == 9 || ii == 11 || ii == 14)	x = 1.0;
		else	x = 0.5;

		if(ii < 4 || ii == 15)		y = 0.0;
		else if(ii < 8 || ii == 16)	y = 1.0;
		else	y = 0.5;

		if(ii == 0 || ii == 4 || ii == 8 || ii == 9 || ii == 17)		z = 0.0;
		else if(ii == 2 || ii == 6 || ii == 10 || ii == 11 || ii == 18)	z = 1.0;
		else	z = 0.5;

		func_0 = (1-x)*(1-y)*(1-z)*val[0] + (1-x)*(1-y)*z*val[3] + (1-x)*y*(1-z)*val[4]
					+ x*(1-y)*(1-z)*val[1] + (1-x)*y*z*val[7] + x*(1-y)*z*val[2]
					+ x*y*(1-z)*val[5] + x*y*z*val[6];
		func_1 = val_new[ii];

		if(func_1 > func_0)	temp = func_1 - func_0;
		else	temp = func_0 - func_1;

		f_x = (1-y)*(1-z)*(val[1] - val[0]) + (1-y)*z*(val[2] - val[3])
					+ y*(1-z)*(val[5] - val[4]) + y*z*(val[6] - val[7]);
		f_y = (1-x)*(1-z)*(val[4] - val[0]) + (1-x)*z*(val[7] - val[3])
					+ x*(1-z)*(val[5] - val[1]) + x*z*(val[6] - val[2]);
		f_z = (1-x)*(1-y)*(val[3] - val[0]) + (1-x)*y*(val[7] - val[4])
					+ x*(1-y)*(val[2] - val[1]) + x*y*(val[6] - val[5]);

		sum_grad += temp / sqrt(f_x*f_x + f_y*f_y + f_z*f_z);
	}

	if(!my_bool) {
		return (sum_grad);
	}
	else
		return (-1.0);

}

/*void Octree::traverse_qef(float err_tol)
{
	CellQueue prev_queue,cur_queue;
	int i,oc_id,level;
	int x, y, z, tx, ty, tz, cell_size;
	float r, r_min, temp, center;

	oc_id=0;
	leaf_num=0;
	center = (dim[0]-1.0f)/2.0f;

	memset(oct_array,0,sizeof(octcell)*octcell_num);
	prev_queue.Add(oc_id);
	while (!prev_queue.Empty()) {
		
		//leaf_idx=leaf_num;
		while (prev_queue.Get(oc_id) >= 0)  {
			level=get_level(oc_id);
			
			octcell2xyz(oc_id,x,y,z,level);
			cell_size = (dim[0]-1)/(1<<level);
			tx = x*cell_size;		ty = y*cell_size;		tz = z*cell_size;
			r = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			r_min = r;
			tx = (x+1)*cell_size;	ty = y*cell_size;		tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			tx = x*cell_size;		ty = (y+1)*cell_size;	tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = (y+1)*cell_size;	tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			
			tx = x*cell_size;		ty = y*cell_size;		tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = y*cell_size;		tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			tx = x*cell_size;		ty = (y+1)*cell_size;	tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = (y+1)*cell_size;	tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			tx = x*cell_size;		ty = y*cell_size;		tz = z*cell_size;
			if(tx < center)	tx = (x+1)*cell_size;
			if(ty < center)	ty = (y+1)*cell_size;
			if(tz < center)	tz = (z+1)*cell_size;

			//if (is_skipcell(oc_id)) {
			if(minmax[oc_id].min > iso_val) {
				continue;
			//} else if (level <= 3 ||(get_err_grad_test(oc_id) > err_tol && level != oct_depth)) {

			//} else if (level <= oct_depth-4 || (level <= oct_depth-2 && is_skipcell(oc_id) == 0)
			//	||(get_err_grad(oc_id) > err_tol && level != oct_depth)) {

			//} else if (level < oct_depth) {
			
			// 1a02 65-129
			} else if ((level <= 1 && r > (dim[0]-1)/4.0*1.733) || (level <= 3 && r < (dim[0]-1)/4.0*1.733)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 && r_min < (dim[0]-1)/4.0*1.0) //1.515
				||(get_err_grad(oc_id) > err_tol && level < oct_depth && r < (dim[0]-1)/8.0*1.8-0.2)) { // 14.3
			
			// mAChEH 161/2
			} else if ((level <= 2 && r > 63) || (level <= 3 && r < 40*1.75 && !(abs(tx-center)>41 && abs(ty-center)>41 && abs(tz-center)>41 && is_skipcell(oc_id) == 0)) // 1.4
				|| (level <= oct_depth-2 && is_skipcell(oc_id) >= 0 && minmax[oc_id].max > -0.5 && r < 40*1.75 && abs(tx-center)<40 && abs(ty-center)<40 && abs(tz-center)<40) // r < 47
				||(get_err_grad(oc_id) > err_tol && level < oct_depth && r < 40*1.75 && abs(tx-center)<40 && abs(ty-center)<40 && abs(tz-center)<40)) { // 14.3, 50

			// ion-acc 97-257
			} else if ((level <= 1 && r > (dim[0]-1)/4.0*1.733) || (level <= 2 && r < (dim[0]-1)/4.0*1.733)
				|| (level <= 4 && r < 48*1.5)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 && r_min < 48*1.3) //1.515
				||(get_err_grad(oc_id) > err_tol && level < oct_depth && r < (dim[0]-1)/8.0*2.5)) { // 14.3
			
			// ion-acc 211-257   adaptive cavity
			} else if ((level <= 3 && r > 80*1.733) || (level <= 3 && r < 80*1.733)
				|| (level <= 4 && r < 80*1.5)
				|| (level <= oct_depth-3 && is_skipcell(oc_id) == 0 && r_min < 80*1.42)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 && 
					tx >= 94 && tx <= 172 && ty >= 149 && ty <= 228 && tz >= 81 && tz <=147)	  //211
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 104 && tx <= 162 && ty >= 149 && ty <= 218 && tz >= 91 && tz <=137)) {  //211
			
			// ion-acc 161-257
			} else if ((level <= 2 && r > 120) || (level <= 3 && abs(tx-center)<80 && abs(ty-center)<80 && abs(tz-center)<80)
				|| (level <= 4 && r < 80 && abs(tx-center)<80 && abs(ty-center)<80 && abs(tz-center)<80)
				|| (level <= oct_depth-3 && is_skipcell(oc_id) == 0 && r < 120 && abs(tx-center)<80 && abs(ty-center)<80 && abs(tz-center)<80)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 && 
					tx >= 90 && tx <= 164 && ty >= 135 && ty <= 220 && tz >= 84 && tz <=154)		//161
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					//tx >= 106 && tx <= 148 && ty >= 139 && ty <= 215 && tz >= 98 && tz <=140)) {  //161
					tx >= 100 && tx <= 154 && ty >= 145 && ty <= 210 && tz >= 94 && tz <=144)) { 

			// sphere 97-257  coarse
			} else if ((level <= 1 && r > 64*1.733) || (level <= 2 && r < 64*1.733)
				|| (level <= 4 && r < 48*1.733)
				|| (level <= 5 && r < 48)
				|| (level <= oct_depth-2 && r < 40.1) //35.1
				|| (level <= oct_depth-1 && r < 35.1) //35.1
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && r < 34*1.733)) {
			
			// sphere 97-257  denser
			} else if ((level <= 4 && r > 64) || (level <= 4 && r < 64)
				|| (level <= 4 && r < 48*1.733)
				|| (level <= 5 && r < 48)
				|| (level <= oct_depth-2 && r < 40.1) //35.1
				|| (level <= oct_depth-1 && r < 35.1) //35.1
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && r < 34*1.733)) {
			
			// mAChE4_den_257 
			} else if (level < oct_depth-2 
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 1
					tx >= 42 && tx <= 101 && ty >= 93 && ty <= 145 && tz >= 119 && tz <= 168)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 47 && tx <= 96 && ty >= 98 && ty <= 140 && tz >= 124 && tz <= 163)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 2
					tx >= 45 && tx <= 92 && ty >= 112 && ty <= 161 && tz >= 16 && tz <= 57)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 50 && tx <= 88 && ty >= 117 && ty <= 156 && tz >= 21 && tz <= 52)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 3
					tx >= 99 && tx <= 158 && ty >= 175 && ty <= 223 && tz >= 171 && tz <= 227)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 104 && tx <= 153 && ty >= 180 && ty <= 218 && tz >= 176 && tz <= 222)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 4
					tx >= 133 && tx <= 179 && ty >= 134 && ty <= 182 && tz >= 83 && tz <= 139)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 138 && tx <= 174 && ty >= 139 && ty <= 177 && tz >= 88 && tz <= 135)) {
			
			// mAChE4_den1_161 
			} else if (level <= oct_depth-5 
				|| (level <= oct_depth-3 && is_skipcell(oc_id) == 0 && r < 128.0f)
				|| (level <= oct_depth-4 && is_skipcell(oc_id) == 0 && r > 127.5f &&
					tx >= 43 && tx <= 91 && ty >= 33 && ty <= 82 && tz >= 40 && tz <= 80)
				|| (level <= oct_depth-3 && is_skipcell(oc_id) == 0 &&
					tx >= 58 && tx <= 76 && ty >= 58 && ty <= 77 && tz >= 50 && tz <= 56)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 1
					tx >= 75 && tx <= 115 && ty >= 89 && ty <= 126 && tz >= 117 && tz <= 160)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 80 && tx <= 110 && ty >= 94 && ty <= 121 && tz >= 122 && tz <= 155)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 2
					tx >= 77 && tx <= 113 && ty >= 100 && ty <= 140 && tz >= 54 && tz <= 90)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 82 && tx <= 108 && ty >= 105 && ty <= 135 && tz >= 59 && tz <= 85)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 3
					tx >= 108 && tx <= 147 && ty >= 147 && ty <= 182 && tz >= 158 && tz <= 191)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 113 && tx <= 142 && ty >= 152 && ty <= 177 && tz >= 163 && tz <= 186)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 4
					tx >= 126 && tx <= 165 && ty >= 124 && ty <= 155 && tz >= 97 && tz <= 136)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 131 && tx <= 160 && ty >= 129 && ty <= 150 && tz >= 102 && tz <= 131)) {
			
			// mAChE4_acc_257 
			} else if (level < oct_depth-2 
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 1
					tx >= 61 && tx <= 106 && ty >= 93 && ty <= 149 && tz >= 40 && tz <= 76)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 66 && tx <= 101 && ty >= 98 && ty <= 144 && tz >= 45 && tz <= 71)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 2
					tx >= 66 && tx <= 101 && ty >= 80 && ty <= 124 && tz >= 113 && tz <= 165)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 71 && tx <= 106 && ty >= 85 && ty <= 119 && tz >= 118 && tz <= 160)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 3
					tx >= 103 && tx <= 150 && ty >= 154 && ty <= 187 && tz >= 161 && tz <= 205)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 108 && tx <= 145 && ty >= 159 && ty <= 192 && tz >= 166 && tz <= 200)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 4
					tx >= 132 && tx <= 178 && ty >= 124 && ty <= 172 && tz >= 95 && tz <= 136)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 137 && tx <= 173 && ty >= 129 && ty <= 167 && tz >= 100 && tz <= 131)
					) {
			
			// mAChE4_acc_257_p1 
			} else if (level < oct_depth-2 
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 3
					tx >= 0 && tx <= 105 && ty >= 103 && ty <= 199 && tz >= 117 && tz <= 215)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 5 && tx <= 100 && ty >= 108 && ty <= 194 && tz >= 122 && tz <= 210)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 4
					tx >= 59 && tx <= 161 && ty >= 43 && ty <= 149 && tz >= 0 && tz <= 77)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 64 && tx <= 156 && ty >= 48 && ty <= 144 && tz >= 0 && tz <= 72)) {
			
			// mAChE4_acc_257_p2 
			} else if (level < oct_depth-2 
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&						// cavity 1
					tx >= 47 && tx <= 127 && ty >= 111 && ty <= 213 && tz >= 0 && tz <= 67)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 52 && tx <= 122 && ty >= 116 && ty <= 208 && tz >= 0 && tz <= 62)
				|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 &&					// cavity 2
					tx >= 57 && tx <= 137 && ty >= 85 && ty <= 163 && tz >= 151 && tz <= 245)
				|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
					tx >= 62 && tx <= 132 && ty >= 90 && ty <= 158 && tz >= 156 && tz <= 240)) {

			// 3sphere2_inner
			} else if (level < oct_depth-4 ||
						(level < oct_depth-3 && r > 15.0f) ||
						(level < oct_depth-2 && r > 25.0f) ||
						(level < oct_depth-1 && is_skipcell(oc_id) == 0 && r > 32.0f) ||
						(get_err(oc_id) > err_tol && level < oct_depth-1 && r > 32.0f)) {
			
			// 3sphere3_inner
			} else if (level < oct_depth-5 ||
						(level < oct_depth-4 && r > 15.0f) ||
						(level < oct_depth-3 && r > 25.0f) ||
						(level < oct_depth-2 && is_skipcell(oc_id) == 0 && r > 32.0f) ||
						(get_err(oc_id) > err_tol && level < oct_depth-2 && r > 32.0f)) {
			
			// 3sphere_outer
			} else if (level < oct_depth-4 ||
						(level < oct_depth-3 && r < 105.0f) ||
						(level < oct_depth-2 && r < 86.0f) ||
						(level < oct_depth-1 && r < 78.0f) ||
						(level < oct_depth-1 && is_skipcell(oc_id) == 0 && r < 74.0f) ||
						(get_err(oc_id) > err_tol && level < oct_depth && r < 74.0f)) {
			
			// 1FSS
			//} else if (level < oct_depth-2 
			//	|| (level < oct_depth-1 && is_skipcell(oc_id) == 0 &&					
			//		tx >= 24 && tx <= 76 && ty >= 40 && ty <= 80 && tz >= 39 && tz <= 87)
			//	|| (get_err_grad(oc_id) > err_tol && level < oct_depth && 
			//		tx >= 29 && tx <= 71 && ty >= 45 && ty <= 75 && tz >= 44 && tz <= 82)) {

			// ribosome_sphere 129257
			//} else if ((level <= 3 && r > 64*1.0) || (level <= 5 && r < 64*1.0)
			//	|| (get_err_grad(oc_id) > err_tol && level < oct_depth && r < 64*1.0)) {
			
			// hemoglobin 129
			//} else if (level < oct_depth-2 || (is_skipcell(oc_id) == 0 && level < oct_depth)) {
			//} else if (level < oct_depth-2 || (is_skipcell(oc_id) == 0 && level < oct_depth &&
			//			tx>7 && tx<122 && ty>7 && ty<122 && tz>7 && tz<122)) {

			// spinalcord
			//} else if (level < oct_depth-2 || (is_skipcell(oc_id) == 0 && level < oct_depth &&
			//			r < 163.0f)) {
			// RDV trimer segmented
			//} else if (level < oct_depth || (is_skipcell(oc_id) == 0 && level < oct_depth &&
			//			r < 163.0f)) {
			// bladder
			//} else if (level < oct_depth-3 || (ty > 144.0f && level < oct_depth-1) || 
			//			(ty <= 144.0f && level < oct_depth-2)) {

			// 1I9B_0.1
			//} else if (level < oct_depth-3 || (is_skipcell(oc_id) == 0 && level < oct_depth-1 && r < 63.0f) ||
			//			(get_err(oc_id) > err_tol && level < oct_depth-1 && r < 63.0f)) {

			// Benzhuo
			} else if (level < oct_depth-2 || (is_skipcell(oc_id) == 0 && level < oct_depth-2) ||
				(get_err(oc_id) > err_tol && r < 60.0f && level < oct_depth) ){//||
						//(is_skipcell(oc_id) == 0 && r < 60.0f && level < oct_depth) ||
						//(get_err(oc_id) > err_tol && level < oct_depth)) {

			
			// hemoglobin_sphere 257
			//} else if ((level <= 3 && r > 64*1.4) || (level <= 5 && r < 64*1.4)
			//	|| (get_err_grad(oc_id) > err_tol && level < oct_depth-1 && r < 64*1.4)) {
			
				cur_queue.Add(oc_id);
				oct_array[oc_id].refine_flag=1;
			} else {
				cut_array[leaf_num++]=oc_id;
			}
		}
		
		while (cur_queue.Get(oc_id) >= 0) {
			level=get_level(oc_id);
			for (i=0; i<8; i++)
				prev_queue.Add(child(oc_id,level,i));
		}
		
	}
}*/

/*void Octree::traverse_qef(float err_tol)
{
	CellQueue prev_queue,cur_queue;
	int i,oc_id,level, start_level, end_level;
	int x, y, z, tx, ty, tz, cell_size;
	float r, r_min, temp, center;

	oc_id=0;
	leaf_num=0;
	center = (dim[0]-1.0f)/2.0f;

	start_level = oct_depth - 3;	// - 3
	end_level = oct_depth;
	if(flag_type == 2) {
		start_level = oct_depth - 3;//3
		end_level = start_level + 1;
	}
	if(flag_type == 3) {
		start_level = oct_depth - 3;//3
		end_level = start_level + 1;
	}

	memset(oct_array,0,sizeof(octcell)*octcell_num);
	prev_queue.Add(oc_id);
	while (!prev_queue.Empty()) {
		
		//leaf_idx=leaf_num;
		while (prev_queue.Get(oc_id) >= 0)  {
			level=get_level(oc_id);
			
			octcell2xyz(oc_id,x,y,z,level);
			cell_size = (dim[0]-1)/(1<<level);
			tx = x*cell_size;		ty = y*cell_size;		tz = z*cell_size;
			r = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			r_min = r;
			tx = (x+1)*cell_size;	ty = y*cell_size;		tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			tx = x*cell_size;		ty = (y+1)*cell_size;	tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = (y+1)*cell_size;	tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			
			tx = x*cell_size;		ty = y*cell_size;		tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = y*cell_size;		tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			tx = x*cell_size;		ty = (y+1)*cell_size;	tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = (y+1)*cell_size;	tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			//if (is_skipcell(oc_id)) {
			if(minmax[oc_id].min > iso_val) {
				continue;
			//} else if (level <= 4 ||(get_err_grad(oc_id) > err_tol && level != oct_depth)) { // heart_0
			//} else if (level < oct_depth-2 || (is_skipcell(oc_id) == 0 && level < oct_depth-2) ||
			//	(get_err(oc_id) > err_tol && r < 60.0f && level < oct_depth) ){
			} else if (level <= start_level || (get_err_grad(oc_id) > err_tol && level < end_level)) {    // head

				cur_queue.Add(oc_id);
				oct_array[oc_id].refine_flag=1;
			} else {
				cut_array[leaf_num++]=oc_id;
			}
		}
		
		while (cur_queue.Get(oc_id) >= 0) {
			level=get_level(oc_id);
			for (i=0; i<8; i++)
				prev_queue.Add(child(oc_id,level,i));
		}
		
	}
}*/

/*void Octree::traverse_qef(float err_tol)
{
	CellQueue prev_queue,cur_queue;
	int i,oc_id,level, start_level, end_level;
	int x, y, z, tx, ty, tz, cell_size;
	float r, r_min, temp, center;

	oc_id=0;
	leaf_num=0;
	center = (dim[0]-1.0f)/2.0f;

	start_level = oct_depth - 3;	// - 3
	end_level = oct_depth;
	if(flag_type == 2) {
		start_level = oct_depth - 3;//3
		end_level = start_level + 1;
	}
	if(flag_type == 3) {
		start_level = oct_depth - 3;//3
		end_level = start_level + 1;
	}

	memset(oct_array,0,sizeof(octcell)*octcell_num);
	prev_queue.Add(oc_id);
	while (!prev_queue.Empty()) {
		
		//leaf_idx=leaf_num;
		while (prev_queue.Get(oc_id) >= 0)  {
			level=get_level(oc_id);
			
			octcell2xyz(oc_id,x,y,z,level);
			cell_size = (dim[0]-1)/(1<<level);
			tx = x*cell_size;		ty = y*cell_size;		tz = z*cell_size;
			r = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			r_min = r;
			tx = (x+1)*cell_size;	ty = y*cell_size;		tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			tx = x*cell_size;		ty = (y+1)*cell_size;	tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = (y+1)*cell_size;	tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			
			tx = x*cell_size;		ty = y*cell_size;		tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = y*cell_size;		tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			tx = x*cell_size;		ty = (y+1)*cell_size;	tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = (y+1)*cell_size;	tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			//if (is_skipcell(oc_id)) {
			if(minmax[oc_id].min > iso_val) {
				continue;
			//} else if (level <= 4 ||(get_err_grad(oc_id) > err_tol && level != oct_depth)) { // heart_0
			} else if (level <= start_level || (get_err_grad(oc_id) > err_tol && level < end_level)) {    // head

			//} else if (level <= 3) {  //4 

			//} else if (level <= oct_depth-4 || (level <= oct_depth-3 && is_skipcell(oc_id) == 0)
			//	||(get_err_grad(oc_id) > err_tol && level != oct_depth)) {

			//} else if (level < oct_depth) {

			//} else if ((level <= 1 && r > (dim[0]-1)/4.0*1.733) || (level <= 2 && r < (dim[0]-1)/4.0*1.733)
			//	|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 && r < (dim[0]-1)/4.0*1.515)
			//	||(get_err_grad(oc_id) > err_tol && level < oct_depth && r < 16*1.8)) { // 14.3

			//} else if (level < oct_depth-2 || (is_skipcell(oc_id) == 0 && level < oct_depth-2) ||
			//	(get_err(oc_id) > err_tol && r < 60.0f && level < oct_depth) ){

				cur_queue.Add(oc_id);
				oct_array[oc_id].refine_flag=1;
			} else {
				cut_array[leaf_num++]=oc_id;
			}
		}
		
		while (cur_queue.Get(oc_id) >= 0) {
			level=get_level(oc_id);
			for (i=0; i<8; i++)
				prev_queue.Add(child(oc_id,level,i));
		}
		
	}
}*/

void Octree::traverse_qef(float err_tol)
{
	CellQueue prev_queue,cur_queue;
	int i,oc_id,level, start_level, end_level;
	int x, y, z, tx, ty, tz, cell_size;
	float r, r_min, temp, center;

	oc_id=0;
	leaf_num=0;
	center = (dim[0]-1.0f)/2.0f;

	start_level = oct_depth - 3;	// - 3
	end_level = oct_depth;
	if(flag_type == 2) {
		start_level = oct_depth - 3;//3
		end_level = start_level + 1;
	}
	if(flag_type == 3) {
		start_level = oct_depth - 3;//3
		end_level = start_level + 1;
	}

	//memset(oct_array,0,sizeof(octcell)*octcell_num);
	std::fill(oct_array.begin(),oct_array.end(),octcell());
	prev_queue.Add(oc_id);
	while (!prev_queue.Empty()) {
		
		//leaf_idx=leaf_num;
		while (prev_queue.Get(oc_id) >= 0)  {
			level=get_level(oc_id);
			
			octcell2xyz(oc_id,x,y,z,level);
			cell_size = (dim[0]-1)/(1<<level);
			tx = x*cell_size;		ty = y*cell_size;		tz = z*cell_size;
			r = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			r_min = r;
			tx = (x+1)*cell_size;	ty = y*cell_size;		tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			tx = x*cell_size;		ty = (y+1)*cell_size;	tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = (y+1)*cell_size;	tz = z*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			
			tx = x*cell_size;		ty = y*cell_size;		tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = y*cell_size;		tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			tx = x*cell_size;		ty = (y+1)*cell_size;	tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;
			tx = (x+1)*cell_size;	ty = (y+1)*cell_size;	tz = (z+1)*cell_size;
			temp = (float) sqrt((tx-center)*(tx-center) + (ty-center)*(ty-center) + (tz-center)*(tz-center));
			if(temp > r) r = temp;
			if(temp < r_min) r_min = temp;

			//if (is_skipcell(oc_id)) {
			if(minmax[oc_id].min > iso_val) {
				continue;
			//} else if (level <= 4 ||(get_err_grad(oc_id) > err_tol && level != oct_depth)) { // heart_0
			//} else if (level <= start_level || (get_err_grad(oc_id) > err_tol && level < end_level)) {    // head
				
				
			} else if (level < oct_depth-2 || (is_skipcell(oc_id) == 0 && level < oct_depth-2) ||
				(get_err(oc_id) > err_tol && r < 60.0f && level < oct_depth) ){

			//} else if (level <= 3) {  //4 

			//} else if (level <= oct_depth-4 || (level <= oct_depth-3 && is_skipcell(oc_id) == 0)
			//	||(get_err_grad(oc_id) > err_tol && level != oct_depth)) {

			//} else if (level < oct_depth) {

			//} else if ((level <= 1 && r > (dim[0]-1)/4.0*1.733) || (level <= 2 && r < (dim[0]-1)/4.0*1.733)
			//	|| (level <= oct_depth-2 && is_skipcell(oc_id) == 0 && r < (dim[0]-1)/4.0*1.515)
			//	||(get_err_grad(oc_id) > err_tol && level < oct_depth && r < 16*1.8)) { // 14.3

				cur_queue.Add(oc_id);
				oct_array[oc_id].refine_flag=1;
			} else {
				cut_array[leaf_num++]=oc_id;
			}
		}
		
		while (cur_queue.Get(oc_id) >= 0) {
			level=get_level(oc_id);
			for (i=0; i<8; i++)
				prev_queue.Add(child(oc_id,level,i));
		}
		
	}
}



void Octree::traverse_qef_interval(float err_tol, float err_tol_in)
{
	CellQueue prev_queue,cur_queue;
	int i,oc_id,level;

	oc_id=0;
	leaf_num=0;

	//memset(oct_array,0,sizeof(octcell)*octcell_num);
	std::fill(oct_array.begin(),oct_array.end(),octcell());
	prev_queue.Add(oc_id);
	while (!prev_queue.Empty()) {
		
		//leaf_idx=leaf_num;
		while (prev_queue.Get(oc_id) >= 0)  {
			level=get_level(oc_id);
			
			//octcell2xyz(oc_id,x,y,z,level);
			
			//if (is_skipcell(oc_id)) {
			if(minmax[oc_id].min > iso_val || minmax[oc_id].max < iso_val_in) {
				continue;
			//} else if (level <= 3 ||(get_err_grad(oc_id) > err_tol && level != oct_depth)) {
			//} else if (level <= 3 ||(is_skipcell(oc_id) == 0  && get_err_grad(oc_id) > err_tol && level != oct_depth)
			//	||(minmax[oc_id].max > iso_val_in &&  minmax[oc_id].min < iso_val_in && get_err_grad(oc_id) > err_tol_in && level != oct_depth)) {
			/*
			} else if (level <= oct_depth-4 ||(is_skipcell_interval(oc_id) == 0  && get_err_grad(oc_id) > err_tol && level != oct_depth)
				|| (level <= oct_depth-3 && is_skipcell_interval(oc_id) == 0)
				|| (minmax[oc_id].max > iso_val_in &&  minmax[oc_id].min < iso_val_in && get_err_grad(oc_id) > err_tol_in && level != oct_depth)) {
			*/
			//} else if (level <= oct_depth-2) {
			// hemoglobin_fat 129
			} else if (level <= 5 ||(is_skipcell_interval(oc_id) == 0  && get_err_grad(oc_id) > err_tol && level != oct_depth)
				//|| (level <=4 && is_skipcell_interval(oc_id) == 0)
				|| (minmax[oc_id].max > iso_val_in &&  minmax[oc_id].min < iso_val_in && get_err_grad(oc_id) > err_tol_in && level != oct_depth)) {

				cur_queue.Add(oc_id);
				oct_array[oc_id].refine_flag=1;
			} else {
				cut_array[leaf_num++]=oc_id;
			}
		}
		
		while (cur_queue.Get(oc_id) >= 0) {
			level=get_level(oc_id);
			for (i=0; i<8; i++)
				prev_queue.Add(child(oc_id,level,i));
		}
		
	}
}

int Octree::is_min_edge(int oc_id, int e_id, unsigned int* vtx, int& vtx_num, int intersect_id, geoframe& geofrm)
{
	int x,y,z,level,i;
	unsigned int temp_vtx[4];
	
	level=get_level(oc_id);
	octcell2xyz(oc_id,x,y,z,level);
	
	vtx_num=4;

	//debug
	//for(i=0;i<4;i++)
	//if(vtx[i]==4294967295)
	//printf("boom! vtx[%d] == 4294967295 before possible change\n",i);

	switch (e_id) {
		case 0 : 
			if (is_refined(x,y,z-1,level) || is_refined(x,y-1,z-1,level) || is_refined(x,y-1,z,level)) return 0;
			temp_vtx[1]=min_vtx(x,y,z-1,level,geofrm);
			temp_vtx[2]=min_vtx(x,y-1,z-1,level,geofrm);
			temp_vtx[3]=min_vtx(x,y-1,z,level,geofrm);
			break;

		case 1 : 
			if (is_refined(x,y-1,z,level) || is_refined(x+1,y-1,z,level) || is_refined(x+1,y,z,level)) return 0;
			temp_vtx[1]=min_vtx(x+1,y,z,level,geofrm);
			temp_vtx[2]=min_vtx(x+1,y-1,z,level,geofrm);
			temp_vtx[3]=min_vtx(x,y-1,z,level,geofrm);
			break;

		case 2 : 
			if (is_refined(x,y,z+1,level) || is_refined(x,y-1,z+1,level) || is_refined(x,y-1,z,level)) return 0;
			temp_vtx[1]=min_vtx(x,y,z+1,level,geofrm);
			temp_vtx[2]=min_vtx(x,y-1,z+1,level,geofrm);
			temp_vtx[3]=min_vtx(x,y-1,z,level,geofrm);
			break;

		case 3 : 
			if (is_refined(x,y-1,z,level) || is_refined(x-1,y-1,z,level) || is_refined(x-1,y,z,level)) return 0;
			temp_vtx[1]=min_vtx(x,y-1,z,level,geofrm);
			temp_vtx[2]=min_vtx(x-1,y-1,z,level,geofrm);
			temp_vtx[3]=min_vtx(x-1,y,z,level,geofrm);
			break;

		case 4 : 
			if (is_refined(x,y,z-1,level) || is_refined(x,y+1,z-1,level) || is_refined(x,y+1,z,level)) return 0;
			temp_vtx[1]=min_vtx(x,y+1,z,level,geofrm);
			temp_vtx[2]=min_vtx(x,y+1,z-1,level,geofrm);
			temp_vtx[3]=min_vtx(x,y,z-1,level,geofrm);
			break;

		case 5 :
			if (is_refined(x,y+1,z,level) || is_refined(x+1,y,z,level) || is_refined(x+1,y+1,z,level)) return 0;
			temp_vtx[1]=min_vtx(x,y+1,z,level,geofrm);
			temp_vtx[2]=min_vtx(x+1,y+1,z,level,geofrm);
			temp_vtx[3]=min_vtx(x+1,y,z,level,geofrm);
			break;

		case 6 :
			if (is_refined(x,y+1,z,level) || is_refined(x,y+1,z+1,level) || is_refined(x,y,z+1,level)) return 0;
			temp_vtx[1]=min_vtx(x,y+1,z,level,geofrm);
			temp_vtx[2]=min_vtx(x,y+1,z+1,level,geofrm);
			temp_vtx[3]=min_vtx(x,y,z+1,level,geofrm);
			break;

		case 7 :
			if (is_refined(x-1,y,z,level) || is_refined(x-1,y+1,z,level) || is_refined(x,y+1,z,level)) return 0;
			temp_vtx[1]=min_vtx(x-1,y,z,level,geofrm);
			temp_vtx[2]=min_vtx(x-1,y+1,z,level,geofrm);
			temp_vtx[3]=min_vtx(x,y+1,z,level,geofrm);
			break;

		case 8 :
			if (is_refined(x,y,z-1,level) || is_refined(x-1,y,z-1,level) || is_refined(x-1,y,z,level)) return 0;
			temp_vtx[1]=min_vtx(x-1,y,z,level,geofrm);
			temp_vtx[2]=min_vtx(x-1,y,z-1,level,geofrm);
			temp_vtx[3]=min_vtx(x,y,z-1,level,geofrm);
			break;

		case 9 :
			if (is_refined(x,y,z-1,level) || is_refined(x+1,y,z-1,level) || is_refined(x+1,y,z,level)) return 0;
			temp_vtx[1]=min_vtx(x,y,z-1,level,geofrm);
			temp_vtx[2]=min_vtx(x+1,y,z-1,level,geofrm);
			temp_vtx[3]=min_vtx(x+1,y,z,level,geofrm);
			break;

		case 10 :
			if (is_refined(x,y,z+1,level) || is_refined(x-1,y,z+1,level) || is_refined(x-1,y,z,level)) return 0;
			temp_vtx[1]=min_vtx(x,y,z+1,level,geofrm);
			temp_vtx[2]=min_vtx(x-1,y,z+1,level,geofrm);
			temp_vtx[3]=min_vtx(x-1,y,z,level,geofrm);
			break;

		case 11 :
			if (is_refined(x,y,z+1,level) || is_refined(x+1,y,z+1,level) || is_refined(x+1,y,z,level)) return 0;
			temp_vtx[1]=min_vtx(x+1,y,z,level,geofrm);
			temp_vtx[2]=min_vtx(x+1,y,z+1,level,geofrm);
			temp_vtx[3]=min_vtx(x,y,z+1,level,geofrm);
			break;
	}

	temp_vtx[0] = min_vtx(x,y,z,level,geofrm);

	assert(intersect_id==1 || intersect_id==-1 || intersect_id==3 || intersect_id==-3);

	if (intersect_id==1 || intersect_id==3) 
		for (i=0;i<4;i++) vtx[i]=temp_vtx[i];
	else if (intersect_id==-1 || intersect_id==-3)
		for (i=0;i<4;i++) vtx[i]=temp_vtx[3-i];

	//debug
	//	for(i=0;i<4;i++)
	//	  if(vtx[i]==4294967295)
	//	    printf("boom! vtx[%d] == 4294967295 after change\n",i);

	//min_vtx() can return -1, so lets check for that
	for(i=0;i<4;i++)
	  if(vtx[i]==(unsigned int)(-1))
	    return 0;

	return 1;

}

void Octree::get_solution(int oc_id, float* pos)
{
	int x,y,z;
	float val[8];
	int level = get_level(oc_id);
	
	int cell_size = (dim[0]-1)/(1<<level);
		
	getCellValues(oc_id, level, val);
	octcell2xyz(oc_id, x, y, z, level);

//	assert(minmax[oc_id].max >= iso_val);
	assert(minmax[oc_id].max >= iso_val || minmax[oc_id].min <= iso_val_in);

	//if(in_out == 0) {
	if(is_skipcell(oc_id) == 0) {
	  pos[0]=qef_array[oc_id] ? qef_array[oc_id][9] : 0.0;
	  pos[1]=qef_array[oc_id] ? qef_array[oc_id][10] : 0.0;
	  pos[2]=qef_array[oc_id] ? qef_array[oc_id][11] : 0.0;
	}
	else {
	  pos[0]=qef_array[oc_id] ? qef_array_in[oc_id][9] : 0.0;
	  pos[1]=qef_array[oc_id] ? qef_array_in[oc_id][10] : 0.0;
	  pos[2]=qef_array[oc_id] ? qef_array_in[oc_id][11] : 0.0;
	}

	// to be corrected
	if (!(pos[0] > x*cell_size && pos[0] < (x+1)*cell_size)) {
		pos[0] = x*cell_size + cell_size/2.;
	}
		
	if (!(pos[1] > y*cell_size && pos[1] < (y+1)*cell_size)) {
		pos[1] = y*cell_size + cell_size/2.;
	}
		
	if (!(pos[2] > z*cell_size && pos[2] < (z+1)*cell_size)) {
		pos[2] = z*cell_size + cell_size/2.;
	}
}

void Octree::get_vtx(int x, int y, int z, int level, float* pos)
{
	int oc_id;
	oc_id = xyz2octcell(x,y,z,level);
	get_solution(oc_id,pos);
}

int Octree::min_vtx(int x, int y, int z, int level, geoframe& geofrm)
{
	int tx, ty, tz, vert;
	float vtx[3], norm[3];
	tx = x;	ty = y;	tz = z;

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

//	if(minmax[oc_id].max <= iso_val)
	if(minmax[oc_id].max <= iso_val && minmax[oc_id].min >= iso_val_in)
		return -1;
	else {
		get_vtx(tx, ty, tz, level, vtx);
		getVertGrad(tx*cell_size, ty*cell_size, tz*cell_size, norm);
		if(in_out == 0) {
			if ((vert = vtx_idx_arr[xyz2octcell(tx, ty, tz, level)]) == -1) {
				vert = geofrm.AddVert(vtx, norm);
				geofrm.AddBound(vert, 1);
				vtx_idx_arr[xyz2octcell(tx, ty, tz, level)] = vert;
				return vert;
			}
			else
				return vert;
		}
		else {
			if ((vert = vtx_idx_arr_in[xyz2octcell(tx, ty, tz, level)]) == -1) {
				vert = geofrm.AddVert(vtx, norm);
				geofrm.AddBound(vert, -1);
				vtx_idx_arr_in[xyz2octcell(tx, ty, tz, level)] = vert;
				return vert;
			}
			else
				return vert;
		}
	}
}

void Octree::eflag_clear()
{
  //memset(ebit,0,octcell_num*4/8);
  std::fill(ebit.begin(),ebit.end(),0);
}

int Octree::is_eflag_on(int x, int y, int z, int level, int e)
{
	int idx;
	switch (e) {
	case 0 : 
		idx = 3*xyz2octcell(x,y,z,level) + 0;
		break;
	case 1 :
		idx = 3*xyz2octcell(x+1,y,z,level) + 2;
		break;
	case 2 :
		idx = 3*xyz2octcell(x,y,z+1,level) + 0;
		break;
	case 3 :
		idx = 3*xyz2octcell(x,y,z,level) + 2;
		break;
	case 4 :
		idx = 3*xyz2octcell(x,y+1,z,level) + 0;
		break;
	case 5 :
		idx = 3*xyz2octcell(x+1,y+1,z,level) + 2;
		break;
	case 6 :
		idx = 3*xyz2octcell(x,y+1,z+1,level) + 0;
		break;
	case 7 :
		idx = 3*xyz2octcell(x,y+1,z,level) + 2;
		break;
	case 8 :
		idx = 3*xyz2octcell(x,y,z,level) + 1;
		break;
	case 9 :
		idx = 3*xyz2octcell(x+1,y,z,level) + 1;
		break;
	case 10 :
		idx = 3*xyz2octcell(x,y,z+1,level) + 1;
		break;
	case 11 :
		idx = 3*xyz2octcell(x+1,y,z+1,level) + 1;
		break;
	}
	
	if (ebit[idx/8]&(1<<(idx%8))) return 1;
	else return 0;
	
}

void Octree::eflag_on(int x, int y, int z, int level, int e)
{
	int idx;
	switch (e) {
	case 0 : 
		idx = 3*xyz2octcell(x,y,z,level) + 0;
		break;
	case 1 :
		idx = 3*xyz2octcell(x+1,y,z,level) + 2;
		break;
	case 2 :
		idx = 3*xyz2octcell(x,y,z+1,level) + 0;
		break;
	case 3 :
		idx = 3*xyz2octcell(x,y,z,level) + 2;
		break;
	case 4 :
		idx = 3*xyz2octcell(x,y+1,z,level) + 0;
		break;
	case 5 :
		idx = 3*xyz2octcell(x+1,y+1,z,level) + 2;
		break;
	case 6 :
		idx = 3*xyz2octcell(x,y+1,z+1,level) + 0;
		break;
	case 7 :
		idx = 3*xyz2octcell(x,y+1,z,level) + 2;
		break;
	case 8 :
		idx = 3*xyz2octcell(x,y,z,level) + 1;
		break;
	case 9 :
		idx = 3*xyz2octcell(x+1,y,z,level) + 1;
		break;
	case 10 :
		idx = 3*xyz2octcell(x,y,z+1,level) + 1;
		break;
	case 11 :
		idx = 3*xyz2octcell(x+1,y,z+1,level) + 1;
		break;
	}
	
	ebit[idx/8]|=(1<<(idx%8));
	
}

void Octree::tetra_to_4_hexa(geoframe& geofrm) {

	// get a list of leaf cells satisfying the error tolerance.
	// for each cell, get the centroid vertex

	int x, y, z, valid_leaf, cell_size, level;
	int vtx_num, intersect_id, i, j;
	unsigned int vtx[4];
	float val[8];

	int k;
	for (k=0;k<octcell_num;k++) vtx_idx_arr[k] = -1;

	for (i = 0; i < leaf_num; i++ ) {
		valid_leaf = cut_array[i] ;
		level = get_level(valid_leaf) ;
		cell_size = (dim[0]-1)/(1<<level);
		octcell2xyz(valid_leaf, x, y, z, level);
		getCellValues(valid_leaf, level, val);
		for (j = 0 ; j < 12 ; j++ ) {
			if (is_eflag_on(x,y,z,level,j)) continue;
			intersect_id = is_intersect(val , j);

			if (intersect_id == 1 || intersect_id == -1) {
				if (is_min_edge(valid_leaf , j , vtx , vtx_num , intersect_id , geofrm)) {
					eflag_on(x , y , z , level , j);
					geofrm.AddQuad_hexa(vtx , vtx_num);
				}
			}
		}
	}

}

void Octree::polygonize(geoframe& geofrm) {

	// get a list of leaf cells satisfying the error tolerance.
	// for each cell, get the centroid vertex

	int x, y, z, valid_leaf, cell_size, level;
	int vtx_num, intersect_id, i, j;
	unsigned int vtx[4];
	float val[8];

	in_out = 0;
	int k;
	for (k = 0; k < octcell_num; k++) vtx_idx_arr[k] = -1;

	printf("leafnum:%d\n",leaf_num);
	for (i = 0; i < leaf_num; i++ ) {
		valid_leaf = cut_array[i] ;
		level = get_level(valid_leaf) ;
		cell_size = (dim[0]-1)/(1<<level);
		octcell2xyz(valid_leaf, x, y, z, level);
		getCellValues(valid_leaf, level, val);

		for (j = 0 ; j < 12 ; j++ ) {
			if (is_eflag_on(x,y,z,level,j)) continue;
			intersect_id = is_intersect(val , j);

			if (intersect_id == 1 || intersect_id == -1) {
				if (is_min_edge(valid_leaf , j , vtx , vtx_num , intersect_id , geofrm)) {
					eflag_on(x , y , z , level , j);
					geofrm.Add_2_Tri(vtx);
				}
			}
		}
	}

}

void Octree::polygonize_interval(geoframe& geofrm) {

	// get a list of leaf cells satisfying the error tolerance.
	// for each cell, get the centroid vertex

	int x, y, z, valid_leaf, cell_size, level;
	int vtx_num, intersect_id, i, j;
	unsigned int vtx[4];
	float val[8];

	int k;
	for (k = 0; k < octcell_num; k++) {vtx_idx_arr[k] = -1;	vtx_idx_arr_in[k] = -1;}

	for (i = 0; i < leaf_num; i++ ) {
		valid_leaf = cut_array[i] ;
		level = get_level(valid_leaf) ;
		cell_size = (dim[0]-1)/(1<<level);
		octcell2xyz(valid_leaf, x, y, z, level);
		getCellValues(valid_leaf, level, val);

		//if(y*cell_size < 4 || y*cell_size >= 60 || x*cell_size < 12 || x*cell_size >= 52)
		//	continue;

		for (j = 0; j < 12; j++ ) {
			if (is_eflag_on(x, y, z, level, j)) continue;
			intersect_id = is_intersect_interval(val, j);

			if (intersect_id == 1 || intersect_id == -1) {
				if(is_skipcell(valid_leaf) == 0) in_out = 0;
				else in_out = 1;
				if (is_min_edge(valid_leaf, j, vtx, vtx_num, intersect_id, geofrm)) {
					eflag_on(x, y, z, level, j);
					//geofrm.AddQuad(vtx, vtx_num);
					geofrm.Add_2_Tri(vtx);
				}
			}

			if (intersect_id == 3 || intersect_id == -3) {
				in_out = 1;
				if (is_min_edge(valid_leaf, j, vtx, vtx_num, intersect_id, geofrm)) {
					eflag_on(x, y, z, level, j);
					geofrm.AddQuad(vtx, vtx_num);
					in_out = 0;
					is_min_edge(valid_leaf, j, vtx, vtx_num, intersect_id, geofrm);
					//geofrm.AddQuad(vtx, vtx_num);
					geofrm.Add_2_Tri(vtx);
				}
			}
		}
	}

}
/*
// read *.rawiv file to show color, e.g., the potential distribution
void Octree::func_val(geoframe& geofrm) {
	int i, j, x, y, z, oc_id, vtx[8];
	float dx, dy, dz, val[8], func_min, func_max;
	float* func_vol;
	FILE *func_fp, *out_func;

	func_fp = fopen("rawiv/result_solid257.rawiv","rb"); // bladder
	//func_fp = fopen("rawiv/p8-segment-100-137.rawiv","rb");
	
	if (func_fp==NULL) {
		printf("wrong name : %s\n","p8-segment-100-137.rawiv");
		return;
	}
	
	getFloat(minext,3,func_fp);	getFloat(maxext,3,func_fp);
	getInt(&nverts,1,func_fp);	getInt(&ncells,1,func_fp);
	getInt(dim,3,func_fp);		getFloat(orig,3,func_fp);		getFloat(span,3,func_fp);

	func_vol = (float*)malloc(sizeof(float)*dim[0]*dim[1]*dim[2]);
	getFloat(func_vol, dim[0]*dim[1]*dim[2], func_fp);

	fclose(func_fp);
	out_func = fopen("data/out_func", "w");

	func_min = 100.0f;
	func_max = -100.0f;
	int sum_15 = 0;	int sum_10 = 0;	int sum_5 = 0;	int sum_1 = 0;
	for(i = 0; i < geofrm.numverts; i++) {
		x = (int) geofrm.verts[i][0];
		y = (int) geofrm.verts[i][1];
		z = (int) geofrm.verts[i][2];

		dx = geofrm.verts[i][0] - x;
		dy = geofrm.verts[i][1] - y;
		dz = geofrm.verts[i][2] - z;

		oc_id = xyz2octcell(x, y, z, oct_depth);
		idx2vtx(oc_id, oct_depth, vtx);
		for (j = 0; j < 8; j++)	val[j] = func_vol[vtx[j]];

		geofrm.funcs[i][0] = (1-dx)*(1-dy)*(1-dz)*val[0] + (1-dx)*(1-dy)*dz*val[3]
						+ (1-dx)*dy*(1-dz)*val[4] + dx*(1-dy)*(1-dz)*val[1] + (1-dx)*dy*dz*val[7]
						+ dx*(1-dy)*dz*val[2] + dx*dy*(1-dz)*val[5] + dx*dy*dz*val[6];

		// RDV or Bladder
		geofrm.funcs[i][0] = val[0];
		for (j = 1; j < 8; j++) {
			if(geofrm.funcs[i][0] < val[j])	geofrm.funcs[i][0] = val[j];
		}

		if(geofrm.funcs[i][0] < func_min)	func_min = geofrm.funcs[i][0];
		if(geofrm.funcs[i][0] > func_max)	func_max = geofrm.funcs[i][0];
		if(fabs(geofrm.funcs[i][0]) > 15.0) sum_15++;
		if(fabs(geofrm.funcs[i][0]) > 10.0) sum_10++;
		if(fabs(geofrm.funcs[i][0]) > 5.0) sum_5++;
		if(fabs(geofrm.funcs[i][0]) > 1.0) sum_1++;
		fprintf(out_func, "%f %f %f %f %f %f\n", (x-32)*1.578125+1.2735, (y-32)*1.578125-4.0895,
											  (z-32)*1.578125+4.0065, geofrm.funcs[i][0], func_min, func_max);
	}
	fprintf(out_func, "sum_pot 15, 10, 5, 1 = %d %d %d %d\n", sum_15, sum_10, sum_5, sum_1);
	for(i = 0; i < geofrm.numverts; i++) {
		//geofrm.funcs[i][0] = (geofrm.funcs[i][0] - func_min) / (func_max - func_min)*2.0f - 1.0f;
	}

	fclose(out_func);

	free(func_vol);
}
*/
// read *.rawiv file to show color, e.g., the potential distribution
void Octree::func_val(geoframe& geofrm) {

	if(!prop_flag)
		return;
	prop_flag = 0;
	int i, j, x, y, z, oc_id, vtx[8];
	float dx, dy, dz, val[8], func_min, func_max;
	//	float* func_vol;
	/*FILE *func_fp, *out_func;

	//func_fp = fopen("rawiv/test/potential211257.rawiv","rb");
	//func_fp = fopen("rawiv/1C2B_full_pot33.rawiv","rb");
	//func_fp = fopen("C:/Documents and Settings/jessica/Desktop/NMJ_rawiv/1C2B_pot129.rawiv","rb");
	//func_fp = fopen("C:/Documents and Settings/jessica/Desktop/NMJ_rawiv/2BG9_pot97129.rawiv","rb");
	func_fp = fopen(prop_fname,"rb");
	
	if (func_fp==NULL) {
		printf("wrong name : %s\n","potential211257.rawiv");
		return;
	}

	strcpy(prop_fname,"");
	
	getFloat(minext,3,func_fp);	getFloat(maxext,3,func_fp);
	getInt(&nverts,1,func_fp);	getInt(&ncells,1,func_fp);
	getInt(dim,3,func_fp);		getFloat(orig,3,func_fp);		getFloat(span,3,func_fp);

	func_vol = (float*)malloc(sizeof(float)*dim[0]*dim[1]*dim[2]);
	getFloat(func_vol, dim[0]*dim[1]*dim[2], func_fp);

	fclose(func_fp);*/

	//func_vol = (float*)malloc(sizeof(float)*dim[0]*dim[1]*dim[2]);

	VolMagick::Volume propData;
	VolMagick::readVolumeFile(propData,prop_fname);

	int old_flag = interior_flag;
	interior_flag = 1;
	loadData(propData/*,func_vol*/);
	interior_flag = old_flag;
	//out_func = fopen("rawiv/out_func", "w");

	func_min = 100.0f;
	func_max = -100.0f;
	int sum_15 = 0;	int sum_10 = 0;	int sum_5 = 0;	int sum_1 = 0;
	//fprintf(out_func, "%d\n", geofrm.numverts);
	for(i = 0; i < geofrm.numverts; i++) {
		x = (int) geofrm.verts[i][0];
		y = (int) geofrm.verts[i][1];
		z = (int) geofrm.verts[i][2];

		dx = geofrm.verts[i][0] - x;
		dy = geofrm.verts[i][1] - y;
		dz = geofrm.verts[i][2] - z;

		oc_id = xyz2octcell(x, y, z, oct_depth);
		idx2vtx(oc_id, oct_depth, vtx);
		for (j = 0; j < 8; j++)	val[j] = orig_vol[vtx[j]];//func_vol[vtx[j]];

		geofrm.funcs[i][0] = (1-dx)*(1-dy)*(1-dz)*val[0] + (1-dx)*(1-dy)*dz*val[3]
						+ (1-dx)*dy*(1-dz)*val[4] + dx*(1-dy)*(1-dz)*val[1] + (1-dx)*dy*dz*val[7]
						+ dx*(1-dy)*dz*val[2] + dx*dy*(1-dz)*val[5] + dx*dy*dz*val[6];
		if(geofrm.funcs[i][0] < func_min)	func_min = geofrm.funcs[i][0];
		if(geofrm.funcs[i][0] > func_max)	func_max = geofrm.funcs[i][0];
		if(fabs(geofrm.funcs[i][0]) > 15.0) sum_15++;
		if(fabs(geofrm.funcs[i][0]) > 10.0) sum_10++;
		if(fabs(geofrm.funcs[i][0]) > 5.0) sum_5++;
		if(fabs(geofrm.funcs[i][0]) > 1.0) sum_1++;
		//fprintf(out_func, "%f %f %f %f %f %f\n", geofrm.verts[i][0], geofrm.verts[i][1], geofrm.verts[i][2], geofrm.funcs[i][0], func_min, func_max);
	}
	//fprintf(out_func, "sum_pot 15, 10, 5, 1 = %d %d %d %d\n", sum_15, sum_10, sum_5, sum_1);
	for(i = 0; i < geofrm.numverts; i++) {
		//geofrm.funcs[i][0] = (geofrm.funcs[i][0] - func_min) / (func_max - func_min)*2.0f - 1.0f;
	}

	//fclose(out_func);

	//free(func_vol);
}
/*
// read *.rawv file to show the color, e.g., hemoglobin
void Octree::func_val(geoframe& geofrm) {
	int i, j, x, y, z, oc_id, vtx[8];
	float dx, dy, dz, val[8];
	unsigned char* func_r;
	unsigned char* func_g;
	unsigned char* func_b;
	unsigned char* temp;
	FILE *func_fp;

	//func_fp = fopen("rawiv/1FSS.rawv", "rb");
	//func_fp = fopen("J:/AChR/nachr_ha7_hom3_mod10_0.1.rawv", "rb");
	func_fp = fopen("J:/AChR/1I9B_ptn0.1_sphere.rawv", "rb");
	//func_fp = fopen("rawiv/bio_mesh/1JJ2/1JJ2129_00625.rawv", "rb");
	
	if (func_fp==NULL) {
		printf("wrong name : %s\n", "heme129.rawv");
		return;
	}
	
	// skip the header 56+(1+64)*4 = 316 bytes
	temp = (unsigned char*)malloc(sizeof(unsigned char)*316);
	getUnChar(temp, 316, func_fp);

	// read r, g, b
	func_r = (unsigned char*)malloc(sizeof(unsigned char)*dim[0]*dim[1]*dim[2]);
	getUnChar(func_r, dim[0]*dim[1]*dim[2], func_fp);
	func_g = (unsigned char*)malloc(sizeof(unsigned char)*dim[0]*dim[1]*dim[2]);
	getUnChar(func_g, dim[0]*dim[1]*dim[2], func_fp);
	func_b = (unsigned char*)malloc(sizeof(unsigned char)*dim[0]*dim[1]*dim[2]);
	getUnChar(func_b, dim[0]*dim[1]*dim[2], func_fp);

	fclose(func_fp);

	//geofrm.funcs = (float(*)[3])malloc(sizeof(float(*)[3])*10000000);

	for(i = 0; i < geofrm.numverts; i++) {

		if(geofrm.bound_sign[i] != 1) continue;

		x = (int) geofrm.verts[i][0];
		y = (int) geofrm.verts[i][1];
		z = (int) geofrm.verts[i][2];

		dx = geofrm.verts[i][0] - x;
		dy = geofrm.verts[i][1] - y;
		dz = geofrm.verts[i][2] - z;

		oc_id = xyz2octcell(x, y, z, oct_depth);
		idx2vtx(oc_id, oct_depth, vtx);

		for (j = 0; j < 8; j++)	val[j] = (float) func_r[vtx[j]];
		geofrm.funcs[i][0] = (1-dx)*(1-dy)*(1-dz)*val[0] + (1-dx)*(1-dy)*dz*val[3]
						+ (1-dx)*dy*(1-dz)*val[4] + dx*(1-dy)*(1-dz)*val[1] + (1-dx)*dy*dz*val[7]
						+ dx*(1-dy)*dz*val[2] + dx*dy*(1-dz)*val[5] + dx*dy*dz*val[6];

		for (j = 0; j < 8; j++)	val[j] = (float) func_g[vtx[j]];
		geofrm.funcs[i][1] = (1-dx)*(1-dy)*(1-dz)*val[0] + (1-dx)*(1-dy)*dz*val[3]
						+ (1-dx)*dy*(1-dz)*val[4] + dx*(1-dy)*(1-dz)*val[1] + (1-dx)*dy*dz*val[7]
						+ dx*(1-dy)*dz*val[2] + dx*dy*(1-dz)*val[5] + dx*dy*dz*val[6];

		for (j = 0; j < 8; j++)	val[j] = (float) func_b[vtx[j]];
		geofrm.funcs[i][2] = (1-dx)*(1-dy)*(1-dz)*val[0] + (1-dx)*(1-dy)*dz*val[3]
						+ (1-dx)*dy*(1-dz)*val[4] + dx*(1-dy)*(1-dz)*val[1] + (1-dx)*dy*dz*val[7]
						+ dx*(1-dy)*dz*val[2] + dx*dy*(1-dz)*val[5] + dx*dy*dz*val[6];
	}

	free(func_r);
	free(func_g);
	free(func_b);
	free(temp);
}
*/
void Octree::mesh_extract(geoframe& geofrm, float err_tol)
{
	geofrm.Clear();
	eflag_clear();
	vflag_clear();

	//flag_type = 1;
	in_out = 0;
	flag_extend = 1;
	int func = 1;

	//BSplineCoeff = (float*)malloc(sizeof(float)*dim[0]*dim[1]*dim[2]);
	BSplineCoeff.resize(dim[0]*dim[1]*dim[2]);
	LBIE::TransImg2Spline(&(orig_vol[0]), &(BSplineCoeff[0]), dim[0], dim[1], dim[2]);


	switch (flag_type) {
	case 0:							// isosurface
		printf("triangle\n");
		polygonize(geofrm);
		if(func)	func_val(geofrm);
		break;

	case 1:							// tetra mesh
		tetrahedralize(geofrm);
		if(func)	func_val(geofrm);
		break;

	case 2:							// quad mesh
		polygonize_quad(geofrm, err_tol);
		if(func)	func_val(geofrm);
		break;

	case 3:							//hexa_mesh
		hexahedralize(geofrm, err_tol);
		if(func)	func_val(geofrm);
		break;

	case 4:							// interval
		polygonize_interval(geofrm);
		if(func)	func_val(geofrm);
		break;

	case 5:							// interval volume
		tetrahedralize_interval(geofrm);
		if(func)	func_val(geofrm);
		break;
	}
}

void Octree::quality_improve(geoframe& geofrm, int improve_method)
{
	switch (improve_method) {
	        case 0:
			break;
		case 1:
			geometric_flow(geofrm);
			break;
		case 2:
			edge_contraction(geofrm);
			break;
		case 3:
			smoothing_joeliu_volume(geofrm,0);
			break;
		case 4:
			smoothing_joeliu_volume(geofrm,1);
			break;
		case 5:
			optimization(geofrm);
			break;
		default:
			break;
	}
}

int Octree::is_intersect(float* val, int e_id)
{
	float f1,f2;

	f1 = val[cube_eid[e_id][0]];
	f2 = val[cube_eid[e_id][1]];

	if (iso_val <= f1 && iso_val >= f2)
		return -1;
	else if (iso_val <= f2 && iso_val >= f1)
		return 1;
	else if (iso_val >= f1 && f1 >= f2)
		return -2;
	else if (iso_val >= f2 && f2 >= f1)
		return 2;
	else return 0;

}

int Octree::is_intersect_interval(float* val, int e_id)
{
	float f1,f2;

	f1 = val[cube_eid[e_id][0]];
	f2 = val[cube_eid[e_id][1]];

	if ((iso_val <= f1 && iso_val >= f2 && f2 >= iso_val_in) ||
		(iso_val_in <= f1 && iso_val_in >= f2 && f1 <= iso_val))
		return -1;
	else if ((iso_val <= f2 && iso_val >= f1 && f1 >= iso_val_in) ||
			 (iso_val_in <= f2 && iso_val_in >= f1 && f2 <= iso_val))
		return 1;
	else if (iso_val >= f1 && f1 >= f2 && f2 >= iso_val_in)
		return -2;
	else if (iso_val >= f2 && f2 >= f1 && f1 >= iso_val_in)
		return 2;
	else if (f1 >= iso_val && iso_val_in >= f2)
		return -3;
	else if (f2 >= iso_val && iso_val_in >= f1)
		return 3;
	else return 0;

}

/*void Octree::read_header()
{
	getFloat(minext,3,vol_fp);
	getFloat(maxext,3,vol_fp);
	
	getInt(&nverts,1,vol_fp);
	getInt(&ncells,1,vol_fp);
	
	getInt(dim,3,vol_fp);
	getFloat(orig,3,vol_fp);
	getFloat(span,3,vol_fp);

	orig[0] = 0.0f;	orig[1] = 0.0f;	orig[2] = 0.0f;
	span[0] = 1.0f;	span[1] = 1.0f;	span[2] = 1.0f;
	
}*/

void Octree::check_topo() {

	int i, j, k, index[8], t;
	float val[8], my_iso;

	my_iso = 0.5f;
	for(k = 0; k < dim[2]-1; k++)
		for(j = 0; j < dim[1]-1; j++)
			for(i = 0; i < dim[0]-1; i++) {
				index[0] = k*dim[0]*dim[1] + j*dim[0] + i;
				index[1] = k*dim[0]*dim[1] + j*dim[0] + (i+1);
				index[2] = (k+1)*dim[0]*dim[1] + j*dim[0] + (i+1);
				index[3] = (k+1)*dim[0]*dim[1] + j*dim[0] + i;
				index[4] = k*dim[0]*dim[1] + (j+1)*dim[0] + i;
				index[5] = k*dim[0]*dim[1] + (j+1)*dim[0] + (i+1);
				index[6] = (k+1)*dim[0]*dim[1] + (j+1)*dim[0] + (i+1);
				index[7] = (k+1)*dim[0]*dim[1] + (j+1)*dim[0] + i;

				for(t = 0; t < 8; t++) val[t] = orig_vol[index[t]];

				if( val[0] > my_iso && val[6] > my_iso && val[1] < my_iso && val[2] < my_iso &&
					val[3] < my_iso && val[4] < my_iso && val[5] < my_iso && val[7] < my_iso ) {
					orig_vol[index[0]] = my_iso - 0.0002f;	orig_vol[index[6]] = my_iso - 0.0002f;
				}
				else if(val[1] > my_iso && val[7] > my_iso && val[0] < my_iso && val[2] < my_iso &&
						val[3] < my_iso && val[4] < my_iso && val[5] < my_iso && val[6] < my_iso ) {
					orig_vol[index[1]] = my_iso - 0.0002f;	orig_vol[index[7]] = my_iso - 0.0002f;
				}
				else if(val[2] > my_iso && val[4] > my_iso && val[0] < my_iso && val[1] < my_iso &&
						val[3] < my_iso && val[5] < my_iso && val[6] < my_iso && val[7] < my_iso ) {
					orig_vol[index[2]] = my_iso - 0.0002f;	orig_vol[index[4]] = my_iso - 0.0002f;
				}
				else if(val[3] > my_iso && val[5] > my_iso && val[0] < my_iso && val[1] < my_iso &&
						val[2] < my_iso && val[4] < my_iso && val[6] < my_iso && val[7] < my_iso ) {
					orig_vol[index[3]] = my_iso - 0.0002f;	orig_vol[index[5]] = my_iso - 0.0002f;
				}
				if(val[0] > my_iso && val[2] > my_iso && val[1] < my_iso && val[3] < my_iso) {
					orig_vol[index[0]] = my_iso - 0.0002f;	orig_vol[index[2]] = my_iso - 0.0002f;
				}
				else if(val[1] > my_iso && val[3] > my_iso && val[0] < my_iso && val[2] < my_iso) {
					orig_vol[index[1]] = my_iso - 0.0002f;	orig_vol[index[3]] = my_iso - 0.0002f;
				}
				if(val[4] > my_iso && val[6] > my_iso && val[5] < my_iso && val[7] < my_iso) {
					orig_vol[index[4]] = my_iso - 0.0002f;	orig_vol[index[6]] = my_iso - 0.0002f;
				}
				else if(val[5] > my_iso && val[7] > my_iso && val[4] < my_iso && val[6] < my_iso) {
					orig_vol[index[5]] = my_iso - 0.0002f;	orig_vol[index[7]] = my_iso - 0.0002f;
				}
				if(val[0] > my_iso && val[7] > my_iso && val[3] < my_iso && val[4] < my_iso) {
					orig_vol[index[0]] = my_iso - 0.0002f;	orig_vol[index[7]] = my_iso - 0.0002f;
				}
				else if(val[3] > my_iso && val[4] > my_iso && val[0] < my_iso && val[7] < my_iso) {
					orig_vol[index[3]] = my_iso - 0.0002f;	orig_vol[index[4]] = my_iso - 0.0002f;
				}
				if(val[1] > my_iso && val[6] > my_iso && val[2] < my_iso && val[5] < my_iso) {
					orig_vol[index[1]] = my_iso - 0.0002f;	orig_vol[index[6]] = my_iso - 0.0002f;
				}
				else if(val[2] > my_iso && val[5] > my_iso && val[1] < my_iso && val[6] < my_iso) {
					orig_vol[index[2]] = my_iso - 0.0002f;	orig_vol[index[5]] = my_iso - 0.0002f;
				}
				if(val[0] > my_iso && val[5] > my_iso && val[1] < my_iso && val[4] < my_iso) {
					orig_vol[index[0]] = my_iso - 0.0002f;	orig_vol[index[5]] = my_iso - 0.0002f;
				}
				else if(val[1] > my_iso && val[4] > my_iso && val[0] < my_iso && val[5] < my_iso) {
					orig_vol[index[1]] = my_iso - 0.0002f;	orig_vol[index[4]] = my_iso - 0.0002f;
				}
				if(val[3] > my_iso && val[6] > my_iso && val[2] < my_iso && val[7] < my_iso) {
					orig_vol[index[3]] = my_iso - 0.0002f;	orig_vol[index[6]] = my_iso - 0.0002f;
				}
				else if(val[2] > my_iso && val[7] > my_iso && val[3] < my_iso && val[6] < my_iso) {
					orig_vol[index[2]] = my_iso - 0.0002f;	orig_vol[index[7]] = my_iso - 0.0002f;
				}
			}



}


/*void Octree::read_data()
{
	int i;
	int j, k, index;
	//float r, r0, r1;
	float center = (dim[0] - 1.0f)/2.0f;

	// currently support only float data
	getFloat(orig_vol, dim[0]*dim[1]*dim[2], vol_fp);

	//check_topo();

	// interior mesh
	for(i = 0; i < dim[0]*dim[1]*dim[2]; i++) {
		orig_vol[i] = -orig_vol[i];
	}
/*
	// exterior mesh -- 2BG9 add one sphere
	float r0, r;
	r0 = center - 0.001f;
	for(k = 0; k < dim[2]; k++)
		for(j = 0; j < dim[1]; j++)
			for(i = 0; i < dim[0]; i++) {
				index = k*dim[0]*dim[1] + j*dim[0] + i;
				
				orig_vol[index] = orig_vol[index] - 1.0f;
				r = (float) sqrt((i-center)*(i-center) + (j-center)*(j-center) + (k-center)*(k-center));
				if(orig_vol[index] < 0.0) {
					if((r>r0) || ((r0-r) < -orig_vol[index]))	orig_vol[index] = r - r0;
				}
			}


	int j , k, index;
	for(k = 0; k < dim[2]; k++)
		for(j = 0; j < dim[1]; j++)
			for(i = 0; i < dim[0]; i++) {
				index = k*dim[0]*dim[1] + j*dim[0] + i;
				if(orig_vol[index] > 0.0) {
					orig_vol[index] = -1.0;
				}
			}


	// AChR (2) -- membrane
	float r0, r;
	int j, k,index;
	r0 = center - 0.001f;
	for(k = 0; k < dim[2]; k++)
		for(j = 0; j < dim[1]; j++)
			for(i = 0; i < dim[0]; i++) {
				index = k*dim[0]*dim[1] + j*dim[0] + i;
				orig_vol[index] = orig_vol[index] - 40.0f;
				r = (float) sqrt((i-center)*(i-center) + (j-center)*(j-center) + (k-center)*(k-center));
				//if(k < center || r > r0 || (orig_vol[index] > 40.0f && r < r0))	orig_vol[index] = 1.0f;
				//else	orig_vol[index] = -1.0f;
				if(orig_vol[index] < 0.0) {
					if((r>r0) || ((r0-r) < -orig_vol[index]))	orig_vol[index] = r - r0;
				}
				if(k < center)	orig_vol[index] = 1.0f;
			}

	// AChR (3) -- exterior
	float r0, r;
	int j, k,index;
	r0 = center - 0.001f;
	for(k = 0; k < dim[2]; k++)
		for(j = 0; j < dim[1]; j++)
			for(i = 0; i < dim[0]; i++) {
				index = k*dim[0]*dim[1] + j*dim[0] + i;
				orig_vol[index] = orig_vol[index] - 40.0f;
				r = (float) sqrt((i-center)*(i-center) + (j-center)*(j-center) + (k-center)*(k-center));
				if(orig_vol[index] < 0.0) {
					if((r>r0) || ((r0-r) < -orig_vol[index]))	orig_vol[index] = r - r0;
				}
				if(k > center)	orig_vol[index] = 1.0f;
			}

	// add an outer box/sphere
	float r0, r;
	//int j, k,index;
	r0 = center - 0.001f;
	for(k = 0; k < dim[2]; k++)
		for(j = 0; j < dim[1]; j++)
			for(i = 0; i < dim[0]; i++) {
				index = k*dim[0]*dim[1] + j*dim[0] + i;
				
				orig_vol[index] = orig_vol[index] - 12.0f;
				
				// box
				//if(i<3 || i>126 || j<3 || j>126 || k<3 || k>126)  orig_vol[index] = 1.0f;

				// sphere
				r = (float) sqrt((i-center)*(i-center) + (j-center)*(j-center) + (k-center)*(k-center));
				if(orig_vol[index] < 0.0) {
					if((r>r0) || ((r0-r) < -orig_vol[index]))	orig_vol[index] = r - r0;
				}
			}


	// RDV_segmentation
	for(k = 0; k < dim[2]; k++)
		for(j = 0; j < dim[1]; j++)
			for(i = 0; i < dim[0]; i++) {
				index = k*dim[0]*dim[1] + j*dim[0] + i;
				//if(orig_vol[index] > -137.0f && orig_vol[index] < -100.0f)
				if(orig_vol[index] < -170.0f)
					orig_vol[index] = (orig_vol[index] + 170.0f + 7.5f)/(255.0f - 170.0f);
				else if(orig_vol[index] < -100.0f)
					orig_vol[index] = (orig_vol[index] + 100.0f)/(170.0f - 100.0f);
				else
					orig_vol[index] = orig_vol[index]/100.0f;
			}

	
	//mache
	float r0, r;
	int j, k,index;
	r0 = center - 0.001f;
	for(k = 0; k < dim[2]; k++)
		for(j = 0; j < dim[1]; j++)
			for(i = 0; i < dim[0]; i++) {
				index = k*dim[0]*dim[1] + j*dim[0] + i;
				
				// cavity modification
				//if(i>=127 && i<=134 && j>=176 && j<=183 && k>=110 && k<=117)	//211
				//	orig_vol[index] = orig_vol[index] + 0.25;
				//if(i>=127 && i<=133 && j>=170 && j<=177 && k>=111 && k<=119)	//161
				//	orig_vol[index] = orig_vol[index] + 0.455;
				
				// add an outer sphere for mAChE
				//orig_vol[index] = 0.5 - orig_vol[index];
				//r = (float) sqrt((i-center)*(i-center) + (j-center)*(j-center) + (k-center)*(k-center));
				//if(orig_vol[index] < 0.0) {
				//	if((r>r0) || ((r0-r) < -orig_vol[index]))	orig_vol[index] = r - r0;
				//}
				
				// add an outer sphere for mAChE4
				orig_vol[index] = orig_vol[index] - 12.0f;
				r = (float) sqrt((i-center)*(i-center) + (j-center)*(j-center) + (k-center)*(k-center));
				if(orig_vol[index] < 0.0) {
					if((r>r0) || ((r0-r) < -orig_vol[index]))	orig_vol[index] = r - r0;
				}
			}


	// two spheres
	float r0, r1, r;
	int j, k, index;
	r0 = 34.0f;//center - 0.001f;
	r1 = 0.0f;//32.0;
	for(k = 0; k < dim[2]; k++)
		for(j = 0; j < dim[1]; j++)
			for(i = 0; i < dim[0]; i++) {
				index = k*dim[0]*dim[1] + j*dim[0] + i;
				r = (float) sqrt((i-center)*(i-center) + (j-center)*(j-center) + (k-center)*(k-center));
				if(r <= r1) orig_vol[index] = r1 - r;
				else if (r >= r0) orig_vol[index] = r - r0;
				else {
					if(r - r1 > r0 - r)	orig_vol[index] = r - r0;
					else orig_vol[index] = r1 - r;
				}
			}
	
	// one spheres
	float r0, r;
	int j, k, index;
	r0 = 10.0f;//center - 0.001f;
	for(k = 0; k < dim[2]; k++)
		for(j = 0; j < dim[1]; j++)
			for(i = 0; i < dim[0]; i++) {
				index = k*dim[0]*dim[1] + j*dim[0] + i;
				r = (float) sqrt((i-center)*(i-center) + (j-center)*(j-center) + (k-center)*(k-center));
				orig_vol[index] = r - r0;
			}

}*/


int Octree::get_neighbor_bit(int id,int level)
{
	int x,y,z;
	int reg_bitmask=0;
	
	octcell2xyz(id,x,y,z,level);
	reg_bitmask|=(is_refined(x  , y-1, z-1, level)<<0);
	reg_bitmask|=(is_refined(x-1, y  , z-1, level)<<1);
	reg_bitmask|=(is_refined(x  , y  , z-1, level)<<2);
	reg_bitmask|=(is_refined(x+1, y  , z-1, level)<<3);
	reg_bitmask|=(is_refined(x  , y+1, z-1, level)<<4);
	
	reg_bitmask|=(is_refined(x-1, y-1, z, level)<<5);
	reg_bitmask|=(is_refined(x  , y-1, z, level)<<6);
	reg_bitmask|=(is_refined(x+1, y-1, z, level)<<7);
	reg_bitmask|=(is_refined(x-1, y  , z, level)<<8);
	reg_bitmask|=(is_refined(x+1, y  , z, level)<<9);
	reg_bitmask|=(is_refined(x-1, y+1, z, level)<<10);
	reg_bitmask|=(is_refined(x  , y+1, z, level)<<11);
	reg_bitmask|=(is_refined(x+1, y+1, z, level)<<12);
	
	reg_bitmask|=(is_refined(x  , y-1, z+1, level)<<13);
	reg_bitmask|=(is_refined(x-1, y  , z+1, level)<<14);
	reg_bitmask|=(is_refined(x  , y  , z+1, level)<<15);
	reg_bitmask|=(is_refined(x+1, y  , z+1, level)<<16);
	reg_bitmask|=(is_refined(x  , y+1, z+1, level)<<17);
	
	return reg_bitmask;
}


int Octree::is_refined2(int x, int y, int z, int level)
{
	int idx=0;
	int res;
	
	res = (1<<level) ;
	
	if (x<0 || y<0 || z<0) return 0;
	if (x>=res || y>=res || z>=res) return 0;
	
	idx = level_id[level] + x + y*res + z*res*res;
	if (oct_array[idx].refine_flag==1) return 1;
	else return 0;
	//return oct_array[idx].refine_flag;
}

int Octree::is_refined(int x, int y, int z, int level)
{
	int idx=0;
	int res;
	
	res = (1<<level) ;
	
	if (x<0 || y<0 || z<0) return 1;
	if (x>=res || y>=res || z>=res) return 1;
	
	idx = level_id[level] + x + y*res + z*res*res;
	if (oct_array[idx].refine_flag==0) return 0;
	else return 1;
	//return oct_array[idx].refine_flag;
}

int Octree::get_depth(int res)
{
	int i=0;
	while (1) {
		if (res<=(1<<i)+1) break;
		i++;
	}
	if (res!=(1<<i)+1) {
		printf("unsupported resolution : %d\n",res);
		//		exit(0);
	}
	return i;
}

int Octree::get_octcell_num(int depth)
{
	int num=0;
	for (int i=0;i<=depth;i++) {
		num+=(1<<(i*3));
	}
	return num;
}

int Octree::get_level(int oc_id)
{
	int num=0;
	int i=0;
	while (1) {
		num+=(1<<(i*3));
		if (num>oc_id) break;
		i++;
	}
	return i;
}

void Octree::octcell2xyz(int oc_id,int& x,int& y,int& z,int level)
{
	int idx;
	int lres;
	
	idx=oc_id-level_id[level];
	lres=level_res[level];
	x=idx%lres;
	y=(idx/lres)%lres;
	z=idx/(lres*lres);
}

int Octree::xyz2octcell(int x,int y,int z,int level)
{
	int lres;
	
	lres=level_res[level];
	if (x<0 || y<0 || z<0 || x>=lres || y>=lres || z>=lres) return -1;
	return level_id[level]+x+y*lres+z*lres*lres;
}

int Octree::child(int oc_id, int level, int i)
{
	int x,y,z;
	int ret_idx;
	octcell2xyz(oc_id,x,y,z,level);
	
	switch (i) {
	case 0 :
		ret_idx= xyz2octcell(x*2,y*2,z*2,level+1);
		break;
	case 1 :
		ret_idx= xyz2octcell(x*2+1,y*2,z*2,level+1);
		break;
	case 2 :
		ret_idx= xyz2octcell(x*2,y*2+1,z*2,level+1);
		break;
	case 3 :
		ret_idx= xyz2octcell(x*2+1,y*2+1,z*2,level+1);
		break;
	case 4 :
		ret_idx= xyz2octcell(x*2,y*2,z*2+1,level+1);
		break;
	case 5 :
		ret_idx= xyz2octcell(x*2+1,y*2,z*2+1,level+1);
		break;
	case 6 :
		ret_idx= xyz2octcell(x*2,y*2+1,z*2+1,level+1);
		break;
	case 7 :
		ret_idx= xyz2octcell(x*2+1,y*2+1,z*2+1,level+1);
		break;
	}
	return ret_idx;
}


int Octree::xyz2vtx(int x, int y, int z)
{
	return x+y*dim[0]+z*dim[0]*dim[1];
}


void Octree::idx2vtx(int oc_id, int level, int* vtx)
{
	int x,y,z;
	int x0,x2,y0,y2,z0,z2;
	int cell_size;
	
	cell_size = (dim[0]-1)/(1<<level);
	
	octcell2xyz(oc_id,x,y,z,level);
	
	x0=x*cell_size; x2=x0+cell_size;
	y0=y*cell_size; y2=y0+cell_size;
	z0=z*cell_size; z2=z0+cell_size;
	
	vtx[0]=xyz2vtx(x0,y0,z0);
	vtx[1]=xyz2vtx(x2,y0,z0);
	vtx[2]=xyz2vtx(x2,y0,z2);
	vtx[3]=xyz2vtx(x0,y0,z2);
	vtx[4]=xyz2vtx(x0,y2,z0);
	vtx[5]=xyz2vtx(x2,y2,z0);
	vtx[6]=xyz2vtx(x2,y2,z2);
	vtx[7]=xyz2vtx(x0,y2,z2);
	
}


int Octree::is_intersect(int in_idx ,float isovalue , float* in_value 
						 ,int& iv_idx , int i, int j, int k , int level , int faceidx, geoframe& geofrm)
{
	float pt[3], norm[3];
	
	EdgeInfo *ei = &edge_dir[faceidx][in_idx];
	float f1,f2;
	
	f1=in_value[ei->d1];
	f2=in_value[ei->d2];
	
	if (!((f1<=isovalue && f2>=isovalue)||(f1>=isovalue && f2<=isovalue))) {
		return 0; // false : not intersect
	} 
	if (f1==f2) return 0;
	
	switch (ei->dir) {
	case 0:
		interpRect3Dpts_x(2*i+ei->di,2*j+ei->dj,2*k+ei->dk, 
			f1, f2,isovalue,pt, norm, level+1);
		break;
		
	case 1:
		interpRect3Dpts_y(2*i+ei->di,2*j+ei->dj,2*k+ei->dk,
			f1, f2, isovalue,pt, norm, level+1);
		break;
		
	case 2:
		interpRect3Dpts_z(2*i+ei->di,2*j+ei->dj,2*k+ei->dk,
			f1, f2, isovalue,pt, norm, level+1);
		break;
	}
	assert(ei->dir==0 || ei->dir==1 || ei->dir==2);
	
	iv_idx = geofrm.AddVert(pt, norm);
	return 1;
	
}


void Octree::interpRect3Dpts_x(int i1, int j1, int k1, 
							   float d1, float d2, float val, float *pt,float *norm,int level)
{
	double ival, gval;
	int cell_size;
	float g1[3], g2[3];
	
	cell_size = (dim[0]-1)/(1<<level);
	
	ival = (val - d1)/(d2 - d1);
	pt[0] = orig[0] + span[0]*(i1 + ival)*cell_size;
	pt[1] = orig[1] + span[1]*j1*cell_size;
	pt[2] = orig[2] + span[2]*k1*cell_size;
	
	//#ifdef GRAD
	getVertGrad( (i1+ival)*cell_size   , j1*cell_size , k1*cell_size , g1);
	getVertGrad( (i1+ival)*cell_size+1 , j1*cell_size , k1*cell_size , g2);
	gval=(double)((i1+ival)*cell_size)-(double)((int)((i1+ival)*cell_size));
	
	norm[0] = g1[0]*(1-gval) + g2[0]*gval;
	norm[1] = g1[1]*(1-gval) + g2[1]*gval;
	norm[2] = g1[2]*(1-gval) + g2[2]*gval;
	
	float len= sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
	norm[0]/=len;
	norm[1]/=len;
	norm[2]/=len;
	
	/*
	
	  // grad addtion
	  norm[0] = g1[0]*(1.0-ival) + g2[0]*ival;
	  norm[1] = g1[1]*(1.0-ival) + g2[1]*ival;
	  norm[2] = g1[2]*(1.0-ival) + g2[2]*ival;
	*/
	//#endif
	
}

void Octree::interpRect3Dpts_y(int i1, int j1, int k1, 
							   float d1, float d2, float val, float *pt,float *norm,int level)
							   
{
	double ival,gval;
	int cell_size;
	float g1[3],g2[3];
	cell_size = (dim[0]-1)/(1<<level);
	
	ival = (val - d1)/(d2 - d1);
	pt[0] = orig[0] + span[0]*i1*cell_size;
	pt[1] = orig[1] + span[1]*(j1 + ival)*cell_size;
	pt[2] = orig[2] + span[2]*k1*cell_size;
	
	//#ifdef GRAD
	// grad addtion
	
	getVertGrad( i1*cell_size   , (j1+ival)*cell_size , k1*cell_size , g1);
	getVertGrad( i1*cell_size , (j1+ival)*cell_size+1 , k1*cell_size , g2);
	gval=(double)((j1+ival)*cell_size)-(double)((int)((j1+ival)*cell_size));
	
	

	norm[0] = g1[0]*(1-gval) + g2[0]*gval;
	norm[1] = g1[1]*(1-gval) + g2[1]*gval;
	norm[2] = g1[2]*(1-gval) + g2[2]*gval;

	float len= sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
	norm[0]/=len;
	norm[1]/=len;
	norm[2]/=len;

	
	//   norm[0] = g1[0]*(1.0-ival) + g2[0]*ival;
	//   norm[1] = g1[1]*(1.0-ival) + g2[1]*ival;
	//   norm[2] = g1[2]*(1.0-ival) + g2[2]*ival;
	//#endif
}


void Octree::interpRect3Dpts_z(int i1, int j1, int k1,
							   float d1, float d2, float val, float *pt,float *norm,int level)
							   
{
	double ival,gval;
	int cell_size;
	
	float g1[3],g2[3];
	
	cell_size = (dim[0]-1)/(1<<level);
	
	ival = (val - d1)/(d2 - d1);
	pt[0] = orig[0] + span[0]*i1*cell_size;
	pt[1] = orig[1] + span[1]*j1*cell_size;
	pt[2] = orig[2] + span[2]*(k1+ival)*cell_size;
	
	//#ifdef GRAD
	
	getVertGrad( i1*cell_size   , j1*cell_size , (k1+ival)*cell_size , g1);
	getVertGrad( i1*cell_size , j1*cell_size , (k1+ival)*cell_size +1, g2);
	gval=(double)((k1+ival)*cell_size)-(double)((int)((k1+ival)*cell_size));

	
	norm[0] = g1[0]*(1-gval) + g2[0]*gval;
	norm[1] = g1[1]*(1-gval) + g2[1]*gval;
	norm[2] = g1[2]*(1-gval) + g2[2]*gval;

	float len= sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
	norm[0]/=len;
	norm[1]/=len;
	norm[2]/=len;

	//#endif
}

float Octree::compute_error(int oc_idx,int level,float& min, float& max)
{
	int x,y,z,x_orig,y_orig,z_orig;
	float sum=0;
	int vtx[8];
	float temp1,temp2,temp3,temp4,temp5,temp6,diff;
	float orig_val,interp_val;
	int cell_size;
	float val[8];
	float x_ratio,y_ratio,z_ratio;
	min=FLOAT_MAXIMUM;
	max=FLOAT_MINIMUM;
	
	cell_size = (dim[0]-1)/(1<<level);
	octcell2xyz(oc_idx,x_orig,y_orig,z_orig,level);
	
	x_orig=x_orig*cell_size;
	y_orig=y_orig*cell_size;
	z_orig=z_orig*cell_size;
	
	idx2vtx(oc_idx,level,vtx);
	int i;
	for (i=0;i<8;i++) {
		val[i]=orig_vol[vtx[i]];
	}
	
	for (z=z_orig;z<=z_orig+cell_size;z++) {
		for (y=y_orig;y<=y_orig+cell_size;y++) {
			for (x=x_orig;x<=x_orig+cell_size;x++) {
				orig_val=orig_vol[xyz2vtx(x,y,z)];
				
				if (min>orig_val) min=orig_val;
				if (max<orig_val) max=orig_val;
				
				x_ratio=((float)(x-x_orig))/((float)cell_size);
				y_ratio=((float)(y-y_orig))/((float)cell_size);
				z_ratio=((float)(z-z_orig))/((float)cell_size);
				
				temp1 = val[0] + (val[1]-val[0])*x_ratio;
				temp2 = val[2] + (val[3]-val[2])*x_ratio;
				temp3 = val[4] + (val[5]-val[4])*x_ratio;
				temp4 = val[6] + (val[7]-val[6])*x_ratio;
				temp5 = temp1  + (temp2-temp1)*y_ratio;
				temp6 = temp3  + (temp4-temp3)*y_ratio;
				interp_val = temp5  + (temp6-temp5)*z_ratio;
				
				diff = (orig_val>interp_val) ? orig_val-interp_val : interp_val-orig_val;
				sum+=diff*diff;
				//sum+=diff;
			}
		}
	}
	
	if (level==oct_depth) return 0;
	
	//if (sum<=0.01) return 0.01;
	//return oct_depth-get_level(oc_idx);
	//printf("sum : %f\n",sum);
	return sum;
}


void Octree::getCellValues(int oc_id,int level,float* val)
{
	int vtx[8];
	idx2vtx(oc_id,level,vtx);
	int i;
	for (i=0;i<8;i++) {
		//val[i]=orig_vol[vtx[i]];
		val[i] = orig_vol[vtx[i]];
	}
}


float Octree::getValue(int i, int j, int k)
{
	/*
	int cell_size=(dim[0]-1)/(1<<level);
	return orig_vol[i*cell_size + j*cell_size*dim[0] + k*cell_size*dim[0]*dim[1]];
	*/

	return orig_vol[i + j*dim[0] + k*dim[0]*dim[1]];
}

//------------------------------------------------------------------
// compute vertex gradient
// Three method are availble now:
// 1. simple finite difference
// 2. use the limit B-spline approximation
// 3. B-spline convolution
//-------------------------------------------------------------------
void Octree::getVertGrad(int i, int j, int k, float g[3]) 
{

	if(flag_normal == 1) {	// central difference
		int lres=dim[0];
		if (i==0) {
			g[0] = getValue(i+1, j, k) - getValue(i, j, k);
		}
		else if (i>=lres-1) {
			g[0] = getValue(i, j, k) - getValue(i-1, j, k);
		}
		else {
			g[0] = (getValue(i+1, j, k) - getValue(i-1, j, k)) * 0.5;
		}

		if (j==0) {
			g[1] = getValue(i, j+1, k) - getValue(i, j, k);
		}
		else if (j>=lres-1) {
			g[1] = getValue(i, j, k) - getValue(i, j-1, k);
		}
		else {
			g[1] = (getValue(i, j+1, k) - getValue(i, j-1, k)) * 0.5;
		}
		
		if (k==0) {
			g[2] = getValue(i, j, k+1) - getValue(i, j, k);
		}
		else if (k>=lres-1) {
			g[2] = getValue(i, j, k) - getValue(i, j, k-1);
		}
		else {
			g[2] = (getValue(i, j, k+1) - getValue(i, j, k-1)) * 0.5;
		}
	}
	else if(flag_normal == 0) {	// Bspline Convolution
		float v[27];
		int ix[3], iy[3], iz[3];
		int l, m, n, t = 0;

		ix[0] = (i-1 >= 0)? i-1:0;
		ix[1] = i;
		ix[2] = (i+1 < dim[0])? i+1:i;
		iy[0] = (j-1 >= 0)? j-1:0;
		iy[1] = j;
		iy[2] = (j+1 < dim[1])? j+1:j;
		iz[0] = (k-1 >= 0)? k-1:0;
		iz[1] = k;
		iz[2] = (k+1 < dim[2])? k+1:k;

		for (n = 0; n < 3; n++) {
			for (m = 0; m < 3; m++) {
				for (l = 0; l < 3; l++) {
					v[t] = getValue(ix[l], iy[m], iz[n]);
					t++;
				}
			}
		}

		g[0] = g[1] = g[2] = 0;
		for (l = 0; l < 27; l++) {
			g[0] += x_grad_mask[l]*v[l];
			g[1] += y_grad_mask[l]*v[l];
			g[2] += z_grad_mask[l]*v[l];
		}

		g[0] /= span[0];
		g[1] /= span[1];
		g[2] /= span[2];
	}
	else {	// Bspline Interpolation
	  LBIE::GradientAtPoint(&(BSplineCoeff[0]), (float)i, (float)j, (float)k, dim[0], dim[1], dim[2], g);
	  g[0]=g[0]/span[0];
	  g[1]=g[1]/span[1];
	  g[2]=g[2]/span[2];
	}

}


int Octree::is_skipcell(int oc_id)
{
	if ((minmax[oc_id].max > iso_val) && (minmax[oc_id].min < iso_val)) 
		return 0;
	else return 1;
}

int Octree::is_skipcell_in(int oc_id)
{
	if ((minmax[oc_id].max > iso_val_in) && (minmax[oc_id].min < iso_val_in)) 
		return 0;
	else return 1;
}

int Octree::is_skipcell_interval(int oc_id)
{
	if ((minmax[oc_id].max > iso_val && minmax[oc_id].min < iso_val) ||
		(minmax[oc_id].max > iso_val_in && minmax[oc_id].min < iso_val_in)) 
		return 0;
	else return 1;
}

void Octree::geometric_flow(geoframe& geofrm) {

	switch (flag_type) {
	case 0:							// triangular mesh
		geometric_flow_tri(geofrm);
		break;

	case 1:							// tetrahedral mesh
		geometric_flow_tet(geofrm);
		break;

	case 2:							// quadrilateral  mesh
		geometric_flow_quad(geofrm);
		break;

	case 3:							// hex mesh
		geometric_flow_hex(geofrm);
		break;

	case 4:							// triangular mesh
		geometric_flow_tri(geofrm);
		break;

	case 5:							// tetraheral mesh
		geometric_flow_tet(geofrm);
		break;
	}

}

// quality improvement with geometric flow -- tri mesh
void Octree::geometric_flow_tri(geoframe& geofrm) {

	int i, j, k, v0, v1, v2, itri, index, maxIndex;
	float mv0[3], mv1[3], mv2[3], p[3], sum_0, sum_1[3], delta_t, area, t;
	int **neighbor;
	//	FILE *output;

	//	output = fopen("qimprove_temp", "w");

	//	fprintf(output, "geometric flow tri\n");

	delta_t = 0.01f;
	maxIndex = 100;

	neighbor = 0;
	neighbor = new int*[geofrm.numverts];

	for(i = 0; i < geofrm.numverts; i++)	neighbor[i] = 0;
	for(i = 0; i < geofrm.numverts; i++)	neighbor[i] = new int[10];
	for(i = 0; i < geofrm.numverts; i++) {
		for (j = 0; j < 10; j++)	neighbor[i][j] = -1;
	}

	for (i = 0; i < geofrm.numtris; i++) {
		v0 = geofrm.triangles[i][0];
		v1 = geofrm.triangles[i][1];
		v2 = geofrm.triangles[i][2];

		for(j = 0; j < 10; j++) {
			if(neighbor[v0][j] == -1) {neighbor[v0][j] = i;	break;}
		}
		for(j = 0; j < 10; j++) {
			if(neighbor[v1][j] == -1) {neighbor[v1][j] = i;	break;}
		}
		for(j = 0; j < 10; j++) {
			if(neighbor[v2][j] == -1) {neighbor[v2][j] = i;	break;}
		}
	}

	//	fprintf(output, "%d %d\n", geofrm.numverts, geofrm.numtris);
	for(index = 0; index < maxIndex; index++) {

	  //		fprintf(output, "%d ", index);
	  //		if(index % 20 == 0 && index > 0) fprintf(output, "\n");

		for(i = 0; i < geofrm.numverts; i++) {

			// calculate the mass center p[3]
			sum_0 = 0.0f;
			for(j = 0; j < 3; j++) {p[j] = 0.0f; sum_1[j] = 0.0f;}
			for(j = 0; j < 10; j++) {
				itri = neighbor[i][j];
				if(itri == -1) break;
				v0 = geofrm.triangles[itri][0];	v1 = geofrm.triangles[itri][1];
				v2 = geofrm.triangles[itri][2];	

				for(k = 0; k < 3; k++) {
					mv0[k] = geofrm.verts[v0][k];	mv1[k] = geofrm.verts[v1][k];
					mv2[k] = geofrm.verts[v2][k];	
				}

				area = area_tri(mv0, mv1, mv2);
				sum_0 += area;
				for(k = 0; k < 3; k++) {
					p[k] += (mv0[k] + mv1[k] + mv2[k])*area/3.0f;
				}

			}
			for(j = 0; j < 3; j++) p[j] /= sum_0;

			// calculate the new position in tangent direction
			// xi+1 = xi + delta_t*((m-xi) - (n, m-xi)n))
			t = 0.0f;
			for(j = 0; j < 3; j++) {
				p[j] -= geofrm.verts[i][j];
				t += p[j]*geofrm.normals[i][j];
			}
			//fprintf(output, "\n%f %f %f\n", geofrm.verts[i][0], geofrm.verts[i][1], geofrm.verts[i][2]);
			for(j = 0; j < 3; j++) {
				geofrm.verts[i][j] += delta_t*p[j];								// mass center
				//geofrm.verts[i][j] += delta_t*(p[j] - t*geofrm.normals[i][j]);		// tangent movement
			}
			//fprintf(output, "%f %f %f\n", geofrm.verts[i][0], geofrm.verts[i][1], geofrm.verts[i][2]);
		}

	} // end of index loop

	for(i = 0; i < geofrm.numverts; i++)	delete [] neighbor[i];
	delete [] neighbor;

	//	fclose(output);
}

// quality improvement with geometric flow -- tet mesh
void Octree::geometric_flow_tet(geoframe& geofrm) {

	int i, j, k, v0, v1, v2, v3, **neighbor;
	int itet, index, maxIndex, bool_0, bool_1, bool_2, bool_3;
	float mv0[3], mv1[3], mv2[3], mv3[3], p[3], volume;
	float sum_0, delta_t, area, t;
	//	FILE *output;

	//	output = fopen("qimprove_temp", "w");

	//	fprintf(output, "geometric flow tet\n");

	delta_t = 0.01f;
	maxIndex = 100;

	neighbor = 0;
	neighbor = new int*[geofrm.numverts];

	for(i = 0; i < geofrm.numverts; i++)	neighbor[i] = 0;
	for(i = 0; i < geofrm.numverts; i++)	neighbor[i] = new int[50];
	for(i = 0; i < geofrm.numverts; i++) {
		for (j = 0; j < 50; j++)	neighbor[i][j] = -1;
	}

	for (i = 0; i < geofrm.numtris/4; i++) {
		v0 = geofrm.triangles[4*i][0];
		v2 = geofrm.triangles[4*i][1];
		v1 = geofrm.triangles[4*i][2];
		v3 = geofrm.triangles[4*i+1][2];

		for(j = 0; j < 50; j++) {
			if(neighbor[v0][j] == -1) {neighbor[v0][j] = i;	break;}
		}
		for(j = 0; j < 50; j++) {
			if(neighbor[v1][j] == -1) {neighbor[v1][j] = i;	break;}
		}
		for(j = 0; j < 50; j++) {
			if(neighbor[v2][j] == -1) {neighbor[v2][j] = i;	break;}
		}
		for(j = 0; j < 50; j++) {
			if(neighbor[v3][j] == -1) {neighbor[v3][j] = i;	break;}
		}
	}

	//	fprintf(output, "%d %d\n", geofrm.numverts, geofrm.numtris/4);
	for(index = 0; index < maxIndex; index++) {

	  //		fprintf(output, "%d ", index);
	  //		if(index % 20 == 0 && index > 0) fprintf(output, "\n");

		for(i = 0; i < geofrm.numverts; i++) {

			if(geofrm.bound_sign[i] == 1) {					// boundary vertex

				// calculate the mass center p[3]
				sum_0 = 0.0f;
				for(j = 0; j < 3; j++) p[j] = 0.0f;
				for(j = 0; j < 50; j++) {
					itet = neighbor[i][j];
					if(itet == -1) break;
					v0 = geofrm.triangles[4*itet][0];	v2 = geofrm.triangles[4*itet][1];
					v1 = geofrm.triangles[4*itet][2];	v3 = geofrm.triangles[4*itet+1][2];

					bool_0 = geofrm.bound_sign[v0] + geofrm.bound_sign[v2] + geofrm.bound_sign[v1];
					bool_1 = geofrm.bound_sign[v0] + geofrm.bound_sign[v1] + geofrm.bound_sign[v3];
					bool_2 = geofrm.bound_sign[v1] + geofrm.bound_sign[v2] + geofrm.bound_sign[v3];
					bool_3 = geofrm.bound_sign[v2] + geofrm.bound_sign[v0] + geofrm.bound_sign[v3];

					if((i == v0 || i == v1 || i == v2) && bool_0 == 3) {
						for(k = 0; k < 3; k++) {
							mv0[k] = geofrm.verts[v0][k];	mv1[k] = geofrm.verts[v2][k];
							mv2[k] = geofrm.verts[v1][k];		
						}
						area = area_tri(mv0, mv1, mv2);
						sum_0 += area;
						for(k = 0; k < 3; k++) {
							p[k] += (mv0[k] + mv1[k] + mv2[k])*area/3.0f;
						}
					}
					if((i == v0 || i == v1 || i == v3) && bool_1 == 3) { 
						for(k = 0; k < 3; k++) {
							mv0[k] = geofrm.verts[v0][k];	mv1[k] = geofrm.verts[v1][k];
							mv2[k] = geofrm.verts[v3][k];	
						}
						area = area_tri(mv0, mv1, mv2);
						sum_0 += area;
						for(k = 0; k < 3; k++) {
							p[k] += (mv0[k] + mv1[k] + mv2[k])*area/3.0f;
						}
					}
					if((i == v1 || i == v2 || i == v3) && bool_2 == 3) {
						for(k = 0; k < 3; k++) {
							mv0[k] = geofrm.verts[v1][k];	mv1[k] = geofrm.verts[v2][k];
							mv2[k] = geofrm.verts[v3][k];	
						}
						area = area_tri(mv0, mv1, mv2);
						sum_0 += area;
						for(k = 0; k < 3; k++) {
							p[k] += (mv0[k] + mv1[k] + mv2[k])*area/3.0f;
						}
					}
					if((i == v2 || i == v0 || i == v3) && bool_3 == 3) {
						for(k = 0; k < 3; k++) {
							mv0[k] = geofrm.verts[v2][k];	mv1[k] = geofrm.verts[v0][k];
							mv2[k] = geofrm.verts[v3][k];	
						}
						area = area_tri(mv0, mv1, mv2);
						sum_0 += area;
						for(k = 0; k < 3; k++) {
							p[k] += (mv0[k] + mv1[k] + mv2[k])*area/3.0f;
						}
					}

				}
				
				for(j = 0; j < 3; j++) p[j] /= sum_0;
				if(sum_0 == 0.0f) {printf("s%d ", i); continue;}

				// calculate the new position in tangent direction
				// xi+1 = xi + delta_t*((m-xi) - (n, m-xi)n))
				t = 0.0f;
				for(j = 0; j < 3; j++) {
					p[j] -= geofrm.verts[i][j];
					t += p[j]*geofrm.normals[i][j];
				}
				for(j = 0; j < 3; j++) {
					//geofrm.verts[i][j] += delta_t*p[j];							// mass center
					geofrm.verts[i][j] += delta_t*(p[j] - t*geofrm.normals[i][j]);		// tangent movement
				}
				
			}
			else {									// interior vertex

				// calculate the volume center p[3]
				sum_0 = 0.0f;
				for(j = 0; j < 3; j++) p[j] = 0.0f;
				for(j = 0; j < 50; j++) {
					itet = neighbor[i][j];
					if(itet == -1) break;
					v0 = geofrm.triangles[4*itet][0];	v2 = geofrm.triangles[4*itet][1];
					v1 = geofrm.triangles[4*itet][2];	v3 = geofrm.triangles[4*itet+1][2];

					for(k = 0; k < 3; k++) {
						mv0[k] = geofrm.verts[v0][k];	mv1[k] = geofrm.verts[v1][k];
						mv2[k] = geofrm.verts[v2][k];	mv3[k] = geofrm.verts[v3][k];
					}
					volume = volume_tet(mv0, mv1, mv2, mv3);
					if(volume < 0.00001f) continue;
					sum_0 += (float) volume;
					for(k = 0; k < 3; k++) {
						p[k] += (mv0[k] + mv1[k] + mv2[k] + mv3[k])*volume/4.0f;
					}
				}

				for(j = 0; j < 3; j++) p[j] /= sum_0;
				//if(sum_0 <= 0.0f || sum_0 > 10000.0f) {printf("v%d %f\n", i, sum_0); continue;}
				if(sum_0 < 0.00001f) continue;
				for(j = 0; j < 3; j++) {
					geofrm.verts[i][j] += delta_t*(p[j] - geofrm.verts[i][j]);					// volume center
				}

			}

		} // end of i loop

	} // end of index loop

	for(i = 0; i < geofrm.numverts; i++)	delete [] neighbor[i];
	delete [] neighbor;

	//	fclose(output);
}

// quality improvement with geometric flow -- quad mesh
void Octree::geometric_flow_quad(geoframe& geofrm) {

	int i, j, k, v0, v1, v2, v3;
	int iquad, index, maxIndex, **neighbor;
	float mv0[3], mv1[3], mv2[3], mv3[3], p[3];
	float sum_0, sum_1[3], delta_t, area, t;
	//	FILE *output;

	//	output = fopen("qimprove_temp", "w");

	//	fprintf(output, "geometric flow quad\n");

	delta_t = 0.01f;
	maxIndex = 100;

	neighbor = 0;
	neighbor = new int*[geofrm.numverts];

	for(i = 0; i < geofrm.numverts; i++)	neighbor[i] = 0;
	for(i = 0; i < geofrm.numverts; i++)	neighbor[i] = new int[10];
	for(i = 0; i < geofrm.numverts; i++) {
		for (j = 0; j < 10; j++)	neighbor[i][j] = -1;
	}

	for (i = 0; i < geofrm.numquads; i++) {
		v0 = geofrm.quads[i][0];
		v1 = geofrm.quads[i][1];
		v2 = geofrm.quads[i][2];
		v3 = geofrm.quads[i][3];

		for(j = 0; j < 10; j++) {
			if(neighbor[v0][j] == -1) {neighbor[v0][j] = i;	break;}
		}
		for(j = 0; j < 10; j++) {
			if(neighbor[v1][j] == -1) {neighbor[v1][j] = i;	break;}
		}
		for(j = 0; j < 10; j++) {
			if(neighbor[v2][j] == -1) {neighbor[v2][j] = i;	break;}
		}
		for(j = 0; j < 10; j++) {
			if(neighbor[v3][j] == -1) {neighbor[v3][j] = i;	break;}
		}
	}

	//	fprintf(output, "%d %d\n", geofrm.numverts, geofrm.numquads);
	for(index = 0; index < maxIndex; index++) {

	  //		fprintf(output, "%d ", index);
	  //		if(index % 20 == 0 && index > 0) fprintf(output, "\n");

		for(i = 0; i < geofrm.numverts; i++) {

			// calculate the mass center p[3]
			sum_0 = 0.0f;
			for(j = 0; j < 3; j++) {p[j] = 0.0f; sum_1[j] = 0.0f;}
			for(j = 0; j < 10; j++) {
				iquad = neighbor[i][j];
				if(iquad == -1) break;
				v0 = geofrm.quads[iquad][0];	v1 = geofrm.quads[iquad][1];
				v2 = geofrm.quads[iquad][2];	v3 = geofrm.quads[iquad][3];

				for(k = 0; k < 3; k++) {
					mv0[k] = geofrm.verts[v0][k];	mv1[k] = geofrm.verts[v1][k];
					mv2[k] = geofrm.verts[v2][k];	mv3[k] = geofrm.verts[v3][k];
				}

				area = area_quad(mv0, mv1, mv2, mv3);
				sum_0 += area;
				for(k = 0; k < 3; k++) {
					p[k] += (mv0[k] + mv1[k] + mv2[k] + mv3[k])*area/4.0f;
				}

			}
			for(j = 0; j < 3; j++) p[j] /= sum_0;

			// calculate the new position in tangent direction
			// xi+1 = xi + delta_t*((m-xi) - (n, m-xi)n))
			t = 0.0f;
			for(j = 0; j < 3; j++) {
				p[j] -= geofrm.verts[i][j];
				//t += p[j]*geofrm.normals[i][j];
			}
			for(j = 0; j < 3; j++) {
				geofrm.verts[i][j] += delta_t*p[j];								// mass center
				//geofrm.verts[i][j] += delta_t*(p[j] - t*geofrm.normals[i][j]);		// tangent movement
			}
		}

	} // end of index loop

	for(i = 0; i < geofrm.numverts; i++)	delete [] neighbor[i];
	delete [] neighbor;

	//	fclose(output);
}

// quality improvement with geometric flow -- hex mesh
void Octree::geometric_flow_hex(geoframe& geofrm) {

	int i, j, k, v0, v1, v2, v3, v4, v5, v6, v7, **neighbor;
	int ihex, index, maxIndex, bool_0, bool_1, bool_2, bool_3, bool_4, bool_5;
	float mv0[3], mv1[3], mv2[3], mv3[3], p[3];
	float mv4[3], mv5[3], mv6[3], mv7[3], volume;
	float sum_0, delta_t, area, t;
	//FILE *output;

	//output = fopen("qimprove_temp", "w");

	//fprintf(output, "geometric flow hex\n");

	delta_t = 0.01f;
	maxIndex = 100;

	neighbor = 0;
	neighbor = new int*[geofrm.numverts];

	for(i = 0; i < geofrm.numverts; i++)	neighbor[i] = 0;
	for(i = 0; i < geofrm.numverts; i++)	neighbor[i] = new int[50];
	for(i = 0; i < geofrm.numverts; i++) {
		for (j = 0; j < 50; j++)	neighbor[i][j] = -1;
	}

	for (i = 0; i < geofrm.numquads/6; i++) {
		v0 = geofrm.quads[6*i][0];		v1 = geofrm.quads[6*i][1];
		v2 = geofrm.quads[6*i][2];		v3 = geofrm.quads[6*i][3];
		v4 = geofrm.quads[6*i+1][1];	v5 = geofrm.quads[6*i+1][0];
		v6 = geofrm.quads[6*i+1][3];	v7 = geofrm.quads[6*i+1][2];

		for(j = 0; j < 50; j++) {
			if(neighbor[v0][j] == -1) {neighbor[v0][j] = i;	break;}
		}
		for(j = 0; j < 50; j++) {
			if(neighbor[v1][j] == -1) {neighbor[v1][j] = i;	break;}
		}
		for(j = 0; j < 50; j++) {
			if(neighbor[v2][j] == -1) {neighbor[v2][j] = i;	break;}
		}
		for(j = 0; j < 50; j++) {
			if(neighbor[v3][j] == -1) {neighbor[v3][j] = i;	break;}
		}
		for(j = 0; j < 50; j++) {
			if(neighbor[v4][j] == -1) {neighbor[v4][j] = i;	break;}
		}
		for(j = 0; j < 50; j++) {
			if(neighbor[v5][j] == -1) {neighbor[v5][j] = i;	break;}
		}
		for(j = 0; j < 50; j++) {
			if(neighbor[v6][j] == -1) {neighbor[v6][j] = i;	break;}
		}
		for(j = 0; j < 50; j++) {
			if(neighbor[v7][j] == -1) {neighbor[v7][j] = i;	break;}
		}
	}

	//fprintf(output, "%d %d\n", geofrm.numverts, geofrm.numquads/6);
	for(index = 0; index < maxIndex; index++) {

	  //fprintf(output, "%d ", index);
	  //		if(index % 20 == 0 && index > 0) fprintf(output, "\n");

		for(i = 0; i < geofrm.numverts; i++) {

			if(geofrm.bound_sign[i] == 1) {					// boundary vertex

				// calculate the mass center p[3]
				sum_0 = 0.0f;
				for(j = 0; j < 3; j++) p[j] = 0.0f;
				for(j = 0; j < 50; j++) {
					ihex = neighbor[i][j];
					if(ihex == -1) break;
					v0 = geofrm.quads[6*ihex][0];		v1 = geofrm.quads[6*ihex][1];
					v2 = geofrm.quads[6*ihex][2];		v3 = geofrm.quads[6*ihex][3];
					v4 = geofrm.quads[6*ihex+1][1];		v5 = geofrm.quads[6*ihex+1][0];
					v6 = geofrm.quads[6*ihex+1][3];		v7 = geofrm.quads[6*ihex+1][2];

					bool_0 = geofrm.bound_sign[v0] + geofrm.bound_sign[v1] + geofrm.bound_sign[v2] + geofrm.bound_sign[v3];
					bool_1 = geofrm.bound_sign[v4] + geofrm.bound_sign[v5] + geofrm.bound_sign[v6] + geofrm.bound_sign[v7];
					bool_2 = geofrm.bound_sign[v0] + geofrm.bound_sign[v3] + geofrm.bound_sign[v7] + geofrm.bound_sign[v4];
					bool_3 = geofrm.bound_sign[v1] + geofrm.bound_sign[v2] + geofrm.bound_sign[v6] + geofrm.bound_sign[v5];
					bool_4 = geofrm.bound_sign[v0] + geofrm.bound_sign[v4] + geofrm.bound_sign[v5] + geofrm.bound_sign[v1];
					bool_5 = geofrm.bound_sign[v3] + geofrm.bound_sign[v7] + geofrm.bound_sign[v6] + geofrm.bound_sign[v2];

					if((i == v0 || i == v1 || i == v2 || i == v3) && bool_0 == 4) { // down
						for(k = 0; k < 3; k++) {
							mv0[k] = geofrm.verts[v0][k];	mv1[k] = geofrm.verts[v1][k];
							mv2[k] = geofrm.verts[v2][k];	mv3[k] = geofrm.verts[v3][k];
						}
						area = area_quad(mv0, mv1, mv2, mv3);
						sum_0 += area;
						for(k = 0; k < 3; k++) {
							p[k] += (mv0[k] + mv1[k] + mv2[k] + mv3[k])*area/4.0f;
						}
					}
					if((i == v4 || i == v5 || i == v6 || i == v7) && bool_1 == 4) { // up
						for(k = 0; k < 3; k++) {
							mv0[k] = geofrm.verts[v4][k];	mv1[k] = geofrm.verts[v5][k];
							mv2[k] = geofrm.verts[v6][k];	mv3[k] = geofrm.verts[v7][k];
						}
						area = area_quad(mv0, mv1, mv2, mv3);
						sum_0 += area;
						for(k = 0; k < 3; k++) {
							p[k] += (mv0[k] + mv1[k] + mv2[k] + mv3[k])*area/4.0f;
						}
					}
					if((i == v0 || i == v3 || i == v7 || i == v4) && bool_2 == 4) {	// left
						for(k = 0; k < 3; k++) {
							mv0[k] = geofrm.verts[v0][k];	mv1[k] = geofrm.verts[v3][k];
							mv2[k] = geofrm.verts[v7][k];	mv3[k] = geofrm.verts[v4][k];
						}
						area = area_quad(mv0, mv1, mv2, mv3);
						sum_0 += area;
						for(k = 0; k < 3; k++) {
							p[k] += (mv0[k] + mv1[k] + mv2[k] + mv3[k])*area/4.0f;
						}
					}
					if((i == v1 || i == v2 || i == v6 || i == v5) && bool_3 == 4) {	// right
						for(k = 0; k < 3; k++) {
							mv0[k] = geofrm.verts[v1][k];	mv1[k] = geofrm.verts[v2][k];
							mv2[k] = geofrm.verts[v6][k];	mv3[k] = geofrm.verts[v5][k];
						}
						area = area_quad(mv0, mv1, mv2, mv3);
						sum_0 += area;
						for(k = 0; k < 3; k++) {
							p[k] += (mv0[k] + mv1[k] + mv2[k] + mv3[k])*area/4.0f;
						}
					}
					if((i == v0 || i == v4 || i == v5 || i == v1) && bool_4 == 4) {	// front
						for(k = 0; k < 3; k++) {
							mv0[k] = geofrm.verts[v0][k];	mv1[k] = geofrm.verts[v4][k];
							mv2[k] = geofrm.verts[v5][k];	mv3[k] = geofrm.verts[v1][k];
						}
						area = area_quad(mv0, mv1, mv2, mv3);
						sum_0 += area;
						for(k = 0; k < 3; k++) {
							p[k] += (mv0[k] + mv1[k] + mv2[k] + mv3[k])*area/4.0f;
						}
					}
					if((i == v3 || i == v7 || i == v6 || i == v2) && bool_5 == 4) {	// back
						for(k = 0; k < 3; k++) {
							mv0[k] = geofrm.verts[v3][k];	mv1[k] = geofrm.verts[v7][k];
							mv2[k] = geofrm.verts[v6][k];	mv3[k] = geofrm.verts[v2][k];
						}
						area = area_quad(mv0, mv1, mv2, mv3);
						sum_0 += area;
						for(k = 0; k < 3; k++) {
							p[k] += (mv0[k] + mv1[k] + mv2[k] + mv3[k])*area/4.0f;
						}
					}

				}
				
				for(j = 0; j < 3; j++) p[j] /= sum_0;
				if(sum_0 == 0.0f) {printf("s%d ", i); continue;}

				// calculate the new position in tangent direction
				// xi+1 = xi + delta_t*((m-xi) - (n, m-xi)n))
				t = 0.0f;
				for(j = 0; j < 3; j++) {
					p[j] -= geofrm.verts[i][j];
					t += p[j]*geofrm.normals[i][j];
				}
				for(j = 0; j < 3; j++) {
					//geofrm.verts[i][j] += delta_t*p[j];							// mass center
					geofrm.verts[i][j] += delta_t*(p[j] - t*geofrm.normals[i][j]);		// tangent movement
				}
				
			}
			else {									// interior vertex

				// calculate the volume center p[3]
				sum_0 = 0.0f;
				for(j = 0; j < 3; j++) p[j] = 0.0f;
				for(j = 0; j < 50; j++) {
					ihex = neighbor[i][j];
					if(ihex == -1) break;
					v0 = geofrm.quads[6*ihex][0];		v1 = geofrm.quads[6*ihex][1];
					v2 = geofrm.quads[6*ihex][2];		v3 = geofrm.quads[6*ihex][3];
					v4 = geofrm.quads[6*ihex+1][1];		v5 = geofrm.quads[6*ihex+1][0];
					v6 = geofrm.quads[6*ihex+1][3];		v7 = geofrm.quads[6*ihex+1][2];

					for(k = 0; k < 3; k++) {
						mv0[k] = geofrm.verts[v0][k];	mv1[k] = geofrm.verts[v1][k];
						mv2[k] = geofrm.verts[v2][k];	mv3[k] = geofrm.verts[v3][k];
						mv4[k] = geofrm.verts[v4][k];	mv5[k] = geofrm.verts[v5][k];
						mv6[k] = geofrm.verts[v6][k];	mv7[k] = geofrm.verts[v7][k];
					}
					volume = volume_hex(mv0, mv1, mv2, mv3, mv4, mv5, mv6, mv7);
					if(volume < 0.00001f) continue;
					sum_0 += (float) volume;
					for(k = 0; k < 3; k++) {
						p[k] += (mv0[k] + mv1[k] + mv2[k] + mv3[k] + mv4[k] + mv5[k] + mv6[k] + mv7[k])*volume/8.0f;
					}
				}

				//if(sum_0 < 0.00001f || sum_0 > 10000.0f) {fprintf(output, "v%d %f\n", i, sum_0); continue;}
				if(sum_0 < 0.00001f) continue;
				for(j = 0; j < 3; j++) p[j] /= sum_0;
				for(j = 0; j < 3; j++) {
					geofrm.verts[i][j] += delta_t*(p[j] - geofrm.verts[i][j]);					// volume center
				}

			}

		} // end of i loop

	} // end of index loop

	for(i = 0; i < geofrm.numverts; i++)	delete [] neighbor[i];
	delete [] neighbor;

	//fclose(output);
}

// calculate normal at v0
void Octree::crossproduct(float v0[3], float v1[3], float v2[3], float* normal) {

	float v01[3], v02[3], g;

	int i;
	for(i = 0; i < 3; i++) {
		v01[i] = v1[i] - v0[i];
		v02[i] = v2[i] - v0[i];
	}

	normal[0] = v01[1]*v02[2] - v02[1]*v01[2];
	normal[1] = v01[2]*v02[0] - v02[2]*v01[0];
	normal[2] = v01[0]*v02[1] - v02[0]*v01[1];


	g = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]; 
	for(i = 0; i < 3; i++) {
		normal[i] /= (float) sqrt(g);
	}

}

// calculate the area for a triangle
float Octree::area_tri(float v0[3], float v1[3], float v2[3]) {

	float a, b, c, p, area;
	int i;

	a = 0.0;		b = 0.0;		c = 0.0;
	for(i = 0; i < 3; i++) {
		a += (v1[i] - v0[i])*(v1[i] - v0[i]);
		b += (v2[i] - v1[i])*(v2[i] - v1[i]);
		c += (v0[i] - v2[i])*(v0[i] - v2[i]);
	}
	a = (float)sqrt(a);		b = (float)sqrt(b);		c = (float)sqrt(c);
	p = (a + b + c) / 2.0f;
	area = (float)sqrt(p * (p - a) * (p - b) * (p - c));

	return area;

}

// calculate the area for a quad
float Octree::area_quad(float v0[3], float v1[3], float v2 [3], float v3[3]) {

	int i, j;
	float u, v, su[3], sv[3], length2_su, length2_sv, suv, area;

	area = 0.0f;

	for(j = 0; j < 4; j++) {

		u = 0.5f - (float)sqrt(3)/6.0f;
		v = u;
		if(j == 1 || j == 3) u = 0.5f + (float)sqrt(3)/6.0f;
		if(j == 2 || j == 3) v = 0.5f + (float)sqrt(3)/6.0f;

		length2_su = 0.0f;	length2_sv = 0.0f;	suv = 0.0f;
		for(i = 0; i < 3; i++) {
			su[i] = (1-v)*(v1[i] - v0[i]) + v*(v2[i] - v3[i]);
			sv[i] = (1-u)*(v3[i] - v0[i]) + u*(v2[i] - v1[i]);
			length2_su += su[i]*su[i];
			length2_sv += sv[i]*sv[i];
			suv += su[i]*sv[i];
		}

		area += (float) sqrt(length2_su*length2_sv - suv*suv)/4.0f;

	}

	return area;

}

// calculate the volume for a tet
float Octree::volume_tet(float v0[3], float v1[3], float v2[3], float v3[3]) {

	float t, a, b, c, p, area, volume;
	float v01[3], v02[3], n[3], temp;

	int i;
	for(i = 0; i < 3; i++) {
		v01[i] = v1[i] - v0[i];
		v02[i] = v2[i] - v0[i];
	}

	n[0] = v01[1]*v02[2] - v01[2]*v02[1];
	n[1] = v01[2]*v02[0] - v01[0]*v02[2];
	n[2] = v01[0]*v02[1] - v01[1]*v02[0];

	temp = (float)sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
	n[0] = - n[0] / temp;
	n[1] = - n[1] / temp;
	n[2] = - n[2] / temp;

	t = (v2[0] - v3[0])*n[0] + (v2[1] - v3[1])*n[1] + (v2[2] - v3[2])*n[2];

	a = (float)sqrt((v1[0] - v0[0])*(v1[0] - v0[0]) + (v1[1] - v0[1])*(v1[1] - v0[1])
						+ (v1[2] - v0[2])*(v1[2] - v0[2]));
	b = (float)sqrt((v2[0] - v1[0])*(v2[0] - v1[0]) + (v2[1] - v1[1])*(v2[1] - v1[1])
						+ (v2[2] - v1[2])*(v2[2] - v1[2]));
	c = (float)sqrt((v2[0] - v0[0])*(v2[0] - v0[0]) + (v2[1] - v0[1])*(v2[1] - v0[1])
						+ (v2[2] - v0[2])*(v2[2] - v0[2]));
	p = a + b + c;
	area = 0.25f*(float)sqrt(p*(p - 2*a)*(p - 2*b)*(p - 2*c));
	volume = area * t / 3.0f;
	
	return volume;

}

// calculate the volume for a hex
float Octree::volume_hex(float v0[3], float v1[3], float v2[3], float v3[3], float v4[3], float v5[3], float v6[3], float v7[3]) {

	int i, j;
	float u, v, w, su[3], sv[3], sw[3], cuv[3], volume;

	volume = 0.0f;

	for(j = 0; j < 8; j++) {

		u = 0.5f - (float)sqrt(3)/6.0f;
		v = u;
		w = u;
		if(j == 1 || j == 3 || j == 5 || j == 7) u = 0.5f + (float)sqrt(3)/6.0f;
		if(j == 2 || j == 3 || j == 6 || j == 7) v = 0.5f + (float)sqrt(3)/6.0f;
		if(j == 4 || j == 5 || j == 6 || j == 7) w = 0.5f + (float)sqrt(3)/6.0f;
		for(i = 0; i < 3; i++) {
			su[i] = (1-v)*(1-w)*(v1[i] - v0[i]) + v*(1-w)*(v2[i] - v3[i]) + (1-v)*w*(v5[i] - v4[i]) + v*w*(v6[i] - v7[i]);
			sv[i] = (1-u)*(1-w)*(v3[i] - v0[i]) + u*(1-w)*(v2[i] - v1[i]) + (1-u)*w*(v7[i] - v4[i]) + u*w*(v6[i] - v5[i]);
			sw[i] = (1-u)*(1-v)*(v4[i] - v0[i]) + u*(1-v)*(v5[i] - v1[i]) + (1-u)*v*(v7[i] - v3[i]) + u*v*(v6[i] - v2[i]);
		}
		cuv[0] = su[1]*sv[2] - su[2]*sv[1];
		cuv[1] = su[2]*sv[0] - su[0]*sv[2];
		cuv[2] = su[0]*sv[1] - su[1]*sv[0];

		volume += (cuv[0]*sw[0]+cuv[1]*sw[1]+cuv[2]*sw[2])/8.0f;
	}

	return volume;

}

void Octree::edge_contraction(geoframe& geofrm) {

	int nv, ntet, new_nv, new_ntet, i, j, v0, v1, v2, v3;
	int *bound_sign, id_min, vv0, vv1, vv2, vv3, b_sign;
	int my_bool, index, maxIndex, **rep, nrep;
	int my_list[70], nlist, r_vert, bool_r, a_vert;
	int sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7;
	int sum8, sum9, sum10, sum11, sum12, sum13, sum14, sum15;
	float **vtx, **normal, vx, vy, vz, nx, ny, nz, sum_0, sum_1, sum_2;
	float mv0[3], mv1[3], mv2[3], mv3[3];
	double edge_ratio, new_edge_ratio, s;
	double e[6], e_min, e_max, aspect, aspect_min, aspect_max;
	FILE *input, /**output,*/ *fn2;

	input = fopen("input.raw", "w");
	fprintf(input, "%d %d\n", geofrm.numverts, geofrm.numtris/4);
	for (i = 0; i < geofrm.numverts; i++)
		fprintf(input,"%f %f %f %f %f %f %d\n", geofrm.verts[i][0], geofrm.verts[i][1], geofrm.verts[i][2], 
			geofrm.normals[i][0], geofrm.normals[i][1], geofrm.normals[i][2], geofrm.bound_sign[i]);
	for (i = 0; i < geofrm.numtris/4; i++)
		fprintf(input,"%d %d %d %d\n", geofrm.triangles[4*i][0], geofrm.triangles[4*i][1],
			geofrm.triangles[4*i][2], geofrm.triangles[4*i+1][2]);
	fclose(input);
	
	//output = fopen("qimprove_temp", "w");

	//fprintf(output, "edge contraction\n");

	maxIndex = 20;

	edge_ratio = 15.0f;				// the final target (e.g., 8.5)
	new_edge_ratio = 65.0f;	    // the current maximum
	bool_r = 0; r_vert = 108637;	// remove vertex 0/1
	a_vert = -1;					// average all neighbor's coord. -1  59157

	for (index = 0; index < maxIndex; index++) {

	  //fprintf(output, "Loop %d\n", index);

	input = fopen("input.raw", "r");

	fscanf(input, "%d %d\n", &nv, &ntet);
	bound_sign = new int[nv];
	rep = new int*[10];
	for(i = 0; i < 10; i++) rep[i] = new int[2];
	vtx = new float*[nv];
	for(i = 0; i < nv; i++)	vtx[i] = new float[3];
	normal = new float*[nv];
	for(i = 0; i < nv; i++)	normal[i] = new float[3];

	for (i = 0; i < nv; i++) {
		for (j = 0; j < 3; j++) {vtx[i][j] = 0.0f; normal[i][j] = 0.0f;}
		bound_sign[i] = 0;
	}
	for (i = 0; i < 10; i++) {
		rep[i][0] = -1;	rep[i][1] = -1;
	}
	for (i = 0; i < 70; i++) {
		my_list[i] = -1;
	}

	for (i = 0; i < nv; i++) {
		//fscanf(input,"%f %f %f\n", &vx , &vy , &vz);
		fscanf(input,"%f %f %f %f %f %f %d\n", &vx, &vy, &vz, &nx, &ny, &nz, &b_sign);
		vtx[i][0] = vx;
		vtx[i][1] = vy;
		vtx[i][2] = vz;

		normal[i][0] = nx;
		normal[i][1] = ny;
		normal[i][2] = nz;

		bound_sign[i] = b_sign;
	}

	aspect_min = 10.0;
	aspect_max = 0.0;

	//	fprintf(output, "tet read begin\n");
	//	fprintf(output, "%d %d\n", nv, ntet);

	fn2 = fopen("new_mesh2.raw", "w");

	sum0 = 0;	sum1 = 0;	sum2 = 0;	sum3 = 0;	sum4 = 0;	sum5 = 0;
	sum6 = 0;	sum7 = 0;	sum8 = 0;	sum9 = 0;	sum10 = 0;	sum11 = 0;
	sum12 = 0;	sum13 = 0;	sum14 = 0;	sum15 = 0;

	nrep = 0;
	new_ntet = 0;
	nlist = 0;

	for(i = 0; i < ntet; i++) {

	  //		if(i%10000 == 0) fprintf(output, "%d ", i);
		fscanf(input, "%d %d %d %d\n", &v0, &v1, &v2, &v3);

		if(v0 == a_vert || v1 == a_vert || v2 == a_vert || v3 == a_vert) {
			if(v0 != a_vert) {
				my_bool = 0;
				for(j = 0; j < nlist; j++) {
					if(v0 == my_list[j]) {my_bool++; break;}
				}
				if(my_bool == 0) {my_list[nlist] = v0; nlist++;}
			}

			if(v1 != a_vert) {
				my_bool = 0;
				for(j = 0; j < nlist; j++) {
					if(v1 == my_list[j]) {my_bool++; break;}
				}
				if(my_bool == 0) {my_list[nlist] = v1; nlist++;}
			}

			if(v2 != a_vert) {
				my_bool = 0;
				for(j = 0; j < nlist; j++) {
					if(v2 == my_list[j]) {my_bool++; break;}
				}
				if(my_bool == 0) {my_list[nlist] = v2; nlist++;}
			}

			if(v3 != a_vert) {
				my_bool = 0;
				for(j = 0; j < nlist; j++) {
					if(v3 == my_list[j]) {my_bool++; break;}
				}
				if(my_bool == 0) {my_list[nlist] = v3; nlist++;}
			}

		}
  
		for(j = 0; j < 6; j++) e[j] = 0.0;

		for(j = 0; j < 3; j++) {
			e[0] += (vtx[v1][j] - vtx[v0][j])*(vtx[v1][j] - vtx[v0][j]);//01
			e[1] += (vtx[v2][j] - vtx[v1][j])*(vtx[v2][j] - vtx[v1][j]);//12
			e[2] += (vtx[v2][j] - vtx[v0][j])*(vtx[v2][j] - vtx[v0][j]);//20
			e[3] += (vtx[v3][j] - vtx[v0][j])*(vtx[v3][j] - vtx[v0][j]);//30
			e[4] += (vtx[v3][j] - vtx[v1][j])*(vtx[v3][j] - vtx[v1][j]);//31
			e[5] += (vtx[v3][j] - vtx[v2][j])*(vtx[v3][j] - vtx[v2][j]);//32
		}

		e_min = e[0];
		e_max = e[0];
		id_min = 0;
		for(j = 1; j < 6; j++) {
			if(e[j] < e_min) {e_min = e[j]; id_min = j;}
			if(e[j] > e_max) e_max = e[j];
		}

		if(sqrt(e_min) == 0.0) {aspect = -1.0; printf("0000\n");}
		else { 
			aspect = sqrt(e_max/e_min);
			if(aspect < aspect_min) aspect_min = aspect;
			if(aspect > aspect_max) aspect_max = aspect;

			if(aspect < 2.5f) sum0++;
			else if(aspect < 5.0f) sum1++;
			else if(aspect < 7.5f) sum2++;
			else if(aspect < 8.5f) sum3++;
			else if(aspect < 12.5f) sum4++;
			else if(aspect < 15.0f) sum5++;
			else if(aspect < 17.5f) sum6++;
			else if(aspect < 20.0f) sum7++;
			else if(aspect < 22.5f) sum8++;
			else if(aspect < 25.0f) sum9++;
			else if(aspect < 27.5f) sum10++;
			else if(aspect < 30.0f) sum11++;
			else if(aspect < 32.5f) sum12++;
			else if(aspect < 45.0f) sum13++;
			else if(aspect < 50.0f) sum14++;
			else sum15++;
		}

		if(index == 0 && aspect_min > edge_ratio) new_edge_ratio = aspect_max - 0.00001f;

		if(aspect > new_edge_ratio && nrep == 0) {
		//if(i == 170502) {
		  //			fprintf(output, "i net = %d\n", i);
		  //			fprintf(output, "%d %d %d %d %f %f %f %d\n", bound_sign[v0], 
		  //bound_sign[v1], bound_sign[v2], bound_sign[v3], 
		  //				sqrt(e_min), sqrt(e_max), aspect, id_min);
		//			fprintf(output, "%f %f %f %f %f %f\n", sqrt(e[0]), sqrt(e[1]), sqrt(e[2]), sqrt(e[3]), sqrt(e[4]), sqrt(e[5]));
			if(id_min == 0) {vv0 = v0;  vv1 = v1;}
			else if(id_min == 1) {vv0 = v1;  vv1 = v2;}
			else if(id_min == 2) {vv0 = v2;  vv1 = v0;}
			else if(id_min == 3) {vv0 = v3;  vv1 = v0;}
			else if(id_min == 4) {vv0 = v3;  vv1 = v1;}
			else {vv0 = v3;  vv1 = v2;}

			//if(i == 168115) {vv0 = v0; vv1 = v3;}

			if(bound_sign[vv0] == 1 && bound_sign[vv1] == 0) {
				rep[nrep][0] = vv0;   rep[nrep][1] = vv1;}
			else if(bound_sign[vv0] == 0 && bound_sign[vv1] == 1) {
				rep[nrep][0] = vv1;   rep[nrep][1] = vv0;}
			//else if(bound_sign[v0] == 1 && bound_sign[v1] == 1 && bound_sign[v2] == 1 && bound_sign[v3] == 1){
			//else if(bound_sign[vv0] == 1 && bound_sign[vv1] == 1){
			else if(bound_sign[vv0] == bound_sign[vv1]){
				rep[nrep][0] = vv0;   rep[nrep][1] = vv1;
				if(vv0 > vv1) {rep[nrep][0] = vv1;   rep[nrep][1] = vv0;}}
			//else if(i == 50563) {rep[nrep][0] = v2;   rep[nrep][1] = v3;}
			else {	//if(bound_sign[vv0] == 0 && bound_sign[vv1] == 0) 
				fprintf(fn2, "%d %d %d %d\n", v0, v1, v2, v3);
				new_ntet++;
			}
			//else {
			//  rep[nrep][0] = vv0;   rep[nrep][1] = vv1;
			//  if(vv0 > vv1) {rep[nrep][0] = vv1;   rep[nrep][1] = vv0;}}

			//			fprintf(output, "%d %d %d %d %d\n", v0, v1, v2, v3, rep[nrep][0], rep[nrep][1]);
			if(rep[nrep][0] != -1 && rep[nrep][1] != -1) {
				for (j = 0; j < 3; j++) {
					vtx[rep[nrep][1]][j] = vtx[rep[nrep][0]][j];
					normal[rep[nrep][1]][j] = normal[rep[nrep][0]][j];
				}
				nrep++;
			}
		}
		else {
			if(aspect > 0.0) {
				fprintf(fn2, "%d %d %d %d\n", v0, v1, v2, v3);
				new_ntet++;
			}
		}
	}
	fclose(fn2);
	fclose(input);

	//	fprintf(output, "\ntet read end\n");

	//	fprintf(output, "%f %f\n", aspect_min, aspect_max);
	//	fprintf(output, "a_vert, nlist = %d %d\n", a_vert, nlist);

	input = fopen("input.raw", "w");
	fn2 = fopen("new_mesh2.raw", "r");

	new_nv = nv - nrep;
	if(bool_r == 1) new_nv--;

	fprintf(input, "%d %d\n", new_nv, new_ntet);
	for(i = 0; i < nv; i++) {
		my_bool = 0;
		for(j = 0; j < nrep; j++) {
			if(i == rep[j][1]) {my_bool++; break;}
		}
		if(bool_r == 1 && i == r_vert) my_bool++;

		if(i == a_vert) {
		  //			fprintf(output, "average %d %f %f %f %d\n", i, vtx[i][0], vtx[i][1], vtx[i][2], bound_sign[i]);
			sum_0 = 0.0;	sum_1 = 0.0;	sum_2 = 0.0;
			for(j = 0; j < nlist; j++) {
				sum_0 += vtx[my_list[j]][0];
				sum_1 += vtx[my_list[j]][1];
				sum_2 += vtx[my_list[j]][2];
			}
			vtx[i][0] = sum_0/nlist;
			vtx[i][1] = sum_1/nlist;
			vtx[i][2] = sum_2/nlist;
			//			fprintf(output, "average %d %f %f %f %d\n", i, vtx[i][0], vtx[i][1], vtx[i][2], bound_sign[i]);

			sum_0 = 0.0;	sum_1 = 0.0;	sum_2 = 0.0;
			for(j = 0; j < nlist; j++) {
				sum_0 += normal[my_list[j]][0];
				sum_1 += normal[my_list[j]][1];
				sum_2 += normal[my_list[j]][2];
			}
			//normal[i][0] = sum_0/nlist;
			//normal[i][1] = sum_1/nlist;
			//normal[i][2] = sum_2/nlist;
		}

		if(my_bool == 0) {
			fprintf(input, "%f %f %f %f %f %f %d\n", vtx[i][0], vtx[i][1], vtx[i][2], 
				normal[i][0], normal[i][1], normal[i][2], bound_sign[i]);
			//fprintf(input, "%f %f %f\n", vtx[i][0], vtx[i][1], vtx[i][2]);
		}
	}

	aspect_min = 10.0;
	aspect_max = 0.0;

	int new_ntet0 = new_ntet;
	for(i = 0; i < new_ntet; i++) {
		if(i%10000 == 0) printf("%d ", i);
		fscanf(fn2, "%d %d %d %d\n", &v0, &v1, &v2, &v3);
		
		for(j = 0; j < nrep; j++) {
			if(v0 == rep[j][1]) v0 = rep[j][0];
			if(v1 == rep[j][1]) v1 = rep[j][0];
			if(v2 == rep[j][1]) v2 = rep[j][0];
			if(v3 == rep[j][1]) v3 = rep[j][0];
		}
		  
		for(j = 0; j < 6; j++) e[j] = 0.0;

		for(j = 0; j < 3; j++) {
			e[0] += (vtx[v1][j] - vtx[v0][j])*(vtx[v1][j] - vtx[v0][j]);//01
			e[1] += (vtx[v2][j] - vtx[v1][j])*(vtx[v2][j] - vtx[v1][j]);//12
			e[2] += (vtx[v2][j] - vtx[v0][j])*(vtx[v2][j] - vtx[v0][j]);//20
			e[3] += (vtx[v3][j] - vtx[v0][j])*(vtx[v3][j] - vtx[v0][j]);//30
			e[4] += (vtx[v3][j] - vtx[v1][j])*(vtx[v3][j] - vtx[v1][j]);//31
			e[5] += (vtx[v3][j] - vtx[v2][j])*(vtx[v3][j] - vtx[v2][j]);//32
		}

		e_min = e[0];
		e_max = e[0];
		id_min = 0;
		for(j = 1; j < 6; j++) {
			if(e[j] < e_min)	e_min = e[j];
			if(e[j] > e_max)	{e_max = e[j]; id_min = j;}
		}

		if(sqrt(e_min) == 0.0) {aspect = -1.0;	printf("1111\n");	new_ntet0 = new_ntet0--;}
		else { 
			aspect = sqrt(e_max/e_min);
			if(aspect < aspect_min) aspect_min = aspect;
			if(aspect > aspect_max) aspect_max = aspect;
			
			vv0 = v0;	vv1 = v1;	vv2 = v2;	vv3 = v3;
			for(j = 0; j < nrep; j++) {
				if(v0 > rep[j][1]) vv0--;
				if(v1 > rep[j][1]) vv1--;
				if(v2 > rep[j][1]) vv2--;
				if(v3 > rep[j][1]) vv3--;
			}

			if(bool_r == 1) {
				if(v0 > r_vert) vv0--;
				if(v1 > r_vert) vv1--;
				if(v2 > r_vert) vv2--;
				if(v3 > r_vert) vv3--;
			}

			for(j = 0; j < 3; j++) {
				mv0[j] = vtx[vv0][j];	mv1[j] = vtx[vv1][j];
				mv2[j] = vtx[vv2][j];	mv3[j] = vtx[vv3][j];
			}
			s = geofrm.testTetrahedron1(mv0, mv1, mv2, mv3);
			if(6.0f*s*1.0e8 < 0.0f)
				fprintf(input, "%d %d %d %d\n", vv0, vv2, vv1, vv3);
			else
				fprintf(input, "%d %d %d %d\n", vv0, vv1, vv2, vv3);
		}
	}

	//	fprintf(output, "\n%f %f\n", aspect_min, aspect_max);
	//	fprintf(output, "\n%d %d %d\n", new_nv, new_ntet, new_ntet0);


	fclose(fn2);
	fclose(input);

	for(i = 0; i < 10; i++)	delete [] rep[i];
	delete [] rep;
	delete [] bound_sign;
	for(i = 0; i < nv; i++)	delete [] vtx[i];
	delete [] vtx;
	for(i = 0; i < nv; i++)	delete [] normal[i];
	delete [] normal;

	//	fprintf(output, "nrep = %d\n", nrep);
	if(nrep == 0) break;
	else new_edge_ratio = aspect_max - 1.0e-6;
	if(new_edge_ratio < edge_ratio) break;

	}	// end loop index

//	fprintf(output, "sum = %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7,
//		sum8, sum9, sum10, sum11, sum12, sum13, sum14, sum15, 
//		sum0+sum1+sum2+sum3+sum4+sum5+sum6+sum7+sum8+sum9+sum10+sum11+sum12+sum13+sum14+sum15);

	input = fopen("input.raw", "r");
	fscanf(input, "%d %d\n", &nv, &ntet);
	geofrm.numverts = nv;
	geofrm.numtris = ntet*4;
	for (i = 0; i < nv; i++) {
		fscanf(input,"%f %f %f %f %f %f %d\n", &vx, &vy, &vz, &nx, &ny, &nz, &b_sign);
		geofrm.verts[i][0] = vx;
		geofrm.verts[i][1] = vy;
		geofrm.verts[i][2] = vz;

		geofrm.normals[i][0] = nx;
		geofrm.normals[i][1] = ny;
		geofrm.normals[i][2] = nz;

		geofrm.bound_sign[i] = b_sign;
	}
	for (i = 0; i < ntet; i++) {
		fscanf(input,"%d %d %d %d\n", &v0, &v1, &v2, &v3);
		for(j = 0; j < 3; j++) {
			mv0[j] = geofrm.verts[v0][j];	mv1[j] = geofrm.verts[v1][j];
			mv2[j] = geofrm.verts[v2][j];	mv3[j] = geofrm.verts[v3][j];
		}
		s = geofrm.testTetrahedron1(mv0, mv1, mv2, mv3);
		if(6.0f*s*1.0e8 < 0.0f) {
			index = v1;	v1 = v2; v2 = index; 
		}
		geofrm.triangles[4*i][0] = v0;		geofrm.triangles[4*i][1] = v1;		geofrm.triangles[4*i][2] = v2; 
		geofrm.triangles[4*i+1][0] = v2;	geofrm.triangles[4*i+1][1] = v1;	geofrm.triangles[4*i+1][2] = v3; 
		geofrm.triangles[4*i+2][0] = v0;	geofrm.triangles[4*i+2][1] = v2;	geofrm.triangles[4*i+2][2] = v3; 
		geofrm.triangles[4*i+3][0] = v0;	geofrm.triangles[4*i+3][1] = v3;	geofrm.triangles[4*i+3][2] = v1; 
	}
	fclose(input);

//	fclose(output);

}

// sign == 0 -- improve the Joe-Liu parameter
// sign == 1 -- improve the minimal volume
void Octree::smoothing_joeliu_volume(geoframe& geofrm, int sign) {

	int maxIndex = 20;

	int nv, ntet, i, j, v0, v1, v2, v3;
	int my_bool, index;
	int my_list[200], nlist, a_vert, id_min;
	int sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7;
	int sum8, sum9, sum10, sum11, sum12, sum13, sum14, sum15;
	float sum_0, sum_1, sum_2, mv0[3], mv1[3], mv2[3], mv3[3];
	double e[6], sum_e, e_min, e_max, *aspectArr;
	double aspect_min, aspect_max, s, biggest_aspect_min = 0.0;
	int optRunNum;

	//	FILE *output;

	//	output = fopen("qimprove_temp", "w");

	nv = geofrm.numverts;
	ntet = geofrm.numtris/4;

	aspectArr = new double[ntet];

	sum0 = 0;	sum1 = 0;	sum2 = 0;	sum3 = 0;
	sum4 = 0;	sum5 = 0;	sum6 = 0;	sum7 = 0;
	sum8 = 0;	sum9 = 0;	sum10 = 0;	sum11 = 0;
	sum12 = 0;	sum13 = 0;	sum14 = 0;	sum15 = 0;

	for(i = 0; i < ntet; i++) {
		v0 = geofrm.triangles[4*i][0];
		v1 = geofrm.triangles[4*i][1];
		v2 = geofrm.triangles[4*i][2];
		v3 = geofrm.triangles[4*i+1][2];

		for(j = 0; j < 6; j++) e[j] = 0.0;

		for(j = 0; j < 3; j++) {
			e[0] += (geofrm.verts[v1][j] - geofrm.verts[v0][j])*(geofrm.verts[v1][j] - geofrm.verts[v0][j]);//01
			e[1] += (geofrm.verts[v2][j] - geofrm.verts[v1][j])*(geofrm.verts[v2][j] - geofrm.verts[v1][j]);//12
			e[2] += (geofrm.verts[v2][j] - geofrm.verts[v0][j])*(geofrm.verts[v2][j] - geofrm.verts[v0][j]);//20
			e[3] += (geofrm.verts[v3][j] - geofrm.verts[v0][j])*(geofrm.verts[v3][j] - geofrm.verts[v0][j]);//30
			e[4] += (geofrm.verts[v3][j] - geofrm.verts[v1][j])*(geofrm.verts[v3][j] - geofrm.verts[v1][j]);//31
			e[5] += (geofrm.verts[v3][j] - geofrm.verts[v2][j])*(geofrm.verts[v3][j] - geofrm.verts[v2][j]);//32
		}

		sum_e = 0.0f;
		for(j = 0; j < 6; j++) sum_e += e[j];

		for(j = 0; j < 3; j++) {
			mv0[j] = geofrm.verts[v0][j];	mv1[j] = geofrm.verts[v1][j];
			mv2[j] = geofrm.verts[v2][j];	mv3[j] = geofrm.verts[v3][j];
		}
		e_min = e[0];
		e_max = e[0];
		id_min = 0;
		for(j = 1; j < 6; j++) {
			if(e[j] < e_min) {e_min = e[j]; id_min = j;}
			if(e[j] > e_max) e_max = e[j];
		}

		s = geofrm.testTetrahedron1(mv0, mv1, mv2, mv3);

		if(sign == 0) 
			aspectArr[i] = 3.0*pow(4.0*fabs(s*6.0), 0.666667)/sum_e;
		else
			aspectArr[i] = s*6.0;

		if(aspectArr[i] > pow(10., -0.25)) sum15++;
		else if(aspectArr[i] > pow(10., -0.50)) sum14++;
		else if(aspectArr[i] > pow(10., -0.75)) sum13++;
		else if(aspectArr[i] > pow(10., -1.00)) sum12++;
		else if(aspectArr[i] > pow(10., -1.25)) sum11++;
		else if(aspectArr[i] > pow(10., -1.50)) sum10++;
		else if(aspectArr[i] > pow(10., -1.75)) sum9++;
		else if(aspectArr[i] > pow(10., -2.00)) sum8++;
		else if(aspectArr[i] > pow(10., -2.25)) sum7++;
		else if(aspectArr[i] > pow(10., -2.50)) sum6++;
		else if(aspectArr[i] > pow(10., -2.75)) sum5++;
		else if(aspectArr[i] > pow(10., -3.00)) sum4++;
		else if(aspectArr[i] > pow(10., -3.25)) sum3++;
		else if(aspectArr[i] > pow(10., -3.50)) sum2++;
		else if(aspectArr[i] > pow(10., -3.75)) sum1++;
		else sum0++;

		//if(aspectArr[i] < pow(10., -2.0)) sum1++;
		//else sum0++;
	}

	//	fprintf(output, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7, 
	//		sum8, sum9, sum10, sum11, sum12, sum13, sum14, sum15, 
	//		sum0+sum1+sum2+sum3+sum4+sum5+sum6+sum7+sum8+sum9+sum10+sum11+sum12+sum13+sum14+sum15);

	//fprintf(output, "%d %d %d %f\n", sum0, sum1, sum0+sum1, (float)sum1/(float)(sum0+sum1));

	srand( (unsigned)time( NULL ) );
	int minAspIndx;
	for(i = 0; i < ntet; i++) 
	{
		if(i == 0) {
			aspect_min = aspectArr[i];
			aspect_max = aspectArr[i];
			minAspIndx = 0;
		}
		if(aspectArr[i] < aspect_min) {
			aspect_min = aspectArr[i];
			minAspIndx = i;
		}
		if(aspectArr[i] > aspect_max) 	
			aspect_max = aspectArr[i];
	}
	a_vert = -1;
	int av0, av1, av2, av3;
	av0 = geofrm.triangles[4*minAspIndx][0];
	av1 = geofrm.triangles[4*minAspIndx][1];
	av2 = geofrm.triangles[4*minAspIndx][2];
	av3 = geofrm.triangles[4*minAspIndx+1][2];
	//	fprintf(output, "%d %d %d %d %d ", minAspIndx, av0, av1, av2, av3);
	//	fprintf(output, "%d %d %d %d\n", geofrm.bound_sign[av0], geofrm.bound_sign[av1], geofrm.bound_sign[av2], 
	//				geofrm.bound_sign[av3]);

	int randomInt = rand();
	if(randomInt % 4 == 0) 
	{
		if(geofrm.bound_sign[av0] == 0 ) {a_vert = av0;}
		else if(geofrm.bound_sign[av1] == 0) {a_vert = av1;}
		else if(geofrm.bound_sign[av2] == 0 ) {a_vert = av2;}
		else if(geofrm.bound_sign[av3] == 0) {a_vert = av3;}
		else {a_vert = av0;}//-1
	}
	else if(randomInt % 4 == 1) 
	{
		if(geofrm.bound_sign[av1] == 0 ) {a_vert = av1;}
		else if(geofrm.bound_sign[av2] == 0) {a_vert = av2;}
		else if(geofrm.bound_sign[av3] == 0 ) {a_vert = av3;}
		else if(geofrm.bound_sign[av0] == 0) {a_vert = av0;}
		else {a_vert = av1;}//-1
	}
	else if(randomInt % 4 == 2) 
	{
		if(geofrm.bound_sign[av2] == 0 ) {a_vert = av2;}
		else if(geofrm.bound_sign[av3] == 0) {a_vert = av3;}
		else if(geofrm.bound_sign[av0] == 0 ) {a_vert = av0;}
		else if(geofrm.bound_sign[av1] == 0) {a_vert = av1;}
		else {a_vert = av2;}//-1
	}
	else  
	{
		if(geofrm.bound_sign[av3] == 0 ) {a_vert = av3;}
		else if(geofrm.bound_sign[av0] == 0) {a_vert = av0;}
		else if(geofrm.bound_sign[av1] == 0 ) {a_vert = av1;}
		else if(geofrm.bound_sign[av2] == 0) {a_vert = av2;}
		else {a_vert = av3;}//-1
	}
	//a_vert =  av2;
//	fprintf(output, "aspect_min, aspect_max, a_vert = %e %e %d\n", aspect_min, aspect_max, a_vert);

	for(index = 0; index < maxIndex; index++)
	{
	  //		fprintf(output, "Loop %d\n", index);

		for (i = 0; i < 200; i++) {
			my_list[i] = -1;
		}
		nlist = 0;
		for(i = 0; i < ntet; i++) {
			v0 = geofrm.triangles[4*i][0];
			v1 = geofrm.triangles[4*i][1];
			v2 = geofrm.triangles[4*i][2];
			v3 = geofrm.triangles[4*i+1][2];

			if(v0 == a_vert || v1 == a_vert || v2 == a_vert || v3 == a_vert) {
				if(v0 != a_vert) {
					my_bool = 0;
					for(j = 0; j < nlist; j++) {
						if(v0 == my_list[j]) {my_bool++; break;}
					}
					if(my_bool == 0) {my_list[nlist] = v0; nlist++;}
				}

				if(v1 != a_vert) {
					my_bool = 0;
					for(j = 0; j < nlist; j++) {
						if(v1 == my_list[j]) {my_bool++; break;}
					}
					if(my_bool == 0) {my_list[nlist] = v1; nlist++;}
				}

				if(v2 != a_vert) {
					my_bool = 0;
					for(j = 0; j < nlist; j++) {
						if(v2 == my_list[j]) {my_bool++; break;}
					}
					if(my_bool == 0) {my_list[nlist] = v2; nlist++;}
				}

				if(v3 != a_vert) {
					my_bool = 0;
					for(j = 0; j < nlist; j++) {
						if(v3 == my_list[j]) {my_bool++; break;}
					}
					if(my_bool == 0) {my_list[nlist] = v3; nlist++;}
				}
			}
		}

		//		fprintf(output, "nlist = %d\n", nlist);

		i = a_vert;
		//		fprintf(output, "average %d %f %f %f %d\n", i, geofrm.verts[i][0], geofrm.verts[i][1], geofrm.verts[i][2], geofrm.bound_sign[i]);
		sum_0 = 0.0;	sum_1 = 0.0;	sum_2 = 0.0;
		for(j = 0; j < nlist; j++) {
			sum_0 += geofrm.verts[my_list[j]][0];
			sum_1 += geofrm.verts[my_list[j]][1];
			sum_2 += geofrm.verts[my_list[j]][2];
		}
		//if(sign == 0) {
		geofrm.verts[i][0] = sum_0/nlist;
		geofrm.verts[i][1] = sum_1/nlist;
		geofrm.verts[i][2] = sum_2/nlist;
		//}
		//		fprintf(output, "average %d %f %f %f %d\n", i, geofrm.verts[i][0], geofrm.verts[i][1], geofrm.verts[i][2], geofrm.bound_sign[i]);

		for(i = 0; i < ntet; i++) {
			if(geofrm.triangles[4*i][0] == a_vert || geofrm.triangles[4*i][1] == a_vert ||geofrm.triangles[4*i][2] == a_vert ||geofrm.triangles[4*i+1][2] == a_vert)
			{
				v0 = geofrm.triangles[4*i][0];
				v1 = geofrm.triangles[4*i][1];
				v2 = geofrm.triangles[4*i][2];
				v3 = geofrm.triangles[4*i+1][2];
				for(j = 0; j < 6; j++) e[j] = 0.0;

				for(j = 0; j < 3; j++) {
					e[0] += (geofrm.verts[v1][j] - geofrm.verts[v0][j])*(geofrm.verts[v1][j] - geofrm.verts[v0][j]);//01
					e[1] += (geofrm.verts[v2][j] - geofrm.verts[v1][j])*(geofrm.verts[v2][j] - geofrm.verts[v1][j]);//12
					e[2] += (geofrm.verts[v2][j] - geofrm.verts[v0][j])*(geofrm.verts[v2][j] - geofrm.verts[v0][j]);//20
					e[3] += (geofrm.verts[v3][j] - geofrm.verts[v0][j])*(geofrm.verts[v3][j] - geofrm.verts[v0][j]);//30
					e[4] += (geofrm.verts[v3][j] - geofrm.verts[v1][j])*(geofrm.verts[v3][j] - geofrm.verts[v1][j]);//31
					e[5] += (geofrm.verts[v3][j] - geofrm.verts[v2][j])*(geofrm.verts[v3][j] - geofrm.verts[v2][j]);//32
				}

				sum_e = 0.0f;
				for(j = 0; j < 6; j++) sum_e += e[j];

				for(j = 0; j < 3; j++) {
					mv0[j] = geofrm.verts[v0][j];	mv1[j] = geofrm.verts[v1][j];
					mv2[j] = geofrm.verts[v2][j];	mv3[j] = geofrm.verts[v3][j];
				}
				e_min = e[0];
				e_max = e[0];
				id_min = 0;
				for(j = 1; j < 6; j++) {
					if(e[j] < e_min) {e_min = e[j]; id_min = j;}
					if(e[j] > e_max) e_max = e[j];
				}

				s = geofrm.testTetrahedron1(mv0, mv1, mv2, mv3);

				if(sign == 0) 
					aspectArr[i] = 3.0*pow(4.0*fabs(s*6.0), 0.666667)/sum_e;
				else
					aspectArr[i] = s*6.0;

				if(s*6.0*1.0e8 < 0.0) 
				{
					geofrm.triangles[4*i][2] = v1;		geofrm.triangles[4*i][1] = v2;
					geofrm.triangles[4*i+1][0] = v1;	geofrm.triangles[4*i+1][1] = v2;
					geofrm.triangles[4*i+2][1] = v1;	geofrm.triangles[4*i+3][2] = v2; 
				}
			}
		}
		for(i = 0; i < ntet; i++) 
		{
			if(i == 0) {
				aspect_min = aspectArr[i];
				aspect_max = aspectArr[i];
				minAspIndx = 0;
			}
			if(aspectArr[i] < aspect_min) {
				aspect_min = aspectArr[i];
				minAspIndx = i;
			}
			if(aspectArr[i] > aspect_max) 	
				aspect_max = aspectArr[i];
		}
		a_vert = -1;
		av0 = geofrm.triangles[4*minAspIndx][0];
		av1 = geofrm.triangles[4*minAspIndx][1];
		av2 = geofrm.triangles[4*minAspIndx][2];
		av3 = geofrm.triangles[4*minAspIndx+1][2];
		//		fprintf(output, "%d %d %d %d ", av0, av1, av2, av3);
		//		fprintf(output, "%d %d %d %d\n", geofrm.bound_sign[av0], geofrm.bound_sign[av1], geofrm.bound_sign[av2], 
		//				geofrm.bound_sign[av3]);

		randomInt = rand();
		if(randomInt % 4 == 0) 
		{
			if(geofrm.bound_sign[av0] == 0 ) {a_vert = av0;}
			else if(geofrm.bound_sign[av1] == 0) {a_vert = av1;}
			else if(geofrm.bound_sign[av2] == 0 ) {a_vert = av2;}
			else if(geofrm.bound_sign[av3] == 0) {a_vert = av3;}
			else {a_vert = av0;}//-1
		}
		else if(randomInt % 4 == 1) 
		{
			if(geofrm.bound_sign[av1] == 0 ) {a_vert = av1;}
			else if(geofrm.bound_sign[av2] == 0) {a_vert = av2;}
			else if(geofrm.bound_sign[av3] == 0 ) {a_vert = av3;}
			else if(geofrm.bound_sign[av0] == 0) {a_vert = av0;}
			else {a_vert = av1;}//-1
		}
		else if(randomInt % 4 == 2) 
		{
			if(geofrm.bound_sign[av2] == 0 ) {a_vert = av2;}
			else if(geofrm.bound_sign[av3] == 0) {a_vert = av3;}
			else if(geofrm.bound_sign[av0] == 0 ) {a_vert = av0;}
			else if(geofrm.bound_sign[av1] == 0) {a_vert = av1;}
			else {a_vert = av2;}//-1
		}
		else  
		{
			if(geofrm.bound_sign[av3] == 0 ) {a_vert = av3;}
			else if(geofrm.bound_sign[av0] == 0) {a_vert = av0;}
			else if(geofrm.bound_sign[av1] == 0 ) {a_vert = av1;}
			else if(geofrm.bound_sign[av2] == 0) {a_vert = av2;}
			else {a_vert = av3;}//-1
		}
		//		fprintf(output, "aspect_min, aspect_max, a_vert = %e %e %d\n", aspect_min, aspect_max, a_vert);
		if(biggest_aspect_min < aspect_min)
		{
			biggest_aspect_min = aspect_min;
			optRunNum = index;
		}
	}

//loop 2

//	fprintf(output, "biggest aspect_min = %e \n", biggest_aspect_min);
	for(i = 0; i < ntet; i++) {
		v0 = geofrm.triangles[4*i][0];
		v1 = geofrm.triangles[4*i][1];
		v2 = geofrm.triangles[4*i][2];
		v3 = geofrm.triangles[4*i+1][2];

		for(j = 0; j < 6; j++) e[j] = 0.0;

		for(j = 0; j < 3; j++) {
			e[0] += (geofrm.verts[v1][j] - geofrm.verts[v0][j])*(geofrm.verts[v1][j] - geofrm.verts[v0][j]);//01
			e[1] += (geofrm.verts[v2][j] - geofrm.verts[v1][j])*(geofrm.verts[v2][j] - geofrm.verts[v1][j]);//12
			e[2] += (geofrm.verts[v2][j] - geofrm.verts[v0][j])*(geofrm.verts[v2][j] - geofrm.verts[v0][j]);//20
			e[3] += (geofrm.verts[v3][j] - geofrm.verts[v0][j])*(geofrm.verts[v3][j] - geofrm.verts[v0][j]);//30
			e[4] += (geofrm.verts[v3][j] - geofrm.verts[v1][j])*(geofrm.verts[v3][j] - geofrm.verts[v1][j]);//31
			e[5] += (geofrm.verts[v3][j] - geofrm.verts[v2][j])*(geofrm.verts[v3][j] - geofrm.verts[v2][j]);//32
		}

		sum_e = 0.0f;
		for(j = 0; j < 6; j++) sum_e += e[j];

		for(j = 0; j < 3; j++) {
			mv0[j] = geofrm.verts[v0][j];	mv1[j] = geofrm.verts[v1][j];
			mv2[j] = geofrm.verts[v2][j];	mv3[j] = geofrm.verts[v3][j];
		}
		e_min = e[0];
		e_max = e[0];
		id_min = 0;
		for(j = 1; j < 6; j++) {
			if(e[j] < e_min) {e_min = e[j]; id_min = j;}
			if(e[j] > e_max) e_max = e[j];
		}

		s = geofrm.testTetrahedron1(mv0, mv1, mv2, mv3);

		if(sign == 0) 
			aspectArr[i] = 3.0*pow(4.0*fabs(s*6.0), 0.666667)/sum_e;
		else
			aspectArr[i] = s*6.0;
	}


	srand( (unsigned)time( NULL ) );
	for(i = 0; i < ntet; i++) 
	{
		if(i == 0) {
			aspect_min = aspectArr[i];
			aspect_max = aspectArr[i];
			minAspIndx = 0;
		}
		if(aspectArr[i] < aspect_min) {
			aspect_min = aspectArr[i];
			minAspIndx = i;
		}
		if(aspectArr[i] > aspect_max) 	
			aspect_max = aspectArr[i];
	}
	a_vert = -1;
	av0 = geofrm.triangles[4*minAspIndx][0];
	av1 = geofrm.triangles[4*minAspIndx][1];
	av2 = geofrm.triangles[4*minAspIndx][2];
	av3 = geofrm.triangles[4*minAspIndx+1][2];
//	fprintf(output, "%d %d %d %d ", av0, av1, av2, av3);
//	fprintf(output, "%d %d %d %d\n", geofrm.bound_sign[av0], geofrm.bound_sign[av1], geofrm.bound_sign[av2], 
//				geofrm.bound_sign[av3]);

	randomInt = rand();
	if(randomInt % 4 == 0) 
	{
		if(geofrm.bound_sign[av0] == 0 ) {a_vert = av0;}
		else if(geofrm.bound_sign[av1] == 0) {a_vert = av1;}
		else if(geofrm.bound_sign[av2] == 0 ) {a_vert = av2;}
		else if(geofrm.bound_sign[av3] == 0) {a_vert = av3;}
		else {a_vert = av0;}//-1
	}
	else if(randomInt % 4 == 1) 
	{
		if(geofrm.bound_sign[av1] == 0 ) {a_vert = av1;}
		else if(geofrm.bound_sign[av2] == 0) {a_vert = av2;}
		else if(geofrm.bound_sign[av3] == 0 ) {a_vert = av3;}
		else if(geofrm.bound_sign[av0] == 0) {a_vert = av0;}
		else {a_vert = av1;}//-1
	}
	else if(randomInt % 4 == 2) 
	{
		if(geofrm.bound_sign[av2] == 0 ) {a_vert = av2;}
		else if(geofrm.bound_sign[av3] == 0) {a_vert = av3;}
		else if(geofrm.bound_sign[av0] == 0 ) {a_vert = av0;}
		else if(geofrm.bound_sign[av1] == 0) {a_vert = av1;}
		else {a_vert = av2;}//-1
	}
	else  
	{
		if(geofrm.bound_sign[av3] == 0 ) {a_vert = av3;}
		else if(geofrm.bound_sign[av0] == 0) {a_vert = av0;}
		else if(geofrm.bound_sign[av1] == 0 ) {a_vert = av1;}
		else if(geofrm.bound_sign[av2] == 0) {a_vert = av2;}
		else {a_vert = av3;}//-1
	}
//	fprintf(output, "aspect_min, aspect_max, a_vert = %e %e %d\n", aspect_min, aspect_max, a_vert);

	for(index = 0; index < maxIndex; index++)
	{
		if(aspect_min > biggest_aspect_min * 0.999)
			break;
		//		fprintf(output, "Loop %d\n", index);

		for (i = 0; i < 200; i++) {
			my_list[i] = -1;
		}
		nlist = 0;
		for(i = 0; i < ntet; i++) {
			v0 = geofrm.triangles[4*i][0];
			v1 = geofrm.triangles[4*i][1];
			v2 = geofrm.triangles[4*i][2];
			v3 = geofrm.triangles[4*i+1][2];

			if(v0 == a_vert || v1 == a_vert || v2 == a_vert || v3 == a_vert) {
				if(v0 != a_vert) {
					my_bool = 0;
					for(j = 0; j < nlist; j++) {
						if(v0 == my_list[j]) {my_bool++; break;}
					}
					if(my_bool == 0) {my_list[nlist] = v0; nlist++;}
				}

				if(v1 != a_vert) {
					my_bool = 0;
					for(j = 0; j < nlist; j++) {
						if(v1 == my_list[j]) {my_bool++; break;}
					}
					if(my_bool == 0) {my_list[nlist] = v1; nlist++;}
				}

				if(v2 != a_vert) {
					my_bool = 0;
					for(j = 0; j < nlist; j++) {
						if(v2 == my_list[j]) {my_bool++; break;}
					}
					if(my_bool == 0) {my_list[nlist] = v2; nlist++;}
				}

				if(v3 != a_vert) {
					my_bool = 0;
					for(j = 0; j < nlist; j++) {
						if(v3 == my_list[j]) {my_bool++; break;}
					}
					if(my_bool == 0) {my_list[nlist] = v3; nlist++;}
				}
			}
		}

		//		fprintf(output, "nlist = %d\n", nlist);

		i = a_vert;
		//		fprintf(output, "average %d %f %f %f %d\n", i, geofrm.verts[i][0], geofrm.verts[i][1], geofrm.verts[i][2], geofrm.bound_sign[i]);
		sum_0 = 0.0;	sum_1 = 0.0;	sum_2 = 0.0;
		for(j = 0; j < nlist; j++) {
			sum_0 += geofrm.verts[my_list[j]][0];
			sum_1 += geofrm.verts[my_list[j]][1];
			sum_2 += geofrm.verts[my_list[j]][2];
		}
		//if(sign == 0) {
		geofrm.verts[i][0] = sum_0/nlist;
		geofrm.verts[i][1] = sum_1/nlist;
		geofrm.verts[i][2] = sum_2/nlist;
		//}
		//		fprintf(output, "average %d %f %f %f %d\n", i, geofrm.verts[i][0], geofrm.verts[i][1], geofrm.verts[i][2], geofrm.bound_sign[i]);

		for(i = 0; i < ntet; i++) {
			if(geofrm.triangles[4*i][0] == a_vert || geofrm.triangles[4*i][1] == a_vert ||geofrm.triangles[4*i][2] == a_vert ||geofrm.triangles[4*i+1][2] == a_vert)
			{
				v0 = geofrm.triangles[4*i][0];
				v1 = geofrm.triangles[4*i][1];
				v2 = geofrm.triangles[4*i][2];
				v3 = geofrm.triangles[4*i+1][2];
				for(j = 0; j < 6; j++) e[j] = 0.0;

				for(j = 0; j < 3; j++) {
					e[0] += (geofrm.verts[v1][j] - geofrm.verts[v0][j])*(geofrm.verts[v1][j] - geofrm.verts[v0][j]);//01
					e[1] += (geofrm.verts[v2][j] - geofrm.verts[v1][j])*(geofrm.verts[v2][j] - geofrm.verts[v1][j]);//12
					e[2] += (geofrm.verts[v2][j] - geofrm.verts[v0][j])*(geofrm.verts[v2][j] - geofrm.verts[v0][j]);//20
					e[3] += (geofrm.verts[v3][j] - geofrm.verts[v0][j])*(geofrm.verts[v3][j] - geofrm.verts[v0][j]);//30
					e[4] += (geofrm.verts[v3][j] - geofrm.verts[v1][j])*(geofrm.verts[v3][j] - geofrm.verts[v1][j]);//31
					e[5] += (geofrm.verts[v3][j] - geofrm.verts[v2][j])*(geofrm.verts[v3][j] - geofrm.verts[v2][j]);//32
				}

				sum_e = 0.0f;
				for(j = 0; j < 6; j++) sum_e += e[j];

				for(j = 0; j < 3; j++) {
					mv0[j] = geofrm.verts[v0][j];	mv1[j] = geofrm.verts[v1][j];
					mv2[j] = geofrm.verts[v2][j];	mv3[j] = geofrm.verts[v3][j];
				}
				e_min = e[0];
				e_max = e[0];
				id_min = 0;
				for(j = 1; j < 6; j++) {
					if(e[j] < e_min) {e_min = e[j]; id_min = j;}
					if(e[j] > e_max) e_max = e[j];
				}

				s = geofrm.testTetrahedron1(mv0, mv1, mv2, mv3);

				if(sign == 0) 
					aspectArr[i] = 3.0*pow(4.0*fabs(s*6.0), 0.666667)/sum_e;
				else
					aspectArr[i] = s*6.0;

				if(s*6.0*1.0e8 < 0.0) 
				{
					geofrm.triangles[4*i][2] = v1;		geofrm.triangles[4*i][1] = v2;
					geofrm.triangles[4*i+1][0] = v1;	geofrm.triangles[4*i+1][1] = v2;
					geofrm.triangles[4*i+2][1] = v1;	geofrm.triangles[4*i+3][2] = v2; 
				}
			}
		}
		for(i = 0; i < ntet; i++) 
		{
			if(i == 0) {
				aspect_min = aspectArr[i];
				aspect_max = aspectArr[i];
				minAspIndx = 0;
			}
			if(aspectArr[i] < aspect_min) {
				aspect_min = aspectArr[i];
				minAspIndx = i;
			}
			if(aspectArr[i] > aspect_max) 	
				aspect_max = aspectArr[i];
		}
		a_vert = -1;
		av0 = geofrm.triangles[4*minAspIndx][0];
		av1 = geofrm.triangles[4*minAspIndx][1];
		av2 = geofrm.triangles[4*minAspIndx][2];
		av3 = geofrm.triangles[4*minAspIndx+1][3];
		//		fprintf(output, "%d %d %d %d ", av0, av1, av2, av3);
		//		fprintf(output, "%d %d %d %d\n", geofrm.bound_sign[av0], geofrm.bound_sign[av1], geofrm.bound_sign[av2], 
		//				geofrm.bound_sign[av3]);

		randomInt = rand();
		if(randomInt % 4 == 0) 
		{
			if(geofrm.bound_sign[av0] == 0 ) {a_vert = av0;}
			else if(geofrm.bound_sign[av1] == 0) {a_vert = av1;}
			else if(geofrm.bound_sign[av2] == 0 ) {a_vert = av2;}
			else if(geofrm.bound_sign[av3] == 0) {a_vert = av3;}
			else {a_vert = av0;}//-1
		}
		else if(randomInt % 4 == 1) 
		{
			if(geofrm.bound_sign[av1] == 0 ) {a_vert = av1;}
			else if(geofrm.bound_sign[av2] == 0) {a_vert = av2;}
			else if(geofrm.bound_sign[av3] == 0 ) {a_vert = av3;}
			else if(geofrm.bound_sign[av0] == 0) {a_vert = av0;}
			else {a_vert = av1;}//-1
		}
		else if(randomInt % 4 == 2) 
		{
			if(geofrm.bound_sign[av2] == 0 ) {a_vert = av2;}
			else if(geofrm.bound_sign[av3] == 0) {a_vert = av3;}
			else if(geofrm.bound_sign[av0] == 0 ) {a_vert = av0;}
			else if(geofrm.bound_sign[av1] == 0) {a_vert = av1;}
			else {a_vert = av2;}//-1
		}
		else  
		{
			if(geofrm.bound_sign[av3] == 0 ) {a_vert = av3;}
			else if(geofrm.bound_sign[av0] == 0) {a_vert = av0;}
			else if(geofrm.bound_sign[av1] == 0 ) {a_vert = av1;}
			else if(geofrm.bound_sign[av2] == 0) {a_vert = av2;}
			else {a_vert = av3;}//-1
		}
		//		fprintf(output, "aspect_min, aspect_max, a_vert = %e %e %d\n", aspect_min, aspect_max, a_vert);
	}

	delete [] aspectArr;

//	fclose(output);
}

void Octree::optimization(geoframe& geofrm) {

	int i, j, v0, v1, v2, v3, v4, v5, v6, v7;
	int index, maxIndex, minAspIndx, hexaIndx, a_vert;
	//int my_list[200], nlist;
	int sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, temp_i, temp_j;
	float grad[3], step;
	float mv0[3], mv1[3], mv2[3], mv3[3], mv4[3], mv5[3], mv6[3], mv7[3], cond_max1;
	float mmv0[3], mmv1[3], mmv2[3], mmv3[3], mmv4[3];
	float jacobian, oddy, condition, jaco_min, jaco_max, oddy_min, oddy_max, cond_min, cond_max;
	float *jacobianArr, *conditionArr, *oddyArr;

	//	FILE *output;

	maxIndex = 20;

	conditionArr=0;
	conditionArr = new float[geofrm.numverts];
	jacobianArr=0;
	jacobianArr = new float[geofrm.numverts];
	oddyArr=0;
	oddyArr = new float[geofrm.numverts];
	for(i = 0; i < geofrm.numverts; i++)	{
		conditionArr[i] = 1.0f;
		jacobianArr[i] = 1.0f;
		oddyArr[i] = 0.0f;
	}

	//	output = fopen("qimprove_temp", "w");

	//	fprintf(output, "optimization\n");
	
	//	fprintf(output, "%d %d\n", geofrm.numverts, geofrm.numquads/6);

	jaco_min = 1.0f;		jaco_max = -1.0f;
	oddy_min = 1000.0f;		oddy_max = -1000.0f;
	cond_min = 1000.0f;		cond_max = -1000.0f;
	float temp = 0.0f;
	float temp0 = 0.0f;
	float temp1 = 0.0f;
	int num = 0;
	for(i = 0; i < geofrm.numquads/6; i++) {
		v0 = geofrm.quads[6*i][0];		v1 = geofrm.quads[6*i][1];
		v2 = geofrm.quads[6*i][2];		v3 = geofrm.quads[6*i][3];
		v4 = geofrm.quads[6*i+1][1];	v5 = geofrm.quads[6*i+1][0];
		v6 = geofrm.quads[6*i+1][3];	v7 = geofrm.quads[6*i+1][2];

		for(j = 0; j < 3; j++) {
			mv0[j] = geofrm.verts[v0][j];	mv1[j] = geofrm.verts[v1][j];
			mv2[j] = geofrm.verts[v2][j];	mv3[j] = geofrm.verts[v3][j];	
			mv4[j] = geofrm.verts[v4][j];	mv5[j] = geofrm.verts[v5][j];
			mv6[j] = geofrm.verts[v6][j];	mv7[j] = geofrm.verts[v7][j];	
		}

		for(j = 0; j < 8; j++) {
			if(j == 0) testHexa(mv0, mv3, mv1, mv4, jacobian, condition, oddy);
			else if(j == 1) testHexa(mv1, mv0, mv2, mv5, jacobian, condition, oddy);
			else if(j == 2) testHexa(mv2, mv1, mv3, mv6, jacobian, condition, oddy);
			else if(j == 3) testHexa(mv3, mv2, mv0, mv7, jacobian, condition, oddy);
			else if(j == 4) testHexa(mv4, mv5, mv7, mv0, jacobian, condition, oddy);
			else if(j == 5) testHexa(mv5, mv6, mv4, mv1, jacobian, condition, oddy);
			else if(j == 6) testHexa(mv6, mv7, mv5, mv2, jacobian, condition, oddy);
			else testHexa(mv7, mv4, mv6, mv3, jacobian, condition, oddy);

			if(jacobian < jaco_min) jaco_min = jacobian;
			if(jacobian > jaco_max) jaco_max = jacobian;
			if(oddy < oddy_min) oddy_min = oddy;
			if(oddy > oddy_max) oddy_max = oddy;
			if(condition < cond_min) {
				cond_min = condition;	
				if(condition < 0.0f) {minAspIndx = j; hexaIndx = i;}
			}
			if(condition > 0.0f && condition > cond_max) {
			//if(condition > 0.0f && condition > 13.0f) {
				cond_max = condition;
				if(cond_min > 0.0f) {minAspIndx = j; hexaIndx = i;}
				//if(cond_min > 0.0f) {minAspIndx = 2; hexaIndx = i;}
			}
			temp_i = 6*i;	temp_j = j;
			if(j > 4) {temp_i = 6*i+1; temp_j = j -4;}
			if(conditionArr[geofrm.quads[temp_i][temp_j]] > 0.0f && condition > conditionArr[geofrm.quads[temp_i][temp_j]]) 
				conditionArr[geofrm.quads[temp_i][temp_j]] = condition;
			if(condition < 0.0f) conditionArr[geofrm.quads[temp_i][temp_j]] = condition;
			if(jacobian < jacobianArr[geofrm.quads[temp_i][temp_j]]) jacobianArr[geofrm.quads[temp_i][temp_j]] = jacobian;
			if(oddy > oddyArr[geofrm.quads[temp_i][temp_j]]) oddyArr[geofrm.quads[temp_i][temp_j]] = oddy;

			//int a = bound_sign[v0]+bound_sign[v1]+bound_sign[v2]+bound_sign[v3]+
			//	bound_sign[v4]+bound_sign[v5]+bound_sign[v6]+bound_sign[v7];
			//if(condition > 13.0f && a < 9) 
			//if(condition < 0.0f) 
			//	printf("bound_sign=%5.1f %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
			//	condition, v0, v1, v2, v3, v4, v5, v6, v7,
			//	bound_sign[v0], bound_sign[v1], bound_sign[v2], bound_sign[v3],
			//	bound_sign[v4], bound_sign[v5], bound_sign[v6], bound_sign[v7], a, i, j);

			temp += jacobian;	temp0 += condition;		temp1 += oddy;
			num++;
		}
	}

	sum0 = 0;	sum1 = 0;	sum2 = 0;	sum3 = 0;
	sum4 = 0;	sum5 = 0;	sum6 = 0;	sum7 = 0;	sum8 = 0;
	for(i = 0; i < geofrm.numverts; i++) {
		if(conditionArr[i] < 0.0f) sum0++;
		else if(conditionArr[i] < 25.0f) sum1++;
		else if(conditionArr[i] < 50.0f) sum2++;
		else if(conditionArr[i] < 75.0f) sum3++;
		else if(conditionArr[i] < 100.0f) sum4++;
		else if(conditionArr[i] < 125.0f) sum5++;
		else if(conditionArr[i] < 150.0f) sum6++;
		else if(conditionArr[i] < 200.0f) sum7++;
		else sum8++;
		//temp += conditionArr[i];
		//temp0 += jacobianArr[i];
		//temp1 += oddyArr[i];
	}
	//	fprintf(output, "statistics = %d %d %d %d %d %d %d %d %d %d\n", sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8,
	//		sum0+sum1+sum2+sum3+sum4+sum5+sum6+sum7+sum8);

	if(minAspIndx < 4)
		a_vert = geofrm.quads[6*hexaIndx][minAspIndx];
	else
		a_vert = geofrm.quads[6*hexaIndx+1][minAspIndx-4];
	
	//	fprintf(output, "num temp = %d %f\n", num, temp);
	//fprintf(output, "jacob oddy cond = (%f, %f, %f) (%f, %f, %f) (%f, %f, %f) %d\n", jaco_max, temp/num, jaco_min, 
	//cond_min, temp0/num, cond_max, oddy_min, temp1/num, oddy_max, a_vert);

	cond_max1 = cond_max + 1.0f;
	for(index = 0; index < maxIndex; index++)
	{
		printf("Loop %d\n", index);
/*
		for (i = 0; i < 200; i++) {
			my_list[i] = -1;
		}
		nlist = 0;
		for(i = 0; i < geofrm.numquads/6; i++) {
			v0 = geofrm.quads[6*i][0];		v1 = geofrm.quads[6*i][1];
			v2 = geofrm.quads[6*i][2];		v3 = geofrm.quads[6*i][3];
			v4 = geofrm.quads[6*i+1][1];	v5 = geofrm.quads[6*i+1][0];
			v6 = geofrm.quads[6*i+1][3];	v7 = geofrm.quads[6*i+1][2];

			vv0 = -1;	vv1 = -1;	vv2 = -1;
			if(v0 == a_vert) {vv0 = v1, vv1 = v3; vv2 = v4;}
			if(v1 == a_vert) {vv0 = v0, vv1 = v2; vv2 = v5;}
			if(v2 == a_vert) {vv0 = v1, vv1 = v3; vv2 = v6;}
			if(v3 == a_vert) {vv0 = v0, vv1 = v2; vv2 = v7;}
			if(v4 == a_vert) {vv0 = v0, vv1 = v5; vv2 = v7;}
			if(v5 == a_vert) {vv0 = v1, vv1 = v4; vv2 = v6;}
			if(v6 == a_vert) {vv0 = v2, vv1 = v5; vv2 = v7;}
			if(v7 == a_vert) {vv0 = v3, vv1 = v4; vv2 = v6;}

			if(vv0 != -1 && vv1 != -1 && vv2 != -1) {
				my_bool = 0;
				for(j = 0; j < nlist; j++) {
					if(vv0 == my_list[j]) {my_bool++; break;}
				}
				sign_bool = (bound_sign[a_vert] == 0) ||
							(bound_sign[a_vert] == 1 && bound_sign[vv0] == 1);
				if(my_bool == 0 && sign_bool == 1) {my_list[nlist] = vv0; nlist++;}

				my_bool = 0;
				for(j = 0; j < nlist; j++) {
					if(vv1 == my_list[j]) {my_bool++; break;}
				}
				sign_bool = (bound_sign[a_vert] == 0) ||
							(bound_sign[a_vert] == 1 && bound_sign[vv1] == 1);
				if(my_bool == 0 && sign_bool == 1) {my_list[nlist] = vv1; nlist++;}

				my_bool = 0;
				for(j = 0; j < nlist; j++) {
					if(vv2 == my_list[j]) {my_bool++; break;}
				}
				sign_bool = (bound_sign[a_vert] == 0) ||
							(bound_sign[a_vert] == 1 && bound_sign[vv2] == 1);
				if(my_bool == 0 && sign_bool == 1) {my_list[nlist] = vv2; nlist++;}
			}

		}

		printf("nlist = %d\n", nlist);

		i = a_vert;
		printf("average %d %f %f %f %d\n", i, vtx[i][0], vtx[i][1], vtx[i][2], bound_sign[i]);
		sum_0 = 0.0f;	sum_1 = 0.0f;	sum_2 = 0.0f;
		for(j = 0; j < nlist; j++) {
			sum_0 += vtx[my_list[j]][0];
			sum_1 += vtx[my_list[j]][1];
			sum_2 += vtx[my_list[j]][2];
		}
		//vtx[i][0] = sum_0/nlist;
		//vtx[i][1] = sum_1/nlist;
		//vtx[i][2] = sum_2/nlist;
		vtx[i][0] = (3.0f*vtx[i][0] - sum_0/nlist)/2.0f;
		vtx[i][1] = (3.0f*vtx[i][1] - sum_1/nlist)/2.0f;
		vtx[i][2] = (3.0f*vtx[i][2] - sum_2/nlist)/2.0f;
		printf("average %d %f %f %f\n", i, vtx[i][0], vtx[i][1], vtx[i][2]);
*/
//		fprintf(output, "average %d %f %f %f %d\n", a_vert, geofrm.verts[a_vert][0], geofrm.verts[a_vert][1], geofrm.verts[a_vert][2], geofrm.bound_sign[a_vert]);
		v0 = geofrm.quads[6*hexaIndx][0];		v1 = geofrm.quads[6*hexaIndx][1];
		v2 = geofrm.quads[6*hexaIndx][2];		v3 = geofrm.quads[6*hexaIndx][3];
		v4 = geofrm.quads[6*hexaIndx+1][1];		v5 = geofrm.quads[6*hexaIndx+1][0];
		v6 = geofrm.quads[6*hexaIndx+1][3];		v7 = geofrm.quads[6*hexaIndx+1][2];

		for(j = 0; j < 3; j++) {
			if(minAspIndx == 0) {
				mmv0[j] = geofrm.verts[v0][j];	mmv1[j] = geofrm.verts[v3][j];	mmv2[j] = geofrm.verts[v1][j];	mmv3[j] = geofrm.verts[v4][j];}
			else if(minAspIndx == 1){
				mmv0[j] = geofrm.verts[v1][j];	mmv1[j] = geofrm.verts[v0][j];	mmv2[j] = geofrm.verts[v2][j];	mmv3[j] = geofrm.verts[v5][j];}
			else if(minAspIndx == 2){
				mmv0[j] = geofrm.verts[v2][j];	mmv1[j] = geofrm.verts[v1][j];	mmv2[j] = geofrm.verts[v3][j];	mmv3[j] = geofrm.verts[v6][j];}
			else if(minAspIndx == 3){
				mmv0[j] = geofrm.verts[v3][j];	mmv1[j] = geofrm.verts[v2][j];	mmv2[j] = geofrm.verts[v0][j];	mmv3[j] = geofrm.verts[v7][j];}
			else if(minAspIndx == 4){
				mmv0[j] = geofrm.verts[v1][j];	mmv1[j] = geofrm.verts[v0][j];	mmv2[j] = geofrm.verts[v2][j];	mmv3[j] = geofrm.verts[v5][j];}
			else if(minAspIndx == 5){
				mmv0[j] = geofrm.verts[v5][j];	mmv1[j] = geofrm.verts[v6][j];	mmv2[j] = geofrm.verts[v4][j];	mmv3[j] = geofrm.verts[v1][j];}
			else if(minAspIndx == 6){
				mmv0[j] = geofrm.verts[v6][j];	mmv1[j] = geofrm.verts[v7][j];	mmv2[j] = geofrm.verts[v5][j];	mmv3[j] = geofrm.verts[v2][j];}
			else {
				mmv0[j] = geofrm.verts[v7][j];	mmv1[j] = geofrm.verts[v4][j];	mmv2[j] = geofrm.verts[v6][j];	mmv3[j] = geofrm.verts[v3][j];}
		}

		getGrad(mmv0, mmv1, mmv2, mmv3, condition, grad);

		condition = cond_max + 1.0f;
		step = 0.5f;
		while(condition - cond_max > 1.0e-6) {
			for(j = 0; j < 3; j++) {
				mmv4[j] = mmv0[j] - step*grad[j];
			}
			getGrad(mmv4, mmv1, mmv2, mmv3, condition, grad);
			if(condition > cond_max) step /= 2.0f;
		}

		geofrm.verts[a_vert][0] = mmv4[0];
		geofrm.verts[a_vert][1] = mmv4[1];
		geofrm.verts[a_vert][2] = mmv4[2];

		//		fprintf(output, "average %d %f %f %f %f\n", a_vert, geofrm.verts[a_vert][0], geofrm.verts[a_vert][1], geofrm.verts[a_vert][2], step);

		for(i = 0; i < geofrm.numverts ; i++) {conditionArr[i] = 1.0f;	jacobianArr[i] = 1.0f;	oddyArr[i] = 0.0f;}
		jaco_min = 1.0f;		jaco_max = -1.0f;
		oddy_min = 1000.0f;		oddy_max = -1000.0f;
		cond_min = 1000.0f;		cond_max = -1000.0f;
		temp = 0.0f;	temp0 = 0.0f;	temp1 = 0.0f;	num = 0;
		for(i = 0; i < geofrm.numquads/6; i++) {

			v0 = geofrm.quads[6*i][0];		v1 = geofrm.quads[6*i][1];
			v2 = geofrm.quads[6*i][2];		v3 = geofrm.quads[6*i][3];
			v4 = geofrm.quads[6*i+1][1];	v5 = geofrm.quads[6*i+1][0];
			v6 = geofrm.quads[6*i+1][3];	v7 = geofrm.quads[6*i+1][2];

			for(j = 0; j < 3; j++) {
				mv0[j] = geofrm.verts[v0][j];	mv1[j] = geofrm.verts[v1][j];
				mv2[j] = geofrm.verts[v2][j];	mv3[j] = geofrm.verts[v3][j];	
				mv4[j] = geofrm.verts[v4][j];	mv5[j] = geofrm.verts[v5][j];
				mv6[j] = geofrm.verts[v6][j];	mv7[j] = geofrm.verts[v7][j];	
			}

			for(j = 0; j < 8; j++) {
				if(j == 0) testHexa(mv0, mv3, mv1, mv4, jacobian, condition, oddy);
				else if(j == 1) testHexa(mv1, mv0, mv2, mv5, jacobian, condition, oddy);
				else if(j == 2) testHexa(mv2, mv1, mv3, mv6, jacobian, condition, oddy);
				else if(j == 3) testHexa(mv3, mv2, mv0, mv7, jacobian, condition, oddy);
				else if(j == 4) testHexa(mv4, mv5, mv7, mv0, jacobian, condition, oddy);
				else if(j == 5) testHexa(mv5, mv6, mv4, mv1, jacobian, condition, oddy);
				else if(j == 6) testHexa(mv6, mv7, mv5, mv2, jacobian, condition, oddy);
				else testHexa(mv7, mv4, mv6, mv3, jacobian, condition, oddy);

				if(jacobian < jaco_min) jaco_min = jacobian;
				if(jacobian > jaco_max) jaco_max = jacobian;
				if(oddy < oddy_min) oddy_min = oddy;
				if(oddy > oddy_max) oddy_max = oddy;
				//if(condition < cond_min) cond_min = condition;
				//if(condition > cond_max) {
				//	cond_max = condition;
				//	minAspIndx = j; hexaIndx = i;
				//}
				if(condition < cond_min) {
					cond_min = condition;	
					if(condition < 0.0f) {minAspIndx = j; hexaIndx = i;}
				}
				if(condition > 0.0f && condition > cond_max) {
				//if(condition > 0.0f && condition > 50.0f) {
					cond_max = condition;
					if(cond_min > 0.0f) {minAspIndx = j; hexaIndx = i;}
					//if(cond_min > 0.0f) {minAspIndx = 4; hexaIndx = i;}
				}
				//if(condition > conditionArr[vertexArr[i][j]]) conditionArr[vertexArr[i][j]] = condition;
				temp_i = 6*i;	temp_j = j;
				if(j > 4) {temp_i = 6*i+1; temp_j = j -4;}
				if(conditionArr[geofrm.quads[temp_i][temp_j]] > 0.0f && condition > conditionArr[geofrm.quads[temp_i][temp_j]]) 
					conditionArr[geofrm.quads[temp_i][temp_j]] = condition;
				if(condition < 0.0f) conditionArr[geofrm.quads[temp_i][temp_j]] = condition;
				if(jacobian < jacobianArr[geofrm.quads[temp_i][temp_j]]) jacobianArr[geofrm.quads[temp_i][temp_j]] = jacobian;
				if(oddy > oddyArr[geofrm.quads[temp_i][temp_j]]) oddyArr[geofrm.quads[temp_i][temp_j]] = oddy;

				temp += jacobian;	temp0 += condition;		temp1 += oddy;	num++;
			}

		}
		
		if(minAspIndx < 4)
			a_vert = geofrm.quads[6*hexaIndx][minAspIndx];
		else
			a_vert = geofrm.quads[6*hexaIndx+1][minAspIndx-4];

		//		fprintf(output, "min, max, a_vert = %f %f %f %f %f %f %d %d\n", jaco_min, jaco_max, 
		//				oddy_min, oddy_max, cond_min, cond_max, a_vert, geofrm.bound_sign[a_vert]);

		if(fabs(cond_max1 - cond_max) < 1.0e-4 && cond_min > 0.0f && cond_max < 10.0f) break;
		//if(cond_max < 216.94f) break;
		cond_max1 = cond_max;
	}

	sum0 = 0;	sum1 = 0;	sum2 = 0;	sum3 = 0;
	sum4 = 0;	sum5 = 0;	sum6 = 0;	sum7 = 0;	sum8 = 0;
	for(i = 0; i < geofrm.numverts; i++) {
		if(conditionArr[i] < 0.0f) sum0++;
		else if(conditionArr[i] < 25.0f) sum1++;
		else if(conditionArr[i] < 50.0f) sum2++;
		else if(conditionArr[i] < 75.0f) sum3++;
		else if(conditionArr[i] < 100.0f) sum4++;
		else if(conditionArr[i] < 125.0f) sum5++;
		else if(conditionArr[i] < 150.0f) sum6++;
		else if(conditionArr[i] < 200.0f) sum7++;
		else sum8++;
		//temp += conditionArr[i];
		//temp0 += jacobianArr[i];
		//temp1 += oddyArr[i];
	}
	//	fprintf(output, "statistics = %d %d %d %d %d %d %d %d %d %d\n", sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8,
	//		sum0+sum1+sum2+sum3+sum4+sum5+sum6+sum7+sum8);

	//	fprintf(output, "jacob oddy cond = (%f, %f, %f) (%f, %f, %f) (%f, %f, %f)\n", jaco_max, temp/num, jaco_min, 
	//			cond_min, temp0/num, cond_max, oddy_min, temp1/num, oddy_max);

	delete [] conditionArr;
	delete [] jacobianArr;
	delete [] oddyArr;

	//	fclose(output);
}

// calculate Jacobian, Condition_no and Oddy metric for a hexa
void Octree::testHexa(float v0[3], float v1[3], float v2[3], float v3[3], float& jacob, float& condition_no, float& oddy) {

	float v01[3], v02[3], v03[3], g11, g12, g13, g22, g23, g33, temp;

	int i;
	for(i = 0; i < 3; i++) {
		v01[i] = v1[i] - v0[i];
		v02[i] = v2[i] - v0[i];
		v03[i] = v3[i] - v0[i];
	}

	g11 = v01[0]*v01[0] + v01[1]*v01[1] + v01[2]*v01[2]; 
	g22 = v02[0]*v02[0] + v02[1]*v02[1] + v02[2]*v02[2];
	g33 = v03[0]*v03[0] + v03[1]*v03[1] + v03[2]*v03[2];

	for(i = 0; i < 3; i++) {
		v01[i] /= (float) sqrt(g11);
		v02[i] /= (float) sqrt(g22);
		v03[i] /= (float) sqrt(g33);
	}

	g11 = 1.0f;		g22 = 1.0f;		g33 = 1.0f;
	g12 = v01[0]*v02[0] + v01[1]*v02[1] + v01[2]*v02[2]; 
	g13 = v01[0]*v03[0] + v01[1]*v03[1] + v01[2]*v03[2]; 
	g23 = v02[0]*v03[0] + v02[1]*v03[1] + v02[2]*v03[2]; 

	jacob = g11*g22*g33 + 2.0f*g12*g13*g23 - g22*g13*g13 - g11*g23*g23 - g12*g12*g33;
	if(jacob < 0.0f) jacob = -(float) sqrt(fabs(jacob));
	else jacob = (float) sqrt(jacob);
	condition_no = (g11 + g22 + g33) / (3.0f*jacob);
	temp = g11*g11 + g22*g22 + g33*g33 + 2.0f*(g12*g12 + g13*g13 + g23*g23);
	oddy = (temp - (g11+g22+g33)*(g11+g22+g33)/3.0f)/(float)pow(fabs(jacob), 4.0f/3.0f);

}

// calculate the gradient of 1/Condition_no function for a hexa
void Octree::getGrad(float v0[3], float v1[3], float v2[3], float v3[3], float& condition_no, float* grad) {

	float v01[3], v02[3], v03[3], g11, g12, g13, g22, g23, g33;
	float jacob, t0, t1;

	int i;
	for(i = 0; i < 3; i++) {
		v01[i] = v1[i] - v0[i];
		v02[i] = v2[i] - v0[i];
		v03[i] = v3[i] - v0[i];
	}

	g11 = v01[0]*v01[0] + v01[1]*v01[1] + v01[2]*v01[2]; 
	g22 = v02[0]*v02[0] + v02[1]*v02[1] + v02[2]*v02[2];
	g33 = v03[0]*v03[0] + v03[1]*v03[1] + v03[2]*v03[2];

	for(i = 0; i < 3; i++) {
		v01[i] /= (float) sqrt(g11);
		v02[i] /= (float) sqrt(g22);
		v03[i] /= (float) sqrt(g33);
	}

	g11 = 1.0f;		g22 = 1.0f;		g33 = 1.0f;
	g12 = v01[0]*v02[0] + v01[1]*v02[1] + v01[2]*v02[2]; 
	g13 = v01[0]*v03[0] + v01[1]*v03[1] + v01[2]*v03[2]; 
	g23 = v02[0]*v03[0] + v02[1]*v03[1] + v02[2]*v03[2]; 

	//jacob = (float) sqrt(g11*g22*g33 + 2.0f*g12*g13*g23 - g22*g13*g13 - g11*g23*g23 - g12*g12*g33);
	jacob = g11*g22*g33 + 2.0f*g12*g13*g23 - g22*g13*g13 - g11*g23*g23 - g12*g12*g33;
	if(jacob < 0.0f) jacob = -(float) sqrt(fabs(jacob));
	else jacob = (float) sqrt(jacob);
	condition_no = (g11 + g22 + g33) / (3.0f*jacob);

	t0 = 2.0f*(v1[0] + v2[0] + v3[0] - 3.0f*v0[0]);
	t1 = (v1[1]-v0[1])*(v3[3]-v2[3]) + (v2[1]-v0[1])*(v1[3]-v3[3]) + (v3[1]-v0[1])*(v2[3]-v1[3]);
	grad[0] = t0/jacob - (g11+g22+g33)*t1/(jacob*jacob);
	//grad[0] = t1/(g11+g22+g33) - jacob*t0/((g11+g22+g33)*(g11+g22+g33));

	t0 = 2.0f*(v1[1] + v2[1] + v3[1] - 3.0f*v0[1]);
	t1 = (v1[0]-v0[0])*(v2[3]-v3[3]) + (v2[0]-v0[0])*(v3[3]-v1[3]) + (v3[0]-v0[0])*(v1[3]-v2[3]);
	grad[1] = t0/jacob - (g11+g22+g33)*t1/(jacob*jacob);
	//grad[1] = t1/(g11+g22+g33) - jacob*t0/((g11+g22+g33)*(g11+g22+g33));

	t0 = 2.0f*(v1[2] + v2[2] + v3[2] - 3.0f*v0[2]);
	t1 = (v1[0]-v0[0])*(v3[2]-v2[2]) + (v2[0]-v0[0])*(v1[2]-v3[2]) + (v3[0]-v0[0])*(v2[2]-v1[2]);
	grad[2] = t0/jacob - (g11+g22+g33)*t1/(jacob*jacob);
	//grad[2] = t1/(g11+g22+g33) - jacob*t0/((g11+g22+g33)*(g11+g22+g33));

	t0 = grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2];
	for(i = 0; i < 3; i++) grad[i] /= (float) sqrt(t0);
}


void Octree::polygonize_quad(geoframe& geofrm, float err_tol) {

	int x, y, z, valid_leaf, level;
	int vtx_num, intersect_id, i, j, k, oc_id[4], flag_method;
	unsigned int vtx[4];
	float val[8];

	in_out = 0;
	for (k = 0; k < octcell_num; k++) {vtx_idx_arr[k] = -1;}

	flag_method = 5;	// 1, 2, 3, 4, 5
	if(flag_method == 5) assign_refine_sign_quad(geofrm, err_tol);

	for (i = 0; i < leaf_num; i++ ) {
		valid_leaf = cut_array[i] ;
		level = get_level(valid_leaf) ;
		octcell2xyz(valid_leaf, x, y, z, level);
		getCellValues(valid_leaf, level, val);

		for (j = 0 ; j < 12 ; j++ ) {
			if (is_eflag_on(x,y,z,level,j)) continue;
			intersect_id = is_intersect(val, j);

			if (intersect_id == 1 || intersect_id == -1) {
				if (is_min_edge(valid_leaf, j, vtx, vtx_num, intersect_id, geofrm)) {
					eflag_on(x, y, z, level, j);
					find_oc_id(x, y, z, level, j, intersect_id, oc_id);
					//find_vtx_new(geofrm, x, y, z, level, j, intersect_id, vtx_new);
					quad_adaptive(geofrm, oc_id, err_tol, vtx, flag_method);
				}
			}
		}	// end j
	}		// end i
}

void Octree::get_VtxNorm(float* vtx, float* norm) 
{
	// finest resolution
	int oc_id, x, y, z;
	float val[8], tx, ty, tz;

	x = (int) vtx[0];	tx = vtx[0] - (float) x;
	y = (int) vtx[1];	ty = vtx[1] - (float) y;
	z = (int) vtx[2];	tz = vtx[2] - (float) z;

	oc_id = xyz2octcell(x, y, z, oct_depth);
	getCellValues(oc_id, oct_depth, val);

	norm[0] = (1-ty)*(1-tz)*(val[1] - val[0]) + (1-ty)*tz*(val[2] - val[3])
				+ ty*(1-tz)*(val[5] - val[4]) + ty*tz*(val[6] - val[7]);
	norm[1] = (1-tx)*(1-tz)*(val[4] - val[0]) + (1-tx)*tz*(val[7] - val[3])
				+ tx*(1-tz)*(val[5] - val[1]) + tx*tz*(val[6] - val[2]);
	norm[2] = (1-tx)*(1-ty)*(val[3] - val[0]) + (1-tx)*ty*(val[7] - val[4])
				+ tx*(1-ty)*(val[2] - val[1]) + tx*ty*(val[6] - val[5]);
}


void Octree::quad_adaptive(geoframe& geofrm, int* oc_id, float err_tol, unsigned int* vtx, int flag_method) {

	switch (flag_method) {
	case 1:							
		quad_adaptive_method1(geofrm, oc_id, err_tol, vtx);
		break;

	case 2:						
		quad_adaptive_method2(geofrm, oc_id, err_tol, vtx);
		break;

	case 3:						
		quad_adaptive_method3(geofrm, oc_id, err_tol, vtx, flag_method);
		break;

	case 4:						
		quad_adaptive_method3(geofrm, oc_id, err_tol, vtx, flag_method);
		break;

	case 5:						
		quad_adaptive_method5(geofrm, oc_id, err_tol, vtx);
		break;
	}

}

void Octree::quad_adaptive_method1(geoframe& geofrm, int* oc_id, float err_tol, unsigned int* vtx) {

	int i, j, vtx_num, xx, yy, zz, level, cell_size, vtx_new_num;
	unsigned int vtx_new[12];
	float x, y, z;
	vtx_num = 4;
	
	if (get_err_grad(oc_id[0]) > err_tol || get_err_grad(oc_id[1]) > err_tol ||
		get_err_grad(oc_id[2]) > err_tol || get_err_grad(oc_id[3]) > err_tol) {

		vtx_new_num = 4;
		geofrm.AddVert_adaptive(vtx, vtx_new);
		for(i = 0; i < vtx_new_num; i++) {
			for(j = 0; j < 4; j++) {
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
		}
		for(i = 0; i < 4; i++) get_vtx_new(geofrm, oc_id[i], vtx[i]);
		geofrm.AddQuad_adaptive(vtx, vtx_new, vtx_num);
	}
	else
		geofrm.AddQuad(vtx, vtx_num);
	
}

void Octree::quad_adaptive_method2(geoframe& geofrm, int* oc_id, float err_tol, unsigned int* vtx) {

	int num_id, i, j, vtx_num, xx, yy, zz, level, cell_size, vtx_new_num;
	unsigned int temp_vtx[4], vtx_new[12];
	float x, y, z;

	vtx_num = 4;

	num_id = 0;
	if(get_err_grad(oc_id[0]) > err_tol) num_id++;
	if(get_err_grad(oc_id[1]) > err_tol) num_id++;
	if(get_err_grad(oc_id[2]) > err_tol) num_id++;
	if(get_err_grad(oc_id[3]) > err_tol) num_id++;

	for(i = 0; i < 4; i++) {
		//if(get_err_grad(oc_id[i]) > err_tol) get_vtx_new(geofrm, oc_id[i], vtx[i]);
		get_vtx_new(geofrm, oc_id[i], vtx[i]);

		//if(geofrm.vtxnew_sign[vtx[i]] == 0) {
		//	get_vtx_new(geofrm, oc_id[i], vtx[i]);
		//	geofrm.AddVtxNew(vtx[i], 1);
		//}
	}

	for(i = 0; i < 4; i++)	temp_vtx[i] = vtx[i];
	if(num_id == 0)
		geofrm.AddQuad(vtx, vtx_num);
	else if(num_id == 1) {
		if(get_err_grad(oc_id[1]) > err_tol) {
			vtx[0] = temp_vtx[1];	vtx[1] = temp_vtx[2];
			vtx[2] = temp_vtx[3];	vtx[3] = temp_vtx[0];
		}
		else if(get_err_grad(oc_id[2]) > err_tol) {
			vtx[0] = temp_vtx[2];	vtx[1] = temp_vtx[3];
			vtx[2] = temp_vtx[0];	vtx[3] = temp_vtx[1];
		}
		else if(get_err_grad(oc_id[3]) > err_tol) {
			vtx[0] = temp_vtx[3];	vtx[1] = temp_vtx[0];
			vtx[2] = temp_vtx[1];	vtx[3] = temp_vtx[2];
		}
		geofrm.AddVert_adaptive_2_1(vtx, vtx_new);
		vtx_new_num = 6;
		for(i = 0; i < vtx_new_num; i++) {
			for(j = 0; j < 4; j++) {
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
			//if(geofrm.vtxnew_sign[vtx_new[i]] == 0) {
			//	get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
			//	geofrm.AddVtxNew(vtx_new[i], 1);
			//}
		}
		geofrm.AddQuad_adaptive_2_1(vtx, vtx_new, vtx_num);
	}
	else if(num_id == 2) {
		if( (get_err_grad(oc_id[0]) > err_tol && get_err_grad(oc_id[2]) > err_tol) ||
			(get_err_grad(oc_id[1]) > err_tol && get_err_grad(oc_id[3]) > err_tol) ) {
			geofrm.AddVert_adaptive_4(vtx, vtx_new);
			vtx_new_num = 12;
			for(i = 0; i < vtx_new_num; i++) {
				for(j = 0; j < 4; j++) {
					level = get_level(oc_id[j]) ;
					cell_size = (dim[0]-1)/(1<<level);
					octcell2xyz(oc_id[j], xx, yy, zz, level);
					x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
					y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
					z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
					if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
				}
				if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
				//if(geofrm.vtxnew_sign[vtx_new[i]] == 0) {
				//	get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
				//	geofrm.AddVtxNew(vtx_new[i], 1);
				//}
			}
			geofrm.AddQuad_adaptive_4(vtx, vtx_new, vtx_num);
		}
		else {
			if(get_err_grad(oc_id[1]) > err_tol && get_err_grad(oc_id[2]) > err_tol) {
				vtx[0] = temp_vtx[1];	vtx[1] = temp_vtx[2];
				vtx[2] = temp_vtx[3];	vtx[3] = temp_vtx[0];
			}
			else if(get_err_grad(oc_id[2]) > err_tol && get_err_grad(oc_id[3]) > err_tol) {
				vtx[0] = temp_vtx[2];	vtx[1] = temp_vtx[3];
				vtx[2] = temp_vtx[0];	vtx[3] = temp_vtx[1];
			}
			else if(get_err_grad(oc_id[3]) > err_tol && get_err_grad(oc_id[0]) > err_tol) {
				vtx[0] = temp_vtx[3];	vtx[1] = temp_vtx[0];
				vtx[2] = temp_vtx[1];	vtx[3] = temp_vtx[2];
			}
			geofrm.AddVert_adaptive_2_3(vtx, vtx_new);
			vtx_new_num = 8;
			for(i = 0; i < vtx_new_num; i++) {
				for(j = 0; j < 4; j++) {
					level = get_level(oc_id[j]) ;
					cell_size = (dim[0]-1)/(1<<level);
					octcell2xyz(oc_id[j], xx, yy, zz, level);
					x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
					y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
					z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
					if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
				}
				if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
				//if(geofrm.vtxnew_sign[vtx_new[i]] == 0) {
				//	get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
				//	geofrm.AddVtxNew(vtx_new[i], 1);
				//}
			}
			geofrm.AddQuad_adaptive_2_3(vtx, vtx_new, vtx_num);
		}
	}
	else {
		geofrm.AddVert_adaptive_4(vtx, vtx_new);
		vtx_new_num = 12;
		for(i = 0; i < vtx_new_num; i++) {
			for(j = 0; j < 4; j++) {
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
			//if(geofrm.vtxnew_sign[vtx_new[i]] == 0) {
			//	get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
			//	geofrm.AddVtxNew(vtx_new[i], 1);
			//}
		}
		geofrm.AddQuad_adaptive_4(vtx, vtx_new, vtx_num);
	}

}

void Octree::quad_adaptive_method3(geoframe& geofrm, int* oc_id, float err_tol, unsigned int* vtx, int flag_method) {

	int num_id, i, j, vtx_num, xx, yy, zz, level, cell_size, vtx_new_num;
	unsigned int temp_vtx[4], vtx_new[12];
	float x, y, z;

	vtx_num = 4;

	num_id = 0;
	if(get_err_grad(oc_id[0]) > err_tol) num_id++;
	if(get_err_grad(oc_id[1]) > err_tol) num_id++;
	if(get_err_grad(oc_id[2]) > err_tol) num_id++;
	if(get_err_grad(oc_id[3]) > err_tol) num_id++;

	for(i = 0; i < 4; i++) {
		//if(get_err_grad(oc_id[i]) > err_tol) get_vtx_new(geofrm, oc_id[i], vtx[i]);
		get_vtx_new(geofrm, oc_id[i], vtx[i]);
	}

	for(i = 0; i < 4; i++)	temp_vtx[i] = vtx[i];
	if(num_id == 0)
		geofrm.AddQuad(vtx, vtx_num);
	else if(num_id == 1) {
		if(get_err_grad(oc_id[1]) > err_tol) {
			vtx[0] = temp_vtx[1];	vtx[1] = temp_vtx[2];
			vtx[2] = temp_vtx[3];	vtx[3] = temp_vtx[0];
		}
		else if(get_err_grad(oc_id[2]) > err_tol) {
			vtx[0] = temp_vtx[2];	vtx[1] = temp_vtx[3];
			vtx[2] = temp_vtx[0];	vtx[3] = temp_vtx[1];
		}
		else if(get_err_grad(oc_id[3]) > err_tol) {
			vtx[0] = temp_vtx[3];	vtx[1] = temp_vtx[0];
			vtx[2] = temp_vtx[1];	vtx[3] = temp_vtx[2];
		}
		geofrm.AddVert_adaptive_3_1(vtx, vtx_new);
		vtx_new_num = 3;
		for(i = 0; i < vtx_new_num; i++) {
			for(j = 0; j < 4; j++) {
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
		}
		geofrm.AddQuad_adaptive_3_1(vtx, vtx_new, vtx_num);
	}
	else if(num_id == 2) {
		if( (get_err_grad(oc_id[0]) > err_tol && get_err_grad(oc_id[2]) > err_tol) ||
			(get_err_grad(oc_id[1]) > err_tol && get_err_grad(oc_id[3]) > err_tol) ) {
			if(get_err_grad(oc_id[1]) > err_tol && get_err_grad(oc_id[3]) > err_tol) {
				vtx[0] = temp_vtx[1];	vtx[1] = temp_vtx[2];
				vtx[2] = temp_vtx[3];	vtx[3] = temp_vtx[0];
			}
			if(flag_method == 3) {
				geofrm.AddVert_adaptive_3_2b(vtx, vtx_new);
				vtx_new_num = 5;
			}
			else {
				geofrm.AddVert_adaptive_4_2b(vtx, vtx_new);
				vtx_new_num = 8;
			}
			for(i = 0; i < vtx_new_num; i++) {
				for(j = 0; j < 4; j++) {
					level = get_level(oc_id[j]) ;
					cell_size = (dim[0]-1)/(1<<level);
					octcell2xyz(oc_id[j], xx, yy, zz, level);
					x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
					y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
					z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
					if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
				}
				if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
			}
			if(flag_method == 3) geofrm.AddQuad_adaptive_3_2b(vtx, vtx_new, vtx_num);
			else geofrm.AddQuad_adaptive_4_2b(vtx, vtx_new, vtx_num);
		}
		else {
			if(get_err_grad(oc_id[1]) > err_tol && get_err_grad(oc_id[2]) > err_tol) {
				vtx[0] = temp_vtx[1];	vtx[1] = temp_vtx[2];
				vtx[2] = temp_vtx[3];	vtx[3] = temp_vtx[0];
			}
			else if(get_err_grad(oc_id[2]) > err_tol && get_err_grad(oc_id[3]) > err_tol) {
				vtx[0] = temp_vtx[2];	vtx[1] = temp_vtx[3];
				vtx[2] = temp_vtx[0];	vtx[3] = temp_vtx[1];
			}
			else if(get_err_grad(oc_id[3]) > err_tol && get_err_grad(oc_id[0]) > err_tol) {
				vtx[0] = temp_vtx[3];	vtx[1] = temp_vtx[0];
				vtx[2] = temp_vtx[1];	vtx[3] = temp_vtx[2];
			}
			geofrm.AddVert_adaptive_3_2a(vtx, vtx_new);
			vtx_new_num = 8;
			for(i = 0; i < vtx_new_num; i++) {
				for(j = 0; j < 4; j++) {
					level = get_level(oc_id[j]) ;
					cell_size = (dim[0]-1)/(1<<level);
					octcell2xyz(oc_id[j], xx, yy, zz, level);
					x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
					y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
					z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
					if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
				}
				if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
			}
			geofrm.AddQuad_adaptive_3_2a(vtx, vtx_new, vtx_num);
		}
	}
	else if(num_id == 3) {
		if(get_err_grad(oc_id[0]) <= err_tol) {
			vtx[0] = temp_vtx[1];	vtx[1] = temp_vtx[2];
			vtx[2] = temp_vtx[3];	vtx[3] = temp_vtx[0];
		}
		else if(get_err_grad(oc_id[1]) <= err_tol) {
			vtx[0] = temp_vtx[2];	vtx[1] = temp_vtx[3];
			vtx[2] = temp_vtx[0];	vtx[3] = temp_vtx[1];
		}
		else if(get_err_grad(oc_id[2]) <= err_tol) {
			vtx[0] = temp_vtx[3];	vtx[1] = temp_vtx[0];
			vtx[2] = temp_vtx[1];	vtx[3] = temp_vtx[2];
		}
		geofrm.AddVert_adaptive_3_3(vtx, vtx_new);
		vtx_new_num = 10;
		for(i = 0; i < vtx_new_num; i++) {
			for(j = 0; j < 4; j++) {
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
		}
		geofrm.AddQuad_adaptive_3_3(vtx, vtx_new, vtx_num);
	}
	else {
		geofrm.AddVert_adaptive_4(vtx, vtx_new);
		vtx_new_num = 12;
		for(i = 0; i < vtx_new_num; i++) {
			for(j = 0; j < 4; j++) {
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
		}
		geofrm.AddQuad_adaptive_4(vtx, vtx_new, vtx_num);
	}

}

void Octree::quad_adaptive_method5(geoframe& geofrm, int* oc_id, float err_tol, unsigned int* vtx) {

	int num_id, i, j, vtx_num, xx, yy, zz, level, cell_size, vtx_new_num;
	unsigned int temp_vtx[4], vtx_new[12];
	float x, y, z;

	vtx_num = 4;

	num_id = 0;
	if(vtx_idx_arr_refine[oc_id[0]] == 1) num_id++;
	if(vtx_idx_arr_refine[oc_id[1]] == 1) num_id++;
	if(vtx_idx_arr_refine[oc_id[2]] == 1) num_id++;
	if(vtx_idx_arr_refine[oc_id[3]] == 1) num_id++;

	for(i = 0; i < 4; i++) {
		//if(get_err_grad(oc_id[i]) > err_tol) get_vtx_new(geofrm, oc_id[i], vtx[i]);
		get_vtx_new(geofrm, oc_id[i], vtx[i]);
	}

	for(i = 0; i < 4; i++)	temp_vtx[i] = vtx[i];
	if(num_id == 0)
		geofrm.AddQuad(vtx, vtx_num);
	else if(num_id == 1) {
		if(vtx_idx_arr_refine[oc_id[1]] == 1) {
			vtx[0] = temp_vtx[1];	vtx[1] = temp_vtx[2];
			vtx[2] = temp_vtx[3];	vtx[3] = temp_vtx[0];
		}
		else if(vtx_idx_arr_refine[oc_id[2]] == 1) {
			vtx[0] = temp_vtx[2];	vtx[1] = temp_vtx[3];
			vtx[2] = temp_vtx[0];	vtx[3] = temp_vtx[1];
		}
		else if(vtx_idx_arr_refine[oc_id[3]] == 1) {
			vtx[0] = temp_vtx[3];	vtx[1] = temp_vtx[0];
			vtx[2] = temp_vtx[1];	vtx[3] = temp_vtx[2];
		}
		geofrm.AddVert_adaptive_3_1(vtx, vtx_new);
		vtx_new_num = 3;
		for(i = 0; i < vtx_new_num; i++) {
			for(j = 0; j < 4; j++) {
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
		}
		geofrm.AddQuad_adaptive_3_1(vtx, vtx_new, vtx_num);
	}
	else if(num_id == 2) {
		if(vtx_idx_arr_refine[oc_id[1]] == 1 && vtx_idx_arr_refine[oc_id[2]] == 1) {
			vtx[0] = temp_vtx[1];	vtx[1] = temp_vtx[2];
			vtx[2] = temp_vtx[3];	vtx[3] = temp_vtx[0];
		}
		else if(vtx_idx_arr_refine[oc_id[2]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1) {
			vtx[0] = temp_vtx[2];	vtx[1] = temp_vtx[3];
			vtx[2] = temp_vtx[0];	vtx[3] = temp_vtx[1];
		}
		else if(vtx_idx_arr_refine[oc_id[3]] == 1 && vtx_idx_arr_refine[oc_id[0]] == 1) {
			vtx[0] = temp_vtx[3];	vtx[1] = temp_vtx[0];
			vtx[2] = temp_vtx[1];	vtx[3] = temp_vtx[2];
		}
		geofrm.AddVert_adaptive_3_2a(vtx, vtx_new);
		vtx_new_num = 8;
		for(i = 0; i < vtx_new_num; i++) {
			for(j = 0; j < 4; j++) {
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
		}
		geofrm.AddQuad_adaptive_3_2a(vtx, vtx_new, vtx_num);
	}
	else if(num_id == 4) {
		geofrm.AddVert_adaptive_4(vtx, vtx_new);
		vtx_new_num = 12;
		for(i = 0; i < vtx_new_num; i++) {
			for(j = 0; j < 4; j++) {
				level = get_level(oc_id[j]) ;
				cell_size = (dim[0]-1)/(1<<level);
				octcell2xyz(oc_id[j], xx, yy, zz, level);
				x = geofrm.verts[vtx_new[i]][0]/cell_size - (float)xx;
				y = geofrm.verts[vtx_new[i]][1]/cell_size - (float)yy;
				z = geofrm.verts[vtx_new[i]][2]/cell_size - (float)zz;
				if(x >= 0.0f && x <= 1.0f && y >= 0.0f && y <= 1.0f && z >= 0.0f && z <= 1.0f) break;
			}
			if(j < 4) get_vtx_new(geofrm, oc_id[j], vtx_new[i]);
		}
		geofrm.AddQuad_adaptive_4(vtx, vtx_new, vtx_num);
	}

}

void Octree::assign_refine_sign_quad(geoframe& geofrm, float err_tol) {

	int x, y, z, valid_leaf, level, num_id, refine_sign;
	int intersect_id, i, j, k, oc_id[4];
	float val[8];

	for (k = 0; k < octcell_num; k++) vtx_idx_arr_refine[k] = -1;

	for (i = 0; i < leaf_num; i++ ) {
		valid_leaf = cut_array[i] ;
		level = get_level(valid_leaf) ;
		octcell2xyz(valid_leaf, x, y, z, level);
		getCellValues(valid_leaf, level, val);

		for (j = 0 ; j < 12 ; j++ ) {
			if (is_eflag_on(x,y,z,level,j)) continue;
			intersect_id = is_intersect(val, j);

			if (intersect_id == 1 || intersect_id == -1) {
				eflag_on(x, y, z, level, j);
				find_oc_id(x, y, z, level, j, intersect_id, oc_id);
				
				num_id = 0;
				if(get_err_grad(oc_id[0]) > err_tol) num_id++;
				if(get_err_grad(oc_id[1]) > err_tol) num_id++;
				if(get_err_grad(oc_id[2]) > err_tol) num_id++;
				if(get_err_grad(oc_id[3]) > err_tol) num_id++;
				
				//if( num_id > 2 || 
				//	(num_id == 2 && get_err_grad(oc_id[0]) > err_tol && get_err_grad(oc_id[2]) > err_tol) ||
				//	(num_id == 2 && get_err_grad(oc_id[1]) > err_tol && get_err_grad(oc_id[3]) > err_tol) ) {
				if(num_id == 4) {
					vtx_idx_arr_refine[oc_id[0]] = 1;
					vtx_idx_arr_refine[oc_id[1]] = 1;
					vtx_idx_arr_refine[oc_id[2]] = 1;
					vtx_idx_arr_refine[oc_id[3]] = 1;
				}

			}
		}	// end j
	}		// end i

	eflag_clear();

	refine_sign = 1;

	while(refine_sign == 1) {

		refine_sign = 0;

		for (i = 0; i < leaf_num; i++ ) {
			valid_leaf = cut_array[i] ;
			level = get_level(valid_leaf) ;
			octcell2xyz(valid_leaf, x, y, z, level);
			getCellValues(valid_leaf, level, val);

			for (j = 0 ; j < 12 ; j++ ) {
				if (is_eflag_on(x,y,z,level,j)) continue;
				intersect_id = is_intersect(val, j);

				if (intersect_id == 1 || intersect_id == -1) {
					eflag_on(x, y, z, level, j);
					find_oc_id(x, y, z, level, j, intersect_id, oc_id);

					num_id = 0;
					if(vtx_idx_arr_refine[oc_id[0]] == 1) num_id++;
					if(vtx_idx_arr_refine[oc_id[1]] == 1) num_id++;
					if(vtx_idx_arr_refine[oc_id[2]] == 1) num_id++;
					if(vtx_idx_arr_refine[oc_id[3]] == 1) num_id++;

					if( num_id > 2 || 
						(num_id == 2 && vtx_idx_arr_refine[oc_id[0]] == 1 && vtx_idx_arr_refine[oc_id[2]] == 1) ||
						(num_id == 2 && vtx_idx_arr_refine[oc_id[1]] == 1 && vtx_idx_arr_refine[oc_id[3]] == 1) ) {
						if(vtx_idx_arr_refine[oc_id[0]] != 1) {vtx_idx_arr_refine[oc_id[0]] = 1; refine_sign = 1;}
						if(vtx_idx_arr_refine[oc_id[1]] != 1) {vtx_idx_arr_refine[oc_id[1]] = 1; refine_sign = 1;}
						if(vtx_idx_arr_refine[oc_id[2]] != 1) {vtx_idx_arr_refine[oc_id[2]] = 1; refine_sign = 1;}
						if(vtx_idx_arr_refine[oc_id[3]] != 1) {vtx_idx_arr_refine[oc_id[3]] = 1; refine_sign = 1;}
					}

				}
			}	// end j
		}		// end i

		eflag_clear();

	}			// end while

}

}
