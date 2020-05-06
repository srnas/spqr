#ifndef MCUTILS_H
#define MCUTILS_H
#include "mc_global.h"



extern double *idum;
extern long idum2;
extern long iy;
extern long iv[NTAB];
float nrran2();

static inline double dist_1d(double x1, double x2, int dir)
{
  double dist=x1-x2;
#ifdef PBC
  if(dist>box_l[dir]/2.0){
    dist=dist-box_l[dir];
  }
  else if(dist<-box_l[dir]/2.0){
    dist=dist+box_l[dir];
  }
#endif
  return dist;
};



static inline double vec_norm(double x, double y, double z){
  return sqrt(x*x+y*y+z*z);
}

static inline double fvec_norm(double *vec){
  return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}
static inline double fcalc_min_dist(double *a, double *b){
  double c[DIM];
  c[0]=a[0]-b[0];  c[1]=a[1]-b[1];  c[2]=a[2]-b[2];
  return fvec_norm(c);
}
static inline void calc_min_vec(double x1, double y1, double z1, double x2, double y2, double z2, double *rvec, double *r){
  rvec[0]=dist_1d(x1,x2,0);
  rvec[1]=dist_1d(y1,y2,1);
  rvec[2]=dist_1d(z1,z2,2);
  //*r=rvec[0]*rvec[0]+rvec[1]*rvec[1]+rvec[2]*rvec[2];
  *r=vec_norm(rvec[0], rvec[1], rvec[2]);
}

static inline double calc_min_dist(double x1, double y1, double z1, double x2, double y2, double z2){
  return vec_norm(dist_1d(x1,x2,0),dist_1d(y1,y2,1),dist_1d(z1,z2,2));
}

static inline double calc_min_dist_sq(double x1, double y1, double z1, double x2, double y2, double z2){
  return SQ(dist_1d(x1,x2,0))+SQ(dist_1d(y1,y2,1))+SQ(dist_1d(z1,z2,2));
}

static inline double rand_d(double a){
  return a*(double)nrran2();
}

static inline double boltzmann_dist(double average, double std_dev){
  /* we create a random gaussian number */
  double a,b;
  double r=2, vel;
  while(r>1){
    a=2.0*rand_d(1)-1.0;
    b=2.0*rand_d(1)-1.0;
    r=a*a+b*b;
  }
  vel=a*sqrt(-2*log(r)/r);
  vel=average+vel*sqrt(std_dev);
  
  return vel;
}

static inline void fold_coordinate(double *folded, double eff_box, double origin){
  /* the coordinate will be folded to the interval [0, box] */
  double coord=*folded;
  while(coord<origin)
    coord+=eff_box;
  while(coord > eff_box+origin)
    coord-=eff_box;
  
  *folded = coord;
}

static inline double max_2(double a, double b){
  if(a>b) return a;
  else return b;
}


static inline int rand_i(int n){
  return  (int)(n*nrran2());
}

static inline double dot_prod(double *a, double *b){
  int i;
  double prod=0.0;
  for(i=0;i<DIM;i++){
    prod+=a[i]*b[i];
  }
  return prod;
}

static inline void svec_prod( double *a, double *b, double *ret, double *norm){
  ret[0]=a[1]*b[2]-a[2]*b[1];
  ret[1]=a[2]*b[0]-a[0]*b[2];
  ret[2]=a[0]*b[1]-a[1]*b[0];
  
  *norm=sqrt(ret[0]*ret[0]+ret[1]*ret[1]+ret[2]*ret[2]);
}


static inline void vec_prod( double *a, double *b, double *ret){
  ret[0]=a[1]*b[2]-a[2]*b[1];
  ret[1]=a[2]*b[0]-a[0]*b[2];
  ret[2]=a[0]*b[1]-a[1]*b[0];
}

static inline void norex_vec_prod(double a0, double a1, double a2,double b0, double b1, double b2, double *ret){
  ret[0]=a1*b2-a2*b1;
  ret[1]=a2*b0-a0*b2;
  ret[2]=a0*b1-a1*b0;
  double norm = sqrt(dot_prod(ret, ret));
  ret[0]/=norm;
  ret[1]/=norm;
  ret[2]/=norm;
}


static inline double get_unf_coo_temp_x(int part){
#ifdef PBC 
  return mc_temp_x[part]+mc_temp_pbox[part][0]*box_l[0];
#else
  return mc_temp_x[part];
#endif
}
static inline double get_unf_coo_temp_y(int part){
#ifdef PBC 
  return mc_temp_y[part]+mc_temp_pbox[part][1]*box_l[1];
#else
  return mc_temp_y[part];
#endif
}
static inline double get_unf_coo_temp_z(int part){
#ifdef PBC 
  return mc_temp_z[part]+mc_temp_pbox[part][2]*box_l[2];
#else
  return mc_temp_z[part];
#endif
}

static inline double get_unf_coo_x(double *rx, int coord){
#ifdef PBC 
  return rx[coord]+mc_pbox[coord][0]*box_l[0];
#else
  return rx[coord];
#endif
}

static inline  double get_unf_coo_y(double *ry, int coord){
#ifdef PBC 
  return ry[coord]+mc_pbox[coord][1]*box_l[1];
#else
  return ry[coord];
#endif
}



static inline  double get_unf_coo_z(double *rz, int coord){
#ifdef PBC 
  return rz[coord]+mc_pbox[coord][2]*box_l[2];
#else
  return rz[coord];
#endif
}

static inline void proj_on_nt(double *d_vec, double *rx, double *ry, double *rz, int nt_c, double *r_vec){
  double pbase_x[DIM], pbase_y[DIM], pbase_z[DIM];//, pbasey[N_PARTS_PER_NT], pbasez[N_PARTS_PER_NT];
  //here we use that the first auxiliar particle is x, and the second, y
  int at_c=N_PARTS_PER_NT*nt_c;
  pbase_x[0]=dist_1d(rx[at_c+1],rx[at_c],0);
  pbase_x[1]=dist_1d(ry[at_c+1],ry[at_c],0);
  pbase_x[2]=dist_1d(rz[at_c+1],rz[at_c],0);
  pbase_y[0]=dist_1d(rx[at_c+2],rx[at_c],0);
  pbase_y[1]=dist_1d(ry[at_c+2],ry[at_c],0);
  pbase_y[2]=dist_1d(rz[at_c+2],rz[at_c],0);
  vec_prod(pbase_x, pbase_y, pbase_z);
  r_vec[0]=(d_vec[0]*pbase_x[0]+d_vec[1]*pbase_x[1]+d_vec[2]*pbase_x[2]);
  r_vec[1]=(d_vec[0]*pbase_y[0]+d_vec[1]*pbase_y[1]+d_vec[2]*pbase_y[2]);
  r_vec[2]=(d_vec[0]*pbase_z[0]+d_vec[1]*pbase_z[1]+d_vec[2]*pbase_z[2]);
}

static inline void proj_on_nt_inv(double *d_vec, double *rx, double *ry, double *rz, int nt_c, double *r_vec){
  double pbase_x[DIM], pbase_y[DIM], pbase_z[DIM];//, pbasey[N_PARTS_PER_NT], pbasez[N_PARTS_PER_NT];
  //here we use that the first auxiliar particle is x, and the second, y
  int at_c=N_PARTS_PER_NT*nt_c;
  pbase_x[0]=dist_1d(rx[at_c+1],rx[at_c],0);
  pbase_x[1]=dist_1d(ry[at_c+1],ry[at_c],0);
  pbase_x[2]=dist_1d(rz[at_c+1],rz[at_c],0);
  pbase_y[0]=dist_1d(rx[at_c+2],rx[at_c],0);
  pbase_y[1]=dist_1d(ry[at_c+2],ry[at_c],0);
  pbase_y[2]=dist_1d(rz[at_c+2],rz[at_c],0);
  vec_prod(pbase_x, pbase_y, pbase_z);
  d_vec[0]*=-1;d_vec[1]*=-1;d_vec[2]*=-1;
  r_vec[0]=(d_vec[0]*pbase_x[0]+d_vec[1]*pbase_x[1]+d_vec[2]*pbase_x[2]);
  r_vec[1]=(d_vec[0]*pbase_y[0]+d_vec[1]*pbase_y[1]+d_vec[2]*pbase_y[2]);
  r_vec[2]=(d_vec[0]*pbase_z[0]+d_vec[1]*pbase_z[1]+d_vec[2]*pbase_z[2]);
}

#endif

