#include "mc_base_pairing.h"

void MC_init_wc_arrays(int nt_n, double *rx, double *ry, double *rz){
  int i, nt_c, face;
  wc_face_ind_trial=(int **)malloc(sizeof(int *)*nt_n);
  wc_face_ene_trial=(double **)malloc(sizeof(double *)*nt_n);
  wc_is_paired_trial=(int **)malloc(sizeof(int *)*nt_n);
  wc_face_ind=(int **)malloc(sizeof(int *)*nt_n);
  wc_face_ene=(double **)malloc(sizeof(double *)*nt_n);
  wc_is_paired=(int **)malloc(sizeof(int *)*nt_n);
  
  for(i=0;i<nt_n;i++){
    wc_face_ind_trial[i]=(int *)malloc(sizeof(int)*WC_FACES);
    wc_face_ene_trial[i]=(double *)malloc(sizeof(double)*WC_FACES);
    wc_is_paired_trial[i]=(int *)malloc(sizeof(int)*WC_FACES);
    wc_face_ind[i]=(int *)malloc(sizeof(int)*WC_FACES);
    wc_face_ene[i]=(double *)malloc(sizeof(double)*WC_FACES);
    wc_is_paired[i]=(int *)malloc(sizeof(int)*WC_FACES);
  }
  
  for(face=0;face<WC_FACES;face++)
    for(nt_c=0;nt_c<nt_n;nt_c++){
      wc_face_ene[nt_c][face]=0;
      wc_face_ind[nt_c][face]=-1;
      wc_is_paired[nt_c][face]=-1;
      wc_face_ene_trial[nt_c][face]=0;
      wc_face_ind_trial[nt_c][face]=-1;
      wc_is_paired_trial[nt_c][face]=-1;
    }
  
  bp_face_ind_trial=(int *)malloc(sizeof(int )*nt_n);
  bp_face_ene_trial=(double *)malloc(sizeof(double )*nt_n);
  bp_is_paired_trial=(int *)malloc(sizeof(int )*nt_n);
  bp_face_ind=(int *)malloc(sizeof(int )*nt_n);
  bp_face_ene=(double *)malloc(sizeof(double )*nt_n);
  bp_is_paired=(int *)malloc(sizeof(int )*nt_n);
  for(nt_c=0;nt_c<nt_n;nt_c++){
    bp_face_ene[nt_c]=0;
    bp_face_ind[nt_c]=-1;
    bp_is_paired[nt_c]=-1;
    
    bp_face_ene_trial[nt_c]=0;
    bp_face_ind_trial[nt_c]=-1;
    bp_is_paired_trial[nt_c]=-1;


    
  }
  
  
  /** FACES **/
  wc_face_angle[0][0]=WCFACEANG_00_0u*M_PI/180;
  wc_face_angle[0][1]=WCFACEANG_00_1u*M_PI/180;
  wc_face_angle[0][2]=WCFACEANG_00_2u*M_PI/180;
  wc_face_angle[0][3]=WCFACEANG_00_3u*M_PI/180;
  wc_face_angle[1][0]=WCFACEANG_01_0u*M_PI/180;
  wc_face_angle[1][1]=WCFACEANG_01_1u*M_PI/180;
  wc_face_angle[1][2]=WCFACEANG_01_2u*M_PI/180;
  wc_face_angle[1][3]=WCFACEANG_01_3u*M_PI/180;
  wc_face_angle[2][0]=WCFACEANG_02_0u*M_PI/180;
  wc_face_angle[2][1]=WCFACEANG_02_1u*M_PI/180;
  wc_face_angle[2][2]=WCFACEANG_02_2u*M_PI/180;
  wc_face_angle[2][3]=WCFACEANG_02_3u*M_PI/180;
  wc_face_angle[3][0]=WCFACEANG_03_0u*M_PI/180;
  wc_face_angle[3][1]=WCFACEANG_03_1u*M_PI/180;
  wc_face_angle[3][2]=WCFACEANG_03_2u*M_PI/180;
  wc_face_angle[3][3]=WCFACEANG_03_3u*M_PI/180;
  
  wc_face_angle[4][0]=WCFACEANG_10_0u*M_PI/180;
  wc_face_angle[4][1]=WCFACEANG_10_1u*M_PI/180;
  wc_face_angle[4][2]=WCFACEANG_10_2u*M_PI/180;
  wc_face_angle[4][3]=WCFACEANG_10_3u*M_PI/180;
  wc_face_angle[5][0]=WCFACEANG_11_0u*M_PI/180;
  wc_face_angle[5][1]=WCFACEANG_11_1u*M_PI/180;
  wc_face_angle[5][2]=WCFACEANG_11_2u*M_PI/180;
  wc_face_angle[5][3]=WCFACEANG_11_3u*M_PI/180;
  wc_face_angle[6][0]=WCFACEANG_12_0u*M_PI/180;
  wc_face_angle[6][1]=WCFACEANG_12_1u*M_PI/180;
  wc_face_angle[6][2]=WCFACEANG_12_2u*M_PI/180;
  wc_face_angle[6][3]=WCFACEANG_12_3u*M_PI/180;
  wc_face_angle[7][0]=WCFACEANG_13_0u*M_PI/180;
  wc_face_angle[7][1]=WCFACEANG_13_1u*M_PI/180;
  wc_face_angle[7][2]=WCFACEANG_13_2u*M_PI/180;
  wc_face_angle[7][3]=WCFACEANG_13_3u*M_PI/180;
  
  wc_face_angle[8][0]=WCFACEANG_20_0u*M_PI/180;
  wc_face_angle[8][1]=WCFACEANG_20_1u*M_PI/180;
  wc_face_angle[8][2]=WCFACEANG_20_2u*M_PI/180;
  wc_face_angle[8][3]=WCFACEANG_20_3u*M_PI/180;
  wc_face_angle[9][0]=WCFACEANG_21_0u*M_PI/180;
  wc_face_angle[9][1]=WCFACEANG_21_1u*M_PI/180;
  wc_face_angle[9][2]=WCFACEANG_21_2u*M_PI/180;
  wc_face_angle[9][3]=WCFACEANG_21_3u*M_PI/180;
  wc_face_angle[10][0]=WCFACEANG_22_0u*M_PI/180;
  wc_face_angle[10][1]=WCFACEANG_22_1u*M_PI/180;
  wc_face_angle[10][2]=WCFACEANG_22_2u*M_PI/180;
  wc_face_angle[10][3]=WCFACEANG_22_3u*M_PI/180;
  wc_face_angle[11][0]=WCFACEANG_23_0u*M_PI/180;
  wc_face_angle[11][1]=WCFACEANG_23_1u*M_PI/180;
  wc_face_angle[11][2]=WCFACEANG_23_2u*M_PI/180;
  wc_face_angle[11][3]=WCFACEANG_23_3u*M_PI/180;
  
  wc_face_angle[12][0]=WCFACEANG_30_0u*M_PI/180;
  wc_face_angle[12][1]=WCFACEANG_30_1u*M_PI/180;
  wc_face_angle[12][2]=WCFACEANG_30_2u*M_PI/180;
  wc_face_angle[12][3]=WCFACEANG_30_3u*M_PI/180;
  wc_face_angle[13][0]=WCFACEANG_31_0u*M_PI/180;
  wc_face_angle[13][1]=WCFACEANG_31_1u*M_PI/180;
  wc_face_angle[13][2]=WCFACEANG_31_2u*M_PI/180;
  wc_face_angle[13][3]=WCFACEANG_31_3u*M_PI/180;
  wc_face_angle[14][0]=WCFACEANG_32_0u*M_PI/180;
  wc_face_angle[14][1]=WCFACEANG_32_1u*M_PI/180;
  wc_face_angle[14][2]=WCFACEANG_32_2u*M_PI/180;
  wc_face_angle[14][3]=WCFACEANG_32_3u*M_PI/180;
  wc_face_angle[15][0]=WCFACEANG_33_0u*M_PI/180;
  wc_face_angle[15][1]=WCFACEANG_33_1u*M_PI/180;
  wc_face_angle[15][2]=WCFACEANG_33_2u*M_PI/180;
  wc_face_angle[15][3]=WCFACEANG_33_3u*M_PI/180;
  
  /* "flipped" */
  
  wc_face_angle_F[0][0]=WCFACEANG_00_0f*M_PI/180;
  wc_face_angle_F[0][1]=WCFACEANG_00_1f*M_PI/180;
  wc_face_angle_F[0][2]=WCFACEANG_00_2f*M_PI/180;
  wc_face_angle_F[0][3]=WCFACEANG_00_3f*M_PI/180;
  wc_face_angle_F[1][0]=WCFACEANG_01_0f*M_PI/180;
  wc_face_angle_F[1][1]=WCFACEANG_01_1f*M_PI/180;
  wc_face_angle_F[1][2]=WCFACEANG_01_2f*M_PI/180;
  wc_face_angle_F[1][3]=WCFACEANG_01_3f*M_PI/180;
  wc_face_angle_F[2][0]=WCFACEANG_02_0f*M_PI/180;
  wc_face_angle_F[2][1]=WCFACEANG_02_1f*M_PI/180;
  wc_face_angle_F[2][2]=WCFACEANG_02_2f*M_PI/180;
  wc_face_angle_F[2][3]=WCFACEANG_02_3f*M_PI/180;
  wc_face_angle_F[3][0]=WCFACEANG_03_0f*M_PI/180;
  wc_face_angle_F[3][1]=WCFACEANG_03_1f*M_PI/180;
  wc_face_angle_F[3][2]=WCFACEANG_03_2f*M_PI/180;
  wc_face_angle_F[3][3]=WCFACEANG_03_3f*M_PI/180;
  
  wc_face_angle_F[4][0]=WCFACEANG_10_0f*M_PI/180;
  wc_face_angle_F[4][1]=WCFACEANG_10_1f*M_PI/180;
  wc_face_angle_F[4][2]=WCFACEANG_10_2f*M_PI/180;
  wc_face_angle_F[4][3]=WCFACEANG_10_3f*M_PI/180;
  wc_face_angle_F[5][0]=WCFACEANG_11_0f*M_PI/180;
  wc_face_angle_F[5][1]=WCFACEANG_11_1f*M_PI/180;
  wc_face_angle_F[5][2]=WCFACEANG_11_2f*M_PI/180;
  wc_face_angle_F[5][3]=WCFACEANG_11_3f*M_PI/180;
  wc_face_angle_F[6][0]=WCFACEANG_12_0f*M_PI/180;
  wc_face_angle_F[6][1]=WCFACEANG_12_1f*M_PI/180;
  wc_face_angle_F[6][2]=WCFACEANG_12_2f*M_PI/180;
  wc_face_angle_F[6][3]=WCFACEANG_12_3f*M_PI/180;
  wc_face_angle_F[7][0]=WCFACEANG_13_0f*M_PI/180;
  wc_face_angle_F[7][1]=WCFACEANG_13_1f*M_PI/180;
  wc_face_angle_F[7][2]=WCFACEANG_13_2f*M_PI/180;
  wc_face_angle_F[7][3]=WCFACEANG_13_3f*M_PI/180;
  
  wc_face_angle_F[8][0]=WCFACEANG_20_0f*M_PI/180;
  wc_face_angle_F[8][1]=WCFACEANG_20_1f*M_PI/180;
  wc_face_angle_F[8][2]=WCFACEANG_20_2f*M_PI/180;
  wc_face_angle_F[8][3]=WCFACEANG_20_3f*M_PI/180;
  wc_face_angle_F[9][0]=WCFACEANG_21_0f*M_PI/180;
  wc_face_angle_F[9][1]=WCFACEANG_21_1f*M_PI/180;
  wc_face_angle_F[9][2]=WCFACEANG_21_2f*M_PI/180;
  wc_face_angle_F[9][3]=WCFACEANG_21_3f*M_PI/180;
  wc_face_angle_F[10][0]=WCFACEANG_22_0f*M_PI/180;
  wc_face_angle_F[10][1]=WCFACEANG_22_1f*M_PI/180;
  wc_face_angle_F[10][2]=WCFACEANG_22_2f*M_PI/180;
  wc_face_angle_F[10][3]=WCFACEANG_22_3f*M_PI/180;
  wc_face_angle_F[11][0]=WCFACEANG_23_0f*M_PI/180;
  wc_face_angle_F[11][1]=WCFACEANG_23_1f*M_PI/180;
  wc_face_angle_F[11][2]=WCFACEANG_23_2f*M_PI/180;
  wc_face_angle_F[11][3]=WCFACEANG_23_3f*M_PI/180;
  
  wc_face_angle_F[12][0]=WCFACEANG_30_0f*M_PI/180;
  wc_face_angle_F[12][1]=WCFACEANG_30_1f*M_PI/180;
  wc_face_angle_F[12][2]=WCFACEANG_30_2f*M_PI/180;
  wc_face_angle_F[12][3]=WCFACEANG_30_3f*M_PI/180;
  wc_face_angle_F[13][0]=WCFACEANG_31_0f*M_PI/180;
  wc_face_angle_F[13][1]=WCFACEANG_31_1f*M_PI/180;
  wc_face_angle_F[13][2]=WCFACEANG_31_2f*M_PI/180;
  wc_face_angle_F[13][3]=WCFACEANG_31_3f*M_PI/180;
  wc_face_angle_F[14][0]=WCFACEANG_32_0f*M_PI/180;
  wc_face_angle_F[14][1]=WCFACEANG_32_1f*M_PI/180;
  wc_face_angle_F[14][2]=WCFACEANG_32_2f*M_PI/180;
  wc_face_angle_F[14][3]=WCFACEANG_32_3f*M_PI/180;
  wc_face_angle_F[15][0]=WCFACEANG_33_0f*M_PI/180;
  wc_face_angle_F[15][1]=WCFACEANG_33_1f*M_PI/180;
  wc_face_angle_F[15][2]=WCFACEANG_33_2f*M_PI/180;
  wc_face_angle_F[15][3]=WCFACEANG_33_3f*M_PI/180;
  
  
  /* wc_face_angle_F[0][0]=WCFACEANG_0_0f*M_PI/180; */
  /* wc_face_angle_F[0][1]=WCFACEANG_0_1f*M_PI/180; */
  /* wc_face_angle_F[0][2]=WCFACEANG_0_2f*M_PI/180; */
  /* wc_face_angle_F[0][3]=WCFACEANG_0_3f*M_PI/180; */
  /* wc_face_angle_F[1][0]=WCFACEANG_1_0f*M_PI/180; */
  /* wc_face_angle_F[1][1]=WCFACEANG_1_1f*M_PI/180; */
  /* wc_face_angle_F[1][2]=WCFACEANG_1_2f*M_PI/180; */
  /* wc_face_angle_F[1][3]=WCFACEANG_1_3f*M_PI/180; */
  /* wc_face_angle_F[2][0]=WCFACEANG_2_0f*M_PI/180; */
  /* wc_face_angle_F[2][1]=WCFACEANG_2_1f*M_PI/180; */
  /* wc_face_angle_F[2][2]=WCFACEANG_2_2f*M_PI/180; */
  /* wc_face_angle_F[2][3]=WCFACEANG_2_3f*M_PI/180; */
  /* wc_face_angle_F[3][0]=WCFACEANG_3_0f*M_PI/180; */
  /* wc_face_angle_F[3][1]=WCFACEANG_3_1f*M_PI/180; */
  /* wc_face_angle_F[3][2]=WCFACEANG_3_2f*M_PI/180; */
  /* wc_face_angle_F[3][3]=WCFACEANG_3_3f*M_PI/180; */
  
  
  bp_face_angle[0][0]=BPFACEANG_0_0*M_PI/180;
  bp_face_angle[0][1]=BPFACEANG_0_1*M_PI/180;
  bp_face_angle[0][2]=BPFACEANG_0_2*M_PI/180;
  bp_face_angle[0][3]=BPFACEANG_0_3*M_PI/180;
  
  bp_face_angle[1][0]=BPFACEANG_1_0*M_PI/180;
  bp_face_angle[1][1]=BPFACEANG_1_1*M_PI/180;
  bp_face_angle[1][2]=BPFACEANG_1_2*M_PI/180;
  bp_face_angle[1][3]=BPFACEANG_1_3*M_PI/180;
  
  bp_face_angle[2][0]=BPFACEANG_2_0*M_PI/180;
  bp_face_angle[2][1]=BPFACEANG_2_1*M_PI/180;
  bp_face_angle[2][2]=BPFACEANG_2_2*M_PI/180;
  bp_face_angle[2][3]=BPFACEANG_2_3*M_PI/180;
  
  bp_face_angle[3][0]=BPFACEANG_3_0*M_PI/180;
  bp_face_angle[3][1]=BPFACEANG_3_1*M_PI/180;
  bp_face_angle[3][2]=BPFACEANG_3_2*M_PI/180;
  bp_face_angle[3][3]=BPFACEANG_3_3*M_PI/180;

  bp_spec_ang[0][0]=0;
  bp_spec_ang[0][1]=0;
  bp_spec_ang[1][0]=0;
  bp_spec_ang[1][1]=0;
  bp_spec_ang[2][0]=BPSPECANG_2_0*M_PI/180;
  bp_spec_ang[2][1]=BPSPECANG_2_1*M_PI/180;
  bp_spec_ang[3][0]=BPSPECANG_3_0*M_PI/180;
  bp_spec_ang[3][1]=BPSPECANG_3_1*M_PI/180;
  
  
  int typ1, typ2;
  for(typ1=0;typ1<N_BASES;typ1++){
    for(typ2=0;typ2<N_BASES;typ2++){
      mc_sameface_min_aWW[typ1][typ2]=0;
      mc_sameface_max_aWW[typ1][typ2]=0;
      mc_sameface_min_aSS[typ1][typ2]=0;
      mc_sameface_max_aSS[typ1][typ2]=0;
      mc_sameface_min_pSS[typ1][typ2]=0;
      mc_sameface_max_pSS[typ1][typ2]=0;
    }
  }
  //THE CHOSEN ONES
  //aWW (cWW)
  mc_sameface_min_aWW[TYP_ADENINE][TYP_ADENINE]    = 1.14153254365;
  mc_sameface_max_aWW[TYP_ADENINE][TYP_ADENINE]    = 4.5851124186;
  mc_sameface_min_aWW[TYP_ADENINE][TYP_CYTOSINE]   = 2.43519039392;
  mc_sameface_max_aWW[TYP_ADENINE][TYP_CYTOSINE]   = 4.51699272869;
  mc_sameface_min_aWW[TYP_CYTOSINE][TYP_ADENINE]   = 2.0439432003 ;
  mc_sameface_max_aWW[TYP_CYTOSINE][TYP_ADENINE]   = 4.51563985998;
  mc_sameface_min_aWW[TYP_CYTOSINE][TYP_CYTOSINE]  = 1.54315635608 ;
  mc_sameface_max_aWW[TYP_CYTOSINE][TYP_CYTOSINE]  = 4.70037958612;
  mc_sameface_min_aWW[TYP_CYTOSINE][TYP_URACIL]    = 2.37848706341;
  mc_sameface_max_aWW[TYP_CYTOSINE][TYP_URACIL]    = 2.4420986618;
  mc_sameface_min_aWW[TYP_URACIL][TYP_CYTOSINE]    = 1.48160065086;
  mc_sameface_max_aWW[TYP_URACIL][TYP_CYTOSINE]    = 1.57630378601;
  mc_sameface_min_aWW[TYP_URACIL][TYP_URACIL]      = 2.38535590415;
  mc_sameface_max_aWW[TYP_URACIL][TYP_URACIL]      = 3.33030576258;

  //aSS (cSS) 
  mc_sameface_min_aSS[TYP_ADENINE][TYP_ADENINE]    = -3.58571575903;
  mc_sameface_max_aSS[TYP_ADENINE][TYP_ADENINE]    = -1.48206652553;
  mc_sameface_min_aSS[TYP_ADENINE][TYP_CYTOSINE]   =-6.11625805136;
  mc_sameface_max_aSS[TYP_ADENINE][TYP_CYTOSINE]   = -2.94407618885;
  mc_sameface_min_aSS[TYP_ADENINE][TYP_GUANINE]    =-3.08299390844;
  mc_sameface_max_aSS[TYP_ADENINE][TYP_GUANINE]    = -1.20686066605;
  mc_sameface_min_aSS[TYP_CYTOSINE][TYP_ADENINE]   =-0.441259159053;
  mc_sameface_max_aSS[TYP_CYTOSINE][TYP_ADENINE]   = 1.2933659403;
  mc_sameface_min_aSS[TYP_GUANINE][TYP_ADENINE]    =-3.44677113957;
  mc_sameface_max_aSS[TYP_GUANINE][TYP_ADENINE]    = -1.84310388525;
  mc_sameface_min_aSS[TYP_GUANINE][TYP_GUANINE]    =-4.71924230423;
  mc_sameface_max_aSS[TYP_GUANINE][TYP_GUANINE]    = -2.22795662094;
  
  
  //pSS (tSS)
  mc_sameface_min_pSS[TYP_ADENINE][TYP_ADENINE]    =-2.77779413743;
  mc_sameface_max_pSS[TYP_ADENINE][TYP_ADENINE]    = -1.49428583788;
  
  /* wc_face_angle[0][1]=WCFACEANG_0_1; */
  /* wc_face_angle[0][2]=; */
  /* wc_face_angle[1][0]=WCFACEANG_1_0; */
  /* wc_face_angle[1][1]=WCFACEANG_1_1+WCFACEANG_1_0; */
  /* wc_face_angle[1][2]=-2.35619; */
  /* wc_face_angle[2][0]=WCFACEANG_2_0; */
  /* wc_face_angle[2][1]=WCFACEANG_2_1+WCFACEANG_2_0; */
  /* wc_face_angle[2][2]=-2.35619; */
  /* wc_face_angle[3][0]=WCFACEANG_3_0; */
  /* wc_face_angle[3][1]=WCFACEANG_3_1+WCFACEANG_3_0; */
  /* wc_face_angle[3][2]=-2.35619; */
  
  
  
  
  double wc_ener=MC_calculate_total_wc_energy(nt_n,-1,  rx, ry, rz);
  MC_update_wc_lists(nt_n);
}

/* int MC_get_wc_face(int typ_ind, double *r_vec, int face){ */
/*   //int basetyp=typ_ind/N_BASES; */
/*   int RET=-1; */
/*   double angle=0.0; */
/*   double a0, a1, a2, a3; */
/*   if(face==NT_UNFLIPPED)    {a0=wc_face_angle[typ_ind][0]  ;a1=wc_face_angle[typ_ind][1]  ;a2=wc_face_angle[typ_ind][2]  ;a3=wc_face_angle[typ_ind][3]  ;} */
/*   else if(face==NT_FLIPPED) {a0=wc_face_angle_F[typ_ind][0];a1=wc_face_angle_F[typ_ind][1];a2=wc_face_angle_F[typ_ind][2];a3=wc_face_angle_F[typ_ind][3];} */
  
/*   angle=atan2(r_vec[1], r_vec[0]); */
/*   if(angle<-M_PI/2.0) angle+=2.0*M_PI; */
  
/*   if(angle < a0 || angle >= a3) RET=-1; */
/*   if(angle > a0 && angle <= a1) RET= WC_FACE_SUGAR; */
/*   if(angle > a1 && angle <= a2) RET= WC_FACE_WATSCRICK; */
/*   if(angle > a2 && angle <= a3) RET= WC_FACE_HOOGSTEEN; */
  
/*   return RET; */
/* } */

int MC_get_wc_face(int typ_ind, double *r_vec, int face){
  //int basetyp=typ_ind/N_BASES;
  int RET=-1;
  double angle=0.0;
  double a0, a1, a2, a3;
  if(face==NT_UNFLIPPED)    {a0=wc_face_angle[typ_ind][0]  ;a1=wc_face_angle[typ_ind][1]  ;a2=wc_face_angle[typ_ind][2]  ;a3=wc_face_angle[typ_ind][3]  ;}
  else if(face==NT_FLIPPED) {a0=wc_face_angle_F[typ_ind][0];a1=wc_face_angle_F[typ_ind][1];a2=wc_face_angle_F[typ_ind][2];a3=wc_face_angle_F[typ_ind][3];}
  
  angle=atan2(r_vec[1], r_vec[0]);
  if(angle<-M_PI/2.0) angle+=2.0*M_PI;
  
  if(angle < a0 || angle >= a3) RET=-1;
  if(angle > a0 && angle <= a1) RET= WC_FACE_SUGAR;
  if(angle > a1 && angle <= a2) RET= WC_FACE_WATSCRICK;
  if(angle > a2 && angle <= a3) RET= WC_FACE_HOOGSTEEN;
  
  return RET;
}

int MC_get_bp_face(int typ_ind, double *r_vec){
  //int basetyp=typ_ind/N_BASES;
  int RET=-1;
  double angle=0.0;
  double a0, a1, a2, a3;
  a0=bp_face_angle[typ_ind][0];a1=bp_face_angle[typ_ind][1];a2=bp_face_angle[typ_ind][2];a3=bp_face_angle[typ_ind][3];
  //  else if(face==NT_FLIPPED) {a0=wc_face_angle_F[typ_ind][0];a1=wc_face_angle_F[typ_ind][1];a2=wc_face_angle_F[typ_ind][2];a3=wc_face_angle_F[typ_ind][3];}
  
  angle=atan2(r_vec[1], r_vec[0]);
  if(angle<-M_PI/2.0) angle+=2.0*M_PI;
  
  if(angle < a0 || angle >= a3) RET=-1;
  if(angle > a0 && angle <= a1) RET= WC_FACE_SUGAR;
  if(angle > a1 && angle <= a2) RET= WC_FACE_WATSCRICK;
  if(angle > a2 && angle <= a3) RET= WC_FACE_HOOGSTEEN;
  
  return RET;
}

void MC_update_wc_lists(int nt_n){
  int n,f;
  for(n=0;n<nt_n;n++){
    for(f=0;f<WC_FACES;f++){
      wc_face_ene[n][f]=wc_face_ene_trial[n][f];
      wc_is_paired[n][f]=wc_is_paired_trial[n][f];
      wc_face_ind[n][f]=wc_face_ind_trial[n][f];
    }
  bp_face_ene[n]  = bp_face_ene_trial[n];
  bp_is_paired[n] = bp_is_paired_trial[n];
  bp_face_ind[n]  = bp_face_ind_trial[n];
  }
}

void MC_check_total_wc_list(int nt_n){
  //here i do the first loop over bases, and the second, over bases and phosphates
  
  int nt_c, nt_ne, n;
  int face_c, face_ne; 
  int ret=0;
  for(nt_c=0;nt_c<nt_n;nt_c++){
    //for(n=0;n<vl_n_pairs[nt_c];n++){
    //nt_ne=vl_neighbor_tables[nt_c][n];
    for(face_c=0;face_c<WC_FACES;face_c++){
      //CASE BASE-BASE
      if(wc_face_ind_trial[nt_c][face_c]>=0  &&  wc_face_ind_trial[nt_c][face_c] < nt_n){  
	nt_ne=wc_face_ind_trial[nt_c][face_c];
	for(face_ne=0;face_ne<WC_FACES;face_ne++){
	  if(wc_face_ind_trial[nt_ne][face_ne]==nt_c){
	    wc_is_paired_trial[nt_c][face_c]=1;
	    wc_is_paired_trial[nt_ne][face_ne]=1;
	    face_ne=WC_FACES+1;
	  }
	}
      }
      //CASE BASE-PHOSPHATE
      else if(wc_face_ind_trial[nt_c][face_c] >= nt_n){
      	nt_ne=wc_face_ind_trial[nt_c][face_c]-nt_n;  //this is the nt of the phosphate
      	//printf("%d   %d %d\n", nt_ne, nt_c, face_c);
	if(bp_face_ind_trial[nt_ne]==nt_c){
      	  wc_is_paired_trial[nt_c][face_c]=2;
      	  bp_is_paired_trial[nt_ne]=1;
      	}
      }
    }
  }
}

void MC_is_wc_corresp(int nt1, int nt2, int nt_n, int *ret_array){
  int a,b;
  int retbb=-1, retbp=-1, retpb=-1;
  
  //BASE-BASE
  
  //nt1 and nt2 are paired through the face b of nt2
  for(a=0;a<WC_FACES;a++){
    if(wc_is_paired[nt1][a]==1){
      for(b=0;b<WC_FACES;b++){
	if(wc_face_ind[nt2][b]==nt1){
	  retbb=b;
	  b=(int)WC_FACES;
	  a=(int)WC_FACES;
	}
      }
    }
  }
  
  //BASE-PHOSPHATE
  
  //nt2 has the phosphate, through face a of nt1
  for(a=0;a<WC_FACES;a++){
    if(wc_is_paired[nt1][a]==2){
      if(bp_face_ind[nt2]==nt1){
	retbp=WC_FACE_PHOSPHATE;
	a=(int)WC_FACES;
      }
    }
  }
  
  //PHOSPHATE-BASE
  
  //nt1 has the phosphate, nt2 the base through face b
  if(bp_is_paired[nt1]==1){
    for(b=0;b<WC_FACES;b++){
      if(wc_is_paired[nt2][b]==2 && wc_face_ind[nt2][b]==nt1+nt_n){
	retpb=b;
	b=(int)WC_FACES;
      }
    }
  }
  
  ret_array[0]=retbb;
  ret_array[1]=retbp;
  ret_array[2]=retpb;
}



double MC_calculate_local_wc_energy(int nt_n, int nt_c,  double *rx, double *ry, double *rz){
  int n,face,nt_t, face2;
  double sugdistsq;
  int flag=0;
  int up_flag=0;
  
  //printf("Local energ, nt_c = %d  %lf %lf %lf\n", nt_c, mc_temp_x[N_PARTS_PER_NT*nt_c],mc_temp_y[N_PARTS_PER_NT*nt_c],mc_temp_z[N_PARTS_PER_NT*nt_c] );
  for(nt_t=0;nt_t<nt_n;nt_t++){
    for(face=0;face<WC_FACES;face++){
      wc_is_paired_trial[nt_t][face]=-1;// wc_is_paired[nt_t][face];
      wc_face_ene_trial[nt_t][face] = wc_face_ene[nt_t][face];
      wc_face_ind_trial[nt_t][face] = wc_face_ind[nt_t][face];
      if(nt_t==nt_c){
	wc_face_ene_trial[nt_t][face] = 0;
	wc_face_ind_trial[nt_t][face] = -1;
      }
    }
    bp_is_paired_trial[nt_t]=-1;// wc_is_paired[nt_t][face];
    bp_face_ene_trial[nt_t] = bp_face_ene[nt_t];
    bp_face_ind_trial[nt_t] = bp_face_ind[nt_t];
    if(nt_t==nt_c){
      bp_face_ene_trial[nt_t] = 0;
      bp_face_ind_trial[nt_t] = -1;
    }
  }
  
  /* printf("LOCAL (zero)\n"); */
  /* for(n=0;n<nt_n;n++){ */
  /*   for(face=0;face<WC_FACES;face++) */
  /*     printf("%d   %d  %d   %lf\t", n, face, wc_face_ind[n][face], wc_face_ene[n][face]); */
  /*   printf("  - %d   %d    %lf\n", n,bp_face_ind[n], bp_face_ene[n]); */
  /* } */
  
  //for(nt_t=0;nt_t<nt_n;nt_t++){
  // THIS LOOP HAS TO BE DONE IN THE NEIGHBOR LIST OF nt_c
  for(n=0;n<vl_n_pairs[nt_c];n++){
    nt_t=vl_neighbor_tables[nt_c][n]; 
    up_flag=MC_update_min_wc(nt_t, nt_c, rx, ry, rz, nt_n);
    if(up_flag==1)
      MC_find_min_wc(nt_t, nt_c, rx, ry, rz, nt_n); //if the nt was pointing to nt_c but its not anymore, it can point somewhere else
  }  
  
  
  MC_check_total_wc_list(nt_n);
  
  /* printf("LOCAL (old)\n"); */
  /* for(n=0;n<nt_n;n++){ */
  /*   for(face=0;face<WC_FACES;face++) */
  /*     printf("%d   %d  %d   %lf\t", n, face, wc_face_ind[n][face], wc_face_ene[n][face]); */
  /*   printf("  - %d   %d    %lf\n", n,bp_face_ind[n], bp_face_ene[n]); */
    
  /* } */
  /*   printf("LOCAL (new)\n"); */
  /* for(n=0;n<nt_n;n++){ */
  /*   for(face=0;face<WC_FACES;face++) */
  /*     printf("%d   %d  %d   %lf\t", n, face, wc_face_ind_trial[n][face], wc_face_ene_trial[n][face]); */
  /*   printf("  - %d   %d    %lf\n", n,bp_face_ind_trial[n], bp_face_ene_trial[n]); */
    
  /* } */
  
  //collect energies
  // THIS IS TEMPORARY, ONLY FOR COMPARING WITH THE OTHER IMPLEMENTATION
  double wc_tot_energ=0, bp_tot_energ=0;
  double wc_tot_energ_old=0, bp_tot_energ_old=0;
  
  for(n=0;n<nt_n;n++)
    for(face=0;face<WC_FACES;face++){
      if(wc_is_paired_trial[n][face]==1){
	wc_tot_energ+=wc_face_ene_trial[n][face]; 
      } 
      else if(wc_is_paired_trial[n][face]==2){
	bp_tot_energ+=wc_face_ene_trial[n][face];
	//printf("LTRIA  %d and %d,  %lf\n", n, wc_face_ind_trial[n][face]-nt_n, wc_face_ene_trial[n][face] ); 
      }
      if(wc_is_paired[n][face]==1){ 
	wc_tot_energ_old+=wc_face_ene[n][face]; 
      } 
      else if(wc_is_paired[n][face]==2){
	bp_tot_energ_old+=wc_face_ene[n][face];
	//printf("LPREV  %d and %d,  %lf\n", n, wc_face_ind[n][face]-nt_n, wc_face_ene[n][face] ); 
      }
    }
  if(flag==1) exit(1);
  return (wc_tot_energ-wc_tot_energ_old)/2.0+(bp_tot_energ-bp_tot_energ_old);
}



int MC_update_min_wc(int nt_ne, int nt_c, double *rx, double *ry, double *rz, int nt_n){
  int at_c=N_PARTS_PER_NT*nt_c, at_ne, typ_ind, tf, ret=0;
  double r_vec[DIM]={0,0,0}, r_vec_inv[DIM]={0,0,0},bas_vec[DIM], r;
  double basedistsq=0, theta, basephosdistsq=0, phosbasedistsq=0, sphi;
    
  at_ne=N_PARTS_PER_NT*nt_ne;
  basedistsq=calc_min_dist_sq(mc_temp_x[at_c],mc_temp_y[at_c], mc_temp_z[at_c], rx[at_ne], ry[at_ne], rz[at_ne]);
  basephosdistsq = calc_min_dist_sq(mc_temp_x[at_c],mc_temp_y[at_c], mc_temp_z[at_c], rx[at_ne+IPHO], ry[at_ne+IPHO], rz[at_ne+IPHO]); //this base, neighboring phosphate
  phosbasedistsq = calc_min_dist_sq(mc_temp_x[at_c+IPHO],mc_temp_y[at_c+IPHO], mc_temp_z[at_c+IPHO], rx[at_ne], ry[at_ne], rz[at_ne]); //this phosphate, neighboring base
  int temp_ind_wc[WC_FACES], temp_ind_bp, temp_ind_pb;
  double temp_ene_wc[WC_FACES], temp_ene_bp, temp_ene_pb;
  /* if((nt_c==4 && nt_ne==7) || (nt_c==7 && nt_ne==4)) printf("bd %d  %d   -   %lf %lf\n", nt_c, nt_ne, basephosdistsq, phosbasedistsq); */
  for(tf=0;tf<WC_FACES;tf++){
    if(wc_face_ind[nt_ne][tf]==nt_c || wc_face_ind[nt_ne][tf]==(nt_c+nt_n))
      wc_face_ene_trial[nt_ne][tf]=0;
  }
  if(bp_face_ind[nt_ne]==nt_c)
    bp_face_ene_trial[nt_ne]=0;
  /* for(tf=0;tf<WC_FACES;tf++) */
  /* if((nt_c==4 && nt_ne==7) || (nt_c==7 && nt_ne==4)) printf("pf %d  %d   -   %lf %lf\n", nt_c, nt_ne,wc_face_ene_trial[nt_ne][tf], wc_face_ene[nt_ne][tf]); */
  //BASE-BASE CHECK
  //if(basedistsq<mc_r_cut_sq){
  if(basedistsq<mc_wc_rcut_sq){
    calc_min_vec(rx[at_ne],ry[at_ne],rz[at_ne],mc_temp_x[at_c],mc_temp_y[at_c],mc_temp_z[at_c],bas_vec, &r);
    proj_on_nt(bas_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(bas_vec, rx,ry,rz,nt_ne, r_vec_inv);
    theta=calc_wc_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z,rx,ry,rz,nt_c,nt_ne); 
    sphi =calc_wc_secdihed(mc_temp_x, mc_temp_y, mc_temp_z,rx,ry,rz,nt_c,nt_ne);
    MC_calc_nnN_watscric(mc_types[at_c], mc_types[at_ne], r_vec, r_vec_inv, theta, sphi, nt_c, nt_ne);
  }
  //BASE-PHOSPHATE CHECK
  if(basephosdistsq<mc_bph_rcut_sq){
    calc_min_vec(rx[at_ne+IPHO],ry[at_ne+IPHO],rz[at_ne+IPHO],mc_temp_x[at_c],mc_temp_y[at_c],mc_temp_z[at_c],bas_vec, &r);
    proj_on_nt(bas_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    MC_calc_npN_phosbase(mc_types[at_c], r_vec, nt_c, nt_ne, nt_n); // the second argument is the owner of the phosphate
  }
  //PHOSPHATE-BASE CHECK
  if(phosbasedistsq<mc_bph_rcut_sq){
    calc_min_vec(mc_temp_x[at_c+IPHO], mc_temp_y[at_c+IPHO], mc_temp_z[at_c+IPHO], rx[at_ne], ry[at_ne], rz[at_ne],bas_vec, &r);
    proj_on_nt(bas_vec, rx, ry, rz, nt_ne, r_vec);
    MC_calc_npN_phosbase(mc_types[at_ne], r_vec, nt_ne, nt_c, nt_n); // the second argument is the owner of the phosphate
  }
  
  //IF NT POINTED TO NT_C BUT :
  //A) IT DOESNT POINT ANYMORE
  //B) ITS ENERGY IS LARGER THAN THE PREVIOUS ONE - IT CAN NOW POINT SOMEWHERE ELSE - THIS INCLUDES CASE A)
  //IF NT DIDN'T POINT TO NT_C : IF NOW INTERACTS WITH NT_C
  for(tf=0;tf<WC_FACES;tf++){
    if( ((wc_face_ind[nt_ne][tf]==nt_c || wc_face_ind[nt_ne][tf]==(nt_c+nt_n)) && wc_face_ene_trial[nt_ne][tf] > wc_face_ene[nt_ne][tf]) || (wc_face_ind_trial[nt_ne][tf]==nt_c && wc_face_ene_trial[nt_ne][tf]< wc_face_ene[nt_ne][tf]))
      ret=1;
  }
  if((bp_face_ind[nt_ne]==nt_c  && bp_face_ene_trial[nt_ne] > bp_face_ene[nt_ne]) || (bp_face_ind_trial[nt_ne]==nt_c &&  bp_face_ene_trial[nt_ne]<bp_face_ene[nt_ne]) )
    ret=1;
    
  return ret;
}

double MC_calculate_total_wc_energy(int nt_n, int nt_c,  double *rx, double *ry, double *rz){
  int n,face, nt_i;
  double sugdistsq;
  double wc_tot_energ=0, bp_tot_energ=0;
  for(nt_i=0;nt_i<nt_n;nt_i++){
    //printf("%d  ", nt_i);
    for(face=0;face<WC_FACES;face++){
      wc_is_paired_trial[nt_i][face]=-1;
    }
    bp_is_paired_trial[nt_i]=-1;
  }
  
  for(nt_i=0;nt_i<nt_n;nt_i++){
    MC_find_min_wc(nt_i, nt_c, rx, ry, rz, nt_n);
  }
  MC_check_total_wc_list(nt_n);
  
  /* printf("TOTAL\n"); */
  /* for(n=0;n<nt_n;n++){ */
  /*   for(face=0;face<WC_FACES;face++) */
  /*     printf("%d   %d  %d   %lf\t", n, face, wc_face_ind_trial[n][face], wc_face_ene_trial[n][face]); */
  /*   printf("  - %d   %d    %lf\n", n,bp_face_ind_trial[n], bp_face_ene_trial[n]); */
    
  /* } */

  
  //printf("TOTAL\n");
  for(nt_i=0;nt_i<nt_n;nt_i++)
    for(face=0;face<WC_FACES;face++){
      //printf("%d  %d ", nt_i, face);// < SEGMENTATION FAULT!
      
      //printf("   %lf (%d)   ", wc_face_ene_trial[nt_i][face], wc_face_ind_trial[nt_i][face]);
      //if(nt_i>nt_n)printf("   %lf (%d)   ", bp_face_ene_trial[nt_i-nt_n], bp_face_ind_trial[nt_i]);
      if(wc_is_paired_trial[nt_i][face]==1){
	wc_tot_energ+=wc_face_ene_trial[nt_i][face];
      }
      
      else if(wc_is_paired_trial[nt_i][face]==2){
      	bp_tot_energ+=wc_face_ene_trial[nt_i][face]; 
	
      } 
    }
  //printf("\n");
  /* //printf("bp tot energ  %lf\n", bp_tot_energ); */
  return wc_tot_energ/2.0 + bp_tot_energ;
}


void MC_find_min_wc(int nt_c, int nt_trial, double *rx, double *ry, double *rz, int nt_n){
  int n, tf, nt_ne, fl, ss, inv_fl, uu;
  int ind_wc=1, ind_phos=1;
  for(tf=0;tf<WC_FACES;tf++){
    wc_face_ene_trial[nt_c][tf]=0;
    wc_face_ind_trial[nt_c][tf]=-1;
  }
  bp_face_ene_trial[nt_c]=0;
  bp_face_ind_trial[nt_c]=-1;
  
  /* printf("%d   \n", nt_c); */
  for(n=0;n<vl_n_pairs[nt_c];n++){
    nt_ne=vl_neighbor_tables[nt_c][n];
    ind_wc=1;
    ind_phos=1;
#ifdef SECONDSTC
    if(wc_sstruct_N[nt_c]==0 && wc_sstruct_N[nt_ne]==0){
	ind_wc=1;
	ind_phos=1;
	inv_fl=1;
	fl=1;
    }
    else if(wc_sstruct_N[nt_c]==-1 || wc_sstruct_N[nt_ne]==-1){
      ind_wc=0;
      ind_phos=0;
      inv_fl=0;
      fl=0;
    }
    else{
      ind_wc=0;
      ind_phos=0;
      inv_fl=0;
      fl=0;
      if(wc_sstruct_N[nt_ne]==0){
	//printf("easy inv\n");
	inv_fl=1;}
      else{
	for(uu=0;uu<wc_sstruct_N[nt_ne];uu++){
	  if(wc_sstruct_neigh[nt_ne][uu*2] == nt_c){
	    inv_fl=1;
	    break;
	  }
	}
      }
      if(wc_sstruct_N[nt_c]==0){
	//printf("easy \n");
	fl=1;
      }
      else{
	for(ss=0;ss<wc_sstruct_N[nt_c];ss++){
	  if(wc_sstruct_neigh[nt_c][ss*2] == nt_ne){
	    fl=1;
	    break;
	  }	}      }
      if(inv_fl==1 && fl==1){
	//printf("CHOSEN!\n");
	if(wc_sstruct_neigh[nt_c][ss*2+1]< 2){
	  ind_wc=1;
	}
	
	//else if(wc_sstruct_neigh[nt_c][ss*2+1]==1){
	//  ind_phos=1;
	//}
	if(wc_sstruct_neigh[nt_c][ss*2+1]>=1)
	  ind_phos=wc_sstruct_neigh[nt_c][ss*2+1];
      }    }
    
    //printf("%d  %d   %d %d\t%d %d\n", nt_c, nt_ne, ind_wc, ind_phos, fl, inv_fl);
    if(ind_wc!=0 && ind_phos!=0)
#endif
      //MC_find_min_wc_per_nt(nt_c, nt_ne, nt_trial, rx, ry, rz, nt_n, 1, 1);
    
      MC_find_min_wc_per_nt(nt_c, nt_ne, nt_trial, rx, ry, rz, nt_n, ind_wc, ind_phos);
  }
  
  //THE ADDITIONAL NEIGHBOR PAIRS
  ind_phos=0;
  ind_wc=0;
  nt_ne=nt_c-1;
  if(nt_ne>=0){
    ind_phos=2;
    ind_wc=1;
    MC_find_min_wc_per_nt(nt_c, nt_ne, nt_trial, rx, ry, rz, nt_n, ind_wc, ind_phos);
  }
  nt_ne=nt_c+1;
  ind_phos=0;
  ind_wc=0;
  if(nt_ne<nt_n){
    ind_phos=3;
    ind_wc=1;
    MC_find_min_wc_per_nt(nt_c, nt_ne, nt_trial, rx, ry, rz, nt_n, ind_wc, ind_phos);
  }
}

void MC_find_min_wc_per_nt(int nt_c, int nt_ne, int nt_trial, double *rx, double *ry, double *rz, int nt_n, int wc_flag, int bp_flag){
  int at_c=N_PARTS_PER_NT*nt_c, at_ne=N_PARTS_PER_NT*nt_ne, typ_ind;
  double r_vec[DIM]={0,0,0}, r_vec_inv[DIM]={0,0,0},bas_vec[DIM], r;
  double basedistsq=0, basephosdistsq=0, phosbasedistsq=0, theta=0, sphi;
  //at_ne=N_PARTS_PER_NT*nt_ne;
  
  if(nt_trial == nt_c) {
    if(wc_flag==1)
      basedistsq     = calc_min_dist_sq(mc_temp_x[at_c],mc_temp_y[at_c], mc_temp_z[at_c], rx[at_ne], ry[at_ne], rz[at_ne]);
    if(bp_flag==1 || bp_flag==2)
      basephosdistsq = calc_min_dist_sq(mc_temp_x[at_c],mc_temp_y[at_c], mc_temp_z[at_c], rx[at_ne+IPHO], ry[at_ne+IPHO], rz[at_ne+IPHO]); //this base, neighboring phosphate
    if(bp_flag==1 || bp_flag==3)
      phosbasedistsq = calc_min_dist_sq(mc_temp_x[at_c+IPHO],mc_temp_y[at_c+IPHO], mc_temp_z[at_c+IPHO], rx[at_ne], ry[at_ne], rz[at_ne]); //this phosphate, neighboring base
  }
  else if(nt_trial==nt_ne) {
    if(wc_flag==1)
      basedistsq     = calc_min_dist_sq(mc_temp_x[at_ne],mc_temp_y[at_ne], mc_temp_z[at_ne], rx[at_c], ry[at_c], rz[at_c]);
    if(bp_flag==1 || bp_flag==2)
      basephosdistsq = calc_min_dist_sq(mc_temp_x[at_ne+IPHO],mc_temp_y[at_ne+IPHO], mc_temp_z[at_ne+IPHO], rx[at_c], ry[at_c], rz[at_c]);//this base, neighboring phosphate
    if(bp_flag==1 || bp_flag==3)
      phosbasedistsq = calc_min_dist_sq(mc_temp_x[at_ne],mc_temp_y[at_ne], mc_temp_z[at_ne], rx[at_c+IPHO], ry[at_c+IPHO], rz[at_c+IPHO]);//this phosphate, neighboring base
  }
  else {
    if(wc_flag==1)
      basedistsq     = calc_min_dist_sq(rx[at_c], ry[at_c], rz[at_c], rx[at_ne], ry[at_ne], rz[at_ne]);
    if(bp_flag==1 || bp_flag==2)
      basephosdistsq = calc_min_dist_sq(rx[at_c], ry[at_c], rz[at_c], rx[at_ne+IPHO], ry[at_ne+IPHO], rz[at_ne+IPHO]);//this base, neighboring phosphate
    if(bp_flag==1 || bp_flag==3)
      phosbasedistsq = calc_min_dist_sq(rx[at_c+IPHO], ry[at_c+IPHO], rz[at_c+IPHO], rx[at_ne], ry[at_ne], rz[at_ne]);//this phosphate, neighboring base
    
  }
  
  if(basedistsq<mc_wc_rcut_sq && wc_flag==1){
    if(nt_trial == nt_c) {
      calc_min_vec(rx[at_ne],ry[at_ne],rz[at_ne],mc_temp_x[at_c],mc_temp_y[at_c],mc_temp_z[at_c],bas_vec, &r);
      proj_on_nt(bas_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
      proj_on_nt_inv(bas_vec, rx,ry,rz,nt_ne, r_vec_inv);
      theta=calc_wc_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z,rx,ry,rz,nt_c,nt_ne); 
      sphi =calc_wc_secdihed(mc_temp_x, mc_temp_y, mc_temp_z,rx,ry,rz,nt_c,nt_ne); 
    } else if(nt_trial==nt_ne){
      calc_min_vec(mc_temp_x[at_ne],mc_temp_y[at_ne], mc_temp_z[at_ne], rx[at_c],ry[at_c],rz[at_c], bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_c, r_vec);
      proj_on_nt_inv(bas_vec,mc_temp_x, mc_temp_y, mc_temp_z,nt_ne, r_vec_inv);
      theta=calc_wc_psdihedr(rx,ry,rz,mc_temp_x, mc_temp_y, mc_temp_z,nt_c,nt_ne); 
      sphi =calc_wc_secdihed(rx,ry,rz,mc_temp_x, mc_temp_y, mc_temp_z,nt_c,nt_ne); 
    } else {
      calc_min_vec(rx[at_ne],ry[at_ne],rz[at_ne],rx[at_c],ry[at_c],rz[at_c],  bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_c, r_vec);
      proj_on_nt_inv(bas_vec,rx, ry, rz,nt_ne, r_vec_inv);
      theta=calc_wc_psdihedr(rx,ry,rz,rx,ry,rz,nt_c,nt_ne); 
      sphi =calc_wc_secdihed(rx,ry,rz,rx,ry,rz,nt_c,nt_ne); 
    }
    MC_calc_nnN_watscric(mc_types[at_c], mc_types[at_ne], r_vec, r_vec_inv, theta, sphi, nt_c, nt_ne);
  }
  
  if(basephosdistsq<mc_bph_rcut_sq && (bp_flag ==1 || bp_flag==2)){ // ne-phosphate   c-base
    if(nt_trial == nt_c) {//typical case
      calc_min_vec(rx[at_ne+IPHO],ry[at_ne+IPHO],rz[at_ne+IPHO],mc_temp_x[at_c],mc_temp_y[at_c],mc_temp_z[at_c],bas_vec, &r);
      proj_on_nt(bas_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    } else if(nt_trial==nt_ne){
      calc_min_vec(mc_temp_x[at_ne+IPHO],mc_temp_y[at_ne+IPHO], mc_temp_z[at_ne+IPHO], rx[at_c],ry[at_c],rz[at_c], bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_c, r_vec);
    } else {
      calc_min_vec(rx[at_ne+IPHO],ry[at_ne+IPHO],rz[at_ne+IPHO],rx[at_c],ry[at_c],rz[at_c],  bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_c, r_vec);
    }
    if(nt_ne!=0)
      MC_calc_npN_phosbase(mc_types[at_c], r_vec, nt_c, nt_ne, nt_n); // the second argument is the owner of the phosphate
  }
  
  if(phosbasedistsq<mc_bph_rcut_sq && (bp_flag==1 || bp_flag==3)){// ne-base    c-phosphate
    if(nt_trial == nt_c) {//typical case
      calc_min_vec(mc_temp_x[at_c+IPHO], mc_temp_y[at_c+IPHO], mc_temp_z[at_c+IPHO], rx[at_ne], ry[at_ne], rz[at_ne],bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_ne, r_vec);
    } else if(nt_trial==nt_ne){
      calc_min_vec(rx[at_c+IPHO],ry[at_c+IPHO], rz[at_c+IPHO], mc_temp_x[at_ne],mc_temp_y[at_ne],mc_temp_z[at_ne], bas_vec, &r);
      proj_on_nt(bas_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_ne, r_vec);
    } else {
      calc_min_vec(rx[at_c+IPHO],ry[at_c+IPHO],rz[at_c+IPHO],rx[at_ne],ry[at_ne],rz[at_ne],  bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_ne, r_vec);
    }
    if(nt_c!=0)
      MC_calc_npN_phosbase(mc_types[at_ne], r_vec, nt_ne, nt_c, nt_n); // the second argument is the owner of the phosphate
  }
  /* printf("\n"); */
}

void MC_calc_npN_phosbase(int typ_a,  double *r_vec, int nt_a, int nt_b, int nt_n){
  //classify type of pair
  int face_a =MC_get_bp_face(typ_a , r_vec);
   
  if(face_a <0){
    return;
  }
  int typ_ind = typ_a; //for the moment, the interaction depends only on the case type
  int face_ind = face_a;
  double nbe=0;
  if(table_npN_N_0[typ_ind][0]>0)
    nbe=calc_npN_tab_phosbase(r_vec, typ_ind); // base-phosphate
  double tot_nb_bp=0;
  
  double well=nb_bp_well[typ_ind][face_ind];
  //calculate subfaces
  int subf=MC_get_BPh_subface(typ_a, r_vec, face_a);
  
  if(subf!=-1){
    //here we modify the potential well
    well=nb_bp_spec_well[typ_ind][subf];
  }
  
  if(nbe!=0 && well!=0 ){
    //printf("%d %d %d  %lf\n", nt_a, nt_b, subf, well);
    tot_nb_bp=nbe;
    tot_nb_bp-=well;
#ifdef BI_ANNEALING
    tot_nb_bp*=bia_lambda;
    if(tot_nb_bp<bia_cap)
      tot_nb_bp=bia_cap;
#endif
 
    
    if(tot_nb_bp < wc_face_ene_trial[nt_a][face_a]){
      wc_face_ene_trial[nt_a][face_a]=tot_nb_bp;
      wc_face_ind_trial[nt_a][face_a]=nt_b+nt_n;
    }
    if(tot_nb_bp < bp_face_ene_trial[nt_b]){
      bp_face_ene_trial[nt_b]=tot_nb_bp;
      bp_face_ind_trial[nt_b]=nt_a;
    }
  }
}
int MC_get_BPh_subface(int type, double *r_vec, int face){
  int ret=-1;
  double angle;
  if(type==TYP_GUANINE && face==WC_FACE_WATSCRICK){
    angle=atan2(r_vec[1], r_vec[0]);
    if(angle<-M_PI/2.0) angle+=2.0*M_PI;
    
    if(angle<bp_spec_ang[type][0]) ret=0;
    else{
      if(angle<bp_spec_ang[type][1]) ret=1;
      else ret=2;
    }
  }
  else if(type==TYP_CYTOSINE && face==WC_FACE_HOOGSTEEN){
    angle=atan2(r_vec[1], r_vec[0]);
    if(angle<-M_PI/2.0) angle+=2.0*M_PI;
    if(angle<bp_spec_ang[type][0]) ret=0;
    else{
      if(angle<bp_spec_ang[type][1]) ret=1;
      else ret=2;
    }
  }
  return ret;
}

double calc_npN_tab_phosbase(double *r_vec, int typ){
  int xin, yin, zin;
  xin=(int)((r_vec[0]-table_npN_params_0[typ][0][1]+0.5*table_npN_params_0[typ][0][0])/table_npN_params_0[typ][0][0]);
  yin=(int)((r_vec[1]-table_npN_params_0[typ][1][1]+0.5*table_npN_params_0[typ][1][0])/table_npN_params_0[typ][1][0]);
  zin=(int)((r_vec[2]-table_npN_params_0[typ][2][1]+0.5*table_npN_params_0[typ][2][0])/table_npN_params_0[typ][2][0]);
  int pos_ind = zin+table_npN_N_0[typ][2]*yin+table_npN_N_0[typ][1]*table_npN_N_0[typ][2]*xin;
  if(xin<0 || xin>=table_npN_N_0[typ][0])
    return 0;
  if(yin<0 || yin>=table_npN_N_0[typ][1])
    return 0;
  if(zin<0 || zin>=table_npN_N_0[typ][2])
    return 0;
  
  return table_npN_0[typ][pos_ind];

}

int eval_wc_secdihed(double sdih, int typ_ind, int face_ind, int ori){
  int ret=1;
  
  double sdmin=0;
  double sdmax=0;

  if(ori==NT_UNFLIPPED){
    if(sdih<0) sdih+=2*M_PI;
    sdmin=wc_secdih_min[typ_ind][face_ind];
    sdmax=wc_secdih_max[typ_ind][face_ind];
    if(sdih < sdmin || sdih > sdmax)
      ret=0;
  }
  else{
   sdmin=wc_secdih_min_F[typ_ind][face_ind];
   sdmax=wc_secdih_max_F[typ_ind][face_ind];
   if(sdih < sdmin || sdih > sdmax)
     ret=0; 
  }
  return ret;
}




int eval_sameface_excl(int typ_c, int typ_ne, int face_c, int face_ne, double *r_vec, double *r_vec_inv, int ori){
  int ret=0;
  if(face_c==WC_FACE_WATSCRICK && face_ne==WC_FACE_WATSCRICK){
    if(ori==NT_UNFLIPPED){
      if((typ_c==TYP_ADENINE && (typ_ne==TYP_ADENINE || typ_ne==TYP_CYTOSINE)) || 
	 (typ_c==TYP_CYTOSINE && (typ_ne==TYP_ADENINE || typ_ne==TYP_CYTOSINE || typ_ne==TYP_URACIL)) || 
	 (typ_c==TYP_URACIL && (typ_ne==TYP_CYTOSINE || typ_ne==TYP_URACIL))){
	ret=1;
	if((r_vec[0] < mc_sameface_min_aWW[typ_c][typ_ne] && r_vec_inv[0] > mc_sameface_max_aWW[typ_ne][typ_c]) || 
	   (r_vec[0] > mc_sameface_max_aWW[typ_c][typ_ne] && r_vec_inv[0] < mc_sameface_min_aWW[typ_ne][typ_c]))
	  ret=0;
      }
    }
  }
  else if(face_c==WC_FACE_SUGAR && face_ne==WC_FACE_SUGAR){
    if(ori==NT_UNFLIPPED){
      if((typ_c==TYP_ADENINE && (typ_ne==TYP_ADENINE || typ_ne==TYP_CYTOSINE || typ_ne==TYP_GUANINE)) || (typ_c==TYP_CYTOSINE && typ_ne==TYP_ADENINE) || 
	 (typ_c==TYP_GUANINE && (typ_ne==TYP_ADENINE || typ_ne==TYP_GUANINE))){
	ret=1;
	if((r_vec[1] < mc_sameface_min_aSS[typ_c][typ_ne] && r_vec_inv[1] > mc_sameface_max_aSS[typ_ne][typ_c]) || 
	   (r_vec[1] > mc_sameface_max_aSS[typ_c][typ_ne] && r_vec_inv[1] < mc_sameface_min_aSS[typ_ne][typ_c]))
	  ret=0;
      }
    }
    else{
      if(typ_c==TYP_ADENINE && typ_ne==TYP_ADENINE){
	ret=1;
	if((r_vec[1] < mc_sameface_min_pSS[typ_c][typ_ne] && r_vec_inv[1] > mc_sameface_max_pSS[typ_ne][typ_c]) || 
	   (r_vec[1] > mc_sameface_max_pSS[typ_c][typ_ne] && r_vec_inv[1] < mc_sameface_min_pSS[typ_ne][typ_c]))
	  ret=0;
      }
    }
    
  }

  return ret;
}

void MC_calc_nnN_watscric(int typ_c, int typ_ne, double *r_vec, double *r_vec_inv, double theta, double sphi, int nt_c, int nt_ne){
  //fflush(stdout);
  //classify type of pair
  //int face_c =MC_get_wc_face(typ_c , r_vec);
  //int face_ne=MC_get_wc_face(typ_ne, r_vec_inv);
  int face_c ;
  int face_ne;
  
  //if(face_c <0 || face_ne <0)
  //  return;
  int typ_ind     = N_BASES*typ_c  + typ_ne;
  int typ_ind_inv = N_BASES*typ_ne + typ_c;
  //int face_ind = WC_FACES*face_c + face_ne;
  int face_ind;
  double nbe=0, nbei=0, nbed=0;
  double tnbed_F=0, tnbed=0, nb_well=0;
  int tshift=0;
  int tshift_i;
  if(table_nnN_N_3[typ_ind]>0) 
    tnbed  =calc_nnN_tab_wc_psdihedr(theta,typ_ind,NT_UNFLIPPED); // watson-crick pseudodihedral
  if(table_nnN_N_3_F[typ_ind]>0) 
    tnbed_F=calc_nnN_tab_wc_psdihedr(theta,typ_ind,  NT_FLIPPED); // watson-crick pseudodihedral
  
  if(tnbed!=0 && tnbed_F==0){
    //CASE UNFLIPPED
    //face_c =MC_get_wc_face(typ_c , r_vec, NT_UNFLIPPED);
    //face_ne=MC_get_wc_face(typ_ne, r_vec_inv,NT_UNFLIPPED);
    face_c =MC_get_wc_face(typ_ind    , r_vec    , NT_UNFLIPPED);
    face_ne=MC_get_wc_face(typ_ind_inv, r_vec_inv, NT_UNFLIPPED);
    if(face_c <0 || face_ne <0)
      return;
    
    //HERE APPLY THE SAMEFACE EXCLUSION
    if(eval_sameface_excl(typ_c, typ_ne, face_c, face_ne, r_vec, r_vec_inv, NT_UNFLIPPED)==1) return;
    
    face_ind = WC_FACES*face_c + face_ne;
    //SECOND DIHEDRAL
    if(eval_wc_secdihed(sphi,  typ_ind, face_ind, NT_UNFLIPPED)==0) return;
    tshift  =table_nnN_N_2[typ_ind][3]*face_ne;
    tshift_i=table_nnN_N_2[typ_ind_inv][3]*face_c;
    //printf("shift  = %d,   %d   %d\n", tshift, table_nnN_N_2[typ_ind][3], face_ne);
    //printf("shiftt = %d,   %d   %d\n", tshift_i, table_nnN_N_2[typ_ind_inv][3], face_c);
    if(table_nnN_N_2[typ_ind][0]>0)
      nbe=calc_nnN_tab_watscric(r_vec, typ_ind, tshift, NT_UNFLIPPED); // watson-crick
    if(table_nnN_N_2_inv[typ_ind][0]>0)
      nbei=calc_nnN_tab_watscric_inv(r_vec_inv, typ_ind, tshift_i, NT_UNFLIPPED); // watson-crick
    nbed=tnbed;
    nb_well=nb_wc_well[typ_ind][face_ind];
  }
  else if(tnbed==0 && tnbed_F!=0){
    //CASE FLIPPED
    face_c =MC_get_wc_face(typ_ind    , r_vec    , NT_FLIPPED);
    face_ne=MC_get_wc_face(typ_ind_inv, r_vec_inv, NT_FLIPPED);
    if(face_c <0 || face_ne <0)
      return;
    
    
    
    //HERE APPLY THE SAMEFACE EXCLUSION
    if(eval_sameface_excl(typ_c, typ_ne, face_c, face_ne, r_vec, r_vec_inv, NT_FLIPPED)==1) return;
    
    
    face_ind = WC_FACES*face_c + face_ne;
    //SECOND DIHEDRAL
    if(eval_wc_secdihed(sphi,  typ_ind, face_ind, NT_FLIPPED)==0) return;
    
    tshift  =table_nnN_N_2[typ_ind][3]*face_ne;
    tshift_i=table_nnN_N_2[typ_ind_inv][3]*face_c;
    
    if(table_nnN_N_2_F[typ_ind][0]>0)
      nbe=calc_nnN_tab_watscric(r_vec, typ_ind, tshift, NT_FLIPPED); // watson-crick
    if(table_nnN_N_2_inv_F[typ_ind][0]>0)
      nbei=calc_nnN_tab_watscric_inv(r_vec_inv, typ_ind, tshift_i, NT_FLIPPED); // watson-crick
    nbed=tnbed_F;
    nb_well=nb_wc_well_F[typ_ind][face_ind];
  }
  double tot_nb_wc=0;
  //printf("PAIR %d %d    %lf %lf %lf   %lf   SHIFT %d %d\t\t%lf %lf %lf   %lf %lf %lf\n", nt_c, nt_ne, nbe, nbei, nbed, nb_well, tshift, tshift_i,r_vec[0], r_vec[1], r_vec[2], r_vec_inv[0], r_vec_inv[1], r_vec_inv[2]);
  if(nbe!=0 && nbei!=0 && nbed!=0 && nb_well!=0 ){
    tot_nb_wc=nbe+nbei+nbed;
    tot_nb_wc-=nb_well;
#ifdef BI_ANNEALING
    tot_nb_wc*=bia_lambda;
    if(tot_nb_wc<bia_cap)
      tot_nb_wc=bia_cap;
#endif
    if(tot_nb_wc < wc_face_ene_trial[nt_c][face_c]){
      wc_face_ene_trial[nt_c][face_c]=tot_nb_wc;
      wc_face_ind_trial[nt_c][face_c]=nt_ne;
    }
    if(tot_nb_wc < wc_face_ene_trial[nt_ne][face_ne]){
      wc_face_ene_trial[nt_ne][face_ne]=tot_nb_wc;
      wc_face_ind_trial[nt_ne][face_ne]=nt_c;
    }
  }
}


double calc_wc_psdihedr(double *rx1, double *ry1, double *rz1, double *rx2, double *ry2, double *rz2, int nt_c, int nt_neigh){
  double cent1[DIM], cent2[DIM];
  double x1[DIM], x2[DIM], y1[DIM], y2[DIM];
  double z1[DIM], z2[DIM];
  double cos_dih_wc=1;
  int at_c=N_PARTS_PER_NT*nt_c;
  x1[0]=dist_1d(rx1[at_c+1],rx1[at_c],0);  x1[1]=dist_1d(ry1[at_c+1],ry1[at_c],1);  x1[2]=dist_1d(rz1[at_c+1],rz1[at_c],2);
  y1[0]=dist_1d(rx1[at_c+2],rx1[at_c],0);  y1[1]=dist_1d(ry1[at_c+2],ry1[at_c],1);  y1[2]=dist_1d(rz1[at_c+2],rz1[at_c],2);
    
  int np=nt_neigh*N_PARTS_PER_NT;
  x2[0]=dist_1d(rx2[np+1],rx2[np],0);  x2[1]=dist_1d(ry2[np+1],ry2[np],1);  x2[2]=dist_1d(rz2[np+1],rz2[np],2);
  y2[0]=dist_1d(rx2[np+2],rx2[np],0);  y2[1]=dist_1d(ry2[np+2],ry2[np],1);  y2[2]=dist_1d(rz2[np+2],rz2[np],2);
  
  vec_prod(x1,y1,z1);
  vec_prod(x2,y2,z2);
  //if(dot_prod(z1, z1)*dot_prod(z2, z2)<0.0001) cos_dih_wc=1;
  //else cos_dih_wc=dot_prod(z1,z2)/sqrt(dot_prod(z1, z1)*dot_prod(z2, z2));
  cos_dih_wc=dot_prod(z1,z2)/sqrt(dot_prod(z1, z1)*dot_prod(z2, z2));
  return acos(cos_dih_wc);
}


double calc_wc_secdihed(double *rx1, double *ry1, double *rz1, double *rx2, double *ry2, double *rz2, int nt_c, int nt_neigh){
  double cent1[DIM], cent2[DIM];
  double x1[DIM], x2[DIM], y1[DIM], y2[DIM];
  double z1[DIM], z2[DIM];
  int at_c=N_PARTS_PER_NT*nt_c;
  x1[0]=dist_1d(rx1[at_c+1],rx1[at_c],0);  x1[1]=dist_1d(ry1[at_c+1],ry1[at_c],1);  x1[2]=dist_1d(rz1[at_c+1],rz1[at_c],2);
  y1[0]=dist_1d(rx1[at_c+2],rx1[at_c],0);  y1[1]=dist_1d(ry1[at_c+2],ry1[at_c],1);  y1[2]=dist_1d(rz1[at_c+2],rz1[at_c],2);
    
  int np=nt_neigh*N_PARTS_PER_NT;
  x2[0]=dist_1d(rx2[np+1],rx2[np],0);  x2[1]=dist_1d(ry2[np+1],ry2[np],1);  x2[2]=dist_1d(rz2[np+1],rz2[np],2);
  y2[0]=dist_1d(rx2[np+2],rx2[np],0);  y2[1]=dist_1d(ry2[np+2],ry2[np],1);  y2[2]=dist_1d(rz2[np+2],rz2[np],2);
  
  vec_prod(x1,y1,z1);
  vec_prod(x2,y2,z2);
  
  /******* HERE WE CALCULATE THE OTHER DIHEDRAL ********/

  double r12[DIM]={rx1[at_c]-rx2[np],ry1[at_c]-ry2[np],rz1[at_c]-rz2[np]};
  double modr12=sqrt(dot_prod(r12, r12));
  r12[0]=r12[0]/modr12;r12[1]=r12[1]/modr12;r12[2]=r12[2]/modr12;
  double t1[DIM], t2[DIM], t3[DIM];
  vec_prod(r12, z1, t1);
  vec_prod(r12, z2, t2);
  vec_prod(t1, t2, t3);
  double phi=atan2(dot_prod(t3, r12), dot_prod(t1,t2));
  /*****************************************************/
  return phi;
}

double calc_nnN_tab_wc_psdihedr(double theta, int typ_ind, int FLFLAG){
  int theta_in=0;
  double ret=0;
  if(FLFLAG==NT_UNFLIPPED){
    theta_in=(int)((theta-table_nnN_params_3[typ_ind][1]+0.5*table_nnN_params_3[typ_ind][0])/table_nnN_params_3[typ_ind][0]);
    if(theta_in>=table_nnN_N_3[typ_ind] || theta_in <0) return 0;
    ret= table_nnN_3[typ_ind][theta_in];
  }
  else if(FLFLAG==NT_FLIPPED){
    theta_in=(int)((theta-table_nnN_params_3_F[typ_ind][1]+0.5*table_nnN_params_3_F[typ_ind][0])/table_nnN_params_3_F[typ_ind][0]);
    if(theta_in>=table_nnN_N_3_F[typ_ind] || theta_in <0) return 0;
    ret= table_nnN_3_F[typ_ind][theta_in];
  }
  return ret;
}




double calc_nnN_tab_watscric(double *r_vec, int typ_ind, int tshift, int FLFLAG){
  int xin, yin, zin;
  double ret=0;
  //printf("tshift=%d \n", tshift);
  if(FLFLAG==NT_UNFLIPPED){
    xin=(int)((r_vec[0]-table_nnN_params_2[typ_ind][0][1]+0.5*table_nnN_params_2[typ_ind][0][0])/table_nnN_params_2[typ_ind][0][0]);
    yin=(int)((r_vec[1]-table_nnN_params_2[typ_ind][1][1]+0.5*table_nnN_params_2[typ_ind][1][0])/table_nnN_params_2[typ_ind][1][0]);
    zin=(int)((r_vec[2]-table_nnN_params_2[typ_ind][2][1]+0.5*table_nnN_params_2[typ_ind][2][0])/table_nnN_params_2[typ_ind][2][0]);
    int pos_ind = zin+table_nnN_N_2[typ_ind][2]*yin+table_nnN_N_2[typ_ind][1]*table_nnN_N_2[typ_ind][2]*xin + tshift;
    //printf("u %d %d %d    %d\n", xin, yin, zin, pos_ind);
    if(xin<0 || xin>=table_nnN_N_2[typ_ind][0])
      return 0;
    if(yin<0 || yin>=table_nnN_N_2[typ_ind][1])
      return 0;
    if(zin<0 || zin>=table_nnN_N_2[typ_ind][2])
      return 0;
    
    ret= table_nnN_2[typ_ind][pos_ind];
  }
  else if(FLFLAG==NT_FLIPPED){
    xin=(int)((r_vec[0]-table_nnN_params_2_F[typ_ind][0][1]+0.5*table_nnN_params_2_F[typ_ind][0][0])/table_nnN_params_2_F[typ_ind][0][0]);
    yin=(int)((r_vec[1]-table_nnN_params_2_F[typ_ind][1][1]+0.5*table_nnN_params_2_F[typ_ind][1][0])/table_nnN_params_2_F[typ_ind][1][0]);
    zin=(int)((r_vec[2]-table_nnN_params_2_F[typ_ind][2][1]+0.5*table_nnN_params_2_F[typ_ind][2][0])/table_nnN_params_2_F[typ_ind][2][0]);
    int pos_ind = zin+table_nnN_N_2_F[typ_ind][2]*yin+table_nnN_N_2_F[typ_ind][1]*table_nnN_N_2_F[typ_ind][2]*xin + tshift;
    //printf("f %d %d %d    %d\n", xin, yin, zin, pos_ind);
    if(xin<0 || xin>=table_nnN_N_2_F[typ_ind][0])
      return 0;
    if(yin<0 || yin>=table_nnN_N_2_F[typ_ind][1])
      return 0;
    if(zin<0 || zin>=table_nnN_N_2_F[typ_ind][2])
      return 0;
    ret= table_nnN_2_F[typ_ind][pos_ind];
  }
  return ret;
}

double calc_nnN_tab_watscric_inv(double *r_vec, int typ_ind, int tshift, int FLFLAG){
  int xin, yin, zin;
  double ret=0;
  if(FLFLAG==NT_UNFLIPPED){
    xin=(int)((r_vec[0]-table_nnN_params_2_inv[typ_ind][0][1]+0.5*table_nnN_params_2_inv[typ_ind][0][0])/table_nnN_params_2_inv[typ_ind][0][0]);
    yin=(int)((r_vec[1]-table_nnN_params_2_inv[typ_ind][1][1]+0.5*table_nnN_params_2_inv[typ_ind][1][0])/table_nnN_params_2_inv[typ_ind][1][0]);
    zin=(int)((r_vec[2]-table_nnN_params_2_inv[typ_ind][2][1]+0.5*table_nnN_params_2_inv[typ_ind][2][0])/table_nnN_params_2_inv[typ_ind][2][0]);
    int pos_ind = zin+table_nnN_N_2_inv[typ_ind][2]*yin+table_nnN_N_2_inv[typ_ind][1]*table_nnN_N_2_inv[typ_ind][2]*xin + tshift;
    if(xin<0 || xin>=table_nnN_N_2_inv[typ_ind][0]) return 0;
    //printf("x out of range, %d max %d\t\t%lf   %lf\n", xin, table_nnN_N_0[typ_ind][0], r_vec[0], table_nnN_N_0[typ_ind][0]*table_nnN_params_0[typ_ind][0][0]+table_nnN_params_0[typ_ind][0][1]);
    if(yin<0 || yin>=table_nnN_N_2_inv[typ_ind][1]) return 0;
    //printf("y out of range, %d max %d\t\t%lf   %lf\n", yin, table_nnN_N_0[typ_ind][1], r_vec[1], table_nnN_N_0[typ_ind][1]*table_nnN_params_0[typ_ind][1][0]+table_nnN_params_0[typ_ind][1][1]);
    if(zin<0 || zin>=table_nnN_N_2_inv[typ_ind][2]) return 0;
    //printf("z out of range, %d max %d\t\t%lf   %lf\n", zin, table_nnN_N_0[typ_ind][2], r_vec[2], table_nnN_N_0[typ_ind][2]*table_nnN_params_0[typ_ind][2][0]+table_nnN_params_0[typ_ind][2][1]);
    ret= table_nnN_2_inv[typ_ind][pos_ind];
  }
  else if(FLFLAG==NT_FLIPPED){
    xin=(int)((r_vec[0]-table_nnN_params_2_inv_F[typ_ind][0][1]+0.5*table_nnN_params_2_inv_F[typ_ind][0][0])/table_nnN_params_2_inv_F[typ_ind][0][0]);
    yin=(int)((r_vec[1]-table_nnN_params_2_inv_F[typ_ind][1][1]+0.5*table_nnN_params_2_inv_F[typ_ind][1][0])/table_nnN_params_2_inv_F[typ_ind][1][0]);
    zin=(int)((r_vec[2]-table_nnN_params_2_inv_F[typ_ind][2][1]+0.5*table_nnN_params_2_inv_F[typ_ind][2][0])/table_nnN_params_2_inv_F[typ_ind][2][0]);
    int pos_ind = zin+table_nnN_N_2_inv_F[typ_ind][2]*yin+table_nnN_N_2_inv_F[typ_ind][1]*table_nnN_N_2_inv_F[typ_ind][2]*xin + tshift;
    if(xin<0 || xin>=table_nnN_N_2_inv_F[typ_ind][0]) return 0;
    //printf("x out of range, %d max %d\t\t%lf   %lf\n", xin, table_nnN_N_0[typ_ind][0], r_vec[0], table_nnN_N_0[typ_ind][0]*table_nnN_params_0[typ_ind][0][0]+table_nnN_params_0[typ_ind][0][1]);
    if(yin<0 || yin>=table_nnN_N_2_inv_F[typ_ind][1]) return 0;
    //printf("y out of range, %d max %d\t\t%lf   %lf\n", yin, table_nnN_N_0[typ_ind][1], r_vec[1], table_nnN_N_0[typ_ind][1]*table_nnN_params_0[typ_ind][1][0]+table_nnN_params_0[typ_ind][1][1]);
    if(zin<0 || zin>=table_nnN_N_2_inv_F[typ_ind][2]) return 0;
    //printf("z out of range, %d max %d\t\t%lf   %lf\n", zin, table_nnN_N_0[typ_ind][2], r_vec[2], table_nnN_N_0[typ_ind][2]*table_nnN_params_0[typ_ind][2][0]+table_nnN_params_0[typ_ind][2][1]);
    ret= table_nnN_2_inv_F[typ_ind][pos_ind];
  }
  return ret;
}

