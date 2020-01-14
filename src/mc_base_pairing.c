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
  
  bp_contacts[TYP_ADENINE][WC_FACE_SUGAR][0]=0   ;                bp_contacts[TYP_ADENINE][WC_FACE_SUGAR][1]=0;                   bp_contacts[TYP_ADENINE][WC_FACE_SUGAR][2]=0;
  bp_contacts[TYP_ADENINE][WC_FACE_WATSCRICK][0]=BP_CONTACT_A_W_X;bp_contacts[TYP_ADENINE][WC_FACE_WATSCRICK][1]=BP_CONTACT_A_W_Y;bp_contacts[TYP_ADENINE][WC_FACE_WATSCRICK][2]=BP_CONTACT_A_W_Z;
  bp_contacts[TYP_ADENINE][WC_FACE_HOOGSTEEN][0]=BP_CONTACT_A_H_X;bp_contacts[TYP_ADENINE][WC_FACE_HOOGSTEEN][1]=BP_CONTACT_A_H_Y;bp_contacts[TYP_ADENINE][WC_FACE_HOOGSTEEN][2]=BP_CONTACT_A_H_Z;
  
  bp_contacts[TYP_CYTOSINE][WC_FACE_SUGAR][0]=0                    ;bp_contacts[TYP_CYTOSINE][WC_FACE_SUGAR][1]=0                    ;bp_contacts[TYP_CYTOSINE][WC_FACE_SUGAR][2]=0;
  bp_contacts[TYP_CYTOSINE][WC_FACE_WATSCRICK][0]=BP_CONTACT_C_W_X ;bp_contacts[TYP_CYTOSINE][WC_FACE_WATSCRICK][1]=BP_CONTACT_C_W_Y ;bp_contacts[TYP_CYTOSINE][WC_FACE_WATSCRICK][2]=BP_CONTACT_C_W_Z;
  bp_contacts[TYP_CYTOSINE][WC_FACE_HOOGSTEEN][0]=BP_CONTACT_C_H_X ;bp_contacts[TYP_CYTOSINE][WC_FACE_HOOGSTEEN][1]=BP_CONTACT_C_H_Y ;bp_contacts[TYP_CYTOSINE][WC_FACE_HOOGSTEEN][2]=BP_CONTACT_C_H_Z;
  
  bp_contacts[TYP_GUANINE][WC_FACE_SUGAR][0]=    BP_CONTACT_G_S_X ;bp_contacts[TYP_GUANINE][WC_FACE_SUGAR][1]=    BP_CONTACT_G_S_Y   ;bp_contacts[TYP_GUANINE][WC_FACE_SUGAR][2]=    BP_CONTACT_G_S_Z;
  bp_contacts[TYP_GUANINE][WC_FACE_WATSCRICK][0]=BP_CONTACT_G_W1_X;bp_contacts[TYP_GUANINE][WC_FACE_WATSCRICK][1]=BP_CONTACT_G_W1_Y  ;bp_contacts[TYP_GUANINE][WC_FACE_WATSCRICK][2]=BP_CONTACT_G_W1_Z;
  bp_contacts[TYP_GUANINE][WC_FACE_HOOGSTEEN][0]=0                ;bp_contacts[TYP_GUANINE][WC_FACE_HOOGSTEEN][1]=0                  ;bp_contacts[TYP_GUANINE][WC_FACE_HOOGSTEEN][2]=0               ;

  bp_contacts[TYP_URACIL][WC_FACE_SUGAR][0]=0                    ;bp_contacts[TYP_URACIL][WC_FACE_SUGAR][1]=0                    ;bp_contacts[TYP_URACIL][WC_FACE_SUGAR][2]=0;
  bp_contacts[TYP_URACIL][WC_FACE_WATSCRICK][0]=BP_CONTACT_U_W_X ;bp_contacts[TYP_URACIL][WC_FACE_WATSCRICK][1]=BP_CONTACT_U_W_Y ;bp_contacts[TYP_URACIL][WC_FACE_WATSCRICK][2]=BP_CONTACT_U_W_Z;
  bp_contacts[TYP_URACIL][WC_FACE_HOOGSTEEN][0]=0                ;bp_contacts[TYP_URACIL][WC_FACE_HOOGSTEEN][1]=0                ;bp_contacts[TYP_URACIL][WC_FACE_HOOGSTEEN][2]=0               ;
    
  bp_hydros[TYP_ADENINE][WC_FACE_SUGAR][0]=0;                 bp_hydros[TYP_ADENINE][WC_FACE_SUGAR][1]=0;                 bp_hydros[TYP_ADENINE][WC_FACE_SUGAR][2]=0;
  bp_hydros[TYP_ADENINE][WC_FACE_WATSCRICK][0]=BP_HYDRO_A_W_X;bp_hydros[TYP_ADENINE][WC_FACE_WATSCRICK][1]=BP_HYDRO_A_W_Y;bp_hydros[TYP_ADENINE][WC_FACE_WATSCRICK][2]=BP_HYDRO_A_W_Z;
  bp_hydros[TYP_ADENINE][WC_FACE_HOOGSTEEN][0]=BP_HYDRO_A_H_X;bp_hydros[TYP_ADENINE][WC_FACE_HOOGSTEEN][1]=BP_HYDRO_A_H_Y;bp_hydros[TYP_ADENINE][WC_FACE_HOOGSTEEN][2]=BP_HYDRO_A_H_Z;
  
  bp_hydros[TYP_CYTOSINE][WC_FACE_SUGAR][0]=0                    ;bp_hydros[TYP_CYTOSINE][WC_FACE_SUGAR][1]=0                    ;bp_hydros[TYP_CYTOSINE][WC_FACE_SUGAR][2]=0;
  bp_hydros[TYP_CYTOSINE][WC_FACE_WATSCRICK][0]=BP_HYDRO_C_W_X ;bp_hydros[TYP_CYTOSINE][WC_FACE_WATSCRICK][1]=BP_HYDRO_C_W_Y ;bp_hydros[TYP_CYTOSINE][WC_FACE_WATSCRICK][2]=BP_HYDRO_C_W_Z;
  bp_hydros[TYP_CYTOSINE][WC_FACE_HOOGSTEEN][0]=BP_HYDRO_C_H_X ;bp_hydros[TYP_CYTOSINE][WC_FACE_HOOGSTEEN][1]=BP_HYDRO_C_H_Y ;bp_hydros[TYP_CYTOSINE][WC_FACE_HOOGSTEEN][2]=BP_HYDRO_C_H_Z;
  
  bp_hydros[TYP_GUANINE][WC_FACE_SUGAR][0]=    BP_HYDRO_G_S_X ;bp_hydros[TYP_GUANINE][WC_FACE_SUGAR][1]=    BP_HYDRO_G_S_Y   ;bp_hydros[TYP_GUANINE][WC_FACE_SUGAR][2]=    BP_HYDRO_G_S_Z;
  bp_hydros[TYP_GUANINE][WC_FACE_WATSCRICK][0]=BP_HYDRO_G_W1_X;bp_hydros[TYP_GUANINE][WC_FACE_WATSCRICK][1]=BP_HYDRO_G_W1_Y  ;bp_hydros[TYP_GUANINE][WC_FACE_WATSCRICK][2]=BP_HYDRO_G_W1_Z;
  bp_hydros[TYP_GUANINE][WC_FACE_HOOGSTEEN][0]=0              ;bp_hydros[TYP_GUANINE][WC_FACE_HOOGSTEEN][1]=0                ;bp_hydros[TYP_GUANINE][WC_FACE_HOOGSTEEN][2]=0             ;

  bp_hydros[TYP_URACIL][WC_FACE_SUGAR][0]=0                    ;bp_hydros[TYP_URACIL][WC_FACE_SUGAR][1]=0                    ;bp_hydros[TYP_URACIL][WC_FACE_SUGAR][2]=0;
  bp_hydros[TYP_URACIL][WC_FACE_WATSCRICK][0]=BP_HYDRO_U_W_X ;bp_hydros[TYP_URACIL][WC_FACE_WATSCRICK][1]=BP_HYDRO_U_W_Y ;bp_hydros[TYP_URACIL][WC_FACE_WATSCRICK][2]=BP_HYDRO_U_W_Z;
  bp_hydros[TYP_URACIL][WC_FACE_HOOGSTEEN][0]=0              ;bp_hydros[TYP_URACIL][WC_FACE_HOOGSTEEN][1]=0              ;bp_hydros[TYP_URACIL][WC_FACE_HOOGSTEEN][2]=0             ;
  bp_G_cont_W2[0]=BP_CONTACT_G_W2_X ;bp_G_cont_W2[1]=BP_CONTACT_G_W2_Y ;bp_G_cont_W2[2]=BP_CONTACT_G_W2_Z;
  bp_G_hydr_W2[0]=BP_HYDRO_G_W2_X   ;bp_G_hydr_W2[1]=BP_HYDRO_G_W2_Y   ;bp_G_hydr_W2[2]=BP_HYDRO_G_W2_Z;
  
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
  //printf("updating\n");
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
    if(nt_t==nt_c || nt_t==nt_c+1){
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
  
  for(nt_t=0;nt_t<nt_n;nt_t++){
    if(nt_t!=nt_c){
      // THIS LOOP HAS TO BE DONE IN THE NEIGHBOR LIST OF nt_c
      //for(n=0;n<vl_n_pairs[nt_c];n++){
      //nt_t=vl_neighbor_tables[nt_c][n]; 
      up_flag=MC_update_min_wc(nt_t, nt_c, rx, ry, rz, nt_n);
      //up_flag=1;
      if(up_flag==1 || up_flag==3)
	MC_find_min_wc(nt_t, nt_c, rx, ry, rz, nt_n); //if the nt was pointing to nt_c but its not anymore, it can point somewhere else
      //if(up_flag==2 || up_flag==3)
    }
  } 
  //this should be only in the case where nt_c and nt_c+1 are bonded!
  if(nt_c+1<nt_n)
    if(MC_are_neighbors(nt_c, nt_c+1))
      MC_find_min_wc(nt_c+1, nt_c, rx, ry, rz, nt_n); //if the nt was pointing to nt_c but its not anymore, it can point somewhere else
  
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
  /* printf("LOCAL BP   %lf    %lf   %lf    ", bp_tot_energ, bp_tot_energ_old, (bp_tot_energ-bp_tot_energ_old)); */
  /* //printf("L(o)  - %d %lf  (%d %lf , %d %lf , %d %lf)\t",nt_c, bp_face_ind[nt_c], bp_face_ene[nt_c], wc_face_ind[nt_c][0], wc_face_ene[nt_c][0],wc_face_ind[nt_c][1], wc_face_ene[nt_c][1],wc_face_ind[nt_c][2], wc_face_ene[nt_c][2]); printf("L  -  %d %lf  (%d %lf , %d %lf , %d %lf)\n",bp_face_ind_trial[nt_c], bp_face_ene_trial[nt_c], wc_face_ind_trial[nt_c][0], wc_face_ene_trial[nt_c][0],wc_face_ind_trial[nt_c][1], wc_face_ene_trial[nt_c][1],wc_face_ind_trial[nt_c][2], wc_face_ene_trial[nt_c][2]); */
  /* printf("L(o)   %d  -   %d %lf  (%d %lf , %d %lf , %d %lf)\t",nt_c, bp_face_ind[nt_c], bp_face_ene[nt_c], wc_face_ind[nt_c][0], wc_face_ene[nt_c][0],wc_face_ind[nt_c][1], wc_face_ene[nt_c][1],wc_face_ind[nt_c][2], wc_face_ene[nt_c][2]); printf("L  -  %d %lf  (%d %lf , %d %lf , %d %lf)\n",bp_face_ind_trial[nt_c], bp_face_ene_trial[nt_c], wc_face_ind_trial[nt_c][0], wc_face_ene_trial[nt_c][0],wc_face_ind_trial[nt_c][1], wc_face_ene_trial[nt_c][1],wc_face_ind_trial[nt_c][2], wc_face_ene_trial[nt_c][2]); */

  /* printf("L(o)   %d  -   %d %lf  (%d %lf , %d %lf , %d %lf)\t",6, bp_face_ind[6], bp_face_ene[6], wc_face_ind[6][0], wc_face_ene[6][0],wc_face_ind[6][1], wc_face_ene[6][1],wc_face_ind[6][2], wc_face_ene[6][2]); printf("L  -  %d %lf  (%d %lf , %d %lf , %d %lf)\n",bp_face_ind_trial[6], bp_face_ene_trial[6], wc_face_ind_trial[6][0], wc_face_ene_trial[6][0],wc_face_ind_trial[6][1], wc_face_ene_trial[6][1],wc_face_ind_trial[6][2], wc_face_ene_trial[6][2]); */
  /* printf("L(o)   %d  -   %d %lf  (%d %lf , %d %lf , %d %lf)\t",7, bp_face_ind[7], bp_face_ene[7], wc_face_ind[7][0], wc_face_ene[7][0],wc_face_ind[7][1], wc_face_ene[7][1],wc_face_ind[7][2], wc_face_ene[7][2]); printf("L  -  %d %lf  (%d %lf , %d %lf , %d %lf)\n",bp_face_ind_trial[7], bp_face_ene_trial[7], wc_face_ind_trial[7][0], wc_face_ene_trial[7][0],wc_face_ind_trial[7][1], wc_face_ene_trial[7][1],wc_face_ind_trial[7][2], wc_face_ene_trial[7][2]); */
  return (wc_tot_energ-wc_tot_energ_old)/2.0+(bp_tot_energ-bp_tot_energ_old);
  
}



int MC_update_min_wc(int nt_ne, int nt_c, double *rx, double *ry, double *rz, int nt_n){
   //here it seems that we see what happens between nt_ne and nt_c, which was just moved
  int at_c=N_PARTS_PER_NT*nt_c, at_ne, typ_ind, tf, ret=0, nt_add=nt_c+1,at_add=(nt_c+1)*N_PARTS_PER_NT;
  //int at_c=N_PARTS_PER_NT*nt_c, at_ne, typ_ind, tf, ret=0;
  double r_vec[DIM]={0,0,0}, r_vec_inv[DIM]={0,0,0},bas_vec[DIM], r;
  double basedistsq=0, theta, basephosdistsq=0, phosbasedistsq=0, sphi;
  double BAS[DIM], SUG1[DIM], SUG2[DIM], PHO[DIM], op1[DIM], op2[DIM], bas_op1[DIM], bas_op2[DIM];
  int bp_case=-1;
  int glp1, glp2;
  at_ne=N_PARTS_PER_NT*nt_ne;
  basedistsq=calc_min_dist_sq(mc_temp_x[at_c],mc_temp_y[at_c], mc_temp_z[at_c], rx[at_ne], ry[at_ne], rz[at_ne]);
  basephosdistsq = calc_min_dist_sq(mc_temp_x[at_c],mc_temp_y[at_c], mc_temp_z[at_c], rx[at_ne+IPHO], ry[at_ne+IPHO], rz[at_ne+IPHO]); //this base, neighboring phosphate
  phosbasedistsq = calc_min_dist_sq(mc_temp_x[at_c+IPHO],mc_temp_y[at_c+IPHO], mc_temp_z[at_c+IPHO], rx[at_ne], ry[at_ne], rz[at_ne]); //this phosphate, neighboring base
  //note that we also have to calculate neighboring phosphate - this base
  double addphosbasedsq;
  // addphosbasedsq = calc_min_dist_sq(rx[at_add+IPHO],ry[at_add+IPHO], rz[at_add+IPHO], rx[at_ne], ry[at_ne], rz[at_ne]); //next phosphate, this base
  
  if(nt_c<nt_n-1) addphosbasedsq = calc_min_dist_sq(rx[at_add+IPHO],ry[at_add+IPHO], rz[at_add+IPHO], rx[at_ne], ry[at_ne], rz[at_ne]); //next phosphate, this base
  else {
    addphosbasedsq=mc_bph_rcut_sq+100;
    //printf("d=%lf\n", calc_min_dist_sq(rx[at_add+IPHO],ry[at_add+IPHO], rz[at_add+IPHO], rx[at_ne], ry[at_ne], rz[at_ne]));
  }
  
  int temp_ind_wc[WC_FACES], temp_ind_bp, temp_ind_pb;
  double temp_ene_wc[WC_FACES], temp_ene_bp, temp_ene_pb;
  /* if((nt_c==4 && nt_ne==7) || (nt_c==7 && nt_ne==4)) printf("bd %d  %d   -   %lf %lf\n", nt_c, nt_ne, basephosdistsq, phosbasedistsq); */
  for(tf=0;tf<WC_FACES;tf++){
    if(wc_face_ind[nt_ne][tf]==nt_c || wc_face_ind[nt_ne][tf]==(nt_c+nt_n))
      wc_face_ene_trial[nt_ne][tf]=0;
  }
  if(bp_face_ind[nt_ne]==nt_c)
    bp_face_ene_trial[nt_ne]=0;
  //BASE-BASE CHECK
  if(basedistsq<mc_wc_rcut_sq){
    calc_min_vec(rx[at_ne],ry[at_ne],rz[at_ne],mc_temp_x[at_c],mc_temp_y[at_c],mc_temp_z[at_c],bas_vec, &r);
    proj_on_nt(bas_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    proj_on_nt_inv(bas_vec, rx,ry,rz,nt_ne, r_vec_inv);
    theta=calc_wc_psdihedr(mc_temp_x, mc_temp_y, mc_temp_z,rx,ry,rz,nt_c,nt_ne); 
    sphi =calc_wc_secdihed(mc_temp_x, mc_temp_y, mc_temp_z,rx,ry,rz,nt_c,nt_ne);
    glp1 =mc_temp_glyc[nt_c] + mc_temp_puck[nt_c];
    glp2 =mc_glyc[nt_ne] + mc_puck[nt_ne];
    MC_calc_nnN_watscric(mc_types[at_c], mc_types[at_ne], r_vec, r_vec_inv, theta, sphi, nt_c, nt_ne, glp1, glp2);
  }
  //BASE-PHOSPHATE CHECK
  if(basephosdistsq<mc_bph_rcut_sq){
    calc_min_vec(rx[at_ne+IPHO],ry[at_ne+IPHO],rz[at_ne+IPHO],mc_temp_x[at_c],mc_temp_y[at_c],mc_temp_z[at_c],bas_vec, &r);
    proj_on_nt(bas_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
    
    //ne has the phosphate!, so ne-1 has the SUG1
    if(nt_ne-1==nt_c) //case 3, rx, temp, temp
      bp_case=3;
    else  // case 2, rx, temp, rx
      bp_case=2;
    MC_calc_npN_phosbase(mc_types[at_c], r_vec, nt_c, nt_ne, nt_n, rx, ry, rz, bp_case );// the second argument is the owner of the phosphate
  }
  //PHOSPHATE-BASE CHECK
  if(phosbasedistsq<mc_bph_rcut_sq){
    calc_min_vec(mc_temp_x[at_c+IPHO], mc_temp_y[at_c+IPHO], mc_temp_z[at_c+IPHO], rx[at_ne], ry[at_ne], rz[at_ne],bas_vec, &r);
    proj_on_nt(bas_vec, rx, ry, rz, nt_ne, r_vec);
    
    //nt_c has the phosphate!!, so nt_c -1 has the SUG1, and this is case 4, temp, rx, rx
    bp_case=4;
    MC_calc_npN_phosbase(mc_types[at_ne], r_vec, nt_ne, nt_c, nt_n, rx, ry, rz,bp_case);  // the second argument is the owner of the phosphate
  }
  
  //FINAL BASE-PHOSPHATE CHECK
  if(addphosbasedsq<mc_bph_rcut_sq){
    //at_xt=at_c+N_PARTS_PER_NT;
    calc_min_vec(rx[at_add+IPHO],ry[at_add+IPHO],rz[at_add+IPHO],rx[at_ne],ry[at_ne],rz[at_ne],bas_vec, &r);
    proj_on_nt(bas_vec, rx, ry, rz, nt_ne, r_vec);
    //ne has the phosphate!, so ne-1 has the SUG1
    bp_case=2; //check this!!!
    MC_calc_npN_phosbase(mc_types[at_ne], r_vec, nt_ne, nt_add, nt_n, rx, ry, rz, bp_case );// the second argument is the owner of the phosphate
  }
  
  //IF NT POINTED TO NT_C BUT :
  //A) IT DOESNT POINT ANYMORE
  //B) ITS ENERGY IS LARGER THAN THE PREVIOUS ONE - IT CAN NOW POINT SOMEWHERE ELSE - THIS INCLUDES CASE A)
  //IF NT DIDN'T POINT TO NT_C : IF NOW INTERACTS WITH NT_C
  int ret2=0;
  for(tf=0;tf<WC_FACES;tf++){
    ////if( ((wc_face_ind[nt_ne][tf]==nt_c || wc_face_ind[nt_ne][tf]==(nt_c+nt_n)) && wc_face_ene_trial[nt_ne][tf] > wc_face_ene[nt_ne][tf]) || (wc_face_ind_trial[nt_ne][tf]==nt_c && wc_face_ene_trial[nt_ne][tf]< wc_face_ene[nt_ne][tf]))
    //if( ((wc_face_ind[nt_ne][tf]==nt_c || wc_face_ind[nt_ne][tf]==(nt_c+nt_n)) && wc_face_ene_trial[nt_ne][tf] > wc_face_ene[nt_ne][tf]) || ((wc_face_ind_trial[nt_ne][tf]==nt_c || wc_face_ind_trial[nt_ne][tf]==(nt_c+nt_n)) && wc_face_ene_trial[nt_ne][tf]< wc_face_ene[nt_ne][tf])
    //|| (wc_face_ind[nt_ne][tf]==(nt_add+nt_n) && wc_face_ene_trial[nt_ne][tf] > wc_face_ene[nt_ne][tf]) || (wc_face_ind_trial[nt_ne][tf]==nt_add && wc_face_ene_trial[nt_ne][tf]<wc_face_ene[nt_ne][tf])	)
    if( ((wc_face_ind[nt_ne][tf]==nt_c || wc_face_ind[nt_ne][tf]==(nt_c+nt_n)) && wc_face_ene_trial[nt_ne][tf] > wc_face_ene[nt_ne][tf])
	|| ( ( (wc_face_ind_trial[nt_ne][tf]==nt_c && wc_face_ind[nt_ne][tf]!=nt_c) || ( wc_face_ind_trial[nt_ne][tf]==(nt_c+nt_n) && wc_face_ind[nt_ne][tf]!=(nt_c+nt_n)) )    && wc_face_ene_trial[nt_ne][tf]< wc_face_ene[nt_ne][tf]))
      ret=1;
    if( (   wc_face_ind[nt_ne][tf]==(nt_add+nt_n) && wc_face_ene_trial[nt_ne][tf] > wc_face_ene[nt_ne][tf]) 
	|| (   wc_face_ind_trial[nt_ne][tf]==nt_add && wc_face_ene_trial[nt_ne][tf]  < wc_face_ene[nt_ne][tf])	)
      ret2=2;
  }
  //if((bp_face_ind[nt_ne]==nt_c  && bp_face_ene_trial[nt_ne] > bp_face_ene[nt_ne]) || (bp_face_ind_trial[nt_ne]==nt_c &&  bp_face_ene_trial[nt_ne]<bp_face_ene[nt_ne]) )
  if(   (bp_face_ind[nt_ne]==nt_c  && bp_face_ene_trial[nt_ne] > bp_face_ene[nt_ne]) || (bp_face_ind_trial[nt_ne]==nt_c && bp_face_ind[nt_ne]!=nt_c &&  bp_face_ene_trial[nt_ne]<bp_face_ene[nt_ne]) )
    ret=1;
  if(nt_add<nt_n)
    if(   (bp_face_ind[nt_add]==nt_ne &&  bp_face_ene_trial[nt_add]<bp_face_ene[nt_add]) ||  (bp_face_ind_trial[nt_add]==nt_ne &&  bp_face_ind[nt_add]!=nt_ne && bp_face_ene_trial[nt_add]<bp_face_ene[nt_ne]) )
      ret2=2;
  
  //return 1;
  return ret;//+ret2;
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
  
  //printf("TOTAL\n");
  for(nt_i=0;nt_i<nt_n;nt_i++)
    for(face=0;face<WC_FACES;face++){
      if(wc_is_paired_trial[nt_i][face]==1){
	wc_tot_energ+=wc_face_ene_trial[nt_i][face];
      }
      
      else if(wc_is_paired_trial[nt_i][face]==2){
      	bp_tot_energ+=wc_face_ene_trial[nt_i][face]; 
	
      } 
    }
  
  //printf("total  BP   %lf\n", bp_tot_energ);
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
    /* if(wc_sstruct_N[nt_c]==0 && wc_sstruct_N[nt_ne]==0){ */
    /* 	ind_wc=1; */
    /* 	ind_phos=1; */
    /* 	inv_fl=1; */
    /* 	fl=1; */
    /* } */
    /* else if(wc_sstruct_N[nt_c]==-1 || wc_sstruct_N[nt_ne]==-1){ */
    /*   ind_wc=0; */
    /*   ind_phos=0; */
    /*   inv_fl=0; */
    /*   fl=0; */
    /* } */
    /* else{ */
    /*   ind_wc=0; */
    /*   ind_phos=0; */
    /*   inv_fl=0; */
    /*   fl=0; */
    /*   if(wc_sstruct_N[nt_ne]==0){ */
    /* 	//printf("easy inv\n"); */
    /* 	inv_fl=1;} */
    /*   else{ */
    /* 	for(uu=0;uu<wc_sstruct_N[nt_ne];uu++){ */
    /* 	  if(wc_sstruct_neigh[nt_ne][uu*2] == nt_c){ */
    /* 	    inv_fl=1; */
    /* 	    break; */
    /* 	  } */
    /* 	} */
    /*   } */
    /*   if(wc_sstruct_N[nt_c]==0){ */
    /* 	//printf("easy \n"); */
    /* 	fl=1; */
    /*   } */
    /*   else{ */
    /* 	for(ss=0;ss<wc_sstruct_N[nt_c];ss++){ */
    /* 	  if(wc_sstruct_neigh[nt_c][ss*2] == nt_ne){ */
    /* 	    fl=1; */
    /* 	    break; */
    /* 	  }	}      } */
    /*   if(inv_fl==1 && fl==1){ */
    /* 	//printf("CHOSEN!\n"); */
    /* 	if(wc_sstruct_neigh[nt_c][ss*2+1]< 2){ */
    /* 	  ind_wc=1; */
    /* 	} */
	
    /* 	//else if(wc_sstruct_neigh[nt_c][ss*2+1]==1){ */
    /* 	//  ind_phos=1; */
    /* 	//} */
    /* 	if(wc_sstruct_neigh[nt_c][ss*2+1]>=1) */
    /* 	  ind_phos=wc_sstruct_neigh[nt_c][ss*2+1]; */
    /*   }    } */
    
    /* //printf("%d  %d   %d %d\t%d %d\n", nt_c, nt_ne, ind_wc, ind_phos, fl, inv_fl); */
    /* if(ind_wc!=0 || ind_phos!=0) */
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
  double BAS[DIM], SUG1[DIM], SUG2[DIM], PHO[DIM], op1[DIM], op2[DIM], bas_op1[DIM], bas_op2[DIM];
  int bp_case=-1;
  int glp1, glp2;
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
      glp1 =mc_temp_glyc[nt_c] + mc_temp_puck[nt_c];
      glp2 =mc_glyc[nt_ne] + mc_puck[nt_ne];
    } else if(nt_trial==nt_ne){
      calc_min_vec(mc_temp_x[at_ne],mc_temp_y[at_ne], mc_temp_z[at_ne], rx[at_c],ry[at_c],rz[at_c], bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_c, r_vec);
      proj_on_nt_inv(bas_vec,mc_temp_x, mc_temp_y, mc_temp_z,nt_ne, r_vec_inv);
      theta=calc_wc_psdihedr(rx,ry,rz,mc_temp_x, mc_temp_y, mc_temp_z,nt_c,nt_ne); 
      sphi =calc_wc_secdihed(rx,ry,rz,mc_temp_x, mc_temp_y, mc_temp_z,nt_c,nt_ne); 
      glp1 =mc_glyc[nt_c] + mc_puck[nt_c];
      glp2 =mc_temp_glyc[nt_ne] + mc_temp_puck[nt_ne];
    } else {
      calc_min_vec(rx[at_ne],ry[at_ne],rz[at_ne],rx[at_c],ry[at_c],rz[at_c],  bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_c, r_vec);
      proj_on_nt_inv(bas_vec,rx, ry, rz,nt_ne, r_vec_inv);
      theta=calc_wc_psdihedr(rx,ry,rz,rx,ry,rz,nt_c,nt_ne); 
      sphi =calc_wc_secdihed(rx,ry,rz,rx,ry,rz,nt_c,nt_ne); 
      glp1 =mc_glyc[nt_c]  + mc_puck[nt_c];
      glp2 =mc_glyc[nt_ne] + mc_puck[nt_ne];
    }
    MC_calc_nnN_watscric(mc_types[at_c], mc_types[at_ne], r_vec, r_vec_inv, theta, sphi, nt_c, nt_ne, glp1, glp2);
  }
  
  if(basephosdistsq<mc_bph_rcut_sq && (bp_flag ==1 || bp_flag==2)){ // ne-phosphate   c-base
    if(nt_trial == nt_c) {//typical case
      calc_min_vec(rx[at_ne+IPHO],ry[at_ne+IPHO],rz[at_ne+IPHO],mc_temp_x[at_c],mc_temp_y[at_c],mc_temp_z[at_c],bas_vec, &r);
      proj_on_nt(bas_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, r_vec);
      //phosphate is in rx array (nt_ne), base is in temp array (nt_c)
      if(nt_trial==nt_ne-1) //here, sug1 must be in temp (case 3), rx, temp, temp
	bp_case=3;
      else //otherwise, it is in rx (case 2), rx, temp, rx
	bp_case=2;
    } else if(nt_trial==nt_ne){
      calc_min_vec(mc_temp_x[at_ne+IPHO],mc_temp_y[at_ne+IPHO], mc_temp_z[at_ne+IPHO], rx[at_c],ry[at_c],rz[at_c], bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_c, r_vec);
      //phosphate is in temp array (nt_ne), base is in rx array (nt_c); case 4 by force, temp, rx, rx
      bp_case=4;
    } else {
      calc_min_vec(rx[at_ne+IPHO],ry[at_ne+IPHO],rz[at_ne+IPHO],rx[at_c],ry[at_c],rz[at_c],  bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_c, r_vec);
      //phosphate is in rx array (nt_ne), base is in rx array (nt_c)
      if(nt_trial==nt_ne-1) //here, sug1 must be in temp (case 1), rx, rx, temp
	bp_case=1;
      else //otherwise, it is in rx (case 0), rx, rx, rx
	bp_case=0;
    }
    if(nt_ne!=0){
      //MC_calc_npN_phosbase(mc_types[at_c], r_vec, nt_c, nt_ne, nt_n); // the second argument is the owner of the phosphate
      /** HERE WE DO THE HB CHECK **/
    /*   BAS[0]=mc_temp_x[at_c];    BAS[1]=mc_temp_y[at_c];    BAS[2]=mc_temp_z[at_c]; */
    /*   SUG1[0]=rx[at_ne+ISUG-N_PARTS_PER_NT];SUG1[1]=ry[at_ne+ISUG-N_PARTS_PER_NT];SUG1[2]=rz[at_ne+ISUG-N_PARTS_PER_NT]; */
    /*   SUG2[0]=rx[at_ne+ISUG               ];SUG2[1]=ry[at_ne+ISUG               ];SUG2[2]=rz[at_ne+ISUG               ]; */
    /*   PHO[0]=rx[at_ne+IPHO];PHO[1]=ry[at_ne+IPHO];PHO[2]=rz[at_ne+IPHO]; */
      
    /* MC_get_op1op2(mc_types[at_c], MC_get_bp_face(mc_types[at_c], r_vec), BAS, PHO, SUG1,SUG2,op1,op2); */
    /* proj_on_nt(op1, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, bas_op1); */
    /* proj_on_nt(op2, mc_temp_x, mc_temp_y, mc_temp_z, nt_c, bas_op2); */
    /* if(MC_check_bph_HB(bas_op1, bas_op2, mc_types[at_c], MC_get_bp_face(mc_types[at_c], r_vec))==1) */
      
      MC_calc_npN_phosbase(mc_types[at_c], r_vec, nt_c, nt_ne, nt_n, rx, ry, rz, bp_case); //the second argument is the owner of the phosphate
      
    }
  }
  
  
  if(phosbasedistsq<mc_bph_rcut_sq && (bp_flag==1 || bp_flag==3)){// ne-base    c-phosphate
     if(nt_trial == nt_c) {//typical case
      calc_min_vec(mc_temp_x[at_c+IPHO], mc_temp_y[at_c+IPHO], mc_temp_z[at_c+IPHO], rx[at_ne], ry[at_ne], rz[at_ne],bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_ne, r_vec);
      bp_case=4;//here, sug1 must be in rx by force (case 4), temp, rx, rx
    } else if(nt_trial==nt_ne){
      calc_min_vec(rx[at_c+IPHO],ry[at_c+IPHO], rz[at_c+IPHO], mc_temp_x[at_ne],mc_temp_y[at_ne],mc_temp_z[at_ne], bas_vec, &r);
      proj_on_nt(bas_vec, mc_temp_x, mc_temp_y, mc_temp_z, nt_ne, r_vec);
      //here, phsphate is in rx, base in temp
      if(nt_trial==nt_c-1) // sug1 is temp, case 3, rx, temp, temp
	bp_case=3;
      else //otherwise, it is in rx (case 2), rx, temp, rx
	bp_case=2;
    } else {
      calc_min_vec(rx[at_c+IPHO],ry[at_c+IPHO],rz[at_c+IPHO],rx[at_ne],ry[at_ne],rz[at_ne],  bas_vec, &r);
      proj_on_nt(bas_vec, rx, ry, rz, nt_ne, r_vec);
      //here, phosphate is in rx, base in rx
      if(nt_trial==nt_c-1) // sug1 is temp, case 1, rx, rx, temp
	bp_case=1;
      else //otherwise, it is in rx (case 0), rx, rx, rx
	bp_case=0;
    }
    if(nt_c!=0){
      //MC_calc_npN_phosbase(mc_types[at_ne], r_vec, nt_ne, nt_c, nt_n); // the second argument is the owner of the phosphate
      /** HERE WE DO THE HB CHECK **/
    /* BAS[0]=rx[at_ne];    BAS[1]=ry[at_ne];    BAS[2]=rz[at_ne]; */
    /* SUG1[0]=rx[at_c+ISUG-N_PARTS_PER_NT];SUG1[1]=ry[at_c+ISUG-N_PARTS_PER_NT];SUG1[2]=rz[at_c+ISUG-N_PARTS_PER_NT]; */
    /* SUG2[0]=mc_temp_x[at_c+ISUG];SUG2[1]=mc_temp_y[at_c+ISUG];SUG2[2]=mc_temp_z[at_c+ISUG]; */
    /* PHO[0]=mc_temp_x[at_c+IPHO];PHO[1]=mc_temp_y[at_c+IPHO];PHO[2]=mc_temp_z[at_c+IPHO]; */
    
    /* MC_get_op1op2(mc_types[at_ne], MC_get_bp_face(mc_types[at_ne], r_vec), BAS, PHO, SUG1,SUG2,op1,op2); */
    /* proj_on_nt(op1, rx, ry, rz, nt_ne, bas_op1); */
    /* proj_on_nt(op2, rx, ry, rz, nt_ne, bas_op2); */
    /* if(MC_check_bph_HB(bas_op1, bas_op2, mc_types[at_ne], MC_get_bp_face(mc_types[at_ne], r_vec))==1) */
      MC_calc_npN_phosbase(mc_types[at_ne], r_vec, nt_ne, nt_c, nt_n, rx, ry, rz, bp_case); 
    
    }
  }
  /* printf("\n"); */
}

int MC_check_bph_HB(double *OP1, double *OP2, int type, int face){
  int ret=0,d;
  double CONTACT[DIM], HYDRO[DIM], va[DIM], vb1[DIM], vb2[DIM], c1, c2;
  //printf("checking type %d  face %d\n", type, face);
  CONTACT[0]=bp_contacts[type][face][0];HYDRO[0]=bp_hydros[type][face][0];
  CONTACT[1]=bp_contacts[type][face][1];HYDRO[1]=bp_hydros[type][face][1];
  CONTACT[2]=bp_contacts[type][face][2];HYDRO[2]=bp_hydros[type][face][2];
  //REMEMBER! GUANINE IS A SPECIAL CASE
  for(d=0;d<DIM;d++){
    va[d]=CONTACT[d]-HYDRO[d];
    vb1[d]=OP1[d]-HYDRO[d];
    vb2[d]=OP2[d]-HYDRO[d];
  }
  c1=dot_prod(vb1, va)/(fvec_norm(vb1)*fvec_norm(va));
  c2=dot_prod(vb2, va)/(fvec_norm(vb2)*fvec_norm(va));
  //printf("OP1: %lf  %lf\tOP2: %lf  %lf\n", fcalc_min_dist(CONTACT, OP1), c1, fcalc_min_dist(CONTACT, OP2), c2);
  if((fcalc_min_dist(CONTACT, OP1)<RHB &&  c1<CHB ) || (fcalc_min_dist(CONTACT, OP2)<RHB &&  c2<CHB ))
    ret=1;
  
  if(type==TYP_GUANINE && face==WC_FACE_WATSCRICK){
    CONTACT[0]=bp_G_cont_W2[0];HYDRO[0]=bp_G_hydr_W2[0];
    CONTACT[1]=bp_G_cont_W2[1];HYDRO[1]=bp_G_hydr_W2[1];
    CONTACT[2]=bp_G_cont_W2[2];HYDRO[2]=bp_G_hydr_W2[2];
    for(d=0;d<DIM;d++){
    va[d]=CONTACT[d]-HYDRO[d];
    vb1[d]=OP1[d]-HYDRO[d];
    vb2[d]=OP2[d]-HYDRO[d];
    }
    c1=dot_prod(vb1, va)/(fvec_norm(vb1)*fvec_norm(va));
    c2=dot_prod(vb2, va)/(fvec_norm(vb2)*fvec_norm(va));
    if((fcalc_min_dist(CONTACT, OP1)<RHB &&  c1<CHB ) || (fcalc_min_dist(CONTACT, OP2)<RHB &&  c2<CHB ))
      ret=1;
    //here, we should state that when two HB are formed, the energy is larger!
  }
  return ret;
}
//void MC_calc_npN_phosbase(int typ_a,  double *r_vec, int nt_a, int nt_b, int nt_n){
void MC_calc_npN_phosbase(int typ_a,  double *r_vec, int nt_a, int nt_b, int nt_n, double *rx, double *ry, double *rz, int bp_flag){
  //NT_A IS THE BASE, NT_B IS THE PHOSPHATE
  //classify type of pair
  int face_a =MC_get_bp_face(typ_a , r_vec);
  int at_a, at_b, at_c, hb_exists;
  double BAS[DIM], PHO[DIM], SUG1[DIM], SUG2[DIM], op1[DIM], op2[DIM], bas_op1[DIM], bas_op2[DIM];
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
  /* //calculate subfaces */
  /* int subf=MC_get_BPh_subface(typ_a, r_vec, face_a); */
  
  /* if(subf!=-1){ */
  /*   //here we modify the potential well */
  /*   well=nb_bp_spec_well[typ_ind][subf]; */
  /* } */
  //if(nt_a==6 || nt_b==6) printf("inside  %lf\n", nbe);
  if(nbe!=0 && well!=0 ){
    /*  printf("%d %d %d  %lf %lf \t %lf %lf %lf\t", nt_a, nt_b, face_a, nbe, well, r_vec[0], r_vec[1], r_vec[2]); */
    /* if(typ_a==0) printf("A "); */
    /* if(typ_a==1) printf("U "); */
    /* if(typ_a==2) printf("G "); */
    /* if(typ_a==3) printf("C "); */
    /* if(face_ind==WC_FACE_SUGAR) printf("S\t"); */
    /* if(face_ind==WC_FACE_HOOGSTEEN) printf("H\t"); */
    /* if(face_ind==WC_FACE_WATSCRICK) printf("W\t"); */
    
    /* printf("P %lf %lf %lf\tS1  %lf %lf %lf\tS2 %lf %lf %lf\tB %lf %lf %lf\n", rx[N_PARTS_PER_NT*nt_b+IPHO], ry[N_PARTS_PER_NT*nt_b+IPHO], rz[N_PARTS_PER_NT*nt_b+IPHO], */
    /* 	   rx[N_PARTS_PER_NT*(nt_b-1)+ISUG], ry[N_PARTS_PER_NT*(nt_b-1)+ISUG], rz[N_PARTS_PER_NT*(nt_b-1)+ISUG], */
    /* 	   rx[N_PARTS_PER_NT*nt_b+ISUG], ry[N_PARTS_PER_NT*nt_b+ISUG], rz[N_PARTS_PER_NT*nt_b+ISUG], */
    /* 	   rx[N_PARTS_PER_NT*nt_a+IBAS], ry[N_PARTS_PER_NT*nt_a+IBAS], rz[N_PARTS_PER_NT*nt_a+IBAS] */
    /* 	   ); */
    //HERE, OUR REFINED TEST CASE, depending on the case
    at_a=nt_a*N_PARTS_PER_NT;at_b=nt_b*N_PARTS_PER_NT;at_c=(nt_b-1)*N_PARTS_PER_NT;//remember, AT_C MUST EXIST TO PERFORM THE CHECK!
    //printf("bp flag = %d  between %d  %d  %d\n", bp_flag, nt_a, nt_b, nt_b-1);
    if(nt_b-1<0) {hb_exists=0; return;}
    if(!MC_are_neighbors(nt_b, nt_b-1)){hb_exists=0; return;}
    switch(bp_flag){
    case 0:
      BAS[0]=rx[at_a];BAS[1]=ry[at_a];BAS[2]=rz[at_a];
      PHO[0]=rx[at_b+IPHO];PHO[1]=ry[at_b+IPHO];PHO[2]=rz[at_b+IPHO];
      SUG2[0]=rx[at_b+ISUG];SUG2[1]=ry[at_b+ISUG];SUG2[2]=rz[at_b+ISUG];
      SUG1[0]=rx[at_c+ISUG];SUG1[1]=ry[at_c+ISUG];SUG1[2]=rz[at_c+ISUG];
      MC_get_op1op2(BAS, PHO, SUG1,SUG2,op1,op2); //and now we project on the base, which is in rx
      proj_on_nt(op1, rx, ry, rz, nt_a, bas_op1);
      proj_on_nt(op2, rx, ry, rz, nt_a, bas_op2);
      //printf("op1onbas   %lf %lf %lf    %lf\n", bas_op1[0], bas_op1[1], bas_op1[2], sqrt(SQ(bas_op1[0])+SQ(bas_op1[1])+SQ(bas_op1[2])));
      //printf("op2onbas   %lf %lf %lf    %lf\n", bas_op2[0], bas_op2[1], bas_op2[2], sqrt(SQ(bas_op2[0])+SQ(bas_op2[1])+SQ(bas_op2[2])));
      hb_exists=MC_check_bph_HB(bas_op1, bas_op2, typ_a, face_a);
      //to get all op1op2
      /* BAS[0]=rx[at_a];BAS[1]=ry[at_a];BAS[2]=rz[at_a]; */
      /* int ph;  */
      /* int at_ph , at_s1; */
      /* for(ph=1;ph<nt_n;ph++){ */
      /* 	at_ph=ph*N_PARTS_PER_NT;at_s1=(ph-1)*N_PARTS_PER_NT; */
      /* 	PHO[0]=rx[at_ph+IPHO];PHO[1]=ry[at_ph+IPHO];PHO[2]=rz[at_ph+IPHO]; */
      /* 	SUG2[0]=rx[at_ph+ISUG];SUG2[1]=ry[at_ph+ISUG];SUG2[2]=rz[at_ph+ISUG]; */
      /* 	SUG1[0]=rx[at_s1+ISUG];SUG1[1]=ry[at_s1+ISUG];SUG1[2]=rz[at_s1+ISUG]; */
      /* 	printf("NT %d   ", ph); */
      /* 	MC_get_op1op2(BAS, PHO, SUG1,SUG2,op1,op2); //and now we project on the base, which is in rx */
      /* } */
      /* exit(10); */
      //if(nt_b-1<0) {printf("case 0 : %lf %lf %lf    %lf %lf %lf \n", bas_op1, bas_op2);}
      break;
    case 1:
      BAS[0]=rx[at_a];BAS[1]=ry[at_a];BAS[2]=rz[at_a];
      PHO[0]=rx[at_b+IPHO];PHO[1]=ry[at_b+IPHO];PHO[2]=rz[at_b+IPHO];
      SUG2[0]=rx[at_b+ISUG];SUG2[1]=ry[at_b+ISUG];SUG2[2]=rz[at_b+ISUG];
      SUG1[0]=mc_temp_x[at_c+ISUG];SUG1[1]=mc_temp_y[at_c+ISUG];SUG1[2]=mc_temp_z[at_c+ISUG];
      MC_get_op1op2(BAS, PHO, SUG1,SUG2,op1,op2); //and now we project on the base, which is in rx
      proj_on_nt(op1, rx, ry, rz, nt_a, bas_op1);
      proj_on_nt(op2, rx, ry, rz, nt_a, bas_op2); 
      hb_exists=MC_check_bph_HB(bas_op1, bas_op2, typ_a, face_a);
      //if(nt_b-1<0) {printf("case 1 : %lf %lf %lf    %lf %lf %lf    %d\n", bas_op1, bas_op2, hb_exists);}
      break;
    case 2:
      BAS[0]=mc_temp_x[at_a];BAS[1]=mc_temp_y[at_a];BAS[2]=mc_temp_z[at_a];
      PHO[0]=rx[at_b+IPHO];PHO[1]=ry[at_b+IPHO];PHO[2]=rz[at_b+IPHO];
      SUG2[0]=rx[at_b+ISUG];SUG2[1]=ry[at_b+ISUG];SUG2[2]=rz[at_b+ISUG];
      SUG1[0]=rx[at_c+ISUG];SUG1[1]=ry[at_c+ISUG];SUG1[2]=rz[at_c+ISUG];
      MC_get_op1op2(BAS, PHO, SUG1,SUG2,op1,op2); //and now we project on the base, which is in mc_temp
      proj_on_nt(op1, mc_temp_x, mc_temp_y, mc_temp_z, nt_a, bas_op1);
      proj_on_nt(op2, mc_temp_x, mc_temp_y, mc_temp_z, nt_a, bas_op2); 
      hb_exists=MC_check_bph_HB(bas_op1, bas_op2, typ_a, face_a);
      //if(nt_b-1<0) {printf("case 2 : %lf %lf %lf    %lf %lf %lf    %d\n", bas_op1, bas_op2, hb_exists);}
      break; 
    case 3:
      BAS[0]=mc_temp_x[at_a];BAS[1]=mc_temp_y[at_a];BAS[2]=mc_temp_z[at_a];
      PHO[0]=rx[at_b+IPHO];PHO[1]=ry[at_b+IPHO];PHO[2]=rz[at_b+IPHO];
      SUG2[0]=rx[at_b+ISUG];SUG2[1]=ry[at_b+ISUG];SUG2[2]=rz[at_b+ISUG];
      SUG1[0]=mc_temp_x[at_c+ISUG];SUG1[1]=mc_temp_y[at_c+ISUG];SUG1[2]=mc_temp_z[at_c+ISUG];
      MC_get_op1op2(BAS, PHO, SUG1,SUG2,op1,op2); //and now we project on the base, which is in mc_temp
      proj_on_nt(op1, mc_temp_x, mc_temp_y, mc_temp_z, nt_a, bas_op1);
      proj_on_nt(op2, mc_temp_x, mc_temp_y, mc_temp_z, nt_a, bas_op2); 
      hb_exists=MC_check_bph_HB(bas_op1, bas_op2, typ_a, face_a);
      //if(nt_b-1<0) {printf("case 3 : %lf %lf %lf    %lf %lf %lf    %d\n", bas_op1, bas_op2, hb_exists);}
      break; 
    case 4:
      BAS[0]=rx[at_a];BAS[1]=ry[at_a];BAS[2]=rz[at_a];
      PHO[0]=mc_temp_x[at_b+IPHO];PHO[1]=mc_temp_y[at_b+IPHO];PHO[2]=mc_temp_z[at_b+IPHO];
      SUG2[0]=mc_temp_x[at_b+ISUG];SUG2[1]=mc_temp_y[at_b+ISUG];SUG2[2]=mc_temp_z[at_b+ISUG];
      SUG1[0]=rx[at_c+ISUG];SUG1[1]=ry[at_c+ISUG];SUG1[2]=rz[at_c+ISUG];
      MC_get_op1op2(BAS, PHO, SUG1,SUG2,op1,op2); //and now we project on the base, which is in rx
      proj_on_nt(op1, rx, ry, rz, nt_a, bas_op1);
      proj_on_nt(op2, rx, ry, rz, nt_a, bas_op2); 
      hb_exists=MC_check_bph_HB(bas_op1, bas_op2, typ_a, face_a);
      //if(nt_b-1<0) {printf("case 4 : %lf %lf %lf    %lf %lf %lf    %d\n", bas_op1, bas_op2, hb_exists);}
      break;
    default: hb_exists=0;
    }
    if(hb_exists==1){ 
      tot_nb_bp=nbe;
      tot_nb_bp-=well;
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
}

void MC_get_op1op2(double *BAS, double *PHO, double *SUG1, double *SUG2, double *op1, double *op2){
  double bbx[DIM], bby[DIM], bbz[DIM], vn, zdotx;
  int d;
  //we calculate the position of the contact atom in the base reference frame
  //and we calculate the positions of the OP1 and OP2 in the backbone reference frame
  for(d=0;d<DIM;d++){
    bbz[d]=SUG2[d]-PHO[d];
    bbx[d]=SUG1[d]-PHO[d];
  }
  vn=vec_norm(bbz[0], bbz[1], bbz[2]);
  for(d=0;d<DIM;d++)
    bbz[d]/=vn;
  zdotx=dot_prod(bbx, bbz);
  for(d=0;d<DIM;d++)
    bbx[d]-=zdotx*bbz[d];
  vn=vec_norm(bbx[0], bbx[1], bbx[2]);
  for(d=0;d<DIM;d++)
    bbx[d]/=vn;
  vec_prod(bbz, bbx, bby);
    
  for(d=0;d<DIM;d++){
    op1[d]=BP_OP1X*bbx[d]+BP_OP1Y*bby[d]+BP_OP1Z*bbz[d]+PHO[d]-BAS[d];
    op2[d]=BP_OP2X*bbx[d]+BP_OP2Y*bby[d]+BP_OP2Z*bbz[d]+PHO[d]-BAS[d];
  }
  //to finally project the difference in the base reference frame
  //printf("OP1 : %lf %lf %lf\nOP2: %lf %lf %lf\n", op1[0]+BAS[0], op1[1]+BAS[1], op1[2]+BAS[2], op2[0]+BAS[0], op2[1]+BAS[1], op2[2]+BAS[2]);
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
  /* else if(type==TYP_CYTOSINE && face==WC_FACE_HOOGSTEEN){ */
  /*   angle=atan2(r_vec[1], r_vec[0]); */
  /*   if(angle<-M_PI/2.0) angle+=2.0*M_PI; */
  /*   if(angle<bp_spec_ang[type][0]) ret=0; */
  /*   else{ */
  /*     if(angle<bp_spec_ang[type][1]) ret=1; */
  /*     else ret=2; */
  /*   } */
  /* } */
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

double MC_sugar_penalization(int face_c, int face_ne, int typ_c, int typ_ne, int glp1, int glp2, int flip){
  double shift=0;
  if((glp1!=0 || glp2!=0) && (face_c==WC_FACE_SUGAR || face_ne==WC_FACE_SUGAR  )){ //someone is flipped AND one of the faces involved is sugar
    double sub1=0, sub2=0, subT=0;
    if(face_c ==WC_FACE_SUGAR  && glp1!=0) sub1=1;
    if(face_ne==WC_FACE_SUGAR && glp2!=0) sub2=1;
    //case antiparallel
    if(flip==0){
      if(face_c==WC_FACE_SUGAR && face_ne==WC_FACE_SUGAR){
	//SS
	shift=0.5*(sub1+sub2)+1;
      }
      else if(face_c==WC_FACE_HOOGSTEEN && face_ne == WC_FACE_SUGAR){
	//HS
	if( !(typ_c==TYP_GUANINE && typ_ne==TYP_GUANINE) && !(typ_c==TYP_URACIL && typ_ne==TYP_ADENINE) && !(typ_c==TYP_URACIL && typ_ne==TYP_GUANINE)){
	  shift=sub2;
	}
      }
      else if(face_ne==WC_FACE_HOOGSTEEN && face_c == WC_FACE_SUGAR){
	//SH
	if( !(typ_ne==TYP_GUANINE && typ_c==TYP_GUANINE) && !(typ_ne==TYP_URACIL && typ_c==TYP_ADENINE) && !(typ_ne==TYP_URACIL && typ_c==TYP_GUANINE)){
	  shift=sub1;
	}
      }
      else if(face_c==WC_FACE_WATSCRICK && face_ne == WC_FACE_SUGAR){
	//WS
	if( !(typ_c==TYP_URACIL && typ_ne==TYP_GUANINE)){
	  shift=sub2;
	}
      }
      else if(face_ne==WC_FACE_WATSCRICK && face_c == WC_FACE_SUGAR){
	//SW
	if( !(typ_ne==TYP_URACIL && typ_c==TYP_GUANINE)){
	  shift=sub1;
	}
      }
    }
    //case parallel
    else if(flip==1){
      if(face_c==WC_FACE_SUGAR && face_ne==WC_FACE_SUGAR){
	//SS
	shift=0.5*(sub1+sub2);
	if(typ_c==TYP_GUANINE && typ_ne==TYP_GUANINE)
	  shift=shift+1;
      }
      else if(face_c==WC_FACE_WATSCRICK && face_ne == WC_FACE_SUGAR){
	//WS
	if( !(typ_c==TYP_ADENINE && typ_ne==TYP_ADENINE) && !(typ_c==TYP_URACIL && typ_ne==TYP_ADENINE)&& !(typ_c==TYP_URACIL && typ_ne==TYP_GUANINE)	    ){
	  shift=sub2;
	}
      }
      else if(face_ne==WC_FACE_WATSCRICK && face_c == WC_FACE_SUGAR){
	//SW
	if( !(typ_ne==TYP_ADENINE && typ_c==TYP_ADENINE) && !(typ_ne==TYP_URACIL && typ_c==TYP_ADENINE)&& !(typ_ne==TYP_URACIL && typ_c==TYP_GUANINE)	    ){
	  shift=sub1;
	}
      }
    }
    else{
      printf("WRONG FLIP FLAG!\n");
      exit(1);
    }
    
  }
  return shift*E_HB;
}

void MC_calc_nnN_watscric(int typ_c, int typ_ne, double *r_vec, double *r_vec_inv, double theta, double sphi, int nt_c, int nt_ne, int glp1, int glp2){
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
  int flipflag=-1;
  if(table_nnN_N_3[typ_ind]>0) 
    tnbed  =calc_nnN_tab_wc_psdihedr(theta,typ_ind,NT_UNFLIPPED); // watson-crick pseudodihedral
  if(table_nnN_N_3_F[typ_ind]>0) 
    tnbed_F=calc_nnN_tab_wc_psdihedr(theta,typ_ind,  NT_FLIPPED); // watson-crick pseudodihedral
  
  if(tnbed!=0 && tnbed_F==0){
    //CASE UNFLIPPED
    flipflag=0;
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
    flipflag=1;
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
    nb_well=nb_well-MC_sugar_penalization(face_c, face_ne, typ_c, typ_ne, glp1, glp2, flipflag);
  }
  double tot_nb_wc=0;
  //printf("PAIR %d %d    %lf %lf %lf   %lf   SHIFT %d %d\t\t%lf %lf %lf   %lf %lf %lf\n", nt_c, nt_ne, nbe, nbei, nbed, nb_well, tshift, tshift_i,r_vec[0], r_vec[1], r_vec[2], r_vec_inv[0], r_vec_inv[1], r_vec_inv[2]);
  if(nbe!=0 && nbei!=0 && nbed!=0 && nb_well!=0 ){
    tot_nb_wc=nbe+nbei+nbed;
    tot_nb_wc-=nb_well;
    // HERE WE APPLY THE SECSTRC
#ifdef SECONDSTC
    if(wc_sstruct_N[nt_c]>0 && wc_sstruct_N[nt_ne]>0){
      if(wc_sstruct_neigh[nt_c][0]==nt_ne && wc_sstruct_neigh[nt_ne][0]==nt_c && 
	 face_ind==(WC_FACES*WC_FACE_WATSCRICK + WC_FACE_WATSCRICK) && (tnbed!=0 && tnbed_F==0) ){
	//so, it makes a cWW pair which is constrained in the secstrc list
	tot_nb_wc-=1000;
      }
    }
#endif
    /**************************/
    
    if(tot_nb_wc>0) tot_nb_wc=0; //FOR SAFETY!! 
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

