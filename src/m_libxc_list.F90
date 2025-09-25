module m_libxc_list
!
! List of selected functionals in libxc, wrapped
! in a derived type record.
!
! Copyright (c) Alberto Garcia, 2014, 2015
!

   type, public :: libxc_kind_t
      character(len=24) :: str
   end type libxc_kind_t

   type(libxc_kind_t), parameter :: EXCH = libxc_kind_t("XC_EXCHANGE")
   type(libxc_kind_t), parameter :: CORR = libxc_kind_t("XC_CORRELATION")
   type(libxc_kind_t), parameter :: EXCH_CORR =  &
      libxc_kind_t("XC_EXCHANGE_CORRELATION")

!character(len=*), parameter :: EXCH = "XC_EXCHANGE"
!character(len=*), parameter :: CORR = "XC_CORRELATION"
!character(len=*), parameter :: EXCH_CORR = "XC_EXCHANGE_CORRELATION"

   type, public :: libxc_t
      character(len=40)  :: name
      integer            :: code
      type(libxc_kind_t) :: xc_kind
!   character(len=40)  :: xc_kind
   end type libxc_t

   type(libxc_t), parameter, public  ::                  &
      XC_EMPTY = libxc_t("XC_EMPTY", 0 , EXCH_CORR),  &
      XC_NOT_IMPL = libxc_t("XC_NOT_IMPL", -1 , EXCH_CORR)

! Exchange

   type(libxc_t), parameter  ::                  &
      XC_LDA_X = libxc_t("XC_LDA_X", 1 , EXCH),  &
      XC_LDA_X_2D = libxc_t("XC_LDA_X_2D", 19 , EXCH),  &
      XC_LDA_X_1D = libxc_t("XC_LDA_X_1D", 21 , EXCH),  &
      XC_GGA_X_Q2D = libxc_t("XC_GGA_X_Q2D", 48 , EXCH),  &
      XC_GGA_X_PBE_MOL = libxc_t("XC_GGA_X_PBE_MOL", 49 , EXCH),  &
      XC_GGA_X_AK13 = libxc_t("XC_GGA_X_AK13", 56 , EXCH),  &
      XC_GGA_X_LV_RPW86 = libxc_t("XC_GGA_X_LV_RPW86", 58 , EXCH),  &
      XC_GGA_X_PBE_TCA = libxc_t("XC_GGA_X_PBE_TCA", 59 , EXCH),  &
      XC_GGA_X_PBEINT = libxc_t("XC_GGA_X_PBEINT", 60 , EXCH),  &
      XC_GGA_X_VMT84_GE = libxc_t("XC_GGA_X_VMT84_GE", 68 , EXCH),  &
      XC_GGA_X_VMT84_PBE = libxc_t("XC_GGA_X_VMT84_PBE", 69 , EXCH),  &
      XC_GGA_X_VMT_GE = libxc_t("XC_GGA_X_VMT_GE", 70 , EXCH),  &
      XC_GGA_X_VMT_PBE = libxc_t("XC_GGA_X_VMT_PBE", 71 , EXCH),  &
      XC_GGA_X_N12 = libxc_t("XC_GGA_X_N12", 82 , EXCH),  &
      XC_GGA_X_SSB_SW = libxc_t("XC_GGA_X_SSB_SW", 90 , EXCH),  &
      XC_GGA_X_SSB = libxc_t("XC_GGA_X_SSB", 91 , EXCH),  &
      XC_GGA_X_SSB_D = libxc_t("XC_GGA_X_SSB_D", 92 , EXCH),  &
      XC_GGA_X_BPCCAC = libxc_t("XC_GGA_X_BPCCAC", 98 , EXCH),  &
      XC_GGA_X_PBE = libxc_t("XC_GGA_X_PBE", 101 , EXCH),  &
      XC_GGA_X_PBE_R = libxc_t("XC_GGA_X_PBE_R", 102 , EXCH),  &
      XC_GGA_X_B86 = libxc_t("XC_GGA_X_B86", 103 , EXCH),  &
      XC_GGA_X_HERMAN = libxc_t("XC_GGA_X_HERMAN", 104 , EXCH),  &
      XC_GGA_X_B86_MGC = libxc_t("XC_GGA_X_B86_MGC", 105 , EXCH),  &
      XC_GGA_X_B88 = libxc_t("XC_GGA_X_B88", 106 , EXCH),  &
      XC_GGA_X_G96 = libxc_t("XC_GGA_X_G96", 107 , EXCH),  &
      XC_GGA_X_PW86 = libxc_t("XC_GGA_X_PW86", 108 , EXCH),  &
      XC_GGA_X_PW91 = libxc_t("XC_GGA_X_PW91", 109 , EXCH),  &
      XC_GGA_X_OPTX = libxc_t("XC_GGA_X_OPTX", 110 , EXCH),  &
      XC_GGA_X_DK87_R1 = libxc_t("XC_GGA_X_DK87_R1", 111 , EXCH),  &
      XC_GGA_X_DK87_R2 = libxc_t("XC_GGA_X_DK87_R2", 112 , EXCH),  &
      XC_GGA_X_LG93 = libxc_t("XC_GGA_X_LG93", 113 , EXCH),  &
      XC_GGA_X_FT97_A = libxc_t("XC_GGA_X_FT97_A", 114 , EXCH),  &
      XC_GGA_X_FT97_B = libxc_t("XC_GGA_X_FT97_B", 115 , EXCH),  &
      XC_GGA_X_PBE_SOL = libxc_t("XC_GGA_X_PBE_SOL", 116 , EXCH),  &
      XC_GGA_X_RPBE = libxc_t("XC_GGA_X_RPBE", 117 , EXCH),  &
      XC_GGA_X_WC = libxc_t("XC_GGA_X_WC", 118 , EXCH),  &
      XC_GGA_X_MPW91 = libxc_t("XC_GGA_X_MPW91", 119 , EXCH),  &
      XC_GGA_X_AM05 = libxc_t("XC_GGA_X_AM05", 120 , EXCH),  &
      XC_GGA_X_PBEA = libxc_t("XC_GGA_X_PBEA", 121 , EXCH),  &
      XC_GGA_X_MPBE = libxc_t("XC_GGA_X_MPBE", 122 , EXCH),  &
      XC_GGA_X_XPBE = libxc_t("XC_GGA_X_XPBE", 123 , EXCH),  &
      XC_GGA_X_2D_B86_MGC = libxc_t("XC_GGA_X_2D_B86_MGC", 124 , EXCH),  &
      XC_GGA_X_BAYESIAN = libxc_t("XC_GGA_X_BAYESIAN", 125 , EXCH),  &
      XC_GGA_X_PBE_JSJR = libxc_t("XC_GGA_X_PBE_JSJR", 126 , EXCH),  &
      XC_GGA_X_2D_B88 = libxc_t("XC_GGA_X_2D_B88", 127 , EXCH),  &
      XC_GGA_X_2D_B86 = libxc_t("XC_GGA_X_2D_B86", 128 , EXCH),  &
      XC_GGA_X_2D_PBE = libxc_t("XC_GGA_X_2D_PBE", 129 , EXCH),  &
      XC_GGA_X_OPTB88_VDW = libxc_t("XC_GGA_X_OPTB88_VDW", 139 , EXCH),  &
      XC_GGA_X_PBEK1_VDW = libxc_t("XC_GGA_X_PBEK1_VDW", 140 , EXCH),  &
      XC_GGA_X_OPTPBE_VDW = libxc_t("XC_GGA_X_OPTPBE_VDW", 141 , EXCH),  &
      XC_GGA_X_RGE2 = libxc_t("XC_GGA_X_RGE2", 142 , EXCH),  &
      XC_GGA_X_RPW86 = libxc_t("XC_GGA_X_RPW86", 144 , EXCH),  &
      XC_GGA_X_KT1 = libxc_t("XC_GGA_X_KT1", 145 , EXCH),  &
      XC_GGA_X_MB88 = libxc_t("XC_GGA_X_MB88", 149 , EXCH),  &
      XC_GGA_X_SOGGA = libxc_t("XC_GGA_X_SOGGA", 150 , EXCH),  &
      XC_GGA_X_SOGGA11 = libxc_t("XC_GGA_X_SOGGA11", 151 , EXCH),  &
      XC_GGA_X_C09X = libxc_t("XC_GGA_X_C09X", 158 , EXCH),  &
      XC_GGA_X_LB = libxc_t("XC_GGA_X_LB", 160 , EXCH),  &
      XC_GGA_X_LBM = libxc_t("XC_GGA_X_LBM", 182 , EXCH),  &
      XC_GGA_X_OL2 = libxc_t("XC_GGA_X_OL2", 183 , EXCH),  &
      XC_GGA_X_APBE = libxc_t("XC_GGA_X_APBE", 184 , EXCH),  &
      XC_GGA_X_HTBS = libxc_t("XC_GGA_X_HTBS", 191 , EXCH),  &
      XC_GGA_X_AIRY = libxc_t("XC_GGA_X_AIRY", 192 , EXCH),  &
      XC_GGA_X_LAG = libxc_t("XC_GGA_X_LAG", 193 , EXCH),  &
      XC_GGA_X_WPBEH = libxc_t("XC_GGA_X_WPBEH", 524 , EXCH),  &
      XC_GGA_X_HJS_PBE = libxc_t("XC_GGA_X_HJS_PBE", 525 , EXCH),  &
      XC_GGA_X_HJS_PBE_SOL = libxc_t("XC_GGA_X_HJS_PBE_SOL", 526 , EXCH),  &
      XC_GGA_X_HJS_B88 = libxc_t("XC_GGA_X_HJS_B88", 527 , EXCH),  &
      XC_GGA_X_HJS_B97X = libxc_t("XC_GGA_X_HJS_B97X", 528 , EXCH),  &
      XC_GGA_X_ITYH = libxc_t("XC_GGA_X_ITYH", 529 , EXCH),  &
      XC_GGA_X_SFAT = libxc_t("XC_GGA_X_SFAT", 530 , EXCH)

! Correlation

   type(libxc_t), parameter  ::                  &
      XC_LDA_C_WIGNER = libxc_t("XC_LDA_C_WIGNER", 2 , CORR),  &
      XC_LDA_C_RPA = libxc_t("XC_LDA_C_RPA", 3 , CORR),  &
      XC_LDA_C_HL = libxc_t("XC_LDA_C_HL", 4 , CORR),  &
      XC_LDA_C_GL = libxc_t("XC_LDA_C_GL", 5 , CORR),  &
      XC_LDA_C_XALPHA = libxc_t("XC_LDA_C_XALPHA", 6 , CORR),  &
      XC_LDA_C_VWN = libxc_t("XC_LDA_C_VWN", 7 , CORR),  &
      XC_LDA_C_VWN_RPA = libxc_t("XC_LDA_C_VWN_RPA", 8 , CORR),  &
      XC_LDA_C_PZ = libxc_t("XC_LDA_C_PZ", 9 , CORR),  &
      XC_LDA_C_PZ_MOD = libxc_t("XC_LDA_C_PZ_MOD", 10 , CORR),  &
      XC_LDA_C_OB_PZ = libxc_t("XC_LDA_C_OB_PZ", 11 , CORR),  &
      XC_LDA_C_PW = libxc_t("XC_LDA_C_PW", 12 , CORR),  &
      XC_LDA_C_PW_MOD = libxc_t("XC_LDA_C_PW_MOD", 13 , CORR),  &
      XC_LDA_C_OB_PW = libxc_t("XC_LDA_C_OB_PW", 14 , CORR),  &
      XC_LDA_C_2D_AMGB = libxc_t("XC_LDA_C_2D_AMGB", 15 , CORR),  &
      XC_LDA_C_2D_PRM = libxc_t("XC_LDA_C_2D_PRM", 16 , CORR),  &
      XC_LDA_C_vBH = libxc_t("XC_LDA_C_vBH", 17 , CORR),  &
      XC_LDA_C_1D_CSC = libxc_t("XC_LDA_C_1D_CSC", 18 , CORR),  &
      XC_LDA_C_ML1 = libxc_t("XC_LDA_C_ML1", 22 , CORR),  &
      XC_LDA_C_ML2 = libxc_t("XC_LDA_C_ML2", 23 , CORR),  &
      XC_LDA_C_GOMBAS = libxc_t("XC_LDA_C_GOMBAS", 24 , CORR),  &
      XC_LDA_C_PW_RPA = libxc_t("XC_LDA_C_PW_RPA", 25 , CORR),  &
      XC_LDA_C_1D_LOOS = libxc_t("XC_LDA_C_1D_LOOS", 26 , CORR),  &
      XC_LDA_C_RC04 = libxc_t("XC_LDA_C_RC04", 27 , CORR),  &
      XC_LDA_C_VWN_1 = libxc_t("XC_LDA_C_VWN_1", 28 , CORR),  &
      XC_LDA_C_VWN_2 = libxc_t("XC_LDA_C_VWN_2", 29 , CORR),  &
      XC_LDA_C_VWN_3 = libxc_t("XC_LDA_C_VWN_3", 30 , CORR),  &
      XC_LDA_C_VWN_4 = libxc_t("XC_LDA_C_VWN_4", 31 , CORR),  &
      XC_GGA_C_Q2D = libxc_t("XC_GGA_C_Q2D", 47 , CORR),  &
      XC_GGA_C_ZPBEINT = libxc_t("XC_GGA_C_ZPBEINT", 61 , CORR),  &
      XC_GGA_C_PBEINT = libxc_t("XC_GGA_C_PBEINT", 62 , CORR),  &
      XC_GGA_C_ZPBESOL = libxc_t("XC_GGA_C_ZPBESOL", 63 , CORR),  &
      XC_GGA_C_N12_SX = libxc_t("XC_GGA_C_N12_SX", 79 , CORR),  &
      XC_GGA_C_N12 = libxc_t("XC_GGA_C_N12", 80 , CORR),  &
      XC_GGA_C_VPBE = libxc_t("XC_GGA_C_VPBE", 83 , CORR),  &
      XC_GGA_C_OP_XALPHA = libxc_t("XC_GGA_C_OP_XALPHA", 84 , CORR),  &
      XC_GGA_C_OP_G96 = libxc_t("XC_GGA_C_OP_G96", 85 , CORR),  &
      XC_GGA_C_OP_PBE = libxc_t("XC_GGA_C_OP_PBE", 86 , CORR),  &
      XC_GGA_C_OP_B88 = libxc_t("XC_GGA_C_OP_B88", 87 , CORR),  &
      XC_GGA_C_FT97 = libxc_t("XC_GGA_C_FT97", 88 , CORR),  &
      XC_GGA_C_SPBE = libxc_t("XC_GGA_C_SPBE", 89 , CORR),  &
      XC_GGA_C_REVTCA = libxc_t("XC_GGA_C_REVTCA", 99 , CORR),  &
      XC_GGA_C_TCA = libxc_t("XC_GGA_C_TCA", 100 , CORR),  &
      XC_GGA_C_PBE = libxc_t("XC_GGA_C_PBE", 130 , CORR),  &
      XC_GGA_C_LYP = libxc_t("XC_GGA_C_LYP", 131 , CORR),  &
      XC_GGA_C_P86 = libxc_t("XC_GGA_C_P86", 132 , CORR),  &
      XC_GGA_C_PBE_SOL = libxc_t("XC_GGA_C_PBE_SOL", 133 , CORR),  &
      XC_GGA_C_PW91 = libxc_t("XC_GGA_C_PW91", 134 , CORR),  &
      XC_GGA_C_AM05 = libxc_t("XC_GGA_C_AM05", 135 , CORR),  &
      XC_GGA_C_XPBE = libxc_t("XC_GGA_C_XPBE", 136 , CORR),  &
      XC_GGA_C_LM = libxc_t("XC_GGA_C_LM", 137 , CORR),  &
      XC_GGA_C_PBE_JRGX = libxc_t("XC_GGA_C_PBE_JRGX", 138 , CORR),  &
      XC_GGA_C_RGE2 = libxc_t("XC_GGA_C_RGE2", 143 , CORR),  &
      XC_GGA_C_WL = libxc_t("XC_GGA_C_WL", 147 , CORR),  &
      XC_GGA_C_WI = libxc_t("XC_GGA_C_WI", 148 , CORR),  &
      XC_GGA_C_SOGGA11 = libxc_t("XC_GGA_C_SOGGA11", 152 , CORR),  &
      XC_GGA_C_WI0 = libxc_t("XC_GGA_C_WI0", 153 , CORR),  &
      XC_GGA_C_SOGGA11_X = libxc_t("XC_GGA_C_SOGGA11_X", 159 , CORR),  &
      XC_GGA_C_APBE = libxc_t("XC_GGA_C_APBE", 186 , CORR),  &
      XC_GGA_C_OPTC = libxc_t("XC_GGA_C_OPTC", 200 , CORR)

! Exchange and correlation
   type(libxc_t), parameter  ::                  &
      XC_LDA_XC_TETER93 = libxc_t("XC_LDA_XC_TETER93", 20 , EXCH_CORR),  &
      XC_GGA_XC_OPBE_D = libxc_t("XC_GGA_XC_OPBE_D", 65 , EXCH_CORR),  &
      XC_GGA_XC_OPWLYP_D = libxc_t("XC_GGA_XC_OPWLYP_D", 66 , EXCH_CORR),  &
      XC_GGA_XC_OBLYP_D = libxc_t("XC_GGA_XC_OBLYP_D", 67 , EXCH_CORR),  &
      XC_GGA_XC_HCTH_407P = libxc_t("XC_GGA_XC_HCTH_407P", 93 , EXCH_CORR),  &
      XC_GGA_XC_HCTH_P76 = libxc_t("XC_GGA_XC_HCTH_P76", 94 , EXCH_CORR),  &
      XC_GGA_XC_HCTH_P14 = libxc_t("XC_GGA_XC_HCTH_P14", 95 , EXCH_CORR),  &
      XC_GGA_XC_B97_GGA1 = libxc_t("XC_GGA_XC_B97_GGA1", 96 , EXCH_CORR),  &
      XC_GGA_XC_HCTH_A = libxc_t("XC_GGA_XC_HCTH_A", 97 , EXCH_CORR),  &
      XC_GGA_XC_KT2 = libxc_t("XC_GGA_XC_KT2", 146 , EXCH_CORR),  &
      XC_GGA_XC_TH1 = libxc_t("XC_GGA_XC_TH1", 154 , EXCH_CORR),  &
      XC_GGA_XC_TH2 = libxc_t("XC_GGA_XC_TH2", 155 , EXCH_CORR),  &
      XC_GGA_XC_TH3 = libxc_t("XC_GGA_XC_TH3", 156 , EXCH_CORR),  &
      XC_GGA_XC_TH4 = libxc_t("XC_GGA_XC_TH4", 157 , EXCH_CORR),  &
      XC_GGA_XC_HCTH_93 = libxc_t("XC_GGA_XC_HCTH_93", 161 , EXCH_CORR),  &
      XC_GGA_XC_HCTH_120 = libxc_t("XC_GGA_XC_HCTH_120", 162 , EXCH_CORR),  &
      XC_GGA_XC_HCTH_147 = libxc_t("XC_GGA_XC_HCTH_147", 163 , EXCH_CORR),  &
      XC_GGA_XC_HCTH_407 = libxc_t("XC_GGA_XC_HCTH_407", 164 , EXCH_CORR),  &
      XC_GGA_XC_EDF1 = libxc_t("XC_GGA_XC_EDF1", 165 , EXCH_CORR),  &
      XC_GGA_XC_XLYP = libxc_t("XC_GGA_XC_XLYP", 166 , EXCH_CORR),  &
      XC_GGA_XC_B97 = libxc_t("XC_GGA_XC_B97", 167 , EXCH_CORR),  &
      XC_GGA_XC_B97_1 = libxc_t("XC_GGA_XC_B97_1", 168 , EXCH_CORR),  &
      XC_GGA_XC_B97_2 = libxc_t("XC_GGA_XC_B97_2", 169 , EXCH_CORR),  &
      XC_GGA_XC_B97_D = libxc_t("XC_GGA_XC_B97_D", 170 , EXCH_CORR),  &
      XC_GGA_XC_B97_K = libxc_t("XC_GGA_XC_B97_K", 171 , EXCH_CORR),  &
      XC_GGA_XC_B97_3 = libxc_t("XC_GGA_XC_B97_3", 172 , EXCH_CORR),  &
      XC_GGA_XC_PBE1W = libxc_t("XC_GGA_XC_PBE1W", 173 , EXCH_CORR),  &
      XC_GGA_XC_MPWLYP1W = libxc_t("XC_GGA_XC_MPWLYP1W", 174 , EXCH_CORR),  &
      XC_GGA_XC_PBELYP1W = libxc_t("XC_GGA_XC_PBELYP1W", 175 , EXCH_CORR),  &
      XC_GGA_XC_SB98_1a = libxc_t("XC_GGA_XC_SB98_1a", 176 , EXCH_CORR),  &
      XC_GGA_XC_SB98_1b = libxc_t("XC_GGA_XC_SB98_1b", 177 , EXCH_CORR),  &
      XC_GGA_XC_SB98_1c = libxc_t("XC_GGA_XC_SB98_1c", 178 , EXCH_CORR),  &
      XC_GGA_XC_SB98_2a = libxc_t("XC_GGA_XC_SB98_2a", 179 , EXCH_CORR),  &
      XC_GGA_XC_SB98_2b = libxc_t("XC_GGA_XC_SB98_2b", 180 , EXCH_CORR),  &
      XC_GGA_XC_SB98_2c = libxc_t("XC_GGA_XC_SB98_2c", 181 , EXCH_CORR),  &
      XC_GGA_XC_MOHLYP = libxc_t("XC_GGA_XC_MOHLYP", 194 , EXCH_CORR),  &
      XC_GGA_XC_MOHLYP2 = libxc_t("XC_GGA_XC_MOHLYP2", 195 , EXCH_CORR),  &
      XC_GGA_XC_TH_FL = libxc_t("XC_GGA_XC_TH_FL", 196 , EXCH_CORR),  &
      XC_GGA_XC_TH_FC = libxc_t("XC_GGA_XC_TH_FC", 197 , EXCH_CORR),  &
      XC_GGA_XC_TH_FCFO = libxc_t("XC_GGA_XC_TH_FCFO", 198 , EXCH_CORR),  &
      XC_GGA_XC_TH_FCO = libxc_t("XC_GGA_XC_TH_FCO", 199 , EXCH_CORR)

end module m_libxc_list
