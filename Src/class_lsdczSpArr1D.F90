! This file holds all the different modules for sparse vectors
! This makes editing easier as we do not need to consider *all* files

! A logical sparse array
module class_lSpArr1D
  use class_lArray1D
!========================
#define TYPE_NAME  lSpArr1D
#define STR_TYPE_NAME "lSpArr1D"
#define TYPE_NAME_ lSpArr1D_
#define NEW_TYPE newlSpArr1D
#define VAR_TYPE lArray1D
#define VAR_NEW_TYPE newlArray1D
#define VAR_TYPE_TYPE logical
! Do not define precision
#include "class_SpArr1D.T90"
!========================
end module class_lSpArr1D

module class_iSpArr1D
  use class_iArray1D
!========================
#define TYPE_NAME  iSpArr1D
#define STR_TYPE_NAME "iSpArr1D"
#define TYPE_NAME_ iSpArr1D_
#define NEW_TYPE newiSpArr1D
#define VAR_TYPE iArray1D
#define VAR_NEW_TYPE newiArray1D
#define VAR_TYPE_TYPE integer
! Do not define precision
#include "class_SpArr1D.T90"
!========================
end module class_iSpArr1D

module class_sSpArr1D
  use class_sArray1D
!========================
#define TYPE_NAME  sSpArr1D
#define STR_TYPE_NAME "sSpArr1D"
#define TYPE_NAME_ sSpArr1D_
#define NEW_TYPE newsSpArr1D
#define VAR_TYPE sArray1D
#define VAR_NEW_TYPE newsArray1D
#define VAR_TYPE_TYPE real
#define PREC sp
#include "class_SpArr1D.T90"
!========================
end module class_sSpArr1D

module class_dSpArr1D
  use class_dArray1D
!========================
#define TYPE_NAME  dSpArr1D
#define STR_TYPE_NAME "dSpArr1D"
#define TYPE_NAME_ dSpArr1D_
#define NEW_TYPE newdSpArr1D
#define VAR_TYPE dArray1D
#define VAR_NEW_TYPE newdArray1D
#define VAR_TYPE_TYPE real
#define PREC dp
#include "class_SpArr1D.T90"
!========================
end module class_dSpArr1D

module class_cSpArr1D
  use class_cArray1D
!========================
#define TYPE_NAME  cSpArr1D
#define STR_TYPE_NAME "cSpArr1D"
#define TYPE_NAME_ cSpArr1D_
#define NEW_TYPE newcSpArr1D
#define VAR_TYPE cArray1D
#define VAR_NEW_TYPE newcArray1D
#define VAR_TYPE_TYPE complex
#define PREC sp
#include "class_SpArr1D.T90"
!========================
end module class_cSpArr1D

module class_zSpArr1D
  use class_zArray1D
!========================
#define TYPE_NAME  zSpArr1D
#define STR_TYPE_NAME "zSpArr1D"
#define TYPE_NAME_ zSpArr1D_
#define NEW_TYPE newzSpArr1D
#define VAR_TYPE zArray1D
#define VAR_NEW_TYPE newzArray1D
#define VAR_TYPE_TYPE complex
#define PREC dp
#include "class_SpArr1D.T90"
!========================
end module class_zSpArr1D

