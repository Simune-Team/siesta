! This file holds all the different modules for sparse vectors
! This makes editing easier as we do not need to consider *all* files

module class_sSpArr2D
  use class_sArray2D
!========================
#define TYPE_NAME  sSpArr2D
#define STR_TYPE_NAME  "sSpArr2D"
#define TYPE_NAME_ sSpArr2D_
#define NEW_TYPE newsSpArr2D
#define VAR_TYPE sArray2D
#define VAR_NEW_TYPE newsArray2D
#define VAR_TYPE_TYPE real
#define PREC sp
#include "class_SpArr2D.T90"
!========================
end module class_sSpArr2D

module class_dSpArr2D
  use class_dArray2D
!========================
#define TYPE_NAME  dSpArr2D
#define STR_TYPE_NAME  "dSpArr2D"
#define TYPE_NAME_ dSpArr2D_
#define NEW_TYPE newdSpArr2D
#define VAR_TYPE dArray2D
#define VAR_NEW_TYPE newdArray2D
#define VAR_TYPE_TYPE real
#define PREC dp
#include "class_SpArr2D.T90"
!========================
end module class_dSpArr2D

module class_cSpArr2D
  use class_cArray2D
!========================
#define TYPE_NAME  cSpArr2D
#define STR_TYPE_NAME  "cSpArr2D"
#define TYPE_NAME_ cSpArr2D_
#define NEW_TYPE newcSpArr2D
#define VAR_TYPE cArray2D
#define VAR_NEW_TYPE newcArray2D
#define VAR_TYPE_TYPE complex
#define PREC sp
#include "class_SpArr2D.T90"
!========================
end module class_cSpArr2D

module class_zSpArr2D
  use class_zArray2D
!========================
#define TYPE_NAME  zSpArr2D
#define STR_TYPE_NAME  "zSpArr2D"
#define TYPE_NAME_ zSpArr2D_
#define NEW_TYPE newzSpArr2D
#define VAR_TYPE zArray2D
#define VAR_NEW_TYPE newzArray2D
#define VAR_TYPE_TYPE complex
#define PREC dp
#include "class_SpArr2D.T90"
!========================
end module class_zSpArr2D

