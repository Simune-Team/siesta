! We define all variables in this one
! As the entries are so short this makes more sense

module class_lArray1D
!========================
#define TYPE_NAME lArray1D
#define STR_TYPE_NAME "lArray1D"
#define TYPE_NAME_ lArray1D_
#define NEW_TYPE newlArray1D
#define VAR_TYPE logical
! DO NOT DEFINE PREC
#include "class_Array1D.T90"
!========================
end module class_lArray1D

module class_iArray1D
!========================
#define TYPE_NAME iArray1D
#define STR_TYPE_NAME "iArray1D"
#define TYPE_NAME_ iArray1D_
#define NEW_TYPE newiArray1D
#define VAR_TYPE integer
! DO NOT DEFINE PREC
#include "class_Array1D.T90"
!========================
end module class_iArray1D


module class_sArray1D
!========================
#define TYPE_NAME sArray1D
#define STR_TYPE_NAME "sArray1D"
#define TYPE_NAME_ sArray1D_
#define NEW_TYPE newsArray1D
#define VAR_TYPE real
#define PREC sp
#include "class_Array1D.T90"
!========================
end module class_sArray1D

module class_dArray1D
!========================
#define TYPE_NAME dArray1D
#define STR_TYPE_NAME "dArray1D"
#define TYPE_NAME_ dArray1D_
#define NEW_TYPE newdArray1D
#define VAR_TYPE real
#define PREC dp
#include "class_Array1D.T90"
!========================
end module class_dArray1D

module class_cArray1D
!========================
#define TYPE_NAME cArray1D
#define STR_TYPE_NAME "cArray1D"
#define TYPE_NAME_ cArray1D_
#define NEW_TYPE newcArray1D
#define VAR_TYPE complex
#define PREC sp
#include "class_Array1D.T90"
!========================
end module class_cArray1D

module class_zArray1D
!========================
#define TYPE_NAME zArray1D
#define STR_TYPE_NAME "zArray1D"
#define TYPE_NAME_ zArray1D_
#define NEW_TYPE newzArray1D
#define VAR_TYPE complex
#define PREC dp
#include "class_Array1D.T90"
!========================
end module class_zArray1D


