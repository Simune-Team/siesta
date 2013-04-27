! We define all variables in this one
! As the entries are so short this makes more sense
module class_lArray2D
!========================
#define TYPE_NAME lArray2D
#define STR_TYPE_NAME "lArray2D"
#define TYPE_NAME_ lArray2D_
#define NEW_TYPE newlArray2D
#define VAR_TYPE logical
#define VAR_INIT .false.
#include "class_Array2D.T90"
!========================
end module class_lArray2D

module class_iArray2D
!========================
#define TYPE_NAME iArray2D
#define STR_TYPE_NAME "iArray2D"
#define TYPE_NAME_ iArray2D_
#define NEW_TYPE newiArray2D
#define VAR_TYPE integer
#define VAR_INIT 0
#include "class_Array2D.T90"
!========================
end module class_iArray2D


module class_sArray2D
!========================
#define TYPE_NAME sArray2D
#define STR_TYPE_NAME "sArray2D"
#define TYPE_NAME_ sArray2D_
#define NEW_TYPE newsArray2D
#define VAR_TYPE real
#define PREC sp
#define VAR_INIT 0._sp
#include "class_Array2D.T90"
!========================
end module class_sArray2D

module class_dArray2D
!========================
#define TYPE_NAME dArray2D
#define STR_TYPE_NAME "dArray2D"
#define TYPE_NAME_ dArray2D_
#define NEW_TYPE newdArray2D
#define VAR_TYPE real
#define PREC dp
#define VAR_INIT 0._dp
#include "class_Array2D.T90"
!========================
end module class_dArray2D

module class_cArray2D
!========================
#define TYPE_NAME cArray2D
#define STR_TYPE_NAME "cArray2D"
#define TYPE_NAME_ cArray2D_
#define NEW_TYPE newcArray2D
#define VAR_TYPE complex
#define PREC sp
#define VAR_INIT cmplx(0._sp,0._sp)
#include "class_Array2D.T90"
!========================
end module class_cArray2D

module class_zArray2D
!========================
#define TYPE_NAME zArray2D
#define STR_TYPE_NAME "zArray2D"
#define TYPE_NAME_ zArray2D_
#define NEW_TYPE newzArray2D
#define VAR_TYPE complex
#define PREC dp
#define VAR_INIT dcmplx(0._dp,0._dp)
#include "class_Array2D.T90"
!========================
end module class_zArray2D


