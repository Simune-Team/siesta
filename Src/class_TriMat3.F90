! We define all variables in this one
! As the entries are so short this makes more sense

! Currently there is no reason to generate all the modules.
! They are not used...
#ifdef __CORRECT_FOR_TRIMAT_USAGE
module class_lTriMat3
!========================
#define TYPE_NAME lTriMat3
#define STR_TYPE_NAME "lTriMat3"
#define TYPE_NAME_ lTriMat3_
#define NEW_TYPE newlTriMat3
#define VAR_TYPE logical
! DO NOT DEFINE PREC
#define VAR_INIT .false.
#include "class_TriMat3.T90"
!========================
end module class_lTriMat3

module class_iTriMat3
!========================
#define TYPE_NAME iTriMat3
#define STR_TYPE_NAME "iTriMat3"
#define TYPE_NAME_ iTriMat3_
#define NEW_TYPE newiTriMat3
#define VAR_TYPE integer
! DO NOT DEFINE PREC
#define VAR_INIT 0
#include "class_TriMat3.T90"
!========================
end module class_iTriMat3


module class_sTriMat3
!========================
#define TYPE_NAME sTriMat3
#define STR_TYPE_NAME "sTriMat3"
#define TYPE_NAME_ sTriMat3_
#define NEW_TYPE newsTriMat3
#define VAR_TYPE real
#define PREC sp
#define VAR_INIT 0._sp
#include "class_TriMat3.T90"
!========================
end module class_sTriMat3

module class_dTriMat3
!========================
#define TYPE_NAME dTriMat3
#define STR_TYPE_NAME "dTriMat3"
#define TYPE_NAME_ dTriMat3_
#define NEW_TYPE newdTriMat3
#define VAR_TYPE real
#define PREC dp
#define VAR_INIT 0._dp
#include "class_TriMat3.T90"
!========================
end module class_dTriMat3

module class_cTriMat3
!========================
#define TYPE_NAME cTriMat3
#define STR_TYPE_NAME "cTriMat3"
#define TYPE_NAME_ cTriMat3_
#define NEW_TYPE newcTriMat3
#define VAR_TYPE complex
#define PREC sp
#define VAR_INIT cmplx(0._sp,0._sp)
#include "class_TriMat3.T90"
!========================
end module class_cTriMat3

#endif

module class_zTriMat3
!========================
#define TYPE_NAME zTriMat3
#define STR_TYPE_NAME "zTriMat3"
#define TYPE_NAME_ zTriMat3_
#define NEW_TYPE newzTriMat3
#define VAR_TYPE complex
#define PREC dp
#define VAR_INIT dcmplx(0._dp,0._dp)
#include "class_TriMat3.T90"
!========================
end module class_zTriMat3


