module m_siestaxc_list

private

type, public :: siestaxc_t
   character(len=10)  :: family
   character(len=20)  :: authors
!   integer            :: code
end type siestaxc_t

! LDA

type(siestaxc_t), public, parameter  ::                  &
         SXC_LDA_PZ   = siestaxc_t("LDA", "PZ"),   &
         SXC_LDA_CA   = siestaxc_t("LDA", "CA"),   &
         SXC_LDA_PW92 = siestaxc_t("LDA", "PW92"),  &
         SXC_LDA_WIGNER = siestaxc_t("LDA", "Wigner"), &
         SXC_LDA_HL     = siestaxc_t("LDA", "Hedin-Lundqvist"), &
         SXC_LDA_GL     = siestaxc_t("LDA", "Gunnarson-Lundqvist"), &
         SXC_LDA_VBH    = siestaxc_t("LDA", "von Barth-Hedin")

! GGA
type(siestaxc_t), public, parameter  ::                  &
         SXC_GGA_PW91   = siestaxc_t("GGA", "PW91"), &
         SXC_GGA_PBE    = siestaxc_t("GGA", "PBE"),  &
!"RPBE - Hammer et al"
         SXC_GGA_RPBE   = siestaxc_t("GGA", "RPBE"), &
!"revPBE Zhang+Yang"
         SXC_GGA_REVPBE = siestaxc_t("GGA", "revPBE"), &
!"Becke-Lee-Yang-Parr"
         SXC_GGA_LYP    = siestaxc_t("GGA", "LYP"),    &
!"Wu-Cohen"
         SXC_GGA_WC     = siestaxc_t("GGA", "WC"),     &
!"Perdew-Burke-Ernzerhof-solid"
         SXC_GGA_PBESOL = siestaxc_t("GGA", "PBEsol"), &
!
         SXC_GGA_PBEJsJrLO  = siestaxc_t("GGA", "PBEJsJrLO"), &
         SXC_GGA_PBEJsJrHEG = siestaxc_t("GGA", "PBEJsJrHEG"), &
         SXC_GGA_PBEGcGxLO  = siestaxc_t("GGA", "PBEGcGxLO"),  &
         SXC_GGA_PBEGcGxHEG = siestaxc_t("GGA", "PBEGcGxHEG"), &
!"Armiento-Mattsson-05"
         SXC_GGA_AM05       = siestaxc_t("GGA", "AM05")

! VDW
type(siestaxc_t), public, parameter  ::                  &
         SXC_VDW_DRSLL   = siestaxc_t("VDW", "DRSLL"), &
         SXC_VDW_LMKLL    = siestaxc_t("VDW", "LMKLL"),  &
         SXC_VDW_KKBM    = siestaxc_t("VDW", "KKBM"),  &
         SXC_VDW_C09    = siestaxc_t("VDW", "C09"),  &
         SXC_VDW_BH    = siestaxc_t("VDW", "BH"),  &
         SXC_VDW_VV    = siestaxc_t("VDW", "VV")

end module m_siestaxc_list
