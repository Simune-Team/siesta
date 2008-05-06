MODULE atom_options

  implicit none
  PUBLIC

  logical :: write_ion_plot_files    ! Write small auxiliary files?

CONTAINS

  subroutine get_atom_options()
    use fdf
    write_ion_plot_files = fdf_boolean('WriteIonPlotFiles',.false.)
  end subroutine get_atom_options

END MODULE atom_options
