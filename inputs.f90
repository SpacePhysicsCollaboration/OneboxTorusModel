MODULE INPUTS

  USE GLOBAL
!a module to read and store inputs from input.dat
!should only be used in onebox.f90
!variables will be global to eliminate the need for passing variables

real                  ::source, o_to_s, transport, tehot
real                  ::neutral_amp, neutral_t0, neutral_width
real                  ::hote_amp, hote_t0, hote_width

CONTAINS

subroutine readInputs()
  open(unit=100, file='inputs.dat', status='old')

  read(100,*) source
  source=source*ROOTPI*.5*Rj*1.0e5
  print *, Rj, ROOTPI, source
  read(100,*) o_to_s
  read(100,*) fehot_const  !declared in global module
  read(100,*) transport
  read(100,*) tehot
  read(100,*) lag_const  !declared in global module
  read(100,*) neutral_amp
  read(100,*) neutral_t0
  read(100,*) neutral_width
  read(100,*) hote_amp
  read(100,*) hote_t0
  read(100,*) hote_width

  close(100)
end subroutine readInputs

END MODULE
