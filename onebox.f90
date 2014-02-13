!#########################################
!           Matthew Copper
! One box model of io plasma torus
! Based on IDL program by Andrew Steffl 
! Started on 2/25/2013
!#########################################

PROGRAM Onebox

  USE DEBUG
  USE FUNCTIONS
  USE TIMESTEP 
  USE ReadEmis
  USE INPUTS
  USE FTMIX

  IMPLICIT NONE

  call model()

CONTAINS 

subroutine model()
  integer             ::nit
  real                ::lontemp  
  real                ::tm, tm0
  type(density)       ::n, ni, np 
  real                ::const
  type(temp)          ::T, Ti, Tp
  real                ::Te0, Ti0, Teh0
  type(height)        ::h, hi
  type(r_ind)         ::ind
  character(len=20)   ::loc
  type(nT)            ::nrg, nTi, nTp
  integer             ::i, j
  real                ::var, n_height
  type(nu)            ::v, vi
  type(r_dep)         ::dep, depi
  type(lat_dist)      ::lat, lati
  type(ft_int)        ::ft
  type(energy)        ::nrgy
  type(ft_mix)        ::plot
  real                ::output_it

!  call initNu(v)

  do i=1, LAT_SIZE
    lat%z(i)= (i-1) * Rj / 5.0  !Initializing lat%z
  end do

  call readInputs()  !call to input.f90 module to read initial variables from 'input.dat'
  call read_rec_tables()

!set trans_ variables (user prompt or formatted file migh be used in the future)
  trans_exp=1.0
  trans_type=.false.

!set dt (2000)
  dt=2000.0

!set run time
  runt=500.0*8.64e4 !500 days in second

  nit=(runt/dt)+1 ! number of iterations to reach 500 days

!set radial distance
  rdist= 6   !in Rj

!set sys3 longitude of box
  lon3=110.0
  lontemp=lon3

!set zoff
  zoff= 6.4 * cos((lon3-200.0) * dTOr) * dTOr * rdist * Rj !in km 
  if (zoff<0) then
    zoff= zoff*(-1.0)
  end if

  n_height = Rj/2

  tm0=0.01

!set density values
  const=1800.0
  n%sp = 0.060 * const
  n%s2p= 0.212 * const
  n%s3p= 0.034 * const
  n%op = 0.242 * const
  n%o2p= 0.123 * n%op

  n%s=25.0
  n%o=50.0

  Te0 = 5.0
  Ti0 = 70.0
  Teh0= tehot!49.0
!  fehot_const= .0022
  trans = 4.62963e-7
  net_source = source ! ~6.3e6

  n%elec = (n%sp + n%op + 2 * (n%s2p + n%o2p) + 3 * n%s3p) * (1.0 - n%protons)
  n%elecHot = n%fh * n%elec / (1.0-n%fh)
  n%fc = 1.0 - n%fh

!set temp values
  T%sp      = Ti0
  T%s2p     = Ti0
  T%s3p     = Ti0
  T%op      = Ti0
  T%o2p     = Ti0
  T%elec    = Te0
  T%elecHot = Teh0

!get scale heights 
  call get_scale_heights(h, T, n)

  if (protons > 0.0) then
    n%protons = protons
  endif

  Teh0 = tehot
  ind%o_to_s= o_to_s !fix make this read from file 
  ind%o2s_spike=2.0

  tau0=transport
  net_source0=net_source
  !fh0 = fehot_const

  h%s=n_height
  h%o=n_height

  call ReadIndices(ind%emis_temp, ind%emis_dens) !reads in the temp and density tables for inerpolate the emmission tables
! The above function should be unneccesary since values are now calculated rather than searched in a table. Calculation is faster thant searching.

!Read in all emission tables. These tables are used to determine power radiated by the torus
  loc='emisSp.dat'
  call ReadEmisTable(loc, ind%emis_sp)
  loc='emisS2p.dat'
  call ReadEmisTable(loc, ind%emis_s2p)
  loc='emisS3p.dat'
  call ReadEmisTable(loc, ind%emis_s3p)
  loc='emisOp.dat'
  call ReadEmisTable(loc, ind%emis_op)
  loc='emisO2p.dat'
  call ReadEmisTable(loc, ind%emis_o2p)

  T%pu_s = Tpu(32.0, rdist*1.0)
  T%pu_o = Tpu(16.0, rdist*1.0)

  T%elecHot=Teh0
 print *, T%elecHot 

  call independent_rates(ind, T, h)

  net_source0=net_source

  n%fc= 1.0 - n%fh   

  n%elec = ((n%sp + n%op) + 2.0*(n%s2p + n%o2p) + 3.0 * n%s3p)/(1.0-n%protons)
  print *, n%elec
  n%elecHot = n%elec * n%fh/n%fc

  nrg%elec = n%elec * T%elec
  nrg%elecHot = n%elecHot * T%elecHot
  nrg%sp = n%sp * T%sp
  nrg%s2p = n%s2p * T%s2p
  nrg%s3p = n%s3p * T%s3p
  nrg%op = n%op * T%op
  nrg%o2p = n%o2p * T%o2p

  ni=n
  np=n

  Ti=T
  Tp=T

  hi=h
 
  nTi=nrg
  nTp=nrg

  vi=v

  lati=lat

  call get_scale_heights(h, T, n)

!open files for outputting data points. These files will be plotted with gplot (or other plotting software)
  open(unit=101, file='ftsp.dat', status='unknown')
  open(unit=102, file='fts2p.dat', status='unknown')
  open(unit=103, file='fts3p.dat', status='unknown')
  open(unit=104, file='ftop.dat', status='unknown')
  open(unit=105, file='fto2p.dat', status='unknown')

  output_it = 0 !This variable determine when data is output. 

  do i=1, nit
    tm = tm0 + (i-1) * dt / 86400

    var =exp(-((tm-neutral_t0)/neutral_width)**2)

    net_source = net_source0*(1.0 + neutral_amp*var)
    ind%o_to_s = o_to_s
!    ind%o_to_s = (otos + o2s_spike * neutral_amp * var) & !o2s_spike
!               / (1.0 + neutral_amp * var)

    n%fh  = fehot_const * (1.0 + hote_amp * var)
    ni%fh = n%fh
    np%fh = n%fh

    n%fc  = 1.0 - n%fh
    ni%fc = n%fc
    np%fc = n%fc

    n%elecHot = n%elec * n%fh/n%fc
    nrg%elecHot = n%elecHot * T%elecHot

    do j=1, LAT_SIZE
      lat%elecHot(j) = n%elecHot
      lati%elecHot(j) = n%elecHot
    end do
    if ( DEBUG_GEN ) then !this variable set in debug.f90
      print *,  "||||||||||||||||||||||||||||||||||||||||||||||"
      print *,  "i = ", i-1
      print *,  "||||||||||||||||||||||||||||||||||||||||||||||"
    endif

    if ( DEBUG_GEN) then !this variable set in debug.f90. 
      print *, "~~~~~~~~~~~~~DENSITY~~~~~~~~~~~~~"
      call output(n)
      print *, "~~~~~~~~~~~~~HEIGHT~~~~~~~~~~~~~~"
      call output(h)
      print *, "~~~~~~~~~~~TEMPERATURE~~~~~~~~~~~"
      call output(T)
      print *, "~~~~~~~~~~~~~~~NU~~~~~~~~~~~~~~~~"
      call output(v)
      print *, "~~~~~~~~~~~~~ENERGY~~~~~~~~~~~~~~"
      call output(nrg)
    endif

    call cm3_latavg_model(n, T, nrg, h, v, ni, Ti, nTi, hi &
                         ,vi, np, Tp, nTp, ind, dep, depi, lat, lati, ft, zoff) 

    call update_temp(n, nrg, T)

    call get_scale_heights(h, T, n)

    call energyBudget(n, h, T, dep, ind, ft, lat, v, nrgy)

    if ( OUTPUT_MIX .and. (nint(output_it)+1 .eq. i) ) then !Output at set intervals when OUTPUT_MIX is true (from debug.f90)
      plot = ftint_mix(n, h) !calculate values to be plotted
      write(101,*) plot%sp
      write(102,*) plot%s2p
      write(103,*) plot%s3p
      write(104,*) plot%op
      write(105,*) plot%o2p
      output_it=output_it + (86400/dt) !Determines when data is output. Set for once each run day (86400/dt).
    endif        

  end do

  close(101)
  close(102)
  close(103)
  close(104)
  close(105)

  print *, '$$--------------------------------'
  print *, '$$ IN-CODE ENERGY BUDGET'
  print *, '$$--------------------------------'
  print *, '$$ ionized S............', nrgy%s_ion
  print *, '$$ ionized O............', nrgy%o_ion
  print *, '$$ charge exchange S....', nrgy%s_cx
  print *, '$$ charge exchange O....', nrgy%o_cx
  print *, '$$ equil with ehot......', nrgy%elecHot_eq + nrgy%tot_eq
  print *, '$$ total in.............', nrgy%P_in + nrgy%tot_eq
  print *, '$$ puv..................', nrgy%Puv
  print *, '$$ fast/ena.............', nrgy%pfast - nrgy%tot_eq
  print *, '$$ transport............', nrgy%ptrans + nrgy%ptrans_elecHot
  print *, '$$ total out............', nrgy%P_out - nrgy%tot_eq
  print *, '$$ in/out...............', (nrgy%P_in + nrgy%tot_eq )/(nrgy%P_out - nrgy%tot_eq )
  print *, ""
  print *, '++++++++++++++++++++++++++++++++++++'
  print *, 'Final Variable Values'
  print *, '++++++++++++++++++++++++++++++++++++'
  print *, 'O/S.........................', o_to_s
  print *, 'Fraction of Hot Electrons...', fehot_const
  print *, 'Transport...................', transport
  print *, 'Hot Electron Temp...........', tehot
  print *, 'Lag Constant................', lag_const
  print *, 'Neutral Amplitude...........', neutral_amp
  print *, 'Inital Neutral Temperature..', neutral_t0
  print *, 'Neutral Width...............', neutral_width
  print *, 'Hot Electron Amplitude......', hote_amp
  print *, 'Hot Electron Initial Temp...', hote_t0
  print *, 'Hot Electron Width..........', hote_width

end subroutine model

END PROGRAM Onebox
