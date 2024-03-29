MODULE DEBUG

  USE varTypes

  IMPLICIT NONE

  INTERFACE output
    module procedure output_density
    module procedure output_energy
    module procedure output_temp 
    module procedure output_height
    module procedure output_nu
    module procedure output_ind 
    module procedure output_dep
    module procedure output_ft 
    module procedure output_all
  END INTERFACE output

  logical             ::HUSH = .true. !controls error output (specifically the interpolate function)
                             ! true means silent, false means printing of error functions
  logical             ::DEBUG_GEN = .true. !controls general debugging output
                                  ! true will produce debugging info at every iteration, false will run silently
  logical             ::OUTPUT_MIX = .true. !Variable to control outputs for plotting
                                   ! true means plot, false means silence outputs

  CONTAINS

!Here are some helpful output functions and a function to initialize a "reaction" data type.

  SUBROUTINE output_density(a)
    type(density)  :: a
    print *, "{", a%sp, a%s2p, a%s3p, a%s4p, a%op, a%o2p     &
                , a%elec, a%elecHot, a%s, a%o, a%fc, a%fh, a%protons, "}" 
  end SUBROUTINE output_density

  SUBROUTINE output_energy(a)
    type(nT)     :: a
    print *, "{", a%sp, a%s2p, a%s3p, a%s4p, a%op, a%o2p, a%elec, a%elecHot, "}"

  end SUBROUTINE output_energy

  SUBROUTINE output_temp(a)
    type(temp)     :: a
    print *,"{", a%sp, a%s2p, a%s3p, a%s4p, a%op, a%o2p, a%elec, a%elecHot &
               , a%pu_s, a%pu_o, "}"

  end SUBROUTINE output_temp

  SUBROUTINE output_height(a)
    type(height)   :: a
    print *, "{", a%s, a%sp, a%s2p, a%s3p, a%s4p, a%o, a%op, a%o2p, a%elec, "}"

  end SUBROUTINE output_height

  SUBROUTINE output_nu(a)
    type(nu)      :: a
    print *, "{", a%sp_s2p, a%sp_s3p, a%sp_s4p, a%sp_op, a%sp_o2p         &
                , a%s2p_s3p, a%s2p_s4p, a%s2p_op, a%s2p_o2p               &
                , a%s3p_s4p, a%s3p_op, a%s3p_o2p, a%s4p_op, a%s4p_o2p     &
                , a%op_o2p, a%sp_elec, a%s2p_elec, a%s3p_elec, a%s4p_elec &
                , a%op_elec, a%o2p_elec, a%sp_elecHot, a%s2p_elecHot      &
                , a%s3p_elecHot, a%s4p_elecHot, a%op_elecHot              &
                , a%o2p_elecHot, a%elec_elecHot, "}"    

  end SUBROUTINE output_nu

  SUBROUTINE output_ind(a)
    type(r_ind)    :: a
    print *, "{", a%ish, a%isph, a%is2ph, a%is3ph, a%ioh, a%io2ph &
                , a%s_production, a%o_production, a%o2s_spike     &
                , a%o_to_s, "}" 

  end SUBROUTINE output_ind

  SUBROUTINE output_dep(a)
    type(r_dep)    :: a
    print *, "{", a%is, a%isp, a%is2p, a%is3p, a%io    &
                , a%iop, a%io2p, a%io2p, a%rsp, a%rs2p &
                , a%rs3p, a%rop, a%ro2p, a%transport, "}"  

  end SUBROUTINE output_dep

  SUBROUTINE output_ft(a)
    type(ft_int)   :: a
    integer        :: i
    print *, "{"
    do i=1, 17
      write(*,"(f15.2, ',' )",advance="no")  a%cx(i) 
    end do
    print *, a%is, ",", a%isp, ",", a%is2p, ",", a%is3p, ","  &
           , a%io, ",", a%iop, ",", a%io2p, ",", a%ish, ","   &
           , a%isph, ",", a%is2ph, ",", a%is3ph,  ","         &
           , a%ioh, ",", a%ioph, ",", a%io2ph,  ","           &
           , a%rsp, ",", a%rs2p, ",", a%rs3p,  ","            &
           , a%rop, ",", a%ro2p, "}"
  end SUBROUTINE output_ft

  SUBROUTINE output_all(a)
    type(reac)     :: a
    call output(a%dens)
    call output(a%nrg)
    call output(a%tmp)
    call output(a%h)
    call output(a%v)
    call output(a%ind)
    call output(a%dep)
    call output(a%ft)
  end SUBROUTINE output_all

  SUBROUTINE init(re, z)
    type(reac)     ::re
    real           ::z, a, stuff, m
    integer        ::i,j
    z=42.0
!!==========DENSITY INITIALIZAION=====================  
    a=10.0
    re%dens%sp     =a
    re%dens%s2p    =a
    re%dens%s3p    =a
    re%dens%s4p    =a
    re%dens%op     =a
    re%dens%o2p    =a
    re%dens%elec   =a
    re%dens%elecHot=a
    re%dens%s      =a
    re%dens%o      =a
    re%dens%fc     =a
    re%dens%fh     =a
    re%dens%protons=a
!!==========HEIGHT INITIALIZAION=====================  
!    a=1.00
    re%h%s   =a
    re%h%sp  =a
    re%h%s2p =a
    re%h%s3p =a
    re%h%s4p =a
    re%h%o   =a
    re%h%op  =a
    re%h%o2p =a
    re%h%elec=a
!!==========TEMP INITIALIZAION=====================  
!    a=100.00
    re%tmp%sp       =a
    re%tmp%s2p      =a
    re%tmp%s3p      =a
    re%tmp%s4p      =a
    re%tmp%op       =a
    re%tmp%o2p      =a
    re%tmp%elec     =a
    re%tmp%elecHot  =a
    re%tmp%pu_s     =a
    re%tmp%pu_o     =a
!!==========FLUX TUBE INITIALIZAION=====================  
!    a=1.00
    do i=1, 17
      re%ft%cx(i) =a
    end do
    re%ft%is    =a
    re%ft%isp   =a
    re%ft%is2p  =a
    re%ft%is3p  =a
    re%ft%io    =a
    re%ft%iop   =a
    re%ft%io2p  =a
    re%ft%ish   =a
    re%ft%isph  =a
    re%ft%is2ph =a
    re%ft%is3ph =a
    re%ft%ioh   =a
    re%ft%ioph  =a
    re%ft%io2ph =a
    re%ft%rsp   =a
    re%ft%rs2p  =a
    re%ft%rs3p  =a
    re%ft%rop   =a
    re%ft%ro2p  =a
!!==========INDEPENDENT RATES INITIALIZAION=====================  
!    a=1.00
    do i=1, 17
      re%ind%cx(i) =a
    end do
    re%ind%ish           =a 
    re%ind%isph          =a 
    re%ind%is2ph         =a 
    re%ind%is3ph         =a 
    re%ind%ioh           =a 
    re%ind%ioph          =a 
    re%ind%io2ph         =a 
    re%ind%s_production  =a 
    re%ind%o_production  =a 
    re%ind%o2s_spike     =a 
    re%ind%o_to_s        =a 
    do i=1, EMIS_SIZE
      do j=1, EMIS_SIZE
        stuff = j*j*i / (i*i*j+1.0)
        re%ind%emis_sp(i,j) = stuff
        re%ind%emis_s2p(i,j)= stuff
        re%ind%emis_s3p(i,j)= stuff
        re%ind%emis_op(i,j) = stuff
        re%ind%emis_o2p(i,j)= stuff
      end do
    end do 
!!==========DEPENDENT RATES INITIALIZAION=====================  
!    a=1.00
    re%dep%is         =a 
    re%dep%isp        =a 
    re%dep%is2p       =a 
    re%dep%is3p       =a 
    re%dep%io         =a 
    re%dep%iop        =a 
    re%dep%io2p       =a 
    re%dep%rsp        =a 
    re%dep%rs2p       =a 
    re%dep%rs3p       =a 
    re%dep%rop        =a 
    re%dep%ro2p       =a 
    re%dep%transport  =a 
!!==========NU INITIALIZAION=====================  
!    a=1.00
    re%v%sp_s2p        =a
    re%v%sp_s3p        =a
    re%v%sp_s4p        =a
    re%v%sp_op         =a
    re%v%sp_o2p        =a
    re%v%sp_elec       =a
    re%v%sp_elecHot    =a
    re%v%s2p_s3p       =a
    re%v%s2p_s4p       =a
    re%v%s2p_op        =a
    re%v%s2p_o2p       =a
    re%v%s2p_elec      =a
    re%v%s2p_elecHot   =a
    re%v%s3p_s4p       =a
    re%v%s3p_op        =a
    re%v%s3p_o2p       =a
    re%v%s3p_elec      =a
    re%v%s3p_elecHot   =a
    re%v%s4p_op        =a
    re%v%s4p_o2p       =a
    re%v%s4p_elec      =a
    re%v%s4p_elecHot   =a
    re%v%op_o2p        =a
    re%v%op_elec       =a
    re%v%op_elecHot    =a
    re%v%o2p_elec      =a
    re%v%o2p_elecHot   =a
    re%v%elec_elecHot  =a
!!==========LAT INITIALIZAION=====================  
!    a=1.00
    do i= 1, LAT_SIZE
      re%lat%z(i)       =a
      re%lat%sp(i)      =a
      re%lat%s2p(i)     =a
      re%lat%s3p(i)     =a
      re%lat%op(i)      =a
      re%lat%o2p(i)     =a
      re%lat%elec(i)    =a
      re%lat%elecHot(i) =a
    end do
!!==========nT INITIALIZAION=====================  
!    a=1000.00
    re%nrg%sp      =a
    re%nrg%s2p     =a
    re%nrg%s3p     =a
    re%nrg%s4p     =a
    re%nrg%op      =a
    re%nrg%o2p     =a
    re%nrg%elec    =a
    re%nrg%elecHot =a

  end SUBROUTINE init

  subroutine initNu(v)
    type(nu)          ::v
    real              ::a=0.0
     
    v%sp_s2p        =a
    v%sp_s3p        =a
    v%sp_s4p        =a
    v%sp_op         =a
    v%sp_o2p        =a
    v%sp_elec       =a
    v%sp_elecHot    =a
    v%s2p_s3p       =a
    v%s2p_s4p       =a
    v%s2p_op        =a
    v%s2p_o2p       =a
    v%s2p_elec      =a
    v%s2p_elecHot   =a
    v%s3p_s4p       =a
    v%s3p_op        =a
    v%s3p_o2p       =a
    v%s3p_elec      =a
    v%s3p_elecHot   =a
    v%s4p_op        =a
    v%s4p_o2p       =a
    v%s4p_elec      =a
    v%s4p_elecHot   =a
    v%op_o2p        =a
    v%op_elec       =a
    v%op_elecHot    =a
    v%o2p_elec      =a
    v%o2p_elecHot   =a
    v%elec_elecHot  =a

  end subroutine initNu

end MODULE DEBUG
