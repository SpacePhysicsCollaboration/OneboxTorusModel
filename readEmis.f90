MODULE ReadEmis

  USE varTypes

  IMPLICIT NONE

  CONTAINS

  subroutine ReadIndices(temp, dens)
    real              ::temp(EMIS_SIZE), dens(EMIS_SIZE)
 
    call ReadIndex('emisTemp.dat', temp)  !make sure to change character lengh if filenames change    
    call ReadIndex('emisDens.dat', dens)      

  end subroutine ReadIndices
  
  subroutine ReadIndex(loc, array)
    real              ::array(EMIS_SIZE)
    character(len=12) ::loc !location or filename of data
    integer           ::i

    
    open(unit=10, file=loc, status="old")
    do i=1, EMIS_SIZE
      read(10,*) array(i)
    end do
    close(10)

  end subroutine ReadIndex

  subroutine ReadEmisTable(loc, emis)
    real              ::emis(EMIS_SIZE,EMIS_SIZE)
    character(len=20) ::loc  !location of data (filename)
    integer           ::i,j

    open(unit=10, file=loc, status="old")

    do i=1, EMIS_SIZE
      do j=1, EMIS_SIZE
        read(10,*), emis(j, i)
      end do
    end do
   
    do i=1, EMIS_SIZE
      do j=1, EMIS_SIZE
        emis(j, i)=emis(j,i)/(1.60217646e-12)
      end do
    end do
!    print *, emis(1,1), emis(1,101) , emis(101,1), emis(101, 101) 

  end subroutine ReadEmisTable

END MODULE ReadEmis

