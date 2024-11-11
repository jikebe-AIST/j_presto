
      module PHYCNS

!****************************************************
!
!     COMMON AREA FOR PHISICAL CONSTANT
!
!****************************************************

      ! Circular constant
      real(8),parameter:: pi = 3.141592653589793d0
      ! Radian-degree conversion unit ( pi / 180.0 )
      real(8),parameter:: rad = pi / 180.d0
      ! Speed of light in vacumm (m/s)
      real(8),parameter:: litspd = 2.99792458d+8
      ! Avogadro number (/MOL)
      real(8),parameter:: avogad = 6.0221367d+23
      ! Electron charge (coulomb)
      real(8),parameter:: eleuni = 1.60217733d-19
      ! Gas-constant (J/MOL*K)
      real(8),parameter:: gascst = 8.31451d0
      ! Boltzman constant (J/K)
      real(8),parameter:: boltz = 1.380658d-23
      ! Joule - calorie conversion unit (J/CAL)
      real(8),parameter:: joucal = 4.184d0
      ! Dielectric constant in vacumm (c**2/N*m) or (farad/m)
      real(8),parameter:: epszer = 8.854187817d-12
      ! Electro-energy constant ((KCAL/MOL)*(electron-unit**2/angstroms))
      ! eekcal = (eleuni**2 * avogad) / (4*pi*epszer*1.d-10*joucal*1.0d-)
      real(8),parameter:: eekcal = 332.06378d0    ! presto
!      real(8),parameter:: eekcal = 18.2223d0**2   ! amber
!      real(8),parameter:: eekcal = 18.222615d0**2 ! gromacs & the others
      ! Small value for machine epsilon
      real(8),parameter:: eps = 1.0d-10

      real(8),parameter:: boltzk = boltz / joucal * avogad * 0.001d0

      end module PHYCNS
