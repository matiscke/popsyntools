      program make_list
c
      IMPLICIT NONE
c
      integer n_coralie
      integer n_disk
      parameter (n_coralie=15)
      parameter (n_disk=12)
      integer nsimul
      parameter (nsimul=100000)
      integer CDnumber
      parameter (CDnumber=0)
      integer method
      parameter (method=2)

      integer idum
      real*8 ran2,gasdev
      external ran2,gasdev
      character*8 date
      character*10 time
      real*8 random
      integer idisk,iok,i,j
      real*8 frac,mstar
      character*170 line

      real*8 feh_bin(1:n_coralie)
      real*8 fraction_feh(1:n_coralie)
      real*8 mdisk(1:n_disk)
      real*8 gamma(1:n_disk)
      real*8 rcore(1:n_disk)
      real*8 h(1:n_disk)
      real*8 psi(1:n_disk)
      real*8 r_in(1:n_disk)

      real*8 lifetime,feh,fpg
      real*8 disk_mass,expo,core_radius,inner_radius
      real*8 sigma,mdot_wind
      integer seed

      ! Units
      real*8 pi
      parameter (pi = 3.141592657)
      real*8 r0
      parameter (r0 = 5.2) ! AU
      real*8 Msun
      parameter (Msun = 2.e33) ! g
      real*8 AU
      parameter (AU = 1.5e13) ! cm
      real*8 feh_sun
      parameter (feh_sun = 0.04)

      WRITE(*,*)'Enter Mstar in Msun units'
      READ(*,*)Mstar

      ! Seed of the random number generator
      call date_and_time(date,time)
      write(*,*) date, "  ", time
      read(time(1:6),"(i6)") idum
      idum = -idum
      write(*,*) 'Seed of random number generator = ',idum

      ! Read data files
      open(unit=17,file='data_Santos_coralie')
      read(17,*)
      read(17,*)
      do i=1,n_coralie
         read(17,*) feh_bin(i),fraction_feh(i)
      end do
      close(17)

      if (method == 1) then
         open(unit=18,file='data_Andrews')
         read(18,*)
         read(18,*)
         read(18,*)
         do i=1,n_disk
            read(18,*) mdisk(i),gamma(i),Rcore(i),h(i),psi(i),r_in(i)
         end do
         close(18)
      end if

      open(unit=2,status='new',file='simulation_list.dat')

      iok = 0
      do while(iok.lt.Nsimul)
         if (method == 1) then
           ! Metallicity; [Fe/H]
            feh = ran2(idum)*1.5-0.8
            random = ran2(idum)
            frac = -1.
            do j=1,n_coralie-1
               if((feh.ge.feh_bin(j)-0.05).and.
     1             (feh.lt.feh_bin(j+1)))then
                  frac = fraction_feh(j)
               end if
            end do
            if (random.gt.frac) then
                  ! Not good
               cycle
            end if
         else
            feh=-0.02 + 0.22 * gasdev(idum)   
         end if
         fpg = 10. ** feh * feh_sun

         ! Disc lifetime
         lifetime = -2.5 * log(ran2(idum))
         ! maximum lifetime = 15 Myr
         if (lifetime .gt. 15) then
            cycle
         end if

         ! Other parameters
         if (method == 1) then
            idisk = int(ran2(idum) * n_disk) + 1 ! Line from "data_Andrews" to use
            disk_mass = mdisk(idisk) !*(1.+(ran2(idum)-0.5)*0.0)
            expo = gamma(idisk) !*(1.+(ran2(idum)-0.5)*0.0)
            core_radius = Rcore(idisk) !*(1.+(ran2(idum)-0.5)*0.0)
            inner_radius = r_in(idisk) !*(1.+(ran2(idum)-0.5)*0.0)
         else
            ! According to distribution from Andrews et al. (2010)
            ! i.e. the same paper, but as statistical distribution
            disk_mass = 10. ** (-1.66 + 0.56 * gasdev(idum))
            disk_mass = disk_mass*Mstar**1.2
            expo = 0.9
            !Md=0.002*(Rc/10AU)**1.6
            core_radius = 10. * (disk_mass / 0.002) ** (1./1.6)
            !core_radius = 30.
            inner_radius = 0.03
            ! Maximum disc mass = 0.16 M_sun
            if (disk_mass > 0.16*Mstar**1.2.OR.disk_mass<0.004*Mstar**1.2) then
               cycle
            end if
         end if

         iok = iok + 1

         ! The data from Andrews et al. assumes fpg=1/100; so compute the
         ! gas disk mass with our value of fpg
         !disk_mass = disk_mass / (100. * fpg)
         ! Note: the integral dr of sigma(r)=2*pi*sigma_0*(r/r0)**(-expo)*exp(-(r/r_core)**(2-expo))
         ! is M_disk=2*pi*sigma_0*(1/r0)**(-expo)*(1/r_core)**(2-expo)/(2-expo)
         ! Here we invert this expression to get sigma_0 from M_disk
c         sigma = disk_mass * (2.+expo) / (2.*pi) * r0**expo *
c     1      core_radius**(-2. - expo)
         !sigma = disk_mass*Msun / ((2.*pi*(1./(r0*AU))**(-expo)*(1./(core_radius*AU))**(2-expo))/(2-expo))
         sigma=(r0*AU)**(-expo)*(core_radius*AU)**(-2+expo)*(2-expo)*disk_mass*Msun/(2*PI)
         ! Convert back to CGS
C         sigma = sigma * Msun / AU**2

         ! External photo-evaporation rate
         mdot_wind = lifetime/1.78/(sigma/100.)**0.4
         mdot_wind = mdot_wind**(-0.4) * 1.e-7

         ! New method to get the seed of the random number generator in planete.f
         seed = INT(abs(log10(fpg) + 1) * 1000.) +
     1      1000 * INT(abs(log10(mdot_wind) + 6.6) * 1000.)

         write(20,*) disk_mass,expo,core_radius,inner_radius,
     1      lifetime,feh,sigma,mdot_wind,fpg,disk_mass*fpg,seed

         ! Write line in simulation_list.dat format
         line(1:3) = 'CD_'
         write(line(4:17),"(i14.14)") CDnumber
         line(18:20) = 'FP_'
         write(line(21:34),"(e14.8)") 1. / fpg
         line(35:37) = 'SI_'
         write(line(38:51),"(e14.8)") sigma
         line(52:54) = 'AI_'
         write(line(55:68),"(e14.8)") inner_radius
         line(69:71) = 'AO_'
         write(line(72:85),"(e14.8)") core_radius
         line(86:88) = 'EX_'
         write(line(89:102),"(e14.8)") expo
         line(103:105) = 'MW_'
         write(line(106:119),"(e14.8)") mdot_wind
         line(120:122) = 'SIM'
         write(line(123:136),"(i14.14)") iok
         line(137:170) = 'AS_0.00000000E+00ST_0.00000000E+00'

         write(2,"(A170)") line
      end do

      write(2,'(A3)') 'END'
      close(2)

      end 
