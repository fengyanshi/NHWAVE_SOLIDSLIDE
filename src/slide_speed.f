
      implicit none
      real l_center2water,gamma_rho,Cd,Cf,Cm,Cn
      real grav,uslide,b_slide,T_slide,w_slide,dt,alpha_slide,time,dis
      integer i,j,k
      real F1,F2,F3,FF

!      open(1,file='coulomb_8deg.txt')
      open(1,file='coulomb_0deg.txt')

      grav=9.8
      b_slide=800.
      T_slide=24.0
      l_center2water=-1000.0
      dt=1.0

      alpha_slide=14.0*3.1415926/180.0
!      alpha_slide=34.0*3.1415926/180.0

      gamma_rho=2650.0/1000.0
      Cd=0.3
      Cf=0.00514
      Cm=3.1415926*T_slide/2.0/b_slide


!      Cn=tan(7.0*3.1415926/180.0)
      Cn=tan(0.0*3.1415926/180.0)

      uslide=0.0

       time=0.0
       dis=0.0
      do k=1,100
        FF=b_slide*0.5-l_center2water
        if (FF<0.0)FF=0.0
        if (FF>b_slide)FF=b_slide
        F1=b_slide*gamma_rho+Cm*FF
        F2=-0.5*(Cf*FF/T_slide+Cd)
        F3=(b_slide*gamma_rho-FF)*grav*
     &  (sin(alpha_slide)-Cn*cos(alpha_slide))

        uslide=uslide+dt*1.0/F1*(F2*uslide*uslide+F3)
        l_center2water=l_center2water-uslide*dt
        dis=dis+uslide*dt
        time=time+dt
        write(1,*)time,uslide,dis
        print*,FF,F1,F2,F3
      enddo
      close(1)


      end
