! Constitutive law: Unified Hardening model for overconsolidated clay
! Reference: Zhou X, Lu DC, Timon R, Du XL.  A return-free unconstrained stress updating 
!            algorithm  for the modified Cam-clay model 
! Author: Zhou xin (zhouxin615@126.com)
      
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)
c
      include 'ABA_PARAM.INC'
      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),
     2 ddsddt(ntens),drplde(ntens),
     3 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     4 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3),
     5 jstep(4)
      dimension dstres(6),d(3,3)

! This code is only available for 3D stress state, so NTENS=6
      double precision lambda,kappa,nu,e0,GK,devep,norm_g,tol,e1,
     1 pt,qt,ft,p,q,pc,dphi,Mf,cp,ck,r,Beta,Ftol,cd,f,chi0,chi1,dev,
     2 pold,qold,pcold,BKw,GKw,nsdt,alpha,alphaj,F0,F1,rho,zeta

      double precision s(6,1),dstra(6,1),dstra1(6,1),De(6,6),
     1 Jg(4,4),g(4,1),dx(4,1),II(6,6),IP(6,6),Ivol(6,6),Isym(6,6),
     2 st(6,1),dum(1,1),x(4,1),paras(4,1),sdold(6,1),delta(6,1),
     3 dgamma(6,1),xold(4,1),nw(6,1),xs(4,1),xcd(3,1),IJg(4,4),
     4 xsold(4,1),xcdold(3,1),dgammaold(6,1),nsdtold,pcold1

      integer K1,K2,Iteration,Iteration1
      DATA Maxit,Mmax /100, 10/
!   The Material parameters of MCC model --------------------------------------------------------------
!   props(1) - lambda   ! compression index
!   props(2) - kappa    ! swell index
!   props(3) - Mf       ! critical state stress ratio
!   props(4) - nu       ! Poisson's ratio
!   props(5) - e1       ! initial void ratio  
      Beta = 1.D-12/2.0d0   ! smoothing parameter
      Ftol =sqrt(2.0d0*Beta)    ! error tolerance for the NMTR method
      tol = 1.D-12              ! tolerance that is a threshold to determine when to use the L'hospital rule to avoid a zero denominator
      rho = 1.D-4
      zeta = 0.1d0
      
      delta = 0.0d0
!   initial values of stress, the strain increment and internal variables
      do K1=1,3
	  s(K1,1) = -stress(K1)
        s(K1+3,1) = -stress(K1+3)
        dstra(K1,1) = -dstran(K1)
        dstra(K1+3,1) = -dstran(K1+3)/2.0d0
        dstra1(K1,1) = -dstran(K1)
        dstra1(K1+3,1) = -dstran(K1+3)
        delta(K1,1) = 1.0d0
      end do      

      lambda = props(1)
      kappa = props(2)
      Mf = props(3)
      nu = props(4)
      e1 = props(5)
!      e0 = props(5)
      if (time(2).LT.1.D-7) then
         e1 = props(5) 
         call K0e0(e0,s,e1,lambda,kappa,Mf,coords) 
         statev(6) = e0
      else    
         e0 = statev(6)
      end if   
      
      
!      write(*,*)"void ratio,",e0,statev(6)
      cp = (1.0D0+e0)/(lambda-kappa)
      ck = (1.0D0+e0)/kappa
      r = 3.0d0*(1.0d0-2.0d0*nu)/(2.0d0*(1.0d0+nu))
      paras(:,1) = (/Mf,cp,ck,r/)
      
      call TenMat1(II,IP,Ivol,Isym)     ! compute the fourth and second order unit tensors in matrix form

      dev = dstra(1,1)+dstra(2,1)+dstra(3,1)   ! the total volume strain
      dgamma = dstra - delta*dev/3.0d0         ! the total deviatoric strain

      call pqs(s,pold,qold)              ! calculate hydrostatic stress, the generalized shear stress, and the deviatoric stress at the previous step.
	do  K1 = 1,3
	    sdold(K1,1) = s(K1,1)-pold
          sdold(K1+3,1) = s(K1+3,1)
      end do      
      
      if (time(2).LT.1.D-7) then
         pcold =(Mf*Mf*pold*pold+qold*qold)/(Mf*Mf*pold)    
         pcold1 = dexp((e1-e0-kappa*dlog(pold))/(lambda-kappa))

!      write(*,*)"pcold,e1,e0,lam,kap,pold",pcold,e1,e0,lambda,kappa,pold
!      write(*,*)"coords,noel,npt,qold",coords,noel,npt,qold
         pcold = max(pcold,pcold1)
      else    
         pcold = statev(1)                     ! pre-consolidation pressure at the previous step
      end if
      dphi = 0.0d0                             ! the initial value of plastic multiplier
! four variable: hydrostatic pressure,generalized shear stress,
!                pre-consolidation pressure, and the plastic multiplier   
      x(:,1) =(/pold,qold,pcold,dphi/)         
      xold = x
! four variable: the secant volume modulus(BKw), the secant shear modulus (GKw),
!                the elastic volume strain (deve),and the eta in Eq.(A8) in the paper (zhou et al. 2021)
      xs = 0.0d0                               

      BKw = pold*(dexp(ck*dev)-1.0d0)/dev 
      if (abs(dev).LT.tol) then
         BKw = pold*ck
      end if
      GKw = BKw*r
      
      De = 0.0d0
      De = (3.0d0*BKw - 2.0d0*GKw)*Ivol+2.0d0*GKw*Isym  ! elastic stiffness matrix
      st = matmul(De,dstra1) + s               ! calculate the trial stress
      
      call norms(st,cd)
      cd = max(cd,1.0d0)**3.0d0   ! dimensional parameter in the smoothing function
      xcd = 0.0d0
      xcd(1,1) = cd        ! xcd(:,1)=(/cd,chi0,chi1/)  chi0 and chi1  in Eq.(A11) in the paper (zhou et al. 2021)  
      
      pt = (st(1,1)+st(2,1)+st(3,1))/3.0D0
      qt = sqrt(((st(1,1)-st(2,1))**2.0D0+(st(2,1)-st(3,1))**2.0D0
     1 +(st(3,1)-st(1,1))**2.0D0+6.0D0*(st(4,1)**2.0D0+st(5,1)**2.0D0
     2 +st(6,1)**2.0D0))/2.0D0)
      x(:,1) =(/pt,qt,pcold,dphi/)     ! the trial stress point is used as the initial point
!      x = xold
      Iteration = 0.0d0
      call gnwf(g,nw,paras,sdold,dev,dgamma,pold,pcold,xs,xcd,x,nsdt)
      norm_g=sqrt(g(1,1)*g(1,1)+g(2,1)*g(2,1)+g(3,1)*g(3,1)
     1  + g(4,1)*g(4,1))      
!      call Trust(paras,sdold,dev,dgamma,xold,xs,xcd,x,Iteration)       ! call NMTR method subroutine
      do while ((norm_g.GT.Ftol).and.(Iteration .LT. Maxit))
         Iteration =Iteration + 1
         call Jacob(Jg,nw,paras,dgamma,pold,pcold,xs,xcd,x,nsdt)
         call inverse(Jg,IJg,4)
         dx = -matmul(IJg,g)
         
         alpha = 1.0d0
         F0 = 0.5d0*norm_g
         
         Iteration1 = 0
         xsold = xs
         xcdold = xcd
         dgammaold = dgamma
         nsdtold = nsdt
         
         do while (Iteration1.LT.Mmax) 
            call gnwf(g,nw,paras,sdold,dev,dgamma,pold,pcold,xs,xcd,
     +   x+alpha*dx,nsdt)
              norm_g=sqrt(g(1,1)*g(1,1)+g(2,1)*g(2,1)+g(3,1)*g(3,1)
     1  + g(4,1)*g(4,1))            
            F1 = 0.5d0*norm_g
              
            if (F1.LT.((1.0d0-2.0d0*rho*alpha)*F0)) then
                go to 192
            else
                alpha = max(zeta*alpha,F0/(F0+F1))
            endif
            
            Iteration1 = Iteration1 + 1
        end do        
 
          
192     continue    
        
         xs = xsold
         xcd = xcdold
         dgamma = dgammaold
         nsdt =   nsdtold
!        write(*,*)"alpha",alpha
!        alpha = 1.0d0
        x = x + alpha*dx         

         call gnwf(g,nw,paras,sdold,dev,dgamma,pold,pcold,xs,xcd,x,nsdt)
         norm_g=sqrt(g(1,1)*g(1,1)+g(2,1)*g(2,1)+g(3,1)*g(3,1)
     1  + g(4,1)*g(4,1))         

      end do
      
      
! compute residual function vector (g),the direction of deviatoric stress tensor (nw), xs, xcd, and nsdt
! the nsdt is 2-norm of elastic deviatoric stress     
       call gnwf(g,nw,paras,sdold,dev,dgamma,pold,pcold,xs,xcd,x,nsdt)
!      norm_g=sqrt((g(1,1)*g(1,1)+g(2,1)*g(2,1)+g(3,1)*g(3,1))/(cd*cd)+
!     1g(4,1)*g(4,1))            ! 2-norm of residual function vector
       norm_g=sqrt(g(1,1)*g(1,1)+g(2,1)*g(2,1)+g(3,1)*g(3,1)
     1  + g(4,1)*g(4,1))                  ! 2-norm of residual function vector        
! output the four independent variables in the residual function      
      p = x(1,1)
      q = x(2,1)
      pc = x(3,1)
      dphi = x(4,1)
! update the stress
      s = p*delta + sqrt(2.0d0/3.0d0)*q*nw      
! compute the smoothing consistent tangent operator
      call ConsStiff(nw,IP,Ivol,DDSDDE,paras,dgamma,pold,xs,xcd,x)
! compute the smoothing continuum tangent operator       
!      call ContStiff(Ivol,Isym,DDSDDE,paras,s,pold,xs,xcd,x)
      
! The time step is required to decrease when the solution of NMTR is divergent   
      if ((Iteration .GE. 100).or.isnan(norm_g).or.
     2isnan(norm2(DDSDDE)).or.isnan(norm2(s))) then
         PNEWDT = 0.25d0
         DDSDDE = -1.D60*II
!         write(*,*)"Iteration",Iteration
!         write(*,*)"s",s
!         write(*,*)"xold",xold
!         write(*,*)"dstran",dstran
         s = -1.D10
         pc = -1.D10         
      end if
       
      do K1=1,6
         stress(K1) = -s(K1,1)
      end do
!   statev(1) - pc        ! pre-consolidation pressure 
!   statev(2) - Iteration ! the iteration number of NMTR
!   statev(3) - norm_g    ! 2-norm of residual function vector
!   statev(4) - f         ! the value of yield function 
!   statev(5) - statev(5) ! Counter 
      f = q*q/(Mf*Mf) + p*(p-pc)
      statev(1) = pc
      statev(2) = Iteration
      statev(3) = norm_g
      statev(4) = f
      statev(5) = g(4,1)
      statev(7) = statev(7) + 1
      statev(8) = dphi

      end


! The subroutine computing the smoothing continuum tangent operator  
	subroutine ContStiff(Ivol,Isym,DDSDDE,paras,s,pold,xs,xcd,x)
      include 'ABA_PARAM.INC'
!	implicit none
	double precision fp,h,p,q,pc,dphi,Mf,cp,chi0,chi1,pold,BKw,GKw
      double precision Ivol(6,6),Isym(6,6),DDSDDE(6,6),De(6,6),
     1 fs(6,1),DfsD(6,6),fsDfs(1,1),x(4,1),paras(4,1),sd(6,1),s(6,1),
     2 xs(4,1),xcd(3,1)
      integer K1,K2
      Mf = paras(1,1)
      cp = paras(2,1)
      BKw = xs(1,1)
      GKw = xs(2,1)
      
      chi0 = xcd(2,1)
      chi1 = xcd(3,1)
         
      De = 0.0d0
      fs = 0.0d0
      DfsD = 0.0d0      
      
      p = x(1,1)
      q = x(2,1)
      pc = x(3,1)
      dphi = x(4,1)
	do  K1 = 1,3
	    sd(K1,1) = s(K1,1)-p
          sd(K1+3,1) = s(K1+3,1)
      end do 
      
      De = (3.0d0*BKw - 2.0d0*GKw)*Ivol+2.0d0*GKw*Isym
! the derivative of yield function to the stress     
      do  K1=1,3
         fs(K1,1) = 3.0d0*sd(K1,1)/(Mf*Mf)+(2.0d0*p-pc)/3.0d0
         fs(K1+3,1) = 6.0d0*sd(K1+3,1)/(Mf*Mf)
      end do
      
      fp = -p                                 ! the derivative of yield function(f) to pre-consolidation pressure (pc)
      h = (2.0d0*p-pc)*cp*pc                  ! the derivative of pre-consolidation pressure (pc) to the plastic multiplier
      
      fsDfs = matmul(matmul(transpose(fs),De),fs)
      DfsD = matmul(matmul(De,fs),matmul(transpose(fs),De))
      
	DDSDDE=De-chi0*DfsD/(chi0*(fsDfs(1,1)-fp*h)-chi1)  ! the smoothing continuum tangent operator  

      return
      
      end
      
! The subroutine computing the smoothing consistent tangent operator       
	subroutine ConsStiff(nw,IP,Ivol,DDSDDE,paras,dgamma,pold,xs,xcd,
     1 x)
      include 'ABA_PARAM.INC'
!	implicit none
	double precision nwgam,a0,a1,a2,a3,a4,a5,b0,b1,b2,c0,c1,c2,c3,
     1 c4,c5,c6,rBK,p,q,pc,dphi,Mf,cp,ck,r,chi0,chi1,tol,BK,pold,
     2 BKw,GKw,deve,eta
      double precision IP(6,6),Ivol(6,6),DDSDDE(6,6),paras(4,1),
     1 Inw(6,6),gamnw(6,6),nwI(6,6),nwnw(6,6),gamI(6,6),x(4,1),
     2 delta(6,1),dgamma(6,1),nw(6,1),xs(4,1),xcd(3,1)
      integer K1,K2
      tol = 1.D-12              ! tolerance that is a threshold to determine when to use the L'hospital rule to avoid a zero denominator
      p = x(1,1)
      q = x(2,1)
      pc = x(3,1)
      dphi = x(4,1)
      Mf = paras(1,1)
      cp = paras(2,1)
      ck = paras(3,1)
      r = paras(4,1)
      BKw = xs(1,1)
      GKw = xs(2,1)
      deve = xs(3,1)
      eta = xs(4,1)  
      
      chi0 = xcd(2,1)
      chi1 = xcd(3,1)
      
      delta = 0.0d0    
      nwgam = 0.0d0
      do  K1=1,3
          delta(K1,1) = 1.0d0
          nwgam = nwgam + nw(K1,1)*dgamma(K1,1)
      end do
      do  K1=4,6
          nwgam = nwgam + 2.0d0*nw(K1,1)*dgamma(K1,1)
      end do      
      
      BK = ck*pold*dexp(ck*deve)
      rBK = (r*BK-GKw)/deve               ! the derivative of shear modulus (f) to the volume strain (pc)
      if (abs(deve).LT.tol) then          ! the L'hospital rule is used when the denominator is close to 0
          rBK = r*pold*ck*ck
      end if
   
      a0 = 1.0d0+pc*cp*dphi+2.0d0*BK*dphi  ! the coefficients in Eqs.(B6) and (B11) in the paper (zhou et al. 2021) 
      a1 = (1.0d0+pc*cp*dphi)/a0
      a2 = -(2.0d0*p-pc)/a0
      a3 = 2.0d0*pc*cp*dphi/a0
      a4 = pc*cp*(2.0d0*p-pc)/(a0*BK)
      a5 = eta
      b0 = chi0*(((nwgam/sqrt(6.0d0)-q*dphi/(Mf*Mf))*a2*    ! the coefficients in Eq.(B15) in the paper (zhou et al. 2021) 
     1 rBK-q*GKw/(Mf*Mf))*12.0d0*q*a5/(Mf*Mf)+((2.0d0*a2-a4)*
     2 p-a2*pc)*BK)+chi1  
      b1 = -chi0*((nwgam/sqrt(6.0d0)-q*dphi/(Mf*Mf))*12.0d0*q*a1
     1 *a5/(Mf*Mf)*rBK+(2.0d0*a1-a3)*p*BK-a1*pc*BK)/b0
      b2 = -chi0*2.0d0*sqrt(6.0d0)*q*GKw*a5/(Mf*Mf*b0)
    
      
      c0 = 2.0d0*GKw*a5                    ! the coefficients in Eq.(B19) in the paper (zhou et al. 2021)   
      c1 = (a1+a2*b1)*BK
      c2 = a2*b2*BK
      c3 = 2.0d0*a5*(a1+a2*b1)*rBK
      c4 = 2.0d0*a2*a5*b2*rBK
      c5 = -2.0d0*sqrt(6.0d0)*q*a5/(Mf*Mf)*(b1*GKw+(a1+a2*b1)*rBK
     1 *dphi)
      c6 = -2.0d0*sqrt(6.0d0)*q*a5/(Mf*Mf)*(b2*GKw+a2*b2*rBK*dphi)
      call TenMat2(dgamma,nw,delta,Inw,gamI,gamnw,nwI,nwnw)    ! the fourth order tensors in the consistent tangent operator 

      DDSDDE = 0.0d0                         ! the smoothing consistent tangent operator  
      DDSDDE =c0*IP+3.0d0*c1*Ivol+c2*Inw+c3*gamI+c4*gamnw+c5*nwI+c6*nwnw

      
      return
      
      end

         