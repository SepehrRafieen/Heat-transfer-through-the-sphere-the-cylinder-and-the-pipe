

    
program pipe_flow

implicit none

character(6)::zz
integer::i,j,k,n,s,smax                         
integer,parameter::imax=82
integer,parameter::jmax=82
integer,parameter::kmax=1000
integer,parameter::single=8             
integer,parameter::double=16              
real(kind=single),dimension(imax,jmax)::T_new,T_old,deltaT
real(kind=single),dimension(imax)::a,b,c,d,T              
real(kind=single),dimension(jmax)::r,An,As,de,dw,v,teta,u,dt      
real(kind=single),dimension(kmax)::Tm,Tw_m,Tn,Nu,z,delta
real(kind=single),dimension(imax,jmax,kmax)::T1
 real(kind=single)::dr                 
 real(kind=single)::Ro=1.0              
 real(kind=single)::deltateta          
 real(kind=single)::dn,ds               
 real(kind=single)::Aw,Ae              
 real(kind=single)::kk=1                
 real(kind=single)::h=5                 
 real(kind=single)::Cp=1                
 real(kind=single)::alfa                
 real(kind=single)::density=1                 
 real(kind=single)::dz_exp
 real(kind=single)::dz
 real(kind=single)::dz_max
 real(kind=single)::pie=4.0*atan(1.0)
 real(kind=single)::su,su1,su2,su3
 real(kind=single)::Um            
 real(kind=single)::Ap 
 integer::counter,kkk
 alfa=kk/density*cp
Ap=pie*(Ro**2)
Um=1000.0
!!!!!!! Meshing in area !!!!!!
dr=Ro/(jmax-1.5)
 r(1)=0
 r(2)=dr/2
 do j=3,jmax
 r(j)=r(j-1)+dr
 end do
 deltateta=(2*pie)/((imax-1)-1)
do i=1,imax-2
teta(i)=(i-1)*deltateta
end do
dn=dr;ds=dr
do j=1,jmax
dw(j)=r(j)*deltateta
de(j)=r(j)*deltateta
end do
Aw=dr;Ae=dr
do j=1,jmax
An(j)=(r(j)+dr/2)*deltateta
As(j)=(r(j)-dr/2)*deltateta
end do
As(2)=0
An(2)=dr*deltateta
do j=1,jmax
v(j)=r(j)*dr*deltateta
end do
do j=1,jmax
u(j)=2*Um*(1-(r(j)/Ro)**2)
end do
u(jmax)=0
dz_exp=(u(jmax-1)*(dr)**2)/(2*alfa)
dz_max=100*dz_exp
!!!  initial temperature !!!
 T1(:,:,1)=0.0
 T_old(:,:)=0.0
 T_new(:,:)=0.0
 T1(:,:,:)=0.0
 

500 format(4f15.5)   
open(210,file="tecplot_result.plt",status="replace",action="write")
write(210,*)'variables = "X" "Y" "Z" "T"'
write(210,*) 'Zone T="Time:',0,',Timeset:',0,'"'
write(210,*) 'strandID=1, SolutionTime=',0
write(210,*) 'i=',imax-1, ', j=',jmax,', k=',kmax,', ZONETYPE=Ordered'
write(210,*) 'DATAPACKING=Point'

do k=1,kmax
do j=1,jmax
do i=1,imax-1
write(210,500)r(j)*cos(teta(i)),r(j)*sin(teta(i)),z(k),T1(i,j,k)
end do
end do
end do
 
close(210)

 
 !!!! Meshing in z space !!!
dz=10*dz_exp
z(1)=0
Nu(1)=0
!!! main loop (length) !!! 
counter=0.0
length: do k=2,100000     
		z(k)=z(k-1)+dz
!T1(:,:,1)=0
! T_old(:,:)=0
! T_new(:,:)=0
		do j=2,jmax-1
		dt(j)=dz/u(j)
		end do
!!!   loop Sweep !!!
		Sweep:do n=1,1000000

!!!    teta sweep !!!
			do i=2,imax-1

					do j=2,jmax-1
						a(j)=-alfa*(As(j)/ds)
						b(j)=(v(j)/dt(j))+alfa*((An(j)/dn)+(As(j)/ds)+(Ae/de(j))+(Aw/dw(j)))
						c(j)=-alfa*(An(j)/dn)

						d(j)=(v(j)/dt(j))*T1(i,j,k-1)+((alfa*Ae)/de(j))*T_old(i-1,j)+((alfa*Aw)/dw(j))*T_old(i+1,j)

						end do 


					
					!!!    Boundry condition (arbitrary) !!!

						a(1)=1;b(1)=-1;c(1)=0;d(1)=0.0
                        if(i>=1 .and. i<int(imax/4))then
						a(jmax)=1.0;b(jmax)=-1.0;c(jmax)=0;d(jmax)=0
						else if(i>=int(imax/4) .and. i<int(imax/2))then
						a(jmax)=0;b(jmax)=1.0;c(jmax)=0;d(jmax)=2.0
						else if(i>=int(imax/2) .and. i<int(3*imax/4))then
						a(jmax)=0;b(jmax)=1.0;c(jmax)=0;d(jmax)=0.0
						else
						a(jmax)=0.;b(jmax)=1.0;c(jmax)=0;d(jmax)=1.0
						end if
						
					call TDMA_Solver(a,b,c,d,T,jmax)
						do j=1,jmax
						T_new(i,j)=T(j)
						end do

				end do 
				!!!   r sweep  !!!

				do j=2,jmax-1

					do i=2,imax-1
						a(i)=-alfa*(Ae/de(j))
						b(i)=(v(j)/dt(j))+alfa*((An(j)/dn)+(As(j)/ds)+(Ae/de(j))+(Aw/dw(j)))
						c(i)=-alfa*(Aw/dw(j))
						d(i)=(v(j)/dt(j))*T1(i,j,k-1)+((alfa*An(j))/dn)*T_new(i,j+1)+((alfa*As(j))/ds)*T_new(i,j-1)
						end do 
						!!!!   Boundry condition (arbitrary)  !!!!!!! 
						a(1)=-1;b(1)=1;c(1)=0;d(1)=0
						a(imax)=0;b(imax)=1;c(imax)=-1;d(imax)=0
						 call periodic_2line_overlap(imax,a,b,c,d,T)
						do i=1,imax
						T_new(i,j)=T(i)
						end do
						T_new(imax,jmax)=T_new(2,jmax)          
						T_new(1,jmax)=T_new(imax-1,jmax)
						T_new(:,1)=(sum(T_new(:,2)))/(imax)
				end do 
				
				do i=1,imax
				do j=2,jmax
				deltaT(i,j)=abs(T_new(i,j)-T_old(i,j))
				end do
				end do

				T_old(:,:)=T_new(:,:)
				!!!!	error  Sweep !!!!
				if (maxval(deltaT)<0.00001)exit
		end do Sweep
		!!!!!	end  Sweep !!!!!

		do i=1,imax
		do j=1,jmax
		T1(i,j,k)=T_new(i,j)
		end do
		end do
		!!!!!!!	calculating Nusselt !!!!!!
		su=0
		do i=1,imax-2
		do j=2,jmax-1
		su=su+(u(j)*T1(i,j,k)*r(j)*dr*deltateta)
		end do
		end do
		Tm(k)=su/(Um*Ap)


		su1=0
		do i=1,imax-2
		su1=su1+(T1(i,jmax,k)-T1(i,jmax-1,k))/(dr)
		end do
		Tn(k)=su1/(imax-2)

		su2=0
		do i=1,imax-2
		su2=su2+T1(i,jmax,k)
		end do
		Tw_m(k)=su2/(imax-2)

		Nu(k)=(2*Ro*Tn(k))/(Tw_m(k)-Tm(k))
        counter=counter+1
		!!!!!!	error  Nusselt !!!!!
		if(abs(Nu(k)-Nu(k-1))<0.001)then
		smax=k
		exit
		end if

		dz=dz*1.05
		if(dz>=dz_max) then
            dz=dz_max
        end if    

		print*,"k=",k,"Nu=",Nu(k),"Z=",z(k),"counter=",counter

open(210,file="tecplot_result.plt",status="old",position="append",action="write")
write(210,*)'variables = "X" "Y" "Z" "T"'
write(210,*) 'Zone T="Time:',counter,',Timeset:',counter,'"'
write(210,*) 'strandID=1, SolutionTime=',counter
write(210,*) 'i=',imax-1, ', j=',jmax,', kkk=',kmax,', ZONETYPE=Ordered'
write(210,*) 'DATAPACKING=Point'

do kkk=1,kmax
do j=1,jmax
do i=1,imax-1
write(210,500)r(j)*cos(teta(i)),r(j)*sin(teta(i)),z(kkk),T1(i,j,kkk)
end do
end do
end do

close(210)

end do length        
!!!!!!!!	end  main (do) loop !!!!!!
!!!!!!!!!	contour !!!!!!!!
open(11,file='T.plt',action='write')
write(11,*)'variables="X,"Y,"T"'

write(11,*)'zone,i=',imax,'j=',jmax
do j=1,jmax
do i=1,imax
write(11,*)r(j)*cos(teta(i)),r(j)*sin(teta(i)),T_new(i,j)
end do
end do


open(13,file='Nusselt.txt',action='write')
 do k=1,smax
write(13,*)z(k),Nu(k),Tm(k),Tw_m(k),Tn(k),dz_max,dz_exp,k
end do

open(23, file='1.plt')
write (23,*) 'VARIABLE= "X","Y","T1"'
write (23,*)  "ZONE ,i=",IMAX-2 ,"I=",JMAX,"k=",smax
do k=2,smax 
 do I=1,Imax-2
   do J=1,Jmax
 write (23,*)   R(J)*COS(TETA(I)),R(J)*SIN(TETA(I)),T1(i,j,int(smax/4))
  end do
 end do
end do
close (23)

open(33, file='2.plt')
write (33,*) 'VARIABLE= "X","Y","T1"'
write (33,*)  "ZONE ,i=",IMAX-2 ,"I=",JMAX,"k=",smax
do k=2,smax 
 do I=1,Imax-2
   do J=1,Jmax


 write (33,*)   R(J)*COS(TETA(I)),R(J)*SIN(TETA(I)),T1(i,j,int(smax/3))
  end do
 end do
end do
close (33)


open(43,file="T_distribution.plt")
write(43,*)'variables="X","Y","Z","T"'
write(43,*)'zone,i=',imax-1,'j=',jmax,'k=',smax

do k=1,int(smax)
do j=1,jmax
do i=1,imax-1
write(43,*)r(j)*cos(teta(i)),r(j)*sin(teta(i)),z(k),T1(i,j,k)
end do
end do
end do
close(43)

open(53, file='3.plt')
write (53,*) 'VARIABLE= "z","R","T1"'
write (53,*)  "ZONE ,k=",smax,"j=",JMAX
do k=2,smax 
   do J=1,Jmax
       write (53,*)   z(k),R(J),T1(1,j,k)
   end do
end do
close (53)

open(63,file="for_checking.plt")
write(63,*)'variables="(x/D)/Pe_D","T_mean","T_wall"'
do k=2,smax
write(63,*)z(k)/Um,Tm(k),Tw_m(k)
end do

open(73,file="answer_Nusselt.plt")
write(73,*)'variables="Z","Nusselt"'
do i=2,smax
write (73,*)z(i),Nu(i)
end do

print*,"Done"
pause
end program 

!-----------------------------------------
                          !TDMA solver
!-----------------------------------------
subroutine TDMA_solver(ai,bi,ci,di,Ti,nn)
integer,intent(in)::nn
integer::j         
integer,parameter::single=8             
integer,parameter::double=16              
real(kind=single),dimension(nn),intent(in)::ai,bi,ci,di
real(kind=single),dimension(nn),intent(out)::Ti
real(kind=single),dimension(nn)::d_a,b_a
b_a(1)=bi(1)
d_a(1)=di(1)
do j=2,nn
b_a(j)=bi(j)-(ci(j-1)*ai(j)/b_a(j-1))
d_a(j)=di(j)-(d_a(j-1)*(ai(j)/b_a(j-1)))
end do
Ti(nn)=d_a(nn)/b_a(nn)
do j=nn-1,1,-1
Ti(j)=(d_a(j)-(ci(j)*Ti(j+1)))/(b_a(j))
end do
    end subroutine TDMA_solver
!---------------------------------------
   !!!  TDMA solver PERIODIC OVERLAP  !!!!
! this is a periodic solver for cases with two lines overlap
! notot is the total number of grid points in the periodic direction including
! the overlap ones
! aco, bco, cco, dco are coefficents for the tridiagonal solver
! eco, and fco are working vectors with the same dimension as the
! coefficents for the tridiagonal solver
! tpime is the solution vector
subroutine periodic_2line_overlap(ntot,aco,bco,cco,dco,tpime)
implicit none
integer,parameter::single=8            
integer,parameter::double=16              
integer,intent(in)::ntot
real(kind=single),dimension(ntot)::aco,bco,cco,dco
real(kind=single),dimension(ntot)::tpime
real(kind=single),dimension(ntot)::eco,fco
integer::n2,n1,i1,i
real(kind=single)::dn,fpn,fpn1,a1,b1,fact

n2=ntot-2
n1=ntot-1
dn=0.

! begin with the first row and start the elimination

eco(1)=aco(1)

! now move to the second row

fact=-aco(2)/bco(1)
bco(2)=bco(2)+fact*cco(1)
dco(2)=dco(2)+fact*dco(1)
eco(2)=fact*eco(1)
do  i=3,ntot-3
fact=-aco(i)/bco(i-1)
bco(i)=bco(i)+fact*cco(i-1)
dco(i)=dco(i)+fact*dco(i-1)
eco(i)=fact*eco(i-1)
end do

! now row ntot-2

fact=-aco(n2)/bco(n2-1)
bco(n2)=bco(n2)+fact*cco(n2-1)
eco(n2)=cco(n2)+fact*eco(n2-1)
dco(n2)=dco(n2)+fact*dco(n2-1)

! now ntot-1

fact=-aco(n1)/bco(n2)
bco(n1)=bco(n1)+fact*eco(n2)
dco(n1)=dco(n1)+fact*dco(n2)

 !now the final row where special things happen

fco(2)=cco(ntot)
fpn1=aco(ntot)
fpn=bco(ntot)

! now go across the final row

fact=-fco(2)/bco(2)
fco(3)=fact*cco(2)
fpn1=fpn1+fact*eco(2)
dn=dn+fact*dco(2)
do i=4,n2
fact=-fco(i-1)/bco(i-1)
fco(i)=fact*cco(i-1)
fpn1=fpn1+fact*eco(i-1)
dn=dn+fact*dco(i-1)
end do

! now at n2

fact=-fco(n2)/bco(n2)
dn=dn+fact*dco(n2)
fpn1=fpn1+fact*eco(n2)

!now at n1

a1=fpn-fpn1*cco(n1)/bco(n1)
b1=dn-fpn1*dco(n1)/bco(n1)

 !start the back solve

tpime(ntot)=b1/a1
tpime(2)=tpime(ntot)
tpime(n1)=(dn-fpn*tpime(ntot))/fpn1
tpime(1)=tpime(n1)
tpime(n2)=(dco(n2)-eco(n2)*tpime(n1))/bco(n2)

 !now we are ready for a recursion relationship
do  i1=3,n2-1
i=ntot-i1
tpime(i)=(dco(i)-eco(i)*tpime(n1)-cco(i)*tpime(i+1))/bco(i)
end do
return
    end subroutine periodic_2line_overlap

!!!! TDMA solver PERIODIC OVERLAP !!!!
