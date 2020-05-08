


program sphere

	implicit none

integer::i,j,k,m,imax,jmax,kmax,n,z,w,i1,wmax,wf
real, allocatable ::TSTAR(:,:,:,:)
real, allocatable ::t(:,:,:),told(:,:,:),tzegond(:,:,:),apk(:,:,:),qk(:,:,:),ank(:,:,:),ask(:,:,:),awi(:,:,:),aei(:,:,:),api(:,:,:),qi(:,:,:),apr(:,:,:),qr(:,:,:),aur(:,:,:),abr(:,:,:)
real, allocatable ::se(:,:),sw(:,:),An(:,:),As(:,:),Au(:,:),Ab(:,:),vol(:,:)
real, allocatable ::aco(:),bco(:),cco(:),dco(:),eco(:),fco(:),tpime(:),ar(:),br(:),cr(:),dr(:),srr(:),r(:),rphalf(:),rmhalf(:),teta(:),phi(:),ss(:),sn(:),Aw(:),Ae(:),tav1(:),tavmax(:),ak(:),bk(:),ck(:),dk(:),sk(:)
real::steta,sr,ti,t0,st,epsilonx,epsilont,ro,alpha,resuteta,resuteta1,resur,resur1,resuphi,resuphi1,TJZERO,n1,n2,dn,fact,fpn,fpn1,a1,b1,su,sb,sphi,resuprime,resuprime1,resut,resut1,cprime,czegond,pie,C1,C
imax=50
jmax=25
kmax=30
st=.1
w=1
allocate (TSTAR(IMAX,IMAX,Imax,1000))
allocate (t(imax,Imax,Imax),told(imax,Imax,Imax),tzegond(imax,Imax,Imax),apk(imax,Imax,Imax),qk(imax,Imax,Imax),ank(imax,Imax,Imax),ask(imax,Imax,Imax),awi(imax,Imax,Imax),aei(imax,Imax,Imax),api(imax,Imax,Imax),qi(imax,Imax,Imax),apr(imax,Imax,Imax),qr(imax,Imax,Imax),aur(imax,Imax,Imax),abr(imax,Imax,Imax))
allocate (se(Imax,Imax),sw(Imax,Imax),An(Imax,Imax),As(Imax,Imax),Au(Imax,Imax),Ab(Imax,Imax),vol(Imax,Imax))
allocate (aco(imax),bco(imax),cco(imax),dco(imax),eco(imax),fco(imax),tpime(imax),ar(imax),br(imax),cr(imax),dr(imax),srr(imax),r(imax),rphalf(imax),rmhalf(imax),teta(imax),phi(imax),ss(imax),sn(imax),Aw(imax),Ae(imax),tav1(imax),tavmax(imax),ak(imax),bk(imax),ck(imax),dk(imax),sk(imax))
call properties
call meshgeometry
call initial

do
told(:,:,:)=t(:,:,:)


do i=1,imax-1
 do j=1,jmax
 do k=1,kmax
TSTAR(i,j,k,w)=t(i,j,k)
end do
end do
end do
w=w+1
 do 

 tzegond(:,:,:)=t(:,:,:)
m=m+1
call tetaphisweep
call rphisweep
call tetarsweep
call resualprime

print *, resuprime,resut,w
if(resuprime<.0001) exit
end do
call resualt
if(resut<.001) exit
do k=1,kmax
T(1,JMAX,k)=1.0
T(IMAX-1,JMAX,K)=1.0
end do
end do
C1=0.0
DO I=1,IMAX-1
DO K=1,KMAX
T(1,1,1)=T(1,1,1)+T(I,2,K)
C1=C1+1
END DO
END DO
T(1,1,1)=T(1,1,1)/C1
wmax=w

do i=1,imax-2
 do j=1,jmax
  do k=1,kmax
TSTAR(i,j,k,wmax)=t(i,j,k)
end do
end do
end do

call  ope
call outp

pause
contains
!---------------------------
subroutine properties
ti=50.
t0=100.
alpha=1.0
end subroutine properties
!---------------------------
subroutine meshgeometry
pie=3.141592
sr=1.000/(jmax-2)
steta=2.*pie/(imax-3)
sphi=pie/(kmax-2)
teta(1)=0.0
teta(2)=steta/2.000
do i=3,imax-2
teta(i)=teta(i-1)+steta
end do
teta(imax-1)=2.00*pie
r(1)=0.00
r(2)=sr/2.0
do j=3, jmax-1
r(j)=r(j-1)+sr
end do
r(jmax)=r(jmax-1)+sr/2.
phi(1)=0.00
phi(2)=sphi/2.00
do k=3,kmax-1
phi(k)=phi(k-1)+sphi
end do
phi(kmax)=phi(kmax-1)+sphi/2.0
rmhalf(1)=0.0
rphalf(1)=0.0
rmhalf(2)=0.0
rphalf(2)=sr
do j=3,jmax-1
rphalf(j)=(r(j)+r(j+1))/2.
rmhalf(j)=(r(j)+r(j-1))/2.
end do
rphalf(jmax)=0.0
rmhalf(jmax)=(r(jmax)+r(jmax-1))/2.
do j=1,jmax
sn(j)=r(j)*sphi
ss(j)=r(j)*sphi
end do
do j=1,jmax
do k=1,kmax
se(j,k)=r(j)*sin(phi(k))*steta
sw(j,k)=r(j)*sin(phi(k))*steta
end do
end do
su=sr
sb=sr
do k=1,kmax
Au(1,k)=0.0
Ab(1,k)=0.0
Ab(2,k)=0.0
Au(2,k)=sr**2*sin(phi(k))*steta*sphi
do j=3,jmax
Au(j,k)=rphalf(j)**2*sin(phi(k))*steta*sphi
Ab(j,k)=rmhalf(j)**2*sin(phi(k))*steta*sphi
end do
end do
do j=1,jmax-1
 do k=1,kmax
An(j,k)=r(j)*sin(phi(k))*steta*sr
As(j,k)=r(j)*sin(phi(k))*steta*sr
vol(j,k)=r(j)**2*sin(phi(k))*steta*sphi*sr
 end do
end do

do j=1,jmax-1
Aw(j)=r(j)*sphi*sr
Ae(j)=r(j)*sphi*sr
end do
end subroutine meshgeometry
!---------------------------------
subroutine initial
do i=1,imax
  do j=1,jmax
  do k=1,kmax
 t(i,j,k)=0.0
 end do
end do  
end do
end subroutine  initial    
!------------------------------
subroutine tetarsweep
tav1(:)=0.0
tavmax(:)=0.0
do j=1,jmax
do i=1,imax-2
tav1(j)=tav1(j)+t(i,j,2)
tavmax(j)=tavmax(j)+t(i,j,kmax-1)
end do
tav1(j)=tav1(j)/(imax-2)
tavmax(j)=tavmax(j)/(imax-2)
end do

  do j=2,jmax-1
 do i=1,imax-2
 
ak(:)=0.0
bk(:)=0.0
ck(:)=0.0
dk(:)=0.0
do k=2,kmax-1 
ak(k)=-alpha*As(j,k)/ss(j)
bk(k)=vol(j,k)/st+alpha*(An(j,k)/sn(j)+As(j,k)/ss(j)+Ae(j)/se(j,k)+Aw(j)/sw(j,k)+Au(j,k)/su+Ab(j,k)/sb)
ck(k)=-alpha*An(j,k)/sn(j)
if(i==1) then
dk(k)=(vol(j,k)/st)*told(i,j,k)+alpha*((Ae(j)/se(j,k))*t(imax-2,j,k)+(Aw(j)/sw(j,k))*t(i+1,j,k)+(Au(j,k)/su)*t(i,j+1,k)+(Ab(j,k)/sb)*t(i,j-1,k))
end if
if(i>1) then
dk(k)=(vol(j,k)/st)*told(i,j,k)+alpha*((Ae(j)/se(j,k))*t(i-1,j,k)+(Aw(j)/sw(j,k))*t(i+1,j,k)+(Au(j,k)/su)*t(i,j+1,k)+(Ab(j,k)/sb)*t(i,j-1,k))
end if
end do
ak(1)=0.0
bk(1)=1.0
ck(1)=0.0
dk(1)=tav1(j)
ak(kmax)=0.0
bk(kmax)=1.0
ck(kmax)=0.0
dk(kmax)=tavmax(j)
ask(i,j,:)=ak(:)
ank(i,j,:)=ck(:)
apk(i,j,:)=bk(:)
qk(i,j,:)=dk(:)
call tdmak
t(i,j,:)=sk(:)
end do
end do
end subroutine tetarsweep
!----------------------------------
subroutine tdmak
real ::fk
n=0
do z=2,kmax
fk=ak(z)/bk(z-1)
bk(z)=bk(z)-ck(z-1)*fk
dk(z)=dk(z)-dk(z-1)*fk
end do
sk(kmax)=dk(kmax)/bk(kmax)
do n=kmax-1,1,-1
sk(n)=(dk(n)-ck(n)*sk(n+1))/bk(n)
end do
return
end subroutine tdmak
!-----------------------------------
subroutine rphisweep

do j=2,jmax-1
 do k=2,kmax-1
aco(:)=0.0
bco(:)=0.0
cco(:)=0.0
dco(:)=0.0
 do i=2,imax-1
aco(i)=-alpha*Ae(j)/se(j,k)
bco(i)=vol(j,k)/st+alpha*(An(j,k)/sn(j)+As(j,k)/ss(j)+Ae(j)/se(j,k)+Aw(j)/sw(j,k)+Au(j,k)/su+Ab(j,k)/sb)
cco(i)=-alpha*Aw(j)/sw(j,k)
dco(i)=(vol(j,k)/st)*told(i,j,k)+alpha*((An(j,k)/sn(j))*t(i,j,k+1)+(As(j,k)/ss(j))*t(i,j,k-1)+(Au(j,k)/su)*t(i,j+1,k)+(Ab(j,k)/sb)*t(i,j-1,k))
end do
aco(1)=-1.0
bco(1)=1.0
cco(1)=0.0
dco(1)=0.0
aco(imax)=0.0
bco(imax)=1.0
cco(imax)=-1.0
dco(imax)=0.0
awi(:,j,k)=aco(:)
aei(:,j,k)=cco(:)
api(:,j,k)=bco(:)
qi(:,j,k)=dco(:)
call period1
t(:,j,k)=tpime(:)
end do
end do
end subroutine rphisweep
!-------------------------------
! this is a periodic solver for cases with two lines overlap
! notot is the total number of grid points in the periodic direction including
! the overlap ones
! aco, bco, cco, dco are coefficents for the tridiagonal solver
! eco, and fco are working vectors with the same dimension as the 
! coefficents for the tridiagonal solver
! tpime is the solution vector

      subroutine period1
      
!
      n2=imax-2
      n1=imax-1
      dn=0.
!
! begin with the first row and start the elimination
!                                                   
      eco(1)=aco(1)
!
! now move to the second row
!     
      fact=-aco(2)/bco(1)                      
      bco(2)=bco(2)+fact*cco(1)
      dco(2)=dco(2)+fact*dco(1)
      eco(2)=fact*eco(1)
      do  i=3,imax-3      
         fact=-aco(i)/bco(i-1)
         bco(i)=bco(i)+fact*cco(i-1)
         dco(i)=dco(i)+fact*dco(i-1)
         eco(i)=fact*eco(i-1)
     end do
!
! now row ntot-2
!
      fact=-aco(n2)/bco(n2-1)               
      bco(n2)=bco(n2)+fact*cco(n2-1)
      eco(n2)=cco(n2)+fact*eco(n2-1)
      dco(n2)=dco(n2)+fact*dco(n2-1)
!
! now ntot-1
!           
      fact=-aco(n1)/bco(n2)
      bco(n1)=bco(n1)+fact*eco(n2)
      dco(n1)=dco(n1)+fact*dco(n2)
!
! now the final row where special things happen
!                                              
      fco(2)=cco(imax)
      fpn1=aco(imax)
      fpn=bco(imax)
!
! now go across the final row
!                            
      fact=-fco(2)/bco(2)
      fco(3)=fact*cco(2)
      fpn1=fpn1+fact*eco(2)
      dn=dn+fact*dco(2)
      do  i=4,n2
         fact=-fco(i-1)/bco(i-1)
         fco(i)=fact*cco(i-1)
         fpn1=fpn1+fact*eco(i-1)
         dn=dn+fact*dco(i-1)
      end do
!
! now at n2
!          
      fact=-fco(n2)/bco(n2)
      dn=dn+fact*dco(n2)
      fpn1=fpn1+fact*eco(n2)
!
! now at n1
!          
      a1=fpn-fpn1*cco(n1)/bco(n1)
      b1=dn-fpn1*dco(n1)/bco(n1)
!
! start the back solve
!      
      tpime(imax)=b1/a1
      tpime(2)=tpime(imax)
      tpime(n1)=(dn-fpn*tpime(imax))/fpn1               
      tpime(1)=tpime(n1)
      tpime(n2)=(dco(n2)-eco(n2)*tpime(n1))/bco(n2)
!
! now we are ready for a recursion relationship
!          
      do  i1=3,n2-1
         i=imax-i1
         tpime(i)=(dco(i)-eco(i)*tpime(n1)-cco(i)*tpime(i+1))/bco(i)         
     end do
      return
      end subroutine period1
!----------------------------------
subroutine tetaphisweep
TJZERO=0.0
C=0.0
DO I=1,IMAX-2
 do k=1,kmax
TJZERO=(TJZERO+T(I,2,k))
c=c+1
END DO
end do
TJZERO=TJZERO/(c)
do i=1,imax-2
 do k=2,kmax-1
ar(:)=0.0
br(:)=0.0
cr(:)=0.0
dr(:)=0.0
 do j=2,jmax-1
ar(j)=-alpha*Ab(j,k)/sb
br(j)=vol(j,k)/st+alpha*(An(j,k)/sn(j)+As(j,k)/ss(j)+Ae(j)/se(j,k)+Aw(j)/sw(j,k)+Au(j,k)/su+Ab(j,k)/sb)
cr(j)=-alpha*Au(j,k)/su
if(i==1) then
dr(j)=(vol(j,k)/st)*told(i,j,k)+alpha*((Ae(j)/se(j,k))*t(imax-2,j,k)+(Aw(j)/sw(j,k))*t(i+1,j,k)+(An(j,k)/sn(j))*t(i,j,k+1)+(As(j,k)/ss(j))*t(i,j,k-1))
end if
if(i>1) then
dr(j)=(vol(j,k)/st)*told(i,j,k)+alpha*((Ae(j)/se(j,k))*t(i-1,j,k)+(Aw(j)/sw(j,k))*t(i+1,j,k)+(An(j,k)/sn(j))*t(i,j,k+1)+(As(j,k)/ss(j))*t(i,j,k-1))
end if
end do


ar(1)=0.0
br(1)=1.0
cr(1)=0.0
dr(1)=TJZERO
ar(jmax)=0.0
if(teta(i)>=0.0.and.teta(i)<pie/2.) then
br(jmax)=1.0
cr(jmax)=0.0
dr(jmax)=1.0
end if
if(teta(i)>=pie/2..and.teta(i)<pie) then
br(jmax)=1.0
cr(jmax)=0.0
dr(jmax)=0.0
end if
if(teta(i)>=pie.and.teta(i)<3.*pie/2.) then
ar(jmax)=-1.0
br(jmax)=1.0
cr(jmax)=0.0
dr(jmax)=0.0
end if
if(teta(i)>=3.*pie/2..and.teta(i)<2.*pie) then
br(jmax)=1.0
cr(jmax)=0.0
dr(jmax)=3.0
end if
abr(i,:,k)=ar(:)
aur(i,:,k)=cr(:)
apr(i,:,k)=br(:)
qr(i,:,k)=dr(:)
call tdmar
t(i,:,k)=srr(:)
end do
end do
end subroutine tetaphisweep
!--------------------------------------
subroutine tdmar
real ::fr
n=0
do z=2,jmax
fr=ar(z)/br(z-1)
br(z)=br(z)-cr(z-1)*fr
dr(z)=dr(z)-dr(z-1)*fr
end do
srr(jmax)=dr(jmax)/br(jmax)
do n=jmax-1,1,-1
srr(n)=(dr(n)-cr(n)*srr(n+1))/br(n)
end do
return
end subroutine tdmar
!--------------------------------------
subroutine resualt
resut=0.0
czegond=0.0
do j=2,jmax-1
 do k=1,kmax
 do i=2,imax-1
  
 resut=resut+abs(t(i,j,k)-told(i,j,k))
 czegond=czegond+1
  end do
end do
end do

resut=resut/czegond
end subroutine resualt
!----------------------------------
subroutine resualprime
resuprime=0.0
cprime=0.0
do j=1,jmax-1
do k=1,kmax
do i=1,imax-2
  resuprime=resuprime+abs(t(i,j,k)-tzegond(i,j,k))
  cprime=cprime+1
 end do
 end do
 end do
resuprime=resuprime/cprime
end subroutine resualprime
!-----------------------------------
subroutine ope
open(1, file='t1.txt')
open(2, file='t5.txt')
open(3, file='t10.txt')
open(4, file='t15.txt')
open(5, file='t20.txt')
open(6, file='t.txt')
open(7, file='ts.txt')
open(8, file='t90.txt')
open(9, file='PHI.txt')
end subroutine ope
!----------------------------------
subroutine outp
 
write (1,*) 'VARIABLE= "X","Y","T1"'
write (1,*)  "ZONE ,j=",JMAX ,"I=",IMAX-1,"k=",kmax
write (2,*) 'VARIABLE= "X","Y","T5"'
write (2,*)  "ZONE ,j=",JMAX ,"I=",IMAX-1,"k=",kmax
write (3,*) 'VARIABLE= "X","Y","T10"'
write (3,*)  "ZONE ,j=",JMAX ,"I=",IMAX-1,"k=",kmax
write (4,*) 'VARIABLE= "X","Y","T15"'
write (4,*)  "ZONE ,j=",JMAX ,"I=",IMAX-1,"k=",kmax
write (5,*)  'VARIABLE= "X","Y","T20"'
write (5,*)   "ZONE ,j=",JMAX ,"I=",IMAX-1,"k=",kmax
write (6,*)  'VARIABLE= "X","Y","Ts"'
write (6,*)   "ZONE ,J=",JMAX ,"I=",IMAX-1,"K=",Kmax
write (7,*)   "ZONE ,I=",IMAX ,"J=",JMAX
do k=1,kmax 
do J=1,Jmax

do I=1,Imax-1
   
 write (1,*)    R(J)*sin(phi(k))*COS(TETA(I)),R(J)*sin(phi(k))*SIN(TETA(I)),R(J)*COS(phi(k)),tstar(i,j,k,2)
  end do
end do
end do
do k=1,kmax 
do J=1,Jmax

do I=1,Imax-1
   
 write (2,*)    R(J)*sin(phi(k))*COS(TETA(I)),R(J)*sin(phi(k))*SIN(TETA(I)),R(J)*COS(phi(k)),tstar(i,j,k,3)
  end do
end do
end do
do k=1,kmax 
do J=1,Jmax

do I=1,Imax-1
   
 write (3,*)     R(J)*sin(phi(k))*COS(TETA(I)),R(J)*sin(phi(k))*SIN(TETA(I)),R(J)*COS(phi(k)),tstar(i,j,k,4)
  end do
end do
end do
do k=1,kmax 
do J=1,Jmax

do I=1,Imax-1
   
 write (4,*)     R(J)*sin(phi(k))*COS(TETA(I)),R(J)*sin(phi(k))*SIN(TETA(I)),R(J)*COS(phi(k)),tstar(i,j,k,5)
  end do
end do
end do
do k=1,kmax 
do J=1,Jmax

do I=1,Imax-1
   
 write (5,*)    R(J)*sin(phi(k))*COS(TETA(I)),R(J)*sin(phi(k))*SIN(TETA(I)),R(J)*COS(phi(k)),tstar(i,j,k,20)
  end do
end do
end do

do k=1,kmax 
do J=1,Jmax

do I=1,Imax-1


  
 write (6,*)     R(J)*sin(phi(k))*COS(TETA(I)),R(J)*sin(phi(k))*SIN(TETA(I)),R(J)*COS(phi(k)),t(i,j,k)
  end do
end do
end do
do j=1,jmax
do i=1,imax


 write (7,*)    R(J)*sin(phi(20))*COS(TETA(I)),R(J)*sin(phi(20))*SIN(TETA(I)),tstar(I,J,3,wmax)
  end do
end do
DO K=1,KMAX
write (9,*)   K,SIN(PHI(K))
END DO

do j=1,jmax

write (8,*) R(J),TSTAR(1,J,15,5)
end do

end subroutine outp
!----------------------------------------

print*,"Done"
pause
end program  
    