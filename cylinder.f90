


    program cylinder

    implicit none

    integer::i,imax,j,jmax,p,s,counter,k,kmax       !i,imax,j,jmax,s (is counter of time space)& p is counter
    !integer,parameter::single=4              !Compiler dependept value
    !integer,parameter::double=8              !Compiler dependept value
    double precision,dimension(200,200,200)::T_old_1,T_old,T_1
    !double precision,dimension(200,200,200)::TdtAR
    double precision,dimension(200)::a,b,c,d,T,r,An,As,de,dw,vol,teta,aco,bco,cco,dco,tpime,z,Au,Ad   !a,b,c are three diagonals of TDMA; d=Right hand side of TDMA; r is radius(m)
    !; An&As is area(m^2); de&dw is didtance beetwin the node:vol is volume (m^3)
    !!!! properties !!!!
    double precision::dr,resuprime,cprime,resut,czegond,tempAvgJ1,ncount                  !element of radius(m)
    double precision::Ra                 !Ro is radius(radian)
    double precision::dteta           !element of angle is the area(m^2)
    double precision::dz
    double precision::dn,ds               !dn & ds is the didtance between nod(m)
    double precision::Aw,Ae               !Aw&Ae is area(m^2)
    double precision::kconduct=1                 !Thermal conductivity(W/m.K)
    double precision::Cp=1                !Specific heat capacity at condtant pressure(J/kg.K)
    double precision::alfa
    double precision::density=1           !density(kg/m^3)
    double precision::dt                  !time dtep
    double precision,parameter::pie=3.1415
    !double precision::TJZERO


    print*,"enter number of point kmax and jmax and imax"
    print*,"imax="
    read*,imax
    print*,"jmax="
    read*,jmax
    print*,"kmax="
    read*,kmax
    print*,"enter radius"
    read*,Ra
    print*,"enter time dtep dt"
    read*,dt

    !----------------------------------
    !!!  Geometry Mesh  !!!
    !----------------------------------

    alfa=kconduct/(density*Cp)
    dr=Ra/(jmax-1.5)
    dteta=(2*pie)/(imax-2)
    teta(1)=0.0
    do i=2,imax-2
        teta(i)=teta(i-1)+dteta
    end do
     z(1)=0.0
    z(kmax)=5.0
    dz=5.0/(kmax-1)
    do k=2,kmax-1
        z(k)=z(k-1)+dz
    end do
    teta(imax-1)=2*pie
    dn=dr;ds=dr
    Aw=dr*dz;Ae=dr*dz
    r(1)=0
    r(2)=dr/2.
    do i=3,imax
        r(i)=r(i-1)+dr
    end do
      
    do j=2,jmax-1
        An(j)=(r(j)+r(j+1))/2*dteta*dz
        As(j)=(r(j)+r(j-1))/2*dteta*dz
    end do
   
    
    do i=2,imax
        Au(i)=r(i)*dteta*dr
        Ad(i)=r(i)*dteta*dr
    end do
  

    do i=1,imax
        vol(i)=r(i)*dr*dteta*dz
    end do
    do i=2,imax-1
        dw(i)=r(i)*dteta
        de(i)=r(i)*dteta
    end do
    !--------------------------------------------------------
    !!! initial prediction for temperature  !!!
    !--------------------------------------------------------
    counter=0
    do i=1,imax
        do j=1,jmax
            do k=1,kmax  
                T_1(i,j,k)=0
                T_old(i,j,k)=0
            end do
        end do
    end do
 500 format(4f15.5)   
    open(21,file="tecplot_result.plt",status="replace",action="write")
write(21,*)'variables = "X" "Y" "Z" "T"'
write(21,*) 'Zone T="Time:',counter*dt,',Timeset:',dt,'"'
write(21,*) 'strandID=1, SolutionTime=',0
write(21,*) 'i=',imax-1, ', j=',jmax,', k=',kmax,', ZONETYPE=Ordered'
write(21,*) 'DATAPACKING=Point'

 do k=1,kmax
        do j=1,jmax
            do i=1,imax-1
                write(21,500)r(j)*cos(teta(i)),r(j)*sin(teta(i)),z(k),T_1(i,j,k)
            end do
        end do
 end do
 
close(21)

    !------------------------------------------
    !!!  do Time space  !!!
    !------------------------------------------
    do 
        !do i=1,imax
        !do j=2,jmax
        ! TdtAR(i,j,w)=T_old_1(i,j)
        !end do
        !end do
        counter=counter+1

        !--------------------------------------------------
        !!!  do R_sweep and teta_sweep   !!!
        !--------------------------------------------------

        do 


        ! i-sweep
        do i=2,imax-1
            do j=2,jmax-1
                do k=2,kmax-1
                    a(k)=-alfa*(Ad(k)/dz)
                    b(k)=(vol(j)/dt)+alfa*((Au(k)/dz)+(Ad(k)/dz)+(An(k)/dn)+(As(k)/ds)+(Ae/de(k))+(Aw/dw(k)))
                    c(k)=-alfa*(Au(k)/dz)
                    d(k)=(vol(j)/dt)*T_1(i,j,k)+((alfa*Ae)/de(k))*T_old(i-1,j,k)+((alfa*Aw)/dw(k))*T_old(i+1,j,k)&
                    +((alfa*As(k))/ds)*T_old(i,j-1,k)+((alfa*An(k))/dn)*T_old(i,j+1,k)
                end do

               
                
                a(1)=0;b(1)=1;c(1)=0;d(1)=0
                a(kmax)=0;b(kmax)=1;c(kmax)=0;d(kmax)=1
               
                call TDMA_Solver(a,b,c,d,T,kmax)
                do k=1,kmax
                    T_old_1(i,j,k)=T(k)
                end do
            end do 

            do  k=2,kmax-1
                do j=2,jmax-1
                    a(j)=-alfa*(As(j)/ds)
                    b(j)=(vol(j)/dt)+alfa*((Au(j)/dz)+(Ad(j)/dz)+(An(j)/dn)+(As(j)/ds)+(Ae/de(j))+(Aw/dw(j)))
                    c(j)=-alfa*(An(j)/dn)
                    d(j)=(vol(j)/dt)*T_1(i,j,k)+((alfa*Ae)/de(j))*T_old(i-1,j,k)+((alfa*Aw)/dw(j))*T_old(i+1,j,k)&
                    +((alfa*Ad(j))/dz)*T_old(i,j,k-1)+((alfa*Au(j))/dz)*T_old(i,j,k+1)
                end do

                a(1)=0;b(1)=1;c(1)=0;d(1)=0
                a(jmax)=0;b(jmax)=1;c(jmax)=0;d(jmax)=1

                
                call TDMA_Solver(a,b,c,d,T,jmax)
                do j=1,jmax
                    T_old_1(i,j,k)=T(j)
                end do
            end do
        end do
        !-------------------------------------
        !j-sweep
        do  j=2,jmax-1
            do i=2,imax-1
                do k=2,kmax-1
                    a(k)=-alfa*(Ad(k)/dz)
                    b(k)=(vol(j)/dt)+alfa*((Au(k)/dz)+(Ad(k)/dz)+(An(k)/dn)+(As(k)/ds)+(Ae/de(k))+(Aw/dw(k)))
                    c(k)=-alfa*(Au(k)/dz)
                    d(k)=(vol(j)/dt)*T_1(i,j,k)+((alfa*Ae)/de(k))*T_old(i-1,j,k)+((alfa*Aw)/dw(k))*T_old(i+1,j,k)&
                    +((alfa*As(k))/ds)*T_old(i,j-1,k)+((alfa*An(k))/dn)*T_old(i,j+1,k)
                end do 

                a(1)=0;b(1)=1;c(1)=0;d(1)=0
                a(kmax)=0 ; b(kmax)=1 ; c(kmax)=0 ; d(kmax)=1

                call TDMA_Solver(a,b,c,d,T,kmax)
                do k=1,kmax
                    T_old_1(i,j,k)=T(k)
                end do
            end do

            do k=2,kmax-1
                do i=2,imax-1
                    a(i)=-alfa*(Ae/de(i))
                    b(i)=(vol(i)/dt)+alfa*((Au(i)/dz)+(Ad(i)/dz)+(An(i)/dn)+(As(i)/ds)+(Ae/de(i))+(Aw/dw(i)))
                    c(i)=-alfa*(Aw/dw(i))
                    d(i)=(vol(i)/dt)*T_1(i,j,k)+((alfa*As(i))/ds)*T_old_1(i,j-1,k)+((alfa*An(i))/dn)*T_old_1(i,j+1,k)+&
                    ((alfa*Ad(i))/dz)*T_old_1(i,j,k-1)+((alfa*Au(i))/dz)*T_old_1(i,j,k+1)
                end do 

                a(1)=-1.0
                b(1)=1.0
                c(1)=0.0
                d(1)=0.0
                a(imax)=0.0
                b(imax)=1.0
                c(imax)=-1.0
                d(imax)=0.0
                
                call PERIODIC2(imax,a,b,c,d,Tpime)
                do i=1,imax
                    T_old_1(i,j,k)=Tpime(i)
                end do
            end do 
        end do

        !-----------------------------
        !k-sweep

        do k=2,kmax-1
            do i=2,imax-1
                do j=2,jmax-1
                    a(j)=-alfa*(As(j)/ds)
                    b(j)=(vol(j)/dt)+alfa*((Au(j)/dz)+(Ad(j)/dz)+(An(j)/dn)+(As(j)/ds)+(Ae/de(j))+(Aw/dw(j)))
                    c(j)=-alfa*(An(j)/dn)
                    d(j)=(vol(j)/dt)*T_1(i,j,k)+((alfa*Ae)/de(j))*T_old(i-1,j,k)+((alfa*Aw)/dw(j))*T_old(i+1,j,k)&
                    +((alfa*Ad(j))/dz)*T_old(i,j,k-1)+((alfa*Au(j))/dz)*T_old(i,j,k+1)
                end do 

                a(1)=0;b(1)=1;c(1)=0;d(1)=1


                if(teta(i)>=0.and.teta(i)<pie/2) then
                    a(jmax)=-1;b(jmax)=1;c(jmax)=0;d(jmax)=0
                end if
                if(teta(i)>=pie/2.and.teta(i)<pie) then
                    a(jmax)=0;b(jmax)=1;c(jmax)=0;d(jmax)=3
                end if
                if(teta(i)>=pie.and.teta(i)<3.*pie/2.) then
                    a(jmax)=0;b(jmax)=1;c(jmax)=0;d(jmax)=2
                end if
                if(teta(i)>=3.*pie/2.and.teta(i)<2.*pie) then
                    a(jmax)=0;b(jmax)=1;c(jmax)=0;d(jmax)=1
                end if
                call TDMA_Solver(a,b,c,d,T,jmax)
                do j=1,jmax
                    T_old_1(i,j,k)=T(j)
                end do
            end do

            do j=2,jmax-1
                do i=2,imax-1
                  a(i)=-alfa*(Ae/de(i))
                    b(i)=(vol(i)/dt)+alfa*((Au(i)/dz)+(Ad(i)/dz)+(An(i)/dn)+(As(i)/ds)+(Ae/de(i))+(Aw/dw(i)))
                    c(i)=-alfa*(Aw/dw(i))
                    d(i)=(vol(i)/dt)*T_1(i,j,k)+((alfa*As(i))/ds)*T_old_1(i,j-1,k)+((alfa*An(i))/dn)*T_old_1(i,j+1,k)+&
                    ((alfa*Ad(i))/dz)*T_old_1(i,j,k-1)+((alfa*Au(i))/dz)*T_old_1(i,j,k+1)  
                end do 

                a(1)=-1.0
                b(1)=1.0
                c(1)=0.0
                d(1)=0.0
                a(imax)=0.0
                b(imax)=1.0
                c(imax)=-1.0
                d(imax)=0.0
                
                call PERIODIC2(imax,a,b,c,d,tpime)    
                do i=1,imax
                    T_old_1(i,j,k)=Tpime(i)
                end do
            end do 
        end do 
        T_old_1(1,:,:)=T_old_1(imax-1,:,:)
        T_old_1(imax,:,:)=T_old_1(2,:,:)
        
    ! j=1
                        nCount=1.0
                        tempAvgJ1=0.0
                        do i=1,imax
                            do k=1,kmax
                            nCount=nCount+1
                            tempAvgJ1=T_old_1(i,2,k)+tempAvgJ1
                            end do
                        end do
                        tempAvgJ1=tempAvgJ1/nCount
                        do i=1,imax
                            do k=1,kmax
                            T_old_1(i,1,k)=tempAvgJ1
                            T_1(i,1,k)=tempAvgJ1
                            end do 
                        end do
                        T_old_1(:,:,1)=0
                        T_1(:,:,1)=0
                        T_old_1(:,:,kmax)=1
                        T_1(:,:,kmax)=1

        !-----------------------------------------------------

        !---------------------------------------------
        !!!  end do R_sweep and teta_sweep  !!!
        !---------------------------------------------
        resuprime=0.0
       ! cprime=0.0                    
        do i=1,imax
            do j=1,jmax
                do k=1,kmax  
                    resuprime=resuprime+Abs(T_old_1(i,j,k)-T_old(i,j,k))
                  !  cprime=cprime+1
                end do
            end do
        end do
        ! resuprime=resuprime/cprime
        do i=1,imax
            do j=1,jmax
                do k=1,kmax  
                    T_old(i,j,k)=T_old_1(i,j,k)
                end do
            end do
        end do 

        if (resuprime<0.0001)exit
        print*,'resuprime',resuprime
        end do

        !----------------------------------------
        !!!  end do time step  !!!
        !----------------------------------------
        resut=0.0
        czegond=0.0
        do i=1,imax
            do j=1,jmax
                do k=1,kmax 
                    ! delta_1(i,j)=Ads(T_old_1(i,j)-T_1(i,j))
                    resut=resut+Abs(T_old_1(i,j,k)-T_1(i,j,k))
                   ! czegond=czegond+1	
                end do
            end do
        end do
        ! resut=resut/czegond
        do i=1,imax
            do j=1,jmax
                do k=1,kmax  
                    T_1(i,j,k)=T_old_1(i,j,k)

                end do
            end do
        end do
        print*,"counter=",counter

        if (resut<0.001) exit
        print*,'resut',resut
        
open(21,file="tecplot_result.plt",status="old",position="append",action="write")
write(21,*)'variables = "X" "Y" "Z" "T"'
write(21,*) 'Zone T="Time:',counter*dt,',Timeset:',dt,'"'
write(21,*) 'strandID=1, SolutionTime=',counter*dt
write(21,*) 'i=',imax-1, ', j=',jmax,', k=',kmax,', ZONETYPE=Ordered'
write(21,*) 'DATAPACKING=Point'

 do k=1,kmax
        do j=1,jmax
            do i=1,imax-1
                write(21,500)r(j)*cos(teta(i)),r(j)*sin(teta(i)),z(k),T_1(i,j,k)
            end do
        end do
 end do
 
close(21)

    end do

    print*,"Done"
    pause
    end program

    !------------------------------------------------
    !!!  subroutine TDMA_solver  !!!
    !------------------------------------------------
    subroutine TDMA_solver(ai,bi,ci,di,Ti,nn)
    integer,intent(in)::nn
    integer::j         !counter
    !integer,parameter::single=4              !Compiler dependept value
    !integer,parameter::double=8              !Compiler dependept value
    double precision,dimension(nn),intent(in)::ai,bi,ci,di
    double precision,dimension(nn),intent(out)::Ti
    ! double precision,dimension(nn),intent(out)::Ti
    double precision,dimension(nn)::d_a,b_a
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

    !---------------------------------------------
    !!! subroutine two over LAP !!!
    !---------------------------------------------

    subroutine PERIODIC2(ntot,aco,bco,cco,dco,tpime)
    double precision,dimension(ntot) :: aco,bco,cco,dco,eco,fco,tpime
    integer ntot
    !
    n2=ntot-2
    n1=ntot-1
    dn=0.
    !
    ! begin with the firdt row and dtart the elimination
    !                                                   
    eco(1)=aco(1)
    !
    ! now move to the second row
    !     
    fact=-aco(2)/bco(1)                      
    bco(2)=bco(2)+fact*cco(1)
    dco(2)=dco(2)+fact*dco(1)
    eco(2)=fact*eco(1)
    do 1 i=3,ntot-3      
        fact=-aco(i)/bco(i-1)
        bco(i)=bco(i)+fact*cco(i-1)
        dco(i)=dco(i)+fact*dco(i-1)
        eco(i)=fact*eco(i-1)
1   continue
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
    fco(2)=cco(ntot)
    fpn1=aco(ntot)
    fpn=bco(ntot)
    !
    ! now go across the final row
    !                            
    fact=-fco(2)/bco(2)
    fco(3)=fact*cco(2)
    fpn1=fpn1+fact*eco(2)
    dn=dn+fact*dco(2)
    do 2 i=4,n2
        fact=-fco(i-1)/bco(i-1)
        fco(i)=fact*cco(i-1)
        fpn1=fpn1+fact*eco(i-1)
        dn=dn+fact*dco(i-1)
2   continue
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
    ! dtart the back solve
    !      
    tpime(ntot)=b1/a1
    tpime(2)=tpime(ntot)
    tpime(n1)=(dn-fpn*tpime(ntot))/fpn1               
    tpime(1)=tpime(n1)
    tpime(n2)=(dco(n2)-eco(n2)*tpime(n1))/bco(n2)
    !
    ! now we are ready for a recursion relationship
    !          
    do 3 i1=3,n2-1
        i=ntot-i1
        tpime(i)=(dco(i)-eco(i)*tpime(n1)-cco(i)*tpime(i+1))/bco(i)         
3   continue  
    return
    end subroutine PERIODIC2
    