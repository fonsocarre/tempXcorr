!
!   interpolationModule.f90
!   Stitchutils
!
!   Created by Alfonso del Carre on 08/06/2015.
!   Copyright 2015 Alfonso del Carre. All rights reserved.
!
module interpolationModule
use iso_c_binding
use IEEE_EXCEPTIONS
use IEEE_ARITHMETIC
implicit none

contains

subroutine interpolate(ixinL, iyinL, ivarin1L, ivarin2L, ixinR,&
                       iyinR, ivarin1R, ivarin2R,&
                       nRowsin, nColsin, ixout, iyout, varout1, varout2,&
                       nRowsout, nColsout, xdisp, ydisp) bind(c)
    real(c_double), intent(IN)          :: ixinL(nRowsin*nColsin)
    real(c_double), intent(IN)          :: iyinL(nRowsin*nColsin)
    real(c_double), intent(IN)          :: ivarin1L(nRowsin*nColsin)
    real(c_double), intent(IN)          :: ivarin2L(nRowsin*nColsin)
    real(c_double), intent(IN)          :: ixinR(nRowsin*nColsin)
    real(c_double), intent(IN)          :: iyinR(nRowsin*nColsin)
    real(c_double), intent(IN)          :: ivarin1R(nRowsin*nColsin)
    real(c_double), intent(IN)          :: ivarin2R(nRowsin*nColsin)
    integer(c_int), intent(IN)          :: nRowsin
    integer(c_int), intent(IN)          :: nColsin
    real(c_double), intent(IN)          :: ixout(nRowsout*nColsout)
    real(c_double), intent(IN)          :: iyout(nRowsout*nColsout)
    integer(c_int), intent(IN)          :: nRowsout
    integer(c_int), intent(IN)          :: nColsout
    real(c_double), intent(IN)          :: xdisp
    real(c_double), intent(IN)          :: ydisp

    ! varout will be a 1D vector by rows
    real(c_double), intent(INOUT)       :: varout1(nRowsout*nColsout,1)
    real(c_double), intent(INOUT)       :: varout2(nRowsout*nColsout,1)

    real(c_double), allocatable         :: xinL(:,:)
    real(c_double), allocatable         :: yinL(:,:)
    real(c_double), allocatable         :: varin1L(:, :)
    real(c_double), allocatable         :: varin2L(:, :)
    real(c_double), allocatable         :: xinR(:, :)
    real(c_double), allocatable         :: yinR(:, :)
    real(c_double), allocatable         :: varin1R(:, :)
    real(c_double), allocatable         :: varin2R(:, :)
    real(c_double), allocatable         :: xout(:, :)
    real(c_double), allocatable         :: yout(:, :)

    real(c_double), allocatable         :: tvarout1(:, :)
    real(c_double), allocatable         :: tvarout2(:, :)
    real(c_double), allocatable         :: ttvarout1(:, :)
    real(c_double), allocatable         :: ttvarout2(:, :)


    integer                             :: i
    integer                             :: j

    real(c_double)                      :: xminL, xmaxL, xminR, xmaxR
    real(c_double)                      :: yminL, ymaxL, yminR, ymaxR
    real(c_double)                      :: dxL, dxR, dxO
    real(c_double), parameter           :: eps = 1.0e-4

    integer                             :: jmaxL, jminR
    integer                             :: jminRout, jmaxLout
    real(c_double)                      :: xyout(2)

    ! problem with heap (or stack?) size in fortran
    allocate(xinL(nRowsin, nColsin))
    allocate(yinL(nRowsin, nColsin))
    allocate(varin1L(nRowsin, nColsin))
    allocate(varin2L(nRowsin, nColsin))
    allocate(xinR(nRowsin, nColsin))
    allocate(yinR(nRowsin, nColsin))
    allocate(varin1R(nRowsin, nColsin))
    allocate(varin2R(nRowsin, nColsin))
    allocate(xout(nRowsout, nColsout))
    allocate(yout(nRowsout, nColsout))

    allocate(tvarout1(nRowsout, nColsout))
    allocate(tvarout2(nRowsout, nColsout))

    xinL = transpose(reshape(ixinL, [nColsin, nRowsin]))
    yinL = transpose(reshape(iyinL, [nColsin, nRowsin]))
    varin1L = transpose(reshape(ivarin1L, [nColsin, nRowsin]))
    varin2L = transpose(reshape(ivarin2L, [nColsin, nRowsin]))
    xinR = transpose(reshape(ixinR, [nColsin, nRowsin]))
    yinR = transpose(reshape(iyinR, [nColsin, nRowsin]))
    varin1R = transpose(reshape(ivarin1R, [nColsin, nRowsin]))
    varin2R = transpose(reshape(ivarin2R, [nColsin, nRowsin]))
    xout = transpose(reshape(ixout, [nColsout, nRowsout]))
    yout = transpose(reshape(iyout, [nColsout, nRowsout]))

    dxL = abs(xinL(1,2) - xinL(1,1))
    dxR = abs(xinR(1,2) - xinR(1,1))
    dxO = abs(xout(1,2) - xout(1,1))

    if (abs(dxL-dxR) >= eps) then
        print*, "dx are different, check algorithm"
        stop
    end if
    if (abs(dxO-dxR) >= eps) then
        print*, "dx are different, check algorithm"
        stop
    end if

    jmaxL = 0
    jminR = 0
    jminRout = 0
    jmaxLout = 0
    do j=1, nColsin
        if (xinL(1,j) < minval(xinR)) jmaxL = j
    end do
    do j=1, nColsin
        if (xinR(1,j) > maxval(xinL)) then
            jminR = j
            exit
        end if
    end do
    do j=1, nColsout
        if (xout(1,j) > xinR(1,jminR)) then
            jminRout = j 
            exit
        end if
    end do
    do j=1, nColsout
        if (xout(1,j) > xinL(1,jmaxL)) then
            jmaxLout = j - 1
            exit
        end if
    end do

    do i=1, nRowsout
        do j=1, nColsout
        if (j <= nColsin) then
            !if (j <= jmaxL) then
            ! only first frame
                !if (abs(xinL(i,j) - xout(i,j)) < eps) then
                ! same initial grid
                    tvarout1(i,j) = varin1L(i,j)
                    tvarout2(i,j) = varin2L(i,j)
                !end if
            !end if
        end if
        end do
    end do



!do j=1, nColsOut
!    if (abs(xout(1,j) - 95.788) < 1.0e-3) then
!        print*, "  sdv"
!    end if    
!end do









! right frame
!call transferMesh1(xinR(:,jminR-3:), yinR(:,jminR-3:),&
!                   varin1R(:,jminR-3:), varin2R(:,jminR-3:),&
!                   xout(:,jminRout-2:), yout(:,jminRout-2:),&
!                   tvarout1(:,jminRout-2:), tvarout2(:,jminRout-2:))


call middleInterpolation(xinR(:,1:jminR+1),&
                         yinR(:,1:jminR+1),&
                         xinL(:,jmaxL-2:),&
                         yinL(:,jmaxL-2:),&
                         varin1R(:,1:jminR+1),&
                         varin1L(:,jmaxL-2:),&
                         varin2R(:,1:jminR+1),&
                         varin2L(:,jmaxL-2:),&
                         xout(:,jmaxLout-1:jminRout+2),&
                         yout(:,jmaxLout-1:jminRout+2),&
                         tvarout1(:,jmaxLout-1:jminRout+2),&
                         tvarout2(:,jmaxLout-1:jminRout+2))

call transferMesh1(xinR, yinR,&
                   varin1R, varin2R,&
                   xout(:,jminRout-2:), yout(:,jminRout-2:),&
                   tvarout1(:,jminRout-2:), tvarout2(:,jminRout-2:))



!call transferMesh1(xinR, yinR,&
!varin1R, varin2R,&
!xout(:,jminRout:), yout(:,jminRout:),&
!tvarout1(:,jminRout:), tvarout2(:,jminRout:))


!! middle zone, two transferMesh1 and then weighting with the distance
!print*, "Transf 2"
!call transferMesh1(xinR(:,1:jminR+1), yinR(:,1:jminR+1),&
!                    varin1R(:,1:jminR+1), varin2R(:,1:jminR+1),&
!                    xout(:,jmaxLout:jminRout-1), yout(:,jmaxLout:jminRout-1),&
!                    tvarout1(:,jmaxLout:jminRout-1), tvarout2(:,jmaxLout:jminRout-1))
!
!!call transferMesh1(xinR, yinR,&
!!varin1R, varin2R,&
!!xout(:,jmaxLout-1:jminRout), yout(:,jmaxLout-1:jminRout),&
!!tvarout1(:,jmaxLout-1:jminRout), tvarout2(:,jmaxLout-1:jminRout))
!
!allocate(ttvarout1(size(tvarout1(:,1)), size(tvarout1(1,:))))
!allocate(ttvarout2(size(tvarout1(:,1)), size(tvarout1(1,:))))
!ttvarout1 = 0.0d0
!ttvarout2 = 0.0d0
!
!print*, "maxval(xin) ", maxval(xinL)
!print*, "maxval(xout(:,jmaxLout:jminRout-1)) ", maxval(xout(:,jmaxLout:jminRout-1))
!
!print*, "Transf 3"
!call transferMesh1(xinL(:,jmaxL-1:), yinL(:,jmaxL-1:),&
!                varin1L(:,jmaxL-1:), varin2L(:,jmaxL-1:),&
!                xout(:,jmaxLout:jminRout-1), yout(:,jmaxLout:jminRout-1),&
!                ttvarout1(:,jmaxLout:jminRout-1), ttvarout2(:,jmaxLout:jminRout-1))

!call transferMesh1(xinL, yinL,&
!varin1L, varin2L,&
!xout(:,jmaxLout:jminRout-1), yout(:,jmaxLout:jminRout-1),&
!ttvarout1(:,jmaxLout:jminRout-1), ttvarout2(:,jmaxLout:jminRout-1))
!
!tvarout1(:,jmaxLout:jminRout-1) = (tvarout1(:,jmaxLout:jminRout-1) + ttvarout1)*0.5
!tvarout2(:,jmaxLout:jminRout-1) = (tvarout2(:,jmaxLout:jminRout-1) + ttvarout2)*0.5

    deallocate(xinL)
    deallocate(yinL)
    deallocate(varin1L)
    deallocate(varin2L)
    deallocate(xinR)
    deallocate(yinR)
    deallocate(varin1R)
    deallocate(varin2R)
    deallocate(xout)
    deallocate(yout)



    varout1 = reshape(transpose(tvarout1), [nRowsout*nColsOut,1])
    varout2 = reshape(transpose(tvarout2), [nRowsout*nColsOut,1 ])
    deallocate(tvarout1)
    deallocate(tvarout2)
end subroutine interpolate

subroutine transferMesh1(xin, yin, varin1, varin2, xout, yout, varout1, varout2)
    real(c_double), intent(IN)          :: xin(:,:)
    real(c_double), intent(IN)          :: yin(:,:)
    real(c_double), intent(IN)          :: varin1(:,:)
    real(c_double), intent(IN)          :: varin2(:,:)
    real(c_double), intent(IN)          :: xout(:,:)
    real(c_double), intent(IN)          :: yout(:,:)
    real(c_double), intent(OUT)         :: varout1(:,:)
    real(c_double), intent(OUT)         :: varout2(:,:)

    integer                             :: ijfirst(2)
    integer                             :: ijlast(2)
    integer                             :: i,j
    integer                             :: iin, jin
    integer                             :: temp
    integer                             :: nColsin, nRowsin
    integer                             :: nColsout, nRowsout

    real(c_double)                      :: coeffVec(4)
    real(c_double)                      :: valueVec(4)
    real(c_double)                      :: xyout(2)

    real(c_double), parameter           :: eps = 1.0e-6

    nColsin = size(xin(1,:))
    nRowsin = size(yin(:,1))
    nColsout = size(xout(1,:))
    nRowsout = size(yout(:,1))

    if (size(varout1(1,:)) /= nColsout .OR. size(varout1(:,1)) /= nRowsout) then
        print*, "Sizes not consistent"
        stop
    end if


!    print*, "xin(1,1) = ", xin(1,1)
!    print*, "xin(1,2) = ", xin(1,2)
!    print*, "xout(1,1) = ", xout(1,1)
    ! find first square that corresponds
    ijfirst = 0
    ijlast = 0
    do j=1, nColsin-1
        if (ijfirst(2) == 0) then
            if (xin(1,j) < xout(1,1) .AND. xin(1,j+1) > xout(1,1)) then
                ijfirst(2) = j
            end if
            if (j == nColsin-1) ijfirst(2) = -1
        end if
    end do
    do i=1, nRowsin-1
        if (ijfirst(1) == 0) then
            if (yin(i,1) > yout(1,1) .AND. yin(i+1,1) < yout(1,1)) then
                ! ijfirst(1) will be -1 if the out mesh is higher than the in
                ijfirst(1) = i
            end if
            if (i == nRowsin-1) ijfirst(1) = -1
        end if
    end do

do j=nColsin, 2, -1
    if (ijlast(2) == 0) then
        if (xin(1,j) > xout(1,nColsOut) .AND. xin(1,j-1) < xout(1,nColsOut)) then
            ijlast(2) = j
        end if
        if (j == 2) ijlast(2) = -1
    end if
end do
do i=nRowsin, 2, -1
    if (ijlast(1) == 0) then
        if (yin(i,1) < yout(nRowsout, 1) .AND. yin(i-1, 1) > yout(1,nRowsout)) then
            ijlast(1) = i
        end if
        if (i == 2) ijlast(1) = -1
    end if
end do

!    print*, "yin(ijfirst,1) = ", yin(ijfirst(1),1)
!    print*, "yin(ijfirst+1,1) = ", yin(ijfirst(1)+1,1)
!    print*, "yout(1,1) = ", yout(1,1)

    if (minval(ijfirst) == 0) then
        print*, "ERROR, ijfirst not found in interpolationModule:transferMesh1"
        stop
    end if


    do j=1, nColsout-1
        jin = j + ijfirst(2) - 1
        do i=1, nRowsout
            iin = i + ijfirst(1) - 1
            xyout = [xout(i,j), yout(i,j)]
            coeffVec(1) = dist(xyout, [xin(iin, jin), yin(iin, jin)])
            coeffVec(2) = dist(xyout, [xin(iin, jin+1), yin(iin, jin+1)])
            coeffVec(3) = dist(xyout, [xin(iin+1, jin+1), yin(iin+1, jin+1)])
            coeffVec(4) = dist(xyout, [xin(iin+1, jin), yin(iin+1, jin)])

            ! normalisation
            coeffVec = 1.0d0/coeffVec
            coeffVec = coeffVec/sum(coeffVec)

            ! first var
            valueVec(1) = varin1(iin, jin)
            valueVec(2) = varin1(iin, jin+1)
            valueVec(3) = varin1(iin+1, jin+1)
            valueVec(4) = varin1(iin+1, jin)

            varout1(i,j) = dot_product(coeffVec, valueVec)

            ! second var
            valueVec(1) = varin2(iin, jin)
            valueVec(2) = varin2(iin, jin+1)
            valueVec(3) = varin2(iin+1, jin+1)
            valueVec(4) = varin2(iin+1, jin)

            varout2(i,j) = dot_product(coeffVec, valueVec)
        end do
    end do


end subroutine transferMesh1






subroutine middleInterpolation(xinR,&
                                 yinR,&
                                 xinL,&
                                 yinL,&
                                 varin1R,&
                                 varin1L,&
                                 varin2R,&
                                 varin2L,&
                                 xout,&
                                 yout,&
                                 tvarout1,&
                                 tvarout2)

    real(c_double), intent(IN)              :: xinR(:,:)
    real(c_double), intent(IN)              :: yinR(:,:)
    real(c_double), intent(IN)              :: xinL(:,:)
    real(c_double), intent(IN)              :: yinL(:,:)
    real(c_double), intent(IN)              :: varin1R(:,:)
    real(c_double), intent(IN)              :: varin1L(:,:)
    real(c_double), intent(IN)              :: varin2R(:,:)
    real(c_double), intent(IN)              :: varin2L(:,:)
    real(c_double), intent(IN)              :: xout(:,:)
    real(c_double), intent(IN)              :: yout(:,:)
    real(c_double), intent(INOUT)           :: tvarout1(:,:)
    real(c_double), intent(INOUT)           :: tvarout2(:,:)

    integer                                 :: i, j, ii, jj, irad
    integer                                 :: nColsin, nRowsin
    integer                                 :: nColsout, nRowsout
    real(c_double)                          :: minRad
    real(c_double)                          :: incRad
    integer, parameter                      :: maxIncRad = 8
    integer, parameter                      :: minPoints = 4
    integer, parameter                      :: maxPoints = 4
    real(c_double)                          :: radius
    real(c_double)                          :: dx
    integer                                 :: counterL, counterR
    real(c_double)                          :: distL(maxPoints/2), distR(maxPoints/2)
    real(c_double)                          :: value1L(maxPoints/2), value2L(maxPoints/2)
    real(c_double)                          :: value1R(maxPoints/2), value2R(maxPoints/2)
    real(c_double)                          :: factor, xminR, xmaxL, val1R, val1L, val2L, val2R
    real(c_double), parameter               :: eps = 1.0e-5
    real(c_double)                          :: temp
    logical                                 :: only1point


    nColsin = size(xinR(1,:))
    nRowsin = size(yinR(:,1))
    nColsout = size(xout(1,:))
    nRowsout = size(yout(:,1))

    xminR = minval(xinR)
    xmaxL = maxval(xinL)

    dx = abs(xout(1,2) - xout(1,1))
    minRad = 0.1*dx
    incRad = 0.5*dx

    do j=1, nColsout
        do i=1, nRowsout
            counterL = 0
            counterR = 0
            distL = 0
            distR = 0
            value1L = 0.0d0
            value2L = 0.0d0
            value1R = 0.0d0
            value2R = 0.0d0
            do irad=1, maxIncRad
                if (counterL > minPoints/2 .AND. counterR > minPoints/2) exit
                radius = minRad + (irad-1)*incRad
                do jj=1, nColsin
                    do ii=1, nRowsin
                        ! left frame
                        if (counterL < maxPoints/2) then
                            if (isInside(xout(i,j), yout(i,j),&
                                         xinL(ii,jj), yinL(ii,jj), radius)) then
                                counterL = counterL + 1
                                if (counterL <= maxPoints/2) then
                                    distL(counterL) = dist([xout(i,j), yout(i,j)],&
                                                          [xinL(ii,jj), yinL(ii,jj)])
                                    value1L(counterL) = varin1L(ii,jj)
                                    value2L(counterL) = varin2L(ii,jj)
                                end if
                            end if
                        end if

                        ! right frame
                        if (counterR < maxPoints/2) then
                            if (isInside(xout(i,j), yout(i,j),&
                                         xinR(ii,jj), yinR(ii,jj), radius)) then
                                counterR = counterR + 1
                                if (counterR <= maxPoints/2) then
                                    distR(counterR) = dist([xout(i,j), yout(i,j)],&
                                                          [xinR(ii,jj), yinR(ii,jj)])
                                    value1R(counterR) = varin1R(ii,jj)
                                    value2R(counterR) = varin2R(ii,jj)
                                end if
                            end if
                        end if
                    end do
                end do
            end do

            only1point = .false.
            temp = 0.0d0
            do ii=1, size(distL)
                if (distL(ii) > eps) then
                    ! DEBUG
                    temp = temp + distL(ii)

                else if (distL(ii) > 0.0d0) then
                    val1L = checkNaN(value1l(ii))
                    val2L = checkNaN(value2l(ii))
                    only1point = .true.
                    exit
                end if
            end do
    
            if (.not. only1point) then
                temp = sum(distR)
                if (temp == 0.0d0) then
                    stop "Error"
                end if
                distL = distL / temp

                val1L = checkNaN(dot_product(distL, value1L))
                val2L = checkNaN(dot_product(distL, value2L))
            end if

            only1point = .false.
            temp = 0.0d0
            do ii=1, size(distR)
                if (distR(ii) > eps) then
                    ! DEBUG
                    temp = temp + distR(ii)
                else  if (distR(ii) > 0.0d0) then
                    val1R = checkNaN(value1R(ii))
                    val2R = checkNaN(value2R(ii))
                    only1point = .true.
                    exit
                end if
            end do
            if (.not. only1point) then
                temp = sum(distR)
                if (temp == 0.0d0) then
                    stop "Error"
                end if
                !print*, "temp = ", temp
                distR = distR / temp
                !print*, "newdist = ", distR

                val1R = checkNaN(dot_product(distR, value1R))
                val2R = checkNaN(dot_product(distR, value2R))
            end if   

            factor = 1.0d0/(xmaxL - xminR) * (xout(i,j) - xminR)
            if (factor > 1) factor = 1.0d0
            if (factor < 0) factor = 0.0d0

            tvarout1(i,j) = (1.0d0 - factor)*val1L + factor*val1R
            tvarout2(i,j) = (1.0d0 - factor)*val2L + factor*val2R
        end do
    end do

end subroutine middleInterpolation

function isInside(xc, yc, x, y, radius)
    real(c_double), intent(IN)      :: xc
    real(c_double), intent(IN)      :: yc
    real(c_double), intent(IN)      :: x
    real(c_double), intent(IN)      :: y
    real(c_double), intent(IN)      :: radius

    logical                         :: isInside

    isInside = .false.
    if ((xc-x)**2 + (yc-y)**2 < radius**2) isInside = .true.
end function isInside


function dist(vec1, vec2)
real(c_double), intent(IN)      :: vec1(:)
real(c_double), intent(IN)      :: vec2(:)

real(c_double)                  :: dist

    dist = norm2(vec2(1)-vec1(1), vec2(2)-vec1(2))
end function

function norm2(x,y)
    real(c_double), intent(IN)      :: x
    real(c_double), intent(IN)      :: y
    real(c_double)                  :: norm2

    norm2 = sqrt(x*x + y*y)
end function

function checkNaN(num)
    real(c_double), intent(IN)          :: num
    real(c_double)                      :: checkNaN

    if (ieee_is_nan(num)) then
        print*, "NaN appeared!!!"
        stop
    end if
    checkNaN = num
end function checkNaN

end module interpolationModule











