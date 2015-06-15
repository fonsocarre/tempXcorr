module xCorrModule
    use iso_c_binding
    implicit none
    
    contains
    
    subroutine xcorr(tmatL, tmatR, nRows, nCols, overlap, xdisp, ydisp,&
                            withInterpolation) bind(c)
        integer(c_int), intent(IN)      :: nRows
        integer(c_int), intent(IN)      :: nCols
        real(c_double), intent(IN)      :: tmatL(nRows*nCols)
        real(c_double), intent(IN)      :: tmatR(nRows*nCols)

        integer(c_int), intent(IN)      :: overlap
        real(c_double), intent(OUT)     :: xdisp
        real(c_double), intent(OUT)     :: ydisp
        logical(c_bool), intent(IN)     :: withInterpolation

        real(c_double)                  :: auxL(nRows, overlap)
        real(c_double)                  :: auxR(nRows, overlap)
        
        real(c_double)                  :: corr(nRows, overlap)
        
        integer                         :: i, ii
        integer                         :: j, jj
        integer                         :: imin, imax, jmin, jmax
        integer                         :: nR2
        integer                         :: dispVec(2)

        real(c_double)                  :: matL(nRows,nCols)
        real(c_double)                  :: matR(nRows,nCols)
        real(c_double)                  :: y1, y2, y3, x1, x2, x3, D

        real(c_double)          :: offset, A, B

        matL = transpose(reshape(tmatL, [nCols, nRows]))
        matR = transpose(reshape(tmatR, [nCols, nRows]))


        xdisp = 0.0d0
        ydisp = 0.0d0

        offset = 11.80d0

        auxL = max(matL(:, nCols-overlap+1:nCols),offset) - offset
        auxR = max(matR(:, 1:overlap), offset) - offset
        
        nR2 = int(nRows/2)
        ! Vertical disp
        ! range of motion: -nR2+1 -> nR2
        do i=-nR2+1, nR2
            ! Horizontal disp
            ! range of motion: 1 -> overlap
            do j=1, overlap
                corr(i+nR2, j) = conv(auxL, auxR, nRows, overlap, j-1, i)
            end do
        end do

        dispVec = maxloc(corr)
        xdisp =  dispVec(2) - overlap
        ydisp = -(dispVec(1)- (nR2)) ! CHANGE was nR2+1

        if (withInterpolation) then
                ! subpixel interpolation
                x1 = dispVec(2) - 1
                x2 = dispVec(2)
                x3 = dispVec(2) + 1
                if (x1 <= 0) x1 = x3 + 1
                if (x3 > overlap) x3 = x1 -1

                y1 = log(corr(dispVec(1), x1)+offset)
                y2 = log(corr(dispVec(1), x2)+offset)
                y3 = log(corr(dispVec(1), x3)+offset)
                A = log(y3) - log(y1)
                B = log(y2) - log(y1)
                xdisp = (B*(x1**2 - x3**2) - A*(x1**2 - x2**2))/(2*A*(x2-x1) - 2*B*(x3-x1))

                x1 = dispVec(1) - 1
                x2 = dispVec(1)
                x3 = dispVec(1) + 1
                if (x1 <= 0) x1 = x3 + 1
                if (x3 > nRows) x3 = x1 -1

                y1 = log(corr(x1, dispVec(2))+offset)
                y2 = log(corr(x2, dispVec(2))+offset)
                y3 = log(corr(x3, dispVec(2))+offset)
                A = log(y3) - log(y1)
                B = log(y2) - log(y1)
                ydisp = (B*(x1**2 - x3**2) - A*(x1**2 - x2**2))/(2*A*(x2-x1) - 2*B*(x3-x1))

                xdisp =  xdisp - overlap
                ydisp = -(ydisp - (nR2))! CHANGE was nR2+1 (-1 for sawtooth, 
                                            ! but it gives problems with cylinder
        end if
    end subroutine xCorr
    
    subroutine tempXcorr(tmat0, tmat1, xmin, xmax, ymin, ymax, offset, expon,&
                            nRows, nCols, xdisp, ydisp, corrVal) bind(c)
        integer(c_int), intent(IN)      :: nRows
        integer(c_int), intent(IN)      :: nCols
        integer(c_int), intent(IN)      :: xmin
        integer(c_int), intent(IN)      :: xmax
        integer(c_int), intent(IN)      :: ymin
        integer(c_int), intent(IN)      :: ymax
        real(c_double), intent(IN)      :: offset
        integer(c_int), intent(IN)      :: expon
        real(c_double), intent(IN)      :: tmat0(nRows*nCols)
        real(c_double), intent(IN)      :: tmat1(nRows*nCols)
        integer(c_int), intent(OUT)     :: xdisp(xmax)
        integer(c_int), intent(OUT)     :: ydisp(xmax) ! TODO ydisp not implemented
        real(c_double), intent(OUT)     :: corrVal(xmax)

        real(c_double)                  :: mat0(nRows, nCols)
        real(c_double)                  :: mat1(nRows, nCols)

        integer                         :: i,j

        mat0 = transpose(reshape(tmat0, [nCols, nRows]))
        mat1 = transpose(reshape(tmat1, [nCols, nRows]))

        mat0 = (max(mat0,offset) - offset)**expon
        mat1 = (max(mat1,offset) - offset)**expon

        do j=1, xmax
            xdisp(i) = i
            corrVal = conv(mat0, mat1, nRows, nCols, i, 0_c_int)
        end do

    end subroutine tempXcorr

    function conv(matL, matR, nRows, nCols, xdisp, ydisp)
        real(c_double)              :: conv
        real(c_double), intent(IN)  :: matL(:,:)
        real(c_double), intent(IN)  :: matR(:,:)
        integer(c_int), intent(IN)  :: nRows
        integer(c_int), intent(IN)  :: nCols
        integer(c_int), intent(IN)  :: xdisp
        integer(c_int), intent(IN)  :: ydisp

        integer                     :: iminL, imaxL
        integer                     :: jminL, jmaxL
        integer                     :: iminR, imaxR
        integer                     :: jminR, jmaxR
        
        iminR = max(1,-(ydisp))
        imaxR = min(nRows, nRows-ydisp+1)
        jminR = max(1, -xdisp)
        jmaxR = min(nCols, nCols-xdisp+1)

        iminL = max(1, (ydisp))
        imaxL = min(nRows, nRows+ydisp+1)
        jminL = max(1, xdisp)
        jmaxL = min(nCols, nCols+xdisp+1)

        if ((imaxR-iminR /= imaxL-iminL) .OR. &
            (jmaxR-jminR /= jmaxL-jminL)) then
            print*, "ERROR, sizes are not the same in conv"
            print*, "    imaxR-iminR = ", imaxR-iminR
            print*, "    imaxL-iminL = ", imaxL-iminL
            print*, "    jmaxR-jminR = ", jmaxR-jminR
            print*, "    jmaxL-jminL = ", jmaxL-jminL
            stop
        end if

        conv = sum(matL(iminL:imaxL, jminL:jmaxL)*&
                   matR(iminR:imaxR, jminR:jmaxR))
    end function conv

end module xCorrModule
