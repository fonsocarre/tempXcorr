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

                y1 = (corr(dispVec(1), x1))
                y2 = (corr(dispVec(1), x2))
                y3 = (corr(dispVec(1), x3))
                A = log(y3) - log(y1)
                B = log(y2) - log(y1)
                xdisp = (B*(x1**2 - x3**2) - A*(x1**2 - x2**2))/(2*A*(x2-x1) - 2*B*(x3-x1))

                x1 = dispVec(1) - 1
                x2 = dispVec(1)
                x3 = dispVec(1) + 1
                if (x1 <= 0) x1 = x3 + 1
                if (x3 > nRows) x3 = x1 -1

                y1 = (corr(x1, dispVec(2)))
                y2 = (corr(x2, dispVec(2)))
                y3 = (corr(x3, dispVec(2)))
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
        integer(c_int), intent(INOUT)     :: xdisp(xmax)
        integer(c_int), intent(INOUT)     :: ydisp(xmax) ! TODO ydisp not implemented
        real(c_double), intent(OUT)     :: corrVal(xmax)

        real(c_double)                  :: mat0(nRows, nCols)
        real(c_double)                  :: mat1(nRows, nCols)

        integer                         :: i,j

        mat0 = transpose(reshape(tmat0, [nCols, nRows]))
        mat1 = transpose(reshape(tmat1, [nCols, nRows]))

        mat0 = (max(mat0,offset) - offset)**expon
        mat1 = (max(mat1,offset) - offset)**expon

        do j=1, xmax
            xdisp(j) = j
            corrVal(j) = conv(mat0, mat1, nRows, nCols, j, 0_c_int)
            print*, xdisp(j), corrVal(j)
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
        !print*, "    jmaxR, jminR = ", jmaxR,jminR
        !print*, "    jmaxL,jminL = ", jmaxL,jminL
        !read(*,*)
        conv = sum(matL(iminL:imaxL, jminL:jmaxL)*&
                   matR(iminR:imaxR, jminR:jmaxR))
    end function conv

    function calculatexcorr(tmat0, tmat1, nRows, nCols, xdisp, expon, offset, jmin, x2)&
                                    result(xcorr) bind(c)
        real(c_double), intent(IN)          :: tmat0(nRows*nCols)
        real(c_double), intent(IN)          :: tmat1(nRows*nCols)
        integer(c_int), intent(IN)          :: nRows
        integer(c_int), intent(IN)          :: nCols
        integer(c_int), intent(IN)          :: xdisp
        integer(c_int), intent(IN)          :: expon
        real(c_double), intent(IN)          :: offset
        integer(c_int), intent(IN)          :: jmin
        integer(c_int), intent(IN)          :: x2

        real(c_double)                      :: xcorr

        real(c_double)                      :: mat0(nRows, nCols)
        real(c_double)                      :: mat1(nRows, nCols)

        real(c_double)                      :: mat00(nRows, nCols-jmin)
        real(c_double)                      :: mat11(nRows, nCols-jmin)
        
        integer                             :: nx, ny

        mat0 = transpose(reshape(tmat0, [nCols, nRows]))
        mat1 = transpose(reshape(tmat1, [nCols, nRows]))

        !mat0 = (max(mat0,offset) - offset)**expon
        !mat1 = (max(mat1,offset) - offset)**expon
        !print*, "nCols = ", nCols
        !print*, "size(mat0(1,:))=", size(mat0(1,:))
        mat0 = mat0**expon
        mat1 = mat1**expon
               

        mat00 = mat0(:,1:nCols-jmin)
        mat11 = mat1(:,1:nCols-jmin)
        
        mat00 = mat00(:,nCols-jmin:1:-1)
        mat11 = mat11(:,nCols-jmin:1:-1)
        
        
        nx = size(mat00(1,:))
        ny = size(mat00(:,1))
        
        !print*, "sxcorr"
        !xcorr = conv(mat00, mat11, ny, nx, xdisp, 0)
        xcorr = conv(mat11(x2:x2,:), mat00(x2:x2,:), 1, nx, xdisp, 0)

    end function calculatexcorr

    subroutine PIV_xCorr(mat0, mat1, nRows, nCols, wNx, wNy, &
                         xdisps, ydisps, xcorrvals) bind(c)
        real(c_double), intent(IN)              :: mat0(nRows*nCols)
        real(c_double), intent(IN)              :: mat1(nRows*nCols)
        integer(c_int), intent(IN)              :: nRows
        integer(c_int), intent(IN)              :: nCols
        integer(c_int), intent(IN)              :: wNx
        integer(c_int), intent(IN)              :: wNy
        real(c_double), intent(OUT)             :: xdisps(nRows*nCols)
        real(c_double), intent(OUT)             :: ydisps(nRows*nCols)
        real(c_double), intent(OUT)             :: xcorrvals(nRows*nCols)

        real(c_double)                          :: tmat0(nRows, nCols)
        real(c_double)                          :: tmat1(nRows, nCols)
        real(c_double)                          :: txdisps(nRows, nCols)
        real(c_double)                          :: tydisps(nRows, nCols)
        real(c_double)                          :: txcorrvals(nRows, nCols)
        integer                                 :: nx, ny
        integer                                 :: i, j

        tmat0 = transpose(reshape(mat0, [nCols, nRows]))
        tmat1 = transpose(reshape(mat1, [nCols, nRows]))

        do j=1, nCols
            do i=1, nRows
                call point_mat_xcorr([i,j], tmat0(i,j), tmat1, wNx, wNy,&
                                txdisps(i,j), tydisps(i,j), txcorrvals(i,j))
            end do
        end do

        xdisps = reshape(transpose(txdisps), [nRows*nCols])
        ydisps = reshape(transpose(tydisps), [nRows*nCols])
        xcorrvals = reshape(transpose(txcorrvals), [nRows*nCols])
        
        !xdisps = reshape(txdisps, [nRows*nCols])
        !ydisps = reshape(tydisps, [nRows*nCols])
        !xcorrvals = reshape(txcorrvals, [nRows*nCols])
    end subroutine PIV_xCorr


    subroutine point_mat_xcorr(ijPoint, valPoint, mat, wNx, wNy, xdisp, ydisp, xcorrval)
        integer(c_int), intent(IN)              :: ijPoint(:)
        real(c_double), intent(IN)              :: valPoint
        real(c_double), intent(IN)              :: mat(:,:)
        integer(c_int), intent(IN)              :: wNx
        integer(c_int), intent(IN)              :: wNy
        real(c_double), intent(OUT)             :: xdisp
        real(c_double), intent(OUT)             :: ydisp
        real(c_double), intent(OUT)             :: xcorrval
        
        integer                                 :: i, j
        integer                                 :: iimin, iimax
        integer                                 :: jjmin, jjmax
        
        integer                                 :: nRows, nCols
        real(c_double), allocatable             :: corrMat(:,:)
        integer(c_int)                          :: disps(2)

        nRows = size(mat(:,1))
        nCols = size(mat(1,:))

        xdisp = 0.0d0
        ydisp = 0.0d0
        xcorrval = 0.0
    
        if (ijPoint(1) - wNy >= 1) then
            iimin = -wNy
        else
            !iimin = ijPoint(1) - wNy
            iimin = -ijPoint(1) + 1
        end if

        if (ijPoint(2) - wNx >= 1) then
            jjmin = -wNx
        else
            jjmin = -ijPoint(2) + 1
        end if

        if (nCols - ijPoint(1) < wNy) then
            iimax = nCols - ijPoint(1) - 1
        else
            iimax =  wNy
        end if

        if (nRows - ijPoint(2) < wNx) then
            jjmax = nRows - ijPoint(2) - 1
        else
            jjmax =  wNx
        end if

        allocate(corrMat(iimax-iimin+1, jjmax-jjmin+1))
        corrMat = 0.0d0

        corrMat = valPoint*mat(ijPoint(1)+iimin:ijPoint(1)+iimax,&
                              ijPoint(2)+jjmin:ijPoint(2)+jjmax)

        !do j=1, jjmax - jjmin
            !do i=1, iimax-iimin
                    !print*, "ijPoint(1) = " , ijPoint(1)
                    !print*, "i = ", i
                    !print*, "iimin = ", iimin
                    !print*, "ijPoint(2) = " , ijPoint(2)
                    !print*, "j = ", j
                    !print*, "jjmin = ", jjmin
                !print*, "iindex = ", i-iimin
                !print*, "jindex = ", j-jjmin
                !read(*,*)
                !corrMat(i,j) = valPoint*mat(i-iimin,&
                                            !j-jjmin)
            !end do
        !end do

        xcorrval = maxval(corrMat)
        disps = maxloc(corrMat)
        
        xdisp = disps(1) - iimin + 1
        ydisp = disps(2) - jjmin + 1


        deallocate(corrMat)
    end subroutine point_mat_xcorr
end module xCorrModule















