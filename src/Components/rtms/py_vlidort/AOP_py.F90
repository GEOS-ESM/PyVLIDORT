subroutine AOP_accum(n, nrh, p, ang, rh, drh, q_mass, &
                     lut_bext, lut_bsca, lut_g, lut_pmatrix, &
                     do_vector, aot, sca, g_acc, pmat_acc)
    implicit None
    
    ! !INPUT PARAMETERS:
    integer,          intent(in)            :: n
    integer,          intent(in)            :: nrh
    integer,          intent(in)            :: p
    integer,          intent(in)            :: ang
    integer,          intent(in)            :: do_vector
    
    real*8,           intent(in)            :: drh
    real*8,           intent(in)            :: rh(n)
    real*8,           intent(in)            :: q_mass(n)
    real*8,           intent(in)            :: lut_bext(nrh)
    real*8,           intent(in)            :: lut_bsca(nrh)
    real*8,           intent(in)            :: lut_g(nrh)
    real*8,           intent(in)            :: lut_pmatrix(p,ang,nrh)
    
    ! !OUTPUT PARAMETERS:
    real*8,           intent(inout)         :: aot(n)
    real*8,           intent(inout)         :: sca(n)
    real*8,           intent(inout)         :: g_acc(n)
    real*8,           intent(inout)         :: pmat_acc(p,ang,n)
    
    ! !LOCAL VARIABLES:
    integer                                 :: i, jp, ka, irh
    real*8                                  :: aot_i, sca_i, q
    
    do i = 1, n
        irh = int(rh(i) / drh + 0.5d0) + 1
        if (irh < 1) irh = 1
        if (irh > nrh) irh = nrh
        
        q = q_mass(i)
        
        aot_i = lut_bext(irh) * q
        sca_i = lut_bsca(irh) * q
        
        aot(i) = aot(i) + aot_i
        sca(i) = sca(i) + sca_i
        g_acc(i) = g_acc(i) + sca_i * lut_g(irh)
        
        if (do_vector == 1) then
            do ka = 1, ang
                do jp = 1, p
                    pmat_acc(jp, ka, i) = pmat_acc(jp, ka, i) + &
                                          lut_pmatrix(jp, ka, irh) * sca_i
                end do
            end do
        end if
        
    end do
end subroutine AOP_accum
