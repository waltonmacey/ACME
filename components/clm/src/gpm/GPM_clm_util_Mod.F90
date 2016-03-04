module GPM_clm_util_mod

private

integer, public, save      :: GPM_MW_nchannel = 9             ! number of channels
integer, public, parameter :: GPM_npol = 2                     ! number of polarizations. (H and V)
real,    public, parameter :: GPM_DefaultEmissivity_MW = 0.95  ! default emissivity at microwave frequencies
     

end module GPM_clm_util_mod
