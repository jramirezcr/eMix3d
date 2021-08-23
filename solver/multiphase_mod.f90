      module multiphase
      real, allocatable, dimension (:,:,:) :: phi, dmins, hevi
      real, allocatable, dimension (:,:,:) :: dxmins, dymins, dzmins
      real :: deltamin, deltamax, volphi0, volphi
      real, allocatable, dimension (:,:,:) :: rhoTwoPhase, muTwoPhase
      real, allocatable, dimension (:,:,:) :: rhoInt, muInt, sobject
      endmodule




