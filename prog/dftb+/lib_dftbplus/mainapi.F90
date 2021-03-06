!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> main module for the DFTB+ API
module dftbp_mainapi
  use dftbp_accuracy, only : dp, mc
  use dftbp_assert
  use dftbp_densedescr, only : TDenseDescr
  use dftbp_environment, only : TEnvironment
  use dftbp_initprogram, only : initProgramVariables, destructProgramVariables,                   &
       & energy, derivs, TRefExtPot, refExtPot, orb, sccCalc, chrgForces, qDepExtPot,             &
       & nAtom, nSpin, nEl0, nEl, speciesName, speciesMass, coord0, latVec, species0, mass,       &
       & tCoordsChanged, tLatticeChanged, tExtField, tExtChrg, tForces, tSccCalc, tDFTBU,         &
       & tMulliken, tSpin, tReadChrg, tMixBlockCharges, tRangeSep, t2Component, tRealHS,          &
       & q0, qInput, qOutput, qInpRed, qOutRed, referenceN0, qDiffRed, nrChrg, nrSpinPol,         &
       & setEquivalencyRelations, iEqOrbitals, nIneqOrb, nMixElements, onSiteElements, denseDesc, &
       & parallelKS, HSqrCplx, SSqrCplx, eigVecsCplx, HSqrReal, SSqrReal, eigVecsReal, getDenseDescCommon, &
       & initializeReferenceCharges, setNElectrons, initializeCharges, qBlockIn, qBlockOut,       &
       & qiBlockIn, qiBlockOut, iEqBlockDFTBU, iEqBlockOnSite, iEqBlockDFTBULS, iEqBlockOnSiteLS, &
       & tStress, totalStress
#:if WITH_SCALAPACK
  use dftbp_initprogram,  only : getDenseDescBlacs
#:endif
  use dftbp_main,         only : processGeometry
  use dftbp_message,      only : error
  use dftbp_orbitals,     only : TOrbitals
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy
#:if WITH_SCALAPACK
  use dftbp_scalapackfx,  only : scalafx_getlocalshape
#:endif
  use dftbp_wrappedintr
  implicit none
  private

  public :: initProgramVariables, destructProgramVariables
  public :: setGeometry, setQDepExtPotProxy, setExternalPotential, setExternalCharges
  public :: getEnergy, getGradients, getExtChargeGradients, getGrossCharges, getStressTensor
  public :: nrOfAtoms
  public :: updateDataDependentOnSpeciesOrdering, checkSpeciesNames

contains

  !> Sets up the atomic geometry
  subroutine setGeometry(coords, latVecs)

    !> atom coordinates
    real(dp), intent(in) :: coords(:,:)

    !> lattice vectors, if periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    @:ASSERT(size(coords,1) == 3)
    coord0(:,:) = coords
    tCoordsChanged = .true.
    if (present(latVecs)) then
      latVec(:,:) = latVecs
      tLatticeChanged = .true.
    else
      tLatticeChanged = .false.
    end if

  end subroutine setGeometry


  !> Returns the free energy of the system at finite temperature
  subroutine getEnergy(env, merminEnergy)

    !> Instance
    type(TEnvironment), intent(inout) :: env

    !> Resulting energy
    real(dp), intent(out) :: merminEnergy

    call recalcGeometry(env)
    merminEnergy = energy%EMermin

  end subroutine getEnergy


  !> get forces on atoms
  subroutine getGradients(env, gradients)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> resulting gradients wrt atom positions
    real(dp), intent(out) :: gradients(:,:)

    if (.not. tForces) then
      call error("Forces not available, you must initialise your calculator&
          & with forces enabled.")
    end if

    @:ASSERT(size(gradients,1) == 3)

    call recalcGeometry(env)
    gradients(:,:) = derivs

  end subroutine getGradients


  !> get stress tensor for unit cell
  subroutine getStressTensor(env, stress)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> resulting gradients wrt atom positions
    real(dp), intent(out) :: stress(:,:)

    if (.not. tStress) then
      call error("Stress tensor not available, you must initialise your calculator with&
          & this property enabled.")
    end if

    call recalcGeometry(env)
    stress(:,:) = totalStress

  end subroutine getStressTensor


  !> get the gross (Mulliken projected) charges for atoms wrt neutral atoms
  subroutine getGrossCharges(atomCharges)

    !> resulting charges
    real(dp), intent(out) :: atomCharges(:)

    atomCharges(:) = sum(q0(:, :, 1) - qOutput(:, :, 1), dim=1)

  end subroutine getGrossCharges

  !> Sets up an external population independent electrostatic potential.
  !>
  !> Sign convention: charge of electron is considered to be positive.
  !>
  subroutine setExternalPotential(atomPot, shellPot, potGrad)

    !> Atomic external potential
    real(dp), intent(in), optional :: atomPot(:)

    !> Shell resolved electrostatic potential
    real(dp), intent(in), optional :: shellPot(:,:)

    !> Gradient of the electrostatic potential
    real(dp), intent(in), optional :: potGrad(:,:)

    ! Using explicit allocation instead of F2003 automatic ones in order to stop eventual
    ! shape mismatches already at this point rather than later deep in the main code
    if (present(atomPot)) then
      if (.not. allocated(refExtPot%atomPot)) then
        allocate(refExtPot%atomPot(nAtom, nSpin))
      end if
      @:ASSERT(all(shape(atomPot) == [nAtom]))
      refExtPot%atomPot(:,1) = atomPot
    end if
    if (present(shellPot)) then
      if (.not. allocated(refExtPot%shellPot)) then
        allocate(refExtPot%shellPot(orb%mShell, nAtom, nSpin))
      end if
      @:ASSERT(all(shape(shellPot) == [orb%mShell, nAtom]))
      refExtPot%shellPot(:,:,1) = shellPot
    end if
    if (present(potGrad)) then
      if (.not. allocated(refExtPot%potGrad)) then
        allocate(refExtPot%potGrad(3, nAtom))
      end if
      @:ASSERT(all(shape(potGrad) == [3, nAtom]))
      refExtPot%potGrad(:,:) = potGrad
    end if
    tExtField = .true.

  end subroutine setExternalPotential


  !> Sets up a generator for external population dependant potentials
  subroutine setQDepExtPotProxy(extPotProxy)

    !> Generator for the external population dependant potential
    type(TQDepExtPotProxy), intent(in) :: extPotProxy

    qDepExtPot = extPotProxy

  end subroutine setQDepExtPotProxy


  !> Sets up external point charges
  subroutine setExternalCharges(chargeCoords, chargeQs, blurWidths)

    !> Coordiante of the external charges
    real(dp), intent(in) :: chargeCoords(:,:)

    !> Charges of the external point charges (sign convention: electron is negative)
    real(dp), intent(in) :: chargeQs(:)

    !> Widths of the Gaussian for each charge used for blurring (0.0 = no blurring)
    real(dp), intent(in), optional :: blurWidths(:)

    tExtChrg = .true.
    if (tForces) then
      if (.not. allocated(chrgForces)) then
        allocate(chrgForces(3, size(chargeQs)))
      end if
    end if
    call sccCalc%setExternalCharges(chargeCoords, chargeQs)

  end subroutine setExternalCharges


  !> Returns the gradient acting on the external point charges
  subroutine getExtChargeGradients(chargeGradients)

    !> Gradients
    real(dp), intent(out) :: chargeGradients(:,:)

    @:ASSERT(tForces .and. allocated(chrgForces))

    chargeGradients(:,:) = chrgForces

  end subroutine getExtChargeGradients


  !> Obtains number of atoms in the system
  function nrOfAtoms()
    integer :: nrOfAtoms

    nrOfAtoms = nAtom

  end function nrOfAtoms


  !> Check that the order of speciesName remains constant
  !> Keeping speciesNames constant avoids the need to reset:
  !> atomEigVal, referenceN0, speciesMass and SK parameters
  !
  !  Even if nAtom is not conserved, it should be possible to know the
  !  total number of species types, nTypes, in a simulation and hence
  !  always keep speciesName constant 
  function checkSpeciesNames(inputSpeciesName) result(tSpeciesNameChanged)
   
    !> Labels of atomic species from external program
    character(len=*), intent(in) :: inputSpeciesName(:)
    !> Has speciesName changed?
    logical                   :: tSpeciesNameChanged
    
    tSpeciesNameChanged = any(speciesName /= inputSpeciesName)

  end function checkSpeciesNames


  !> When order of atoms changes, update arrays containing atom type indices,
  !> and all subsequent dependencies.
  !  Updated data returned via module use statements
  subroutine updateDataDependentOnSpeciesOrdering(env, inputSpecies)
    
    !> dftb+ environment 
    type(TEnvironment),   intent(in) :: env
    !> types of the atoms (nAllAtom) 
    integer,              intent(in) :: inputSpecies(:)
    !> Dummy arguments. Won't be used if not allocated
    real(dp), allocatable :: initialCharges(:), initialSpins(:,:)
    type(TWrappedInt1), allocatable :: customOccAtoms(:)
    real(dp), allocatable :: customOccFillings(:,:)

    if(size(inputSpecies) /= nAtom)then
       call error("Number of atoms must be keep constant in simulation")
    endif
    
    species0 = inputSpecies
    mass =  updateAtomicMasses(species0)
    orb%nOrbAtom = updateAtomicOrbitals(species0)
    
    !Used in partial charge initialisation
    call setEquivalencyRelations(species0, sccCalc, orb, onSiteElements, iEqOrbitals, &
         & iEqBlockDFTBU, iEqBlockOnSite, iEqBlockDFTBULS, iEqBlockOnSiteLS, nIneqOrb, nMixElements)
#:if WITH_SCALAPACK
    call updateBLACSDecomposition(env, denseDesc)
    call reallocateHSArrays(env, denseDesc, HSqrCplx, SSqrCplx, eigVecsCplx, HSqrReal, &
         & SSqrReal, eigVecsReal)
#:endif
    !If atomic order changes, partial charges need to be initialised,                                                   
    !else wrong charge will be associated with each atom                                                                
    call initializeReferenceCharges(species0, referenceN0, orb, customOccAtoms, &
         & customOccFillings, q0)
    call setNElectrons(q0, nrChrg, nrSpinPol, nEl, nEl0)
    call initializeCharges(species0, speciesName, orb, nEl, iEqOrbitals, nIneqOrb, &
         & nMixElements, initialSpins, initialCharges, nrChrg, q0, qInput, qOutput, &
         & qInpRed, qOutRed, qDiffRed, qBlockIn, qBlockOut, qiBlockIn, qiBlockOut)
    
  end subroutine updateDataDependentOnSpeciesOrdering


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update order of nr. atomic orbitals for each atom, orb%nOrbAtom
  function updateAtomicOrbitals(species0) result(nOrbAtomReordered)
    !> Type of the atoms (nAtom)
    integer,  intent(in)  :: species0(:)
    !> Nr. of orbitals for each atom (nAtom)  
    integer,  allocatable :: nOrbAtomReordered(:)
    
    @:ASSERT(nAtom == size(species0))
    allocate(nOrbAtomReordered(nAtom))
    nOrbAtomReordered(:) = orb%nOrbSpecies(species0(:))
  end function updateAtomicOrbitals


  !> Update atomic masses
  function updateAtomicMasses(species0) result(massReordered)
    !> Type of the atoms (nAtom)
    integer,  intent(in)  :: species0(:)
    !> List of atomic masses (nAtom)
    real(dp), allocatable :: massReordered(:)

    @:ASSERT(nAtom == size(speciesMass))
    allocate(massReordered(nAtom))
    massReordered = speciesMass(species0)
  end function updateAtomicMasses

    
  !> Update dense matrix descriptor for H and S in BLACS decomposition
  subroutine updateBLACSDecomposition(env, denseDesc)

    !> Environment settings
    type(TEnvironment), intent(in)    :: env
    !> Dense matrix descriptor for H and S
    type(TDenseDescr),  intent(inout) :: denseDesc

  #:if WITH_SCALAPACK
    ! Specificaly, denseDesc uses orb%nOrbAtom
    call getDenseDescCommon(orb, nAtom, t2Component, denseDesc)
    call getDenseDescBlacs(env, env%blacs%rowBlockSize, env%blacs%columnBlockSize, denseDesc)
  #:endif

  end subroutine updateBLACSDecomposition

  
  !> Reassign Hamiltonian, overlap and eigenvector arrays
  !
  ! If nAtom is constant and one is running without BLACS, this should not be required
  ! hence preprocessed out 
  ! May require extending if ((nAtom not constant) and (not BLACS))
  subroutine reallocateHSArrays(env, denseDesc, HSqrCplx, SSqrCplx, eigVecsCplx, &
       & HSqrReal, SSqrReal ,eigVecsReal)

    !> Environment instance
    type(TEnvironment), intent(in) :: env
    !> Dense matrix descriptor for H and S
    type(TDenseDescr),  intent(in) :: denseDesc
    !> Square dense hamiltonian storage for cases with k-points
    complex(dp), allocatable, intent(inout) :: HSqrCplx(:,:)
    !> Square dense overlap storage for cases with k-points
    complex(dp), allocatable, intent(inout) :: SSqrCplx(:,:)
    !> Complex eigenvectors
    complex(dp), allocatable, intent(inout) :: eigvecsCplx(:,:,:)
    !> Square dense hamiltonian storage
    real(dp),    allocatable, intent(inout) :: HSqrReal(:,:)
    !> Square dense overlap storage
    real(dp),    allocatable, intent(inout) :: SSqrReal(:,:)
    !> Real eigenvectors
    real(dp),    allocatable, intent(inout) :: eigvecsReal(:,:,:)
    ! Local variables 
    integer :: nLocalRows, nLocalCols, nLocalKS

#:if WITH_SCALAPACK 
    !Retrieved from index array for spin and k-point index
    nLocalKS = size(parallelKS%localKS, dim=2)
  
    !Get nLocalRows and nLocalCols 
    call scalafx_getlocalshape(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, &
                               nLocalRows, nLocalCols)
    
    !Complex
    if (t2Component .or. .not. tRealHS) then
       if( (size(HSqrCplx,1) /= nLocalRows) .or. (size(HSqrCplx,2) /= nLocalCols) )then
          deallocate(HSqrCplx)
          deallocate(SSqrCplx)
          deallocate(eigVecsCplx)
          allocate(HSqrCplx(nLocalRows, nLocalCols))
          allocate(SSqrCplx(nLocalRows, nLocalCols))
          allocate(eigVecsCplx(nLocalRows, nLocalCols, nLocalKS))
       endif
    !Real   
    else
       if( (size(HSqrReal,1) /= nLocalRows) .or. (size(HSqrReal,2) /= nLocalCols) )then
          deallocate(HSqrReal)
          deallocate(SSqrReal)
          deallocate(eigVecsReal)
          allocate(HSqrReal(nLocalRows, nLocalCols))
          allocate(SSqrReal(nLocalRows, nLocalCols))
          allocate(eigVecsReal(nLocalRows, nLocalCols, nLocalKS))
       endif

    end if

#:endif
  end subroutine reallocateHSArrays


  !> re-evaluate the energy/forces if the geometry changes
  subroutine recalcGeometry(env)

    !> instance
    type(TEnvironment), intent(inout) :: env

    logical :: tStopDriver, tStopScc, tExitGeoOpt

    if (tLatticeChanged .or. tCoordsChanged) then
      call processGeometry(env, 1, 1, .false., tStopDriver, tStopScc, tExitGeoOpt)
      tLatticeChanged = .false.
      tCoordsChanged = .false.
    end if

  end subroutine recalcGeometry


end module dftbp_mainapi
