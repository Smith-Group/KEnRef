/*
 * KEnRefForceProvider.cpp
 *
 *      Author: amr
 */

#include <iostream>
//#include <typeinfo>
#include <cmath>
#include <memory>
#include <Eigen/Core>

#include "gromacs/mdtypes/commrec.h"
#include "mpi.h"

#include "gromacs/pbcutil/pbc.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "core/KEnRef.h"
#include "core/kabsch.h"
#include "gmxinterface/KEnRefForceProvider.h"
#include "gmxinterface/KEnRefMDModule.h"

#include <utility>
#include<unistd.h>

#define VERBOSE false
#define VALIDATE_VECTORS false

static constexpr auto singleStr = "single";

KEnRefForceProvider::KEnRefForceProvider() = default;

KEnRefForceProvider::~KEnRefForceProvider() = default;

[[maybe_unused]] KEnRefForceProvider::KEnRefForceProvider(KEnRefForceProvider &&other) noexcept = default;

[[maybe_unused]] KEnRefForceProvider::KEnRefForceProvider(const KEnRefForceProvider &other) = default;

[[maybe_unused]] void KEnRefForceProvider::setSimulationContext(gmx::SimulationContext *simulationContext) {
    this->simulationContext_ = simulationContext;
}

void KEnRefForceProvider::setGuideAtom0Indices(std::shared_ptr<std::vector<int> const> targetAtoms0Indices) {
    this->guideAtom0Indices_ = std::move(targetAtoms0Indices);
}

void KEnRefForceProvider::setGuideAtomsReferenceCoords(
        std::shared_ptr<const CoordsMatrixType<KEnRef_Real_t>> &guideAtomsReferenceCoords) {
    this->guideAtomsReferenceCoords_ = std::move(guideAtomsReferenceCoords);
    //keep another cashed version after moving its Center of Mass (COM) to the origin for faster processing
    this->guideAtomsReferenceCoordsCentered_ = std::make_shared<const CoordsMatrixType<KEnRef_Real_t>>(
            Kabsch_Umeyama<KEnRef_Real_t>::translateCenterOfMassToOrigin(*this->guideAtomsReferenceCoords_));
    //TODO do we need a guideAtomsReferenceCoordsCOM_ ?
    this->lastFrameGuideAtomsX_ZEROIndexed_ = std::make_shared<CoordsMatrixType<KEnRef_Real_t>>(
            this->guideAtomsReferenceCoords_->rows(), this->guideAtomsReferenceCoords_->cols());
}


void KEnRefForceProvider::calculateForces(const gmx::ForceProviderInput &forceProviderInput,
                                          gmx::ForceProviderOutput *forceProviderOutput) {
    auto begin = std::chrono::high_resolution_clock::now();
#if VERBOSE
    std::cout << "calculateForces() called" << std::endl;
#endif
    const auto homenr = forceProviderInput.homenr_; // total number of atoms in the system (or domain dec ?)
    GMX_ASSERT(homenr >= 0, "number of home atoms must be non-negative.");

//    const auto &x = forceProviderInput.x_;
    const auto &cr = forceProviderInput.cr_;
    const auto &step = forceProviderInput.step_;
    //    const auto& t  = forceProviderInput.t_; // Not needed (at least yet)
    const auto &force = forceProviderOutput->forceWithVirial_.force_;

    bool isMultiSimulation = this->simulationContext_->multiSimulation_ != nullptr;
    int numSimulations = isMultiSimulation ? this->simulationContext_->multiSimulation_->numSimulations_ : 1;
    int simulationIndex = isMultiSimulation ? this->simulationContext_->multiSimulation_->simulationIndex_ : 0;
    MPI_Comm mainRanksComm = isMultiSimulation
                             ? this->simulationContext_->multiSimulation_->mainRanksComm_ : MPI_COMM_NULL;
    if (!paramsInitialized) {
        std::string alt_out_path = IoUtils::strip_enclosing_quotes(IoUtils::getEnvParam("KENREF_ALT_OUT_PATH", std::string("")));
        if (!alt_out_path.empty()) {
            if (isMultiSimulation) {
                auto pos = alt_out_path.find(singleStr);
                std::string simIndexStr = std::to_string(simulationIndex + 1);
                if (pos != std::string::npos)
                    alt_out_path.replace(pos, strlen(singleStr), simIndexStr, 0, simIndexStr.length());
            }
            alt_out_path = std::filesystem::absolute(alt_out_path);
            std::cout << "Attempting to redirect output to " << alt_out_path << std::endl;
            if (std::freopen(alt_out_path.c_str(), "a", stdout)) {
                std::cout << "Redirected stdout stream" << std::endl;
            } else {
                std::cout << "FAILED to redirect output to " << alt_out_path << std::endl;
            }
        }
    }

    if (step % 10 == 0)
        std::cout
                << "--> numSimulations " << numSimulations << "\n"
                << "--> rankInDefaultCommunicator " << cr.rankInDefaultCommunicator << " " << (isMultiSimulation ? simulationIndex : -1) << "\n"
                << "--> simulationIndex " << simulationIndex << "\tstep " << step << std::endl;

    if (!paramsInitialized) {
        volatile bool holdToDebug = false;
        while (/*simulationIndex > 0 &&*/ holdToDebug) {
            sleep(1);
        }
    }


    if (!paramsInitialized) {
        GMX_ASSERT(check_box(PbcType::Unset, forceProviderInput.box_) == nullptr, "Invalid box.");
        std::cout << "Number of atoms = " << homenr << std::endl;
        std::cout << "havePPDomainDecomposition(cr): " << havePPDomainDecomposition(&cr) << std::endl;
        std::cout << "haveDDAtomOrdering(cr): " << haveDDAtomOrdering(cr) << std::endl;
        std::cout << "cr.dd->nnodes: " << cr.dd->nnodes << std::endl;

#if VERBOSE
        std::cout<< "Global to Local Atom number mapping:" << std::endl;
        for(int i = 0; i < homenr; i++){
            const int* aLocal = &i;
            if ((cr.dd == nullptr) || (aLocal = cr.dd->ga2la->findHome(i))){
                std::cout<< i<< " : " << static_cast<size_t>(*aLocal) << std::endl;
            }
        }
#endif
        fillParamsStep0(homenr, numSimulations, forceProviderInput); //TODO Optimize this function via OMP
        paramsInitialized = true;
    }

    // const auto &atomName_to_atomSub0Id_map = *this->atomName_to_atomSub0Id_map_;
    const auto &sub0Id_to_global1Id = *this->sub0Id_to_global1Id_;
    // const auto &experimentalData_table = *this->experimentalData_table_;
    // const auto &atomName_pairs = *this->atomName_pairs_;
    const auto &g0 = *this->g0_;
    CoordsMatrixType<KEnRef_Real_t> &subAtomsX = *this->subAtomsX_;
    CoordsMatrixType<KEnRef_Real_t> &allSimulationsSubAtomsX = *this->allSimulationsSubAtomsX_;
    //setting the coordinate values to ZERO (or ONE) is dangerous because it causes Invalid floating point operation

#if VERBOSE
    int iZero = 0;
    const int *piZero = cr.dd->ga2la->findHome(iZero);
    std::cout << x[*piZero][0] << " " << x[*piZero][1] << " " << x[*piZero][2] << "\t\t" << std::endl << std::endl;
    //Note that the first atom of guideAtomsX (i.e. guideAtomsX[0]) is not the same subAtomsX_[0], and even subAtomsX_[0] ((may)) later not be the first atom in the system.
#endif

    // ================= fit all models to reference ====================
    const CoordsMatrixType<KEnRef_Real_t> &guideAtomsX_ZEROIndexed =
        getGuideAtomsX(*this->guideAtom0Indices_, forceProviderInput, true);
    restoreNoJump(const_cast<CoordsMatrixType<KEnRef_Real_t> &>(guideAtomsX_ZEROIndexed),
                  *this->lastFrameGuideAtomsX_ZEROIndexed_, forceProviderInput.box_, true,
                  gmx_omp_nthreads_get(ModuleMultiThread::Default), (step % 10 == 0));

    //TODO later fix and use guideAtomsReferenceCoordsCentered_ (may need to pass the old center centered matrix + original center).
    //N.B. You must restore no jump BEFORE calling find3DAffineTransform(), applyTransform(), or applyInverseOfTransform().
    const auto &affine = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform(
            guideAtomsX_ZEROIndexed,
            *this->guideAtomsReferenceCoordsCentered_,
            false, false, false);
//    std::cout << "affine: \n" << affine.matrix() << std::endl;

    //Fill needed atoms of subAtomsX with atoms (in the original order).
    // Also, transform to Angstrom to compare/fit against Angstrom from PDB
    //N.P. We chose Angstrom, although we are working in Gromacs (Nanometer-based environment) because the energy
    // was calculated using the PDB, which is Angstrom based.
    // BTW, don't panic about scaling the forces back, because that is handled/included in the force constant.
    fillSubAtomsX(subAtomsX, sub0Id_to_global1Id, forceProviderInput, true);
    restoreNoJump(subAtomsX, *this->lastFrameSubAtomsX_, forceProviderInput.box_, true,
                  gmx_omp_nthreads_get(ModuleMultiThread::Default), (step % 10 == 0));
    (*this->lastFrameSubAtomsX_)(Eigen::indexing::all, Eigen::indexing::all) = subAtomsX(Eigen::indexing::all, Eigen::indexing::all);
    (*this->lastFrameGuideAtomsX_ZEROIndexed_)(Eigen::indexing::all, Eigen::indexing::all) = guideAtomsX_ZEROIndexed(Eigen::indexing::all, Eigen::indexing::all);

    if (haveDDAtomOrdering(cr)) {
        //TODO handle Domain Decomposition
    }

    // Copy all subAtomsXAfterFitting into its corresponding section of allSimulationsSubAtomsX (after fitting)

    //    I don't think this line is important. Only for easy printing
    //    gmx_barrier(mainRanksComm);

//    subAtomsX= Kabsch_Umeyama<KEnRef_Real_t>::translateCenterOfMassToOrigin(subAtomsX);
//    std::cout << "subAtomsX B4 transform: \n" << subAtomsX(Eigen::seqN(0, 5, 10), Eigen::indexing::all) << std::endl;
    //N.B. You must restore no jump BEFORE calling find3DAffineTransform(), applyTransform(), or applyInverseOfTransform().
    subAtomsX = Kabsch_Umeyama<KEnRef_Real_t>::applyTransform(affine, subAtomsX);
//    subAtomsXAfterFitting = subAtomsX;
//    std::cout << "subAtomsX after transform: \n" << subAtomsX(Eigen::seqN(0, 5, 10), Eigen::indexing::all) << std::endl;

    // Gather allSimulationsSubAtomsX to rank 0 (in allSimulationsSubAtomsX)
    if (isMultiSimulation) {
        MPI_Gather(subAtomsX.data(), static_cast<int>(subAtomsX.size()),
                   ((std::is_same_v<KEnRef_Real_t, float>) ? MPI_FLOAT : MPI_DOUBLE),
                   allSimulationsSubAtomsX.data(), static_cast<int>(subAtomsX.size()),
                   ((std::is_same_v<KEnRef_Real_t, float>) ? MPI_FLOAT : MPI_DOUBLE), 0,
                   mainRanksComm);
        //I don't think this line is important. Only for easy printing
        //        gmx_barrier(mainRanksComm);
    } else {
        // NOthing. They are one and the same already
    }
    // ================= end of fit all models to reference ====================

#if VERBOSE
    if(simulationIndex == 0){
        std::cout << "allSimulationsSubAtomsX shape is (" << allSimulationsSubAtomsX.rows() << ", " << allSimulationsSubAtomsX.cols() << ")" << std::endl;
//			std::cout << "allSimulationsSubAtomsX" << std::endl << allSimulationsSubAtomsX << std::endl;
        for (int r = 0; r < allSimulationsSubAtomsX.size()/10; ++r) {
            for(int i = 0; i< numSimulations; i++){
                for(int j = 0; j < 3; j++){
                    std::cout << std::setw(6) << allSimulationsSubAtomsX(subAtomsX.rows() * i + r , j) << " ";
                }
                std::cout << "\t| ";
            }
            std::cout << std::endl;
        }
    }
#endif

    KEnRef_Real_t energy = 0;
    std::optional<std::vector<CoordsMatrixType<KEnRef_Real_t> > > all_derivatives_vector_optional;

    //in the master rank
    if (!isMultiSimulation || simulationIndex == 0) {
        //collect all matrices of all replica into a vector of atom coordinates.
        std::vector<CoordsMatrixType<KEnRef_Real_t> > allSimulationsSubAtomsX_vector;
        allSimulationsSubAtomsX_vector.reserve(numSimulations);
        for (int i = 0; i < numSimulations; i++) {
            //TODO Review (and Simplify?). May need memory alignment
            allSimulationsSubAtomsX_vector.emplace_back(std::move(
                CoordsMapType<KEnRef_Real_t>(&allSimulationsSubAtomsX.data()[i * subAtomsX.size()], subAtomsX.rows(), 3)));
        }

        const auto &atomId_pairs = *this->atomId_pairs_;
        //do force calculations
        switch (selectedEnergyModel) {
            case PLATEAU_VALUES:
                // std::vector<std::vector<std::vector<int> > > &simulated_grouping_list = ;
                std::tie(energy, all_derivatives_vector_optional) =
                    KEnRef<KEnRef_Real_t>::coord_array_to_g_energy(
                        allSimulationsSubAtomsX_vector, atomId_pairs, *this->simulated_grouping_list_, g0, this->k_, this->n_,
                        true, gmx_omp_nthreads_get(ModuleMultiThread::Default));
                break;
            case SIGMA:
                std::tie(energy, all_derivatives_vector_optional) =
                    KEnRef<KEnRef_Real_t>::coord_array_to_sigma_energy(allSimulationsSubAtomsX_vector,
                        this->rates_, *this->SpecDenData_, this->proton_mhz_, this->k_, this->n_, *this->atomName_to_atomSub0Id_map_,
                        true, gmx_omp_nthreads_get(ModuleMultiThread::Default));
                break;
            case UNKNOWN:
                GMX_ASSERT(false, "ERROR: UNKNOWN energy model. please set \"ENERGY_MODEL\" in your params file");
        }
        if (step % 10 == 0)
            std::cout << "Step: " << step << " Energy: " << energy << std::endl;
#if VERBOSE
        std::cout << "energy = " << energy << ", allDerivatives_vector:" << std::endl;
        for (int i = 0; i < all_derivatives_vector_optional.value().size(); i++) {
            std::cout << "model " << i << " shape (" << all_derivatives_vector_optional.value()[i].rows() << " x " << all_derivatives_vector_optional.value()[i].cols()
                      << ")" << std::endl << all_derivatives_vector_optional.value()[i] << std::endl;
        }
#endif
    }

    //Rectify derivatives: Scatter them, then Transform them back
    CoordsMapType<KEnRef_Real_t> derivatives_map(nullptr, 0, 3);
    auto allDerivatives_buffer = this->allDerivatives_buffer_.get();
    auto derivatives_buffer = this->derivatives_buffer_.get();
    // First, prepare the buffer
    if (simulationIndex == 0) { //whether master rank or single simulation
        //I will use the slow method of copying data now, as it is less error-prone. TODO To change it, we need first to make sure the function returns the data contagiously.
        for (int i = 0; i < all_derivatives_vector_optional.value().size(); ++i) {
            const auto &matrix = all_derivatives_vector_optional.value()[i];
            std::copy_n(matrix.data(), subAtomsX.size(), &allDerivatives_buffer[i * subAtomsX.size()]);
        }
    }
    // Use it to distribute all derivatives
    if (isMultiSimulation) {
        MPI_Scatter(allDerivatives_buffer, static_cast<int>(subAtomsX.size()), ((std::is_same<KEnRef_Real_t, float>()) ? MPI_FLOAT : MPI_DOUBLE),
                    derivatives_buffer, static_cast<int>(subAtomsX.size()), ((std::is_same<KEnRef_Real_t, float>()) ? MPI_FLOAT : MPI_DOUBLE), 0, mainRanksComm);
    } else { //if a single simulation
        //Do nothing. allDerivatives_buffer and derivatives_buffer are already the same.
    }

    //once you have the derivatives, retrieve them from the buffer
    new(&derivatives_map) CoordsMapType<KEnRef_Real_t>(derivatives_buffer, subAtomsX.rows(), 3);

    // Transform them back
    //N.B. You must restore no jump BEFORE calling find3DAffineTransform(), applyTransform(), or applyInverseOfTransform().
    CoordsMatrixType<KEnRef_Real_t> derivatives_rectified = Kabsch_Umeyama<KEnRef_Real_t>::applyInverseOfTransform(affine, derivatives_map);
//    CoordsMatrixType<KEnRef_Real_t> derivatives_rectified = derivatives_map;


    KEnRef<KEnRef_Real_t>::saturate(derivatives_rectified, this->maxForceSquared_,
                                    gmx_omp_nthreads_get(ModuleMultiThread::Default));

#if VERBOSE
    std::cout << "computeVirial_ = " << std::boolalpha << forceProviderOutput->forceWithVirial_.computeVirial_ <<
            std::endl;
#endif

    if (isMultiSimulation && haveDDAtomOrdering(cr)) {
        // Note: this assumes that all ranks are hitting this line, which is not generally true.
        // I need to find the right sub-communicator. What I really want is a _scoped_ communicator...
        gmx_barrier(cr.mpi_comm_mygroup);
    }

    //Finally, add them to corresponding atoms
    for (int i = 0; i < subAtomsX.rows(); i++) {
        const int *piGlobal = new int{sub0Id_to_global1Id[i] - 1};
        const int *piLocal = cr.dd->ga2la->findHome(*piGlobal); //TODO Confirm whether it use global or local ID?

#if VERBOSE
        std::cout << sub0Id_to_global1Id[i] << "\t" << *piGlobal << "\t" << *piLocal << "\t force: " << force[*piLocal][0] << "," << force[*piLocal][1] << "," <<force[*piLocal][2] << std::endl;
#endif
        //next line assumes that the basic type of force is **real**
        // TODO optimize this line/process
        force[*piLocal] -= {
                static_cast<real>(derivatives_rectified(i, 0)),
                static_cast<real>(derivatives_rectified(i, 1)),
                static_cast<real>(derivatives_rectified(i, 2))
        };
    }

    /////////////////// print angle between vectors //////////////////////////////////////
    if (derivatives_rectified.rows() == 2) {
        const auto &vec1 = subAtomsX.row(1) - subAtomsX.row(0);
        KEnRef_Real_t norm1 = vec1.norm();
        std::cout << "Step\t" << step << "\tDist\t" << norm1 << "\tEnergy\t" << energy << std::endl;
#if VALIDATE_VECTORS
        bool found0 = false, found180 = false; KEnRef_Real_t targetLength = 7.0;
        for (int i = 0; i < derivatives_rectified.rows(); i++) {
            Eigen::RowVector3<KEnRef_Real_t> vec2 = derivatives_rectified.row(i);
            KEnRef_Real_t v1DotV2 = vec1.dot(vec2);
            KEnRef_Real_t norm2 = vec2.norm();
            std::cout << "Vec1 " << vec1 << " norm1 " << norm1 << " Vec2 " << vec2 << " norm2 " << norm2 << " dot "
                      << v1DotV2 << std::flush;
            KEnRef_Real_t cos = v1DotV2 / (norm1 * norm2);
            if (abs(cos) > 1) cos = round(cos);
            KEnRef_Real_t theta = std::acos(cos);
            const auto thetaInDegrees = theta * 180 / M_PI;
            std::cout << " cos " << cos << " theta (in degrees) " << thetaInDegrees << std::endl;

            if(thetaInDegrees < 2) {
                found0 = true;
                bool wider = norm1 > targetLength;
                bool shrinking = /*thetaInDegrees < 2 &&*/ i == 0;
                //assert either (wider and shrinking) or (smaller and expanding),
                // then negate all of that because the value is subtracted from forces
                assert(!((wider && shrinking) || (!wider && !shrinking) /*(expanding)*/));
            }else if(thetaInDegrees > 178){
                found180 = true;
                //opposite to the above condition. no need to check.
            }else{
                std::cerr << "WARNING: Theta = " << thetaInDegrees << " degrees." << std::endl;
//            assert(false);
            }
        }
//    assert(found0 && found180);
        if (!(found0 && found180)){
            std::cerr << "WARNING: angles are NOT 0 and 180 degrees." << std::endl;
        }
#endif
    }
    /////////////////// End print angle between vectors //////////////////////////////////////

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    this->calculateForces_time += elapsed.count();
    if (!(step % 10) && simulationIndex == 0)
        printf("This iteration (%ld): %.5f seconds. All walltime %.3f seconds\n", step,
               static_cast<double>(elapsed.count()) * 1e-9,
               static_cast<double>(calculateForces_time) * 1e-9);


    //    std::cout << "final force values of simulation # " << simulationIndex << std::endl;
    //    for(int i = 0; i < globalAtomIdFlags_->size() / 10; i++){
    //    	std::cout << ((*globalAtomIdFlags_).at(i) ? "*" : " ") << "\t" << force[i][0] << "\t" << force[i][1] << "\t" << force[i][2] << ((*globalAtomIdFlags_).at(i) ? "\t*" : "") << std::endl;
    //    }

    //I don't think this line is important. Only for easy printing
    //	if(isMultiSimulation) gmx_barrier(mainRanksComm);
    if (!(step % 10) && simulationIndex == 0)
        std::cout << "=================" << std::endl;
}

void KEnRefForceProvider::fillSubAtomsX(CoordsMatrixType<KEnRef_Real_t> &subAtomsX,
                                        const std::vector<int> &sub0Id_to_global1Id,
                                        const gmx::ForceProviderInput &forceProviderInput, const bool toAngstrom) {
    for (int i = 0; i < subAtomsX.rows(); i++) {
        const int *piGlobal = new int{sub0Id_to_global1Id[i] - 1};
        const int *piLocal = forceProviderInput.cr_.dd->ga2la->findHome(*piGlobal);
        GMX_ASSERT(piLocal, "ERROR: Can't find local index of atom");
        const gmx::RVec atom_x = forceProviderInput.x_[*piLocal];
#if VERBOSE
        std::cout << sub0Id_to_global1Id[i] << "\t" << *piGlobal << "\t" << *piLocal << "\t x: " << atom_x[0] << ", " << atom_x[1] << ", " << atom_x[2] << std::endl;
#endif
        if constexpr (std::is_same_v<KEnRef_Real_t, real>) {
            auto subAtomsX_buffer = subAtomsX.data();
            const auto rvec = atom_x.as_vec();
            std::copy_n(rvec, 3, &subAtomsX_buffer[i * 3]);
        } else {
            for (int j = 0; j < 3; ++j) {
                subAtomsX(i, j) = static_cast<KEnRef_Real_t>(atom_x[j]);
            }
        }
    }
    if (toAngstrom)
        subAtomsX *= 10;
}

CoordsMatrixType<KEnRef_Real_t> KEnRefForceProvider::getGuideAtomsX(const std::vector<int> &guideAtom0Indices,
                                                                    const gmx::ForceProviderInput &forceProviderInput,
                                                                    const bool toAngstrom) {
    long guideAtom0IndicesSize = static_cast<long>(guideAtom0Indices.size());
    auto guideAtomsX_ZEROIndexed = CoordsMatrixType<KEnRef_Real_t>(guideAtom0IndicesSize, 3);
    KEnRef_Real_t *guideAtomsX_ZEROIndexed_buffer = guideAtomsX_ZEROIndexed.data();
    for (auto i = 0; i < guideAtom0IndicesSize; i++) {
        const int *pi = &guideAtom0Indices[i];
        const int *piLocal = forceProviderInput.cr_.dd->ga2la->findHome(*pi);
        GMX_ASSERT(piLocal, "ERROR: Can't find local index of atom");
        const gmx::RVec atom_x = forceProviderInput.x_[*piLocal];

        auto rvec = atom_x.as_vec();
        std::copy_n(rvec, 3, &guideAtomsX_ZEROIndexed_buffer[i * 3]);
        //		// for test only
        //		for(auto j = 0; j < 3; j++){
        //			guideAtomsX_buffer[i * 3 + j] = (isMultiSimulation? 100000 * simulationIndex : 0) + 100 * i + j;
        //		}
    }
    //TODO make a unit test to validate that the value coming in rvec is equal to the value in guideAtomsX_ZEROIndexed
#if VERBOSE
    std::cout << "guideAtomsX_ZEROIndexed shape is (" << guideAtomsX_ZEROIndexed.rows() << ", " << guideAtomsX_ZEROIndexed.cols() << ")" << std::endl;
    std::cout << "guideAtomsX_ZEROIndexed" << std::endl << guideAtomsX_ZEROIndexed << std::endl;
#endif
    if (toAngstrom)
        guideAtomsX_ZEROIndexed *= 10;
    return guideAtomsX_ZEROIndexed; //RETURN BY VALUE
}

void KEnRefForceProvider::restoreNoJump(CoordsMatrixType<KEnRef_Real_t> &atoms,
                                        const CoordsMatrixType<KEnRef_Real_t> &reference,
                                        const matrix &box_, const bool toAngstrom, int numOmpThreads, const bool printStatistics) {
#define RESTORE_NO_JUMP_VERBOSE true
    matrix box; // This is real[DIM][DIM] - a 2D array on stack
    const real scale_factor = toAngstrom ? 10.0 : 1.0;
    msmul(box_, scale_factor, box);

    // Precompute all necessary box information
    const KEnRef_Real_t box_half[3] = {box[XX][XX] * 0.5, box[YY][YY] * 0.5, box[ZZ][ZZ] * 0.5};
    // Precompute non-zero dimensions to avoid branching in inner loop
    int active_dims[3];  // Stack allocation instead of vector
    int num_active = 0;


    for (int m = 2; m >= 0; --m) {
        if (box_half[m] != 0) {
            active_dims[num_active++] = m;
        }
    }

#if RESTORE_NO_JUMP_VERBOSE
    Eigen::RowVectorXi updatedLocations;
    std::atomic<bool> updated{false};  //Use atomic for thread safety
    if (printStatistics) {
        updatedLocations = Eigen::RowVectorXi::Zero(atoms.rows());
    }
#endif

    const uint num_atoms = atoms.rows();

    // Modified: Better loop scheduling for load balancing
    // TODO try also schedule(dynamic, 64) or schedule(dynamic, 16) or even smaller
#if RESTORE_NO_JUMP_VERBOSE
    #pragma omp parallel for num_threads(numOmpThreads) schedule(guided)\
        reduction(||:updated) shared(atoms, box, box_half, reference, updatedLocations, printStatistics, active_dims)
#else
    #pragma omp parallel for num_threads(numOmpThreads) schedule(guided)\
        shared(atoms, box, box_half, reference, printStatistics, active_dims)
#endif
    for (int i = 0; i < num_atoms; ++i) {
#if RESTORE_NO_JUMP_VERBOSE
        bool local_updated = false;
#endif
        // Get row references once per atom
        auto atom = atoms.row(i);
        const auto &refAtom = reference.row(i);

        //Iterate only over active dimensions
        for (int dim_idx = 0; dim_idx < num_active; ++dim_idx) {
            const int m = active_dims[dim_idx];
            const KEnRef_Real_t diff = atom[m] - refAtom[m];
            const KEnRef_Real_t half_len = box_half[m];

            if (diff <= -half_len) {
                KEnRef_Real_t correction_count = std::floor((-diff - half_len) / half_len * 2.0) + 1.0;
                // Manual unrolling for DIM=3
                if (m >= 0) atom[0] += correction_count * box[m][0];
                if (m >= 1) atom[1] += correction_count * box[m][1];
                if (m >= 2) atom[2] += correction_count * box[m][2];
#if RESTORE_NO_JUMP_VERBOSE
                if (printStatistics) {
                    updatedLocations(i) = 1;
                    local_updated = true;
                }
#endif
            }
            else if (diff > half_len) {
                KEnRef_Real_t correction_count = std::floor((diff - half_len) / half_len * 2.0) + 1.0;
                if (m >= 0) atom[0] -= correction_count * box[m][0];
                if (m >= 1) atom[1] -= correction_count * box[m][1];
                if (m >= 2) atom[2] -= correction_count * box[m][2];
#if RESTORE_NO_JUMP_VERBOSE
                if (printStatistics) {
                    updatedLocations(i) = 1;
                    local_updated = true;
                }
#endif
            }
        }
#if RESTORE_NO_JUMP_VERBOSE
        if (local_updated) {
            updated.store(true, std::memory_order_relaxed);
        }
#endif
    }

#if RESTORE_NO_JUMP_VERBOSE
    if (printStatistics && updated.load(std::memory_order_relaxed)) {
        // More efficient output for large systems
        const int updated_count = updatedLocations.sum();
        std::cout << "INFO: Restored NoJump in " << updated_count
                  << " atoms (out of " << num_atoms << ")" << std::endl;
        if (updated_count < 50) {  // Only show details for small numbers
            std::cout << "Updated atom indices:\n" << updatedLocations << std::endl;
        }
    }
#endif
    // delete[] box; //no need anymore
#undef RESTORE_NO_JUMP_VERBOSE
}

void KEnRefForceProvider::fillParamsStep0(const size_t homenr, int numSimulations, const gmx::ForceProviderInput &forceProviderInput) {
    auto begin = std::chrono::high_resolution_clock::now();
    bool isMultiSimulation = this->simulationContext_->multiSimulation_ != nullptr;
    std::cout << "Energy model: ";
    switch (selectedEnergyModel) {
        case SIGMA:
            std::cout << "SIGMA" << std::endl;
            break;
        case PLATEAU_VALUES:
            std::cout << "PLATEAU_VALUES" << std::endl;
            break;
        default:
            std::cout << "UNKNOWN" << std::endl;
            break;
    }
    this->atomName_to_atomGlobalId_map_ = std::make_shared<std::map<std::string, int> >(
        IoUtils::getAtomMappingFromPdb<std::string, int>(KEnRefMDModule::ATOMNAME_MAPPING_FILENAME,IoUtils::fill_atomId_to_index_Map));
    GMX_ASSERT(!atomName_to_atomGlobalId_map_->empty(), "No atom mapping found");
    auto &atomName_to_atomGlobalId_map = *this->atomName_to_atomGlobalId_map_;

    this->maxForceSquared_ = std::pow(IoUtils::getEnvParam("KENREF_MAXFORCE", std::sqrt(this->maxForceSquared_)), 2.f);
    this->k_ = IoUtils::getEnvParam("KENREF_K", this->k_);
    this->n_ = IoUtils::getEnvParam("KENREF_N", this->n_);
    int maxAtomPairsToRead = IoUtils::getEnvParam("KENREF_MAX_ATOMPAIRS_TO_READ", -1); //TODO remove

    std::cout << "KEnRef_Real_t type is: " << typeid(KEnRef_Real_t).name() << '\n';

    if (selectedEnergyModel == SIGMA) {
        this->proton_mhz_ = IoUtils::getEnvParam("KENREF_PROTON_MHZ", this->proton_mhz_);

        const bool handleNames = IoUtils::should_handleNames(atomName_to_atomGlobalId_map);

        int maxId0 = -1;
        for (const auto &[atomName, id1]: atomName_to_atomGlobalId_map) {
            std::string tempName = atomName;
            if (const int id0 = id1 - 1; id0 > maxId0) {
                maxId0 = id0;
            }
        }

        const std::vector<std::string> spec_den_data_prefixes  {"1-1", "1-2", "1-3", "2-2", "2-3", "3-3"};
        std::vector<SpecDenData<KEnRef_Real_t>> spec_den_data_vector;
        spec_den_data_vector.reserve(spec_den_data_prefixes.size());
        //TODO unify `this->atomName_pairs_` among model cases
        this->atomName_pairs_ = new std::vector<std::tuple<std::string, std::string> >();
        for (int i = 0; i < spec_den_data_prefixes.size(); ++i) {
            //AtomPairs
            std::string atomPairAndSigmaFileName = KEnRefMDModule::EXPERIMENTAL_DATA_FOLDER + spec_den_data_prefixes[i]+"_atom_pairs.csv";
            std::cout << atomPairAndSigmaFileName << std::endl;
            const Table& atomPairAndSigmaTable = IoUtils::readTable(atomPairAndSigmaFileName,true,false, "\\s*,\\s*", -1, true);
            std::vector<std::tuple<std::string, std::string>> atomPairs(atomPairAndSigmaTable.rowCount());
            for (int j = 0; j < atomPairAndSigmaTable.rowCount(); ++j) {
                auto atom1 = IoUtils::normalizeName(atomPairAndSigmaTable.at(j, 0), handleNames);
                auto atom2 = IoUtils::normalizeName(atomPairAndSigmaTable.at(j, 1), handleNames);
                atomPairs.at(j) = std::move(std::tuple<std::string, std::string>{ atom1, atom2 });
                this->atomName_pairs_->emplace_back(atom1, atom2);
            }
            // sigma
            std::vector<KEnRef_Real_t> sigmasVec{};
            for (int row = 0; row < atomPairAndSigmaTable.rowCount(); ++row) {
                if (! atomPairAndSigmaTable.isRowComplete(row))
                    break;
                const auto &valueStr = atomPairAndSigmaTable.at(row, 2);
                std::istringstream iss(valueStr);
                KEnRef_Real_t value;
                iss >> value;
                sigmasVec.emplace_back(value);
            }
            std::optional<NamedVector<KEnRef_Real_t>> sigma = NamedVector<KEnRef_Real_t>(sigmasVec.size());
            for (int j = 0; j < sigma->rows(); ++j) {
                sigma.value()(j, 0) = sigmasVec[j];
            }
            //multiple_grouping
            std::string multiple_grouping_fileName = KEnRefMDModule::EXPERIMENTAL_DATA_FOLDER + spec_den_data_prefixes[i]+"_groupings.csv";
            const auto &grouping_matrix = NamedMatrix<int>(IoUtils::readTable(multiple_grouping_fileName, false, false).toNamedMatrix<int>().array() - 1);
            const auto &multiple_grouping = IoUtils::grouping_mat_to_subset_idx(grouping_matrix);
            //a_coef
            std::string aCoefFileName = KEnRefMDModule::EXPERIMENTAL_DATA_FOLDER + spec_den_data_prefixes[i]+"_a_coef.csv";
            const auto &a_coef = IoUtils::readTable(aCoefFileName, true,false, "\\s*,\\s*", -1, false).toNamedMatrix<KEnRef_Real_t>();
            //lambda_coef
            std::string lambdaCoefFileName = KEnRefMDModule::EXPERIMENTAL_DATA_FOLDER + spec_den_data_prefixes[i]+"_lambda_coef.csv";
            const auto &lambda_coef = IoUtils::readTable(lambdaCoefFileName, true,true, "\\s*,\\s*", -1, false).toNamedMatrix<KEnRef_Real_t>();

            // spec_den_data_vector.emplace_back(atomPairs, sigma, multiple_grouping, a_coef, lambda_coef);
            spec_den_data_vector.emplace_back(std::move(SpecDenData<KEnRef_Real_t>{atomPairs, sigma, multiple_grouping, a_coef, lambda_coef}));
        }
        this->SpecDenData_ = std::make_shared<std::vector<SpecDenData<KEnRef_Real_t>>>(spec_den_data_vector);
    }else /*if (selectedEnergyModel == PLATEAU_VALUES)*/ {

        this->experimentalData_table_ = std::make_shared<Table>(
            IoUtils::readTable(KEnRefMDModule::EXPERIMENTAL_DATA_FILENAME, true, false,
                "\\s*,\\s*", maxAtomPairsToRead));
        GMX_ASSERT(experimentalData_table_ && !(maxAtomPairsToRead && experimentalData_table_->rowCount() > 0), "No simulated data found");
        const Table &table = *experimentalData_table_;
#if VERBOSE
        const auto& [table_header, table_data] = *experimentalData_table_;
        IoUtils::printVector(table_header);
        for(const auto& record: table_data){
            IoUtils::printVector(record);
        }
#endif
        //First confirm whether we should handle the atom names
        bool handleUnpreparedAtomNames = false;
        for (int row = 0; row <table.rowCount(); row++) {
            //        std::cout << "checking [" << record[1] <<"] and [" << record[2] << "]: ";
            if (IoUtils::isNotPrepared(table(row,"atom1")) || IoUtils::isNotPrepared(table(row, "atom2"))) {
                //            std::cout << "TRUE, Exiting" << std::endl;
                std::cerr << "WARNING: It seems that your data is from an unprepared file. We will try to handle it, but we can not guarantee the results." << std::endl;
                handleUnpreparedAtomNames = true;
                break;
            }
            //        std::cout << "FALSE" << std::endl;
        }
        //TODO if handleUnpreparedAtomNames is false, we need to simplify the code
        this->atomName_pairs_ = new std::vector<std::tuple<std::string, std::string> >();
        for (int row = 0; row < table.rowCount(); row++) {
            std::string atom1 = IoUtils::normalizeName(table(row, "atom1"), handleUnpreparedAtomNames);
            std::string atom2 = IoUtils::normalizeName(table(row, "atom2"), handleUnpreparedAtomNames);
            this->atomName_pairs_->emplace_back(atom1, atom2);
        }
#if VERBOSE
        for(auto [atom1, atom2]: *atomName_pairs_){
            std::cout << "[" << atom1 << "], [" << atom2 << "]" << std::endl;
        }
#endif
        this->g0_ = new Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, Eigen::Dynamic>(table.rowCount(), 2);
        auto &g0 = *g0_;
        for (int i = 0; i < g0.rows(); ++i) {
            std::istringstream temp1(table(i,"g1")), temp2(table(i, "g2"));
            temp1 >> g0(i, 0);
            temp2 >> g0(i, 1);
        }
#if VERBOSE
        std::cout << *g0_ << std::endl;
#endif

    }

#if VERBOSE
    for (const auto& [name, globalId] : atomName_to_atomGlobalId_map){
        std::cout << "[" << name << "]\t:" << globalId << std::endl;
    }
    // for(auto& entry: atomName_to_atomGlobalId_map_){
    //     auto[name, globalId] = entry;
    //     std::cout << "[" << name << "]\t:" << globalId << std::endl;
    // }
#endif
    int maxAtomIdOfInterest = -1; // If you want to use size_t, then you can NOT use -1 as an initial value
    this->globalAtomIdFlags_ = std::make_shared<std::vector<bool> >(homenr, false);
    auto &globalAtomIdFlags = *this->globalAtomIdFlags_; //ONE based
    //scan the atom pairs to do:
    //1) find the highest globalAtomId number of interest
    //2) fill in the subAtomsFilter
    //3) do a quick sanity scan on the availability of all atomname atomID maping
    int tempI;
    for (const auto &[a1, a2]: *this->atomName_pairs_) {
#if VERBOSE
        std::cout << "[" << a1 << "]\t" << atomName_to_atomGlobalId_map.at(a1) << "\t";
        std::cout << "[" << a2 << "]\t" << atomName_to_atomGlobalId_map.at(a2) << std::endl;
#endif
        //In the next lines I use .at() instead of [] deliberately; to throw an exception if unexpected name found
        if ((tempI = atomName_to_atomGlobalId_map.at(a1)) > maxAtomIdOfInterest) maxAtomIdOfInterest = tempI;
        globalAtomIdFlags[tempI] = true; // Here I use [] instead of at(), because I am sure about the boundaries.
        if ((tempI = atomName_to_atomGlobalId_map.at(a2)) > maxAtomIdOfInterest) maxAtomIdOfInterest = tempI;
        globalAtomIdFlags[tempI] = true; // Here I use [] instead of at(), because I am sure about the boundaries.
    }
#if VERBOSE
    IoUtils::printVector(globalAtomIdFlags);
#endif
    globalAtomIdFlags.resize(maxAtomIdOfInterest + 1);

    auto global1Id_to_sub0Id = std::vector<int>(globalAtomIdFlags.size(), -1);
    auto sub0Id_to_global1Id = std::vector<int>(globalAtomIdFlags.size(), -1);

    int localId = 0;
    for (int i = 0; i < globalAtomIdFlags.size(); i++) {
        if (globalAtomIdFlags[i]) {
            global1Id_to_sub0Id[i] = localId;
            sub0Id_to_global1Id[localId] = i;
            localId++;
        }
    }
    sub0Id_to_global1Id.resize(localId);

    this->global1Id_to_sub0Id_ = std::make_shared<std::vector<int> >(std::move(global1Id_to_sub0Id));
    this->sub0Id_to_global1Id_ = std::make_shared<std::vector<int> >(std::move(sub0Id_to_global1Id));

    this->atomName_to_atomSub0Id_map_ = std::make_shared<std::map<std::string, int> >();
    auto &atomName_to_atomSub0Id_map = *this->atomName_to_atomSub0Id_map_;
    for (const auto &[name, globalId]: atomName_to_atomGlobalId_map){
        if (globalId > maxAtomIdOfInterest)
            continue;
        atomName_to_atomSub0Id_map[name] = global1Id_to_sub0Id_->at(globalId);
    }
    this->atomId_pairs_ = std::make_shared<std::vector<std::tuple<int, int> > >(KEnRef<KEnRef_Real_t>::atomNamePairs_2_atomIdPairs(*atomName_pairs_, atomName_to_atomSub0Id_map));
    // Cache the atomId_pairs subgroups in every specDenData object here
    if (selectedEnergyModel == SIGMA) {
        for (int i = 0; i < this->SpecDenData_->size(); ++i) {
            auto & currentSpecDenData = this->SpecDenData_->at(i);
            const auto &atomPairs = currentSpecDenData.get_atom_pairs();
            currentSpecDenData.set_atomIdPairs_to_sub0Atom_id_pairs_cache({
                KEnRef<KEnRef_Real_t>::atomNamePairs_2_atomIdPairs(atomPairs, atomName_to_atomSub0Id_map)
            });
        }
    }

#if VERBOSE
    for(const auto& [name, subId]: atomName_to_atomSub0Id_map){
        std::cout << "[" << name << "]\t:" << subId << std::endl;
    }
#endif

    if (selectedEnergyModel == PLATEAU_VALUES) {
        switch (numSimulations) {
            case 1:
                this->simulated_grouping_list_ = std::make_shared<std::vector<std::vector<std::vector<int>>>>(
                    std::vector<std::vector<std::vector<int> > >{{{0}}, {{0}}});
                break;
            case 2:
                this->simulated_grouping_list_ = std::make_shared<std::vector<std::vector<std::vector<int>>>>(
                    std::vector<std::vector<std::vector<int>>>{{{0, 1}}, {{0}, {1}}});
                break;
            case 3:
                this->simulated_grouping_list_ = std::make_shared<std::vector<std::vector<std::vector<int>>>>(
                    std::vector<std::vector<std::vector<int>>>{{{0, 1, 2}}, {{0}, {1}, {2}}});
                break;
            default:
                GMX_ASSERT(numSimulations <= 3, "I don't know how to handle more than 3 simulations yet");
        }
    }
    this->subAtomsX_ = std::make_shared<CoordsMatrixType<KEnRef_Real_t> >(this->sub0Id_to_global1Id_->size(), 3);//contains needed atoms only
    this->lastFrameSubAtomsX_ = std::make_shared<CoordsMatrixType<KEnRef_Real_t>>(this->subAtomsX_->rows(), this->subAtomsX_->cols());
#if VERBOSE
    auto subAtomsX = *this->subAtomsX_; std::cout << "subAtomsX_ shape is (" << subAtomsX.rows() << ", " <<
        subAtomsX.cols() << ")" << std::endl;
#endif
    this->allSimulationsSubAtomsX_ =
            isMultiSimulation
                ? std::make_shared<CoordsMatrixType<KEnRef_Real_t> >(numSimulations * this->subAtomsX_->rows(), 3)
                : this->subAtomsX_;
#if VERBOSE
    auto allSimulationsSubAtomsX = *this->allSimulationsSubAtomsX_; std::cout << "allSimulationsSubAtomsX_ shape is ("
        << allSimulationsSubAtomsX.rows() << ", " << allSimulationsSubAtomsX.cols() << ")" << std::endl;
#endif

    this->allDerivatives_buffer_ = std::shared_ptr<KEnRef_Real_t[]>(
        new KEnRef_Real_t[this->subAtomsX_->size() * numSimulations]);
    this->derivatives_buffer_ =
        isMultiSimulation
            ? std::shared_ptr<KEnRef_Real_t[]>(new KEnRef_Real_t[this->subAtomsX_->size()])
            : this->allDerivatives_buffer_;

    fillSubAtomsX(*this->lastFrameSubAtomsX_, *this->sub0Id_to_global1Id_, forceProviderInput, true);
    (*this->lastFrameGuideAtomsX_ZEROIndexed_)(Eigen::indexing::all, Eigen::indexing::all) = getGuideAtomsX(*this->guideAtom0Indices_, forceProviderInput, true);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    this->calculateForces_time -= elapsed.count(); //Exclude THIS method from wall time calculation
}

#undef VERBOSE
