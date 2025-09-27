/*
 * KEnRef.h
 *
 *  Created on: May 8, 2023
 *      Author: amr
 */

#ifndef KENREF_H_
#define KENREF_H_

#include <memory>
#include <utility>
#include <vector>
#include <map>
#include <tuple>
#include <string>

#include <Eigen/Dense>

#include "Table.h"
// #include "../config/KEnRefConfig.h"

#ifdef DOUBLE
typedef double KEnRef_Real_t;
#else
typedef float KEnRef_Real_t;
#endif
template<typename KEnRef_Real>
using CoordsMatrixType = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 3, Eigen::RowMajor>;
template<typename KEnRef_Real>
using CoordsMapType = Eigen::Map<CoordsMatrixType<KEnRef_Real> >;
template<typename KEnRef_Real>
using CoordsMapTypeConst = Eigen::Map<const CoordsMatrixType<KEnRef_Real>>; // a read-only map

template<typename KEnRef_Real>
class SpecDenData {
    std::vector<std::tuple<std::string, std::string> > atom_pairs{};
    std::optional<NamedVector<KEnRef_Real>> sigma = std::nullopt;
    std::vector<std::vector<std::vector<int> > > multiple_grouping{};
    NamedMatrix<KEnRef_Real> a_coef;
    NamedMatrix<KEnRef_Real> lambda_coef;
    std::map<std::tuple<std::string, std::string>, size_t> atomNamePairs_to_atomPairIndex;
    std::optional<std::vector<std::tuple<int, int> > > atomIdPairs_to_sub0Atom_id_pairs_cache_ = std::nullopt;
    KEnRef_Real getSigma(size_t atomPairId) const;

public:
    SpecDenData(const std::vector<std::tuple<std::string, std::string> > &atom_pairs,
                const std::optional<NamedVector<KEnRef_Real> > &sigma,
                const std::vector<std::vector<std::vector<int> > > &multiple_grouping,
                const NamedMatrix<KEnRef_Real> &a_coef,
                const NamedMatrix<KEnRef_Real> &lambda_coef)
        : atom_pairs(atom_pairs), sigma(sigma), multiple_grouping(multiple_grouping), a_coef(a_coef),
        lambda_coef(std::move(lambda_coef)) {

        for (int i = 0; i < atom_pairs.size(); ++i) {
            const auto &atomPair = atom_pairs[i];
            atomNamePairs_to_atomPairIndex[atomPair] = i;
        }
    }
    const NamedVector<KEnRef_Real>& sigmas() const {return *this->sigma;}
    KEnRef_Real getSigma(const std::tuple<std::string, std::string>& atomPairs) const;
    void set_sigmas(const NamedVector<KEnRef_Real>& sigmas) {
        this->sigma = sigmas;
    }
    void setAtomPairs(std::vector<std::tuple<std::string, std::string>> &atomPairs);

    [[nodiscard]] const std::vector<std::tuple<std::string, std::string>> & get_atom_pairs() const { return atom_pairs; }
    [[nodiscard]] const NamedMatrix<KEnRef_Real> &get_a_coef() const { return a_coef; }
    [[nodiscard]] const NamedMatrix<KEnRef_Real> &get_lambda_coef() const { return lambda_coef; }
    [[nodiscard]] const std::vector<std::vector<std::vector<int> > > &get_multiple_grouping() const { return multiple_grouping; }
    [[nodiscard]] const std::optional<std::vector<std::tuple<int, int>>> & get_atomIdPairs_to_sub0Atom_id_pairs_cache() const {
        return atomIdPairs_to_sub0Atom_id_pairs_cache_;
    }
    void set_atomIdPairs_to_sub0Atom_id_pairs_cache(const std::optional<std::vector<std::tuple<int, int>>> &atomIdPairs_to_sub0Atom_id_pairs_cache) {
        this->atomIdPairs_to_sub0Atom_id_pairs_cache_ = atomIdPairs_to_sub0Atom_id_pairs_cache;
    }
};


template<typename KEnRef_Real>
class KEnRef final {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	KEnRef();
    ~KEnRef();
//	KEnRef(const KEnRef &other);
//  KEnRef(KEnRef &&other) noexcept ;
//	KEnRef& operator=(const KEnRef &other);
//	KEnRef& operator=(KEnRef &&other);

    enum lossFunction{POWER_SCALED_LOSS_FUNCTION, LOG_ABS_DIFFERENCE_OVER_OPTIMUM_LOSS_FUNCTION};
    inline static const std::string KC = "kc";

    static std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >
    array_shift(const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &arrays_to_shift, uint nShift,
                int numOmpThreads = 0);

    //First element of the returned tuple is (numPairs, <d_1,d_2,d_3,d_4,d_5,>)
    //Second element of the returned tuple is (numPairs, <  x1d1,x2d1,x3d1,
    //                                                      x1d2,x2d2,x3d2,
    //                                                      x1d3,x2d3,x3d3,
    //                                                      x1d4,x2d4,x3d4,
    //                                                      x1d5,x2d5,x3d5>)
	static std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> >
	r_array_to_d_array(const CoordsMatrixType<KEnRef_Real> &Nxyz, bool gradient = false, bool addEpsilon = false, int numOmpThreads = 0);

    //return tuple where item0 is dipole-dipole interaction tensors (model<pairs, 5_tensor_elements>)
    //item1 is derivatives (It is a vector of 2D Matrix (models<pairId, (5_tensor_elements * XYZ)>).
    static std::tuple<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >, std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > >
    r_array_to_d_array(
        const std::vector<CoordsMatrixType<KEnRef_Real> > &models_Nxyz, //model<pairID, XYZ>
        bool gradient = false, bool addEpsilon = false, int numOmpThreads = 0
			);


    /** Calculate group norm squared from dipole-dipole interaction tensors, and optionally their gradients in the 5 tensor dimensions.
	 *
	 * @param  d_arrays vector (models<pairId, interaction tensor_elements>) with interaction tensors
	 * @param  grouping groupings of models to average interaction tensors (per dipole-dipole interaction pair), i.e. outer list for pairId and inner list for modelId
	 * @param  gradient whether to calculate & return the derivates
	 * @param  numOmpThreads number of OpenMP threads
	 * @return tuple containing
	 *			1) vector of norm squared for each atom pair,
	 *			2) derivates in a vector of Matrix (models<pairId, derivates of interaction tensor_elements>) or empty vector
	 */
    static std::tuple<Eigen::VectorX<KEnRef_Real>, std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > >
    d_array_to_g(
        const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &d_arrays,
        const std::vector<std::vector<int> > &grouping,
        bool gradient = false, int numOmpThreads = 0
    );

    // Calculate group norm squared from dipole-dipole interaction tensors
    static std::tuple<std::vector<Eigen::VectorX<KEnRef_Real> >, std::vector<std::vector<Eigen::Matrix<KEnRef_Real,
        Eigen::Dynamic, 5> > > >
    d_array_to_g_multiple_groupings(
        const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &d_array,
        //vector (models<pairId, tensor_elements>) with interaction tensors
        const std::vector<std::vector<std::vector<int> > > &groupings,
        //groupings of models to average interaction tensors
        bool gradient = false, int numOmpThreads = 0
    );

    static std::tuple<std::vector<Eigen::VectorX<KEnRef_Real> >, std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > > >
    d_array_to_g_matrix(
        const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &d_array,
        //vector (models<pairId, tensor_elements>) with interaction tensors
        const Eigen::MatrixX<int> &grouping_mat,
        //groupings of models to average interaction tensors (per dipole-dipole interaction pair), i.e., outer list for pairId and inner list for modelId
        bool gradient = false, int numOmpThreads = 0);


    static Eigen::MatrixX<KEnRef_Real>
    a_matrix_to_sigma_backprop(
        const Eigen::MatrixX<KEnRef_Real> &d_sigma_d_a_matrix,
        const Eigen::MatrixX<KEnRef_Real> &d_energy_d_sigma,
        int numOmpThreads = 0);

    static std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >
    d_array_to_g_matrix_backprop(
        const std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > > &d_g_matrix_d_d_array,
        const Eigen::MatrixX<KEnRef_Real> &d_energy_d_g_matrix,
        /*int num_interactions, int num_ensembleMembers,*/ int numOmpThreads = 0);

    static std::vector<CoordsMatrixType<KEnRef_Real> >
    r_array_to_d_array_backprop(
        const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > &d_d_array_d_r_array,
        const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &d_energy_d_d_vector,
        int numOmpThreads);


    //Calculate internuclear vectors from atomic coordinates
    //!\param coord_array vector of Nx3 Matrix (models<atoms, xyz>) with atomic coordinates
    //!\param atomId_pairs matrix with each row having the indices of an atom pair (first dimension in `coord_array` Matrices)
    //!\param numOmpThreads number of available MpenMP threads
    //!\return a vector of Matrices (models<pairIndex, xyz>) with internuclear vectors.
    //Atom pair names in R code follow the format `resSeq:Atom-resSeq:Atom` (not implemented here).
    static std::vector<CoordsMatrixType<KEnRef_Real> >
    coord_array_to_r_array(
        const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
        const std::vector<std::tuple<int, int> > &atomId_pairs, int numOmpThreads = 0
    );

    static CoordsMatrixType<KEnRef_Real>
    coord_array_to_r_array(
        const CoordsMatrixType<KEnRef_Real> &coord_array,
        const std::vector<std::tuple<int, int> > &atomId_pairs, int numOmpThreads = 0);

    static std::vector<CoordsMatrixType<KEnRef_Real> >
    coord_array_to_r_array_backprop(
        const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
        const std::vector<std::tuple<int, int> > &atomId_pairs,
        const std::vector<CoordsMatrixType<KEnRef_Real> > &d_energy_d_r_array,
        int numOmpThreads=0);

    static std::vector<std::tuple<int, int> >
    atomNamePairs_2_atomIdPairs(const std::vector<std::tuple<std::string, std::string> > &atomName_pairs,
                                const std::map<std::string, int> &atomNames_2_atomIds, int numOmpThreads = 0);

    //Calculate g value restraint energy from atomic coordinates
    static std::tuple<KEnRef_Real, std::optional<std::vector<CoordsMatrixType<KEnRef_Real> > > >
    coord_array_to_g_energy(
        const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
        //Every vector item is Nx3 Matrix representing atom coordinates of a model.
        const std::vector<std::tuple<std::string, std::string> > &atomName_pairs,
        // Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
        const std::vector<std::vector<std::vector<int> > > &grouping_list,
        // list of lists of integer vectors giving groupings of models to average interaction tensors
        const Eigen::MatrixX<KEnRef_Real> &g0, //target group norm squared values
        const std::map<std::string, int> &atomNames_2_atomIds,
        KEnRef_Real k = 1.0, //force constant
        KEnRef_Real n = 0.25,
        bool gradient = false, int numOmpThreads = 0
    );

    //Calculate g value restraint energy from atomic coordinates
    static std::tuple<KEnRef_Real, std::optional<std::vector<CoordsMatrixType<KEnRef_Real> > > >
    coord_array_to_g_energy(
        const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
        //Every vector item is Nx3 Matrix representing atom coordinates of a model.
        const std::vector<std::tuple<int, int> > &atomId_pairs,
        // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
        const std::vector<std::vector<std::vector<int> > > &grouping_list,
        // list of lists of integer vectors giving groupings of models to average interaction tensors
        const Eigen::MatrixX<KEnRef_Real> &g0, //target group norm squared values
        KEnRef_Real k = 1.0, //force constant
        KEnRef_Real n = 0.25,
        bool gradient = false, int numOmpThreads = 0
    );

    //Calculate g value restraint energy from atomic coordinates
    static std::tuple<KEnRef_Real, std::vector<CoordsMatrixType<KEnRef_Real> > >
    coord_array_to_g_energy_refactored(
        const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
        //Every vector item is Nx3 Matrix representing atom coordinates of a model.
        const std::vector<std::tuple<std::string, std::string> > &atomName_pairs,
        // Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
        const std::vector<std::vector<std::vector<int> > > &grouping_list,
        // list of lists of integer vectors giving groupings of models to average interaction tensors
        const Eigen::MatrixX<KEnRef_Real> &g0, //target group norm squared values
        const std::map<std::string, int> &atomNames_2_atomIds,
        KEnRef_Real k = 1.0, //force constant
        KEnRef_Real n = 0.25,
        bool gradient = false, int numOmpThreads = 0
    );

    //Calculate g value restraint energy from atomic coordinates
    static std::tuple<KEnRef_Real, std::vector<CoordsMatrixType<KEnRef_Real> > >
    coord_array_to_g_energy_refactored(
        const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
        //Every vector item is Nx3 Matrix representing atom coordinates of a model.
        const std::vector<std::tuple<int, int> > &atomId_pairs,
        // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
        const std::vector<std::vector<std::vector<int> > > &grouping_list,
        // list of lists of integer vectors giving groupings of models to average interaction tensors
        const Eigen::MatrixX<KEnRef_Real> &g0, //target group norm squared values
        KEnRef_Real k = 1.0, //force constant
        KEnRef_Real n = 0.25,
        bool gradient = false, int numOmpThreads = 0
    );


    static Eigen::MatrixX<KEnRef_Real>
    coord_array_to_g_matrix(
        const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
        //Every vector item is Nx3 Matrix representing atom coordinates of a model.
        const std::vector<std::tuple<int, int> > &atomId_pairs,
        // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
        const std::vector<std::vector<std::vector<int> > > &grouping_list,
        // list of lists of integer vectors giving groupings of models to average interaction tensors
        int numOmpThreads = 0
    );

    static NamedMatrix<KEnRef_Real>
    g_matrix_to_a_matrix(const Eigen::MatrixX<KEnRef_Real> &g_matrix, const NamedMatrix<KEnRef_Real> &a_coef,
        int numOmpThreads = 0);

    // Calculate restraint energy from group norm squared values
    // returns restraint energy calculated using \eqn{k*(g^n -g0^n)^2} //TODO correct the function.
    static std::tuple<Eigen::MatrixX<KEnRef_Real>, std::optional<Eigen::MatrixX<KEnRef_Real> > >
    power_scaled_loss_function(
    const Eigen::MatrixX<KEnRef_Real>& g, // current group norm squared values
    const Eigen::MatrixX<KEnRef_Real>& g0, // target group norm squared values.
    KEnRef_Real k = 1.0, // force constant
    KEnRef_Real n = 0.25, // correction power.
    bool gradient = false,
    int numOmpThreads = 0
    );

    static std::tuple<Eigen::MatrixX<KEnRef_Real>, std::optional<Eigen::MatrixX<KEnRef_Real> > >
    log_abs_diff_over_optimum_loss_function(
    Eigen::MatrixX<KEnRef_Real> g, // current group norm squared values
    Eigen::MatrixX<KEnRef_Real> g0, // target group norm squared values.
    KEnRef_Real k = 1.0, // force constant
    bool gradient = false,
    int numOmpThreads = 0
    );

    //Collects list/vector of norm squared of all groups in a single matrix (num_pairIds, num_models (or num of grouping vectors?))
    static Eigen::MatrixX<KEnRef_Real>
    vectorOfVectors_to_Matrix(std::vector<Eigen::VectorX<KEnRef_Real> > g_vect, int numOmpThreads = 0);

    static void saturate(CoordsMatrixType<KEnRef_Real> &derivatives_rectified, KEnRef_Real thresholdSquared, int numOmpThreads = 0);

    static Eigen::VectorX<KEnRef_Real>
    s2OrderParams(
        const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
        //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
        const std::vector<std::tuple<int, int> > &atomId_pairs,
        // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
        int numOmpThreads);

    // lambda_vector <- -colSums(rates[rownames(spec_den_data_list[[i]][["lambda_coef"]])]*spec_den_data_list[[i]][["lambda_coef"]])
    static NamedRowVector<KEnRef_Real>
    calculateLambdaVector(const SpecDenData<KEnRef_Real> &currentSpecDenData, const NamedRowVector<KEnRef_Real> &rates, int numOmpThreads=0);

    static std::tuple<std::vector<NamedVector<KEnRef_Real>>, std::optional<std::vector<Eigen::MatrixX<KEnRef_Real> > > >
        coord_array_to_sigma(
            const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
            // Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
            const NamedRowVector<KEnRef_Real> &rates,
            const std::vector<SpecDenData<KEnRef_Real> > &spec_den_data_list,
            KEnRef_Real proton_mhz,
            std::map<std::string, int> atomNames_2_atomIds, bool gradient=false, int numOmpThreads=0);


    static Eigen::MatrixX<KEnRef_Real>
    g_matrix_to_a_matrix_backprop(
        const Eigen::MatrixX<KEnRef_Real>& a_coef,
        const Eigen::MatrixX<KEnRef_Real>& d_energy_d_a_matrix,
        int numOmpThreads = 0);

    static std::tuple<KEnRef_Real, std::optional<std::vector<CoordsMatrixType<KEnRef_Real> > > >
    coord_array_to_sigma_energy(
        std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
        // Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
        const NamedRowVector<KEnRef_Real> &rates,
        const std::vector<SpecDenData<KEnRef_Real> > &spec_den_data_list,
        KEnRef_Real proton_mhz,
        KEnRef_Real k,
        KEnRef_Real n, const std::map<std::string, int> &atomNames_2_atomIds, bool gradient, int numOmpThreads = 0 /*,
        lossFunction lossFunc = POWER_SCALED_LOSS_FUNCTION*/);

    /**Calculate sigma from a matrix of a values
     *
     * @param a_matrix matrix of a values with columns associated with eigenvalues
     * @param lambda_prime_vec eigenvalues augmented with tumbling rate
     * @param proton_mhz spectrometer proton field strength in MHz
     * @param gradient a logical value indicating whether to calculate the derivative
     * @return
     */
    [[deprecated("Use a_matrix_to_sigma() instead.")]]
    std::tuple<NamedRowVector<KEnRef_Real>, std::optional<NamedMatrix<KEnRef_Real>>>
    static a_matrix_to_sigma_horizontal(const NamedMatrix<KEnRef_Real>& a_matrix,
        const NamedRowVector<KEnRef_Real>& lambda_prime_vec,
        KEnRef_Real proton_mhz, bool gradient = false);

    /**Calculate sigma from a matrix of a values
     *
     * @param a_matrix matrix of a values with columns associated with eigenvalues
     * @param lambda_prime_vec eigenvalues augmented with tumbling rate
     * @param proton_mhz spectrometer proton field strength in MHz
     * @param gradient a logical value indicating whether to calculate the derivative
     * @return
     */
    std::tuple<NamedVector<KEnRef_Real>, std::optional<NamedMatrix<KEnRef_Real>>>
    static a_matrix_to_sigma(const NamedMatrix<KEnRef_Real>& a_matrix,
        const NamedRowVector<KEnRef_Real>& lambda_prime_vec,
        KEnRef_Real proton_mhz, bool gradient = false);
};

#endif /* KENREF_H_ */
