/*
 * KEnRef.h
 *
 *  Created on: May 8, 2023
 *      Author: amr
 */

#ifndef KENREF_H_
#define KENREF_H_

#include <utility>
#include <vector>
#include <map>
#include <tuple>
#include <string>
#include <optional>
#include <cassert>

#include <Eigen/Dense>

#include "Table.h"
// #include "../config/KEnRefConfig.h"

#ifdef KENREF_DOUBLE
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

/**
 * Shared base for the spectral-density data variants. It holds the parts that are
 * common to both the sigma-based data (see read_spec_den_data() in ke.R) and the
 * relaxation-based data (see read_spec_den_relax_data() in ke.R): the atom pairs,
 * the multiple groupings, and the `a`/`lambda` coefficient matrices, together with
 * the atom-pair-name -> index lookup and the cached atom-id pairs.
 */
template<typename KEnRef_Real>
class SpecDenDataBase {
protected:
    std::vector<std::tuple<std::string, std::string> > atom_pairs{};
    std::vector<std::vector<std::vector<int> > > multiple_grouping{};
    NamedMatrix<KEnRef_Real> a_coef;
    NamedMatrix<KEnRef_Real> lambda_coef;
    std::map<std::tuple<std::string, std::string>, size_t> atomNamePairs_to_atomPairIndex;
    std::optional<std::vector<std::tuple<int, int> > > atomIdPairs_to_sub0Atom_id_pairs_cache_ = std::nullopt;

public:
    SpecDenDataBase(const std::vector<std::tuple<std::string, std::string> > &atom_pairs,
                    const std::vector<std::vector<std::vector<int> > > &multiple_grouping,
                    const NamedMatrix<KEnRef_Real> &a_coef,
                    const NamedMatrix<KEnRef_Real> &lambda_coef)
        : atom_pairs(atom_pairs), multiple_grouping(multiple_grouping), a_coef(a_coef),
        lambda_coef(std::move(lambda_coef)) {

        for (int i = 0; i < atom_pairs.size(); ++i) {
            const auto &atomPair = atom_pairs[i];
            atomNamePairs_to_atomPairIndex[atomPair] = i;
        }
    }
    virtual ~SpecDenDataBase() = default;

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

/**
 * Sigma-based spectral-density data. Mirrors read_spec_den_data() in ke.R: it
 * augments the shared data with the per-atom-pair `sigma` target vector.
 */
template<typename KEnRef_Real>
class SpecDenData : public SpecDenDataBase<KEnRef_Real> {
    std::optional<NamedVector<KEnRef_Real>> sigma = std::nullopt;
    KEnRef_Real getSigma(size_t atomPairId) const;

public:
    SpecDenData(const std::vector<std::tuple<std::string, std::string> > &atom_pairs,
                const std::optional<NamedVector<KEnRef_Real> > &sigma,
                const std::vector<std::vector<std::vector<int> > > &multiple_grouping,
                const NamedMatrix<KEnRef_Real> &a_coef,
                const NamedMatrix<KEnRef_Real> &lambda_coef)
        : SpecDenDataBase<KEnRef_Real>(atom_pairs, multiple_grouping, a_coef, lambda_coef), sigma(sigma) {}

    const NamedVector<KEnRef_Real>& sigmas() const {return *this->sigma;}
    KEnRef_Real getSigma(const std::tuple<std::string, std::string>& atomPairs) const;
    void set_sigmas(const NamedVector<KEnRef_Real>& sigmas) {
        this->sigma = sigmas;
    }
};

/**
 * One spectral-density term array for a single relaxation rate. Mirrors the
 * `spec_den_term_array` member of an R `relax_entry`, but stored as the two 2-D
 * slices the calculation actually consumes: a_matrix_to_relax() in ke.R immediately
 * splits the R (pairs, terms, 2) array into a coef matrix and a freq matrix, so we
 * store them pre-split.
 *
 * Both are (n_data_rows x terms) with column names = the frequency-term labels
 * (e.g. "wHmwN", "wN", "2wH"). Rows are positional: atom-pair identity lives once in
 * SpecDenDataBase and is not duplicated here.
 *
 * Note: freq.array().square() is loop-invariant in a_matrix_to_relax() and may be
 * precomputed once per block at calculation time.
 */
template<typename KEnRef_Real>
struct SpecDenTermArray {
    NamedMatrix<KEnRef_Real> coef;   // n_data_rows x terms
    NamedMatrix<KEnRef_Real> freq;   // n_data_rows x terms (same column labels as coef)
};

/**
 * One entry of `relax_data_list`, for a single relaxation rate. Mirrors an R
 * `relax_entry = list(value, k, spec_den_term_array)`, with the rate name (R's list
 * `names` attribute) carried alongside so the flattened `*_atom_relax.csv` column
 * names can be regenerated as `<rate_name>_<term>_<coef|freq>`, `<rate_name>_value`,
 * and `<rate_name>_k`.
 *
 * `value` and the term arrays have n_data_rows rows; `k` has length 1 or n_data_rows
 * (length 1 broadcasts). All are positional/nameless on the row axis.
 */
template<typename KEnRef_Real>
struct RelaxEntry {
    std::string rate_name;                              // e.g. "r1_400"
    std::optional<NamedVector<KEnRef_Real>> value;      // n_data_rows x 1  (<rate>_value), or nullopt
    std::optional<NamedVector<KEnRef_Real>> k;          // 1 or n_data_rows  (<rate>_k), or nullopt
    SpecDenTermArray<KEnRef_Real> spec_den_term_array;

    // Forwarding accessors to the term-array slices
    // (R: relax_entry[["spec_den_term_array"]][, , "coef"|"freq"]).
    [[nodiscard]] const NamedMatrix<KEnRef_Real> &coef() const { return spec_den_term_array.coef; }
    [[nodiscard]] const NamedMatrix<KEnRef_Real> &freq() const { return spec_den_term_array.freq; }
};

/**
 * Relaxation-based spectral-density data. Mirrors read_spec_den_relax_data() in
 * ke.R (the full `spec_den_relax_data` object). The shared `a_coef`/`lambda_coef`
 * of the base correspond to `a_int_coef`/`lambda_int_coef` there.
 *
 * `unit` is the dipole-dipole unit-vector flag forwarded to r_array_to_d_array()
 * (true for e.g. a fixed backbone N-H bond). On disk it is encoded by the first two
 * `*_atom_relax.csv` columns ending in `_unit`. It is INDEPENDENT of `handleNames`:
 * the atom-pair identifiers are still atom names and are still normalized when
 * handleNames is set.
 *
 * Wrapping: only the first n_data_rows atom pairs carry relaxation data; the full
 * atom_pairs list may be a whole multiple of that (the trailing rows alias back,
 * exactly as sigma does). After the array_shift pipeline the predicted rates collapse
 * to n_data_rows, so each rate's value/k/coef/freq aligns row-for-row with the
 * a-matrix. Use blockRow() to map a full atom-pair index/name to its data row (idx % N).
 */
template<typename KEnRef_Real>
class SpecDenRelaxData : public SpecDenDataBase<KEnRef_Real> {
    bool unit_ = false;
    std::vector<RelaxEntry<KEnRef_Real> > relax_data_list_{};
    std::map<std::string, size_t> rateName_to_index_{};
    size_t n_data_rows_ = 0;
    std::string atom_relax_header_{};   // verbatim *_atom_relax.csv header line (byte-exact round-trip)

public:
    SpecDenRelaxData(const std::vector<std::tuple<std::string, std::string> > &atom_pairs,
                     const bool unit,
                     std::vector<RelaxEntry<KEnRef_Real> > relax_data_list,
                     const size_t n_data_rows,
                     std::string atom_relax_header,
                     const std::vector<std::vector<std::vector<int> > > &multiple_grouping,
                     const NamedMatrix<KEnRef_Real> &a_int_coef,
                     const NamedMatrix<KEnRef_Real> &lambda_int_coef)
        : SpecDenDataBase<KEnRef_Real>(atom_pairs, multiple_grouping, a_int_coef, lambda_int_coef),
          unit_(unit),
          relax_data_list_(std::move(relax_data_list)),
          n_data_rows_(n_data_rows),
          atom_relax_header_(std::move(atom_relax_header)) {
        for (size_t i = 0; i < relax_data_list_.size(); ++i) {
            rateName_to_index_[relax_data_list_[i].rate_name] = i;
        }
        // Wrapping invariant: the full atom-pair list is a whole multiple of the data rows.
        assert((n_data_rows_ == 0 || atom_pairs.size() % n_data_rows_ == 0)
               && "atom_pairs count must be a whole multiple of n_data_rows");
    }

    [[nodiscard]] bool is_unit() const { return unit_; }
    [[nodiscard]] size_t n_data_rows() const { return n_data_rows_; }
    [[nodiscard]] const std::string &get_atom_relax_header() const { return atom_relax_header_; }
    [[nodiscard]] const std::vector<RelaxEntry<KEnRef_Real> > &get_relax_data_list() const { return relax_data_list_; }

    // R: relax_data_list[[rate_name]]
    [[nodiscard]] bool has_rate(const std::string &rate_name) const {
        return rateName_to_index_.find(rate_name) != rateName_to_index_.end();
    }
    [[nodiscard]] const RelaxEntry<KEnRef_Real> &rate(const std::string &rate_name) const {
        return relax_data_list_.at(rateName_to_index_.at(rate_name));
    }

    // Map a full atom-pair index/name to its data row, mirroring SpecDenData::getSigma's `% size` wrapping.
    [[nodiscard]] size_t blockRow(const size_t pairIndex) const { return pairIndex % n_data_rows_; }
    [[nodiscard]] size_t blockRow(const std::tuple<std::string, std::string> &pair) const {
        return this->atomNamePairs_to_atomPairIndex.at(pair) % n_data_rows_;
    }

    // Debugging convenience: attach atom-pair identifiers as row names to a positional
    // n_data_rows-row matrix (e.g. a coef/freq slice). Names are regenerated on demand
    // from the base atom_pairs rather than stored per substructure.
    [[nodiscard]] NamedMatrix<KEnRef_Real> withPairNames(const NamedMatrix<KEnRef_Real> &m) const {
        NamedMatrix<KEnRef_Real> out(m);
        std::vector<std::string> rowNames;
        rowNames.reserve(n_data_rows_);
        for (size_t r = 0; r < n_data_rows_; ++r) {
            const auto &p = this->atom_pairs.at(r);
            rowNames.emplace_back(std::get<0>(p) + " / " + std::get<1>(p));
        }
        out.setRowNames(rowNames);
        return out;
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

    enum class lossFunction: int {POWER_SCALED_LOSS_FUNCTION, LOG_ABS_DIFFERENCE_OVER_OPTIMUM_LOSS_FUNCTION};
    enum class energyModel: int {UNKNOWN, SIGMA, PLATEAUS};
    const static inline std::map<std::string, energyModel> energyModelMap{{"UNKNOWN", energyModel::UNKNOWN}, {"SIGMA", energyModel::SIGMA}, {"PLATEAUS", energyModel::PLATEAUS}};
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
	// Computes the d_array (and optional gradient) for a contiguous row-range [startRow, startRow+numRows)
	// of Nxyz. The defaults (startRow=0, numRows=-1) cover the whole matrix, reproducing the original
	// whole-array behavior exactly; passing a sub-range processes just that batch of pair-rows (bit-for-bit
	// identical, since every operation is per-row). The returned matrices have `numRows` (the resolved
	// block length) rows. The vector overload below uses this to fill its result in independent row-blocks.
	static std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> >
	r_array_to_d_array(const CoordsMatrixType<KEnRef_Real> &Nxyz, bool unit = false, bool gradient = false, bool addEpsilon = false, int numOmpThreads = 0,
	                   Eigen::Index startRow = 0, Eigen::Index numRows = -1);

    //return tuple where item0 is dipole-dipole interaction tensors (model<pairs, 5_tensor_elements>)
    //item1 is derivatives (It is a vector of 2D Matrix (models<pairId, (5_tensor_elements * XYZ)>).
    static std::tuple<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >, std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > >
    r_array_to_d_array(
        const std::vector<CoordsMatrixType<KEnRef_Real> > &models_Nxyz, //model<pairID, XYZ>
        bool unit = false, bool gradient = false, bool addEpsilon = false, int numOmpThreads = 0
			);

    // Buffer-reuse variant of the above: fills the caller-provided ret1/ret2 (resized in place) instead
    // of returning fresh vectors. A caller that keeps ret1/ret2 alive across calls (same model count and
    // pair count) avoids re-allocating / re-faulting the [N×5] and [N×15] results every call — the hot
    // per-step path for MD refinement. ret2 is filled only when `gradient` (else sized to 0 rows).
    static void
    r_array_to_d_array_into(
        const std::vector<CoordsMatrixType<KEnRef_Real> > &models_Nxyz,
        bool unit, bool gradient, bool addEpsilon, int numOmpThreads,
        std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &ret1,
        std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > &ret2);


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
    calculateLambdaVector(const SpecDenDataBase<KEnRef_Real> &currentSpecDenData, const NamedRowVector<KEnRef_Real> &rates, int numOmpThreads=0);

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

    /** Calculate relaxation rates from internal-motion amplitudes/eigenvalues and overall-tumbling modes.
     *
     * Port of a_matrix_to_relax() in ke.R. For each overall mode i and internal eigenvalue j, with
     * lambda_prime = lambda_int_vec[j] + lambda_overall_vec[i], accumulates over pairs:
     *   value += a_int_matrix[:,j] * ( -a_overall_matrix[:,i] * lambda_prime
     *                                  * rowSum_t( coef[:,t] / (lambda_prime^2 + freq[:,t]^2) ) )
     * The two spectral-density slices coef/freq come pre-split in `spec_den_term_array`.
     *
     * @return tuple of (relaxation rate per pair [N], optional gradient d(relax)/d(a_int) [N x n_int]).
     */
    static std::tuple<Eigen::VectorX<KEnRef_Real>, std::optional<Eigen::MatrixX<KEnRef_Real>>>
    a_matrix_to_relax(
        const Eigen::MatrixX<KEnRef_Real>& a_int_matrix,
        const Eigen::RowVectorX<KEnRef_Real>& lambda_int_vec,
        const Eigen::MatrixX<KEnRef_Real>& a_overall_matrix,
        const Eigen::RowVectorX<KEnRef_Real>& lambda_overall_vec,
        const SpecDenTermArray<KEnRef_Real>& spec_den_term_array,
        bool gradient = false, int numOmpThreads = 0);

    /** Back-propagate the relaxation rate to the internal-motion amplitude matrix.
     *
     * Port of a_matrix_to_relax_backprop() in ke.R: elementwise
     *   d_energy_d_a_int[:,j] = d_relax_d_a_matrix[:,j] * d_energy_d_relax[:]
     * (the per-pair upstream gradient broadcasts across the n_int columns).
     *
     * @param d_relax_d_a_matrix [pairs x n_int] gradient returned by a_matrix_to_relax().
     * @param d_energy_d_relax   [pairs] upstream gradient of the energy w.r.t. each relaxation rate.
     * @return [pairs x n_int] contribution to d(energy)/d(a_int_matrix).
     */
    static Eigen::MatrixX<KEnRef_Real>
    a_matrix_to_relax_backprop(
        const Eigen::MatrixX<KEnRef_Real>& d_relax_d_a_matrix,
        const Eigen::VectorX<KEnRef_Real>& d_energy_d_relax,
        int numOmpThreads = 0);

    /** Calculate overall rotational-diffusion tumbling modes from the diffusion tensor.
     *
     * Port of dxyz_dunit_to_overall_modes() in ke.R, 3-D-array autocorrelation path only
     * (dunit_b == dunit_a). The coord pipeline always passes the array_shift()-ed unit-tensor
     * array, so the regularized matrix path is not needed here.
     *
     * `dunit_a_array_shifted` is models<pairs, 5> (the shifted unit rank-2 tensors). For each
     * pair the overall-mode amplitudes are averaged over models:
     *   a_overall[p,i] = mean_m ( (d_pm . v_i)^2 ),  v_i = rank-2 diffusion eigenvectors.
     * The number of columns collapses to the diffusion symmetry: isotropic (1), axially
     * symmetric Dx==Dy (3), or fully anisotropic (5). lambda_overall_vec holds the matching
     * negative decay rates.
     *
     * @return tuple (a_overall_matrix [pairs x nOverall], lambda_overall_vec [nOverall]).
     */
    static std::tuple<Eigen::MatrixX<KEnRef_Real>, Eigen::RowVectorX<KEnRef_Real>>
    dxyz_dunit_to_overall_modes(
        const Eigen::Matrix<KEnRef_Real, 3, 1>& dxyz_vec,
        const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>>& dunit_a_array_shifted,
        KEnRef_Real tol = static_cast<KEnRef_Real>(1.4901161193847656e-08), // sqrt(DBL_EPSILON)
        int numOmpThreads = 0);

    /** Calculate per-pair relaxation rates from atomic coordinates.
     *
     * Port of coord_array_to_relax() in ke.R. For each SpecDenRelaxData substructure it runs the
     * shared front-half pipeline (r_array -> d_array -> array_shift -> g_matrix -> a_int_matrix,
     * lambda_int via calculateLambdaVector), derives the overall tumbling modes from the unit
     * dipole tensors, then evaluates a_matrix_to_relax() for every rate in relax_data_list.
     *
     * @return vector (one per substructure) of [n_data_rows x n_rates] relaxation-rate matrices,
     *         column order matching relax_data_list.
     */
    static std::vector<NamedMatrix<KEnRef_Real>>
    coord_array_to_relax(
        const std::vector<CoordsMatrixType<KEnRef_Real>>& coord_array,
        const NamedRowVector<KEnRef_Real>& rates,
        const std::vector<SpecDenRelaxData<KEnRef_Real>>& spec_den_relax_data_list,
        const std::map<std::string, int>& atomNames_2_atomIds,
        int numOmpThreads = 0);

    /** Calculate the relaxation-rate restraint energy (and optional gradient) from coordinates.
     *
     * Port of coord_array_to_relax_energy() in ke.R using power_scaled_loss_function. Flattens the
     * predicted and target (value) relaxation rates plus per-rate force constants k, applies the
     * loss, and sums. The gradient (w.r.t. coord_array, scaled by 1e-10) propagates through the
     * internal-motion amplitudes only — overall tumbling modes are treated as constant, exactly as
     * documented for the R reference.
     *
     * @return tuple (scalar energy, optional gradient as models<atoms,3>).
     */
    static std::tuple<KEnRef_Real, std::optional<std::vector<CoordsMatrixType<KEnRef_Real>>>>
    coord_array_to_relax_energy(
        std::vector<CoordsMatrixType<KEnRef_Real>>& coord_array,
        const NamedRowVector<KEnRef_Real>& rates,
        const std::vector<SpecDenRelaxData<KEnRef_Real>>& spec_den_relax_data_list,
        KEnRef_Real k, KEnRef_Real n,
        const std::map<std::string, int>& atomNames_2_atomIds,
        bool gradient = false, int numOmpThreads = 0);
};

#endif /* KENREF_H_ */
