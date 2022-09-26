#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "instruments.hpp"

/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> discard(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_discard, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> discard(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_discard, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi, const std::vector<idx>& subsys,
///       const std::vector<idx>& dims)
TEST(qpp_ip, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       ip(const Eigen::MatrixBase<Derived>& phi,
///       const Eigen::MatrixBase<Derived>& psi, const std::vector<idx>& subsys,
///       idx d = 2)
TEST(qpp_ip, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(
///       const Eigen::MatrixBase<Derived>& A, const cmat& U)
TEST(qpp_measure, Orthonormal) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(
///       const Eigen::MatrixBase<Derived>& A, const cmat& V,
///       const std::vector<idx>& target, const std::vector<idx>& dims,
///       bool destructive = true)
TEST(qpp_measure, RankOneQudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>> measure(
///       const Eigen::MatrixBase<Derived>& A, const cmat& V,
///       const std::vector<idx>& target, idx d = 2, bool destructive = true)
TEST(qpp_measure, RankOneQubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks)
TEST(qpp_measure, KrausInitList) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks, const std::vector<idx>& target,
///       const std::vector<idx>& dims, bool destructive = true)
TEST(qpp_measure, SubsysKrausInitList) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::initializer_list<cmat>& Ks, const std::vector<idx>& target,
///       idx d = 2, bool destructive = true)
TEST(qpp_measure_DeathTest, SubsysKrausInitListQubits)
{
    using namespace qpp::literals;

    auto const psi = Eigen::Vector4cd::Unit(0).eval();

    EXPECT_EXIT({ qpp::measure(psi, { 00_prj, 01_prj, 10_prj, 11_prj }, { 0, 1 }, 2, false); exit(0);},
        ::testing::ExitedWithCode(0), "");
}
TEST(qpp_measure_DeathTest, SubsysKrausInitListQubitsDestructive)
{
    using namespace qpp::literals;

    auto const psi = Eigen::Vector4cd::Unit(0).eval();

    EXPECT_EXIT({ qpp::measure(psi, { 00_prj, 01_prj, 10_prj, 11_prj }, { 0, 1 }, 2, true); exit(0);},
        ::testing::ExitedWithCode(0), "");
}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks)
TEST(qpp_measure, KrausVector) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks, const std::vector<idx>& target,
///       const std::vector<idx>& dims, bool destructive = true)
TEST(qpp_measure, SubsysKrausVector) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       std::tuple<idx, std::vector<double>, std::vector<cmat>>
///       measure(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<cmat>& Ks, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_measure, SubsysKrausVectorQubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::tuple<std::vector<idx>, double, cmat>
///       measure_seq(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2, bool destructive = true)
TEST(qpp_measure_seq, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::tuple<std::vector<idx>, double, cmat>
///       measure_seq(const Eigen::MatrixBase<Derived>& A,
///       std::vector<idx> target, std::vector<idx> dims,
///       bool destructive = true)
TEST(qpp_measure_seq, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> reset(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       const std::vector<idx>& dims)
TEST(qpp_reset, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived>
///       dyn_mat<typename Derived::Scalar> reset(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
///       idx d = 2)
TEST(qpp_reset, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::vector<idx>
///       sample(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_sample, Qudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::vector<idx>
///       sample(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_sample, Qubits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::map<std::vector<idx>, idx>
///       sample(idx num_samples, const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, const std::vector<idx>& dims)
TEST(qpp_sample, MultipleQudits) {}
/******************************************************************************/
/// BEGIN template <typename Derived> std::map<std::vector<idx>, idx>
///       sample(idx num_samples,const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& target, idx d = 2)
TEST(qpp_sample, MultipleQubits) {}
/******************************************************************************/
