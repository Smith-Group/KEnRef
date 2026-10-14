/* KEnRef must subscribe to GROMACS's simulation-setup notifiers exactly ONCE.
 *
 * mdrun.cpp subscribes the module itself, because GROMACS 2025 never does it for a module added from
 * outside. GROMACS 2026 does, so there the module is subscribed a second time on the same notifiers,
 * and every callback would be registered (and run) twice. These tests count the subscriptions through
 * the diagnostic printed at the subscribe call, so they fail if any path subscribes more than once.
 *
 * The multisim subscription needs a multi-simulation record, which is built by hand below (no MPI
 * calls are made, so no MPI_Init is needed). */
#include "config.h" // GMX_LIB_MPI

#include <filesystem>
#include <memory>
#include <string>

#include "gtest/gtest.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdrun/simulationcontext.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/version.h"
#include "gmxinterface/KEnRefMDModule.h"
#include "gmxinterface/gmxwrapper.h"

namespace {

const std::string kSubscribed = "subscribed to the multisim record";

int occurrences(const std::string& text, const std::string& what) {
    int n = 0;
    for (auto pos = text.find(what); pos != std::string::npos; pos = text.find(what, pos + what.size())) {
        ++n;
    }
    return n;
}

// One replica of a two-replica run, set up the way mdrun.cpp sets it up.
class SubscribeOnceTest : public ::testing::Test {
protected:
    gmx::SimulationContext context_{GMX_LIB_MPI ? MPI_COMM_WORLD : MPI_COMM_NULL, {}};
    std::shared_ptr<gmx::IMDModule> module_;
    KEnRefMDModule* raw_ = nullptr;

    void SetUp() override {
        const auto res = std::filesystem::path(__FILE__).parent_path().parent_path() / "res" / "10nsstart+fitting";
        Settings::indexFileName = (res / "KEnRefAtomIndex.ndx").string();
        Settings::refFileName = (res / "00ns.pdb").string();
        context_.multiSimulation_ = std::make_unique<gmx_multisim_t>(2, 0, MPI_COMM_NULL, MPI_COMM_NULL);
        module_ = KEnRefModuleInfo::create();
        raw_ = dynamic_cast<KEnRefMDModule*>(module_.get());
        raw_->setSimulationContext(&context_);
    }
};

} // namespace

// What GROMACS 2026 does to mdrun.cpp's module, without depending on the GROMACS version.
TEST_F(SubscribeOnceTest, RepeatOnTheSameNotifiersIsIgnored) {
    gmx::MDModulesNotifiers notifiers;
    testing::internal::CaptureStdout();
    raw_->subscribeToSimulationSetupNotifications(&notifiers);
    raw_->subscribeToSimulationSetupNotifications(&notifiers);
    EXPECT_EQ(occurrences(testing::internal::GetCapturedStdout(), kSubscribed), 1);
}

// The count that matters: a callback registered twice RUNS twice. The topology callback is
// unconditional, so one notification must reach it exactly once.
TEST_F(SubscribeOnceTest, OneNotificationRunsEachCallbackOnce) {
    gmx::MDModulesNotifiers notifiers;
    raw_->subscribeToSimulationSetupNotifications(&notifiers);
    raw_->subscribeToSimulationSetupNotifications(&notifiers);
    const gmx_mtop_t mtop;
    testing::internal::CaptureStdout();
    notifiers.simulationSetupNotifier_.notify(mtop);
    EXPECT_EQ(occurrences(testing::internal::GetCapturedStdout(), "KEnRef: topology received"), 1);
}

// The guard is per notifiers object, not a blanket "subscribe once ever".
TEST_F(SubscribeOnceTest, DifferentNotifiersAreEachSubscribed) {
    gmx::MDModulesNotifiers first;
    gmx::MDModulesNotifiers second;
    testing::internal::CaptureStdout();
    raw_->subscribeToSimulationSetupNotifications(&first);
    raw_->subscribeToSimulationSetupNotifications(&second);
    EXPECT_EQ(occurrences(testing::internal::GetCapturedStdout(), kSubscribed), 2);
}

// mdrun.cpp's own call followed by the runner's dispatch, through the real GROMACS MDModules: one
// subscription on 2025 (GROMACS skips added modules) and, with the guard, one on 2026 as well.
TEST_F(SubscribeOnceTest, MdrunAndGromacsTogetherSubscribeOnce) {
    gmx::MDModules mdModules;
    testing::internal::CaptureStdout();
    raw_->subscribeToSimulationSetupNotifications(const_cast<gmx::MDModulesNotifiers*>(&mdModules.notifiers()));
#if GMX_VERSION >= 20260000
    mdModules.add("KEnRef", module_);
#else
    mdModules.add(module_);
#endif
    mdModules.subscribeToSimulationSetupNotifications();
    EXPECT_EQ(occurrences(testing::internal::GetCapturedStdout(), kSubscribed), 1);
}
