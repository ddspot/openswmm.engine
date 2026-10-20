// SPDX-License-Identifier: Apache-2.0
/**
 * @file MeshVegetation.hpp
 * @brief The vegetated mesh cell (fork `urbaneva/stage1-vegetat`): UrbanEVA's
 *        vegetation law as the top boundary of a `[2D_AQUIFER]` column.
 *
 * Per `[2D_VEGETATION]` row (scope `*` / TAG / CELL) a cell carries a canopy
 * (Braden I-3 interception on the cell's rain, Ei-2 canopy evaporation) and
 * its evaporation demand becomes `Ep · [SCF·(1 − f) + (1 − SCF)·g]` with
 * `Ep = E_pot · Ks(season) · Kc · Kc_pattern[month]`, `f` the wet-canopy
 * fraction and `g` the Es-2 stress on the root zone's water content
 * (`SubsurfaceSolver::rootZoneTheta`, `ROOT_DEPTH`). The marcher then pays
 * ponded-water evaporation from that demand and the subsurface takes the rest
 * at the column top (`GW_ET BOUNDARY_ET / BOTH`) — nothing downstream changes.
 * Intercepted rain never reaches the surface; the canopy terms are reported
 * beside the 2D continuity block, outside its balance.
 */
#ifndef OPENSWMM_ENGINE_2D_MESH_VEGETATION_HPP
#define OPENSWMM_ENGINE_2D_MESH_VEGETATION_HPP
#include <cstdint>
#include <string>
#include <vector>

namespace openswmm { struct SimulationContext; }

namespace openswmm::twoD {
struct MeshData;
struct SurfaceStateData;
struct SubsurfaceConfig;
struct GwUnitFactors;
class SubsurfaceSolver;

class MeshVegetation {
public:
    /// Resolve the rows onto cells (needs the aquifer: FC defaults to the
    /// column's drainage field capacity, WP to θ_r). Empty return = OK.
    std::string initialize(const MeshData& mesh, const SubsurfaceConfig& cfg, const GwUnitFactors& uf,
                           const SimulationContext& ctx, const SubsurfaceSolver& gw,
                           std::vector<std::string>& warnings);
    bool active() const noexcept { return active_; }

    /// Forcing-cadence refresh: Ep, Smax, SCF, root-zone stress, demand →
    /// `state.evap_rate` on vegetated cells that are not forced through the API.
    void refresh(const SimulationContext& ctx, const SubsurfaceSolver& gw, SurfaceStateData& state, bool evap_on);
    /// Before the marcher advances: intercept from the gross rain, write the net
    /// rain the surface sees. `rain_gross` is the refreshed per-cell rate (m/s).
    void beforeAdvance(SurfaceStateData& state, const std::vector<double>& rain_gross, double dt);
    /// After the advance: the canopy evaporates Ep · f · SCF over `dt`.
    void afterAdvance(double dt);

    double interceptedVolume() const noexcept { return intercepted_m3_; }   ///< cumulative, m³
    double canopyEvapVolume() const noexcept { return canopy_evap_m3_; }    ///< cumulative, m³
    double canopyStorageVolume() const noexcept;                            ///< now, m³

private:
    bool active_ = false;
    int n_ = 0;
    std::vector<int8_t> enabled_;
    std::vector<double> area_, laimean_, kc_, sl_m_, nfc_, fc_, wp_, ks0_, ks1_, ks2_;
    std::vector<int> lai_idx_, kc_idx_;
    // state
    std::vector<double> canopy_m_, smax_m_, scf_, ep_, demand_, stress_, wet_;
    double intercepted_m3_ = 0.0, canopy_evap_m3_ = 0.0;
};

}  // namespace openswmm::twoD
#endif
