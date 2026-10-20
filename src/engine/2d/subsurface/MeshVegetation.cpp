// SPDX-License-Identifier: Apache-2.0
#include "MeshVegetation.hpp"
#include <algorithm>
#include <cmath>
#include "SubsurfaceData.hpp"
#include "SubsurfaceSections.hpp"
#include "SubsurfaceSolver.hpp"
#include "../data/MeshData.hpp"
#include "../data/SurfaceStateData.hpp"
#include "../../core/SimulationContext.hpp"
#include "../../core/DateTime.hpp"
#include "../../hydrology/VegetationLaw.hpp"

namespace openswmm::twoD {

std::string MeshVegetation::initialize(const MeshData& mesh, const SubsurfaceConfig& cfg, const GwUnitFactors& uf,
                                       const SimulationContext& ctx, const SubsurfaceSolver& gw,
                                       std::vector<std::string>& warnings) {
    active_ = false;
    if (cfg.vegetation.empty()) return {};
    if (!gw.active()) return "[2D_VEGETATION] needs a [2D_AQUIFER] column under the cells (the root zone is the aquifer's unsaturated zone).";
    n_ = mesh.n_cells();
    const auto n = static_cast<std::size_t>(n_);
    enabled_.assign(n, 0);
    for (auto* v : {&area_, &laimean_, &kc_, &sl_m_, &nfc_, &fc_, &wp_, &ks0_, &ks1_, &ks2_,
                    &canopy_m_, &smax_m_, &scf_, &ep_, &demand_, &stress_, &wet_})
        v->assign(n, 0.0);
    lai_idx_.assign(n, -1); kc_idx_.assign(n, -1);
    const double rain_depth = uf.rate * 3600.0;        // authored mm (SI) or in (US) → m
    for (int pass = 0; pass < 3; ++pass) {
        for (const auto& r : cfg.vegetation) {
            if (r.scope != pass) continue;
            const int lai_idx = ctx.patterns.find(r.lai_pattern);
            const int kc_idx = ctx.patterns.find(r.kc_pattern);
            if (lai_idx < 0) return "[2D_VEGETATION] LAI pattern '" + r.lai_pattern + "' is not a [PATTERNS] entry";
            if (kc_idx < 0) return "[2D_VEGETATION] Kc pattern '" + r.kc_pattern + "' is not a [PATTERNS] entry";
            for (int i = 0; i < n_; ++i) {
                if (pass == 1) {
                    if (i >= static_cast<int>(mesh.tri_tag.size()) || mesh.tri_tag[static_cast<std::size_t>(i)] != r.tag) continue;
                } else if (pass == 2) {
                    if (r.cell != i) continue;
                }
                const auto u = static_cast<std::size_t>(i);
                enabled_[u] = 1;
                laimean_[u] = r.laimean; kc_[u] = r.kc; sl_m_[u] = r.sl * rain_depth; nfc_[u] = r.nfc;
                ks0_[u] = r.ks[0]; ks1_[u] = r.ks[1]; ks2_[u] = r.ks[2];
                lai_idx_[u] = lai_idx; kc_idx_[u] = kc_idx;
                const double theta_r = gw.state().theta_r[u], theta_s = gw.state().theta_s[u];
                wp_[u] = r.wp >= 0.0 ? r.wp : theta_r;
                fc_[u] = r.fc >= 0.0 ? r.fc : gw.drainageFieldCapacity(i);
                if (!(fc_[u] > wp_[u])) {
                    fc_[u] = wp_[u] + 0.05 * std::max(theta_s - wp_[u], 0.0);
                    warnings.push_back("[2D_VEGETATION] cell " + std::to_string(i + 1) +
                                       ": the column's drainage field capacity is at or below the wilting point; "
                                       "the stress law uses FC = WP + 5 % of the free pore space. Author FC explicitly or use a steeper soil law.");
                }
            }
        }
    }
    for (int i = 0; i < n_; ++i) area_[static_cast<std::size_t>(i)] = mesh.tri_area[static_cast<std::size_t>(i)];
    active_ = std::any_of(enabled_.begin(), enabled_.end(), [](int8_t e) { return e != 0; });
    intercepted_m3_ = canopy_evap_m3_ = 0.0;
    return {};
}

void MeshVegetation::refresh(const SimulationContext& ctx, const SubsurfaceSolver& gw, SurfaceStateData& state, bool evap_on) {
    if (!active_) return;
    const int month0 = datetime::monthOfYear(ctx.current_date) - 1;
    const int doy = datetime::dayOfYear(ctx.current_date);
    const double e_pot = std::max(0.0, ctx.climate_state.evap_rate) * 0.3048;   // internal ft/s → m/s
    const double root = gw.rootDepth();
    for (int i = 0; i < n_; ++i) {
        const auto u = static_cast<std::size_t>(i);
        if (!enabled_[u]) continue;
        const double lai = laimean_[u] * veglaw::monthlyPatternFactor(ctx.patterns, lai_idx_[u], month0);
        const double kc = kc_[u] * veglaw::monthlyPatternFactor(ctx.patterns, kc_idx_[u], month0);
        ep_[u] = e_pot * veglaw::seasonalKs({ks0_[u], ks1_[u], ks2_[u]}, doy) * kc;
        veglaw::canopyFactors(lai, sl_m_[u], smax_m_[u], scf_[u]);
        if (smax_m_[u] <= 0.0 && canopy_m_[u] > 0.0) {
            // leafless: the store goes to the surface as rain over the next step
            canopy_m_[u] = 0.0;      // (small; booked nowhere — the canopy is outside the 2D balance)
        }
        const double theta = gw.rootZoneTheta(i, root);
        stress_[u] = veglaw::es2Stress(theta, fc_[u], wp_[u], nfc_[u]);
        wet_[u] = veglaw::canopyWet(canopy_m_[u], smax_m_[u]);
        demand_[u] = veglaw::demandRate(ep_[u], scf_[u], stress_[u], wet_[u]);
        if (evap_on && state.evap_forced[u] == 0) state.evap_rate[u] = demand_[u];
    }
}

void MeshVegetation::beforeAdvance(SurfaceStateData& state, const std::vector<double>& rain_gross, double dt) {
    if (!active_ || rain_gross.size() != static_cast<std::size_t>(n_)) return;
    for (int i = 0; i < n_; ++i) {
        const auto u = static_cast<std::size_t>(i);
        if (!enabled_[u]) continue;
        const double gross = std::max(0.0, rain_gross[u]);
        const double inter = veglaw::bradenInterception(gross, dt, smax_m_[u], scf_[u], canopy_m_[u]);
        canopy_m_[u] += inter;
        intercepted_m3_ += inter * area_[u];
        state.rainfall[u] = gross - inter / dt;
    }
}

void MeshVegetation::afterAdvance(double dt) {
    if (!active_) return;
    for (int i = 0; i < n_; ++i) {
        const auto u = static_cast<std::size_t>(i);
        if (!enabled_[u] || canopy_m_[u] <= 0.0 || smax_m_[u] <= 0.0) continue;
        const double wet = veglaw::canopyWet(canopy_m_[u], smax_m_[u]);
        const double ei = std::min(canopy_m_[u], ep_[u] * dt * wet * scf_[u]);
        canopy_m_[u] -= ei;
        canopy_evap_m3_ += ei * area_[u];
    }
}

double MeshVegetation::canopyStorageVolume() const noexcept {
    double v = 0.0;
    for (int i = 0; i < n_; ++i) {
        const auto u = static_cast<std::size_t>(i);
        if (enabled_[u]) v += canopy_m_[u] * area_[u];
    }
    return v;
}

}  // namespace openswmm::twoD
