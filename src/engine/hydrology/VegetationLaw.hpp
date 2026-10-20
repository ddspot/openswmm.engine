// SPDX-License-Identifier: Apache-2.0
// UrbanEVA vegetation law — the pieces shared by the storage-node LID
// (hydrology/LidNode.cpp, Stage 1 of the port) and the vegetated mesh cell
// (2d/subsurface/MeshVegetation.cpp). All depths in one length unit, rates
// per second; the Braden law is written per hour as in lidproc.c.
#ifndef OPENSWMM_HYDROLOGY_VEGETATION_LAW_HPP
#define OPENSWMM_HYDROLOGY_VEGETATION_LAW_HPP
#include <algorithm>
#include <array>
#include <cmath>
#include "../data/InflowData.hpp"

namespace openswmm::veglaw {

/// UrbanEVA's seasonal Ks: anchored at day 79 / 172 / 263 / 355 (winter,
/// spring, summer, winter) and interpolated linearly.
inline double seasonalKs(const std::array<double, 3>& ks, int doy) {
    static const int days[4] = {79, 172, 263, 355};
    const double v[4] = {ks[0], ks[1], ks[2], ks[0]};
    if (doy <= days[0] || doy >= days[3]) {
        const double d0 = days[3] - 365.0, d1 = days[0];
        const double x = doy <= days[0] ? doy : doy - 365.0;
        return v[3] + (v[0] - v[3]) * (x - d0) / (d1 - d0);
    }
    for (int k = 0; k < 3; ++k)
        if (doy >= days[k] && doy <= days[k + 1])
            return v[k] + (v[k + 1] - v[k]) * (doy - days[k]) / double(days[k + 1] - days[k]);
    return ks[0];
}

/// A MONTHLY [PATTERNS] factor (1 when the pattern is missing or not MONTHLY).
inline double monthlyPatternFactor(const PatternData& patterns, int idx, int month0) {
    if (idx < 0 || idx >= patterns.count()) return 1.0;
    const auto& f = patterns.factors[static_cast<std::size_t>(idx)];
    if (patterns.types[static_cast<std::size_t>(idx)] != 0 || f.empty()) return 1.0;
    return f[static_cast<std::size_t>(month0 % static_cast<int>(f.size()))];
}

/// Braden canopy: Smax = SL · LAI per unit canopy area, SCF = 1 − 0.7^LAI.
inline void canopyFactors(double lai, double sl, double& smax, double& scf) {
    smax = std::max(0.0, sl * lai);
    scf = lai > 0.0 ? 1.0 - std::pow(0.7, lai) : 0.0;
}

/// Wet-canopy fraction f = (S/Smax)^(2/3) on the whole-area store (lidproc.c).
inline double canopyWet(double store, double smax) {
    return (smax > 0.0 && store > 0.0) ? std::pow(std::min(1.0, store / smax), 2.0 / 3.0) : 0.0;
}

/// Braden (I-3) interception over `dt` from a rain rate (depth/s): the canopy
/// takes I = Smax · (1 − 1 / (1 + SCF · P / Smax)) per hour (I ≤ SCF · P) and
/// its whole-area store fills to Smax · SCF. Returns the depth intercepted.
inline double bradenInterception(double rain_rate, double dt, double smax, double scf, double store) {
    if (!(rain_rate > 0.0) || !(smax > 0.0) || !(scf > 0.0) || !(dt > 0.0)) return 0.0;
    const double p_h = rain_rate * 3600.0;
    const double i_h = smax * (1.0 - 1.0 / (1.0 + scf * p_h / smax));
    const double headroom = std::max(0.0, smax * scf - store);
    return std::min({rain_rate * dt, i_h / 3600.0 * dt, headroom});
}

/// Es-2 soil-moisture stress sqrt(min(w / nFC, 1)), w the relative plant-
/// available water (θ − WP) / (FC − WP).
inline double es2Stress(double theta, double fc, double wp, double nfc) {
    if (!(fc > wp) || !(nfc > 0.0)) return 1.0;
    const double w = (theta - wp) / (fc - wp);
    return std::sqrt(std::clamp(w / nfc, 0.0, 1.0));
}

/// The demand the soil sees: Es (bare share, stressed) + Et (canopy share
/// less the wet fraction that evaporates as interception, Et-1).
inline double demandRate(double ep, double scf, double g_soil, double wet) {
    return std::max(0.0, ep) * ((1.0 - scf) * g_soil + scf * (1.0 - std::clamp(wet, 0.0, 1.0)));
}

}  // namespace openswmm::veglaw
#endif
