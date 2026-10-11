// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Caleb Buahin
#ifndef OPENSWMM_LID_NODE_DATA_HPP
#define OPENSWMM_LID_NODE_DATA_HPP
#include <array>
#include <vector>
#include <string>
#include "../hydrology/RichardsColumn.hpp"
#include "../hydrology/Infiltration.hpp"

namespace openswmm {
// Authored, ordered layers. Thickness/suction are in rain-depth units (in/mm),
// conductivity in in/hr or mm/hr. BOTTOM is a boundary, not a numbered layer.
enum class LidNodeLayerKind { Surface = 0, Media = 1, Aggregate = 2, Bottom = 3 };
struct LidLayerTreatment {
    std::string pollutant;
    double removal = 0.0; // fraction of outgoing mass
    double decay = 0.0; // first-order rate, per day
    std::string expression; // optional R= or C=, applied after fixed removal
};
struct LidNodeLayer {
    LidNodeLayerKind kind = LidNodeLayerKind::Media;
    std::array<double, 7> params{};
    std::vector<LidLayerTreatment> treatment;
    richards::Material retention; // explicit inverse-metre parameters
    richards::Options flow; // control options, stored on first layer only
};
/// UrbanEVA vegetation on a storage-node LID (``NAME VEGETAT …`` row, the
/// SWMM-UrbanEVA grammar): crop factor and leaf-area index with MONTHLY
/// patterns, Braden interception (I-3) with wetted-fraction canopy
/// evaporation (Ei-2), unstressed transpiration (Et-1) and square-root
/// soil-moisture stress (Es-2) on the MEDIA water, seasonal Ks evapfactors
/// (``NAME EVAPFACTORS winter spring summer``). The storage node's Fevap stays
/// a multiplier on the climate rate.
struct LidNodeVegetation {
    bool enabled = false;
    double laimean = 0.0;          ///< mean leaf-area index (scaled by lai_pattern)
    double cropfactor = 1.0;       ///< Kc (scaled by kc_pattern)
    double sl = 0.0;               ///< interception storage per unit LAI, authored mm
    double nfc = 0.6;              ///< stress onset as a fraction of plant-available water
    std::string lai_pattern;       ///< [PATTERNS] MONTHLY name
    std::string kc_pattern;        ///< [PATTERNS] MONTHLY name
    int lai_pattern_idx = -1;      ///< resolved at initialize()
    int kc_pattern_idx = -1;
    std::array<double, 3> ks{0.78, 0.82, 0.85};  ///< seasonal evapfactors (winter, spring, summer)
    // UrbanEVA model selectors, kept for round trip; this build implements
    // I-3 / Ei-2 / Et-1 / Es-2 regardless of the authored codes.
    int soilevapmodel = 2, intercepmodel = 3, ievapmodel = 2, soiltransmodel = 1;
    int evap_while_ponding = 1, infilmodel = 1;
};

/// UrbanEVA irrigation on a storage-node LID (``NAME IRRIGATION DepletionFact
/// IrrigationCoeff Pattern SourceNode [TriggerHour]``): once a day at the
/// trigger hour, if the MEDIA water is below WP + d·(FC − WP) and the MONTHLY
/// pattern is > 0, a dose (c·FC − θ_thr)·thickness·area·pattern is drawn from
/// the source storage node and applied at the top of the column. Booked:
/// the source loses the volume the LID gains, so routing continuity closes.
struct LidNodeIrrigation {
    bool enabled = false;
    double depletion_fact = 0.0;   ///< trigger: fraction of plant-available water drawn down
    double irrigation_coeff = 0.0; ///< target: fraction of field capacity to refill to
    std::string pattern;           ///< [PATTERNS] MONTHLY name
    int pattern_idx = -1;
    std::string source;            ///< storage node supplying the water
    int source_idx = -1;
    int trigger_hour = 0;          ///< local hour the daily check runs (0-23)
};

/// Capillary rise from a `[2D_AQUIFER]` cell into a storage-node LID
/// (``NAME CAPILLARY ET_DEPTH Reach`` row). Stage 3a of the UrbanEVA port,
/// the SWMM 5 ``gwater.c`` lower-zone-ET mechanism carried over: after the
/// vegetation law has drawn what the MEDIA can give, the unmet potential ET
/// is supplied from the saturated store of the cell under the node's bed
/// when the water table lies within ``reach`` below the LID bottom,
/// ``rise = D · max(0, 1 − gap / reach)``. Booked on the aquifer side as
/// subsurface ET (the node's 1D books are untouched); the LID keeps the
/// cumulative volume for the vegetation summary.
/// Two laws share the row. ``ET_DEPTH`` (mode 0, Stage 3a) is the one above.
/// ``DARCY`` (mode 1, Stage 3b) is the process form: an unsaturated Darcy
/// flux between the aquifer soil halfway down the gap and the lowest MEDIA
/// cell, ``q = K̄ · ((ψ_m − ψ_a)/Δz − 1)`` upward when the media is the drier
/// side, bounded by the media's free pore space and the cell's drainable
/// water; the water enters the media (θ rises) through the ordinary node
/// exchange channel and ET then follows from the vegetation law. The media's
/// retention is Clapp–Hornberger on the authored suction head,
/// ``ψ(θ) = ψ_b · (φ/θ)^b`` (``b`` optional, default 4), its conductivity the
/// Clapp–Hornberger conductivity ``K = Ks · (θ/φ)^(2b+3)`` (the LID drainage
/// law ``exp(−HCO·(φ−θ))`` keeps far too much K in the dry range).
struct LidNodeCapillary {
    bool enabled = false;
    int mode = 0;                  ///< 0 = ET_DEPTH, 1 = DARCY
    double reach = 0.0;            ///< authored rain-depth units (mm / in) below the LID bottom
    double ch_b = 4.0;             ///< DARCY: Clapp–Hornberger retention exponent of the media
};

struct LidNodeCell {
    double bottom = 0.0, top = 0.0; // feet above invert
    double geometric_volume = 0.0, area = 0.0;
    double porosity = 1.0, theta = 0.0, field_capacity = 0.0, wilting_point = 0.0;
    double conductivity = 0.0, conductivity_slope = 0.0, suction = 0.0;
    LidNodeLayerKind kind = LidNodeLayerKind::Aggregate;
    int layer = 0;
};
struct LidWaterTransfer { int from, to; double volume; }; // -1 mobile, -2 lateral, -3 evaporation
struct LidPortTransfer { int link, cell; double volume; }; // positive into held water
struct LidNodeState {
    bool richards = false;
    std::vector<double> richards_water; // ft3, full porous stores; surface row is unused
    std::vector<double> richards_pressure; // ft, negative means suction
    richards::Report richards_report; // SI diagnostics for last accepted interval
    double richards_bottom_loss = 0.0; // ft3, outside-mesh system loss
    // One modified Green-Ampt history for the media directly below SURFACE.
    // Other layer/port wetting uses the physical cell and mobile inventories.
    GreenAmptState surface_infil;
    int infiltration_cell = -1;
    double infiltration_head = 0.0; // last reconciled mobile water-table depth

    std::vector<double> quality_mass; // [cell * pollutant + pollutant], retained stores
    std::vector<double> quality_old_water;
    std::vector<LidWaterTransfer> quality_transfers;
    std::vector<LidPortTransfer> quality_ports;
    std::vector<double> quality_outlet_conc; // [link * pollutant + pollutant], NaN means mobile
    double quality_old_mobile = 0.0;
    double mobile_delta = 0.0; // accepted port balance in the current hydraulic iteration
    std::vector<double> port_delta; // accepted link exchange volumes, reset each hydraulic iteration
    std::vector<LidNodeCell> cells; // top to bottom; MEDIA subdivided dynamically
    double held_volume = 0.0;      ///< cells + canopy
    double evap_volume = 0.0;      ///< soil/media + canopy evaporation this step
    double captured_flow = 0.0;
    double treated_volume = 0.0;
    double canopy_volume = 0.0;    ///< ft3, intercepted water on the vegetation
    double canopy_evap_total = 0.0;///< ft3, cumulative canopy (interception) evaporation
    double transp_total = 0.0;     ///< ft3, cumulative media evaporation under the vegetation law
    double irrigation_total = 0.0; ///< ft3, cumulative irrigation applied
    int irrigation_events = 0;     ///< number of doses
    long irrigation_last_day = -1; ///< whole-day index of the last daily check
    double caprise_total = 0.0;    ///< ft3, cumulative capillary rise from the aquifer (ET_DEPTH: as ET; DARCY: into the media)
    double caprise_step = 0.0;     ///< ft3, ET_DEPTH capillary ET this step (not in evap_volume: the aquifer books it)
};
struct LidNodeConfig {
    int control = -1;
    double initial_saturation = 0.0; // percent
};
struct LidNodeOutlet {
    int link = -1;
    int layer = 1; // one-based, top to bottom
    bool top = false;
};
} // namespace openswmm
#endif
