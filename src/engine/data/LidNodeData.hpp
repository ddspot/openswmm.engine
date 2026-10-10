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
/// (``NAME EVAPFACTORS winter spring summer``). Replaces the storage node's
/// Fevap scaling of the climate rate when present.
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
