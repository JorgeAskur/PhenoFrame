#include "maize_c_api.h"

#include "Maize.h"
#include "Descriptor.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <new>

struct MaizeHandle {
    MaizeModel model;
};

namespace {
    constexpr float kPi = 3.14159265358979323846f;

    bool isValid(MaizeHandle* handle) {
        return handle != nullptr;
    }

    bool hasTiller(MaizeHandle* handle, int index) {
        return handle && index >= 0 && index < (int)handle->model.plant().tillers.size();
    }

    bool hasLeaf(MaizeHandle* handle, int tiller_index, int leaf_index) {
        if (!handle || tiller_index < 0 || tiller_index >= (int)handle->model.plant().tillers.size()) return false;
        return leaf_index >= 0 && leaf_index < (int)handle->model.plant().tillers[tiller_index].leaves.size();
    }

    void apply_tiller(const MaizeTillerDesc& src, TillerDesc& dst) {
        dst.radius = src.radius;
        dst.alphaDeg = src.alpha_deg;
        dst.stemShrink = src.stem_shrink;
        dst.azimuthDeg = src.azimuth_deg;
        dst.azimuthNoise = src.azimuth_noise;
        dst.randomSeed = src.random_seed;
        if (src.type && *src.type) dst.type = src.type;
    }

    void apply_leaf(const MaizeLeafDesc& src, LeafDesc& dst) {
        dst.id = src.id;
        dst.distance = src.distance;
        dst.leafLength = src.leaf_length;
        dst.leafWidth = src.leaf_width;
        dst.azimuthDeg = src.azimuth_deg;
        dst.leafAngle = src.leaf_angle;
        dst.droopiness = src.droopiness;
        dst.stemInclinationDeg = src.stem_inclination_deg;
        dst.splinePoints = src.spline_points;
        dst.widthTaper = src.width_taper;
        dst.leafTwist = src.leaf_twist;
        dst.leafCurl = src.leaf_curl;
        dst.waveLAmp = src.wave_l_amp;
        dst.waveLFreq = src.wave_l_freq;
        dst.waveLPhase = src.wave_l_phase;
        dst.waveRAmp = src.wave_r_amp;
        dst.waveRFreq = src.wave_r_freq;
        dst.waveRPhase = src.wave_r_phase;
        dst.surfaceNoiseAmp = src.surface_noise_amp;
        dst.surfaceNoiseFreq = src.surface_noise_freq;
        dst.midribTipTaperStart = src.midrib_tip_taper_start;
        dst.midribTextureStrength = src.midrib_texture_strength;
        dst.midribWidth = src.midrib_width;
        dst.liguleWrapLengthScale = src.ligule_wrap_length_scale;
        dst.liguleUnfoldSharpness = src.ligule_unfold_sharpness;
        dst.sheathOuterScale = src.sheath_outer_scale;
        dst.useCtrlOverrides = (src.use_ctrl_overrides != 0);
    }

    float clamp01(float u) {
        return (u < 0.0f) ? 0.0f : ((u > 1.0f) ? 1.0f : u);
    }

    Vect3d sample_center_spline(const std::vector<Vect3d>& ctrl, float u) {
        if (ctrl.empty()) return Vect3d(0.0f, 0.0f, 0.0f);
        if (ctrl.size() == 1) return ctrl.front();
        return Maize::deBoor(clamp01(u), ctrl);
    }

    float spline_arc_length(const std::vector<Vect3d>& ctrl, int samples) {
        if (ctrl.size() < 2) return 0.0f;
        Vect3d prev = sample_center_spline(ctrl, 0.0f);
        float length = 0.0f;
        for (int i = 1; i <= samples; ++i) {
            const float u = (float)i / (float)samples;
            Vect3d p = sample_center_spline(ctrl, u);
            length += (p - prev).Length();
            prev = p;
        }
        return length;
    }

    int leaf_geom_count(const MaizeModel& model) {
        int count = 0;
        for (const auto& g : model.leafGeoms()) {
            if (g.type != LeafGeom::GEOM_SHEATH) ++count;
        }
        return count;
    }

    const LeafGeom* find_leaf_geom(const MaizeModel& model, int leaf_index) {
        if (leaf_index < 0) return nullptr;
        int idx = 0;
        for (const auto& g : model.leafGeoms()) {
            if (g.type == LeafGeom::GEOM_SHEATH) continue;
            if (idx == leaf_index) return &g;
            ++idx;
        }
        return nullptr;
    }
}

MaizeHandle* maize_create(void) {
    return new (std::nothrow) MaizeHandle();
}

void maize_destroy(MaizeHandle* handle) {
    if (!handle) return;
    delete handle;
}

int maize_set_species(MaizeHandle* handle, const char* species) {
    if (!isValid(handle) || !species) return 0;
    handle->model.plant().species = species;
    return 1;
}

void maize_reset_plant(MaizeHandle* handle) {
    if (!isValid(handle)) return;
    handle->model.plant() = PlantDesc{};
}

int maize_set_global_settings(
    MaizeHandle* handle,
    float density_v,
    float density_u,
    float stem_density_scale,
    int stem_row_override,
    float stem_ribbon_spacing) {
    if (!isValid(handle)) return 0;
    handle->model.m_densityV = density_v;
    handle->model.m_densityU = density_u;
    handle->model.m_stemDensityScale = stem_density_scale;
    return 1;
}

int maize_add_tiller(MaizeHandle* handle, const MaizeTillerDesc* desc) {
    if (!isValid(handle)) return -1;
    TillerDesc td;
    if (desc) apply_tiller(*desc, td);
    handle->model.plant().tillers.push_back(td);
    return (int)handle->model.plant().tillers.size() - 1;
}

int maize_set_tiller(MaizeHandle* handle, int tiller_index, const MaizeTillerDesc* desc) {
    if (!isValid(handle) || !desc) return 0;
    if (!hasTiller(handle, tiller_index)) return 0;
    apply_tiller(*desc, handle->model.plant().tillers[tiller_index]);
    return 1;
}

int maize_tiller_count(MaizeHandle* handle) {
    if (!isValid(handle)) return 0;
    return (int)handle->model.plant().tillers.size();
}

int maize_add_leaf(MaizeHandle* handle, int tiller_index, const MaizeLeafDesc* desc) {
    if (!isValid(handle) || !desc) return -1;
    if (!hasTiller(handle, tiller_index)) return -1;
    LeafDesc ld;
    apply_leaf(*desc, ld);
    handle->model.plant().tillers[tiller_index].leaves.push_back(ld);
    return (int)handle->model.plant().tillers[tiller_index].leaves.size() - 1;
}

int maize_set_leaf(MaizeHandle* handle, int tiller_index, int leaf_index, const MaizeLeafDesc* desc) {
    if (!isValid(handle) || !desc) return 0;
    if (!hasLeaf(handle, tiller_index, leaf_index)) return 0;
    apply_leaf(*desc, handle->model.plant().tillers[tiller_index].leaves[leaf_index]);
    return 1;
}

int maize_leaf_count(MaizeHandle* handle, int tiller_index) {
    if (!hasTiller(handle, tiller_index)) return 0;
    return (int)handle->model.plant().tillers[tiller_index].leaves.size();
}

int maize_rebuild_geometry(MaizeHandle* handle) {
    if (!isValid(handle)) return 0;
    handle->model.rebuildGeometry();
    return 1;
}

int maize_get_triangle_count(MaizeHandle* handle) {
    if (!isValid(handle)) return -1;
    return handle->model.getTriangleCount();
}

int maize_leaf_spline_trait_count(MaizeHandle* handle) {
    if (!isValid(handle)) return 0;
    handle->model.rebuildGeometry();
    return leaf_geom_count(handle->model);
}

int maize_get_leaf_spline_traits(
    MaizeHandle* handle,
    int leaf_index,
    MaizeLeafSplineTraits* out_traits) {
    if (!isValid(handle) || !out_traits || leaf_index < 0) return 0;

    handle->model.rebuildGeometry();
    const LeafGeom* leaf = find_leaf_geom(handle->model, leaf_index);
    if (!leaf || leaf->ctrlCenter.empty()) return 0;

    const auto& ctrl = leaf->ctrlCenter;
    const Vect3d connection = sample_center_spline(ctrl, 0.0f);
    const Vect3d tip = sample_center_spline(ctrl, 1.0f);
    const float du = (ctrl.size() > 1) ? (1.0f / (float)std::max(16, (int)ctrl.size() * 8)) : 1.0f;

    Vect3d tangent = sample_center_spline(ctrl, du) - connection;
    if (tangent.Length() < 1e-8f && ctrl.size() > 1) {
        tangent = ctrl[1] - ctrl[0];
    }

    const float tx = tangent.x();
    const float ty = tangent.y();
    const float tz = tangent.z();
    float azimuth = std::atan2(tz, tx) * (180.0f / kPi);
    if (azimuth < 0.0f) azimuth += 360.0f;
    const float horizontal = std::sqrt(tx * tx + tz * tz);
    const float inclination = std::atan2(ty, horizontal) * (180.0f / kPi);

    std::memset(out_traits, 0, sizeof(MaizeLeafSplineTraits));
    out_traits->connection_x = connection.x();
    out_traits->connection_y = connection.y();
    out_traits->connection_z = connection.z();
    out_traits->tip_x = tip.x();
    out_traits->tip_y = tip.y();
    out_traits->tip_z = tip.z();
    out_traits->base_tangent_x = tx;
    out_traits->base_tangent_y = ty;
    out_traits->base_tangent_z = tz;
    out_traits->leaf_length = spline_arc_length(ctrl, 128);
    out_traits->azimuth_deg = azimuth;
    out_traits->inclination_deg = inclination;
    return 1;
}

int maize_save_obj(
    MaizeHandle* handle,
    const char* path,
    const char* leaf_texture_path,
    const char* stem_texture_path,
    int include_stem,
    int include_leaves,
    int separate_leaves) {
    if (!isValid(handle) || !path) return 0;
    return handle->model.savePlantOBJ(
        path,
        leaf_texture_path,
        stem_texture_path,
        include_stem != 0,
        include_leaves != 0,
        separate_leaves != 0) ? 1 : 0;
}

int maize_load_xml(MaizeHandle* handle, const char* path) {
    if (!isValid(handle) || !path) return 0;

    PlantDesc global_backup = Maize::plant();
    bool ok = Descriptor::loadPlantXML(path);
    if (ok) {
        handle->model.plant() = Maize::plant();
    }
    Maize::plant() = global_backup;
    return ok ? 1 : 0;
}

int maize_save_xml(MaizeHandle* handle, const char* path) {
    if (!isValid(handle) || !path) return 0;

    PlantDesc global_backup = Maize::plant();
    Maize::plant() = handle->model.plant();
    bool ok = Descriptor::savePlantXML(path);
    Maize::plant() = global_backup;
    return ok ? 1 : 0;
}

void maize_default_tiller_desc(MaizeTillerDesc* out_desc) {
    if (!out_desc) return;
    std::memset(out_desc, 0, sizeof(MaizeTillerDesc));
    out_desc->type = "main";
    out_desc->radius = 0.02f;
    out_desc->alpha_deg = 0.0f;
    out_desc->stem_shrink = 0.001f;
    out_desc->azimuth_deg = 180.0f;
    out_desc->azimuth_noise = 45.0f;
    out_desc->random_seed = 1337;
}

void maize_default_leaf_desc(MaizeLeafDesc* out_desc) {
    if (!out_desc) return;
    std::memset(out_desc, 0, sizeof(MaizeLeafDesc));
    out_desc->id = 0;
    out_desc->distance = 0.0f;
    out_desc->leaf_length = 0.7f;
    out_desc->leaf_width = 0.08f;
    out_desc->azimuth_deg = 0.0f;
    out_desc->leaf_angle = 45.0f;
    out_desc->droopiness = 0.5f;
    out_desc->stem_inclination_deg = 0.0f;
    out_desc->spline_points = 4;
    out_desc->width_taper = 1.0f;
    out_desc->leaf_twist = 0.0f;
    out_desc->leaf_curl = 0.0f;
    out_desc->wave_l_amp = 0.0f;
    out_desc->wave_l_freq = 0.0f;
    out_desc->wave_l_phase = 0.0f;
    out_desc->wave_r_amp = 0.0f;
    out_desc->wave_r_freq = 0.0f;
    out_desc->wave_r_phase = 0.0f;
    out_desc->surface_noise_amp = 0.01f;
    out_desc->surface_noise_freq = 8.0f;
    out_desc->midrib_tip_taper_start = 0.75f;
    out_desc->midrib_texture_strength = 0.35f;
    out_desc->midrib_width = 0.075f;
    out_desc->ligule_wrap_length_scale = 6.8f;
    out_desc->ligule_unfold_sharpness = 2.2f;
    out_desc->sheath_outer_scale = 1.12f;
    out_desc->use_ctrl_overrides = 0;
}

int maize_c_api_version_major(void) { return MAIZE_C_API_VERSION_MAJOR; }
int maize_c_api_version_minor(void) { return MAIZE_C_API_VERSION_MINOR; }
int maize_c_api_version_patch(void) { return MAIZE_C_API_VERSION_PATCH; }
int maize_leaf_desc_size(void) { return (int)sizeof(MaizeLeafDesc); }
int maize_tiller_desc_size(void) { return (int)sizeof(MaizeTillerDesc); }
int maize_leaf_spline_traits_size(void) { return (int)sizeof(MaizeLeafSplineTraits); }
