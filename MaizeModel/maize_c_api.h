#pragma once

#ifdef _WIN32
#if defined(MAIZE_C_API_BUILD)
#define MAIZE_C_API __declspec(dllexport)
#else
#define MAIZE_C_API __declspec(dllimport)
#endif
#else
#define MAIZE_C_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct MaizeHandle MaizeHandle;

typedef struct {
    const char* type;
    float radius;
    float alpha_deg;
    float stem_shrink;
    float azimuth_deg;
    float azimuth_noise;
    unsigned int random_seed;
} MaizeTillerDesc;

typedef struct {
    int id;
    float distance;
    float leaf_length;
    float leaf_width;
    float azimuth_deg;
    float leaf_angle;
    float droopiness;
    float stem_inclination_deg;
    int spline_points;
    float width_taper;
    float leaf_twist;
    float leaf_curl;
    float wave_l_amp;
    float wave_l_freq;
    float wave_l_phase;
    float wave_r_amp;
    float wave_r_freq;
    float wave_r_phase;
    float surface_noise_amp;
    float surface_noise_freq;
    float midrib_tip_taper_start;
    float midrib_texture_strength;
    float midrib_width;
    float ligule_wrap_length_scale;
    float ligule_unfold_sharpness;
    float sheath_outer_scale;
    int use_ctrl_overrides;
} MaizeLeafDesc;

typedef struct {
    float connection_x;
    float connection_y;
    float connection_z;
    float tip_x;
    float tip_y;
    float tip_z;
    float base_tangent_x;
    float base_tangent_y;
    float base_tangent_z;
    float leaf_length;
    float azimuth_deg;
    float inclination_deg;
} MaizeLeafSplineTraits;

MAIZE_C_API MaizeHandle* maize_create(void);
MAIZE_C_API void maize_destroy(MaizeHandle* handle);

MAIZE_C_API int maize_set_species(MaizeHandle* handle, const char* species);
MAIZE_C_API void maize_reset_plant(MaizeHandle* handle);

MAIZE_C_API int maize_set_global_settings(
    MaizeHandle* handle,
    float density_v,
    float density_u,
    float stem_density_scale,
    int stem_row_override,
    float stem_ribbon_spacing);

MAIZE_C_API int maize_add_tiller(MaizeHandle* handle, const MaizeTillerDesc* desc);
MAIZE_C_API int maize_set_tiller(MaizeHandle* handle, int tiller_index, const MaizeTillerDesc* desc);
MAIZE_C_API int maize_tiller_count(MaizeHandle* handle);

MAIZE_C_API int maize_add_leaf(MaizeHandle* handle, int tiller_index, const MaizeLeafDesc* desc);
MAIZE_C_API int maize_set_leaf(MaizeHandle* handle, int tiller_index, int leaf_index, const MaizeLeafDesc* desc);
MAIZE_C_API int maize_leaf_count(MaizeHandle* handle, int tiller_index);

MAIZE_C_API int maize_rebuild_geometry(MaizeHandle* handle);
MAIZE_C_API int maize_get_triangle_count(MaizeHandle* handle);
MAIZE_C_API int maize_leaf_spline_trait_count(MaizeHandle* handle);
MAIZE_C_API int maize_get_leaf_spline_traits(
    MaizeHandle* handle,
    int leaf_index,
    MaizeLeafSplineTraits* out_traits);

MAIZE_C_API int maize_save_obj(
    MaizeHandle* handle,
    const char* path,
    const char* leaf_texture_path,
    const char* stem_texture_path,
    int include_stem,
    int include_leaves,
    int separate_leaves);

MAIZE_C_API int maize_load_xml(MaizeHandle* handle, const char* path);
MAIZE_C_API int maize_save_xml(MaizeHandle* handle, const char* path);

MAIZE_C_API void maize_default_tiller_desc(MaizeTillerDesc* out_desc);
MAIZE_C_API void maize_default_leaf_desc(MaizeLeafDesc* out_desc);

/* ABI version introspection — added in C-API 1.1.0 */
#define MAIZE_C_API_VERSION_MAJOR 1
#define MAIZE_C_API_VERSION_MINOR 1
#define MAIZE_C_API_VERSION_PATCH 0

MAIZE_C_API int maize_c_api_version_major(void);
MAIZE_C_API int maize_c_api_version_minor(void);
MAIZE_C_API int maize_c_api_version_patch(void);
MAIZE_C_API int maize_leaf_desc_size(void);
MAIZE_C_API int maize_tiller_desc_size(void);
MAIZE_C_API int maize_leaf_spline_traits_size(void);

#ifdef __cplusplus
}
#endif
