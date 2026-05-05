import MaizeWrapper


def build_plant():
    gen = MaizeWrapper.Maize()
    gen.reset()
    gen.set_species("maize")
    gen.set_global_settings(
        density_v=200.0,
        density_u=200.0,
        stem_density_scale=0.45,
        stem_row_override=0,
        stem_ribbon_spacing=0.3,
    )

    tiller = MaizeWrapper.Tiller(
        type="main",
        radius=0.02,
        alpha_deg=4.0,
        stem_shrink=0.001,
        azimuth_deg=180.0,
        azimuth_noise=20.0,
        random_seed=1337,
    ).to_ctype()
    tiller_index = gen.add_tiller(tiller)

    leaves = [
        MaizeWrapper.Leaf(id=0, distance=0.15, leaf_length=0.35, leaf_width=0.06, leaf_angle=32.0, droopiness=-25.0),
        MaizeWrapper.Leaf(id=1, distance=0.11, leaf_length=0.40, leaf_width=0.075, leaf_angle=38.0, droopiness=-15.0),
        MaizeWrapper.Leaf(id=2, distance=0.18, leaf_length=0.50, leaf_width=0.09, leaf_angle=45.0, droopiness=-10.0),
        MaizeWrapper.Leaf(id=3, distance=0.27, leaf_length=0.42, leaf_width=0.10, leaf_angle=52.0, droopiness=-0.0),
    ]
    for leaf in leaves:
        gen.add_leaf(tiller_index, leaf.to_ctype())

    gen.rebuild()
    gen.save_obj("plants/procedural_from_code")


def main():
    MaizeWrapper.build_example()
    print("Maize model built successfully.")

    build_plant()
    print("Procedural code-built plant saved to plants/procedural_from_code.obj")

    MaizeWrapper.from_xml_to_obj(
        "plants/plant_0.xml",
        "plants/myplant",
        leaf_texture_path="plants/maize_leaf.png",
        stem_texture_path="plants/maize_stem_texture.png",
    )
    print("XML-loaded plant saved to plants/myplant.obj")

    MaizeWrapper.from_xml_to_obj(
        "plants/plant_0.xml",
        "plants/separated",
        leaf_texture_path="plants/maize_leaf.png",
        stem_texture_path="plants/maize_stem_texture.png",
        include_stem=True,
        include_leaves=True,
        separate_leaves=True,
    )
    print("Separated parts saved to plants/separated/")


if __name__ == "__main__":
    main()
