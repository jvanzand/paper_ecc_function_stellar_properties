"""Named stellar subsamples used by the eccentricity occurrence calculations."""


SAMPLE_DEFINITIONS = {
    "allstars": (None, "All Stars"),
    "highMstar": ("Mstar > 1", r"$M_{\star} > 1 M_{\odot}$"),
    "lowMstar": ("Mstar <= 1", r"$M_{\star} \leq 1 M_{\odot}$"),
    "highFeH": ("feh > 0", "[Fe/H] > 0"),
    "lowFeH": ("feh <= 0", r"[Fe/H] $\leq$ 0"),
    "highAct": (
        "logrhk > -4.969 & 0.82 <= Mstar <= 1.21",
        # log R'_HK = -4.969 corresponds to 5.0 Gyr (Mamajek & Hillenbrand 2008).
        "Age < 5.0 Gyr",
    ),
    "lowAct": (
        "logrhk <= -4.969 & 0.82 <= Mstar <= 1.21",
        r"Age $\geq$ 5.0 Gyr",
    ),
}


def _combined_samples():
    definitions = {}
    mass = {"highMstar": "Mstar > 1", "lowMstar": "Mstar <= 1"}
    metallicity = {"highFeH": "feh > 0", "lowFeH": "feh <= 0"}
    activity = {
        "highAct": "logrhk > -4.969 & 0.82 <= Mstar <= 1.21",
        "lowAct": "logrhk <= -4.969 & 0.82 <= Mstar <= 1.21",
    }

    for mass_name, mass_query in mass.items():
        for feh_name, feh_query in metallicity.items():
            name = mass_name + feh_name
            definitions[name] = (
                "({}) & ({})".format(mass_query, feh_query), name,
            )
            for activity_name, activity_query in activity.items():
                three_parameter_name = name + activity_name
                definitions[three_parameter_name] = (
                    "({}) & ({}) & ({})".format(
                        mass_query, feh_query, activity_query
                    ),
                    three_parameter_name,
                )
    return definitions


SAMPLE_DEFINITIONS.update(_combined_samples())

# Public shape expected by occurrence.run.run_multiple().
tier2_df_cuts_dict = {
    name: [{"star_df_query": query}, title]
    for name, (query, title) in SAMPLE_DEFINITIONS.items()
}
