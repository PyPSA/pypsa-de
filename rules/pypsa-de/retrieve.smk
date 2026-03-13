# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC BY 4.0

if config["pypsa-de"]["use_internal_db"]:
    if (ARIADNE_DATABASE_INTERNAL := dataset_version("ariadne_database_internal"))[
        "source"
    ] in ["primary"]:

        rule retrieve_ariadne_database_internal:
            params:
                source="internal",
            output:
                data="data/ariadne_database.csv",
            log:
                "logs/retrieve_ariadne_database_internal_primary.log",
            resources:
                mem_mb=1000,
            script:
                scripts("pypsa-de/retrieve_ariadne_database.py")

else:
    if (ARIADNE_DATABASE := dataset_version("ariadne_database"))["source"] in [
        "primary"
    ]:

        rule retrieve_ariadne_database:
            params:
                source="primary",
            output:
                data="data/ariadne_database.csv",
            log:
                "logs/retrieve_ariadne_database_primary.log",
            resources:
                mem_mb=1000,
            script:
                scripts("pypsa-de/retrieve_ariadne_database.py")

    if (ARIADNE_DATABASE := dataset_version("ariadne_database"))["source"] in [
        "archive"
    ]:

        rule retrieve_ariadne_database:
            params:
                source="archive",
                version=ARIADNE_DATABASE["version"],
            input:
                raw_xlsx=storage(ARIADNE_DATABASE["url"]),
            output:
                data="data/ariadne_database.csv",
            log:
                "logs/retrieve_ariadne_database_archive.log",
            resources:
                mem_mb=1000,
            script:
                scripts("pypsa-de/retrieve_ariadne_database.py")


if (ARIADNE_TEMPLATE := dataset_version("ariadne_template"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_ariadne_template:
        input:
            storage(ARIADNE_TEMPLATE["url"]),
        output:
            "data/template_ariadne_database.xlsx",
        run:
            move(input[0], output[0])


if (OPEN_MASTR := dataset_version("open_mastr"))["source"] in ["primary", "archive"]:

    rule retrieve_open_mastr:
        input:
            storage(OPEN_MASTR["url"]),
        params:
            "data/mastr",
        output:
            "data/mastr/bnetza_open_mastr_2023-08-08_B_biomass.csv",
            "data/mastr/bnetza_open_mastr_2023-08-08_B_combustion.csv",
        run:
            unpack_archive(input[0], params[0])


if (EGON := dataset_version("egon"))["source"] in ["build"]:

    rule retrieve_egon_data:
        params:
            url=EGON["url"],
            folder=EGON["folder"],
        output:
            spatial=f"{EGON['folder']}/demandregio_spatial_2018.json",
            mapping=f"{EGON['folder']}/mapping_technologies.json",
        shell:
            """
            mkdir -p {params.folder}
            curl -o {output.spatial} "{params.url}?id_spatial=5&year=2018"
            curl -o {output.mapping} "{params.url}_description?id_spatial=5"
            """


if (EGON := dataset_version("egon"))["source"] in ["archive"]:

    rule retrieve_egon_data:
        input:
            spatial=storage(f"{EGON['url']}/demandregio_spatial_2018.json"),
            mapping=storage(f"{EGON['url']}/mapping_technologies.json"),
        output:
            spatial=f"{EGON['folder']}/demandregio_spatial_2018.json",
            mapping=f"{EGON['folder']}/mapping_technologies.json",
        run:
            copy2(input["spatial"], output["spatial"])
            copy2(input["mapping"], output["mapping"])
