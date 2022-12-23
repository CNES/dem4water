"""."""
import json
import os
import unicodedata

from dem4water.tools.generate_list_from_DB import create_dam_list_from_db


def create_folder(path_in):
    """."""
    if not os.path.exists(path_in):
        os.mkdir(path_in)


def ensure_log_name(name):
    """."""
    out_name = unicodedata.normalize("NFKD", name)
    out_name = out_name.encode("ascii", "ignore")
    out_name = out_name.decode("utf-8")
    return out_name


def genrerate_folders(output_path):
    """."""
    create_folder(output_path)

    camp_path = os.path.join(output_path, "camp")
    create_folder(camp_path)

    extract_path = os.path.join(output_path, "extracts")
    create_folder(extract_path)

    log_path = os.path.join(output_path, "log")
    create_folder(log_path)


def write_json(global_config_json):
    """."""
    with open(global_config_json, encoding="utf-8") as config_in:
        config = json.load(config_in)
        # Get functionnal parameters
        # Generate folder
        output_path = config["campaign"]["output_path"]
        genrerate_folders(output_path)

        database = config["campaign"]["database"]
        reference = config["campaign"]["reference"]
        watermap = config["campaign"]["watermap"]
        dem = config["campaign"]["dem"]
        customs_files = config["campaign"]["customs_files"]
        dam_id_column = config["campaign"]["id_dam_column"]
        dam_name_column = config["campaign"]["dam_name_column"]

        # Ensure output path exists
        output_list = os.path.join(output_path, "dam_list.txt")
        dict_all_dams = create_dam_list_from_db(
            database, dam_id_column, dam_name_column, output_list
        )
        for dam, id_dam in dict_all_dams.items():
            dict_dam = {}
            output_dam_camp_path = os.path.join(output_path, "camp", dam)
            create_folder(output_dam_camp_path)
            output_dam_tmp = os.path.join(output_path, "camp", dam, "tmp")
            create_folder(output_dam_tmp)
            output_dam_extract_path = os.path.join(output_path, "extracts", dam)
            create_folder(output_dam_extract_path)
            extract_dem = os.path.join(
                output_dam_extract_path, f"dem_extract_{dam}.tif"
            )

            extract_wmap = os.path.join(
                output_dam_extract_path, f"wmap_extract_{dam}.tif"
            )

            elevsamp = config["cut_contourlines"]["elevsampling"]
            contourline_file = os.path.join(
                output_dam_extract_path, f"{dam}_contourlines{elevsamp}m.json"
            )
            if os.path.exists(contourline_file):
                cached_cont_file = contourline_file
            else:
                cached_cont_file = None
            # TODO: handle customs
            if customs_files is not None:
                print(f"Looking for custom files for {dam}")
            cutline_file = ""
            daminfo_file = ""
            szi_dat_file = ""

            dict_dam["area_mapping"] = {
                "input_database": database,
                "dam_id": id_dam,
                "id_column_name": dam_id_column,
                "watermap": watermap,
                "dem": dem,
                "out_wmap": extract_wmap,
                "out_dem": extract_dem,
                **config["area_mapping"],
            }
            dict_dam["find_pdb_and_cutline"] = {
                "dam_id": id_dam,
                "id_db": dam_id_column,
                "infile": database,
                "out": output_dam_camp_path,
                "out_cutline": cutline_file,
                "out_daminfo": daminfo_file,
                "tmp": output_dam_tmp,
                **config["find_pdb_and_cutline"],
            }
            dict_dam["cutline_score"] = {
                "infile": cutline_file,
                "watermap": extract_wmap,
                "out": output_dam_camp_path,
                **config["cutline_score"],
            }

            dict_dam["cut_countourlines"] = {
                "info": daminfo_file,
                "dem": extract_dem,
                "cut": cutline_file,
                "level": cached_cont_file,
                "cache": output_dam_extract_path,
                "tmp": output_dam_tmp,
                "out": output_dam_tmp,
                **config["cut_contourlines"],
            }

            dict_dam["szi_to_model"] = {
                "szi_file": szi_dat_file,
                "database": database,
                "watermap": extract_wmap,
                "daminfo": daminfo_file,
                **config["szi_to_model"],
            }

            if reference is not None:
                dict_dam["val_report"] = {}
            with open(
                os.path.join(output_dam_camp_path, f"param_{dam}.json"),
                "w",
                encoding="utf-8",
            ) as write_file:
                json.dump(dict_dam, write_file, indent=4)


write_json("/home/btardy/Documents/activites/code/dem4water/global.json")
