"""."""
import glob
import json
import os
import shutil
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


def read_list_file(input_file):
    """Read the output of generate_list_from_DB."""
    dams_dict = {}
    with open(input_file, "r", encoding="utf-8") as file_name:
        dams_to_process = file_name.readlines()
        for dam in dams_to_process:
            dams_dict[dam.split(",")[1].rstrip()] = dam.split(",")[0]
        return dams_dict


def generate_folders(output_path):
    """."""
    create_folder(output_path)

    camp_path = os.path.join(output_path, "camp")
    create_folder(camp_path)

    extract_path = os.path.join(output_path, "extracts")
    create_folder(extract_path)

    download_path = os.path.join(output_path, "download")
    create_folder(download_path)

    log_path = os.path.join(output_path, "log")
    create_folder(log_path)


def copy_customs_files_to_camp_folder(custom_path, dest_path, dam):
    """Copy custom files if exists and rename it to the correct nomenclature."""
    out_dam = None
    out_cut = None
    out_dat = None
    out_json = None
    if custom_path is not None:
        dam_info = glob.glob(os.path.join(custom_path, f"{dam}_daminfo*.*json"))
        if dam_info:
            out_dam = os.path.join(dest_path, f"{dam}_daminfo_custom.json")
            shutil.copy(dam_info[0], out_dam)
        cutline = glob.glob(os.path.join(custom_path, f"{dam}_*cutline*.*json"))
        if cutline:
            out_cut = os.path.join(dest_path, f"{dam}_cutline_custom.geojson")
            shutil.copy(cutline[0], out_cut)
        szi_dat = glob.glob(os.path.join(custom_path, f"{dam}_SZi*.dat"))
        if szi_dat:
            out_dat = os.path.join(dest_path, f"{dam}_SZi_custom.dat")
            shutil.copy(szi_dat[0], out_dat)
        json_cust = glob.glob(os.path.join(custom_path, f"params_{dam}*.json"))
        if json_cust:
            out_json = os.path.join(dest_path, f"params_{dam}_custom.json")
            shutil.copy(json_cust[0], out_json)
    return out_dam, out_cut, out_dat, out_json


def write_json(
    global_config_json,
    output_force_path=None,
    version_name=None,
    concat=False,
    ref_only=False,
    input_force_list=None,
):
    """."""
    generated_json = []
    with open(global_config_json, encoding="utf-8") as config_in:
        config = json.load(config_in)
        # Get functionnal parameters
        # Generate folder
        if output_force_path is not None:
            output_path = os.path.join(output_force_path, version_name)
        else:
            output_path = config["campaign"]["output_path"]
        generate_folders(output_path)
        if version_name:
            with open(
                os.path.join(output_path, "version.txt"), "w", encoding="utf-8"
            ) as version_file:
                version_file.write(version_name)
        database = config["campaign"]["database"]
        reference = config["campaign"]["reference"]
        # watermap = config["campaign"]["watermap"]
        dem = config["campaign"]["dem"]
        customs_files = config["campaign"]["customs_files"]
        dam_id_column = config["campaign"]["id_dam_column"]
        dam_name_column = config["campaign"]["dam_name_column"]
        mode = config["campaign"]["mode"]
        # Ensure output path exists
        output_list = os.path.join(output_path, "dam_list.txt")
        if input_force_list is not None:
            dict_all_dams = read_list_file(input_force_list)
        else:
            dict_all_dams = create_dam_list_from_db(
                database, dam_id_column, dam_name_column, output_list, concat
            )

        keys_ref = []
        if reference is not None:
            with open(reference, encoding="utf-8") as ref_file:
                ref_cont = json.load(ref_file)
                dict_ref = dict(ref_cont)
                keys_ref = dict_ref.keys()
                keys_ref = [int(key) for key in keys_ref]

        for dam, list_id_alt in dict_all_dams.items():
            id_dam = list_id_alt[0]
            maximum_alt = list_id_alt[1]
            if ref_only and keys_ref:
                if id_dam not in keys_ref:
                    print(f"{dam} {id_dam} not in reference")
                    print(keys_ref)
                    continue
            dict_dam = {}
            dam_path_name = dam.replace(" ", "-")
            output_dam_camp_path = os.path.join(output_path, "camp", dam_path_name)
            create_folder(output_dam_camp_path)
            output_dam_tmp = os.path.join(output_path, "camp", dam_path_name, "tmp")
            create_folder(output_dam_tmp)
            output_download_path = os.path.join(output_path, "download", dam_path_name)
            create_folder(output_download_path)
            output_dam_extract_path = os.path.join(
                output_path, "extracts", dam_path_name
            )
            create_folder(output_dam_extract_path)
            extract_dem = os.path.join(
                output_dam_extract_path, f"dem_extract_{dam_path_name}.tif"
            )
            if mode == "GDP":
                extract_wmap = os.path.join(output_dam_camp_path, "waterbody_bin.tif")
            else:
                extract_wmap = os.path.join(
                    output_dam_extract_path, f"wmap_extract_{dam_path_name}.tif"
                )
            extract_db = os.path.join(
                output_dam_extract_path, f"DB_{dam_path_name}.geojson"
            )

            elevsamp = config["cut_contourlines"]["elevsampling"]
            contourline_file = os.path.join(
                output_dam_extract_path, f"{dam_path_name}_contourlines{elevsamp}m.json"
            )
            if os.path.exists(contourline_file):
                cached_cont_file = contourline_file
            else:
                cached_cont_file = None
            # Handle customs files
            # If exists copy to camp folder
            (
                daminfo_file,
                cutline_file,
                szi_dat_file,
                json_file,
            ) = copy_customs_files_to_camp_folder(
                customs_files, output_dam_camp_path, dam_path_name
            )

            if cutline_file is None:
                cutline_file = os.path.join(
                    output_dam_camp_path, f"{dam_path_name}_cutline.geojson"
                )
            if szi_dat_file is None:
                szi_dat_file = os.path.join(
                    output_dam_camp_path, f"{dam_path_name}_SZi.dat"
                )

            dam_log = ensure_log_name(dam_path_name)
            dict_dam["chain"] = {
                "log_folder": os.path.join(output_path, "log"),
                "log_out": os.path.join(
                    output_path, "log", f"{dam_log}_{id_dam}_out.log"
                ),
                "log_err": os.path.join(
                    output_path, "log", f"{dam_log}_{id_dam}_err.log"
                ),
            }
            # dict_dam["area_mapping"] = {
            #     "dam_database": database,
            #     "dam_id": id_dam,
            #     "dam_name": dam,
            #     "id_db": dam_id_column,
            #     "watermap": watermap,
            #     "to_crop_file": dem,
            #     "out_dir": output_dam_extract_path,
            #     "out_wmap": extract_wmap,
            #     "out_dem": extract_dem,
            #     "output_download_path": output_download_path,
            #     "mode": mode,
            #     **config["area_mapping"],
            # }
            dict_dam["area_mapping"] = {
                "dam_name": dam,
                "dam_database": database,
                "out_dir": output_dam_extract_path,
                "dem": dem,
                "dam_id": id_dam,
                **config["area_mapping"],
            }
            if mode == "GDP":
                dict_dam["find_cutline_and_pdb"] = {
                    "database_file": extract_db,
                    "dem_raster": extract_dem,
                    "work_dir": output_dam_camp_path,
                    "id_db": id_dam,
                    "dam_name": dam,
                    "maximum_alt": maximum_alt,
                    **config["find_cutline_and_pdb"],
                }
            else:
                dict_dam["find_pdb_and_cutline"] = {
                    "infile": database,
                    "dam_id": id_dam,
                    "id_db": dam_id_column,
                    "watermap": extract_wmap,
                    "dem": extract_dem,
                    "out": output_dam_camp_path,
                    "info": daminfo_file,
                    "tmp": output_dam_tmp,
                    **config["find_pdb_and_cutline"],
                }

            # If no custom daminfo, it was generated by the previous function
            if daminfo_file is None:
                daminfo_file = os.path.join(
                    output_dam_camp_path, f"{dam_path_name}_daminfo.json"
                )
            dict_dam["cutline_score"] = {
                "infile": cutline_file,
                "watermap": extract_wmap,
                "out": output_dam_camp_path,
                **config["cutline_score"],
            }

            dict_dam["cut_contourlines"] = {
                "info": daminfo_file,
                "dem": extract_dem,
                "cutline": cutline_file,
                "level": cached_cont_file,
                "cache": output_dam_extract_path,
                "tmp": output_dam_tmp,
                "out": output_dam_camp_path,
                "mode": mode,
                **config["cut_contourlines"],
            }

            dict_dam["szi_to_model"] = {
                "szi_file": szi_dat_file,
                "database": database,
                "watermap": extract_wmap,
                "daminfo": daminfo_file,
                "outfile": os.path.join(
                    output_dam_camp_path, f"{dam_path_name}_model.png"
                ),
                # "custom_szi": szi_dat_file,
                **config["szi_to_model"],
            }

            if reference is not None:
                dict_dam["val_report"] = {
                    "infile": os.path.join(
                        output_dam_camp_path, f"{dam_path_name}_model.json"
                    ),
                    "outfile": os.path.join(
                        output_dam_camp_path, f"{dam_path_name}_report.png"
                    ),
                    "reffile": reference,
                    **config["val_report"],
                }
            json_out_file = os.path.join(
                output_dam_camp_path, f"params_{dam_path_name}.json"
            )
            with open(
                json_out_file,
                "w",
                encoding="utf-8",
            ) as write_file:
                json.dump(dict_dam, write_file, indent=4)

            if json_file is not None:
                generated_json.append(json_file)
            else:
                generated_json.append(json_out_file)
    return generated_json
