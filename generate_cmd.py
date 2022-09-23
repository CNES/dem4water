import geopandas as gpd
import argparse
import sys
import os
from datetime import datetime
import shutil
import glob


def save_previous_run(output_path, name):
    current_time = datetime.now()
    date_time = current_time.strftime("%Y%m%dT%Hh%Mm")
    save_path = os.path.join(output_path, "camp", name, f"{name}_{date_time}")
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    else:
        raise OSError(f"{save_path} already exists.""Conflict detected abort the processing.")
    file_of_interest = ["*_daminfo.json", "*_cutline.json", "*.png","*model.json","tmp/*.png"]
    for pred in file_of_interest:
        files = glob.glob(os.path.join(output_path, "camp", name, pred))
        for infile in files:
            shutil.copy(infile, save_path)
            print(infile, " copied to ", save_path)


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--infile", help="Input file")
    parser.add_argument("--dam_name", help="Dam NAME column", default="DAM_NAME")
    parser.add_argument("--id_db", help="Dam id field in database")
    parser.add_argument("-w", "--watermap", help="Input water map file")
    parser.add_argument("-d", "--dem", help="Input DEM")
    parser.add_argument("-r", "--radius", help="Extract radius (m)", default=None)
    parser.add_argument("-e", "--elevoffdam", help="Offset to elevation DAM (m)", default=50)
    
    parser.add_argument("-o", "--out", help="Output directory")
    parser.add_argument("-c", "--corrections_folder", help="Folder containing files to manually correct informations", default=None)
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")
    parser.add_argument("--exec_mode", action="store_true", help="Submit jobs to PBS")
    args = parser.parse_args(arguments)

    gdf_db = gpd.GeoDataFrame().from_file(args.infile)
    cmd_list=[]
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    logs_dir =  os.path.join(args.out, "logs_pbs")
    if not os.path.exists(logs_dir):
        os.mkdir(logs_dir)
    if args.radius is None:
        rad = ""
    else:
        rad = ",RADIUS={args.radius}"
    for name, id_dam in zip(gdf_db[args.dam_name], gdf_db[args.id_db]):
        param_dam = ""
        param_cut = ""
        param_vsurf = ""
        save = False
        if args.corrections_folder is not None:
            dam_info = glob.glob(os.path.join(args.corrections_folder, f"{name}_daminfo.json"))
            cutline = glob.glob(os.path.join(args.corrections_folder, f"{name}_cutline.json"))
            vsurface = glob.glob(os.path.join(args.corrections_folder, f"{name}_vSurfaces.json"))
            if dam_info:
                if len(dam_info) != 1:
                    raise ValueError(f"More than one file found for dam info for {name} dam.")
                param_dam = f",INFO_DAM={dam_info[0]}"
                save = True
            if cutline:
                if len(cutline) != 1:
                    raise ValueError(f"More than one file found for cutline for {name} dam.")
                param_cut = f",CUTLINE={cutline[0]}"
                save = True
            if vsurface:
                if len(vsurface) != 1:
                    raise ValueError(f"More than one file found for vsurface for {name} dam.")
                param_vsurf =f",VSURF={vsurface[0]}"
                save = True
            if save:
                print(f"Save previous results for {name}")
                save_previous_run(args.out, name)
                cmd = (f"qsub -v WD=$PWD,DAM={name},DAM_ID={id_dam},ID_FIELD={args.id_db},"
                       f"DB_PATH={args.infile},DEM_PATH={args.dem},WMAP_PATH={args.watermap},"
                       f"ROOT_DIR={args.out}{rad}{param_dam}{param_cut}{param_vsurf}"
                       f",ELEV_OFF_DAM={args.elevoffdam}"
                       f" -e {os.path.join(logs_dir, name+'_err.log')}"
                       f" -o {os.path.join(logs_dir, name+'_out.log')}"
                       " compute_hsv.pbs")
                cmd_list.append(cmd)
     
        else:
            cmd = (f"qsub -v WD=$PWD,DAM={name},DAM_ID={id_dam},ID_FIELD={args.id_db},"
                   f"DB_PATH={args.infile},DEM_PATH={args.dem},WMAP_PATH={args.watermap},"
                   f"ROOT_DIR={args.out}{rad}"
                   f",ELEV_OFF_DAM={args.elevoffdam}"
                   f" -e {os.path.join(logs_dir, name+'_err.log')}"
                   f" -o {os.path.join(logs_dir, name+'_out.log')}"
                   " compute_hsv.pbs")
            cmd_list.append(cmd)
        
    current_time = datetime.now()
    date_time = current_time.strftime("%Y%m%dT%Hh%Mm")
    with open(os.path.join(args.out, f"command_list_{date_time}.txt"),
              "w", encoding="utf-8") as out_file:
        out_file.write("\n".join(cmd_list))
    for cmd in cmd_list:
        if args.exec_mode:
            print(cmd)
            # barbare
            os.system(cmd)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
