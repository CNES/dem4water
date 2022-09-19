import geopandas as gpd
import argparse
import sys
import os
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
        
        
        cmd = (f"qsub -v WD=$PWD,DAM={name},DAM_ID={id_dam},ID_FIELD={args.id_db},"
               f"DB_PATH={args.infile},DEM_PATH={args.dem},WMAP_PATH={args.watermap},"
               f"ROOT_DIR={args.out}{rad},ELEV_OFF_DAM={args.elevoffdam}"
               f" -e {os.path.join(logs_dir, name+'_err.log')}"
               f" -o {os.path.join(logs_dir, name+'_out.log')}"
               " compute_hsv.pbs")
        cmd_list.append(cmd)
        if args.exec_mode:
            print(cmd)
            # barbare
            os.system(cmd)
            
    with open(os.path.join(args.out, "command_list.txt"), "w", encoding="utf-8") as out_file:
        out_file.write("\n".join(cmd_list))
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
