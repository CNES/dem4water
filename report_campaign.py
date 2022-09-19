import sys
import os
import argparse
import glob
import pandas as pd
import numpy as np
import json
def ensure_combo(folder):
    content = glob.glob(os.path.join(folder,"*combo.png"))
    if len(content) != 1:
        print(f"erreur barrage {folder}")


def check_pbd_delta(folder, name):
    
    szi_file = glob.glob(os.path.join(folder,"*SZi.dat"))
    dam_info = os.path.join(folder,f"{name}_daminfo.json")

        # input(elev)
    if szi_file:
        with open(dam_info, "r") as dam_inf:
            dam = json.load(dam_inf)
            # print(dam)
            # print(dam.keys())
            dam_info_dict = dam["features"][1]["properties"]
            elev = dam_info_dict["elev"]
            szi_file = szi_file[0]
            df = pd.read_csv(szi_file, header=None, names=["Zi","SZi"], sep=" ")
            # print(df.columns)
            # print(df.Zi)
            vals = list(df.Zi)
            delta = vals[-2] - vals[-1]
            delta2 = float(elev) - vals[-1]
            if delta > 150:
                print(name," ", delta)
    else:
        print(f"Erreur not dat file found in {folder}")
        delta = -1
        delta2 = -1
    return name, delta, delta2
    
def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-in", "--infolder", help="Folder containing outputs")

    args = parser.parse_args(arguments)

    list_dir = glob.glob(args.infolder+"/*")
    names = []
    deltas = []
    list_dir.sort()
    for folder in list_dir :
        dam_name = folder.split("/")[-1]
        ensure_combo(folder)
        ret =  check_pbd_delta(folder, dam_name)
        deltas.append(ret)
    
    deltas = sorted(deltas, key=lambda x: x[1])
    # print(deltas)
    for k, v, v2 in deltas:
        print(f"|{k}|{v}|{v2}|")
        
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
