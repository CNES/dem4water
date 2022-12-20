#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""Provide tools to generate json files."""
import argparse
import json
import os
import sys

from dem4water.area_mapping import area_mapping_args
from dem4water.cut_contourlines import cut_countourlines_ars
from dem4water.cutline_score import cutline_score_parameters
from dem4water.szi_from_contourline import szi_from_contourline_parameters
from dem4water.szi_to_model import szi_to_model_parameters


def get_all_parameters(output_path):
    """Create a json file from parameters."""
    all_parameters = {}
    # area mapping
    parser = area_mapping_args()
    area_map_args = parser.parse_args()
    all_parameters["area_mapping"] = {}

    for arg in vars(area_map_args):
        value = getattr(area_map_args, arg)
        if value is not None:
            all_parameters["area_mapping"][arg] = value

    # szi_from_contourline
    parser = szi_from_contourline_parameters()
    szi_cont_args = parser.parse_args()
    all_parameters["szi_from_contourlines"] = {}

    for arg in vars(szi_cont_args):
        value = getattr(szi_cont_args, arg)
        if value is not None:
            all_parameters["szi_from_contourlines"][arg] = value

    # cutline score
    parser = cutline_score_parameters()
    score_args = parser.parse_args()
    all_parameters["cutline_score"] = {}

    for arg in vars(score_args):
        value = getattr(score_args, arg)
        if value is not None:
            all_parameters["cutline_score"][arg] = value

    # cut_contourlines
    parser = cut_countourlines_ars()
    cut_args = parser.parse_args()
    all_parameters["cut_contourlines"] = {}

    for arg in vars(cut_args):
        value = getattr(cut_args, arg)
        if value is not None:
            all_parameters["cut_contourlines"][arg] = value

    # szi_to_model
    parser = szi_to_model_parameters()
    model_args = parser.parse_args()
    all_parameters["szi_to_model"] = {}

    for arg in vars(model_args):
        value = getattr(model_args, arg)
        if value is not None:
            all_parameters["szi_to_model"][arg] = value

    with open(
        os.path.join(output_path, "global.json"), "w", encoding="utf-8"
    ) as write_file:
        json.dump(all_parameters, write_file, indent=4)


def main():
    """Cli function to generate global config json."""
    parser_in = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser_in.add_argument(
        "-o", "--output_path", help="Global json config file", required=True
    )
    args = parser_in.parse_args()
    get_all_parameters(args.output_path)


if __name__ == "__main__":
    sys.exit(main())
