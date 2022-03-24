#!/usr/bin/env python3

from math import atan2, cos, pi, sin, sqrt


def distance(lat1, lng1, lat2, lng2):
    # return distance as meter (if wanted km distance, remove "* 1000")
    radius = 6371 * 1000

    dLat = (lat2 - lat1) * pi / 180
    dLng = (lng2 - lng1) * pi / 180

    lat1 = lat1 * pi / 180
    lat2 = lat2 * pi / 180

    val = sin(dLat / 2) * sin(dLat / 2) + sin(dLng / 2) * sin(dLng / 2) * cos(
        lat1
    ) * cos(lat2)
    ang = 2 * atan2(sqrt(val), sqrt(1 - val))

    return radius * ang
