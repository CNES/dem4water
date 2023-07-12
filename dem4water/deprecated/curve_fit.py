import pickle

import matplotlib.pyplot as plt
import numpy as np

# from scipy.interpolate import interp1d
# from scipy.optimize import minimize
# from scipy.signal import savgol_filter


def nderiv(y_v, x_v):
    """Différence finie, dérivée de la fonction f."""
    nb_v = len(y_v)
    der = np.zeros(nb_v, "d")  # virgule flottante à double précision (double)
    # différences de part et d'autre
    # centrées sur les points intérieurs
    for i in range(1, nb_v - 1):
        der[i] = (y_v[i + 1] - y_v[i - 1]) / (x_v[i + 1] - x_v[i - 1])
    # différences sur un seul côté pour les extrémités
    der[0] = (y_v[1] - y_v[0]) / (x_v[1] - x_v[0])
    der[nb_v - 1] = (y_v[nb_v - 1] - y_v[nb_v - 2]) / (x_v[nb_v - 1] - x_v[nb_v - 2])
    return der


def compute_pdb_dem(dam):

    path = f"/home/btardy/Documents/activites/WATER/test_pdb/water_body/camp/{dam}"
    di = pickle.load(
        open(
            f"{path}/tmp/pdb_profile.pickle",
            "rb",
        )
    )

    found_pdb = False
    rad_l = di["rad_l"]
    der = di["der"]
    alt_l = di["alt_l"]
    alt_l = alt_l[::-1]
    rad_l = rad_l[::-1]
    der = der[::-1]

    fig, axs = plt.subplots(2)
    axs[0].plot(rad_l, alt_l, "o-r")
    axs[0].set(xlabel="Search Area to the Dam (m)", ylabel="Minimum Elevation")
    axs[0].label_outer()
    axs[1].plot(rad_l, abs(der), "o-b")
    axs[1].set(xlabel="Search Area to the Dam (m)", ylabel="d(Minimum Elevation)")
    axs[1].label_outer()

    for i_r, i_a, i_d in zip(rad_l, alt_l, der):
        #  print(i_r, i_a, i_d)
        if abs(i_d) > 0.15:
            #  and abs(i_d) < 0.2 :
            found_pdb = True
            rad_pdb = i_r
            alt_pdb = i_a
            fig.suptitle("PDB profile (1st pass detection)")
            axs[0].plot(rad_pdb, alt_pdb, "x")
            axs[1].plot(rad_pdb, abs(i_d), "x")
            break

    if found_pdb is True:
        print(f"@radius= {rad_pdb}m: local min = {alt_pdb}")
    else:
        for i_r, i_a, i_d in zip(rad_l, alt_l, der):
            #  print(i_r, i_a, i_d)
            if abs(i_d) > 0.10:
                found_pdb = True
                rad_pdb = i_r
                alt_pdb = i_a
                fig.suptitle("PDB profile (2nd pass detection)")
                axs[0].plot(rad_pdb, alt_pdb, "x")
                axs[1].plot(rad_pdb, abs(i_d), "x")
                print("PDB found during 2nd pass - It may not be reliable")
                print(f"@radius= {rad_pdb}m: local min = {alt_pdb}")
                break
            # Only if not if processed

    if found_pdb is True:
        print(f"@radius= {rad_pdb}m: local min = {alt_pdb}")
    else:
        for i_r, i_a, i_d in zip(rad_l, alt_l, der):
            #  print(i_r, i_a, i_d)
            if abs(i_d) > 0.02:
                found_pdb = True
                rad_pdb = i_r
                alt_pdb = i_a
                fig.suptitle("PDB profile (3rd pass detection)")
                axs[0].plot(rad_pdb, alt_pdb, "x")
                axs[1].plot(rad_pdb, abs(i_d), "x")
                print("PDB found during 2nd pass - It may not be reliable")
                print(f"@radius= {rad_pdb}m: local min = {alt_pdb}")
                break
            # Only if not if processed
            print(f"@radius= {rad_pdb}m: elev = {alt_pdb}m: delev = {i_d}")
    fig.savefig(f"{path}/../{dam}_old_pdb_profile.png")


def compute_pdb(dam):
    path = f"/home/btardy/Documents/activites/WATER/test_pdb/water_body/camp/{dam}"
    di = pickle.load(
        open(
            f"{path}/tmp/pdb_profile.pickle",
            "rb",
        )
    )

    rad = di["rad_l"]
    der = di["der"]
    alt = di["alt_l"]

    print(rad)
    # yhat = savgol_filter((rad, alt), 10, 3)
    # f2 = interp1d(rad, alt, kind="cubic")

    poly = np.polyfit(rad, alt, 5)
    poly_y = np.poly1d(poly)(rad)

    der2 = nderiv(poly_y, rad)
    fig, axs = plt.subplots(2)
    axs[0].plot(rad, alt, "r")
    # axs[0].plot(rad, f2(rad), "m")
    axs[0].plot(rad, poly_y, "b")
    for i_r, i_a, i_d in zip(rad, alt, der2):
        if i_d > -0.2:
            axs[0].plot(i_r, i_a, "x")
            break
    axs[1].plot(rad, der, "r")
    axs[1].plot(rad, der2, "b")
    axs[0].label_outer()
    axs[1].label_outer()
    for i_r, i_a, i_d in zip(rad, alt, der2):
        if i_d > -0.2:
            axs[1].plot(i_r, i_d, "x")
            break
    fig.savefig(f"{path}/../{dam}_new_pdb_profile.png")


def find_values_pdb_old(alt, rad, der):

    for i, (i_r, i_a, i_d) in enumerate(zip(rad, alt, der)):
        #  print(i_r, i_a, i_d)
        if abs(i_d) > 0.15:
            return i_r, i_a, i_d
    for i, (i_r, i_a, i_d) in enumerate(zip(rad, alt, der)):

        if abs(i_d) > 0.10:
            return i_r, i_a, i_d
    for i, (i_r, i_a, i_d) in enumerate(zip(rad, alt, der)):
        if abs(i_d) > 0.02:
            return i_r, i_a, i_d
    return 0, 0, 0


def find_values_pdb__pente(alt, rad, der):
    print(alt)
    print(rad)
    print(der)
    size = 0  # len(rad) - 1
    der = abs(der)
    for i, (i_r, i_a, i_d) in enumerate(zip(rad, alt, der)):
        #  print(i_r, i_a, i_d)
        if not i == size:
            if der[i - 1] < der[i]:
                if abs(i_d) > 0.15:
                    return i_r, i_a, i_d
    for i, (i_r, i_a, i_d) in enumerate(zip(rad, alt, der)):

        if not i == size:
            if der[i - 1] < der[i]:
                if abs(i_d) > 0.10:
                    return i_r, i_a, i_d
    for i, (i_r, i_a, i_d) in enumerate(zip(rad, alt, der)):

        if not i == size:
            if der[i - 1] < der[i]:
                if abs(i_d) > 0.02:
                    return i_r, i_a, i_d
    return 0, 0, 0


def find_values_pdb_new(alt, rad):
    poly = np.polyfit(rad, alt, 10)
    poly_y = np.poly1d(poly)(rad)

    der2 = nderiv(poly_y, rad)
    for i_r, i_a, i_d in zip(rad, alt, der2):
        if i_d > -0.04:
            print("sup 0.04")
            return i_r, i_a, i_d, poly_y, der2

    for i_r, i_a, i_d in zip(rad, alt, der2):

        if i_d > -0.2:
            return i_r, i_a, i_d, poly_y, der2

    return 0, 0, 0, poly_y, der2


def find_values_pdb_new2(alt, rad):
    poly = np.polyfit(rad, alt, 10)
    # res = minimize(np.poly1d(poly), 50)
    # print(res)
    poly_y = np.poly1d(poly)(rad)
    der2 = nderiv(poly_y, rad)

    for i_r, i_a, i_d in zip(rad, alt, der2):
        if i_d < -0.2:
            return i_r, i_a, i_d, poly_y, der2

    return 0, 0, 0, poly_y, der2


def draw_graphe(dam):
    path = f"/home/btardy/Documents/activites/WATER/test_pdb/water_body/camp/{dam}"
    di = pickle.load(
        open(
            f"{path}/tmp/pdb_profile.pickle",
            "rb",
        )
    )
    rad_o = di["rad_l"]
    der_o = di["der"]
    alt_o = di["alt_l"]
    ind = np.argwhere(np.array(rad_o) > 400)[0][0]

    rad = rad_o[:ind]
    der = der_o[:ind]
    alt = alt_o[:ind]
    alt_l = alt[::-1]
    rad_l = rad[::-1]
    der = der[::-1]

    r_o, a_o, d_o = find_values_pdb_old(alt_l, rad_l, der)
    r_n2, a_n2, d_n2, poly_y2, der22 = find_values_pdb_new2(alt_l, rad_l)
    r_o2, a_o2, d_o2 = find_values_pdb__pente(alt_l, rad_l, der22)
    r_n, a_n, d_n, poly_y, der2 = find_values_pdb_new(alt, rad)
    fig, axs = plt.subplots(2)
    axs[0].plot(rad_o, alt_o, ".-k", alpha=0.5)
    axs[0].plot(rad_l, alt_l, "r")
    axs[0].set(xlabel="Search Area to the Dam (m)", ylabel="Minimum Elevation")
    axs[0].label_outer()
    # axs[1].plot(rad_l, abs(der), "o-b")
    axs[1].plot(rad_o, (der_o), "k", alpha=0.5)
    axs[1].plot(rad_l, (der), "r")
    axs[1].plot(rad_l, (der22), "c")

    axs[0].plot(rad, poly_y, "m", alpha=0.5)
    axs[1].plot(rad, (der2), "m", alpha=0.5)
    axs[1].set(xlabel="Search Area to the Dam (m)", ylabel="d(Minimum Elevation)")
    axs[1].label_outer()

    if r_o2 > 0:
        axs[0].plot(r_o2, a_o2, "sc", label="origin pente")
        axs[1].plot(r_o2, (d_o2), "sc")
        print(f"{dam} old2 : {a_o2}")

    if r_o > 0:
        axs[0].plot(r_o, a_o, "xk", label="origin")
        axs[1].plot(r_o, (d_o), "xk")
        print(f"{dam} old : {a_o}")

    if r_n > 0:
        axs[0].plot(r_n, a_n, "om", label="lisse")
        axs[1].plot(r_n, (d_n), "om")
        print(f"{dam} new : {a_n}")
    if r_n2 > 0:
        axs[0].plot(r_n2, a_n2, "^b", label="lisse pente")
        axs[1].plot(r_n2, (d_n2), "^b")
        print(f"{dam} new2 : {a_n2}")
    fig.legend()
    fig.savefig(f"{path}/../{dam}_cmp_pdb_profile.png")


for dam in ["Alesani", "Aussoue", "Bissorte", "Bouvante", "Soulages"]:
    # compute_pdb(dam)
    # compute_pdb_dem(dam)
    draw_graphe(dam)
