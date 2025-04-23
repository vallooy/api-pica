# -------------------------------------------------------------------------------
# Name:  algorithme calcul plus court chemin entre ceps à échantillonner
# Author:     © Philippe Vismara - Institut Agro - 2025
# version : 2025-03-20
# Ce logiciel est régi par la licence CeCILL-C soumise au droit français et respectant les principes de diffusion des logiciels libres.
# Vous pouvez utiliser, modifier et/ou redistribuer ce programme sous les conditions de la licence CeCILL-C telle que diffusée sur le site http://www.cecill.info
# -------

from geopy import distance
import pandas as pd
import fiona
from shapely.geometry import shape
import math

import cpmpy as cp

def read_geo_to_dataframe(path):
    records = []
    with fiona.open(path, 'r') as src:
        for feat in src:
            props = feat["properties"]
            geom = shape(feat["geometry"])
            props["geometry"] = geom
            # Ajouter la propriété 'Site' si elle existe dans le GeoJSON
            if "Site" in feat["properties"]:
                props["Site"] = feat["properties"]["Site"]
            else:
                raise KeyError("La propriété 'Site' est absente dans une ou plusieurs entités du fichier GeoJSON.")
            records.append(props)
    df = pd.DataFrame(records)
    return df


class Pica:
    def __init__(self, rows_file, geopandas_file):
        self.__distance_matrix = []
        # read file
        self.__geo_data = read_geo_to_dataframe(geopandas_file)
        print(self.__geo_data)
        print(self.__geo_data.columns)
        # definir l'index de la dataframe
        self.__geo_data.set_index("Site", inplace=True)

        self.__rows_data = pd.read_csv(rows_file)
        # on rajoute au début de data_ext un rang 0 dont les extremites ont le mêmes coordonnées que le premier inter-rang
        ligne_rang1 = self.__rows_data[self.__rows_data["rang"] == 1]
        ligne_rang0 = ligne_rang1.copy()
        ligne_rang0["rang"] = 0
        self.__rows_data = pd.concat([ligne_rang0, self.__rows_data], ignore_index=True)
        self.start_site = -1

    # -- constantes pour les directions (UP = côté du rang en direction du rang inférieur, DOWN = côté du rang en direction du rang supérieur)
    UP = 0
    DOWN = -1
    # -- precision pour les calculs de distance
    PRECISION = 1000

    def distance_matrix(self):
        return self.__distance_matrix

    def rows_data(self):
        return self.__rows_data

    def geo_data(self):
        return self.__geo_data

    def __compute_matrix(self):
        # on contruit une matrice de distance (2*nb_echantillons + 1) x (2*nb_echantillons + 1)
        # le premier site (indice 0) est le point de départ
        # les autres sites sont les sites à échantillonner (2 par cep) correspondant aux positiosn avant (down) et après (up) le rang du cep
        # la matrice contient les distances entre les points
        self.__distance_matrix = []

        # on ajoute les distances au point de départ du parcours (premier_site)
        ligne = [0]
        for cep in self.__sites:
            dist, chemin = self.trajet_vers(self.start_site, Pica.UP, cep, Pica.DOWN)
            ligne.append(math.trunc(Pica.PRECISION * dist))
            dist, chemin = self.trajet_vers(self.start_site, Pica.UP, cep, Pica.UP)
            ligne.append(math.trunc(Pica.PRECISION * dist))
        self.__distance_matrix.append(ligne)

        # on ajoute les distances entre les sites à échantillonner
        for i, cepi in enumerate(self.__sites):
            for cote_i in [Pica.DOWN, Pica.UP]:
                # on initialise la ligne de la matrice avec la distance au premier_point
                ligne = [self.__distance_matrix[0][i * 2 + 2 + cote_i]]
                for j, cepj in enumerate(self.__sites):
                    for cote_j in [Pica.DOWN, Pica.UP]:
                        dist, chemin = self.trajet_vers(cepi, cote_i, cepj, cote_j)
                        ligne.append(math.trunc(Pica.PRECISION * dist))
                self.__distance_matrix.append(ligne)

        # print(matrice_distances)

    # -- fonction pour calculer la longueur d'un chemin
    def __distances(self, listegps):
        dist = 0
        for i in range(1, len(listegps)):
            dist = dist + distance.geodesic(listegps[i - 1], listegps[i]).m
        return dist

    # --- accès aux coordonnées GPS d'un site (latitude, longitude) donc Y, X
    def coordonnes_gps_site(self, num_site):
        # return (self._geo_data["geometry"][num_site].y, self._geo_data["geometry"][num_site].x)
        return (
            self.__geo_data.loc[num_site, "geometry"].y,
            self.__geo_data.loc[num_site, "geometry"].x,
        )

    # -- acces au rang d'un site
    def site_row(self, num_site):
        return self.__geo_data.loc[num_site, "rang"]

    # -- calcul des points intremédiares (extrémitéés d'inter-rangs) entre deux points
    def trajet_vers(self, depart, cote_rang_depart, arrivee, cote_rang_arrivee):
        """calcule le meilleur chemin de depart a arrivee et retourne la distance
        ainsi que les coordonnées GPS du chemin (hors points de départ et d'arrivée)
        en passant éventuellement par les extrémités des deux inter-rangs ET par les extrémités intermédiaires.
        Si le trajet est direct, la liste des coordonnées sera vide sinon elle contiendra les coordonnées de 2 ou plus extrémités de rangs.
        Le trajet est direct si les deux points sont sur le même intrerrang-rang.

        Args:
        depart (int): Point de départ.
        cote_rang_depart (str): Côté du rang de départ (0 ou -1).
        arrivee (int): Point d'arrivée.
        cote_rang_depart (str): Côté du rang d'arrivée (0 ou -1).

        Returns:
            tuple: Une paire contenant la distance totale (float) et une liste des coordonnées GPS (list) du chemin.
        Exemple:
        Par exemple pour : trajet_vers(115, UP, 309, DOWN)
        on obtiendra : (108.27954733465921, [(43.54790492901712, 3.840181709048596), (43.54792220843044, 3.840194712848685), (43.54793948783887, 3.8402207204488614), (43.54796147980605, 3.840242393449003), (43.54797875920314, 3.8402553972490985), (43.547999180302455, 3.8402662337491686), (43.54775327653986, 3.840861608668911)])
        """
        if depart == arrivee and cote_rang_depart == cote_rang_arrivee:
            return 0, []
        distot = 0
        chemin = []
        # on doit optimiser le chemin entre l'échantillon pt -1 et pt
        ir_1 = self.__geo_data.loc[depart, "rang"] + cote_rang_depart
        ir_2 = self.__geo_data.loc[arrivee, "rang"] + cote_rang_arrivee
        # coord gps
        gps_depart = self.coordonnes_gps_site(depart)
        gps_arrivee = self.coordonnes_gps_site(arrivee)
        # on ne passe par les extrémités que si les deux points ne sont pas sur le même rang
        if ir_1 != ir_2:
            # coordonnées des extremités de rang à traverser
            liste_gps_gauche = []
            liste_gps_droite = []
            # attention il faut tester si ir_1 < ir_2 ou l'inverse
            # pour créer la liste des inter-rangs depuis ir1 jusqu'à ir2
            if ir_1 < ir_2:
                liste_ir = range(ir_1, ir_2 + 1)
            else:
                liste_ir = range(ir_1, ir_2 - 1, -1)
            # on parcours la liste des inter-rangs depuis ir1 jusqu'à ir2
            for ir in liste_ir:
                liste_gps_gauche.append(
                    (
                        self.__rows_data.loc[
                            self.__rows_data["rang"] == ir, "left_lat"
                        ].values[0],
                        self.__rows_data.loc[
                            self.__rows_data["rang"] == ir, "left_lon"
                        ].values[0],
                    )
                )
                liste_gps_droite.append(
                    (
                        self.__rows_data.loc[
                            self.__rows_data["rang"] == ir, "right_lat"
                        ].values[0],
                        self.__rows_data.loc[
                            self.__rows_data["rang"] == ir, "right_lon"
                        ].values[0],
                    )
                )
            # chemin à gauche
            dist_g = self.__distances([gps_depart] + liste_gps_gauche + [gps_arrivee])
            # chemin à droite
            dist_d = self.__distances([gps_depart] + liste_gps_droite + [gps_arrivee])
            # on compare pour ajouter les extrémités du cote le plus court
            if dist_g < dist_d:
                chemin += liste_gps_gauche
                distot += dist_g
            else:
                chemin += liste_gps_droite
                distot += dist_d
        else:
            distot += distance.geodesic(gps_arrivee, gps_depart).m

        return distot, chemin

    def optimise_route(self, start_site, sites):
        # si start_site ou les numeros de sites ne sont pas tous des indices de la dataframe geo_data, on déclenche une erreur
        if not (start_site in self.__geo_data.index) or not all(
            [site in self.__geo_data.index for site in sites]
        ):
            raise ValueError(
                "the sites must be part of the values in the Site column of the geo_data dataframe"
            )

        self.start_site = start_site
        self.__sites = sites
        self.__compute_matrix()

        sites = [self.start_site] + [self.__sites[i] for i in range(len(self.__sites))]

        nb_etapes = len(self.__sites) + 1
        nb_sommets = 2 * len(self.__sites) + 1

        matrice_distances = cp.cpm_array(self.__distance_matrix)
        maxi = max(max(ligne) for ligne in matrice_distances)

        model = cp.Model()

        order = cp.intvar(lb=0, ub=(nb_etapes - 1), name="order", shape=nb_etapes)
        cote = cp.intvar(lb=-1, ub=0, name="cote", shape=nb_etapes)
        pos = cp.intvar(lb=0, ub=nb_sommets, name="pos", shape=nb_etapes)

        distances = cp.intvar(lb=0, ub=maxi, name="distance", shape=nb_etapes)

        model += cp.AllDifferent(order)

        model += [order[0] == 0]
        model += [cote[0] == 0]

        cp.circuit

        # casser symétrie
        model += [order[1] > order[nb_etapes - 1]]

        for i in range(0, nb_etapes):
            model += [pos[i] == 2 * order[i] + cote[i]]
            model += distances[i] == matrice_distances[pos[i], pos[(i + 1) % nb_etapes]]

        travel_distance = cp.sum(distances)

        model.minimize(travel_distance)

        solution = [(self.start_site, Pica.UP, 0)]

        # Résoudre le modèle
        if model.solve():

            gps_circuit = [self.coordonnes_gps_site(self.start_site)]

            sol_distances = distances.value()
            sol_order = order.value()
            sol_cote = cote.value()

            disttot = 0
            for idx in range(1, nb_etapes):
                if sol_cote[idx] == -1:
                    cote = Pica.DOWN
                else:
                    cote = Pica.UP
                lesite = sites[sol_order[idx]]
                dist, chemin = self.trajet_vers(
                    solution[-1][0], solution[-1][1], lesite, cote
                )
                disttot += dist
                # print(f"{predi}->{i} {lesite} dist={dist}  matrix={matrice_distances_cp[predi][i]}")
                gps_circuit += chemin + [self.coordonnes_gps_site(lesite)]
                solution.append((lesite, cote, dist))

            # print(self.start_site)
            dist, chemin = self.trajet_vers(
                solution[-1][0], solution[-1][1], self.start_site, Pica.UP
            )
            disttot += dist
            gps_circuit += chemin + [self.coordonnes_gps_site(self.start_site)]
            solution.append((self.start_site, Pica.UP, dist))
            # print("Distance totale =", disttot, travel_distance.value() / Pica.PRECISION)
            #
            return disttot, solution, gps_circuit
        else:
            return None, None, None
