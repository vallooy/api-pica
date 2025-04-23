from fastapi import FastAPI
from algo_pica_cpmpy import Pica
from pathlib import Path
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

class Parcours(BaseModel):
    parcelle: str
    millesime: str
    parcours: list

app = FastAPI()

origins = [
    "http://localhost:8080",
    "http://localhost:4200",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def root():
    return {"message": "Martin rocks!"}

@app.post("/optimise/")
def optimise(parcours: Parcours):
    ## Variable pour test
    # point de depart
    premier_site = 1
    # valeur fixe pour test
    ceps = [116, 310, 316, 528, 701, 921, 1327]
    ##---- Fin de variable de test
    print(parcours)
    #Import des fichiers
    json_path = Path(__file__).parent
    nom_fichier_geopandas = json_path / f"{parcours.parcelle}_point_algo_v2_rang.json"
    nom_fichier_extremites = json_path / f"{parcours.parcelle}_extremites_rangs.csv"

    #Création de l'instance Pica 
    algo = Pica(nom_fichier_extremites, nom_fichier_geopandas)

    #Paramètre du premier site
    if parcours.parcelle == "Larzat":
        premier_site = 9
    elif parcours.parcelle == "Arnel":
        premier_site = 38
    elif parcours.parcelle == "Estagnol":
        premier_site = 732
    
    longueur, solution, gps_circuit = algo.optimise_route(premier_site, parcours.parcours)
    if solution is None:
        print("Pas de solution")
    else:
        print("Distance totale =", longueur)
        print(gps_circuit)
    
    return longueur, solution, gps_circuit
    
    
    

