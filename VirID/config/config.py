import os
from os.path import join
import sys

LOG_TASK = 21
VirID_path = os.path.abspath(os.path.join(os.path.dirname(__file__),os.path.pardir))


Package_PATH = join(VirID_path, 'external')

RdRP_DB_PATH = join(VirID_path, 'data', 'diamond_database', 'RdRP_230330_rmdup')
RdRP_DB_TABLE_PATH = join(VirID_path, 'data', 'diamond_database', 'RdRp_expanded_20230330.csv')
Pathogenicity_VIRUS_CSV = join(VirID_path, 'data', 'Pathogenicity_VIRUS_CSV.csv')
RdRP_ICTV_NR = join(VirID_path, 'data', 'RdRP_ICTV_NR.csv')
ICTV_TAXON_CSV = join(VirID_path, 'data', 'ICTV_TAXON.csv')
INTER_RdRP_CLUSTR_PATH = join(VirID_path, 'data', 'RdRP_tree_230328')
TAXONKIT_DB = join(VirID_path, 'data', 'taxonkit_db')


try:
    # expandvars is required to transform things like $HOME
    VirID_DB_PATH = os.path.expandvars(os.environ['VirID_DB_PATH'])
except KeyError:
    print("The 'VirID_DB_PATH' environment variable is not defined.".center(80) + '\n')
    sys.exit(1)


NR_DB_PATH = VirID_DB_PATH+"/NR/nr"
NT_DB_PATH = VirID_DB_PATH+"/NT/nt"
rRNA_DB_PATH = VirID_DB_PATH+'/rRNA/rRNA_cutout_ref'
PROT_ACC2TAXID = VirID_DB_PATH+'/accession2taxid/prot.accession2taxid'



AA_LENGTH_DICT = {
    "Astro-Poty":200,
    "Birna-Permutotetra":200,
    "Negative_Bunya-Arena":200,
    "Cysto":200,
    "Flavi":200,
    "Hepe-Virga":200,
    "Hypo":200,
    "Luteo-Sobemo":200,
    "Negative_Mono-Chu":200,
    "Narna-Levi":200,
    "Nido":200,
    "Negative_Orthomyxo":200,
    "Partiti-Picobirna":200,
    "Picorna-Calici":400,
    "Weivirus":200,
    "Yanvirus":200,
    "Yuevirus":200,
    "Zhaovirus":200,
    "Qinvirus":200,
    "Reo":200,
    "na-daxi":200, # Quenya

    "other-negative":200,
    "Tombus-Noda":200,
    "Toti-Chryso":200
}




Pathogenicity_VIRUS_GROUP = {
"Astro-Poty":["Astroviridae_Mamastrovirus","Potyviridae","Astroviridae","Astroviridae_Avastrovirus"],
"Birna-Permutotetra":["Birnaviridae"],
"Negative_Bunya-Arena":["Hantaviridae","Fimoviridae","Arenaviridae",
                        "Nairoviridae_Orthonairovirus","Peribunyaviridae",
                        "Phenuiviridae","Tospoviridae","Phenuiviridae_Bandavirus",
                        "Phenuiviridae_Phlebovirus"],
"Flavi":["Orthoflavivirus","Hepacivirus","Pegivirus","Pestivirus"],
"Hepe-Virga":["Alphaflexiviridae","Benyviridae","Betaflexiviridae",
              "Bromoviridae","Closteroviridae",
              "Hepe-Virga","Hepeviridae","Idaeovirus",
              "Togaviridae_Alphavirus","Togaviridae_Rubivirus","Tymoviridae","Virgaviridae"],
"Negative_Mono-Chu":["Bornaviridae","Filoviridae","Paramyxoviridae","Pneumoviridae","Rhabdoviridae_Betarhabdovirinae",
                     "Rhabdoviridae_Alpharhabdovirinae"],
"Narna-Levi":["Ourmiavirus"],
"Luteo-Sobemo":["Solemoviridae_Polerovirus","Solemoviridae_Sobemovirus","Solemoviridae_Enamovirus"],
"Nido":["Coronaviridae","Arteriviridae"],
"Picorna-Calici":["Caliciviridae","Picornaviridae","Secoviridae"],
"Negative_Orthomyxo":["Thogotovirus","Influenzavirus"],
"Partiti-Picobirna":["Partiti-Picobirna","Partitiviridae",
                     "Partitiviridae_Deltapartitivirus",
                     "Picobirnaviridae"],
"Reo":["Orbivirus","Orthoreovirus","Oryzavirus","Rotavirus","Phytoreovirus","Reoviridae",
       "Seadornavirus","Coltivirus","Fijivirus","Aquareovirus"],
"Tombus-Noda":["Tombusviridae","Tombusviridae_Calvusvirinae","Tombusviridae_Procedovirinae",
               "Tombusviridae_Regressovirinae"]
}



Group_Default_Taxon = {
"Astro-Poty":"Astro-Poty",
"Birna-Permutotetra":"Birna-Permutotetra",
"Negative_Bunya-Arena":"Negative_Bunya-Arena",
"Flavi":"Flaviviridae",
"Hepe-Virga":"Hepe-Virga",
"Negative_Mono-Chu":"Negative_Mono-Chu",
"Hypo":"Hypoviridae",
"Narna-Levi":"Narna-Levi",
"Luteo-Sobemo":"Luteo-Sobemo",
"Nido":"Nido",
"Picorna-Calici":"Picorna-Calici",
"Negative_Orthomyxo":"Orthomyxoviridae",
"Partiti-Picobirna":"Partiti-Picobirna",
"Reo":"Reoviridae",
"Tombus-Noda":"Tombus-Noda",
"Toti-Chryso":"Toti-Chryso"
}



