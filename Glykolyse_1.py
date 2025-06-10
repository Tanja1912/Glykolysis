#Importieren der Bibliotheken 
#numpy --> numerische Berechnungen

import numpy as np


# --- Objektorientierte Klassen ---

# -----------------------------
# Klasse für Metabolite
# -----------------------------

class Metabolite:                               #Modelliert Metaboliten in biochemischen Netzwerk 
    def __init__(self, name, initial_conc):     #erstellung Metabolit
        self.name = name                        #Name Molekül
        self.conc = initial_conc                # Ausgangskonzentration (mMol/L)
        self.history = [initial_conc]           #Speicherung der Werte
    
    def update_conc(self, delta):               #Aktualisiert Konzentration der Metaboliten basierend auf delta
        self.conc += delta                      #Veränderung Konzentration um delta
        self.conc = max(self.conc, 0)           #Konzentration darf nicht negativ sein
        self.history.append(self.conc)          #Speichern der Konzentration für spätere Darstellung
        
    def __repr__(self):                         #Darstellung von Objekt als Text
        return f"{self.name}: {self.conc:.3f} mM"  #formatierter String mit Name des Metaboliten mit Konzentration auf 3 Dezimalstellen
# -----------------------------
# Klasse für Enzyme
# -----------------------------
class Enzyme:                                   #Enzym Klasse mit Name,Turnover number, Enzymkonzentration und Michaelis-Menten-Konstante
    def __init__(self, name, kcat, enzyme_conc, km = 1.0):
        self.name = name
        self.enzyme_conc = enzyme_conc         #Schätzung von Enzymaktivität und Zellvolumen in mMol
        self.kcat = kcat                       #kcat in 1/s
        self.vmax = kcat * enzyme_conc         #Berechnung vmax aus kcat und enzym Konzentration in mMol/s
        self.km = km                           #km in mMol

    def rate(self, substrate_conc):            #Berechnung Enzymgeschwindigkeit basierend auf der Substratkonzentration nach Michaelis-Menten
        s = substrate_conc              
        if s > 0:                              #Substrat > 0 --> Berechnung vmax
            return self.vmax * s / (self.km + s)
        else:
            return 0                           #Substrat < 0 --> Geschwindigkeit = o
    def __repr__(self):
        return f"Enzyme({self.name}, kcat={self.kcat}, [E]={self.enzyme_conc}, "f"vmax={self.vmax:.3f}, km={self.km})"

# -----------------------------
# Klasse für Reaktionen
# -----------------------------
class Reaction:
    def __init__(self, name, substrate, product, enzyme):    
        self.name = name
        self.substrate = substrate                             # Substrat 
        self.product = product                                 # entstehendes Produkt
        self.enzyme = enzyme
    def rate(self):
        return self.enzyme.rate(self.substrate.conc)           #Berechnet die aktuelle Reaktionsgeschwindigkeit basierend auf atueller Substratkonzentration
    
    def step(self, dt):                                        #zeitliche Änderung der Konzentration in Zeitabschnitt
        v = self.rate()                                        #Momentane Reaktionsgeschwindigkeit
        delta = v * dt                                         #Substratmenge die umgesetzt wird
        self.substrate.update_conc(-delta)                     #Veränderung der Substrat Konzentration --> sinkt bei -delta 
        self.product.update_conc(delta)                        #Veränderung der Produkt Konzentration --> steigt bei delta

    def __repr__(self):
        return f"Reaction({self.name}, Enzyme={self.enzyme})"
        
# -----------------------------
# Klasse für Splitreaktion
# -----------------------------    
class SplitReaction:                                                
    def __init__(self, name, substrate, product1, product2, enzyme):
        self.name = name
        self.substrate = substrate
        self.product1 = product1                                    #Spaltung der entstehenden Produkte in Produkt 1 und Produkt 2
        self.product2 = product2
        self.enzyme = enzyme

    def rate(self):
        return self.enzyme.rate(self.substrate.conc)                #siehe oben
   
    def step(self, dt):
        v = self.rate()                                             #berechnet Reaktionsgeschwindigkeit
        delta = v * dt                                              #Menge an Substrat die in Zeit dt umgesetzt wird  
        self.substrate.update_conc(-delta)                          #Verringerung Substrat um delta
        self.product1.update_conc(delta / 2)                        #Erhöhung Konzentration der Produkte um je 50%
        self.product2.update_conc(delta / 2)

    def __repr__(self):
        return f"SplitReaction({self.name}, Enzyme={self.enzyme})"

# -----------------------------
# Klasse für den Stoffwechselweg
# -----------------------------
class GlycolysisPathway:
    def __init__(self, glucose_conc=10.0):                                                              #Übergabe Anfangskonzentration von Glukose (10mMol/L); Platzhalter: kann in Streamlit angepasst werden oder beim erstellen eines Objekts der Klasse 
        self.glucose = Metabolite("Glukose", glucose_conc)                                              #definieren der Metabolite mit spezifischen Namen und Ausgangskonzentrationen (alle starten bei 0 außer Glukose)                         
        self.g6p = Metabolite("Glukose-6-phosphat", 0.0)
        self.f6p = Metabolite("Fruktose-6-phosphat", 0.0)
        self.f1_6bp = Metabolite("Fruktose-1,6-bisphosphat", 0.0)
        self.dhap = Metabolite("Dihydroxyacetonphosphat", 0.0)
        self.g3p = Metabolite("Glycerinaldehyd-3-phosphat", 0.0)
        self.bpg_1_3 = Metabolite("1,3-Bisphosphoglycerat", 0.0)
        self.pg_3 = Metabolite("3-Phosphoglycerat", 0.0)
        self.pg_2 = Metabolite("2-Phosphoglycerat", 0.0)
        self.pep = Metabolite("Phosphoenolpyruvat", 0.0)
        self.pyruvate = Metabolite("Pyruvat", 0.0) 
                                                                                                          #definieren der Enzyme mit spezifischen Namen, Turnover Number, typischen Enzymkonzentrationen sowie Michaelis-Menten-Konstante → Werte aus der Literatur 
        self.hexokinase = Enzyme("Hexokinase", kcat=200, enzyme_conc=0.0025, km=0.05)                    
        self.isomerase = Enzyme("Glukose-6-phosphat-Isomerase", kcat=150, enzyme_conc=0.002, km=0.1)      
        self.pfk = Enzyme("Phosphofructokinase", kcat=300, enzyme_conc=0.002, km=0.08)                    
        self.aldolase = Enzyme("Aldolase", kcat=100, enzyme_conc=0.0015, km=0.03)                         
        self.triosephosphat_isomerase = Enzyme("Triosephosphatisomerase", kcat=4300, enzyme_conc=0.002, km=0.6)               
        self.gapdh = Enzyme("Glycerinaldehyd-3-phosphat-Dehydrogenase", kcat=250, enzyme_conc=0.002, km=0.02)                                
        self.pgk = Enzyme("3-Phosphoglyceratkinase", kcat=300, enzyme_conc=0.002, km=0.2)                                   
        self.pgm = Enzyme("Phosphoyglyceratmutase", kcat=100, enzyme_conc=0.002, km=0.15)                                    
        self.enolase = Enzyme("Enolase", kcat=200, enzyme_conc=0.002, km=0.07)                            
        self.pyruvate_kinase = Enzyme("Pyruvat-Kinase", kcat=350, enzyme_conc=0.002, km=0.07)             
        
        self.reactions = [                                                                                #Reaktionen definieren Ablauf aus Literatur; Umsetzung Substrat zu Produkt durch Enzym
            Reaction("Glukose → Glukose-6-phosphat", self.glucose, self.g6p, self.hexokinase),
            Reaction("Glukose-6-phosphat → Fruktose-6-phosphat", self.g6p, self.f6p, self.isomerase),
            Reaction("Fruktose-6-phosphat → Fruktose-1,6-bisphosphat", self.f6p, self.f1_6bp, self.pfk),
            SplitReaction("Fruktose-1,6-bisphosphat → Dihydroxyacetonphosphat + Glycerinaldehyd-3-phosphat", self.f1_6bp, self.dhap, self.g3p, self.aldolase),
            Reaction("Dihydroxyacetonphosphat → Glycerinaldehyd-3-phosphat", self.dhap, self.g3p, self.triosephosphat_isomerase),
            Reaction("Glycerinaldehyd-3-phosphat → 1,3-Bisphosphoglycerat", self.g3p, self.bpg_1_3, self.gapdh),
            Reaction("1,3-Bisphosphoglycerat → 3-Phosphoglycerat", self.bpg_1_3, self.pg_3, self.pgk),
            Reaction("3-Phosphoglycerat → 2-Phosphoglycerat", self.pg_3, self.pg_2, self.pgm),
            Reaction("2-Phosphoglycerat → Phosphoenolpyruvat", self.pg_2, self.pep, self.enolase),
            Reaction("Phosphoenolpyruvat → Pyruvat", self.pep, self.pyruvate, self.pyruvate_kinase),
                    ]
    def simulate(self, steps=100, dt=1.0):                                                                #Verlauf der Konzentrationen der Metabolite wird gespeichert (100 Zeitschritte mit 1 sek Länge → Anpassung über streamlit) über Dictionary
        history = {                            
        "Glukose": [],
        "Glukose-6-phosphat": [],
        "Fruktose-6-phosphat": [],
        "Fruktose-1,6-bisphosphat": [],
        "Dihydroxyacetonphosphat": [],
        "Glycerinaldehyd-3-phosphat": [],
        "1,3-Bisphosphoglycerat": [],
        "3-Phosphoglycerat": [],
        "2-Phosphoglycerat": [],
        "Phosphoenolpyruvat": [],
        "Pyruvat": [],
    }

        for _ in range(steps):                                                                            #durchlaufen der Schleife steps-mal
            for reaction in self.reactions:                                                               #aufrufen der Methode step(dt) für jede Reaktion
                reaction.step(dt)                                                                        
                
            history["Glukose"].append(self.glucose.conc)                                                  #aktuelle Konzentration der Metabolite wird im history Dictionary gespeichert durch append Funktion                                                  
            history["Glukose-6-phosphat"].append(self.g6p.conc)
            history["Fruktose-6-phosphat"].append(self.f6p.conc)
            history["Fruktose-1,6-bisphosphat"].append(self.f1_6bp.conc)
            history["Dihydroxyacetonphosphat"].append(self.dhap.conc)
            history["Glycerinaldehyd-3-phosphat"].append(self.g3p.conc)
            history["1,3-Bisphosphoglycerat"].append(self.bpg_1_3.conc)
            history["3-Phosphoglycerat"].append(self.pg_3.conc)
            history["2-Phosphoglycerat"].append(self.pg_2.conc)
            history["Phosphoenolpyruvat"].append(self.pep.conc)
            history["Pyruvat"].append(self.pyruvate.conc)
        
        return history                                                                                     #Rückgabe des Dictionarys mit den Konzentrationsverläufen über die Zeit                                                             





