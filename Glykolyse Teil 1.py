#Importieren der Bibliotheken 
#numpy --> numerische Berechnungen
#matplotlib.pyplot --> visualisierung von Daten
#networkx --> modellieren und zeichnen von Netzwerken

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# --- Objektorientierte Klassen --- für Metabolite Aufgabe 1

class Metabolite:                               #Modelliert Metaboliten in biochemischen Netzwerk 
    def __init__(self, name, initial_conc):     #erstellung Metabolit
        self.name = name                        #Name Molekül
        self.conc = initial_conc                # Ausgangskonzentration (mM)
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
class Enzyme:                                   #Enzym Klasse mit Name, Maximale Geschwindikeit, Michaelis-Menten Konstante
    def __init__(self, name, kcat, enzyme_conc, km = 1.0):
        self.name = name
        self.enzyme_conc = enzyme_conc         #Schätzung von Enzymaktivität und Zellvolumen in mM
        self.kcat = kcat                       #kcat in 1/s
        self.vmax = kcat * enzyme_conc         #Berechnung vmax aus kcat und enzym Konzentration in mM/s
        self.km = km                           #km in mM

    def rate(self, substrate_conc):            #Berechnung Rate nach Michaelis-Menten
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
        return self.enzyme.rate(self.substrate.conc)           #aktuelle Reaktionsgeschwindigkeit und übergibt Wert an aktuelle Substrat Konz
    
    def step(self, dt):                                        #zeitliche Änderung der KOnzentration in Zeitabschnitt
        v = self.rate()                                        #Momentane Reaktionsv
        delta = v * dt                                         #Substratmenge die umgesetzt wird
        self.substrate.update_conc(-delta)
        self.product.update_conc(delta)

    def __repr__(self):
        return f"Reaction({self.name}, Enzyme={self.enzyme})"
    
class SplitReaction:                                                #chatgpt prompt:wie kann ich eine Reaktion modellieren, bei der ein Substrat in zwei Produkte gespalten wird (z. B. F1,6BP → DHAP + G3P)?
    def __init__(self, name, substrate, product1, product2, enzyme):
        self.name = name
        self.substrate = substrate
        self.product1 = product1
        self.product2 = product2
        self.enzyme = enzyme

    def rate(self):
        return self.enzyme.rate(self.substrate.conc)

    def step(self, dt):
        v = self.rate()
        delta = v * dt
        self.substrate.update_conc(-delta)
        self.product1.update_conc(delta / 2)                #je 50%
        self.product2.update_conc(delta / 2)

    def __repr__(self):
        return f"SplitReaction({self.name}, Enzyme={self.enzyme})"

# -----------------------------
# Klasse für den Stoffwechselweg
# -----------------------------
class GlycolysisPathway:
    def __init__(self):                   #definieren der Metabolite und Ausgangskonzentrationen                 
        self.glucose = Metabolite("Glukose", 10.0)                                               #Ausgangskonz. von Glucose 10 mM
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
                                                                                                          #Enzyme definieren und Reaktionsgeschwindigkeit 
        self.hexokinase = Enzyme("Hexokinase", kcat=200, enzyme_conc=0.0025, km=0.05)                     #Hexokinase angaben aus Hecokinase of Human Erythrocytes G. Gerber et al.
        self.isomerase = Enzyme("Glukose-6-phosphat-Isomerase", kcat=150, enzyme_conc=0.002, km=0.1)      #aus Studies on human triosephosphate isomerase. I. Isolation and properties of the enzyme from erythrocytes E.E. Rozacky et al.
        self.pfk = Enzyme("Phosphofructokinase", kcat=300, enzyme_conc=0.002, km=0.08)                    #aus Type 2 diabetes differentially affects the substrate saturation kinetic attributes of erythrocyte hexokinase and phosphofructokinase S. Katyare et al.
        self.aldolase = Enzyme("Aldolase", kcat=100, enzyme_conc=0.0015, km=0.03)                         #aus Human aldolase A natural mutants: relationship between flexibility of the C-terminal region and enzyme function G. Esposito et al.
        self.triosephosphat_isomerase = Enzyme("Triosephosphatisomerase", kcat=4300, enzyme_conc=0.002, km=0.6)               #aus Lehrbuch der Biochemie D.Voet et al. 
        self.gapdh = Enzyme("Glycerinaldehyd-3-phosphat-Dehydrogenase", kcat=250, enzyme_conc=0.002, km=0.02)                                #aus Immunoaffinity purification and characterization of glyceraldehyde-3-phosphate dehydrogenase from human erythrocytes D. Mountassif
        self.pgk = Enzyme("3-Phosphoglyceratkinase", kcat=300, enzyme_conc=0.002, km=0.2)                                     #aus A study of phosphoglycerate kinase in human erythrocytes. I. Enzyme isolation, purification and assay M. Ali, Y. S. Brownstone
        self.pgm = Enzyme("Phosphoyglyceratmutase", kcat=100, enzyme_conc=0.002, km=0.15)                                    #aus Lehrbuch der Biochemie D.Voet et al. 
        self.enolase = Enzyme("Enolase", kcat=200, enzyme_conc=0.002, km=0.07)                            #aus Riboregulation of Enolase 1 activity controls glycolysis and embryonic stem cell differentiation I. Huppertz et al.
        self.pyruvate_kinase = Enzyme("Pyruvat-Kinase", kcat=350, enzyme_conc=0.002, km=0.07)             #aus Structure and Function of Human Erythrocyte Pyruvate Kinase: MOLECULAR BASIS OF NONSPHEROCYTIC HEMOLYTIC ANEMIA G.Valentini et al.
                                                                                                           
        self.reactions = [                                                                                #Reaktionen definieren
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




    def simulate(self, steps=100, dt=1.0):                  #Simulation wird durchgeführt mit 100 Schritten mit 1.0 Schrittweite
        history = {                                         #Konzentrationswerte von Glucose,G6P, F6P und Pyruvat
            "Glukose": [],
            "Glukose-6-phosphat": [],
            "Fruktose-6-phosphat": [],
            "Fruktose-1,6-bisphosphat": [],
            "Dihydroxyacetonphosphat + Glycerinaldehyd-3-phosphat" : [],
            "Glycerinaldehyd-3-phosphat" : [],
            "1,3-Bisphosphoglycerat" : [],
            "3-Phosphoglycerat" : [],
            "2-Phosphoglycerat" : [],
            "Phosphoenolpyruvat" : [],
            "Pyruvat" :[], 
        }

        for _ in range(steps):                              #Simulation für angegebene Schritte wird ausgeführt
            for reaction in self.reactions:
                reaction.step(dt)
            history["Glukose"].append(self.glucose.conc)    #speichern der Konzentrationen
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

        return history
# -----------------------------
# Ausführung & Visualisierung
# -----------------------------
if __name__ == "__main__":
    model = GlycolysisPathway()                     #Modell wird erstellt mit allen Metaboliten  Enzymen und Reaktionen
    data = model.simulate(steps=100, dt=1.0)         #Simulation über 100 Zeitschritte mit 1 EInheiten --> gespeichert in Data

    for met, values in data.items():                #für jeden Metabolit wird der zeitliche Verlauf als Linie geplottet
        plt.plot(values, label=met)
        
    plt.xlabel("Zeit (a.u.)")                             #x-Wert: Zeit
    plt.ylabel("Konzentration (mM)")                      #y-Wert: Konzentration
    plt.title("Glykolyse: Metabolitenkonzentrationen")   #Diagrammtitel: Glykolyse: Metabolitenkonzentrationen
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


