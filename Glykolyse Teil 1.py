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
    def __init__(self, name, rate, km = 1.0):
        self.name = name
        self.vmax = rate
        self.km = km

    def rate(self, substrate_conc):            #Berechnung Rate nach Michaelis-Menten
        s = substrate_conc              
        if s > 0:                              #Substrat > 0 --> Berechnung vmax
            return self.vmax * s / (self.km + s)
        else:
            return 0                           #Substrat < 0 --> Geschwindigkeit = o
    def __repr__(self):
        return f"Enzyme({self.name}, vmax={self.vmax}, km={self.km})"

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

# -----------------------------
# Klasse für den Stoffwechselweg
# -----------------------------
class GlycolysisPathway:
    def __init__(self):                                    
        self.glucose = Metabolite("Glukose", 10.0)          #Metaboliten definieren und Ausgangskonzentration
        self.g6p = Metabolite("Glucose-6-Phosphat", 0.0)
        self.f6p = Metabolite("Fructose-6-Phosphat", 0.0)
        self.pyruvate = Metabolite("Pyruvat", 0.0)

        self.hexokinase = Enzyme("Hexokinase", rate=0.5)     #Enzyme definieren und Reaktionsgeschwindigkeit
        self.isomerase = Enzyme("Isomerase", rate=0.4)
        self.pyruvate_kinase = Enzyme("Pyruvat-Kinase", rate=0.7)
                                                            #Reaktionen definieren
        self.reactions = [
            Reaction("Glukose → G6P", self.glucose, self.g6p, self.hexokinase),
            Reaction("G6P → F6P", self.g6p, self.f6p, self.isomerase),
            Reaction("F6P → Pyruvat", self.f6p, self.pyruvate, self.pyruvate_kinase),
        ]

    def simulate(self, steps=100, dt=1.0):                  #Simulation wird durchgeführt mit 100 Schritten mit 1.0 Schrittweite
        history = {                                         #Konzentrationswerte von Glucose,G6P, F6P und Pyruvat
            "Glukose": [],
            "G6P": [],
            "F6P": [],
            "Pyruvat": []
        }

        for _ in range(steps):                              #Simulation für angegebene Schritte wird ausgeführt
            for reaction in self.reactions:
                reaction.step(dt)
            history["Glukose"].append(self.glucose.conc)    #Speichern der jeweiligen Werte
            history["G6P"].append(self.g6p.conc)
            history["F6P"].append(self.f6p.conc)
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


