# Importieren der Bibliotheken 

import numpy as np   # numpy --> numerische Berechnungen


# Objektorientierte Klassen

# Klasse für Metabolitesc


class Metabolite:                               # Modelliert Metaboliten in biochemischen Netzwerk 
    
    
    def __init__(self, name, initial_conc):              # Erstellung Metabolit
        self.name = name                        # Name Molekül
        self.conc = initial_conc                # Ausgangskonzentration (mMol/L)
        self.history = [initial_conc]           # Speicherung der Werte
    
    def update_conc(self, delta):               # Aktualisiert Konzentration der Metaboliten basierend auf delta
        self.conc += delta                      # Veränderung Konzentration um delta
        self.conc = max(self.conc, 0)           # Sicherstellen dass Konzentration nicht < 0
        self.history.append(self.conc)          # Speichern der Konzentration für spätere Darstellung
        
    def __repr__(self):                            # Darstellung von Objekt als Text
        return f"{self.name}: {self.conc:.3f} mM"  # formatierter String mit Name des Metaboliten mit Konzentration auf 3 Dezimalstellen


# Klasse für Enzyme


class Enzyme:                                   # Enzym Klasse mit Name,Turnover number (kcat), Enzymkonzentration und Michaelis-Menten-Konstante (km)
   
    
    def __init__(self, name, kcat, enzyme_conc, km=1.0):
        self.name = name
        self.enzyme_conc = enzyme_conc         # Schätzung von Enzymaktivität und Zellvolumen in mMol
        self.kcat = kcat                       # kcat (Turnover Number) in 1/s
        self.vmax = kcat * enzyme_conc         # Berechnung vmax (maximale Reaktionsgeschwindigkeit) aus kcat und Enzymkonzentration in mMol/s
        self.km = km                           # km in mMol

    
    def rate(self, substrate_conc):            # Berechnung Enzymgeschwindigkeit basierend auf der Substratkonzentration nach Michaelis-Menten
        
        s = substrate_conc              
        if s > 0:                              # wenn Substrat > 0 berechne die Geschwindigkeit (vmax) nach Michaelis-Menten 
            return self.vmax * s / (self.km + s)  # Berechnung nach Michaelis-Menten-Gleichung
        else:
            return 0                           # Substrat < 0 --> Geschwindigkeit = 0 --> keine negative Geschwindigkeit
   
    
    def __repr__(self):
        return (f"Enzyme({self.name}, kcat={self.kcat}, [E]={self.enzyme_conc}, "
                f"vmax={self.vmax:.3f}, km={self.km})")


#Klasse für Reaktionen


class Reaction:
    
    
    def __init__(self, name, substrate, product, enzyme):                                      # initialisiert enzymkatalysierte Reaktion
        self.name = name                                       # Name der Reaktion
        self.substrate = substrate                             # Substrat 
        self.product = product                                 # entstehendes Produkt
        self.enzyme = enzyme                                   # verwendetes Enzym
    
    
    def rate(self):
        return self.enzyme.rate(self.substrate.conc)           # Berechnet die aktuelle Reaktionsgeschwindigkeit basierend auf atueller Substratkonzentration
    
   
    def step(self, dt):                                            # führt Zeitschritt der Reaktion durch und aktualisiert die Konzentrationen
        v = self.rate()                                            # Momentane Reaktionsgeschwindigkeit
        delta = v * dt                                             # Substratmenge die umgesetzt wird   
        delta = min(delta, self.substrate.conc)                    # sicherstellen, dass nur so viel Substrat umgesetzt wird wie tatsächlich vorhanden ist
        self.substrate.update_conc(-delta)                         # Abnahme der Substratkonzentration um delta
        self.product.update_conc(delta)                            # Zunahme der Produktkonzentration um delta


    def __repr__(self):
        return f"Reaction({self.name}, Enzyme={self.enzyme})"
        

# Klasse für Splitreaktion --> Spezialfall Spaltung eines Substrats in 2 Produkte


class SplitReaction:                                                
    
    
    def __init__(self, name, substrate, product1, product2, enzyme):
        self.name = name
        self.substrate = substrate
        self.product1 = product1                                    # Spaltung der entstehenden Produkte in Produkt 1 und Produkt 2
        self.product2 = product2
        self.enzyme = enzyme

    
    def rate(self):
        return self.enzyme.rate(self.substrate.conc)                

    
    def step(self, dt):                                             
        v = self.rate()                                            
        delta = v * dt                                             
        delta = min(delta, self.substrate.conc)                    
        self.substrate.update_conc(-delta)                          
        self.product1.update_conc(delta)                            # Zunahme beider Produkte um delta --> korrektes Verhältnis der Umsetzung aus Fruktose-1,6-bisphosphat → Dihydroxyacetonphosphat + Glycerinaldehyd-3-phosphat
        self.product2.update_conc(delta)


    def __repr__(self):
        return f"SplitReaction({self.name}, Enzyme={self.enzyme})"


# Klasse für den Stoffwechselweg, mit allen Metaboliten, Enzymen und den jeweiligen Reaktionen


class GlycolysisPathway:
    
    
    def __init__(self, glucose_conc=10.0):             # Übergabe Anfangskonzentration von Glukose (10mMol/L); Platzhalter: kann in Streamlit angepasst werden oder beim Erstellen eines Objekts der Klasse 
        self.glucose = Metabolite(
            name="Glukose", 
            initial_conc=glucose_conc
        )                                              # definieren der Metabolite mit spezifischen Namen und Ausgangskonzentrationen (alle starten bei 0 außer Glukose)                         
        self.g6p = Metabolite(
            name="Glukose-6-phosphat",
            initial_conc=0.0
        )
        self.f6p = Metabolite(
            name="Fruktose-6-phosphat",
            initial_conc=0.0
        )
        self.f1_6bp = Metabolite(
            name="Fruktose-1,6-bisphosphat",
            initial_conc=0.0
        )
        self.dhap = Metabolite(
            name="Dihydroxyacetonphosphat",
            initial_conc=0.0
        )
        self.g3p = Metabolite(
            name="Glycerinaldehyd-3-phosphat",
            initial_conc=0.0
        )
        self.bpg_1_3 = Metabolite(
            name="1,3-Bisphosphoglycerat", 
            initial_conc=0.0
        )
        self.pg_3 = Metabolite(
            name="3-Phosphoglycerat", 
            initial_conc=0.0
        )
        self.pg_2 = Metabolite(
            name="2-Phosphoglycerat", 
            initial_conc=0.0
        )
        self.pep = Metabolite(
            name="Phosphoenolpyruvat", 
            initial_conc=0.0
        )
        self.pyruvate = Metabolite(
            name="Pyruvat", 
            initial_conc=0.0
        ) 
                                                      # definieren der Enzyme mit Namen, Turnover Number (kcat), typischen Enzymkonzentrationen sowie Michaelis-Menten-Konstante(km) → Annäherungswerte aus der Literatur → Anpassung über streamlit
        self.hexokinase = Enzyme(
            name="Hexokinase", 
            kcat=200, 
            enzyme_conc=0.0025, 
            km=0.05
        )                    
        self.isomerase = Enzyme(
            name="Glukose-6-phosphat-Isomerase",
            kcat=150, 
            enzyme_conc=0.002,
            km=0.1
        )      
        self.pfk = Enzyme(
            name="Phosphofructokinase",
            kcat=300, 
            enzyme_conc=0.002, 
            km=0.08
        )                    
        self.aldolase = Enzyme(
            name="Aldolase", 
            kcat=100, 
            enzyme_conc=0.0015, 
            km=0.03
        )                         
        self.triosephosphat_isomerase = Enzyme(
            name="Triosephosphatisomerase", 
            kcat=4300, 
            enzyme_conc=0.002,
            km=0.6
        )               
        self.gapdh = Enzyme(
            name="Glycerinaldehyd-3-phosphat-Dehydrogenase",
            kcat=250, 
            enzyme_conc=0.002,
            km=0.02
        )                                
        self.pgk = Enzyme(
            name="3-Phosphoglyceratkinase",
            kcat=300, 
            enzyme_conc=0.002,
            km=0.2
        )                                   
        self.pgm = Enzyme(
            name="Phosphoglyceratmutase",
            kcat=100,
            enzyme_conc=0.002, 
            km=0.15
        )                                    
        self.enolase = Enzyme(
            name="Enolase",
            kcat=200, 
            enzyme_conc=0.002, 
            km=0.07)                            
        self.pyruvate_kinase = Enzyme(
            name="Pyruvat-Kinase", 
            kcat=350, 
            enzyme_conc=0.002, 
            km=0.07
        )             
        
        self.reactions = [                                                                # Reaktionen definieren Ablauf aus Literatur; Umsetzung Substrat zu Produkt durch Enzym
            Reaction(
                name="Glukose → Glukose-6-phosphat",
                substrate=self.glucose, 
                product=self.g6p, 
                enzyme=self.hexokinase
            ),
            Reaction(
                name="Glukose-6-phosphat → Fruktose-6-phosphat", 
                substrate=self.g6p,
                product=self.f6p, 
                enzyme=self.isomerase
            ),
            Reaction(
                name="Fruktose-6-phosphat → Fruktose-1,6-bisphosphat", 
                substrate=self.f6p, 
                product=self.f1_6bp,
                enzyme=self.pfk
            ),
            SplitReaction(
                name="Fruktose-1,6-bisphosphat → Dihydroxyacetonphosphat + Glycerinaldehyd-3-phosphat", 
                substrate=self.f1_6bp,
                product1=self.dhap,
                product2=self.g3p, 
                enzyme=self.aldolase
            ),
            Reaction(
                name="Dihydroxyacetonphosphat → Glycerinaldehyd-3-phosphat", 
                substrate=self.dhap, 
                product=self.g3p,
                enzyme=self.triosephosphat_isomerase
            ),
            Reaction(
                name="Glycerinaldehyd-3-phosphat → 1,3-Bisphosphoglycerat", 
                substrate=self.g3p,
                product=self.bpg_1_3, 
                enzyme=self.gapdh
            ),
            Reaction(
                name="1,3-Bisphosphoglycerat → 3-Phosphoglycerat", 
                substrate=self.bpg_1_3, 
                product=self.pg_3, 
                enzyme=self.pgk
            ),
            Reaction(
                name="3-Phosphoglycerat → 2-Phosphoglycerat", 
                substrate=self.pg_3, 
                product=self.pg_2, 
                enzyme=self.pgm
            ),
            Reaction(
                name="2-Phosphoglycerat → Phosphoenolpyruvat", 
                substrate=self.pg_2, 
                product=self.pep,
                enzyme=self.enolase
            ),
            Reaction(
                name="Phosphoenolpyruvat → Pyruvat",
                substrate=self.pep, 
                product=self.pyruvate, 
                enzyme=self.pyruvate_kinase
            ),
        ]

    
    def simulate(self, steps=100, dt=0.1):                                                                # Speichert Verlauf der Metabolitkonzentrationen über Anzahl von Zeitschritten (100 Zeitschritte mit je dt Sekunden → Anpassung über streamlit) 
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
                                                                                                    # Anfangswerte speichern --> sodass Streamlit bei angegebener Konzentration beginnt und nicht bereits das erste Mal die Schleife durchläuft
        history["Glukose"].append(self.glucose.conc)
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

        for _ in range(steps):                                                                            # Durchlaufen der Schleife steps-mal
            for reaction in self.reactions:                                                               # Aufrufen der Methode step(dt) für jede Reaktion
                reaction.step(dt)                                                                        
                
            history["Glukose"].append(self.glucose.conc)                                                  # neue Konzentration der Metabolite wird in einer Liste gespeichert welche im history-Dictionary unter dem jeweiligen Namen gespeichert wird                                                  
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






