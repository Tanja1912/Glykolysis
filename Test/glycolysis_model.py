# glycolysis_model.py

# Importieren der notwendigen Bibliotheken
import numpy as np
import networkx as nx

# --- Klassen für die Modellierung ---

class Metabolite:
    """
    Modelliert einen Metaboliten im biochemischen Netzwerk.
    Behält den Überblick über die Konzentration und deren Verlauf.
    """
    def __init__(self, name, initial_conc=0.0):
        self.name = name
        self.conc = initial_conc  # Anfangskonzentration in mM
        self.history = [initial_conc]  # Speichert den Konzentrationsverlauf

    def update_conc(self, delta):
        """Aktualisiert die Konzentration und stellt sicher, dass sie nicht negativ wird."""
        self.conc += delta
        self.conc = max(self.conc, 0) # Konzentration kann nicht < 0 sein
        self.history.append(self.conc)

    def __repr__(self):
        """Textuelle Darstellung des Metaboliten."""
        return f"{self.name}: {self.conc:.4f} mM"

class Enzyme:
    """
    Modelliert ein Enzym mit Michaelis-Menten-Kinetik.
    """
    def __init__(self, name, kcat, enzyme_conc, km=1.0):
        self.name = name
        self.enzyme_conc = enzyme_conc  # Enzymkonzentration in mM
        self.kcat = kcat  # katalytische Konstante (Substratumsätze pro Sekunde)
        self.vmax = kcat * enzyme_conc  # Maximale Reaktionsgeschwindigkeit in mM/s
        self.km = km  # Michaelis-Konstante in mM

    def rate(self, substrate_conc):
        """Berechnet die Reaktionsgeschwindigkeit basierend auf der Substratkonzentration."""
        s = substrate_conc
        if s > 0:
            return self.vmax * s / (self.km + s)
        return 0.0

    def __repr__(self):
        """Textuelle Darstellung des Enzyms."""
        return f"Enzyme({self.name}, vmax={self.vmax:.4f}, km={self.km})"

class Reaction:
    """
    **Verbesserte, generalisierte Reaktionsklasse.**
    Kann Reaktionen mit mehreren Substraten und Produkten und beliebiger Stöchiometrie modellieren.
    Ersetzt die Notwendigkeit einer separaten 'SplitReaction'-Klasse.
    """
    def __init__(self, name, substrates, products, enzyme):
        self.name = name
        # Substrate/Produkte werden als Listen von Tupeln übergeben: [(metabolit, koeffizient), ...]
        self.substrates = substrates
        self.products = products
        self.enzyme = enzyme

    def rate(self):
        """
        Berechnet die aktuelle Reaktionsgeschwindigkeit.
        Für Einfachheit wird die Rate durch das erste Substrat in der Liste bestimmt.
        """
        # Annahme: Das erste Substrat ist das primäre, rate-bestimmende Substrat.
        primary_substrate = self.substrates[0][0]
        return self.enzyme.rate(primary_substrate.conc)

    def step(self, dt):
        """
        Führt einen Simulationsschritt für die Reaktion durch.
        Aktualisiert die Konzentrationen aller beteiligten Metaboliten gemäß ihrer Stöchiometrie.
        """
        v = self.rate()
        delta = v * dt  # Gesamtumsatz in diesem Zeitschritt

        # Verbrauche Substrate
        for metabolite, coeff in self.substrates:
            metabolite.update_conc(-delta * coeff)

        # Produziere Produkte
        for metabolite, coeff in self.products:
            metabolite.update_conc(delta * coeff)

    def __repr__(self):
        """Textuelle Darstellung der Reaktion."""
        return f"Reaction({self.name})"

# --- Hauptklasse für den Stoffwechselweg ---

class GlycolysisPathway:
    """
    Modelliert den gesamten Glykolyse-Weg.
    Initialisiert Metaboliten, Enzyme und Reaktionen.
    Führt die Simulation durch und kann einen Graphen des Weges erstellen.
    """
    def __init__(self, initial_glucose=10.0, enzyme_configs=None):
        if enzyme_configs is None:
            enzyme_configs = {} # Leeres dict, falls keine Konfiguration übergeben wird

        # Metaboliten definieren
        self.glucose = Metabolite("Glukose", initial_glucose)
        self.g6p = Metabolite("G6P", 0.0)
        self.f6p = Metabolite("F6P", 0.0)
        self.f1_6bp = Metabolite("F1,6BP", 0.0)
        self.dhap = Metabolite("DHAP", 0.0)
        self.g3p = Metabolite("G3P", 0.0)
        self.bpg_1_3 = Metabolite("1,3-BPG", 0.0)
        self.pg_3 = Metabolite("3-PG", 0.0)
        self.pg_2 = Metabolite("2-PG", 0.0)
        self.pep = Metabolite("PEP", 0.0)
        self.pyruvate = Metabolite("Pyruvat", 0.0)
        
        # Co-Faktoren hinzufügen für eine realistischere Modellierung
        self.atp = Metabolite("ATP", 5.0)  # Hohe Anfangskonz.
        self.adp = Metabolite("ADP", 0.5)  # Niedrige Anfangskonz.
        self.nad_plus = Metabolite("NAD+", 2.0)
        self.nadh = Metabolite("NADH", 0.1)

        # Sammle alle Metaboliten in einem Dictionary für einfachen Zugriff
        self.metabolites = {m.name: m for m in [
            self.glucose, self.g6p, self.f6p, self.f1_6bp, self.dhap, self.g3p,
            self.bpg_1_3, self.pg_3, self.pg_2, self.pep, self.pyruvate,
            self.atp, self.adp, self.nad_plus, self.nadh
        ]}

        # Enzyme definieren mit Standardwerten
        # Die Konzentrationen können durch 'enzyme_configs' überschrieben werden
        self.hexokinase = Enzyme("Hexokinase", kcat=200, enzyme_conc=enzyme_configs.get("Hexokinase", 0.0025), km=0.05)
        self.isomerase = Enzyme("Isomerase", kcat=150, enzyme_conc=0.002, km=0.1)
        self.pfk = Enzyme("PFK", kcat=300, enzyme_conc=enzyme_configs.get("PFK", 0.002), km=0.08)
        self.aldolase = Enzyme("Aldolase", kcat=100, enzyme_conc=0.0015, km=0.03)
        self.triosephosphat_isomerase = Enzyme("TIM", kcat=4300, enzyme_conc=0.002, km=0.6)
        self.gapdh = Enzyme("GAPDH", kcat=250, enzyme_conc=0.002, km=0.02)
        self.pgk = Enzyme("PGK", kcat=300, enzyme_conc=0.002, km=0.2)
        self.pgm = Enzyme("PGM", kcat=100, enzyme_conc=0.002, km=0.15)
        self.enolase = Enzyme("Enolase", kcat=200, enzyme_conc=0.002, km=0.07)
        self.pyruvate_kinase = Enzyme("Pyruvat-Kinase", kcat=350, enzyme_conc=enzyme_configs.get("Pyruvat-Kinase", 0.002), km=0.07)

        # Reaktionen mit korrekter Stöchiometrie definieren
        self.reactions = [
            Reaction("Hexokinase", [(self.glucose, 1), (self.atp, 1)], [(self.g6p, 1), (self.adp, 1)], self.hexokinase),
            Reaction("Isomerase", [(self.g6p, 1)], [(self.f6p, 1)], self.isomerase),
            Reaction("PFK", [(self.f6p, 1), (self.atp, 1)], [(self.f1_6bp, 1), (self.adp, 1)], self.pfk),
            Reaction("Aldolase", [(self.f1_6bp, 1)], [(self.dhap, 1), (self.g3p, 1)], self.aldolase),
            Reaction("TIM", [(self.dhap, 1)], [(self.g3p, 1)], self.triosephosphat_isomerase),
            Reaction("GAPDH", [(self.g3p, 1), (self.nad_plus, 1)], [(self.bpg_1_3, 1), (self.nadh, 1)], self.gapdh),
            Reaction("PGK", [(self.bpg_1_3, 1), (self.adp, 1)], [(self.pg_3, 1), (self.atp, 1)], self.pgk),
            Reaction("PGM", [(self.pg_3, 1)], [(self.pg_2, 1)], self.pgm),
            Reaction("Enolase", [(self.pg_2, 1)], [(self.pep, 1)], self.enolase),
            Reaction("Pyruvat-Kinase", [(self.pep, 1), (self.adp, 1)], [(self.pyruvate, 1), (self.atp, 1)], self.pyruvate_kinase),
        ]

    def simulate(self, steps=100, dt=1.0):
        """Führt die Simulation über eine gegebene Anzahl von Zeitschritten durch."""
        for _ in range(steps):
            for reaction in self.reactions:
                reaction.step(dt)
        
        # Sammelt die Verlaufsdaten aller Metaboliten für die Visualisierung
        history = {name: met.history for name, met in self.metabolites.items()}
        return history

    def create_pathway_graph(self):
        """
        **Neue Funktion: Erstellt einen gerichteten Graphen des Stoffwechselweges mit NetworkX.**
        """
        G = nx.DiGraph()
        for reaction in self.reactions:
            # Knoten für die Reaktion selbst hinzufügen
            reaction_node = f"{reaction.enzyme.name}"
            G.add_node(reaction_node, type='enzyme')
            
            # Kanten von Substraten zur Reaktion
            for sub, _ in reaction.substrates:
                G.add_edge(sub.name, reaction_node)

            # Kanten von der Reaktion zu den Produkten
            for prod, _ in reaction.products:
                G.add_edge(reaction_node, prod.name)
        return G