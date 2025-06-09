# app.py
# Importieren der notwendigen Bibliotheken
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import networkx as nx
from glycolysis_model import GlycolysisPathway # Importiert die Logik aus der anderen Datei

st.set_page_config(layout="wide") # Nutzt die volle Seitenbreite
st.title("🔬 Interaktive Glykolyse-Simulation")

st.markdown("""
Dieses Tool modelliert den Fluss von Metaboliten durch die Glykolyse. 
Über die Regler in der Seitenleiste können Sie die Anfangsbedingungen und die Enzymaktivitäten (Konzentrationen) 
der wichtigsten regulatorischen Schritte verändern, um deren Auswirkungen auf den Stoffwechselweg zu untersuchen.
""")

# --- Interaktive Steuerung in der Seitenleiste ---
st.sidebar.header("Simulations-Parameter")

# Globale Simulationsparameter
steps = st.sidebar.slider("Anzahl der Simulationsschritte", 10, 1000, 200, 10)
dt = st.sidebar.slider("Zeitschrittweite (dt)", 0.1, 5.0, 1.0, 0.1)

st.sidebar.header("Anfangsbedingungen & Engpässe")

# Regler für die Substrat-Verfügbarkeit
initial_glucose = st.sidebar.slider(
    "Anfangskonzentration Glukose (mM)", 
    min_value=1.0, max_value=20.0, value=10.0, step=0.5
)

# Regler für die Konzentration der Schlüsselenzyme, um Engpässe zu simulieren
st.sidebar.markdown("**Enzymkonzentrationen (mM)**")
hexokinase_conc = st.sidebar.slider("Hexokinase", 0.0005, 0.005, 0.0025, 0.0005, format="%.4f")
pfk_conc = st.sidebar.slider("Phosphofructokinase (PFK)", 0.0005, 0.005, 0.0020, 0.0005, format="%.4f")
pyruvate_kinase_conc = st.sidebar.slider("Pyruvat-Kinase", 0.0005, 0.005, 0.0020, 0.0005, format="%.4f")

# Button zum Starten der Simulation
if st.button("🚀 Simulation starten und analysieren"):

    # Dictionary mit den vom User gewählten Enzymkonzentrationen erstellen
    enzyme_configs = {
        "Hexokinase": hexokinase_conc,
        "PFK": pfk_conc,
        "Pyruvat-Kinase": pyruvate_kinase_conc
    }

    # Modell mit den gewählten Parametern initialisieren
    model = GlycolysisPathway(initial_glucose=initial_glucose, enzyme_configs=enzyme_configs)
    
    # Simulation durchführen
    data = model.simulate(steps=steps, dt=dt)

    st.header("Simulationsergebnisse")
    
    # Ergebnisse in zwei Spalten anzeigen
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Metaboliten-Konzentrationen über die Zeit")
        
        # Plot erstellen
        fig, ax = plt.subplots(figsize=(10, 7))
        
        # Nur die Haupt-Kohlenstoff-Metaboliten und ATP/ADP für Übersichtlichkeit plotten
        main_metabolites = [
            "Glukose", "G6P", "F1,6BP", "G3P", "PEP", "Pyruvat", "ATP", "ADP"
        ]
        for met, values in data.items():
            if met in main_metabolites:
                ax.plot(np.arange(len(values)) * dt, values, label=met)
        
        ax.set_xlabel(f"Zeit (s)")
        ax.set_ylabel("Konzentration (mM)")
        ax.set_title("Verlauf der Metaboliten-Konzentrationen")
        ax.legend(loc='upper right')
        ax.grid(True)
        plt.tight_layout()
        st.pyplot(fig)

    with col2:
        st.subheader("Visualisierung des Stoffwechselwegs")
        
        # Netzwerk-Graph erstellen und zeichnen
        G = model.create_pathway_graph()
        
        fig_net, ax_net = plt.subplots(figsize=(10, 10))
        
        # Layout für den Graphen
        pos = nx.spring_layout(G, k=0.9, iterations=50)

        # Knotenfarben basierend auf Typ (Metabolit oder Enzym)
        node_colors = ['#1f78b4' if G.nodes[n].get('type') != 'enzyme' else '#e31a1c' for n in G.nodes]

        nx.draw(G, pos, with_labels=True, node_size=2500, node_color=node_colors, 
                font_size=8, font_weight='bold', edge_color='gray', width=1.5, ax=ax_net)
        
        ax_net.set_title("Flussdiagramm der Glykolyse")
        st.pyplot(fig_net)

        st.markdown("""
        **Legende:**
        - <span style='color:#1f78b4'>**Blaue Knoten**</span>: Metaboliten
        - <span style='color:#e31a1c'>**Rote Knoten**</span>: Enzyme
        """, unsafe_allow_html=True)