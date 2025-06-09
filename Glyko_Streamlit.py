import streamlit as st
import matplotlib.pyplot as plt
import graphviz  
from Glykolyse_1 import GlycolysisPathway

st.set_page_config(page_title="Glykolyse-Simulation", layout="wide")
st.title("Glykolyse-Simulation")

# Simulationsparameter
st.sidebar.header("Simulation")
steps = st.sidebar.slider("Anzahl Schritte", 10, 1000, 100, step=10)
dt = st.sidebar.number_input("Zeitschritt (dt)", 0.01, 10.0, 1.0, step=0.1)
glucose_input = st.slider("Glukose-Konzentration (mmol/L)", min_value=0.0, max_value=10.0, step=0.1)

# Modell laden
model = GlycolysisPathway(glucose_conc=glucose_input)

# Mapping von Enzymnamen zu Instanzen
enzym_map = {
    "Hexokinase": model.hexokinase,
    "Glukose-6-phosphat-Isomerase": model.isomerase,
    "Phosphofructokinase": model.pfk,
    "Aldolase": model.aldolase,
    "Triosephosphat-Isomerase": model.triosephosphat_isomerase,
    "Glycerinaldehyd-3-phosphat-Dehydrogenase": model.gapdh,
    "Phosphoglyceratkinase": model.pgk,
    "Phosphoglyceratmutase": model.pgm,
    "Enolase": model.enolase,
    "Pyruvat-Kinase": model.pyruvate_kinase,
}

# Sidebar: Enzymparameter
st.sidebar.header("Enzymparameter")
for name, enzyme in enzym_map.items():
    with st.sidebar.expander(f"{name}", expanded=False):
        enzyme.kcat = st.slider(f"{name} kcat (1/s)", 50, 5000, enzyme.kcat, 10)
        enzyme.enzyme_conc = st.number_input(f"{name} – [E] (mM)", 0.0001, 0.01, enzyme.enzyme_conc, 0.0001)
        enzyme.km = st.number_input(f"{name} Km (mM)", 0.01, 1.0, enzyme.km, 0.01)
        enzyme.vmax = enzyme.kcat * enzyme.enzyme_conc  # vmax aktualisieren

# Simulation starten
if st.button("Simulation starten"):
    data = model.simulate(steps=steps, dt=dt)
    st.success("Simulation abgeschlossen!")

    # Plot
    st.subheader("Konzentrationsverläufe der Metaboliten")
    fig, ax = plt.subplots(figsize=(10, 6))
    for met, values in data.items():
        ax.plot(values, label=met)
    ax.set_xlabel("Zeit (s)")
    ax.set_ylabel("Konzentration (mM)")
    ax.legend()
    ax.grid(True)
    st.pyplot(fig)



st.title("Glykolyse: Fluss der Metaboliten")

st.markdown("""
Diese Visualisierung zeigt den Fluss der Metaboliten durch die Glykolyse. 
Zwischen den Metaboliten sind die katalysierenden Enzyme angegeben.

Reversible Reaktionen = Doppelpfeil  
Irreversible Reaktionen = Einfacher Pfeil
""")

# Graphviz Diagramm
dot = graphviz.Digraph()

# Metaboliten
metabolites = [
    "Glucose",
    "Glucose-6-phosphat",
    "Fructose-6-phosphat",
    "Fructose-1,6-bisphosphat",
    "Dihydroxyacetonphosphat",
    "Glycerinaldehyd-3-phosphat",
    "1,3-Bisphosphoglycerat",
    "3-Phosphoglycerat",
    "2-Phosphoglycerat",
    "Phosphoenolpyruvat",
    "Pyruvat"
]

# Reaktionen: (Start, Ende, Enzym, reversible=True/False)
steps = [
    ("Glucose", "Glucose-6-phosphat", "Hexokinase", False),
    ("Glucose-6-phosphat", "Fructose-6-phosphat", "Glucose-6-phosphat-Isomerase", True),
    ("Fructose-6-phosphat", "Fructose-1,6-bisphosphat", "Phosphofructokinase", False),
    ("Fructose-1,6-bisphosphat", "Dihydroxyacetonphosphat", "Aldolase", True),
    ("Fructose-1,6-bisphosphat", "Glycerinaldehyd-3-phosphat", "Aldolase", True),
    ("Dihydroxyacetonphosphat", "Glycerinaldehyd-3-phosphat", "Triosephosphat-Isomerase", True),
    ("Glycerinaldehyd-3-phosphat", "1,3-Bisphosphoglycerat", "Glycerinaldehyd-3-phosphat-Dehydrogenase", True),
    ("1,3-Bisphosphoglycerat", "3-Phosphoglycerat", "Phosphoglycerat-Kinase", False),
    ("3-Phosphoglycerat", "2-Phosphoglycerat", "Phosphoglycerat-Mutase", True),
    ("2-Phosphoglycerat", "Phosphoenolpyruvat", "Enolase", True),
    ("Phosphoenolpyruvat", "Pyruvat", "Pyruvatkinase", False)
]

# Knoten hinzufügen
for m in metabolites:
    dot.node(m)

# Schritte mit Pfeilrichtung je nach Reversibilität
for start, end, enzyme, reversible in steps:
    direction = "both" if reversible else "forward"
    dot.edge(start, end, label=enzyme, dir=direction)

# Diagramm anzeigen
st.graphviz_chart(dot)


