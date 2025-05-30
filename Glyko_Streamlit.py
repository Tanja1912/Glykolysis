
import streamlit as st
import matplotlib.pyplot as plt
from Glykolyse_1 import GlycolysisPathway

st.title("Glykolyse-Simulation")

# Eingabeparameter
steps = st.slider("Anzahl der Simulationsschritte", min_value=10, max_value=1000, value=100, step=10)
dt = st.number_input("Zeitschrittweite (dt)", min_value=0.01, max_value=10.0, value=1.0, step=0.1)

if st.button("Simulation starten"):
    model = GlycolysisPathway()
    data = model.simulate(steps=steps, dt=dt)

    st.write(f"Simulation mit {steps} Schritten und dt = {dt}")

    # Plotten mit matplotlib und Streamlit
    fig, ax = plt.subplots(figsize=(10, 6))
    for met, values in data.items():
        ax.plot(values, label=met)
    ax.set_xlabel("Zeit (a.u.)")
    ax.set_ylabel("Konzentration (mM)")
    ax.set_title("Glykolyse: Metabolitenkonzentrationen")
    ax.legend()
    ax.grid(True)
    st.pyplot(fig)



