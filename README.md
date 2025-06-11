# Glykolyse-Simulation

Simulierung des biochemischen Prozesses der Glykolyse mit objektorientierter Programmierung. Es wurden reale Enzymparameter verwendet die dem Stoffwechselweg, bei dem Glukose schrittweise durch die Enzyme zu Pyruvat abgebaut wird, folgen.

## Projektinhalte

- `Simulationsmodul`: Berechnet den Konzentrationsverlauf auf Basis der Enzymparameter

- `Webinterface`: Ermöglicht einfache Benutzerinteraktion und Visualisierung

- `Stoffwechselpfaddiagramm`: Veranschaulicht die metabolischen Schritte und Enzymkatalysen der Glykolyse


## Verwendete Bibliotheken
- `numpy` – für numerische Berechnungen  
- `streamlit` – für die Erstellung interaktiver Webanwendungen direkt in Python. Ermöglicht die einfache Integration von Slidern, Buttons, Diagrammen und Layouts ohne Frontend-Kenntnisse
- `matplotlib` – für die Darstellung der zeitlichen Konzentrationsverläufe der Metaboliten in Form von Liniengrafiken
- `graphviz` – für die Visualisierung des Glykolysewegs als gerichteter Graph, einschließlich reversibler und irreversibler Reaktionen
- `Glykolyse_1 (eigenes Modul)`- Enthält das Modell `GlycolysisPathway`, das die mathematische Simulation der Glykolyse übernimmt. Dieses Modell basiert auf Michaelis-Menten-Kinetik und führt die zeitliche Integration der Reaktionsschritte durch
### Struktur und Funktionsweise

#### Klassen

- **`Metabolite`**  
  → Modellierung einzelner Moleküle wie Glukose und Pyruvat mit ihrem Konzentrationsverlauf über die Zeit

- **`Enzyme`**  
  → Implementierung der Michaelis-Menten-Kinetik mit Werten aus der Literatur für kcat (Turnover number), Km (Michaelis-Menten-Konstante), und typischen Enzymkonzentrationen

- **`Reaction`**  
  → Modellierung der enzymatischen Umwandlungen eines Substrats zum jeweiligem Produkt

- **`SplitReaction`**  
  → Modellierung Spaltungsreaktion bei der ein Substrat in zwei Produkte reagiert (hier: Spaltung Fruktose-1,6-bisphosphat → Dihydroxyacetonphosphat + Glycerinaldehyd-3-phosphat) 

- **`GlycolysisPathway`**  
  → Modellierung aller Metabolite, Enzyme und Reaktionen des Glykolysewegs
  → Mit der Methode `simulate()` werden zeitlich aufgelöste Konzentrationsverläufe berechnet


## Beispiel: Simulation starten

Stellen Sie sicher, dass alle Dateien in dem selben Ordner installiert sind. Sie können das Programm mit dem folgenden Terminalbefehl starten:

```python
streamlit run Glyko_Streamlit.py