# Minimalflächenberechnung (Plateau-Problem) – C++ Projekt

**Sprache:** C++  
**Tools:** FreeFem++, automatische Differenzierung  
**Zeitraum:** 2014–2015  
**Studierende:** Yingchun SONG  

---

## 🎯 Projektziel

Ziel des Projekts ist die numerische Lösung des klassischen Plateau-Problems.  
Dabei wird eine Fläche minimaler Fläche gesucht, deren Rand eine vorgegebene geschlossene Kurve im dreidimensionalen Raum bildet.

Die gesuchte Fläche wird als Graph einer Funktion z = f(x, y) über einem zweidimensionalen Gebiet modelliert, diskretisiert und anschließend numerisch minimiert.

---

## 🧠 Mathematischer Hintergrund

Das Minimierungsproblem basiert auf der Summe der Flächen aller Dreieckselemente einer Triangulation.  
Dabei gilt:

- Die Oberfläche wird durch eine Triangulation des Gebiets approximiert.  
- Jedes Dreieck liefert einen positiven Flächenbeitrag.  
- Die Gesamtfläche ergibt sich aus der Summe aller Dreiecksflächen.

Die Fläche eines einzelnen Dreiecks wird über das Kreuzprodukt der Kantenvektoren berechnet (halbe Norm des Kreuzproduktes).  
Dieses Prinzip dient als Grundlage für das gesamte Minimierungsverfahren.

---

## ⚙️ Numerische Umsetzung

Die Implementierung umfasst folgende Schritte:

### **1. Triangulation**
- Zerlegung des Gebiets in Dreiecke (Finite-Elemente-Struktur).
- Definition der Randkurve.

### **2. Formulierung des Minimierungsproblems**
- Berechnung der Dreiecksflächen auf Basis der aktuellen Funktion z(x, y).
- Aufstellen des Gesamtflächenfunktionals.

### **3. Gradientenverfahren**
- Iteratives Update des Oberflächenprofils:
  - z_{n+1} = z_n – Schrittweite * Gradient des Flächenfunktionals  
- Ziel ist die Verringerung der Gesamtfläche bei jedem Schritt.

### **4. Newton-Verfahren**
- Lösung des Gleichungssystems, das aus der Bedingung „Gradient gleich Null“ entsteht.
- Verwendung der Hesse-Matrix (zweite Ableitungen) für schnellere Konvergenz.
- Für die linearen Gleichungssysteme wird das konjugierte Gradientenverfahren eingesetzt.

### **5. Visualisierung**
- Darstellung der triangulierten Minimalfläche.
- Vergleich verschiedener Iterationsstufen zur Analyse der Konvergenz.

## 🧪 Ergebnisse (Beispiele)

| Mesh         | Methode           | Iterationen | Zeit (s) | Fehler |
|--------------|-------------------|-------------|----------|--------|
| `c10.msh`    | Gradient          | 435         | 0        | ~1e-8  |
|              | Newton            | 3           | 0        | ~1e-8  |
| `c50.msh`    | Gradient          | 4254        | 9        | ~1e-5  |
|              | Newton            | >10800      | >?       | ~1e-5  |
| `carre.msh`  | Gradient          | 1           | 0        | ~1e-17 |
| `cercle.msh` | Newton            | 6           | 2        | ~1e-9  |

## 📷 Visualisierung

Bitte die Bilddateien im Ordner `/images` bereitstellen:

- `catenoid_example.jpeg`
- `scherk_example.png`
- `helicoid_example.png`
- `result_comparison.png`

## 📁 Projektstruktur

```
plateau-minimal-surface/
├── src/
├── meshes/
├── results/
├── images/
└── README.md
```

## 🧠 Gelernt

- Numerische Optimierung und Differenzierung
- Umsetzung komplexer Gleichungen in C++
- Anwendung mathematischer Modelle auf physikalisch motivierte Probleme
