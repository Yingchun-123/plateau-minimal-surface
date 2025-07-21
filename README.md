
# 🧩 Minimalflächenberechnung (Plateau-Problem) – C++ Projekt

**Sprache:** C++  
**Tools:** FreeFem++, Automatische Differenzierung  
**Zeitraum:** 2014–2015  
**Studierende:** Yingchun SONG

## 🎯 Projektziel

Das Ziel dieses Projekts war die numerische Lösung des klassischen Plateau-Problems: 
Gesucht ist eine Fläche minimaler Fläche, deren Rand eine vorgegebene geschlossene Kurve im ℝ³ ist.  
Die Fläche wird dabei als Graph einer Funktion \( z = f(x, y) \) über einem 2D-Gebiet modelliert und diskretisiert.

## 🧠 Mathematischer Hintergrund

Die Funktional, die minimiert wird, ist:

\[
J(z) = \sum_{K \in \Sigma_h} \text{Fläche}(K)
\]

Dabei ist \( \Sigma_h \) eine Triangulation des Gebiets und \( K \) ein Dreieck mit Fläche \( |K| = \frac{1}{2} | \vec{AB} \wedge \vec{AC} | \)

## 🔧 Technische Umsetzung

- Definition und Implementierung der Funktionen:
  - `aire()`: Fläche eines Dreiecks im ℝ³
  - `daire()`: Erste Ableitungen der Fläche nach z-Werten
  - `ddaire()`: Zweite Ableitungen (Hesse-Matrix)
- Entwicklung einer `Mesh2`-Struktur zur Verwaltung von Knoten, Dreiecken und Randpunkten
- Einlesen von `.msh`-Dateien, Extraktion von Eckpunkten, Flächen, Kanten
- Implementierung von:
  - `J()`, `dJ()`, `ddJ()` zur Berechnung von Funktional, Gradient und Hesse
  - Gradientenabstieg
  - Newton-Verfahren
  - Konjugierter Gradientenalgorithmus für lineare Gleichungssysteme

## ⚙️ Verwendete Algorithmen

- **Gradientenverfahren**:
  - Iteratives Update \( z_{n+1} = z_n - \rho \nabla J(z_n) \)
- **Newton-Verfahren**:
  - Löst \( \nabla J(z) = 0 \) mit Hesse-Matrix
  - Verwendung des konjugierten Gradientenverfahrens zur Lösung von \( H u = b \)

## ✅ Validierung & Beispiele

Die Algorithmen wurden mit analytisch bekannten Minimalflächen getestet:

- **Caténoïde**: \( z = \pm \text{acosh}(\sqrt{x^2 + y^2}) \)
- **Scherk-Fläche**: \( z = \ln \cos(x) - \ln \cos(y) \)
- **Helikoide**: \( f(x,y) = \tan(a y / x) \)

Vergleiche zwischen numerischer Lösung \( u_h \) und exakter Lösung \( f(x, y) \) wurden durchgeführt.

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

- `catenoid_example.png`
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
