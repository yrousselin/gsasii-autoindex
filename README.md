# Auto-indexing (unit-cell search) with several methods for GSAS-II

*[Version française : LISEZMOI.md](LISEZMOI.md)*

A plugin for **GSAS-II** that looks for the **unit cell of a powder pattern** by running several indexing methods side by side. It compares their solutions with the same figures of merit, and sends the cells you choose to the *Unit Cells List* of your project.

It is the first step of a simple workflow:
**Auto-index** (this plugin) → [**gsasii-cellsearch**](https://github.com/yrousselin/gsasii-cellsearch) (is this cell already known?) → [**gsasii-phasematch**](https://github.com/yrousselin/gsasii-phasematch) (which phases are in the sample?).

## What it does

- Adds the menu entry **Export → Entire project as → (user) Auto-index (find unit cell)** to GSAS-II.
- Reads the peaks of the chosen histogram: the *Index Peak List* if it is filled, otherwise the *Peak List*.
- Runs the methods you tick **in parallel**:

| Method | What it is | Crystal systems |
|---|---|---|
| **Conograph** | external indexing program (recommended) | all, including triclinic |
| **Dichotomy** | Python re-implementation of the DICVOL principle (successive volume dichotomy) | cubic to monoclinic |
| **Trial-and-error** | Python re-implementation of the TREOR principle (trials on low-angle lines) | cubic to monoclinic |
| **GSAS-II native** | the indexing engine built into GSAS-II | all (slow in triclinic) |

- Refines each candidate cell on the peaks, then recomputes **M20 and F20 the same way for every method**. Their solutions can therefore be compared fairly.
- Shows a report with a comparison of the methods and an overall ranking. You tick the cells to keep.
- Writes the ticked cells into the **Unit Cells List** of the open project. Nothing is saved to disk until you save the project.
- The calculation runs in a separate process: GSAS-II stays usable, and a **Cancel** button stops the search.

**The interface is in French.** The dichotomy and trial-and-error engines are re-implementations of the principles of DICVOL and TREOR, **not the original programs**.

## Requirements

| What | Why | Where to get it |
|---|---|---|
| **GSAS-II** | hosts the plugin; provides numpy, scipy and wxPython | [GSAS-II](https://advancedphotonsource.github.io/GSAS-II-tutorials/install.html) |
| **Conograph** *(recommended)* | the most robust method, the only one for triclinic cells besides GSAS-II | [Z-Code downloads](https://z-code-software.com/downloads/) |
| A peak list | the input of every method | *Peak List* or *Index Peak List* in GSAS-II |

Without Conograph, the three other methods still work.

**Tip:** peaks fitted in GSAS-II (*Peak List → Peak fitting*) give more precise cells than peaks that are only picked.

## Installation

### Step 1: copy the files into GSAS-II

**macOS (and Linux)**
1. Download this repository: green **Code** button → **Download ZIP**, then double-click the ZIP to unpack it.
2. Open the **Terminal**. Type `bash `, with a space after it, drag `install_autoindex_plugin.sh` into the window, and press Enter. The two files are copied into `~/.GSASII/exports/`.

**Windows**
1. Download the ZIP as above, then right-click it → **Extract All**.
2. Double-click **`install_windows.bat`**. The two files are copied into `C:\Users\<your name>\.GSASII\exports\`.

### Step 2: install Conograph (recommended)

1. Download Conograph from the [Z-Code downloads page](https://z-code-software.com/downloads/) and unpack it in a folder of your choice. It can also be compiled from its source code.
2. Note the full path of the **Conograph** program (`Conograph.exe` on Windows).
3. The plugin finds it automatically if it is on the system path or if the `CONOGRAPH` environment variable points to it. Otherwise, type its path in the **Exécutable Conograph** field of the plugin window the first time. It is remembered afterwards.

### Step 3: restart GSAS-II

## Usage

1. Open a GSAS-II project with a powder pattern, and fill the *Peak List* (ideally with fitted peaks).
2. Open **Export → Entire project as → (user) Auto-index (find unit cell)**.
3. In the window, choose the histogram and the methods. Adjust the settings if needed:
   - maximum number of peaks used;
   - zero shift;
   - tolerance on 2θ;
   - maximum duration of each search;
   - quick or extended (*poussée*) Conograph search.
4. Click OK. The methods run in parallel; you can cancel at any time.
5. Read the report, tick the cells to keep, and confirm. They appear in the **Unit Cells List** of the histogram.
6. Save the project.

**From the command line** (with GSAS-II's Python). The cells are then written to a **copy** of the project:

```
python g2_autoindex.py project.gpx                       # asks which methods to run
python g2_autoindex.py project.gpx --methods c,d,t       # Conograph, dichotomy, trial-and-error
python g2_autoindex.py project.gpx --methods c --maxpeaks 25
python g2_autoindex.py --help                            # all options
```

## Limitations

- The dichotomy and trial-and-error engines do not handle triclinic cells. Use Conograph or GSAS-II for those.
- Indexing depends on the quality of the peak positions. Impurity lines, or too few lines, can give wrong or no solutions. Always check the proposed cell, for example by a Le Bail fit.
- Developed and used on macOS. Windows should work, but has not been tested yet. Feedback is welcome.

## Credits and references

- **GSAS-II:** B. H. Toby, R. B. Von Dreele, *J. Appl. Cryst.* **2013**, 46, 544–549. The native method uses the GSAS-II indexing module.
- **Conograph** (powder indexing): « New functions and graphical user interface attached to powder indexing software », *J. Appl. Cryst.* **2017**, 50, 651–659.
- **DICVOL** (dichotomy method): D. Louër, M. Louër, *J. Appl. Cryst.* **1972**, 5, 271–275; A. Boultif, D. Louër, *J. Appl. Cryst.* **2004**, 37, 724–731.
- **TREOR** (trial-and-error method): P.-E. Werner, L. Eriksson, M. Westdahl, *J. Appl. Cryst.* **1985**, 18, 367–370.
- **Figures of merit:** M20 (de Wolff, *J. Appl. Cryst.* **1968**, 1, 108–113) and F20 (Smith & Snyder, *J. Appl. Cryst.* **1979**, 12, 60–65).
- Written by Yoann Rousselin with the help of Claude (Anthropic) as a coding assistant.

## Licence

MIT. See [LICENSE](LICENSE). Conograph and GSAS-II are distributed under their own licences.

## Author

**Yoann Rousselin**, CNRS research engineer (X-ray crystallography), ICMUB, Dijon (France).
Questions or bug reports: please open an *Issue* on this page.
