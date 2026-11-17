#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
g2_autoindex.py (v6) - Recherche de maille (indexation) multi-methodes a partir d'un projet GSAS-II (.gpx)

Methodes, au choix (option --methods, ou question posee au lancement) :
  c  conograph   : programme Conograph (executable externe, licence MIT), tous systemes, tres robuste
  d  dichotomie  : reecriture en Python du principe de DICVOL (dichotomie par tranches de volume)
  t  treor       : reecriture en Python du principe de TREOR (essai-erreur sur les raies de bas angle)
  g  gsas        : moteur natif de GSAS-II (GSASIIindex, type SVD-Index) - lent, surtout en triclinique
Les moteurs "dichotomie" et "treor" traitent du cubique au monoclinique (pas le triclinique) et ne
sont pas les programmes originaux DICVOL / TREOR, mais des reecritures de leur principe.

Le script lit les pics du projet (Index Peak List si elle est remplie, sinon Peak List), lance les
methodes choisies en parallele, affine chaque maille candidate sur les pics, recalcule M20 / F20 de
facon uniforme, compare les methodes entre elles, puis vous demande quelles mailles transmettre a
GSAS-II (ecrites dans une COPIE du projet, dans la "Unit Cells List").

Usage :
    python g2_autoindex.py projet.gpx                        (demande les methodes)
    python g2_autoindex.py projet.gpx --methods c,d,t        (tout sauf GSAS-II)
    python g2_autoindex.py projet.gpx --methods c --maxpeaks 25

Dependances : numpy, scipy. GSAS-II n'est necessaire que pour la methode "gsas" et pour l'ecriture
des mailles dans le projet.
Conseil : des pics ajustes dans GSAS-II (Peak List > Peak fitting) donnent des mailles plus precises
que des pics simplement pointes.
"""
import argparse, contextlib, io, json, multiprocessing, os, pickle, re, shutil, subprocess, sys, tempfile, time
import numpy as np

# --------------------------------------------------------------------------
# 1. Lecture du .gpx (sans importer GSAS-II)
# --------------------------------------------------------------------------
class _Stub(object):
    def __init__(self, *a, **k): pass
    def __setstate__(self, s): self.__dict__['state'] = s

_SAFE = ('numpy', 'builtins', 'collections', 'copyreg', '_codecs')
_BAD = ('eval', 'exec', 'open', '__import__', 'getattr', 'compile', 'system')

class _Unpickler(pickle.Unpickler):
    """N'instancie que des objets numpy / types de base ; tout le reste devient un objet inerte."""
    def find_class(self, module, name):
        if module.split('.')[0] in _SAFE and name not in _BAD:
            return pickle.Unpickler.find_class(self, module, name)
        return type(str(name), (_Stub,), {})

def read_gpx(path):
    recs = []
    with open(path, 'rb') as f:
        while True:
            try:
                recs.append(_Unpickler(f, encoding='latin-1').load())
            except EOFError:
                break
    hists, phases = [], []
    for rec in recs:
        name = rec[0][0]
        if name.startswith('PWDR'):
            sub = dict((s[0], s[1]) for s in rec[1:])
            data = rec[0][1][1]
            hists.append(dict(name=name, x=np.array(data[0], float), y=np.array(data[1], float),
                              limits=sub.get('Limits'), inst=sub['Instrument Parameters'][0],
                              peaks=sub.get('Peak List', {}), ipeaks=sub.get('Index Peak List', [[], []])))
        elif name == 'Phases':
            for s in rec[1:]:
                g = s[1]['General']
                phases.append(dict(name=s[0], sg=g['SGData']['SpGrp'], cell=[float(v) for v in g['Cell'][1:8]]))
    return hists, phases

def inst_val(inst, key, default=None):
    v = inst.get(key)
    return default if v is None else v[1]

def get_wavelength(inst):
    lam = inst_val(inst, 'Lam')
    if lam is None:
        lam = inst_val(inst, 'Lam1')
    if lam is None:
        raise ValueError("Longueur d'onde introuvable dans les parametres instrumentaux")
    return float(lam)

def fwhm_from_inst(inst, tth):
    """Largeur a mi-hauteur (deg) d'apres U,V,W,X,Y de GSAS-II (centidegres)."""
    th = np.radians(np.asarray(tth, float) / 2.)
    U, V, W = [float(inst_val(inst, k, 0.)) for k in 'UVW']
    X, Y = [float(inst_val(inst, k, 0.)) for k in 'XY']
    sig2 = np.maximum(U * np.tan(th) ** 2 + V * np.tan(th) + W, 0.5)      # centideg^2
    fg = np.sqrt(8. * np.log(2.) * sig2)
    fl = np.maximum(X / np.cos(th) + Y * np.tan(th), 0.)
    f = (fg ** 5 + 2.69269 * fg ** 4 * fl + 2.42843 * fg ** 3 * fl ** 2 + 4.47163 * fg ** 2 * fl ** 3
         + 0.07842 * fg * fl ** 4 + fl ** 5) ** 0.2
    return np.clip(f / 100., 0.02, 1.0)

def get_peaks(hist):
    """Retourne (2theta, intensite, origine). Priorite a l'Index Peak List si elle contient des pics utilises."""
    ip = hist['ipeaks'][0] if hist['ipeaks'] else []
    rows = [p for p in ip if len(p) > 2 and p[2]]
    if len(rows) >= 5:
        return np.array([p[0] for p in rows], float), np.array([p[1] for p in rows], float), 'Index Peak List'
    pk = hist['peaks'].get('peaks', []) if isinstance(hist['peaks'], dict) else hist['peaks']
    if len(pk) < 5:
        raise ValueError("Moins de 5 pics dans le projet : faites d'abord une recherche de pics dans GSAS-II")
    return np.array([p[0] for p in pk], float), np.array([p[2] for p in pk], float), 'Peak List'

# --------------------------------------------------------------------------
# 2. Outils cristallographiques : positions calculees, M20, F20
# --------------------------------------------------------------------------
def cell_to_recip_metric(cell):
    a, b, c, al, be, ga = cell[:6]
    ca, cb, cg = [np.cos(np.radians(v)) for v in (al, be, ga)]
    G = np.array([[a * a, a * b * cg, a * c * cb], [a * b * cg, b * b, b * c * ca], [a * c * cb, b * c * ca, c * c]])
    return np.linalg.inv(G), np.sqrt(abs(np.linalg.det(G)))

_CENT = {'P': lambda h, k, l: np.ones(h.shape, bool),
         'I': lambda h, k, l: (h + k + l) % 2 == 0,
         'F': lambda h, k, l: ((h + k) % 2 == 0) & ((k + l) % 2 == 0),
         'A': lambda h, k, l: (k + l) % 2 == 0,
         'B': lambda h, k, l: (h + l) % 2 == 0,
         'C': lambda h, k, l: (h + k) % 2 == 0,
         'R': lambda h, k, l: (-h + k + l) % 3 == 0}

def calc_q(cell, qmax, cent='P'):
    """Liste triee des Q = 1/d^2 distincts (a 1e-7 pres) jusqu'a qmax."""
    Gs, vol = cell_to_recip_metric(cell)
    a, b, c = cell[:3]
    dmin = 1. / np.sqrt(qmax)
    hm, km, lm = [int(np.ceil(v / dmin)) + 1 for v in (a, b, c)]
    h, k, l = np.meshgrid(np.arange(0, hm + 1), np.arange(-km, km + 1), np.arange(-lm, lm + 1), indexing='ij')
    h, k, l = h.ravel(), k.ravel(), l.ravel()
    keep = (h > 0) | ((h == 0) & (k > 0)) | ((h == 0) & (k == 0) & (l > 0))
    keep &= _CENT.get(cent, _CENT['P'])(h, k, l)
    h, k, l = h[keep], k[keep], l[keep]
    q = (Gs[0, 0] * h * h + Gs[1, 1] * k * k + Gs[2, 2] * l * l
         + 2 * (Gs[0, 1] * h * k + Gs[0, 2] * h * l + Gs[1, 2] * k * l))
    q = np.sort(q[q <= qmax * 1.0001])
    if len(q) == 0:
        return q, vol
    uniq = q[np.concatenate(([True], np.diff(q) > 1e-7 * max(q[-1], 1e-9)))]
    return uniq, vol

def unindexed_peaks(cell, cent, tth_obs, lam, zero=0., tol_deg=0.05):
    """Positions observees (2theta) sans raie calculee a moins de tol_deg."""
    tth = np.sort(np.asarray(tth_obs, float))
    qmax = (2. * np.sin(np.radians(tth[-1] / 2.)) / lam) ** 2 * 1.03
    qc, _ = calc_q(cell, qmax, cent)
    tc = 2. * np.degrees(np.arcsin(np.clip(lam * np.sqrt(qc) / 2., 0, 1))) + zero
    i = np.clip(np.searchsorted(tc, tth), 1, len(tc) - 1)
    d = np.minimum(np.abs(tth - tc[i - 1]), np.abs(tth - tc[i]))
    return [(float(t), float(x)) for t, x in zip(tth, d) if x > tol_deg]

def figures_of_merit(cell, cent, tth_obs, lam, zero=0., nfom=20, tol_deg=0.05):
    """M20 (de Wolff), F20 (Smith & Snyder), nombre de pics non indexes parmi les nfom premiers."""
    tth = np.sort(np.asarray(tth_obs, float))[:nfom] - zero
    n = len(tth)
    qobs = (2. * np.sin(np.radians(tth / 2.)) / lam) ** 2
    qmax = qobs[-1]
    qc, vol = calc_q(cell, qmax * 1.02, cent)
    if len(qc) == 0:
        return dict(M=0., F=0., nunidx=n, ncalc=0, vol=vol, d2th=9.)
    idx = np.clip(np.searchsorted(qc, qobs), 1, len(qc) - 1) if len(qc) > 1 else np.zeros(n, int)
    lo = qc[np.clip(idx - 1, 0, len(qc) - 1)]; hi = qc[np.clip(idx, 0, len(qc) - 1)]
    qnear = np.where(np.abs(qobs - lo) < np.abs(qobs - hi), lo, hi)
    tthc = 2. * np.degrees(np.arcsin(np.clip(lam * np.sqrt(qnear) / 2., 0, 1)))
    d2 = np.abs(tth - tthc)
    nunidx = int(np.sum(d2 > tol_deg))
    ncalc = max(int(np.sum(qc <= qmax * (1 + 1e-9))), 1)
    dq = np.mean(np.abs(qobs - qnear))
    M = qmax / (2. * max(dq, 1e-12) * ncalc)
    F = n / (max(np.mean(d2), 1e-6) * ncalc)
    return dict(M=M, F=F, nunidx=nunidx, ncalc=ncalc, vol=vol, d2th=float(np.mean(d2)))

# --------------------------------------------------------------------------
# 3. Moteur Conograph
# --------------------------------------------------------------------------
_CONO_BRAVAIS = [  # (nom dans le XML, systeme, centrage, balise de sortie)
    ('Cubic(F)', 'cubique', 'F'), ('Cubic(I)', 'cubique', 'I'), ('Cubic(P)', 'cubique', 'P'),
    ('Hexagonal', 'hexagonal', 'P'), ('Rhombohedral', 'rhomboedrique', 'R'),
    ('Tetragonal(I)', 'quadratique', 'I'), ('Tetragonal(P)', 'quadratique', 'P'),
    ('Orthorhombic(F)', 'orthorhombique', 'F'), ('Orthorhombic(I)', 'orthorhombique', 'I'),
    ('Orthorhombic(C)', 'orthorhombique', 'C'), ('Orthorhombic(P)', 'orthorhombique', 'P'),
    ('Monoclinic(C)', 'monoclinique', 'C'), ('Monoclinic(P)', 'monoclinique', 'P'),
    ('Triclinic', 'triclinique', 'P')]

def find_conograph(user=None):
    cands = [user] if user else []
    cands += [os.environ.get('CONOGRAPH'), shutil.which('Conograph'), shutil.which('Conograph.exe'),
              os.path.expanduser('~/ProjectPowderConograph/Conograph/src/Conograph'),
              os.path.expanduser('~/ProjectPowderConograph/Conograph/src/Conograph.exe'),
              shutil.which('Conograph CUI.exe'),                                   # version Windows de Z-Code
              r'C:\Program Files (x86)\Conograph\Conograph CUI.exe',
              r'C:\Program Files\Conograph\Conograph CUI.exe']
    for c in cands:
        if c and os.path.isfile(os.path.expanduser(c)) and os.access(os.path.expanduser(c), os.X_OK):
            return os.path.abspath(os.path.expanduser(c))
    return None

def write_conograph_inputs(wd, hist, tth, inten, fwhm, lam, zero, search, npeaks, nthreads):
    with open(os.path.join(wd, 'cntl.inp.xml'), 'w') as f:
        f.write("<ZCodeParameters>\n  <ConographInputFile>\n"
                "    <ControlParamFile> g2.inp.xml </ControlParamFile>\n"
                "    <PeakDataFile> g2_pks.histogramIgor </PeakDataFile>\n"
                "    <OutputFile> output/g2.index.xml </OutputFile>\n"
                "  </ConographInputFile>\n</ZCodeParameters>\n")
    os.makedirs(os.path.join(wd, 'output'), exist_ok=True)
    outflags = ''.join('    <Output%s> 1 </Output%s>\n' % (t, t) for t in (
        'Triclinic', 'MonoclinicP', 'MonoclinicB', 'OrthorhombicP', 'OrthorhombicC', 'OrthorhombicI',
        'OrthorhombicF', 'TetragonalP', 'TetragonalI', 'Rhombohedral', 'Hexagonal', 'CubicP', 'CubicI', 'CubicF'))
    with open(os.path.join(wd, 'g2.inp.xml'), 'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8" ?>\n<ZCodeParameters>\n  <ConographParameters>\n'
                '    <IsAngleDispersion> 1 </IsAngleDispersion>\n'
                '    <ConversionParameters>       </ConversionParameters>\n'
                '    <ZeroPointShiftParameter> %.5f </ZeroPointShiftParameter>\n'
                '    <WaveLength> %.6f </WaveLength>\n'
                '    <SearchLevel> %d </SearchLevel>\n'
                '    <MaxNumberOfPeaks> %s </MaxNumberOfPeaks>\n'
                '    <CriticalValueForLinearSum> 1.0 </CriticalValueForLinearSum>\n'
                '    <MinPrimitiveUnitCellVolume> AUTO </MinPrimitiveUnitCellVolume>\n'
                '    <MaxPrimitiveUnitCellVolume> AUTO </MaxPrimitiveUnitCellVolume>\n'
                '    <MaxNumberOfTwoDimTopographs> AUTO </MaxNumberOfTwoDimTopographs>\n'
                '    <MaxNumberOfLatticeCandidates> AUTO </MaxNumberOfLatticeCandidates>\n'
                '%s'
                '    <Resolution> 0.03 </Resolution>\n'
                '    <MaxNumberOfPeaksForFOM> 20 </MaxNumberOfPeaksForFOM>\n'
                '    <MinFOM> 3.0 </MinFOM>\n'
                '    <MaxNumberOfMillerIndicesInRange> AUTO </MaxNumberOfMillerIndicesInRange>\n'
                '    <MinNumberOfMillerIndicesInRange> AUTO </MinNumberOfMillerIndicesInRange>\n'
                '    <MaxUnitCellEdgeABC> 1000.0 </MaxUnitCellEdgeABC>\n'
                '    <MinUnitCellEdgeABC> 0.0 </MinUnitCellEdgeABC>\n'
                '    <NumberOfThreadsToUse> %s </NumberOfThreadsToUse>\n'
                '    <AxisForRhombohedralSymmetry> Hexagonal </AxisForRhombohedralSymmetry>\n'
                '    <AxisForMonoclinicSymmetry> B </AxisForMonoclinicSymmetry>\n'
                '    <ThresholdOnNormM> 1.9 </ThresholdOnNormM>\n'
                '    <ThresholdOnRevM> 1.0 </ThresholdOnRevM>\n'
                '    <MinDistanceBetweenLatticePoints> 2.0 </MinDistanceBetweenLatticePoints>\n'
                '  </ConographParameters>\n</ZCodeParameters>\n'
                % (zero, lam, 1 if search == 'regular' else 0, ('%d' % npeaks) if npeaks else 'AUTO',
                   outflags, nthreads))
    x, y = hist['x'], hist['y']
    with open(os.path.join(wd, 'g2_pks.histogramIgor'), 'w') as f:
        f.write('IGOR\nWAVES/O wave0, wave1, sqrt_of_wave1\nBEGIN\n')
        for a, b in zip(x, y):
            f.write('%12.5f %15.6E %15.6E\n' % (a, b, np.sqrt(max(b, 1.))))
        f.write('END\nWAVES/O peak, peakpos, height, FWHM, Flag\nBEGIN\n')
        for i, (t, h, w) in enumerate(zip(tth, inten, fwhm)):
            f.write('%5d %15.6E %15.6E %15.6E %4d\n' % (i + 1, t, max(h, 1.), w, 1))
        f.write('END\n')

_NUM = r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?'
_SYS = dict((n, (s, c)) for n, s, c in _CONO_BRAVAIS)

def parse_conograph_xml(path):
    """Extrait toutes les mailles candidates du fichier *.index.xml de Conograph."""
    txt = open(path, errors='replace').read()
    cands, seen = [], set()
    for num, blk in re.findall(r'<LatticeCandidate number="\s*(\d+)\s*">(.*?)</LatticeCandidate>', txt, flags=re.S):
        if num in seen:
            continue
        seen.add(num)
        def tag(t):
            m = re.search(r'<%s\b[^>]*>(.*?)</%s>' % (t, t), blk, flags=re.S)
            return m.group(1).strip() if m else None
        def num1(t):
            v = tag(t)
            return float(re.findall(_NUM, v)[0]) if v else float('nan')
        cell = tag('OptimizedLatticeParameters') or tag('ReducedLatticeParameters')
        vals = [float(v) for v in re.findall(_NUM, cell or '')][:6]
        if len(vals) < 6:
            continue
        brav = tag('CrystalSystem') or 'Triclinic'
        system, cent = _SYS.get(brav, ('triclinique', 'P'))
        cands.append(dict(engine='Conograph', label=num, cell=vals, bravais=brav, system=system, cent=cent,
                          M_cono=num1('FigureOfMeritWolff'), Mrev=num1('ReversedFigureOfMeritWolff'),
                          Msym=num1('SymmetricFigureOfMeritWolff')))
    return cands

def run_conograph(exe, hist, tth, inten, lam, zero=0., search='quick', npeaks=None, nthreads='MAX',
                  timeout=900, keep=None, verbose=True):
    wd = keep or tempfile.mkdtemp(prefix='g2cono_')
    os.makedirs(wd, exist_ok=True)
    fw = fwhm_from_inst(hist['inst'], tth)
    write_conograph_inputs(wd, hist, tth, inten, fw, lam, zero, search, npeaks, nthreads)
    t0 = time.time()
    try:
        # "quit" est indispensable : sans lui le programme boucle et remplit le disque
        p = subprocess.run([exe], cwd=wd, input='quit\n' * 3, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, universal_newlines=True, timeout=timeout)
        log = p.stdout[-4000:]
    except subprocess.TimeoutExpired:
        log = 'Conograph arrete apres %d s (timeout)' % timeout
    dt = time.time() - t0
    try:
        os.remove(os.path.join(wd, 'LOG_CONOGRAPH.txt'))
    except OSError:
        pass
    out = os.path.join(wd, 'output', 'g2.index.xml')
    if not os.path.isfile(out):
        raise RuntimeError("Conograph n'a pas produit de resultat.\n" + log)
    return out, dt, wd

# --------------------------------------------------------------------------
# 3b. Moteur natif GSAS-II (GSASIIindex.DoIndexPeaks), un processus par reseau de Bravais
# --------------------------------------------------------------------------
_G2_BRAVAIS = [  # ordre impose par GSAS-II : (nom, systeme, centrage, duree max relative)
    ('Cubic-F', 'cubique', 'F', 0.2), ('Cubic-I', 'cubique', 'I', 0.2), ('Cubic-P', 'cubique', 'P', 0.2),
    ('Trigonal-R', 'rhomboedrique', 'R', 0.3), ('Hexagonal-P', 'hexagonal', 'P', 0.3),
    ('Tetragonal-I', 'quadratique', 'I', 0.3), ('Tetragonal-P', 'quadratique', 'P', 0.3),
    ('Orthorhombic-F', 'orthorhombique', 'F', 0.4), ('Orthorhombic-I', 'orthorhombique', 'I', 0.4),
    ('Orthorhombic-A', 'orthorhombique', 'A', 0.4), ('Orthorhombic-B', 'orthorhombique', 'B', 0.4),
    ('Orthorhombic-C', 'orthorhombique', 'C', 0.4), ('Orthorhombic-P', 'orthorhombique', 'P', 0.4),
    ('Monoclinic-I', 'monoclinique', 'I', 0.7), ('Monoclinic-A', 'monoclinique', 'A', 0.7),
    ('Monoclinic-C', 'monoclinique', 'C', 0.7), ('Monoclinic-P', 'monoclinique', 'P', 0.7),
    ('Triclinic', 'triclinique', 'P', 1.0)]

def find_gsas2(user=None):
    """Retourne le dossier a ajouter a sys.path pour importer GSAS-II, ou None."""
    cands = [user, os.environ.get('GSAS2DIR')]
    for base in ('~/g2main', '~/gsas2main', '~/gsas2full', '~/g2full', '~/GSAS-II', '~/gsas2'):
        cands += [os.path.join(base, 'GSAS-II'), base]
    cands.append(os.path.join(os.path.dirname(os.path.dirname(sys.executable)), 'GSAS-II'))
    for c in cands:
        if not c:
            continue
        c = os.path.abspath(os.path.expanduser(c))
        if os.path.isfile(os.path.join(c, 'GSASII', 'GSASIIindex.py')):
            return c                                  # version recente : paquet GSASII
        if os.path.isfile(os.path.join(c, 'GSASIIindex.py')):
            return c                                  # ancienne disposition : modules a plat
    return None

def _import_g2index(g2dir):
    if g2dir and g2dir not in sys.path:
        sys.path.insert(0, g2dir)
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        try:
            from GSASII import GSASIIindex as G2indx
        except ImportError:
            import GSASIIindex as G2indx
    return G2indx

def _native_worker(args):
    """Recherche pour UN reseau de Bravais. Relance par tranches et s'arrete des qu'une bonne maille sort."""
    g2dir, ibrav, tth, inten, lam, zero, tmax = args
    t0 = time.time()
    try:
        G2indx = _import_g2index(g2dir)
    except Exception as e:
        return ibrav, [], 0., 'import impossible : %s' % e
    cells = []
    peaks = [[float(t), float(i), True, False, 0, 0, 0,
              lam / (2. * np.sin(np.radians((t - zero) / 2.))), 0.] for t, i in zip(tth, inten)]
    bravais = [i == ibrav for i in range(18)]
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            ok, dmin, cc = G2indx.DoIndexPeaks(peaks, [0, zero, 4, 25.0], bravais, None, timeout=tmax)
    except Exception as e:
        return ibrav, cells, time.time() - t0, 'erreur : %s' % e
    # X20 (c[1]) = nombre de pics non indexes parmi les 20 premiers : on ecarte les mailles qui en laissent trop
    cells = [[float(v) for v in c[:10]] for c in (cc if ok else []) if c[1] <= 3]
    cells.sort(key=lambda c: -c[0])
    return ibrav, cells[:8], time.time() - t0, ''

# --------------------------------------------------------------------------
# 3c. Moteurs Python : dichotomie (principe de DICVOL) et essai-erreur (principe de TREOR)
#     Ce sont des reecritures du principe de ces methodes, pas les programmes originaux.
#     Systemes traites : cubique, quadratique, hexagonal/trigonal, orthorhombique, monoclinique.
#     (Le triclinique est laisse a Conograph : trop long en Python pur.)
# --------------------------------------------------------------------------
# Dans chaque systeme, Q = 1/d^2 est LINEAIRE en quelques parametres p :
#   cubique        Q = (h2+k2+l2) A
#   quadratique    Q = (h2+k2) A + l2 C
#   hexagonal      Q = (h2+hk+k2) A + l2 C            (A = 4/(3a2))
#   orthorhombique Q = h2 A + k2 B + l2 C
#   monoclinique   Q = h2 A + k2 B + l2 C + hl D      (D = 2 a* c* cos(beta*), axe unique b)
_PY_SYSTEMS = ['cubique', 'quadratique', 'hexagonal', 'orthorhombique', 'monoclinique']
_SYS_NPAR = {'cubique': 1, 'quadratique': 2, 'hexagonal': 2, 'orthorhombique': 3, 'monoclinique': 4}
_SYS_MULT = {'cubique': 48., 'quadratique': 16., 'hexagonal': 24., 'orthorhombique': 8., 'monoclinique': 4.}
_COSB_MAX = 0.75          # beta entre 90 et ~131 degres
_CENT_OPTIONS = {'cubique': 'PIF', 'quadratique': 'PI', 'hexagonal': 'PR', 'orthorhombique': 'PCABIF',
                 'monoclinique': 'PCAI'}

def coef_matrix(system, hm):
    """Lignes distinctes des coefficients de Q pour tous les hkl jusqu'a l'indice hm."""
    r = np.arange(0, hm + 1)
    if system == 'monoclinique':
        h, k, l = [v.ravel() for v in np.meshgrid(r, r, np.arange(-hm, hm + 1), indexing='ij')]
        C = np.stack([h * h, k * k, l * l, h * l], 1)
    else:
        h, k, l = [v.ravel() for v in np.meshgrid(r, r, r, indexing='ij')]
        if system == 'cubique':
            C = (h * h + k * k + l * l).reshape(-1, 1)
        elif system == 'quadratique':
            C = np.stack([h * h + k * k, l * l], 1)
        elif system == 'hexagonal':
            C = np.stack([h * h + h * k + k * k, l * l], 1)
        else:
            C = np.stack([h * h, k * k, l * l], 1)
    C = np.unique(C, axis=0)
    return C[np.abs(C).sum(1) > 0].astype(float)

def params_to_cells(system, P):
    """Parametres lineaires (N x n) -> mailles (N x 6)."""
    P = np.atleast_2d(np.asarray(P, float))
    N = len(P)
    out = np.empty((N, 6)); out[:, 3:] = 90.
    if system == 'cubique':
        out[:, 0] = out[:, 1] = out[:, 2] = 1. / np.sqrt(P[:, 0])
    elif system == 'quadratique':
        out[:, 0] = out[:, 1] = 1. / np.sqrt(P[:, 0]); out[:, 2] = 1. / np.sqrt(P[:, 1])
    elif system == 'hexagonal':
        out[:, 0] = out[:, 1] = 2. / np.sqrt(3. * P[:, 0]); out[:, 2] = 1. / np.sqrt(P[:, 1]); out[:, 5] = 120.
    elif system == 'orthorhombique':
        out[:, :3] = 1. / np.sqrt(P[:, :3])
    else:
        A, B, C, D = P.T
        cb = np.clip(D / (2. * np.sqrt(A * C)), -0.999, 0.999)       # cos(beta*)
        sb = np.sqrt(1. - cb * cb)
        out[:, 0] = 1. / (np.sqrt(A) * sb); out[:, 1] = 1. / np.sqrt(B); out[:, 2] = 1. / (np.sqrt(C) * sb)
        out[:, 4] = 180. - np.degrees(np.arccos(cb))
    return out

def params_volume(system, P):
    P = np.atleast_2d(P)
    if system == 'cubique':
        return P[:, 0] ** -1.5
    if system == 'quadratique':
        return 1. / (P[:, 0] * np.sqrt(P[:, 1]))
    if system == 'hexagonal':
        return 2. / (np.sqrt(3.) * P[:, 0] * np.sqrt(P[:, 1]))
    if system == 'orthorhombique':
        return 1. / np.sqrt(P[:, 0] * P[:, 1] * P[:, 2])
    return 1. / np.sqrt(np.maximum(P[:, 1] * (P[:, 0] * P[:, 2] - P[:, 3] ** 2 / 4.), 1e-30))

def eval_params(P, C, qobs, eps, max_bytes=6e7):
    """Pour chaque jeu de parametres : nombre de raies observees indexees (a eps pres), ecart moyen en Q,
    nombre de raies calculees jusqu'a la derniere raie observee."""
    P = np.atleast_2d(P)
    N, K, M = len(P), len(qobs), len(C)
    nidx = np.zeros(N, int); dq = np.zeros(N); ncalc = np.zeros(N, int)
    step = max(1, int(max_bytes / (8. * M * K)))
    for s in range(0, N, step):
        Q = P[s:s + step] @ C.T
        d = np.abs(Q[:, :, None] - qobs[None, None, :]).min(1)
        nidx[s:s + step] = (d <= eps).sum(1)
        dq[s:s + step] = np.minimum(d, 3. * eps).mean(1)
        ncalc[s:s + step] = (Q <= qobs[-1] + eps[-1]).sum(1)
    return nidx, dq, ncalc

def _prep_obs(tth, lam, zero, nlines, tol_deg):
    t = np.sort(np.asarray(tth, float))[:nlines] - zero
    th = np.radians(t / 2.)
    qobs = (2. * np.sin(th) / lam) ** 2
    eps = 2. * np.sin(2. * th) / lam ** 2 * np.radians(tol_deg)       # dQ pour une erreur tol_deg sur 2theta
    return qobs, eps

def _vmax_system(system, qobs, vmax_user):
    """Volume maximal raisonnable : au plus ~4 fois plus de raies calculees que de raies observees."""
    est = 4. * len(qobs) * _SYS_MULT[system] / (4.19 * qobs[-1] ** 1.5)
    return min(vmax_user, est)

def _top_candidates(system, P, C, qobs, eps, nimp, engine, nkeep=12, relax=1.):
    if len(P) == 0:
        return []
    # regroupe les solutions quasi identiques
    key = np.round(np.log(np.maximum(np.abs(P), 1e-12)) / 0.004).astype(np.int64)
    _, iu = np.unique(key, axis=0, return_index=True)
    P = P[iu]
    # relax > 1 : le centre d'une boite de dichotomie n'est pas encore affine, on tolere un ecart plus grand
    nidx, dq, ncalc = eval_params(P, C, qobs, eps * relax)
    ok = nidx >= len(qobs) - nimp
    P, dq, ncalc = P[ok], dq[ok], ncalc[ok]
    if len(P) == 0:
        return []
    m20 = qobs[-1] / (2. * np.maximum(dq, 1e-9) * np.maximum(ncalc, 1))
    o = np.argsort(-m20)[:nkeep]
    cells = params_to_cells(system, P[o])
    return [dict(engine=engine, label='', cell=[float(v) for v in c], system=system, cent='P', bravais='',
                 m20v=float(m)) for c, m in zip(cells, m20[o])]

# Arret anticipe (comportement par defaut de DICVOL) : des qu'un moteur a une bonne solution dans un systeme,
# ses recherches en cours dans les systemes de symetrie plus basse sont arretees.
_SYS_LEVEL = {'cubique': 0, 'hexagonal': 1, 'quadratique': 1, 'orthorhombique': 2, 'monoclinique': 3}
_ENG_SLOT = {'dichotomie': 0, 'treor': 1}
_STOP = None

def _init_pool(shared):
    global _STOP
    _STOP = shared

def _should_stop(engine, system):
    if _STOP is None:
        return False
    v = _STOP[_ENG_SLOT[engine]]
    return v > 0 and (v - 1) < _SYS_LEVEL[system]

def _signal_found(engine, system):
    if _STOP is None:
        return
    k, lev = _ENG_SLOT[engine], _SYS_LEVEL[system] + 1
    with _STOP.get_lock():
        if _STOP[k] == 0 or lev < _STOP[k]:
            _STOP[k] = lev

_STOP_MSG = 'arrete : solution deja trouvee en symetrie plus haute'

# ---------------- essai-erreur (principe de TREOR) ----------------
def _trial_rows(system):
    """Petits indices essayes pour les raies de base."""
    if system == 'cubique':
        return coef_matrix(system, 4)[:24]
    if system in ('quadratique', 'hexagonal'):
        C = coef_matrix(system, 4)
        return C[(C[:, 0] <= 13) & (C[:, 1] <= 16)]
    if system == 'orthorhombique':
        return coef_matrix(system, 2)
    r = np.arange(0, 3)
    h, k, l = [v.ravel() for v in np.meshgrid(r, r, np.arange(-2, 3), indexing='ij')]
    m = (h + k + np.abs(l) <= 3) & (h + k + np.abs(l) > 0) & ~((h == 0) & (l < 0))
    return np.unique(np.stack([h * h, k * k, l * l, h * l], 1)[m], axis=0).astype(float)

def treor_system(system, tth, lam, zero, opt):
    """Essai-erreur : on attribue de petits hkl a quelques raies de bas angle, on en deduit la maille
    (systeme lineaire), puis on garde les mailles qui indexent toutes les raies (sauf nimp)."""
    import itertools
    t0 = time.time()
    n = _SYS_NPAR[system]
    qobs, eps = _prep_obs(tth, lam, zero, opt['nlines'], opt['tol'])
    if len(qobs) < n + 3:
        return [], 0., 'trop peu de raies'
    amax, amin = opt['amax'], 2.0
    hm = int(min(10, np.ceil(amax * np.sqrt(qobs[-1])) + 1))
    C = coef_matrix(system, hm)
    pmin = max(1. / amax ** 2, qobs[-1] / hm ** 2); pmax = 1. / amin ** 2
    if system == 'hexagonal':
        pmin_a, pmax_a = pmin * 4. / 3., pmax * 4. / 3.
    vmax = _vmax_system(system, qobs, opt['vmax'])
    T = _trial_rows(system)
    nb = {1: 4, 2: 5, 3: 6, 4: 6}[n]
    nb = min(nb, len(qobs))
    idx = np.indices((len(T),) * n).reshape(n, -1).T
    mats = T[idx]                                        # (T^n, n, n)
    good = np.abs(np.linalg.det(mats)) > 0.5
    mats = mats[good]
    sols = []
    msg = ''
    for combo in itertools.combinations(range(nb), n):
        if time.time() - t0 > opt['budget']:
            msg = 'interrompu (budget de temps)'
            break
        if opt.get('early') and _should_stop('treor', system):
            msg = _STOP_MSG
            break
        qb = qobs[list(combo)]
        p = np.linalg.solve(mats, np.broadcast_to(qb, (len(mats), n))[..., None])[..., 0]
        if system == 'monoclinique':
            p[:, 3] = np.abs(p[:, 3])
            swap = p[:, 0] > p[:, 2]
            p[swap, 0], p[swap, 2] = p[swap, 2].copy(), p[swap, 0].copy()
            ok = (np.all(p[:, :3] >= pmin, 1) & np.all(p[:, :3] <= pmax, 1)
                  & (p[:, 3] <= 2. * _COSB_MAX * np.sqrt(np.abs(p[:, 0] * p[:, 2]))))
        elif system == 'hexagonal':
            ok = (p[:, 0] >= pmin_a) & (p[:, 0] <= pmax_a) & (p[:, 1] >= pmin) & (p[:, 1] <= pmax)
        else:
            if system == 'orthorhombique':
                p = np.sort(p, axis=1)
            ok = np.all(p >= pmin, 1) & np.all(p <= pmax, 1)
        p = p[ok & np.all(np.isfinite(p), 1)]
        if len(p) == 0:
            continue
        p = p[params_volume(system, p) <= vmax]
        if len(p) == 0:
            continue
        key = np.round(np.log(np.maximum(p, 1e-12)) / 0.002).astype(np.int64)
        _, iu = np.unique(key, axis=0, return_index=True)
        p = p[iu]
        nidx, dq, ncalc = eval_params(p, C, qobs, eps)
        keep = nidx >= len(qobs) - opt['nimp']
        if keep.any():
            sols.append(p[keep])
    P = np.concatenate(sols) if sols else np.zeros((0, n))
    return _top_candidates(system, P, C, qobs, eps, opt['nimp'], 'Essai-erreur'), time.time() - t0, msg

# ---------------- dichotomie (principe de DICVOL) ----------------
def _box_volume_bounds(system, lo, hi):
    if system == 'monoclinique':
        vmin = 1. / np.sqrt(hi[:, 1] * np.maximum(hi[:, 0] * hi[:, 2] - lo[:, 3] ** 2 / 4., 1e-30))
        vmax = 1. / np.sqrt(lo[:, 1] * np.maximum(lo[:, 0] * lo[:, 2] - hi[:, 3] ** 2 / 4., 1e-30))
        return vmin, vmax
    return params_volume(system, hi), params_volume(system, lo)

def _box_valid(system, lo, hi):
    if system == 'orthorhombique':                       # ordre A <= B <= C (a >= b >= c)
        return (lo[:, 0] <= hi[:, 1]) & (lo[:, 1] <= hi[:, 2])
    if system == 'monoclinique':                         # A <= C et beta* admissible
        return (lo[:, 0] <= hi[:, 2]) & (lo[:, 3] <= 2. * _COSB_MAX * np.sqrt(hi[:, 0] * hi[:, 2]))
    return np.ones(len(lo), bool)

def _box_test(lo, hi, Cp, Cn, qlo, qhi, nimp, max_bytes=8e7):
    """Une boite est conservee si chaque raie observee tombe dans l'intervalle [Qmin, Qmax] d'au moins un hkl."""
    N, M = len(lo), len(Cp)
    keep = np.zeros(N, bool)
    step = max(1, int(max_bytes / (16. * M)))
    for s in range(0, N, step):
        l, h = lo[s:s + step], hi[s:s + step]
        Qmin = l @ Cp.T + h @ Cn.T
        Qmax = h @ Cp.T + l @ Cn.T
        alive = np.arange(len(l)); miss = np.zeros(len(l), int)
        for i in range(len(qlo)):
            ok = ((Qmin <= qhi[i]) & (Qmax >= qlo[i])).any(1)
            miss[alive[~ok]] += 1
            sel = miss[alive] <= nimp
            if not sel.all():
                alive, Qmin, Qmax = alive[sel], Qmin[sel], Qmax[sel]
            if len(alive) == 0:
                break
        keep[s + alive] = True
    return keep

def _split(lo, hi):
    n = lo.shape[1]
    mid = 0.5 * (lo + hi)
    L, H = [], []
    for bits in range(2 ** n):
        b = np.array([(bits >> j) & 1 for j in range(n)], bool)
        L.append(np.where(b, mid, lo)); H.append(np.where(b, hi, mid))
    return np.concatenate(L), np.concatenate(H)

def dicho_system(system, tth, lam, zero, opt):
    """Dichotomie : l'espace des parametres est decoupe en boites ; on elimine toute boite dans laquelle
    une raie observee ne peut etre indexee par aucun hkl, on coupe les autres en deux, et on recommence.
    Recherche par tranches de volume croissant : on s'arrete a la premiere tranche qui contient des solutions."""
    t0 = time.time()
    n = _SYS_NPAR[system]
    qobs, eps = _prep_obs(tth, lam, zero, opt['nlines'], opt['tol'])
    if len(qobs) < n + 3:
        return [], 0., 'trop peu de raies'
    amax, amin = opt['amax'], 2.0
    hm = int(min(10, np.ceil(amax * np.sqrt(qobs[-1])) + 1))
    C = coef_matrix(system, hm)
    Cp, Cn = np.maximum(C, 0.), np.minimum(C, 0.)
    qlo, qhi = qobs - eps, qobs + eps
    pmin = max(1. / amax ** 2, qobs[-1] / hm ** 2); pmax = 1. / amin ** 2
    vmax = _vmax_system(system, qobs, opt['vmax'])
    ratio = 1.30
    n0 = int(np.ceil(np.log(pmax / pmin) / np.log(ratio)))
    e = pmin * (pmax / pmin) ** (np.arange(n0 + 1) / float(n0))
    nlev = 7
    # boites initiales
    if system == 'cubique':
        lo, hi = e[:-1, None], e[1:, None]
    elif system in ('quadratique', 'hexagonal'):
        f = 4. / 3. if system == 'hexagonal' else 1.
        i, j = [v.ravel() for v in np.meshgrid(np.arange(n0), np.arange(n0), indexing='ij')]
        lo = np.stack([f * e[i], e[j]], 1); hi = np.stack([f * e[i + 1], e[j + 1]], 1)
    else:
        i, j, k = [v.ravel() for v in np.meshgrid(np.arange(n0), np.arange(n0), np.arange(n0), indexing='ij')]
        lo = np.stack([e[i], e[j], e[k]], 1); hi = np.stack([e[i + 1], e[j + 1], e[k + 1]], 1)
        if system == 'monoclinique':
            m = lo[:, 0] <= hi[:, 2]
            lo, hi = lo[m], hi[m]
            nd = 4
            dmax = 2. * _COSB_MAX * np.sqrt(hi[:, 0] * hi[:, 2])
            L, H = [], []
            for q in range(nd):
                L.append(np.column_stack([lo, dmax * q / nd])); H.append(np.column_stack([hi, dmax * (q + 1) / nd]))
            lo, hi = np.concatenate(L), np.concatenate(H)
    m = _box_valid(system, lo, hi)
    lo, hi = lo[m], hi[m]
    bvmin, bvmax = _box_volume_bounds(system, lo, hi)
    shells = vmax * np.array([0., 0.08, 0.16, 0.28, 0.42, 0.60, 0.80, 1.0])
    msg, found = '', []
    nimp = opt['nimp']
    for s in range(len(shells) - 1):
        v1, v2 = shells[s], shells[s + 1]
        m = (bvmax >= v1) & (bvmin <= v2)
        stack = [(lo[m], hi[m], 0)]
        centers = []
        while stack:
            if time.time() - t0 > opt['budget']:
                msg = 'interrompu (budget de temps), volume explore < %.0f A3' % v1
                stack = []
                break
            if opt.get('early') and _should_stop('dichotomie', system):
                msg = _STOP_MSG
                stack = []
                break
            l, h, lev = stack.pop()
            if len(l) == 0:
                continue
            if len(l) > 20000:                           # traitement par paquets (memoire)
                stack.append((l[20000:], h[20000:], lev))
                l, h = l[:20000], h[:20000]
            k = _box_valid(system, l, h)
            a, b = _box_volume_bounds(system, l, h)
            k &= (b >= v1) & (a <= v2)
            l, h = l[k], h[k]
            if len(l) == 0:
                continue
            k = _box_test(l, h, Cp, Cn, qlo, qhi, nimp)
            l, h = l[k], h[k]
            if len(l) == 0:
                continue
            if lev >= nlev:
                centers.append(0.5 * (l + h))
            else:
                l2, h2 = _split(l, h)
                stack.append((l2, h2, lev + 1))
        if centers:
            found = np.concatenate(centers)
            break
        if msg:
            break
    P = found if len(found) else np.zeros((0, n))
    if system == 'monoclinique' and len(P):
        P = P[P[:, 3] < 2. * 0.999 * np.sqrt(P[:, 0] * P[:, 2])]
    return _top_candidates(system, P, C, qobs, eps, nimp, 'Dichotomie', relax=3.), time.time() - t0, msg

def mono_standard(cell):
    """Reduit le plan (a, c) d'une maille monoclinique : a et c les plus courts possibles, beta >= 90, a <= c."""
    a, b, c, be = cell[0], cell[1], cell[2], np.radians(cell[4])
    u = np.array([a, 0.]); v = np.array([c * np.cos(be), c * np.sin(be)])
    for _ in range(50):
        if u @ u > v @ v:
            u, v = v, u
        m = np.round((u @ v) / (u @ u))
        if m == 0:
            break
        v = v - m * u
    if u @ v > 0:
        v = -v
    lu, lv = np.sqrt(u @ u), np.sqrt(v @ v)
    beta = np.degrees(np.arccos(np.clip((u @ v) / (lu * lv), -1, 1)))
    return [float(lu), float(b), float(lv), 90., float(beta), 90.]

def assign_centering(c, tth, lam, zero):
    """Les moteurs Python cherchent une maille P. On l'affine (maille + zero), puis on teste si un mode de
    reseau centre indexe aussi bien toutes les raies : il donne moins de raies calculees, donc un meilleur M20."""
    system = c['system']
    cell = list(c['cell'])
    if system == 'monoclinique':
        cell = mono_standard(cell)
    elif system == 'orthorhombique':                     # convention d'affichage : a <= b <= c
        cell = sorted(cell[:3]) + [90., 90., 90.]
    z = zero
    try:
        cell, z, _ = refine_cell(cell, system, 'P', tth, lam, zero)
    except Exception:
        pass
    best, bm, bn = 'P', None, None
    for cent in _CENT_OPTIONS.get(system, 'P'):
        try:
            f = figures_of_merit(cell, cent, tth, lam, z, 20)
            fa = figures_of_merit(cell, cent, tth, lam, z, len(tth))
        except Exception:
            continue
        if cent == 'P':
            bm, bn = f['M'], fa['nunidx']
        elif bm is not None and fa['nunidx'] <= bn and f['M'] > bm * 1.05:
            best, bm = cent, f['M']
    if system == 'monoclinique' and best == 'A':         # convention : C plutot que A (echange de a et c)
        cell = [cell[2], cell[1], cell[0], 90., cell[4], 90.]; best = 'C'
    if system == 'hexagonal' and best == 'R':
        system = 'rhomboedrique'
    names = dict(((sy, ce), nm) for nm, sy, ce, _ in _G2_BRAVAIS)
    return dict(c, cell=[float(v) for v in cell], cent=best, system=system,
                bravais=names.get((system, best), system))

def _py_worker(args):
    """Une tache = un moteur Python sur un systeme cristallin (execute dans un processus separe).
    Passes successives tant que rien n'est trouve : stricte, puis une raie non indexee toleree,
    puis tolerance angulaire elargie (erreur de zero, pics imprecis)."""
    engine, system, tth, lam, zero, opt = args
    fn = dicho_system if engine == 'dichotomie' else treor_system
    tth = np.asarray(tth, float)
    cands, dt, msg, note = [], 0., '', ''
    # etapes : chaque etape regroupe une ou deux variantes (nimp = raies non indexees tolerees, facteur de tolerance)
    stages = [[(0, 1.)], [(0, 1.7), (1, 1.)], [(1, 1.7)]]
    try:
        for stage in stages:
            notes = []
            for nimp, ftol in stage:
                left = opt['budget'] - dt
                if left < 2.:
                    msg = msg or 'interrompu (budget de temps)'
                    break
                o = dict(opt, nimp=nimp, tol=opt['tol'] * ftol, budget=left)
                cc, d, m = fn(system, tth, lam, zero, o)
                dt += d
                msg = msg or m
                if cc:
                    cands += cc
                    if ftol > 1.:
                        notes.append('tolerance %.3f deg' % o['tol'])
                    if nimp:
                        notes.append('1 raie non indexee toleree')
                if m:
                    break
            if cands or msg:
                note = ('trouve avec ' + ' / '.join(notes)) if notes else ''
                break
        if cands and max(c['m20v'] for c in cands) >= 10.:
            _signal_found(engine, system)
    except Exception as e:
        return engine, system, [], dt, 'erreur : %s' % e
    return engine, system, cands, dt, msg or note

# --------------------------------------------------------------------------
# 4. Affinement de la maille sur les pics, comparaison, classement
# --------------------------------------------------------------------------
def _pack(cell, system):
    a, b, c, al, be, ga = cell[:6]
    if system == 'cubique': return [a]
    if system in ('quadratique', 'hexagonal', 'rhomboedrique'): return [a, c]
    if system == 'orthorhombique': return [a, b, c]
    if system == 'monoclinique': return [a, b, c, be]
    return [a, b, c, al, be, ga]

def _unpack(p, system):
    if system == 'cubique': return [p[0], p[0], p[0], 90., 90., 90.]
    if system == 'quadratique': return [p[0], p[0], p[1], 90., 90., 90.]
    if system in ('hexagonal', 'rhomboedrique'): return [p[0], p[0], p[1], 90., 90., 120.]
    if system == 'orthorhombique': return [p[0], p[1], p[2], 90., 90., 90.]
    if system == 'monoclinique': return [p[0], p[1], p[2], 90., p[3], 90.]
    return list(p[:6])

def refine_cell(cell, system, cent, tth_obs, lam, zero0=0., tol_deg=0.06, refine_zero=True):
    """Moindres carres sur 2theta (maille + zero), en n'utilisant que les pics indexes a tol_deg pres."""
    try:
        from scipy.optimize import least_squares
    except ImportError:
        return list(cell[:6]), zero0, 0
    tth = np.sort(np.asarray(tth_obs, float))
    qmax = (2. * np.sin(np.radians(tth[-1] / 2.)) / lam) ** 2 * 1.03
    def nearest(p):
        c = _unpack(p[:-1], system)
        qc, _ = calc_q(c, qmax, cent)
        tc = 2. * np.degrees(np.arcsin(np.clip(lam * np.sqrt(qc) / 2., 0, 1))) + p[-1]
        i = np.clip(np.searchsorted(tc, tth), 1, len(tc) - 1)
        return np.where(np.abs(tth - tc[i - 1]) < np.abs(tth - tc[i]), tc[i - 1], tc[i])
    p = np.array(_pack(cell, system) + [zero0], float)
    used = 0
    for _ in range(3):
        try:
            sel = np.abs(tth - nearest(p)) < tol_deg
        except Exception:
            break
        used = int(sel.sum())
        if used < len(p) + 2:
            break
        def res(pp):
            if not refine_zero:
                pp = np.concatenate((pp[:-1], [zero0]))
            try:
                return (tth - nearest(pp))[sel]
            except Exception:
                return np.full(used, 1.)
        lo = np.concatenate((p[:-1] * 0.97, [zero0 - 0.15])); hi = np.concatenate((p[:-1] * 1.03, [zero0 + 0.15]))
        r = least_squares(res, p, bounds=(np.minimum(lo, hi), np.maximum(lo, hi)), diff_step=1e-5)
        p = r.x
    if not refine_zero:
        p[-1] = zero0
    return _unpack(p[:-1], system), float(p[-1]), used

def same_lattice(c1, cent1, c2, cent2, tol=0.012):
    """Deux mailles decrivent-elles le meme reseau ? (independant du choix d'axes)"""
    v1, v2 = prim_volume(c1, cent1), prim_volume(c2, cent2)
    if abs(v1 - v2) > 0.05 * v2:                       # exclut les surmailles et sous-mailles
        return False
    return _same_sig(_qsig(c1, cent1), _qsig(c2, cent2), tol)

_CENT_MULT = {'P': 1., 'A': 2., 'B': 2., 'C': 2., 'I': 2., 'F': 4., 'R': 3.}

def prim_volume(cell, cent):
    return cell_to_recip_metric(cell)[1] / _CENT_MULT.get(cent, 1.)

def _qsig(cell, cent, nq=40):
    q, _ = calc_q(cell, 1.0, cent)
    return q[:nq]

def _same_sig(q1, q2, tol=0.006):
    """Meme reseau si, sur le domaine commun, chaque Q calcule de l'un a un voisin chez l'autre (a 2*tol pres
    en relatif). Comparaison par ensembles : insensible au choix d'axes et aux raies dedoublees d'une
    description de symetrie plus basse (ex. monoclinique avec beta = 90.03)."""
    if len(q1) < 8 or len(q2) < 8:
        return False
    qlim = min(q1[-1], q2[-1]) * 0.98
    a, b = q1[q1 <= qlim], q2[q2 <= qlim]
    if len(a) < 6 or len(b) < 6:
        return False
    def covered(x, y):
        i = np.clip(np.searchsorted(y, x), 1, len(y) - 1)
        d = np.minimum(np.abs(x - y[i - 1]), np.abs(x - y[i]))
        return bool(np.all(d / x < 2 * tol))
    return covered(a, b) and covered(b, a)

_SYM_ORDER = ['cubique', 'hexagonal', 'rhomboedrique', 'quadratique', 'orthorhombique', 'monoclinique', 'triclinique']

def rank_candidates(cands, tth, lam, zero, nfom=20, refine=True, nkeep=60):
    out = []
    # pre-tri sans affinement : seules les meilleures candidates sont affinees (gain de temps)
    if refine and len(cands) > nkeep:
        pre = []
        for c in cands:
            try:
                f = figures_of_merit(c['cell'], c['cent'], tth, lam, zero, len(tth))
                f20 = figures_of_merit(c['cell'], c['cent'], tth, lam, zero, nfom)
                pre.append((f20['M'] * (1. - f['nunidx'] / float(len(tth))) ** 2, id(c), c))
            except Exception:
                pass
        pre.sort(key=lambda t: -t[0])
        cands, per = [t[2] for t in pre[:nkeep]], {}
        for t in pre[nkeep:]:                      # et toujours les 5 meilleures de chaque type de reseau,
            k = (t[2]['system'], t[2]['cent'])     # car une erreur de zero penalise surtout les hautes symetries
            if per.get(k, 0) < 5:
                per[k] = per.get(k, 0) + 1
                cands.append(t[2])
    for c in cands:
        cell, z = c['cell'], zero
        try:
            if refine:
                cell, z, _ = refine_cell(c['cell'], c['system'], c['cent'], tth, lam, zero)
            fom = figures_of_merit(cell, c['cent'], tth, lam, z, nfom)
            allp = figures_of_merit(cell, c['cent'], tth, lam, z, len(tth))
        except Exception:
            continue
        d = dict(c); d.update(cell_ref=cell, zero=z, M20=fom['M'], F20=fom['F'], nunidx20=fom['nunidx'],
                             ncalc20=fom['ncalc'], vol=fom['vol'], nunidx_all=allp['nunidx'], nobs=len(tth),
                             sig=_qsig(cell, c['cent']), vprim=prim_volume(cell, c['cent']),
                             engines=set([c['engine']]))
        # score = M20, penalise par les pics non indexes, avec un leger bonus a la symetrie la plus haute
        # (a qualite egale on essaie d'abord la description la plus symetrique)
        d['score'] = (d['M20'] * (1. - d['nunidx_all'] / float(max(d['nobs'], 1))) ** 2
                      * (1. + 0.08 * (len(_SYM_ORDER) - 1 - _SYM_ORDER.index(d['system']))))
        out.append(d)
    # fusion des mailles equivalentes (meme reseau, quel que soit le choix d'axes). La description principale est
    # celle de meilleur score ou, a score voisin, de plus haute symetrie. Les descriptions du meme reseau dans un
    # AUTRE systeme cristallin sont conservees a part ("alts") : en cas de pseudo-symetrie (ex. monoclinique avec
    # a ~ c, qui ressemble a un orthorhombique C), c'est la description de plus basse symetrie qui est la vraie.
    out.sort(key=lambda d: (-d['score']))
    uniq = []
    for d in out:
        for u in uniq:
            if abs(d['vprim'] - u['vprim']) < 0.04 * u['vprim'] and _same_sig(d['sig'], u['sig']):
                eng = u['engines'] | d['engines']
                alts = u.get('alts', {})
                if d['system'] == u['system']:
                    pass
                elif (_SYM_ORDER.index(d['system']) < _SYM_ORDER.index(u['system'])
                        and d['score'] > 0.8 * u['score']):
                    old = dict((k, v) for k, v in u.items() if k != 'alts')
                    old['engines'] = set(old['engines'])
                    alts[old['system']] = old
                    alts.pop(d['system'], None)
                    u.clear(); u.update(d)
                else:
                    cur = alts.get(d['system'])
                    if cur is None:
                        alts[d['system']] = dict(d, engines=set(d['engines']))
                    elif d['cent'] == 'P' and cur['cent'] != 'P':      # on prefere la description primitive
                        alts[d['system']] = dict(d, engines=set(d['engines']) | cur['engines'])
                    else:
                        cur['engines'] |= d['engines']
                u['engines'] = eng
                u['alts'] = alts
                break
        else:
            uniq.append(d)
    uniq.sort(key=lambda d: -d['score'])
    # liste finale : chaque reseau, suivi (pour les 3 premiers) de ses descriptions de symetrie plus basse
    final = []
    for i, u in enumerate(uniq):
        alts = u.pop('alts', {})
        u['equiv_of'] = 0
        final.append(u)
        main_rank = len(final)
        if i < 3:
            if (u['system'] == 'orthorhombique' and u['cent'] in 'ABC'
                    and alts.get('monoclinique', {}).get('cent') != 'P'):
                d = _derive_mono_from_ortho_c(u, tth, lam, zero, nfom)
                if d is not None:
                    alts['monoclinique'] = d
            for sysname in _SYM_ORDER:
                if sysname in alts and _SYM_ORDER.index(sysname) > _SYM_ORDER.index(u['system']):
                    al = alts[sysname]
                    al['equiv_of'] = main_rank
                    final.append(al)
    return final

def _derive_mono_from_ortho_c(u, tth, lam, zero, nfom=20):
    """Un reseau orthorhombique a une face centree est aussi un reseau monoclinique P dont deux parametres
    sont egaux. On construit cette maille monoclinique et on l'affine librement (a, b, c, beta independants)."""
    try:
        a, b, c = u['cell_ref'][:3]
        if u['cent'] == 'A':
            a, b, c = b, c, a            # face centree (b, c) -> on la ramene dans (a, b)
        elif u['cent'] == 'B':
            a, b, c = a, c, b
        e = 0.5 * np.sqrt(a * a + b * b)
        beta = np.degrees(np.arccos((a * a - b * b) / (a * a + b * b)))
        if beta < 90.:
            beta = 180. - beta
        cell0 = [e, c, e, 90., beta, 90.]
        cell, z, _ = refine_cell(cell0, 'monoclinique', 'P', tth, lam, zero)
        cell = mono_standard(cell)
        fom = figures_of_merit(cell, 'P', tth, lam, z, nfom)
        allp = figures_of_merit(cell, 'P', tth, lam, z, len(tth))
        d = dict(engine='derivee', label='', cell=cell, bravais='Monoclinic-P', system='monoclinique', cent='P',
                 cell_ref=cell, zero=z, M20=fom['M'], F20=fom['F'], nunidx20=fom['nunidx'], ncalc20=fom['ncalc'],
                 vol=fom['vol'], nunidx_all=allp['nunidx'], nobs=len(tth), sig=_qsig(cell, 'P'),
                 vprim=prim_volume(cell, 'P'), engines=set(['derivee']))
        d['score'] = d['M20'] * (1. - d['nunidx_all'] / float(max(d['nobs'], 1))) ** 2
        return d
    except Exception:
        return None

def format_table(ranked, nmax=12, ref=None):
    L = []
    L.append('%3s %-15s %8s %8s %8s %7s %7s %7s %9s %6s %6s %7s  %-18s %s' % (
        '#', 'reseau', 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'V(A3)', 'M20', 'F20', 'non-idx', 'trouvee par', ''))
    for i, d in enumerate(ranked[:nmax]):
        c = d['cell_ref']
        flag = ''
        if d.get('equiv_of'):
            flag = '= meme reseau que n. %d, symetrie plus basse' % d['equiv_of']
        if ref is not None and same_lattice(c, d['cent'], ref[0], ref[1]):
            flag = (flag + ' ; ' if flag else '') + '<= phase du projet'
        L.append('%3d %-15s %8.4f %8.4f %8.4f %7.3f %7.3f %7.3f %9.2f %6.1f %6.1f %3d/%-3d  %-18s %s' % (
            i + 1, d['bravais'], c[0], c[1], c[2], c[3], c[4], c[5], d['vol'], d['M20'], d['F20'],
            d['nunidx_all'], d['nobs'], '+'.join(sorted(d['engines'])), flag))
    return '\n'.join(L)

# --------------------------------------------------------------------------
# 4b. Transmission des mailles choisies a GSAS-II (copie du projet, Unit Cells List)
# --------------------------------------------------------------------------
_IBRAV = {('cubique', 'F'): 0, ('cubique', 'I'): 1, ('cubique', 'P'): 2, ('rhomboedrique', 'R'): 3,
          ('hexagonal', 'P'): 4, ('quadratique', 'I'): 5, ('quadratique', 'P'): 6,
          ('orthorhombique', 'F'): 7, ('orthorhombique', 'I'): 8, ('orthorhombique', 'A'): 9,
          ('orthorhombique', 'B'): 10, ('orthorhombique', 'C'): 11, ('orthorhombique', 'P'): 12,
          ('monoclinique', 'I'): 13, ('monoclinique', 'A'): 14, ('monoclinique', 'C'): 15,
          ('monoclinique', 'P'): 16, ('triclinique', 'P'): 17}
_BRAV_SYMB = ['Fm3m', 'Im3m', 'Pm3m', 'R3-H', 'P6/mmm', 'I4/mmm', 'P4/mmm', 'Fmmm', 'Immm', 'Ammm', 'Bmmm',
              'Cmmm', 'Pmmm', 'I2/m', 'A2/m', 'C2/m', 'P2/m', 'P1']
_BRAV_SG = ['F m 3 m', 'I m 3 m', 'P m 3 m', 'R 3 m', 'P 6/m m m', 'I 4/m m m', 'P 4/m m m', 'F m m m',
            'I m m m', 'A m m m', 'B m m m', 'C m m m', 'P m m m', 'I 2/m', 'A 2/m', 'C 2/m', 'P 2/m', 'P -1']

def ask_selection(n, default=5):
    """Demande a l'utilisateur les numeros des mailles a transmettre. Retourne une liste d'indices (base 0)."""
    if not sys.stdin.isatty():
        return []
    try:
        txt = input("\nMailles a transmettre a GSAS-II pour verification visuelle ?\n"
                    "  numeros separes par des virgules (ex. 1,3), 'a' = les %d premieres, Entree = aucune : "
                    % min(default, n)).strip().lower()
    except (EOFError, KeyboardInterrupt):
        return []
    return parse_selection(txt, n, default)

def parse_selection(txt, n, default=5):
    if not txt:
        return []
    if txt in ('a', 'all', 'tout', 'toutes'):
        return list(range(min(default, n)))
    sel = []
    for t in re.split(r'[,; ]+', txt):
        if t.isdigit() and 1 <= int(t) <= n and int(t) - 1 not in sel:
            sel.append(int(t) - 1)
    return sel

def build_unit_cells_data(old_uc, chosen, tth_max, lam):
    """Construit le contenu de l'item 'Unit Cells List' de GSAS-II a partir des mailles choisies.
    chosen : liste de dict avec cell_ref, system, cent, M20, nunidx20, vol, zero."""
    cells = []
    for k, d in enumerate(chosen):
        c = d['cell_ref']
        ib = _IBRAV.get((d['system'], d['cent']), 17)
        cells.append([float(d['M20']), int(d['nunidx20']), ib] + [float(v) for v in c[:6]] + [float(d['vol']),
                     k == 0, True])                  # [M20, X20, ibrav, a..gamma, V, 'use', 'keep']
    uc = old_uc or []
    controls = list(uc[0]) if len(uc) > 0 and len(uc[0]) >= 12 else [0, 0.0, 4, 25.0, 0, 'P1', 1., 1., 1., 90., 90., 90.]
    ssopt = uc[4] if len(uc) > 4 and isinstance(uc[4], dict) else {}
    magcells = uc[5] if len(uc) > 5 else []
    b = chosen[0]
    ib = _IBRAV.get((b['system'], b['cent']), 17)
    controls[1] = float(b['zero'])                   # zero affine avec la meilleure maille
    controls[5] = _BRAV_SYMB[ib]
    controls[6:12] = [float(v) for v in b['cell_ref'][:6]]
    controls = controls[:12] + [float(b['vol']), _BRAV_SG[ib]] + list(controls[14:15] or [[]])
    bravais = [i == ib for i in range(18)]
    dmin = float(lam / (2. * np.sin(np.radians(tth_max / 2.)))) - 0.005
    return [controls, bravais, cells, dmin, ssopt, magcells]

def build_index_peaks(tth, inten, lam):
    rows = [[float(t), float(i), True, False, 0, 0, 0,
             float(lam / (2. * np.sin(np.radians(t / 2.)))), 0.] for t, i in zip(tth, inten)]
    return [rows, [0.] * len(rows)]

def write_cells_to_gpx(gpx, out_gpx, hist_name, chosen, tth, inten, lam, g2dir):
    """Ecrit les mailles choisies dans la 'Unit Cells List' d'une COPIE du projet (l'original n'est pas modifie).
    Dans GSAS-II : ouvrir la copie > Unit Cells List > cliquer sur une ligne pour superposer les raies calculees."""
    if g2dir and g2dir not in sys.path:
        sys.path.insert(0, g2dir)
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        try:
            from GSASII import GSASIIscriptable as G2sc
        except ImportError:
            import GSASIIscriptable as G2sc
        proj = G2sc.G2Project(gpx)
    hist = [h for h in proj.histograms() if h.name == hist_name]
    if not hist:
        raise RuntimeError('histogramme introuvable dans le projet')
    data = hist[0].data
    data['Unit Cells List'] = build_unit_cells_data(data.get('Unit Cells List'), chosen, max(tth), lam)
    ipl = data.get('Index Peak List')
    if not ipl or not len(ipl[0]):                    # Index Peak List vide : on y copie les pics utilises
        data['Index Peak List'] = build_index_peaks(tth, inten, lam)
    with contextlib.redirect_stdout(io.StringIO()):
        proj.save(out_gpx)
    return out_gpx

# --------------------------------------------------------------------------
# 5. Programme principal
# --------------------------------------------------------------------------
_METHODS = [('conograph', 'Conograph'), ('dichotomie', 'Dichotomie'), ('treor', 'Essai-erreur'), ('gsas', 'GSAS-II')]
_ALIASES = {'c': 'conograph', 'cono': 'conograph', 'conograph': 'conograph',
            'd': 'dichotomie', 'dicho': 'dichotomie', 'dichotomie': 'dichotomie', 'dicvol': 'dichotomie',
            't': 'treor', 'treor': 'treor', 'essai': 'treor', 'essai-erreur': 'treor', 'n-treor': 'treor',
            'g': 'gsas', 'gsas': 'gsas', 'gsas2': 'gsas', 'gsas-ii': 'gsas', 'svd': 'gsas'}

def parse_methods(txt):
    if txt.strip().lower() in ('all', 'tout', 'toutes'):
        return [m for m, _ in _METHODS]
    out = []
    for t in re.split(r'[,; +]+', txt.strip().lower()):
        if t in _ALIASES and _ALIASES[t] not in out:
            out.append(_ALIASES[t])
        elif t:
            raise ValueError("methode inconnue : '%s' (choix : c, d, t, g)" % t)
    return out

def ask_methods(default='c,d,t'):
    if not sys.stdin.isatty():
        return parse_methods(default)
    print("Methodes disponibles : c = Conograph, d = dichotomie (type DICVOL), t = essai-erreur (type TREOR),\n"
          "                       g = GSAS-II natif (lent)")
    while True:
        try:
            txt = input("Methodes a lancer, separees par des virgules [Entree = %s] : " % default).strip()
        except (EOFError, KeyboardInterrupt):
            txt = ''
        try:
            return parse_methods(txt or default)
        except ValueError as e:
            print(e)

def _dispatch(task):
    """Point d'entree unique des processus paralleles."""
    if task[0] == 'gsas':
        ibrav, cells, dt, msg = _native_worker(task[1])
        name, system, cent, _ = _G2_BRAVAIS[ibrav]
        cands = [dict(engine='GSAS-II', label='', cell=c[3:9], bravais=name, system=system, cent=cent)
                 for c in cells]
        return 'gsas', name, cands, dt, msg
    engine, system, cands, dt, msg = _py_worker(task[1])
    return engine, system, cands, dt, msg

def main(argv=None):
    ap = argparse.ArgumentParser(description="Recherche de maille multi-methodes a partir d'un projet GSAS-II")
    ap.add_argument('gpx', nargs='?', help='projet GSAS-II (.gpx)')
    ap.add_argument('--input-json', help="donnees fournies par l'interface graphique (remplace le .gpx)")
    ap.add_argument('--json-out', help="ecrit aussi les resultats dans ce fichier JSON (pour l'interface graphique)")
    ap.add_argument('--methods', help="methodes a lancer : c (Conograph), d (dichotomie), t (essai-erreur), "
                                      "g (GSAS-II), ou 'all' ; ex. --methods c,d,t. Sans cette option, la question est posee.")
    ap.add_argument('--hist', type=int, default=0, help="numero de l'histogramme PWDR (defaut 0)")
    ap.add_argument('--maxpeaks', type=int, default=40, help='ne garder que les N premiers pics (defaut 40)')
    ap.add_argument('--zero', type=float, default=None, help='decalage de zero en degres (defaut : celui du projet)')
    ap.add_argument('--budget', type=float, default=90., help='duree maximale (s) de chaque recherche parallele (defaut 90)')
    ap.add_argument('--nproc', type=int, default=0, help='nombre de processus paralleles (defaut : tous les coeurs)')
    g = ap.add_argument_group('dichotomie et essai-erreur')
    g.add_argument('--nlines', type=int, default=20, help='nombre de raies utilisees (defaut 20)')
    g.add_argument('--tol', type=float, default=0.03, help='tolerance sur 2theta en degres (defaut 0.03)')
    g.add_argument('--amax', type=float, default=25., help='parametre de maille maximal en A (defaut 25)')
    g.add_argument('--vmax', type=float, default=4000., help='volume maximal en A3 (defaut 4000)')
    g.add_argument('--exhaustive', action='store_true',
                   help="ne pas arreter les systemes de basse symetrie quand une solution existe en symetrie plus haute")
    g = ap.add_argument_group('Conograph')
    g.add_argument('--conograph', help="chemin de l'executable Conograph (sinon recherche automatique)")
    g.add_argument('--search', choices=['quick', 'regular'], default='quick',
                   help='quick = rapide ; regular = plus pousse (basse symetrie difficile)')
    g.add_argument('--npeaks', type=int, default=0, help='nombre de pics transmis a Conograph (defaut : automatique)')
    g.add_argument('--threads', default='MAX', help='nombre de coeurs pour Conograph (MAX, DEFAULT ou un entier)')
    g.add_argument('--timeout', type=int, default=900, help='duree maximale de Conograph, en secondes')
    g.add_argument('--keep', help='dossier ou conserver les fichiers de travail de Conograph')
    g = ap.add_argument_group('GSAS-II')
    g.add_argument('--gsas2', help="dossier d'installation de GSAS-II (sinon recherche automatique)")
    g.add_argument('--npeaks-native', type=int, default=30, help='nombre de pics transmis au moteur GSAS-II (defaut 30)')
    g.add_argument('--no-triclinic', action='store_true', help='moteur GSAS-II : ne pas chercher en triclinique')
    g = ap.add_argument_group('sorties')
    g.add_argument('--select', help="mailles a transmettre a GSAS-II sans poser la question (ex. 1,3 ou a) ; '' = aucune")
    g.add_argument('--out-gpx', help='nom de la copie du projet recevant les mailles (defaut : <projet>_index.gpx)')
    g.add_argument('--out', help='fichier texte de resultat (defaut : <projet>_index.txt)')
    a = ap.parse_args(argv)

    if a.input_json:
        with open(a.input_json) as f:
            J = json.load(f)
        h = dict(name=J['name'], x=np.array(J['x'], float), y=np.array(J['y'], float), limits=None,
                 inst=dict((k, [v, v, False]) for k, v in J['inst'].items()),
                 peaks={'peaks': [[t, 0, i] for t, i in zip(J['tth'], J['inten'])]}, ipeaks=[[], []])
        phases = J.get('phases', [])
        a.gpx = a.gpx or J.get('project', 'projet.gpx')
        a.select = a.select if a.select is not None else ''
        a.methods = a.methods if a.methods is not None else 'c,d,t'
    elif a.gpx:
        hists, phases = read_gpx(a.gpx)
        hists = [hh for hh in hists if 'Autobkg' not in hh['name']] or hists
        if not hists:
            sys.exit('Aucun histogramme de poudre (PWDR) dans ce projet')
        h = hists[a.hist]
    else:
        ap.error('indiquez un projet .gpx')
    inst = h['inst']
    if not str(inst_val(inst, 'Type', 'PXC')).startswith('P') or 'C' not in str(inst_val(inst, 'Type', 'PXC'))[2:3]:
        sys.exit("Cette version ne traite que les donnees a longueur d'onde constante")
    lam = get_wavelength(inst)
    zero = float(inst_val(inst, 'Zero', 0.)) if a.zero is None else a.zero
    tth, inten, origin = get_peaks(h)
    if a.input_json:
        origin = J.get('origin', origin)
    o = np.argsort(tth); tth, inten = tth[o], inten[o]
    tth, inten = tth[:a.maxpeaks], inten[:a.maxpeaks]
    rep = []
    def say(s=''):
        print(s); rep.append(s)
    say('Projet      : %s' % os.path.basename(a.gpx))
    say('Histogramme : %s' % h['name'])
    say("Pics        : %d (%s), de %.3f a %.3f deg ; lambda = %.5f A ; zero = %.4f deg" % (
        len(tth), origin, tth[0], tth[-1], lam, zero))
    try:
        methods = parse_methods(a.methods) if a.methods is not None else ask_methods()
    except ValueError as e:
        sys.exit(str(e))
    if not methods:
        sys.exit('Aucune methode choisie')
    say('Methodes    : %s' % ', '.join(dict(_METHODS)[m] for m in methods))

    ref = None
    if phases:
        p = phases[0]
        ref = (p['cell'][:6], p['sg'].strip()[0] if p['sg'].strip()[0] in _CENT else 'P')
    t_start = time.time()
    by_method = dict((m, dict(cands=[], dt=0., msgs=[], ran=False)) for m in methods)

    # --- taches paralleles : dichotomie / essai-erreur (une par systeme cristallin), GSAS-II (une par reseau)
    opt = dict(nlines=a.nlines, tol=a.tol, amax=a.amax, vmax=a.vmax, nimp=0, budget=a.budget,
               early=not a.exhaustive)
    tasks = []
    for m in ('dichotomie', 'treor'):
        if m in methods:
            by_method[m]['ran'] = True
            for system in _PY_SYSTEMS:                    # hautes symetries d'abord (arret anticipe)
                tasks.append((m, (m, system, list(tth), lam, zero, opt)))
    if 'gsas' in methods:
        g2dir = find_gsas2(a.gsas2)
        if g2dir:
            by_method['gsas']['ran'] = True
            n_nat = min(len(tth), a.npeaks_native)
            for ib in range(17 if not a.no_triclinic else 16, -1, -1):
                tasks.append(('gsas', (g2dir, ib, list(tth[:n_nat]), list(inten[:n_nat]), lam, zero,
                                       max(3., min(a.budget, 60.) * _G2_BRAVAIS[ib][3]))))
        else:
            by_method['gsas']['msgs'].append('installation introuvable (option --gsas2 ou variable GSAS2DIR)')
    pool = res = None
    if tasks:
        shared = multiprocessing.Array('i', 2)
        pool = multiprocessing.Pool(a.nproc or max(1, min(len(tasks), multiprocessing.cpu_count())),
                                    initializer=_init_pool, initargs=(shared,))
        res = pool.map_async(_dispatch, tasks, chunksize=1)

    # --- Conograph, pendant que les autres tournent
    if 'conograph' in methods:
        exe = find_conograph(a.conograph)
        if exe:
            by_method['conograph']['ran'] = True
            try:
                xml, dt, wd = run_conograph(exe, h, tth, inten, lam, zero, a.search, a.npeaks or None,
                                            a.threads, a.timeout, a.keep)
                by_method['conograph']['cands'] = parse_conograph_xml(xml)
                by_method['conograph']['dt'] = dt
                if not a.keep:
                    shutil.rmtree(wd, ignore_errors=True)
            except Exception as e:
                by_method['conograph']['msgs'].append('echec : %s' % str(e).strip().splitlines()[-1])
        else:
            by_method['conograph']['msgs'].append('executable introuvable (option --conograph ou variable CONOGRAPH)')

    if res is not None:
        try:
            out = res.get(timeout=a.budget * 4 + 300)
            pool.close(); pool.join()
        except multiprocessing.TimeoutError:
            pool.terminate(); out = []
            say('Recherches paralleles interrompues (delai depasse).')
        for kind, what, cands, dt, msg in out:
            d = by_method[kind]
            if kind in ('dichotomie', 'treor'):
                cands = [assign_centering(c, tth, lam, zero) for c in cands]
            d['cands'] += cands
            d['dt'] = max(d['dt'], dt)                   # taches paralleles : duree = la plus longue
            if msg:
                d['msgs'].append('%s : %s' % (what, msg))

    # --- classement par methode, puis classement general (mailles equivalentes fusionnees)
    all_cands = []
    for m in methods:
        d = by_method[m]
        d['ranked'] = rank_candidates(d['cands'], tth, lam, zero) if d['cands'] else []
        all_cands += d['cands']
    if not all_cands:
        for m in methods:
            for msg in by_method[m]['msgs']:
                say('%-12s : %s' % (dict(_METHODS)[m], msg))
        say('\nAucune maille candidate.')
        if a.json_out:
            with open(a.json_out, 'w') as f:
                json.dump(dict(report='\n'.join(rep), ranked=[], lam=lam, tth=[], inten=[], table=[]), f)
        return 1
    ranked = rank_candidates(all_cands, tth, lam, zero)

    say('Duree totale : %.1f s' % (time.time() - t_start))
    if phases:
        say('Phase du projet (pour comparaison) : %s %s  %s  V=%.2f' % (
            p['name'], p['sg'], ' '.join('%.4f' % v for v in p['cell'][:6]), p['cell'][6]))
    say('\nCOMPARAISON DES METHODES (meilleure maille de chacune)')
    say('%-12s %6s %5s  %-15s %8s %8s %8s %7s %7s %7s %9s %6s %7s  %s' % (
        'methode', 'duree', 'cand.', 'reseau', 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'V(A3)', 'M20', 'non-idx',
        'rang general'))
    for m in methods:
        d = by_method[m]; name = dict(_METHODS)[m]
        if not d['ranked']:
            why = ' ; '.join(d['msgs']) if d['msgs'] else ('aucune maille trouvee' if d['ran'] else 'non lancee')
            say('%-12s %5.1fs %5d  %s' % (name, d['dt'], 0, why))
            continue
        bst = d['ranked'][0]; c = bst['cell_ref']
        pos = [i for i, r in enumerate(ranked) if abs(r['vprim'] - bst['vprim']) < 0.04 * r['vprim']
               and _same_sig(r['sig'], bst['sig'])]
        say('%-12s %5.1fs %5d  %-15s %8.4f %8.4f %8.4f %7.3f %7.3f %7.3f %9.2f %6.1f %3d/%-3d  %s' % (
            name, d['dt'], len(d['cands']), bst['bravais'], c[0], c[1], c[2], c[3], c[4], c[5], bst['vol'],
            bst['M20'], bst['nunidx_all'], bst['nobs'], ('n. %d' % (pos[0] + 1)) if pos else '-'))
        for msg in d['msgs']:
            say('%-12s        %s' % ('', msg))
    say("(M20, F20 et pics non indexes sont recalcules par ce script de la meme facon pour toutes les methodes.)")

    say('\nCLASSEMENT GENERAL (mailles equivalentes fusionnees)')
    say(format_table(ranked, 12, ref))
    best = ranked[0]
    say('\nMeilleure proposition : %s  a=%.4f b=%.4f c=%.4f alpha=%.3f beta=%.3f gamma=%.3f  V=%.2f A3  zero=%.4f deg' % (
        (best['bravais'],) + tuple(best['cell_ref']) + (best['vol'], best['zero'])))
    nran = sum(1 for m in methods if by_method[m]['ranked'])
    if len(best['engines']) > 1:
        say('Elle est trouvee independamment par %d methodes sur %d : %s.' % (
            len(best['engines']), nran, ', '.join(sorted(best['engines']))))
    elif nran > 1:
        if best['system'] == 'triclinique':
            say("Elle n'est trouvee que par %s : c'est attendu, la dichotomie et l'essai-erreur ne cherchent pas "
                "en triclinique." % ', '.join(best['engines']))
        else:
            say("Attention : elle n'est trouvee que par une seule methode (%s)." % ', '.join(best['engines']))
    others = [r for r in ranked[1:] if not r.get('equiv_of')]
    if others and others[0]['score'] > 0.6 * best['score']:
        say("Attention : un autre reseau (n. %d) a un score proche ; verifiez les deux dans GSAS-II."
            % (ranked.index(others[0]) + 1))
    eq = [r for r in ranked if r.get('equiv_of') == 1 and r['M20'] >= 0.9 * best['M20']
          and r['nunidx_all'] <= best['nunidx_all']]
    if eq:
        say("ATTENTION, pseudo-symetrie probable : le meme reseau decrit en symetrie plus basse (%s) indexe les pics "
            "aussi bien,\nalors qu'il devrait etre moins bien note s'il etait vraiment %s. La vraie maille est "
            "peut-etre celle de plus basse symetrie.\nLes positions de pics ne suffisent pas a trancher : comparez "
            "les deux par un Le Bail (ou verifiez avec la maille attendue)."
            % (', '.join('n. %d %s' % (ranked.index(r) + 1, r['bravais']) for r in eq), best['bravais']))
    if best['M20'] < 10:
        say("Attention : M20 < 10, l'indexation n'est pas fiable (pics imprecis, impurete ou mauvais zero ?).")
    un = unindexed_peaks(best['cell_ref'], best['cent'], tth, lam, best['zero'])
    if un:
        say('Pics non indexes par la meilleure proposition (2theta, ecart a la raie la plus proche) : '
            + ', '.join('%.3f (%.3f)' % u for u in un))
    out = a.out or os.path.splitext(a.gpx)[0] + '_index.txt'
    try:
        with open(out, 'w') as f:
            f.write('\n'.join(rep) + '\n')
        print('\nResultat enregistre dans : %s' % out)
    except OSError as e:
        print("\nImpossible d'ecrire %s (%s)" % (out, e))
    if a.json_out:
        keys = ('bravais', 'system', 'cent', 'M20', 'F20', 'nunidx20', 'nunidx_all', 'nobs', 'vol', 'zero')
        def _row(d):
            r = dict((k, (float(d[k]) if isinstance(d[k], (float, np.floating)) else
                          int(d[k]) if isinstance(d[k], (int, np.integer)) else d[k])) for k in keys)
            r['cell_ref'] = [float(v) for v in d['cell_ref'][:6]]
            r['engines'] = sorted(d['engines'])
            return r
        with open(a.json_out, 'w') as f:
            json.dump(dict(report='\n'.join(rep), ranked=[_row(d) for d in ranked[:20]], lam=lam,
                           tth=[float(t) for t in tth], inten=[float(i) for i in inten],
                           table=format_table(ranked, 20, ref).splitlines()), f)
    # --- validation par l'utilisateur, puis transmission a GSAS-II
    if a.select is None:
        sel = ask_selection(min(len(ranked), 12))
    else:
        sel = parse_selection(a.select.strip().lower(), min(len(ranked), 12))
    if sel:
        g2dir = find_gsas2(a.gsas2)
        if not g2dir:
            print("GSAS-II introuvable (option --gsas2) : impossible d'ecrire dans le projet.")
            return 0
        out_gpx = a.out_gpx or os.path.splitext(a.gpx)[0] + '_index.gpx'
        try:
            write_cells_to_gpx(a.gpx, out_gpx, h['name'], [ranked[i] for i in sel], tth, inten, lam, g2dir)
            print("\n%d maille(s) ecrite(s) dans une copie du projet : %s" % (len(sel), out_gpx))
            print("Dans GSAS-II : ouvrez cette copie, cliquez sur 'Unit Cells List' de l'histogramme ; chaque ligne du\n"
                  "tableau affiche ses raies calculees sur le diagramme. Puis 'Cell Index/Refine' > 'Refine Cell' et\n"
                  "'Make new phase' pour lancer le Le Bail. Le projet d'origine n'a pas ete modifie.")
        except Exception as e:
            print("Ecriture dans le projet impossible : %s" % e)
    return 0

if __name__ == '__main__':
    sys.exit(main())
