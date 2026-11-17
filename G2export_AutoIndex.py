# -*- coding: utf-8 -*-
'''Entree de menu GSAS-II : Export -> Entire project as -> (user) Auto-index (find unit cell)

Interface graphique du module g2_autoindex.py (recherche de maille multi-methodes).
A placer, avec g2_autoindex.py, dans ~/.GSASII/exports/ puis redemarrer GSAS-II.

Deroulement :
  1. fenetre d'options : histogramme, methodes (Conograph, dichotomie, essai-erreur, GSAS-II natif), reglages ;
  2. le calcul tourne dans un processus separe (GSAS-II reste utilisable, bouton d'annulation) ;
  3. fenetre de resultats : comparaison des methodes, classement general, cases a cocher ;
  4. les mailles cochees sont ecrites dans la "Unit Cells List" de l'histogramme du projet OUVERT
     (rien n'est enregistre sur disque tant que vous ne sauvegardez pas le projet).
'''
from __future__ import division, print_function
import os
import sys
import json
import signal
import tempfile
import subprocess
import importlib.util
import numpy as np
import wx
from .. import GSASIIfiles as G2fil
from .. import GSASIIdataGUI as G2gd

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPT = os.path.join(_HERE, 'g2_autoindex.py')
_PREFS = os.path.join(_HERE, 'g2_autoindex_prefs.json')
_MODULE = None


def _engine():
    '''Charge g2_autoindex.py (une seule fois) pour reutiliser ses fonctions.'''
    global _MODULE
    if _MODULE is None:
        spec = importlib.util.spec_from_file_location('g2_autoindex', _SCRIPT)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _MODULE = mod
    return _MODULE


def _python_exe():
    '''Python a utiliser pour le calcul : celui de GSAS-II, version "console" si elle existe.'''
    for cand in (os.path.join(sys.prefix, 'bin', 'python'), os.path.join(sys.prefix, 'python.exe'), sys.executable):
        if cand and os.path.isfile(cand):
            return cand
    return sys.executable


def _load_prefs():
    try:
        with open(_PREFS) as f:
            return json.load(f)
    except Exception:
        return {}


def _save_prefs(p):
    try:
        with open(_PREFS, 'w') as f:
            json.dump(p, f, indent=1)
    except Exception:
        pass


def _plain(v):
    '''Valeur simple (nombre ou texte) pour le JSON.'''
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return float(v)
    if isinstance(v, str):
        return v
    return None


class OptionsDialog(wx.Dialog):
    '''Choix de l'histogramme, des methodes et des reglages.'''
    def __init__(self, parent, histnames, info, prefs, cono_path):
        wx.Dialog.__init__(self, parent, title='Recherche de maille (auto-indexation)',
                           style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self.info = info
        main = wx.BoxSizer(wx.VERTICAL)

        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(wx.StaticText(self, label='Histogramme : '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hist = wx.Choice(self, choices=histnames)
        self.hist.SetSelection(0)
        self.hist.Bind(wx.EVT_CHOICE, self.OnHist)
        row.Add(self.hist, 1, wx.EXPAND)
        main.Add(row, 0, wx.EXPAND | wx.ALL, 8)
        self.summary = wx.StaticText(self, label='')
        main.Add(self.summary, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)

        box = wx.StaticBoxSizer(wx.StaticBox(self, label='Methodes'), wx.VERTICAL)
        self.cb = {}
        labels = [('c', 'Conograph  (tous systemes, rapide et robuste)'),
                  ('d', 'Dichotomie, type DICVOL  (cubique a monoclinique)'),
                  ('t', 'Essai-erreur, type TREOR  (cubique a monoclinique)'),
                  ('g', 'GSAS-II natif  (tous systemes, lent)')]
        chosen = prefs.get('methods', 'c,d,t').split(',')
        for key, lab in labels:
            c = wx.CheckBox(self, label=lab)
            c.SetValue(key in chosen)
            self.cb[key] = c
            box.Add(c, 0, wx.ALL, 3)
        if not cono_path:
            self.cb['c'].SetValue(False)
            self.cb['c'].SetLabel('Conograph  (executable introuvable : indiquez son chemin ci-dessous)')
        main.Add(box, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 8)

        grid = wx.FlexGridSizer(0, 2, 5, 8)
        grid.AddGrowableCol(1, 1)

        def add(label, ctrl):
            grid.Add(wx.StaticText(self, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(ctrl, 1, wx.EXPAND)
        self.maxpeaks = wx.SpinCtrl(self, min=8, max=200, initial=int(prefs.get('maxpeaks', 40)))
        add('Nombre maximal de pics utilises :', self.maxpeaks)
        self.zero = wx.TextCtrl(self, value='')
        add('Decalage de zero (deg) :', self.zero)
        self.tol = wx.TextCtrl(self, value=str(prefs.get('tol', 0.03)))
        add('Tolerance sur 2theta, dichotomie / essai-erreur (deg) :', self.tol)
        self.budget = wx.SpinCtrl(self, min=10, max=3600, initial=int(prefs.get('budget', 90)))
        add('Duree maximale de chaque recherche (s) :', self.budget)
        self.search = wx.Choice(self, choices=['rapide (quick)', 'poussee (regular)'])
        self.search.SetSelection(1 if prefs.get('search') == 'regular' else 0)
        add('Recherche Conograph :', self.search)
        self.exhaustive = wx.CheckBox(self, label='chercher aussi en basse symetrie quand une solution existe plus haut')
        self.exhaustive.SetValue(bool(prefs.get('exhaustive', False)))
        add('Dichotomie / essai-erreur :', self.exhaustive)
        self.cono = wx.TextCtrl(self, value=cono_path or prefs.get('conograph', ''))
        add('Executable Conograph :', self.cono)
        main.Add(grid, 0, wx.EXPAND | wx.ALL, 8)

        btns = wx.StdDialogButtonSizer()
        ok = wx.Button(self, wx.ID_OK, 'Lancer la recherche')
        ok.SetDefault()
        btns.AddButton(ok)
        btns.AddButton(wx.Button(self, wx.ID_CANCEL, 'Annuler'))
        btns.Realize()
        main.Add(btns, 0, wx.EXPAND | wx.ALL, 8)
        self.SetSizer(main)
        self.OnHist(None)
        main.Fit(self)
        self.SetMinSize(self.GetSize())
        self.CentreOnParent()

    def OnHist(self, event):
        i = self.info[self.hist.GetSelection()]
        self.summary.SetLabel('%d pics (%s)   lambda = %.5f A' % (i['npeaks'], i['origin'], i['lam']))
        self.zero.SetValue('%.4f' % i['zero'])

    def values(self):
        def num(ctrl, default):
            try:
                return float(ctrl.GetValue().replace(',', '.'))
            except ValueError:
                return default
        return dict(hist=self.hist.GetSelection(),
                    methods=','.join(k for k in 'cdtg' if self.cb[k].GetValue()),
                    maxpeaks=self.maxpeaks.GetValue(), zero=num(self.zero, 0.), tol=num(self.tol, 0.03),
                    budget=self.budget.GetValue(),
                    search='regular' if self.search.GetSelection() == 1 else 'quick',
                    exhaustive=self.exhaustive.GetValue(), conograph=self.cono.GetValue().strip())


class ResultsDialog(wx.Dialog):
    '''Comparaison des methodes + choix des mailles a transmettre.'''
    def __init__(self, parent, result):
        wx.Dialog.__init__(self, parent, title='Recherche de maille : resultats',
                           style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER | wx.MAXIMIZE_BOX)
        mono = wx.Font(11, wx.FONTFAMILY_TELETYPE, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
        main = wx.BoxSizer(wx.VERTICAL)
        main.Add(wx.StaticText(self, label='Compte rendu (comparaison des methodes et classement general) :'),
                 0, wx.ALL, 6)
        txt = wx.TextCtrl(self, value=result.get('report', ''),
                          style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_DONTWRAP | wx.HSCROLL)
        txt.SetFont(mono)
        main.Add(txt, 3, wx.EXPAND | wx.LEFT | wx.RIGHT, 6)
        table = result.get('table', [])
        self.rows = table[1:] if len(table) > 1 else []
        if self.rows:
            main.Add(wx.StaticText(self, label='Cochez les mailles a envoyer dans la "Unit Cells List" de GSAS-II :'),
                     0, wx.ALL, 6)
            head = wx.StaticText(self, label='   ' + table[0])
            head.SetFont(mono)
            main.Add(head, 0, wx.LEFT | wx.RIGHT, 6)
            self.clb = wx.CheckListBox(self, choices=self.rows)
            self.clb.SetFont(mono)
            self.clb.Check(0, True)
            main.Add(self.clb, 2, wx.EXPAND | wx.LEFT | wx.RIGHT, 6)
        else:
            self.clb = None
        btns = wx.BoxSizer(wx.HORIZONTAL)
        if self.rows:
            send = wx.Button(self, wx.ID_OK, 'Envoyer les mailles cochees vers Unit Cells List')
            send.SetDefault()
            btns.Add(send, 0, wx.ALL, 6)
        btns.AddStretchSpacer()
        btns.Add(wx.Button(self, wx.ID_CANCEL, 'Fermer'), 0, wx.ALL, 6)
        main.Add(btns, 0, wx.EXPAND)
        self.SetSizer(main)
        self.SetSize((1180, 760))
        self.CentreOnParent()

    def checked(self):
        if self.clb is None:
            return []
        return [i for i in range(len(self.rows)) if self.clb.IsChecked(i)]


class ExportAutoIndex(G2fil.ExportBaseclass):
    '''Recherche de maille multi-methodes, lancee depuis le menu Export.

    :param wx.Frame G2frame: fenetre principale de GSAS-II
    '''
    def __init__(self, G2frame):
        super(self.__class__, self).__init__(
            G2frame=G2frame,
            formatName='Auto-index (find unit cell)',
            extension='.txt',
            longFormatName='Find the unit cell with several indexing methods (Conograph, dichotomy, trial-and-error)'
            )
        self.exporttype = ['project']

    # ------------------------------------------------------------------ lecture de l'arbre
    def _child(self, parent, name):
        item = G2gd.GetGPXtreeItemId(self.G2frame, parent, name)
        return item if item else None

    def _hist_data(self, name):
        '''Rassemble, pour un histogramme PWDR, ce dont le moteur a besoin.'''
        G2frame = self.G2frame
        eng = _engine()
        hid = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, name)
        data = G2frame.GPXtree.GetItemPyData(hid)
        x = np.ma.getdata(np.asarray(data[1][0], float))
        y = np.ma.getdata(np.asarray(data[1][1], float))
        inst = G2frame.GPXtree.GetItemPyData(self._child(hid, 'Instrument Parameters'))[0]
        pk_id = self._child(hid, 'Peak List')
        ip_id = self._child(hid, 'Index Peak List')
        peaks = G2frame.GPXtree.GetItemPyData(pk_id) if pk_id else {}
        ipeaks = G2frame.GPXtree.GetItemPyData(ip_id) if ip_id else [[], []]
        hist = dict(name=name, x=x, y=y, limits=None, inst=inst, peaks=peaks or {}, ipeaks=ipeaks or [[], []])
        tth, inten, origin = eng.get_peaks(hist)
        return dict(hid=hid, hist=hist, tth=tth, inten=inten, origin=origin,
                    lam=eng.get_wavelength(inst), zero=float(eng.inst_val(inst, 'Zero', 0.)))

    def _phases(self):
        out = []
        try:
            G2frame = self.G2frame
            pid = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, 'Phases')
            if not pid:
                return out
            item, cookie = G2frame.GPXtree.GetFirstChild(pid)
            while item:
                gen = G2frame.GPXtree.GetItemPyData(item)['General']
                out.append(dict(name=G2frame.GPXtree.GetItemText(item), sg=gen['SGData']['SpGrp'],
                                cell=[float(v) for v in gen['Cell'][1:8]]))
                item, cookie = G2frame.GPXtree.GetNextChild(pid, cookie)
        except Exception:
            pass
        return out

    # ------------------------------------------------------------------ calcul dans un processus separe
    def _run(self, cmd, logfile):
        '''Lance le calcul et attend sa fin en gardant l'interface vivante. Retourne False si annule.'''
        kw = {}
        if os.name == 'posix':
            kw['start_new_session'] = True                # pour pouvoir arreter aussi les processus enfants
        with open(logfile, 'w') as log:
            proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, **kw)
            dlg = wx.ProgressDialog('Recherche de maille', 'Recherche en cours...\n(les methodes tournent en parallele)',
                                    parent=self.G2frame,
                                    style=wx.PD_APP_MODAL | wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME)
            ok = True
            try:
                while proc.poll() is None:
                    wx.MilliSleep(250)
                    res = dlg.Pulse()
                    cont = res[0] if isinstance(res, tuple) else res
                    if not cont:
                        ok = False
                        try:
                            if os.name == 'posix':
                                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                            else:
                                proc.terminate()
                        except Exception:
                            pass
                        break
            finally:
                dlg.Destroy()
            try:
                proc.wait(timeout=10)
            except Exception:
                pass
        return ok

    # ------------------------------------------------------------------ point d'entree du menu
    def Exporter(self, event=None):
        G2frame = self.G2frame
        try:
            self._exporter()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            wx.MessageBox('Erreur dans la recherche de maille :\n%s' % exc, 'Auto-index', wx.ICON_ERROR, G2frame)

    def _exporter(self):
        G2frame = self.G2frame
        if not os.path.isfile(_SCRIPT):
            wx.MessageBox('Fichier introuvable :\n%s' % _SCRIPT, 'Auto-index', wx.ICON_ERROR, G2frame)
            return
        eng = _engine()
        names = [n for n in G2gd.GetGPXtreeDataNames(G2frame, ['PWDR']) if 'Autobkg' not in n]
        hd, problems = {}, []
        for n in names:
            try:
                hd[n] = self._hist_data(n)
            except Exception as exc:
                problems.append('%s : %s' % (n, exc))
        names = [n for n in names if n in hd]
        if not names:
            wx.MessageBox("Aucun histogramme de poudre avec au moins 5 pics.\n"
                          "Faites d'abord une recherche de pics (Peak List).\n\n" + '\n'.join(problems),
                          'Auto-index', wx.ICON_INFORMATION, G2frame)
            return
        prefs = _load_prefs()
        cono = eng.find_conograph(prefs.get('conograph') or None)
        info = [dict(npeaks=len(hd[n]['tth']), origin=hd[n]['origin'], lam=hd[n]['lam'], zero=hd[n]['zero'])
                for n in names]
        dlg = OptionsDialog(G2frame, names, info, prefs, cono)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            v = dlg.values()
        finally:
            dlg.Destroy()
        if not v['methods']:
            wx.MessageBox('Aucune methode cochee.', 'Auto-index', wx.ICON_INFORMATION, G2frame)
            return
        _save_prefs(dict((k, v[k]) for k in ('methods', 'maxpeaks', 'tol', 'budget', 'search', 'exhaustive', 'conograph')))
        name = names[v['hist']]
        d = hd[name]
        inst = d['hist']['inst']
        project = G2frame.GSASprojectfile or 'projet.gpx'

        wd = tempfile.mkdtemp(prefix='g2autoindex_')
        fin, fout, flog = [os.path.join(wd, f) for f in ('input.json', 'result.json', 'log.txt')]
        ftxt = os.path.join(wd, 'report.txt')
        jinst = {}
        for k, val in inst.items():
            try:
                p = _plain(val[1])
            except Exception:
                p = None
            if p is not None:
                jinst[k] = p
        with open(fin, 'w') as f:
            json.dump(dict(name=name, project=os.path.basename(project), x=d['hist']['x'].tolist(),
                           y=d['hist']['y'].tolist(), inst=jinst, tth=[float(t) for t in d['tth']],
                           inten=[float(i) for i in d['inten']], origin=d['origin'], phases=self._phases()), f)
        g2dir = os.path.dirname(os.path.dirname(os.path.abspath(G2fil.__file__)))
        cmd = [_python_exe(), _SCRIPT, '--input-json', fin, '--json-out', fout, '--out', ftxt, '--select', '',
               '--methods', v['methods'], '--maxpeaks', str(v['maxpeaks']), '--zero', repr(v['zero']),
               '--tol', repr(v['tol']), '--budget', str(v['budget']), '--search', v['search'], '--gsas2', g2dir]
        if v['exhaustive']:
            cmd.append('--exhaustive')
        if v['conograph']:
            cmd += ['--conograph', v['conograph']]
        if not self._run(cmd, flog):
            return
        if not os.path.isfile(fout):
            try:
                with open(flog) as f:
                    log = f.read()[-3000:]
            except Exception:
                log = ''
            wx.MessageBox("Le calcul n'a pas produit de resultat.\n\n" + log, 'Auto-index', wx.ICON_ERROR, G2frame)
            return
        with open(fout) as f:
            result = json.load(f)

        # compte rendu enregistre a cote du projet, comme en ligne de commande
        try:
            if G2frame.GSASprojectfile:
                with open(os.path.splitext(G2frame.GSASprojectfile)[0] + '_index.txt', 'w') as f:
                    f.write(result.get('report', '') + '\n')
        except Exception:
            pass

        dlg = ResultsDialog(G2frame, result)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            sel = dlg.checked()
        finally:
            dlg.Destroy()
        if not sel:
            return
        chosen = [result['ranked'][i] for i in sel]
        hid = d['hid']
        uc_id = self._child(hid, 'Unit Cells List')
        if not uc_id:
            wx.MessageBox("Item 'Unit Cells List' introuvable pour cet histogramme.", 'Auto-index', wx.ICON_ERROR, G2frame)
            return
        old = G2frame.GPXtree.GetItemPyData(uc_id)
        new = eng.build_unit_cells_data(old, chosen, max(result['tth']), result['lam'])
        G2frame.GPXtree.SetItemPyData(uc_id, new)
        ip_id = self._child(hid, 'Index Peak List')
        if ip_id:
            ipl = G2frame.GPXtree.GetItemPyData(ip_id)
            if not ipl or not len(ipl[0]):
                G2frame.GPXtree.SetItemPyData(ip_id, eng.build_index_peaks(result['tth'], result['inten'], result['lam']))
        # affiche la Unit Cells List (on passe par l'histogramme pour forcer le rafraichissement)
        try:
            G2frame.GPXtree.SelectItem(hid)
            wx.CallAfter(G2frame.GPXtree.SelectItem, uc_id)
        except Exception:
            pass
        wx.MessageBox("%d maille(s) ajoutee(s) a la 'Unit Cells List'.\n\n"
                      "Cliquez sur une ligne du tableau pour superposer ses raies calculees au diagramme,\n"
                      "puis utilisez 'Cell Index/Refine' > 'Refine Cell' et 'Make new phase'.\n\n"
                      "Le projet n'est pas enregistre automatiquement." % len(chosen),
                      'Auto-index', wx.ICON_INFORMATION, G2frame)
