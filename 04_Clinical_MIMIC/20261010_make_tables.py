# -*- coding: utf-8 -*-
r"""Table 1 / Table 2 generation and verification

Previously the two main-text tables were **manually transcribed from the results CSVs**, with no generation script
and no audit trail. This script **regenerates both tables from the source data** and compares them cell by cell
with the tables in the manuscript docx, reporting any differences.

Source data
  Table 1  D_bilirubin_phototherapy_robust.csv (7 models)
           + D_bilirubin_pairs.csv (pairs and positive count for the "never phototherapy (ICD-9)" subset)
  Table 2  C_transfusion_robust_results.csv (32 rows = 4 indices x 2 exposures x 4 models)

Output
  04_Results/05_Figure_Data/table1_source_values.csv
  04_Results/05_Figure_Data/table2_source_values.csv
  cell-by-cell verification results printed to the terminal

Run: python 20261010_make_tables.py
"""
import sys, os, csv, io, re
sys.stdout.reconfigure(encoding='utf-8')
from docx import Document

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.abspath(os.path.join(HERE, '..', '..'))


def res_dir(prefix, root):
    """Resolve a 04_Results sub-folder by its numeric prefix.

    The folders carry the original names in the packaged copy and English names in the
    public release; matching on the prefix keeps this script working with either layout.
    """
    if os.path.isdir(root):
        for name in sorted(os.listdir(root)):
            if name.startswith(prefix):
                return os.path.join(root, name)
    return os.path.join(root, prefix)


RES = os.path.join(PKG, '04_Results')
CL = res_dir('04_', RES)
F5 = res_dir('05_', RES)
DOCX = os.path.join(PKG, '01_Manuscript',
                    'Manuscript_IL6_hepcidin_bilirubin_triangulation_SR_v1.0_20261009.docx')

# The saved result tables predate the English release and store the model labels in the
# original form. Map them to the canonical English labels so the lookups below work against
# either the old files or freshly regenerated ones.
MODEL_ALIAS_D = {
    '主模型（阳性 vs 阴性）': 'Main model (positive vs negative)',
    '+ 调整前后时间间隔与基线值': '+ adjusted for interval and baseline value',
    '① 调整窗口内光疗': '(1) adjusted for phototherapy in window',
    '② 仅窗口内无光疗': '(2) no phototherapy in window only',
    '③ 仅该次住院从未光疗(ICD)': '(3) never phototherapy (ICD-9) only',
    '⑤ 确定病原(非CoNS)阳性 vs 阴性': '(5) definite pathogen (non-CoNS) positive vs negative',
    '⑥ 仅CoNS阳性 vs 阴性': '(6) CoNS-only positive vs negative',
}
MODEL_ALIAS_C = {
    '主（全部）': 'Main (all)',
    '① 调整输血': '(1) adjusted for transfusion',
    '② 仅未输血': '(2) never transfused only',
    '③ 仅输血前化验': '(3) pre-transfusion measurements only',
}

def load(f):
    return list(csv.DictReader(io.open(os.path.join(CL, f), encoding='utf-8-sig')))

D = load('D_bilirubin_phototherapy_robust.csv')
PAIRS = load('D_bilirubin_pairs.csv')
C = [r for r in csv.DictReader(io.open(os.path.join(CL, 'C_transfusion_robust_results.csv'), encoding='utf-8-sig'))]

for r in D:
    r['model'] = MODEL_ALIAS_D.get(r['model'], r['model'])
for r in C:
    r['model'] = MODEL_ALIAS_C.get(r['model'], r['model'])

def drow(label):
    for r in D:
        if r['model'] == label:
            return r
    raise KeyError(label)

def crow(outcome, exposure, model):
    for r in C:
        if r['outcome'] == outcome and r['exposure'] == exposure and r['model'] == model:
            return r
    raise KeyError((outcome, exposure, model))

def f2(x):
    """Two decimals, except for a non-zero bound that would otherwise print as 0.00.

    A bound of -0.0015 mg/dL printed as "-0.00" reads as "no effect" and would be
    misquoted; those keep four decimals so the manuscript wording is reproduced.
    """
    v = float(x)
    s = '%+.2f' % v
    return ('%+.4f' % v) if v != 0 and float(s) == 0 else s

def fp(p):
    p = float(p)
    return ('%.4f' % p) if p < 0.01 else ('%.3f' % p)

def fpm(p):          # use absolute values consistently in the table (strip the plus sign)
    return fp(p).lstrip('+')

# ---------------- Table 1 ----------------
ROWS1 = [
    ('Unadjusted',                                'Main model (positive vs negative)'),
    ('Adjusted for interval and baseline value',  '+ adjusted for interval and baseline value'),
    ('+ adjusted for phototherapy in window',     '(1) adjusted for phototherapy in window'),
    ('Restricted to no phototherapy in window',   '(2) no phototherapy in window only'),
    ('__NEVER_PHOTOTHERAPY__',                    None),      # see special handling below (n/case count taken from the paired data)
    ('Definite pathogen (non-coagulase-negative) vs negative', '(5) definite pathogen (non-CoNS) positive vs negative'),
    ('Coagulase-negative staphylococci only vs negative',      '(6) CoNS-only positive vs negative'),
]
t1 = []
for label_en, label_src in ROWS1:
    if label_src is None:      # subset "never phototherapy (ICD-9)": the model cannot be fitted, only report n and the case count
        sub = [x for x in PAIRS if x.get('any_ph_icd') == '0']
        nc = sum(1 for x in sub if x.get('positive') in ('1', 'True', 'true'))
        t1.append(['Restricted to admissions never receiving phototherapy (ICD-9)',
                   '%d (%d)' % (len(sub), nc), 'not estimable (too few cases)',
                   '\u2014', '\u2014', '0'])
        continue
    r = drow(label_src)
    t1.append([label_en, '%s (%s)' % (r['n'], r['ncase']), f2(r['b']),
               '%s to %s' % (f2(r['lo']), f2(r['hi'])), fpm(r['p']),
               r['estimable'] if 'estimable' in r else '1'])

# ---------------- Table 2 ----------------
ROWS2 = [
    ('MCV (fL)', 'mcv',     'pos_noncons', 'Main (all)',            'Main'),
    ('MCV (fL)', 'mcv',     'pos_noncons', '(1) adjusted for transfusion',            '+ adjusted for transfusion'),
    ('MCV (fL)', 'mcv',     'pos_noncons', '(2) never transfused only',            'Never transfused only'),
    ('MCV (fL)', 'mcv(preTx)', 'pos_noncons', '(3) pre-transfusion measurements only',      'Pre-transfusion measurements only'),
    ('MCH (pg)', 'mch',     'pos_noncons', 'Main (all)',            'Main'),
    ('MCH (pg)', 'mch',     'pos_noncons', '(1) adjusted for transfusion',            '+ adjusted for transfusion'),
    ('MCH (pg)', 'mch',     'pos_noncons', '(2) never transfused only',            'Never transfused only'),
    ('MCH (pg)', 'mch(preTx)', 'pos_noncons', '(3) pre-transfusion measurements only',      'Pre-transfusion measurements only'),
    ('MCV (fL)', 'mcv',     'pos_blood',   'Main (all)',            'Main'),
    ('MCV (fL)', 'mcv',     'pos_blood',   '(1) adjusted for transfusion',            '+ adjusted for transfusion'),
    ('MCH (pg)', 'mch',     'pos_blood',   'Main (all)',            'Main'),
    ('MCH (pg)', 'mch',     'pos_blood',   '(1) adjusted for transfusion',            '+ adjusted for transfusion'),
    ('RDW (%)',  'rdw',     'pos_noncons', 'Main (all)',            'Main'),
    ('RDW (%)',  'rdw(preTx)', 'pos_noncons', '(3) pre-transfusion measurements only',      'Pre-transfusion measurements only'),
    ('Hct (%)',  'hct',     'pos_blood',   'Main (all)',            'Main'),
    ('Hct (%)',  'hct',     'pos_blood',   '(1) adjusted for transfusion',            '+ adjusted for transfusion'),
]
t2 = []
for out_en, out_key, exp_key, model_src, model_en in ROWS2:
    r = crow(out_key, exp_key, model_src)
    t2.append([out_en, 'Definite pathogen' if exp_key == 'pos_noncons' else 'Any positive (incl. CoNS)',
               model_en, '%s (%s)' % (r['n'], r['ncase']),
               '%s (%s, %s)' % (f2(r['b']), f2(r['lo']), f2(r['hi'])), fpm(r['p'])])

# ---------------- Write out ----------------
with io.open(os.path.join(F5, 'table1_source_values.csv'), 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Model', 'n (positive)', 'Difference in delta-bilirubin (mg/dL)', '95% CI', 'p', 'estimable'])
    w.writerows(t1)
with io.open(os.path.join(F5, 'table2_source_values.csv'), 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Outcome', 'Exposure', 'Model', 'n (cases)', 'beta (95% CI)', 'p'])
    w.writerows(t2)

# ---------------- Cell-by-cell comparison with the manuscript tables ----------------
doc = Document(DOCX)
print('=' * 108)
print('Comparison with the manuscript docx tables')
print('=' * 108)
bad = 0
for ti, mine in enumerate([t1, t2], 1):
    tbl = doc.tables[ti - 1]
    ref = [[c.text.strip() for c in r.cells] for r in tbl.rows]
    print('\nTABLE %d: manuscript %d data rows | generated %d data rows' % (ti, len(ref) - 1, len(mine)))
    for i, row in enumerate(mine):
        if i + 1 >= len(ref):
            print('   WARNING: manuscript is missing this row:', row[0]); bad += 1; continue
        r = ref[i + 1]
        if ti == 1:
            pairs = [(row[1], r[1]), (row[2], r[2]), (row[3], r[3]), (row[4], r[4])]
        else:
            pairs = [(row[3], r[3]), (row[4], r[4]), (row[5], r[5])]
        for a, b in pairs:
            na, nb = re.findall(r'-?\d+\.?\d*', a.replace('\u2212', '-')), re.findall(r'-?\d+\.?\d*', b.replace('\u2212', '-'))
            ok = True
            if len(na) == len(nb):
                for x, y in zip(na, nb):
                    fx, fy = float(x), float(y)
                    tol = 5e-3 if max(abs(fy), abs(fx)) > 1 else 5e-4   # allow last-digit rounding
                    if abs(fx - fy) > max(tol, abs(fy) * 0.02):
                        ok = False
            elif a.strip() != b.strip():
                ok = False
            if not ok:
                print('   MISMATCH: %-46s generated=%s | manuscript=%s' % (row[0][:46], a, b)); bad += 1
    # rows present in the manuscript but not generated
    for r in ref[1 + len(mine):]:
        print('   WARNING: extra row in the manuscript:', r[0]); bad += 1
print()
print('differing cells:', bad, '->', 'all consistent' if bad == 0 else 'needs manual review')
print('Wrote table1_source_values.csv / table2_source_values.csv to 04_Results/05_Figure_Data/')
