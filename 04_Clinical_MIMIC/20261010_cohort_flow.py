# -*- coding: utf-8 -*-
r"""Cohort flow recomputation and audit trail  --  faithfully reproducing the definitions of 20261008_neo_infection_bilirubin*.py

MIMIC-IV step by step:
  (1) NICU/SCN admissions
  (2) of which those with >=1 total bilirubin
  (3) restricted to the "culture-positive" vs "no infection exposure" groups, with an APR-DRG weight band -> the fully adjusted model sample
     (model: log(bilirubin) ~ infection + postnatal age + weight band, HC3)
  (4) the structure of that sample (measurements / admissions / infants), and the same model re-estimated
     with standard errors clustered by admission and by infant -- an infant contributes on average 2.3
     measurements, so HC3 (robust to heteroscedasticity) does not by itself allow for the correlation
     between repeated measurements of the same infant
  (5) sex of the two compared groups

MIMIC-III step by step (from the packaged nicu_base.csv / nicu_culture.csv, no database needed)

Output: 04_Results/04_Clinical_MIMIC/cohort_flow.csv
"""
import sys, os, csv, io, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import psycopg2
import statsmodels.formula.api as smf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '04_Results'))


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


RES = res_dir('04_', ROOT)
rows = []

con = psycopg2.connect(host='localhost', port=5432, user='postgres', password='1234', dbname='mimiciv10')
con.autocommit = True
q = lambda s: pd.read_sql(s, con)

# (1) Cohort
adm = q("SELECT hadm_id, subject_id FROM work.nicu_adm WHERE in_nicu OR in_scn")
rows.append(('MIMIC-IV', 'NICU/SCN admissions', len(adm), 'work.nicu_adm WHERE in_nicu OR in_scn', 'PostgreSQL mimiciv10'))

# (2) Bilirubin records
bili = q("""SELECT hadm_id, valuenum AS bili, age_hours FROM work.neo_lab
            WHERE itemid=50885 AND valuenum>0 AND valuenum<45 AND age_hours IS NOT NULL""")
bili['age_hours'] = pd.to_numeric(bili.age_hours, errors='coerce')
bili['bili'] = pd.to_numeric(bili.bili, errors='coerce')
bili = bili.dropna(subset=['age_hours', 'bili'])
bili = bili[bili.age_hours >= 0]
rows.append(('MIMIC-IV', 'admissions with at least one total bilirubin',
             len(set(adm.hadm_id) & set(bili.hadm_id)),
             'cohort & neo_lab itemid=50885 (0<value<45)', 'PostgreSQL mimiciv10'))
# The main analysis (models A/B/C) is restricted to within 336 h (14 days) of birth, consistent with 20261008_neo_infection_bilirubin_align.py
bili = bili[bili.age_hours <= 336]
print('  [window] number of bilirubin records within 0-336 h = %d' % len(bili))

# Infection exposure
pos = set(q("""SELECT DISTINCT hadm_id FROM work.neo_cx
               WHERE spec_class='blood' AND culture_positive AND org_name IS NOT NULL""").hadm_id)
sent = set(q("SELECT DISTINCT hadm_id FROM work.neo_cx WHERE spec_class='blood'").hadm_id)

rx = q("""SELECT hadm_id, starttime::timestamp s, stoptime::timestamp e, route
          FROM work.neo_rx WHERE drug_class='antimicrobial'""")
SYS = ("IV", "IM", "PO", "PO/NG", "NG", "PO/OG", "OG", "SC", "SUBCUT",
       "IV DRIP", "IV BOLUS", "PIV", "IVPCA", "PB", "IV/PO", "ID")
rx = rx[rx.route.isin(SYS)].dropna(subset=['s', 'e'])
rx = rx[rx.e > rx.s]
abx = {}
for h, g in rx.groupby('hadm_id'):
    iv = g[['s', 'e']].sort_values('s').values
    mg, (cs, ce) = [], tuple(iv[0])
    for s, e in iv[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            mg.append((cs, ce)); cs, ce = s, e
    mg.append((cs, ce))
    abx[h] = sum((e - s) / np.timedelta64(1, 'h') for s, e in mg) / 24.0

df = adm[adm.hadm_id.isin(set(bili.hadm_id))].copy()
df['bc_pos'] = df.hadm_id.isin(pos)
df['bc_sent'] = df.hadm_id.isin(sent)
df['abx_days'] = df.hadm_id.map(abx)
b = bili.merge(df[['hadm_id', 'subject_id', 'bc_pos', 'bc_sent', 'abx_days']], on='hadm_id', how='inner')
b['group'] = np.where(b.bc_pos, 'pos',
                      np.where(~b.bc_sent & (b.abx_days.isna() | (b.abx_days == 0)), 'no_infection', 'other'))
n_pos = b[b.group == 'pos'].hadm_id.nunique()
n_noinf = b[b.group == 'no_infection'].hadm_id.nunique()
rows.append(('MIMIC-IV', 'admissions in the two compared groups (culture-positive or no infection exposure)',
             n_pos + n_noinf, 'culture-positive %d + no-infection-exposure %d' % (n_pos, n_noinf), 'PostgreSQL mimiciv10'))

# (3) Fully adjusted model
bw = q("""SELECT hadm_id, MAX(drg_code::text) AS apr_bw FROM hosp.drgcodes
          WHERE drg_type='APR' AND drg_code::text IN
          ('593','602','607','612','614','621','622','625','626','633','634','639','640')
          GROUP BY 1""")
BW_ORDER = {'593': 1, '602': 2, '607': 3, '612': 4, '614': 5, '621': 6, '622': 7,
            '625': 8, '626': 9, '633': 10, '634': 11, '639': 12, '640': 13}
g = b[b.group.isin(['pos', 'no_infection'])].copy()
m = g[g.bili > 0].merge(bw, on='hadm_id', how='left')
m['bw_rank'] = pd.to_numeric(m.apr_bw.map(BW_ORDER), errors='coerce')
m['logb'] = np.log(m.bili)
m['age_d'] = m.age_hours / 24.0
m['infect'] = (m.group == 'pos').astype(int)
for spec, lab, need in [('logb ~ infect', 'A infection only', ['logb', 'age_d']),
                        ('logb ~ infect + age_d', 'B + age at measurement', ['logb', 'age_d']),
                        ('logb ~ infect + age_d + bw_rank', 'C + age + APR-DRG weight band',
                         ['logb', 'age_d', 'bw_rank'])]:
    d = m.dropna(subset=need)
    r = smf.ols(spec, data=d).fit(cov_type='HC3')
    rows.append(('MIMIC-IV', 'model %s' % lab, int(r.nobs),
                 'log-scale beta %+.4f (SE %.4f, p=%.3g)' % (r.params['infect'], r.bse['infect'], r.pvalues['infect']),
                 'PostgreSQL mimiciv10'))
    print('  [model %s] n=%d  beta(log)=%+.4f  SE=%.4f  p=%.3g' % (
        lab, int(r.nobs), r.params['infect'], r.bse['infect'], r.pvalues['infect']))

# ---- (4) Sample structure and cluster-robust sensitivity for the fully adjusted model ----
# The fully adjusted model is fitted at the level of the individual measurement, so an infant
# contributes on average 2.3 records. HC3 standard errors are robust to heteroscedasticity but
# not to the correlation between repeated measurements of the same infant; the models below
# therefore re-estimate model C with standard errors clustered by admission and by infant.
SPEC_C = 'logb ~ infect + age_d + bw_rank'
dc = m.dropna(subset=['logb', 'age_d', 'bw_rank'])
n_adm, n_inf = int(dc.hadm_id.nunique()), int(dc.subject_id.nunique())
print()
print('  [sample] %d measurements from %d admissions (%d infants); %.2f measurements per admission'
      % (len(dc), n_adm, n_inf, len(dc) / n_adm))
rows.append(('MIMIC-IV', 'admissions contributing to the fully adjusted model', n_adm,
             '%d measurements from %d infants (%.2f measurements per admission)'
             % (len(dc), n_inf, len(dc) / n_adm), 'PostgreSQL mimiciv10'))

print('  [sensitivity] model C with standard errors clustered by admission / by infant')
for grp, lab in [('hadm_id', 'admission'), ('subject_id', 'infant')]:
    rc = smf.ols(SPEC_C, data=dc).fit(cov_type='cluster', cov_kwds={'groups': dc[grp]})
    ci = rc.conf_int().loc['infect']
    print('  [cluster by %-8s] n=%d  beta(log)=%+.4f  SE=%.4f  95%% CI %+.4f to %+.4f  p=%.3g'
          % (lab, int(rc.nobs), rc.params['infect'], rc.bse['infect'], ci[0], ci[1],
             rc.pvalues['infect']))
    rows.append(('MIMIC-IV', 'model C with standard errors clustered by %s' % lab, int(rc.nobs),
                 'log-scale beta %+.4f (SE %.4f, 95%% CI %+.4f to %+.4f, p=%.3g)'
                 % (rc.params['infect'], rc.bse['infect'], ci[0], ci[1], rc.pvalues['infect']),
                 'PostgreSQL mimiciv10'))

# ---- (5) Sex of the two compared groups ---------------------------------------------------
# MIMIC-IV ships the patient table as hosp.patients; this instance keeps it alongside the other
# core tables, so the schema is resolved from information_schema rather than hard-coded.
pat_tbl = q("""SELECT table_schema || '.' || table_name AS t FROM information_schema.tables
               WHERE table_name = 'patients'""").t.iloc[0]
pat = q('SELECT subject_id, gender FROM %s' % pat_tbl)
sa = adm.merge(pat, on='subject_id', how='left')
sa = sa[sa.hadm_id.isin(set(g.hadm_id))]
nM, nF = int((sa.gender == 'M').sum()), int((sa.gender == 'F').sum())
print('  [sex] two compared groups (n=%d): male %d (%.1f%%), female %d (%.1f%%)'
      % (len(sa), nM, 100.0 * nM / len(sa), nF, 100.0 * nF / len(sa)))
rows.append(('MIMIC-IV', 'sex of the two compared groups', len(sa),
             'male %d (%.1f%%), female %d (%.1f%%)' % (nM, 100.0 * nM / len(sa), nF, 100.0 * nF / len(sa)),
             '%s.gender' % pat_tbl))
con.close()

# ============================ MIMIC-III ============================
base = list(csv.DictReader(io.open(os.path.join(RES, 'nicu_base.csv'), encoding='utf-8-sig')))
cul = list(csv.DictReader(io.open(os.path.join(RES, 'nicu_culture.csv'), encoding='utf-8-sig')))
all_ids = {r['hadm_id'] for r in base}
keep = {r['hadm_id'] for r in base if r['age_days_at_icu'].strip() and float(r['age_days_at_icu']) <= 28}
c = [r for r in cul if r['hadm_id'] in keep]
posr = [r for r in c if (r.get('positive') or '').strip() in ('1', 'True', 'true')]
hpos = {r['hadm_id'] for r in posr}
noncons = {r['hadm_id'] for r in posr
           if (r.get('orgs') or '')
           and 'staphylococcus' not in (r.get('orgs') or '').lower()
           and 'coagulase negative' not in (r.get('orgs') or '').lower()}
rows += [
    ('MIMIC-III', 'NICU stays (first_careunit=NICU), before the age filter', len(all_ids), 'nicu_base.csv', 'nicu_base.csv'),
    ('MIMIC-III', 'admissions with age at ICU entry <= 28 d (analysis cohort)', len(keep), 'age_days_at_icu <= 28', 'nicu_base.csv'),
    ('MIMIC-III', 'with >=1 blood culture', len({r['hadm_id'] for r in c}), 'within the 7,940 cohort', 'nicu_culture.csv'),
    ('MIMIC-III', 'blood-culture episodes', len(c), 'admission x collection time x specimen type', 'nicu_culture.csv'),
    ('MIMIC-III', 'positive episodes', len(posr), '%.1f%% of episodes' % (100.0 * len(posr) / len(c)), 'nicu_culture.csv'),
    ('MIMIC-III', 'admissions with a positive culture', len(hpos), '%.1f%% of the 7,940 cohort' % (100.0 * len(hpos) / len(keep)), 'nicu_culture.csv'),
    ('MIMIC-III', 'admissions with a definite (non-CoNS) pathogen', len(noncons), 'CoNS excluded', 'nicu_culture.csv'),
]

out = os.path.join(RES, 'cohort_flow.csv')
with io.open(out, 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f)
    w.writerow(['layer', 'stage', 'n', 'note', 'source'])
    w.writerows(rows)

print()
print('=' * 104)
for r in rows:
    print('  %-10s %-62s %7s  %s' % (r[0], r[1][:62], r[2], r[3][:60]))
print('=' * 104)
print('Written:', out)
