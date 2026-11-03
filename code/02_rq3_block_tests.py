"""RQ3 supplementary: likelihood-ratio tests for predictor blocks in the NB model (MLE alpha). Run: python code/02_rq3_block_tests.py"""
import sys, warnings, pandas as pd, statsmodels.formula.api as smf
from scipy.stats import chi2, spearmanr
warnings.filterwarnings('ignore')
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent)); from common import load, load_annotators, RESULTS, FIGURES
OUT=str(RESULTS); _, reg = load()
reg['model']=pd.Categorical(reg.model,['Nano Banana Pro','GPT-Image-1','GPT-Image-2','Nano Banana 2'])
reg['morphology']=pd.Categorical(reg.morphology,['Rectangular','Asymmetrical/Irregular','Composite/Sprawling','L-Shape'])
full='flaw_count ~ C(morphology)+C(model)+bedroom_count+floor_area+log_rooms'
fit=lambda f: smf.negativebinomial(f,reg).fit(disp=0,maxiter=500); F=fit(full); rows=[]
for label,f in (('morphology (3 df)',full.replace('C(morphology)+','')),('model cohort (3 df)',full.replace('+C(model)','')),
  ('all prompt size/geometry terms: morphology+bedrooms+area+log rooms (6 df)','flaw_count ~ C(model)'),
  ('log prompted rooms (1 df)',full.replace('+log_rooms','')),('floor area (1 df)',full.replace('+floor_area','')),('bedroom count (1 df)',full.replace('+bedroom_count',''))):
    m=fit(f); s=2*(F.llf-m.llf); df=int(round(F.df_model-m.df_model)); rows.append(dict(block=label,LR_chi2=s,df=df,p=chi2.sf(s,df)))
r=spearmanr(reg.prompted_room_count,reg.flaw_count); rows.append(dict(block='Spearman rho prompted rooms vs errors',LR_chi2=r.statistic,df=None,p=r.pvalue))
pd.DataFrame(rows).to_csv(OUT+'/rq3_block_tests.csv',index=False); print(pd.DataFrame(rows).round(4).to_string())
