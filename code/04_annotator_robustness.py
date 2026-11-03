"""Robustness of RQ2-RQ4 to annotator (Limitations): assignment balance, RQ2 permutations within annotator, RQ3 with annotator, RQ4 within annotator strata.
Run: python code/04_annotator_robustness.py"""
import pandas as pd, numpy as np, itertools, warnings, statsmodels.formula.api as smf
from scipy.stats import chi2, chi2_contingency, rankdata, t as tdist
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent)); from common import load, load_annotators, RESULTS, FIGURES
lab, reg = load(); reg = reg.set_index('sample_id')
am = load_annotators()
grp=am.annotator.map(lambda a: a if a in ('A1','A2','A3','A4') else ('Unknown' if a=='NR' else 'Other')); reg['ann']=grp.reindex(reg.index).fillna('Unknown')
print(reg.ann.value_counts().to_dict())
# balance of assignment
for v in ('model','morphology'): 
    ct=pd.crosstab(reg.ann,reg[v]); c,p,dfree,_=chi2_contingency(ct); print(f'annotator x {v}: chi2={c:.1f}, df={dfree}, p={p:.3f}')
print('annotator x prompted rooms (Kruskal-style means):',reg.groupby('ann').prompted_room_count.mean().round(2).to_dict())
reg['model']=pd.Categorical(reg.model,['Nano Banana Pro','GPT-Image-1','GPT-Image-2','Nano Banana 2']); reg['morphology']=pd.Categorical(reg.morphology,['Rectangular','Asymmetrical/Irregular','Composite/Sprawling','L-Shape'])
reg['ann']=pd.Categorical(reg.ann,['A1','A2','A3','A4','Other','Unknown'])
base='flaw_count ~ C(morphology)+C(model)+bedroom_count+floor_area+log_rooms'
F=smf.negativebinomial(base+'+C(ann)',reg).fit(disp=0,maxiter=500); R=smf.negativebinomial(base,reg).fit(disp=0,maxiter=500)
s=2*(F.llf-R.llf); df=int(round(F.df_model-R.df_model)); print(f'ANNOTATOR effect on errors per layout: LR chi2={s:.2f}, df={df}, p={chi2.sf(s,df):.4f}; alpha {R.params["alpha"]:.3f} -> {F.params["alpha"]:.3f}')
print(reg.groupby('ann').flaw_count.agg(['count','mean','std']).round(2))
print(np.exp(F.params.filter(like='ann')).round(2).to_dict())
# does adding annotator change RQ3 conclusions?
blk=smf.negativebinomial('flaw_count ~ C(model)+C(ann)',reg).fit(disp=0,maxiter=500); s2=2*(F.llf-blk.llf); print(f'RQ3 size/geometry block jointly, adjusting annotator: LR chi2={s2:.2f}, df=6, p={chi2.sf(s2,6):.4f}')
# RQ2: permute model within annotator strata
cnt=pd.crosstab(lab.sample_id,lab.category).reindex(reg.index).fillna(0).values; mod=reg.model.astype(str).values; an=reg.ann.astype(str).values
def stat(l): T=np.array([cnt[l==k].sum(0) for k in np.unique(mod)]); T=T[:,T.sum(0)>0]; return chi2_contingency(T,correction=False)[0]
obs=stat(mod); rng=np.random.default_rng(2026)
def sp(l):
    out=l.copy()
    for a in np.unique(an):
        idx=np.where(an==a)[0]; out[idx]=l[rng.permutation(idx)]
    return out
null=np.array([stat(sp(mod)) for _ in range(5000)]); print(f'RQ2 omnibus chi2={obs:.1f}; permutation within annotator strata p={(1+(null>=obs).sum())/5001:.4f}')
for a,b in (('GPT-Image-1','GPT-Image-2'),('Nano Banana Pro','Nano Banana 2')):
    m=np.isin(mod,[a,b]); cm=cnt[m]; lm=mod[m]; am_=an[m]
    def st2(l): T=np.array([cm[l==k].sum(0) for k in (a,b)]); T=T[:,T.sum(0)>0]; return chi2_contingency(T,correction=False)[0]
    o=st2(lm); nl=[]
    for _ in range(5000):
        l2=lm.copy()
        for g in np.unique(am_):
            ii=np.where(am_==g)[0]; l2[ii]=lm[rng.permutation(ii)]
        nl.append(st2(l2))
    print(f'  {a} vs {b}: chi2={o:.1f}, perm p={(1+(np.array(nl)>=o).sum())/5001:.4f}, Bonferroni(2)={min(1,2*(1+(np.array(nl)>=o).sum())/5001):.4f}')
# RQ4 omnibus within annotator x room-tertile strata
B=(pd.crosstab(lab.sample_id,lab.category).reindex(reg.index).fillna(0)>0).astype(int); cats=list(B.columns); X=B.values; n=len(X); iu=np.triu_indices(13,1)
mp=lambda M: np.nan_to_num(np.corrcoef(M.T))[iu].mean()
q=pd.qcut(reg.prompted_room_count,3,labels=False,duplicates='drop').values; strata=pd.Series(an).astype(str)+'_'+pd.Series(q).astype(str)
o=mp(X); nl=[]
for _ in range(5000):
    M=X.copy()
    for g in strata.unique():
        idx=np.where(strata.values==g)[0]
        for j in range(13): M[idx,j]=X[rng.permutation(idx),j]
    nl.append(mp(M))
print(f'RQ4 omnibus mean phi={o:.4f}; permutation within annotator x room-tertile strata p={(1+(np.array(nl)>=o).sum())/5001:.4f} (null mean {np.mean(nl):.4f})')
# RQ4 count-level partial Spearman with annotator dummies
cn=pd.DataFrame(cnt,index=reg.index,columns=cats); tot=cn.values.sum(1); rk=rankdata
Dm=np.column_stack([pd.get_dummies(reg.model,drop_first=True).astype(float).values,pd.get_dummies(reg.ann,drop_first=True).astype(float).values])
def part(a,b,extra=None):
    Z=np.column_stack([np.ones(n),rk(reg.log_rooms.values),Dm]+([extra] if extra is not None else []))
    ra,rb=rk(cn[a].values),rk(cn[b].values)
    ea=ra-Z@np.linalg.lstsq(Z,ra,rcond=None)[0]; eb=rb-Z@np.linalg.lstsq(Z,rb,rcond=None)[0]
    r=np.corrcoef(ea,eb)[0,1]; dfree=n-Z.shape[1]-2; return r,2*tdist.sf(abs(r*np.sqrt(dfree/(1-r*r))),dfree)
for label,fn in (('rooms+model+annotator',lambda a,b:part(a,b)),('rooms+model+annotator+total errors excl. A,B',lambda a,b:part(a,b,rk(tot-cn[a].values-cn[b].values)))):
    Q=pd.DataFrame([(a,b)+fn(a,b) for a,b in itertools.combinations(cats,2)],columns=['A','B','rho','p']); Q['q_BH']=multipletests(Q.p,method='fdr_bh')[1]
    print(f'{label}: nominal p<.05={(Q.p<.05).sum()}, BH q<.05={(Q.q_BH<.05).sum()}, positive rho={(Q.rho>0).sum()}/78')
    print(Q.sort_values('p').head(8).round(4).to_string(index=False))
    Q.sort_values('p').to_csv(str(RESULTS)+'/rq4_countlevel_annotator_'+('strict' if 'total' in label else 'basic')+'.csv',index=False)
