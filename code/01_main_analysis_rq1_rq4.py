"""Re-run of the reduced statistical set (RQ1-RQ4). Does not modify the manuscript.
Inputs: data/annotations.csv, data/layouts.csv
Outputs: results/ (CSVs, co-occurrence heatmap, results_summary.md).  Seeds fixed."""
import sys, itertools, warnings, json
import numpy as np, pandas as pd, networkx as nx
import statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.stats.proportion import proportion_confint
from statsmodels.stats.multitest import multipletests
from scipy.stats import chi2_contingency, spearmanr
warnings.filterwarnings('ignore')
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent)); from common import load, load_annotators, RESULTS, FIGURES
OUT = str(RESULTS)
SEED=2026
lab, reg = load()
DOM={'Redundant / Unnecessary Elements':'Graphic Drafting','Symbol Errors':'Graphic Drafting','Line & Layer Errors':'Graphic Drafting','Scale Errors':'Graphic Drafting',
 'Missing Openings & Elements':'Spatial Circulation','Misaligned Openings & Elements':'Spatial Circulation','Blocked Pathways & Movement':'Spatial Circulation',
 'Poor Room Layout & Placement':'Ergonomic Layout','Dimension & Measurement Errors':'Ergonomic Layout','Space & Clearance Errors':'Ergonomic Layout','Missing Required Spaces':'Ergonomic Layout',
 'Overlapping / Colliding Elements':'Structural Logic','Unclear Wall Types':'Regulatory Compliance'}
lab['domain']=lab.category.map(DOM); assert lab.domain.notna().all()
reg=reg.set_index('sample_id'); ids=reg.index
cnt=pd.crosstab(lab.sample_id,lab.category).reindex(ids).fillna(0).astype(int)   # layout x category error counts
B=(cnt>0).astype(int); cats=list(cnt.columns); N=len(ids); NERR=int(cnt.values.sum())
summ=[]; S=summ.append
S(f'# Reduced statistical set: re-run results\n\nN layouts = {N}; N errors = {NERR}; categories = {len(cats)}. Seed = {SEED}.\n')

# ---------- RQ1: shares + layout prevalence (Wilson 95%) ----------
rows=[]
for level,col in (('domain','domain'),('category','category')):
    err=lab.groupby(col).size()
    pres=(pd.crosstab(lab.sample_id,lab[col]).reindex(ids).fillna(0)>0).sum()
    for k in err.index:
        lo,hi=proportion_confint(pres[k],N,method='wilson')
        rows.append(dict(level=level,name=k,errors=int(err[k]),error_share_pct=100*err[k]/NERR,layouts_affected=int(pres[k]),prevalence_pct=100*pres[k]/N,wilson_lo_pct=100*lo,wilson_hi_pct=100*hi))
R1=pd.DataFrame(rows).sort_values(['level','errors'],ascending=[True,False]); R1.to_csv(OUT+'/rq1_shares_prevalence.csv',index=False)
S('## RQ1. Error shares and layout prevalence (Wilson 95% CI)\n'); d=R1[R1.level=='domain']
S('| Domain | Errors | Share % | Layouts affected | Prevalence % [95% CI] |\n|---|---|---|---|---|')
for _,r in d.iterrows(): S(f'| {r["name"]} | {r.errors} | {r.error_share_pct:.2f} | {r.layouts_affected} | {r.prevalence_pct:.1f} [{r.wilson_lo_pct:.1f}, {r.wilson_hi_pct:.1f}] |')
S(f'\nLayouts with at least one error: {(cnt.sum(1)>0).sum()}/{N}. Category-level table: rq1_shares_prevalence.csv\n')

# ---------- RQ2: category x model chi2 + Cramer V, layout-level permutation ----------
model=reg.model.values; X=cnt.values; rng=np.random.default_rng(SEED)
def chi_V(mask,labels):
    g=np.unique(labels[mask]); T=np.array([X[mask][labels[mask]==k].sum(0) for k in g]); T=T[:,T.sum(0)>0]
    c=chi2_contingency(T,correction=False)[0]; return c,np.sqrt(c/(T.sum()*(min(T.shape)-1))),T
def perm(mask,P=10000):
    obs,V,T=chi_V(mask,model); lab_=model[mask]; Xm=X[mask]; gs=np.unique(lab_); null=np.empty(P)
    for i in range(P):
        pl=rng.permutation(lab_); Tn=np.array([Xm[pl==k].sum(0) for k in gs]); Tn=Tn[:,Tn.sum(0)>0]
        null[i]=chi2_contingency(Tn,correction=False)[0]
    return obs,V,T.shape,(1+(null>=obs).sum())/(P+1),np.median(null)
allm=np.ones(N,bool); o,V,sh,p,nm=perm(allm)
rows=[dict(test='omnibus (4 models)',chi2=o,V=V,df_naive=(sh[0]-1)*(sh[1]-1),perm_p=p,null_median_chi2=nm,bonferroni_p=np.nan)]
S('## RQ2. Category distribution by model (layout-level permutation, 10,000 permutations)\n')
S(f'Omnibus: chi2 = {o:.1f}, Cramer V = {V:.3f}, permutation p = {p:.4f} (min attainable {1/10001:.4f}); permutation-null median chi2 = {nm:.1f} versus naive df = {(sh[0]-1)*(sh[1]-1)}, i.e. error-level chi2 is inflated by clustering.\n')
S('Planned contrasts (Bonferroni, M = 2):\n\n| Contrast | chi2 | V | perm p | Bonferroni p |\n|---|---|---|---|---|')
for a,b in (('GPT-Image-1','GPT-Image-2'),('Nano Banana Pro','Nano Banana 2')):
    m=np.isin(model,[a,b]); o,V,sh,p,nm=perm(m); pb=min(1,2*p)
    rows.append(dict(test=f'{a} vs {b}',chi2=o,V=V,df_naive=(sh[0]-1)*(sh[1]-1),perm_p=p,null_median_chi2=nm,bonferroni_p=pb))
    S(f'| {a} vs {b} | {o:.1f} | {V:.3f} | {p:.4f} | {pb:.4f} |')
pd.DataFrame(rows).to_csv(OUT+'/rq2_model_category_tests.csv',index=False)
mt=reg.groupby('model').flaw_count.agg(['count','mean','std']).round(2); mt.to_csv(OUT+'/rq2_errors_per_layout_by_model.csv')
S('\nErrors per layout by model:\n\n'+mt.to_markdown()+'\n')

# ---------- RQ3: NB GLM with log(rooms) covariate (MLE alpha, HC1 robust SE) ----------
rho=spearmanr(reg.prompted_room_count,reg.flaw_count)
reg['model']=pd.Categorical(reg.model,['Nano Banana Pro','GPT-Image-1','GPT-Image-2','Nano Banana 2'])
reg['morphology']=pd.Categorical(reg.morphology,['Rectangular','Asymmetrical/Irregular','Composite/Sprawling','L-Shape'])
f='flaw_count ~ C(morphology) + C(model) + bedroom_count + floor_area + log_rooms'
nb=smf.negativebinomial(f,reg).fit(disp=0,maxiter=500,cov_type='HC1')
ci=nb.conf_int(); T3=pd.DataFrame({'term':nb.params.index,'coef':nb.params.values,'IRR':np.exp(nb.params.values),'IRR_lo':np.exp(ci[0].values),'IRR_hi':np.exp(ci[1].values),'p':nb.pvalues.values})
T3=T3[T3.term!='Intercept'].copy(); T3['term']=T3.term.str.replace(r'C\((\w+)\)\[T\.(.*)\]',r'\2',regex=True); T3.loc[T3.term=='alpha','term']='alpha (NB dispersion, MLE)'
T3.to_csv(OUT+'/rq3_negbin_logrooms_covariate.csv',index=False)
# sensitivity: log prompted rooms as an offset instead of a covariate
nb0=smf.negativebinomial(f.replace(' + log_rooms',''),reg,offset=reg.log_rooms).fit(disp=0,maxiter=500,cov_type='HC1')
S('## RQ3. Negative binomial GLM with log(rooms) as a covariate (reference: Rectangular, Nano Banana Pro; robust HC1 SE)\n')
S('| Term | IRR | 95% CI | p |\n|---|---|---|---|')
for _,r in T3.iterrows():
    if r.term.startswith('alpha'): continue
    S(f'| {r.term} | {r.IRR:.3f} | [{r.IRR_lo:.3f}, {r.IRR_hi:.3f}] | {r.p:.3f} |')
lr=T3[T3.term=='log_rooms'].iloc[0]
S(f'\nlog(rooms) coefficient (elasticity) = {lr.coef:.3f} (95% CI {np.log(lr.IRR_lo):.3f} to {np.log(lr.IRR_hi):.3f}); the offset model assumes exactly 1.0.')
S(f'MLE alpha = {nb.params["alpha"]:.3f}. AIC: covariate model {nb.aic:.1f} vs offset model {nb0.aic:.1f}.')
S(f'Spearman rho (prompted rooms vs errors per layout) = {rho.statistic:.3f}, p = {rho.pvalue:.4f}.')
cm=reg[['prompted_room_count','bedroom_count','floor_area']].corr(method='spearman').round(2)
S(f'Predictor collinearity (Spearman): rooms-bedrooms {cm.iloc[0,1]}, rooms-area {cm.iloc[0,2]}, bedrooms-area {cm.iloc[1,2]}.\n')
pd.DataFrame({'spec':['covariate','offset'],'AIC':[nb.aic,nb0.aic],'log_rooms_coef':[lr.coef,1.0]}).to_csv(OUT+'/rq3_offset_vs_covariate.csv',index=False)

# ---------- RQ4a: pairwise logistic (all pairs, controlling log rooms), BH ----------
rows=[]
for a,b in itertools.combinations(cats,2):
    r=dict(A=a,B=b,n_A=int(B[a].sum()),n_B=int(B[b].sum()),n_both=int((B[a]&B[b]).sum()),phi=np.corrcoef(B[a],B[b])[0,1])
    try:
        m=sm.Logit(B[a],sm.add_constant(pd.DataFrame({'b':B[b],'lr':reg.log_rooms}))).fit(disp=0,maxiter=200)
        r.update(OR=np.exp(m.params['b']),OR_lo=np.exp(m.conf_int().loc['b',0]),OR_hi=np.exp(m.conf_int().loc['b',1]),p=m.pvalues['b'])
    except Exception: r.update(OR=np.nan,p=np.nan)
    rows.append(r)
P=pd.DataFrame(rows); P['q_BH']=multipletests(P.p,method='fdr_bh')[1]; P['p_bonferroni']=np.minimum(P.p*len(P),1)
P=P.sort_values('p'); P.to_csv(OUT+'/rq4_pairwise_logistic_all78.csv',index=False)
S(f'## RQ4a. Pairwise co-occurrence: logistic regression of A on B + log(rooms), all {len(P)} pairs, BH-FDR\n')
S(f'Nominal p < 0.05: {(P.p<.05).sum()} pairs (about {0.05*len(P):.1f} expected by chance). BH q < 0.05: {(P.q_BH<.05).sum()}. Bonferroni < 0.05: {(P.p_bonferroni<.05).sum()}. Smallest q = {P.q_BH.min():.3f}.\n')
S('Top 8 by p:\n\n| A | B | both/nA/nB | phi | OR [95% CI] | p | q (BH) |\n|---|---|---|---|---|---|---|')
for _,r in P.head(8).iterrows(): S(f'| {r.A} | {r.B} | {r.n_both}/{r.n_A}/{r.n_B} | {r.phi:.3f} | {r.OR:.2f} [{r.OR_lo:.2f}, {r.OR_hi:.2f}] | {r.p:.4f} | {r.q_BH:.3f} |')
S('')

# ---------- RQ4b/c: network (phi > 0.05) + centrality + bootstrap ----------
THR=0.05
def graph(Bm):
    C=np.nan_to_num(np.corrcoef(Bm.T)); G=nx.Graph(); G.add_nodes_from(range(len(C)))
    for i,j in itertools.combinations(range(len(C)),2):
        if C[i,j]>THR: G.add_edge(i,j,weight=C[i,j],dist=1/C[i,j])
    return G
def cent(G):
    n=G.number_of_nodes(); W=nx.to_numpy_array(G,nodelist=range(n),weight='weight'); v=np.abs(np.linalg.eigh(W)[1][:,-1])
    bw=nx.betweenness_centrality(G,weight='dist'); bu=nx.betweenness_centrality(G); dg=nx.degree_centrality(G)
    return v/np.linalg.norm(v) if False else v, np.array([bw[i] for i in range(n)]),np.array([bu[i] for i in range(n)]),np.array([dg[i] for i in range(n)])
G=graph(B.values); ev,bw,bu,dg=cent(G)
# eigenvector scaled like networkx (unit L2 norm) already
C0=pd.DataFrame({'category':cats,'domain':[DOM[c] for c in cats],'degree':dg,'betweenness_dist_1_over_phi':bw,'betweenness_unweighted':bu,'eigenvector_supplementary':ev}).sort_values('betweenness_dist_1_over_phi',ascending=False)
# alternative betweenness definitions (sensitivity only; the paper uses distance = 1/phi)
bp=nx.betweenness_centrality(G,weight='weight'); C0['betweenness_phi_as_length_(paper)']=[bp[cats.index(c)] for c in C0.category]
nT=1000*2; Bm=B.values; rng=np.random.default_rng(SEED); top=np.zeros((3,len(cats))); rk=[]; used=0
for _ in range(2000):
    bm=Bm[rng.integers(0,N,N)]
    if (bm.std(0)==0).any(): continue
    e,b_,_,d=cent(graph(bm)); used+=1
    top[0,np.argmax(e)]+=1; top[1,np.argmax(b_)]+=1
    mx=d.max(); top[2,np.where(d==mx)[0]]+=1.0/ (d==mx).sum()
    rk.append(pd.Series(e).rank(ascending=False).values)
top/=used; rk=np.array(rk)
C0['P_top_betweenness']=[top[1,cats.index(c)] for c in C0.category]; C0['P_top_degree_ties_split']=[top[2,cats.index(c)] for c in C0.category]; C0['P_top_eigenvector']=[top[0,cats.index(c)] for c in C0.category]
C0['eigen_rank_lo']=[np.percentile(rk[:,cats.index(c)],2.5) for c in C0.category]; C0['eigen_rank_hi']=[np.percentile(rk[:,cats.index(c)],97.5) for c in C0.category]
C0.to_csv(OUT+'/rq4_network_centrality_bootstrap.csv',index=False)
S(f'## RQ4b. Co-occurrence network (phi > {THR}; {G.number_of_edges()} edges) and bootstrap stability ({used} resamples of layouts)\n')
S('| Category | Degree | Betweenness (dist = 1/phi) | P(top betweenness) | P(top degree) | Eigenvector (suppl.) | P(top eigenvector) |\n|---|---|---|---|---|---|---|')
for _,r in C0.head(6).iterrows(): S(f'| {r.category} | {r.degree:.3f} | {r.betweenness_dist_1_over_phi:.3f} | {r.P_top_betweenness:.2f} | {r.P_top_degree_ties_split:.2f} | {r.eigenvector_supplementary:.3f} | {r.P_top_eigenvector:.2f} |')
S('\nFull table (13 categories, all variants): rq4_network_centrality_bootstrap.csv')
sym=C0[C0.category=='Symbol Errors'].iloc[0]
S(f'\nBetweenness of Symbol Errors under alternative definitions: distance = 1/phi (used in the paper) {sym.betweenness_dist_1_over_phi:.3f}; unweighted {sym.betweenness_unweighted:.3f}; phi as length {sym["betweenness_phi_as_length_(paper)"]:.3f}.\n')

# ---------- heatmap ----------
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
order=C0.sort_values(['domain','category']).category.tolist()  # group by domain
n=len(order); M=np.full((n,n),np.nan); Q=np.zeros((n,n),bool)
for _,r in P.iterrows():
    i,j=order.index(r.A),order.index(r.B); M[i,j]=M[j,i]=r.phi; Q[i,j]=Q[j,i]=r.p<.05
fig,ax=plt.subplots(figsize=(9,7.5)); im=ax.imshow(M,cmap='RdBu_r',vmin=-.25,vmax=.25)
for i in range(n):
    for j in range(n):
        if Q[i,j]: ax.text(j,i,'*',ha='center',va='center',fontsize=14,fontweight='bold')
ax.set_xticks(range(n)); ax.set_yticks(range(n)); ax.set_xticklabels(order,rotation=60,ha='right',fontsize=8); ax.set_yticklabels(order,fontsize=8)
cb=plt.colorbar(im,ax=ax,fraction=.04); cb.set_label('phi (layout-level co-occurrence)')
ax.set_title('Co-occurrence of 13 failure categories (300 layouts)\n* = nominal p < 0.05 adjusted for room count; none survive BH-FDR',fontsize=10)
plt.tight_layout(); plt.savefig(OUT+'/rq4_cooccurrence_heatmap.png',dpi=200)
open(OUT+'/results_summary.md','w').write('\n'.join(summ))
print('\n'.join(summ))
