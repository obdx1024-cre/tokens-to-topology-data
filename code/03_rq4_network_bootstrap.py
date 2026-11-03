"""RQ4 Figure 5: co-occurrence network (phi > 0.05) + betweenness with bootstrap 95% intervals. Run: python code/03_rq4_network_bootstrap.py"""
import sys, itertools, warnings, numpy as np, pandas as pd, networkx as nx
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent)); from common import load, load_annotators, RESULTS, FIGURES
OUT=str(RESULTS); lab, reg = load()
B=(pd.crosstab(lab.sample_id,lab.category)>0).astype(int).reindex(reg.sample_id).astype(int); cats=list(B.columns); Bm=B.values; N=len(Bm)
DOM={'Redundant / Unnecessary Elements':'Graphic Drafting','Symbol Errors':'Graphic Drafting','Line & Layer Errors':'Graphic Drafting','Scale Errors':'Graphic Drafting',
 'Missing Openings & Elements':'Spatial Circulation','Misaligned Openings & Elements':'Spatial Circulation','Blocked Pathways & Movement':'Spatial Circulation',
 'Poor Room Layout & Placement':'Ergonomic Layout','Dimension & Measurement Errors':'Ergonomic Layout','Space & Clearance Errors':'Ergonomic Layout','Missing Required Spaces':'Ergonomic Layout',
 'Overlapping / Colliding Elements':'Structural Logic','Unclear Wall Types':'Regulatory Compliance'}
SHORT={'Redundant / Unnecessary Elements':'Redundant Elements','Symbol Errors':'Symbol Errors','Line & Layer Errors':'Line & Layer','Scale Errors':'Scale Errors','Missing Openings & Elements':'Missing Openings',
 'Misaligned Openings & Elements':'Misaligned Openings','Blocked Pathways & Movement':'Blocked Pathways','Poor Room Layout & Placement':'Poor Room Layout','Dimension & Measurement Errors':'Dimension Errors',
 'Space & Clearance Errors':'Space & Clearance','Missing Required Spaces':'Missing Spaces','Overlapping / Colliding Elements':'Overlapping Elements','Unclear Wall Types':'Unclear Wall Types'}
COL={'Graphic Drafting':'#2f5f98','Spatial Circulation':'#d9822b','Ergonomic Layout':'#3f9a72','Structural Logic':'#8a5aa8','Regulatory Compliance':'#7a7a7a'}
def graph(M,thr=0.05):
    C=np.nan_to_num(np.corrcoef(M.T)); G=nx.Graph(); G.add_nodes_from(range(len(C)))
    for i,j in itertools.combinations(range(len(C)),2):
        if C[i,j]>thr: G.add_edge(i,j,weight=C[i,j],dist=1/C[i,j])
    return G
def met(G):
    n=G.number_of_nodes(); bw=nx.betweenness_centrality(G,weight='dist'); dg=nx.degree_centrality(G); return np.array([bw[i] for i in range(n)]),np.array([dg[i] for i in range(n)])
G=graph(Bm); bw,dg=met(G); rng=np.random.default_rng(2026); BW=[];DG=[]
for _ in range(2000):
    m=Bm[rng.integers(0,N,N)]
    if (m.std(0)==0).any(): continue
    b,d=met(graph(m)); BW.append(b); DG.append(d)
BW=np.array(BW); DG=np.array(DG)
tab=pd.DataFrame({'category':cats,'domain':[DOM[c] for c in cats],'degree':dg,'betweenness':bw,
  'betweenness_lo':np.percentile(BW,2.5,0),'betweenness_hi':np.percentile(BW,97.5,0),'degree_lo':np.percentile(DG,2.5,0),'degree_hi':np.percentile(DG,97.5,0)}).sort_values('betweenness',ascending=False)
tab.to_csv(OUT+'/figure5_centrality_bootstrap_intervals.csv',index=False)
s=cats.index('Symbol Errors'); nb=[cats[j] for j in G.neighbors(s)]
cross=sum(DOM[cats[u]]!=DOM[cats[v]] for u,v in G.edges())
info=f'Symbol Errors neighbours ({len(nb)}): '+'; '.join(f'{SHORT[c]} [{DOM[c]}]' for c in nb)+f'\nCross-domain edges: {cross}/{G.number_of_edges()}; Symbol Errors cross-domain edges: {sum(DOM[c]!="Graphic Drafting" for c in nb)}/{len(nb)}\n'
open(OUT+'/figure5_network_notes.txt','w').write(info); print(info); print(tab.round(3).to_string())
# ---- figure ----
order=sorted(range(len(cats)),key=lambda i:(list(COL).index(DOM[cats[i]]),cats[i]))
ang={i:np.pi/2-2*np.pi*k/len(cats) for k,i in enumerate(order)}; pos={i:np.array([np.cos(a),np.sin(a)]) for i,a in ang.items()}
fig=plt.figure(figsize=(7.6,4.1)); axA=fig.add_axes([0.0,0.10,0.60,0.86]); axB=fig.add_axes([0.745,0.22,0.25,0.68])
for u,v,dd in G.edges(data=True): axA.plot(*zip(pos[u],pos[v]),color='#9aa5b1' if s not in (u,v) else '#222',lw=0.6+dd['weight']*9,alpha=.55 if s not in (u,v) else .8,zorder=1,solid_capstyle='round')
for i in range(len(cats)):
    axA.scatter(*pos[i],s=60+dg[i]*520,color=COL[DOM[cats[i]]],edgecolor='#222' if i==s else 'white',linewidth=2.0 if i==s else 1.2,zorder=3)
    a=ang[i]; ha='left' if np.cos(a)>0.15 else ('right' if np.cos(a)<-0.15 else 'center'); va='bottom' if np.sin(a)>0.9 else ('top' if np.sin(a)<-0.9 else 'center')
    axA.text(1.22*np.cos(a),1.22*np.sin(a),SHORT[cats[i]],ha=ha,va=va,fontsize=7.6,fontweight='bold' if i==s else 'normal')
axA.set_xlim(-2.45,2.15); axA.set_ylim(-1.45,1.45); axA.set_aspect('equal'); axA.axis('off')
for d,c in COL.items(): axA.scatter([],[],s=40,color=c,label=d)
fig.legend(*axA.get_legend_handles_labels(),fontsize=7.4,frameon=False,loc='lower left',bbox_to_anchor=(0.02,0.0),ncol=5,handletextpad=.3,columnspacing=1.2)
t=tab.iloc[::-1].reset_index(drop=True); y=np.arange(len(t))
axB.hlines(y,t.betweenness_lo,t.betweenness_hi,color='#555',lw=1.1,zorder=2); axB.scatter(t.betweenness,y,color=[COL[d] for d in t.domain],s=34,zorder=3,edgecolor='white',linewidth=.6)
axB.set_yticks(y); axB.set_yticklabels([SHORT[c] for c in t.category],fontsize=7.4); axB.set_xlabel('Betweenness centrality\n(95% bootstrap interval)',fontsize=8); axB.tick_params(axis='x',labelsize=7.5); axB.tick_params(axis='y',length=0)
axB.grid(axis='x',color='#e3e3e3',lw=.6); axB.set_axisbelow(True); axB.set_xlim(left=-0.02)
for sp in ('top','right','left'): axB.spines[sp].set_visible(False)
fig.text(0.02,0.95,'A',fontsize=10,fontweight='bold'); fig.text(0.63,0.95,'B',fontsize=10,fontweight='bold')
plt.savefig(OUT+'/figure5_network.png',dpi=300); plt.savefig(OUT+'/figure5_network.pdf')
