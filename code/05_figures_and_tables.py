"""One consistent style for Figures 3-5 and Tables 1-2.
Every image is drawn at the same physical width (6.7 in) and exported at the same pixel width (4,468 px, ~667 dpi),
with the same font (Liberation Sans = Arial metrics) and the same point sizes, so text appears the same size on the page
when each image is inserted at full text width.
Run after 01 and 03: python code/05_figures_and_tables.py  (writes PNG + PDF to figures/)"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent)); from common import load, RESULTS, FIGURES
import sys, itertools, textwrap, numpy as np, pandas as pd, networkx as nx
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT=str(FIGURES); lab_all, reg_all = load()
W=6.7; WPX=4468; DPI=WPX/W
FS=8.0; FS_S=7.5; FS_H=8.5; FS_P=10.0                     # base, small, header, panel letter (pt)
plt.rcParams.update({'font.family':'Liberation Sans','font.size':FS,'axes.labelsize':FS,'xtick.labelsize':FS_S,'ytick.labelsize':FS,
                     'legend.fontsize':FS_S,'axes.linewidth':0.8,'xtick.major.width':0.8,'ytick.major.width':0.8,'pdf.fonttype':42})
def save(fig,name):
    fig.savefig(f'{OUT}/{name}.png',dpi=DPI,facecolor='white'); fig.savefig(f'{OUT}/{name}.pdf',facecolor='white'); plt.close(fig)
    from PIL import Image; print(name, Image.open(f'{OUT}/{name}.png').size)
BLUE,ORANGE='#2f5f98','#d9822b'

# ---------------- Figure 3 ----------------
d=pd.read_csv(RESULTS/'rq1_shares_prevalence.csv'); d=d[d.level=='domain'].sort_values('prevalence_pct')
fig,ax=plt.subplots(figsize=(W,2.15)); fig.subplots_adjust(left=0.2,right=0.97,top=0.96,bottom=0.22)
y=np.arange(len(d)); p=d.prevalence_pct.values
ax.barh(y,p,height=0.52,color=BLUE,zorder=2)
ax.errorbar(p,y,xerr=[p-d.wilson_lo_pct.values,d.wilson_hi_pct.values-p],fmt='none',ecolor='#222',elinewidth=0.9,capsize=3,capthick=0.9,zorder=3)
for yi,(hi,pv,n) in enumerate(zip(d.wilson_hi_pct,p,d.layouts_affected)): ax.text(hi+1.2,yi,f'{pv:.1f}%  (n = {n})',va='center',fontsize=FS_S,color='#222')
ax.set_yticks(y); ax.set_yticklabels(d.name); ax.tick_params(axis='y',length=0)
ax.set_xlim(0,118); ax.set_xticks(range(0,101,20)); ax.set_xticklabels([f'{v}%' for v in range(0,101,20)])
ax.xaxis.grid(True,color='#ddd',lw=.6,zorder=0); ax.set_axisbelow(True)
for s in ('top','right','left'): ax.spines[s].set_visible(False)
ax.set_xlabel('Layouts containing at least one error (%), N = 300; whiskers show 95% Wilson CI',loc='right')
save(fig,'figure3_domain_prevalence')

# ---------------- Figure 4 ----------------
sev=pd.crosstab(lab_all.model,lab_all.severity,normalize='index')*100
op=pd.crosstab(lab_all.sample_id,lab_all.category).reindex(reg_all.sample_id).fillna(0).groupby(reg_all.set_index('sample_id').model).mean()
models=['GPT-Image-1','GPT-Image-2','Nano Banana Pro','Nano Banana 2']; y=np.arange(4)[::-1]
fig=plt.figure(figsize=(W,2.55)); axA=fig.add_axes([0.165,0.18,0.405,0.62]); axB=fig.add_axes([0.62,0.18,0.355,0.62])
cols={'Critical':'#1f3f6a','Functional':'#5f8fc4','Cosmetic':'#c6d9ef'}; left=np.zeros(4)
for k,c in cols.items():
    v=sev.loc[models,k].values; axA.barh(y,v,left=left,color=c,height=0.6,edgecolor='white',lw=1,label=k)
    for yi,l,vv in zip(y,left,v): axA.text(l+vv/2,yi,f'{vv:.0f}%',ha='center',va='center',fontsize=FS_S,color='white' if k!='Cosmetic' else '#222')
    left+=v
axA.set_yticks(y); axA.set_yticklabels(models); axA.tick_params(axis='y',length=0); axA.set_xlim(0,100); axA.set_xticks([0,25,50,75,100])
axA.set_xlabel("Share of each model's errors (%)")
axA.legend(ncol=3,frameon=False,loc='lower left',bbox_to_anchor=(-0.01,1.0),handlelength=1,handleheight=.7,columnspacing=1.2,handletextpad=0.4)
h=0.3; a=op.loc[models,'Missing Openings & Elements'].values; b=op.loc[models,'Misaligned Openings & Elements'].values
axB.barh(y+h/2,a,height=h,color=BLUE,label='Missing Openings'); axB.barh(y-h/2,b,height=h,color=ORANGE,label='Misaligned Openings')
for yy,vv in list(zip(y+h/2,a))+list(zip(y-h/2,b)): axB.text(vv+0.05,yy,f'{vv:.2f}',va='center',fontsize=FS_S)
axB.set_yticks([]); axB.set_xlim(0,3.2); axB.set_xticks([0,1,2,3]); axB.set_xlabel('Mean errors per layout')
axB.xaxis.grid(True,color='#ddd',lw=.6); axB.set_axisbelow(True)
axB.legend(ncol=2,frameon=False,loc='lower left',bbox_to_anchor=(-0.02,1.0),handlelength=1,handleheight=.7,columnspacing=1.0,handletextpad=0.4)
for ax in (axA,axB):
    for s in ('top','right'): ax.spines[s].set_visible(False)
fig.text(0.012,0.93,'A',fontsize=FS_P,fontweight='bold'); fig.text(0.6,0.93,'B',fontsize=FS_P,fontweight='bold')
save(fig,'figure4_model_comparison')

# ---------------- Figure 5 (shorter) ----------------
lab, reg = lab_all, reg_all
Bm=(pd.crosstab(lab.sample_id,lab.category)>0).astype(int).reindex(reg.sample_id).astype(int); cats=list(Bm.columns); Bm=Bm.values
DOM={'Redundant / Unnecessary Elements':'Graphic Drafting','Symbol Errors':'Graphic Drafting','Line & Layer Errors':'Graphic Drafting','Scale Errors':'Graphic Drafting',
 'Missing Openings & Elements':'Spatial Circulation','Misaligned Openings & Elements':'Spatial Circulation','Blocked Pathways & Movement':'Spatial Circulation',
 'Poor Room Layout & Placement':'Ergonomic Layout','Dimension & Measurement Errors':'Ergonomic Layout','Space & Clearance Errors':'Ergonomic Layout','Missing Required Spaces':'Ergonomic Layout',
 'Overlapping / Colliding Elements':'Structural Logic','Unclear Wall Types':'Regulatory Compliance'}
SHORT={'Redundant / Unnecessary Elements':'Redundant Elements','Symbol Errors':'Symbol Errors','Line & Layer Errors':'Line & Layer','Scale Errors':'Scale Errors','Missing Openings & Elements':'Missing Openings',
 'Misaligned Openings & Elements':'Misaligned Openings','Blocked Pathways & Movement':'Blocked Pathways','Poor Room Layout & Placement':'Poor Room Layout','Dimension & Measurement Errors':'Dimension Errors',
 'Space & Clearance Errors':'Space & Clearance','Missing Required Spaces':'Missing Spaces','Overlapping / Colliding Elements':'Overlapping Elements','Unclear Wall Types':'Unclear Wall Types'}
COL={'Graphic Drafting':BLUE,'Spatial Circulation':ORANGE,'Ergonomic Layout':'#3f9a72','Structural Logic':'#8a5aa8','Regulatory Compliance':'#7a7a7a'}
C=np.nan_to_num(np.corrcoef(Bm.T)); G=nx.Graph(); G.add_nodes_from(range(len(cats)))
for i,j in itertools.combinations(range(len(cats)),2):
    if C[i,j]>0.05: G.add_edge(i,j,weight=C[i,j])
dg=nx.degree_centrality(G); s=cats.index('Symbol Errors')
tab=pd.read_csv(RESULTS/'figure5_centrality_bootstrap_intervals.csv').sort_values('betweenness',ascending=False)
assert abs(dg[s]-0.75)<1e-9 and G.number_of_edges()==34
order=sorted(range(len(cats)),key=lambda i:(list(COL).index(DOM[cats[i]]),cats[i]))
ang={i:np.pi/2-2*np.pi*k/len(cats) for k,i in enumerate(order)}; pos={i:np.array([np.cos(a),np.sin(a)]) for i,a in ang.items()}
H5=2.75; fig=plt.figure(figsize=(W,H5))
axA=fig.add_axes([0.0,0.02,0.6,0.86]); axB=fig.add_axes([0.765,0.2,0.225,0.68])
for u,v,dd in G.edges(data=True):
    hub=s in (u,v); axA.plot(*zip(pos[u],pos[v]),color='#222' if hub else '#9aa5b1',lw=0.4+dd['weight']*6,alpha=.8 if hub else .5,zorder=1,solid_capstyle='round')
for i in range(len(cats)):
    axA.scatter(*pos[i],s=30+dg[i]*260,color=COL[DOM[cats[i]]],edgecolor='#222' if i==s else 'white',linewidth=1.5 if i==s else 0.8,zorder=3)
    a=ang[i]; ha='left' if np.cos(a)>0.15 else ('right' if np.cos(a)<-0.15 else 'center'); va='bottom' if np.sin(a)>0.9 else ('top' if np.sin(a)<-0.9 else 'center')
    axA.text(1.17*np.cos(a),1.17*np.sin(a),SHORT[cats[i]],ha=ha,va=va,fontsize=FS_S,fontweight='bold' if i==s else 'normal')
axA.set_xlim(-2.3,2.3); axA.set_ylim(-1.32,1.32); axA.set_aspect('equal'); axA.axis('off')
for dname,c in COL.items(): axA.scatter([],[],s=28,color=c,label=dname)
fig.legend(*axA.get_legend_handles_labels(),frameon=False,loc='upper left',bbox_to_anchor=(0.04,1.0),ncol=5,handletextpad=.2,columnspacing=1.4)
t=tab.iloc[::-1].reset_index(drop=True); yy=np.arange(len(t))
axB.hlines(yy,t.betweenness_lo,t.betweenness_hi,color='#555',lw=0.9,zorder=2)
axB.scatter(t.betweenness,yy,color=[COL[x] for x in t.domain],s=18,zorder=3,edgecolor='white',linewidth=.5)
axB.set_yticks(yy); axB.set_yticklabels([SHORT[c] for c in t.category],fontsize=FS_S); axB.tick_params(axis='y',length=0); axB.set_ylim(-0.6,len(t)-0.4)
axB.set_xlabel('Betweenness (95% CI)'); axB.set_xticks([0,0.2,0.4]); axB.set_xlim(-0.02,0.45)
axB.grid(axis='x',color='#e3e3e3',lw=.5); axB.set_axisbelow(True)
for sp in ('top','right','left'): axB.spines[sp].set_visible(False)
fig.text(0.012,0.86,'A',fontsize=FS_P,fontweight='bold'); fig.text(0.6,0.86,'B',fontsize=FS_P,fontweight='bold')
save(fig,'figure5_network')

# ---------------- Tables (column layout, same frame as the user's design) ----------------
def table(name,corner,side,groups,leftw=0.95,top_pad=0.06,line=0.135):
    fig=plt.figure(figsize=(W,3)); r=fig.canvas.get_renderer()
    colw=(W-leftw)/len(groups)
    def width(txt,fs,bold):
        t=fig.text(0,0,txt,fontsize=fs,fontweight='bold' if bold else 'normal'); bb=t.get_window_extent(renderer=r); t.remove(); return bb.width/fig.dpi
    def wrap(txt,fs,bold,maxw):
        words=txt.split(); lines=['']
        for w_ in words:
            trial=(lines[-1]+' '+w_).strip()
            if width(trial,fs,bold)<=maxw or not lines[-1]: lines[-1]=trial
            else: lines.append(w_)
        return lines
    nrows=max(len(g[2]) for g in groups)
    cell=[[ (wrap(it[0],FS,False,colw-0.1),it[1]) if k<len(g[2]) else None for g in groups for it in [g[2][k] if k<len(g[2]) else ('','')]] for k in range(nrows)]
    rowlines=[max(len(c[0])+1 if c else 1 for c in row) for row in cell]
    titles=[wrap(g[0],FS_H,True,colw-0.08) for g in groups]; hl=max(len(t) for t in titles)
    hdr=(hl+1)*line+0.16; body=sum(rowlines)*line+0.11*nrows+0.10
    H=top_pad*2+hdr+body; fig.set_size_inches(W,H)
    ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,W); ax.set_ylim(H,0); ax.axis('off')
    T=top_pad; Hd=T+hdr; Bt=H-top_pad
    for yv,lw in ((T,1.3),(Hd,0.7),(Bt,1.3)): ax.plot([0,W],[yv,yv],color='k',lw=lw,solid_capstyle='butt')
    ax.plot([leftw,leftw],[T,Bt],color='#bbb',lw=0.5)
    hy=(T+Hd)/2
    ax.text(leftw/2,hy-line*0.45,corner[0],ha='center',va='center',fontsize=FS_H,fontweight='bold'); ax.text(leftw/2,hy+line*0.55,corner[1],ha='center',va='center',fontsize=FS_S)
    by=(Hd+Bt)/2
    ax.text(leftw/2,by-line*0.45,side[0],ha='center',va='center',fontsize=FS_H,fontweight='bold'); ax.text(leftw/2,by+line*0.55,side[1],ha='center',va='center',fontsize=FS_S)
    for gi,(title,sub,items) in enumerate(groups):
        cx=leftw+colw*(gi+0.5)
        tl=titles[gi]; ty=hy-(hl+1)*line/2+line/2+(hl-len(tl))*line/2
        for L in tl: ax.text(cx,ty,L,ha='center',va='center',fontsize=FS_H,fontweight='bold'); ty+=line
        ax.text(cx,ty,sub,ha='center',va='center',fontsize=FS_S)
    y0=Hd+0.05
    for k,row in enumerate(cell):
        hgt=rowlines[k]*line+0.11
        for gi,c in enumerate(row):
            cx=leftw+colw*(gi+0.5); cy=y0+hgt/2
            if c is None: ax.text(cx,cy,'–',ha='center',va='center',fontsize=FS,color='#666'); continue
            nm,val=c; blk=(len(nm)+1)*line; ty=cy-blk/2+line/2
            for L in nm: ax.text(cx,ty,L,ha='center',va='center',fontsize=FS); ty+=line
            ax.text(cx,ty,val,ha='center',va='center',fontsize=FS_S,color='#444')
        y0+=hgt
    save(fig,name)
T1=[('Graphic Drafting','Σ = 1,080 (39.9%)',[('Redundant / Unnecessary Elements','n = 443 (16.4%)'),('Symbol Errors','n = 289 (10.7%)'),('Line & Layer Errors','n = 263 (9.7%)'),('Scale Errors','n = 85 (3.1%)')]),
    ('Spatial Circulation','Σ = 879 (32.5%)',[('Missing Openings & Elements','n = 519 (19.2%)'),('Misaligned Openings & Elements','n = 208 (7.7%)'),('Blocked Pathways & Movement','n = 152 (5.6%)')]),
    ('Ergonomic Layout','Σ = 672 (24.8%)',[('Poor Room Layout & Placement','n = 344 (12.7%)'),('Dimension & Measurement Errors','n = 147 (5.4%)'),('Space & Clearance Errors','n = 141 (5.2%)'),('Missing Required Spaces','n = 40 (1.5%)')]),
    ('Structural Logic','Σ = 48 (1.8%)',[('Overlapping / Colliding Elements','n = 48 (1.8%)')]),
    ('Regulatory Compliance','Σ = 29 (1.1%)',[('Unclear Wall Types','n = 29 (1.1%)')])]
def columns_fitted(name,groups,line=0.118,gap=0.085,hgap=0.11,pad=0.045,margin=0.05):
    fig=plt.figure(figsize=(W,3)); r=fig.canvas.get_renderer()
    def wd(t,fs,b=False):
        x=fig.text(0,0,t,fontsize=fs,fontweight='bold' if b else 'normal'); w=x.get_window_extent(r).width/fig.dpi; x.remove(); return w
    def best2(t):                      # best split of a name into one or two balanced lines
        ws=t.split(); opts=[[t]]+[[' '.join(ws[:k]),' '.join(ws[k:])] for k in range(1,len(ws))]
        return min(opts,key=lambda o:(max(wd(l,FS) for l in o),len(o)))
    labels=[('Domains','macro-level'),('Categories','13 types')]
    leftw=max(max(wd(a,FS,True),wd(b,FS_S)) for a,b in labels)+2*margin
    need=[]; split={}
    nr=max(len(g[2]) for g in groups)
    def row_single(k,skip):
        return all(wd(g[2][k][0],FS)<=1.0 for gi,g in enumerate(groups) if gi!=skip and k<len(g[2]))
    for gi,(title,sub,items) in enumerate(groups):
        w=max(wd(title,FS,True),wd(sub,FS_S))
        for nm,v in items:
            one=wd(nm,FS); two=best2(nm); k=[x[0] for x in items].index(nm)
            if row_single(k,gi): two=[nm]
            split[nm]=two
            w=max(w,wd(v,FS_S),max(wd(l,FS) for l in two))
        need.append(w+margin)
    free=W-leftw-sum(need); assert free>=0, f'does not fit: {free:.3f} in short'
    # use the spare width to keep short names on one line where possible, then spread the rest evenly
    colw=[n+free/len(need) for n in need]
    lines={}
    for (title,sub,items),cw in zip(groups,colw):
        for nm,v in items: lines[nm]=[nm] if wd(nm,FS)<=cw-margin else split[nm]
    nrows=max(len(g[2]) for g in groups)
    cells=[[(lines[g[2][k][0]],g[2][k][1]) if k<len(g[2]) else None for g in groups] for k in range(nrows)]
    rl=[max(len(c[0])+1 if c else 1 for c in row) for row in cells]
    hdr=2*line+2*hgap; body=sum(rl)*line+gap*(nrows+1); H=2*pad+hdr+body; fig.set_size_inches(W,H)
    ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,W); ax.set_ylim(H,0); ax.axis('off')
    T=pad; Hd=T+hdr; B=H-pad
    for yv,lw in ((T,1.3),(Hd,0.7),(B,1.3)): ax.plot([0,W],[yv,yv],color='k',lw=lw,solid_capstyle='butt')
    ax.plot([leftw,leftw],[T,B],color='#bbb',lw=0.5)
    for (a,b),yc in zip(labels,((T+Hd)/2,(Hd+B)/2)):
        ax.text(leftw/2,yc-line/2,a,ha='center',va='center',fontsize=FS,fontweight='bold'); ax.text(leftw/2,yc+line/2,b,ha='center',va='center',fontsize=FS_S)
    xs=np.cumsum([leftw]+colw[:-1]); cx=[x+w/2 for x,w in zip(xs,colw)]
    for (title,sub,_),c in zip(groups,cx):
        yc=(T+Hd)/2; ax.text(c,yc-line/2,title,ha='center',va='center',fontsize=FS,fontweight='bold'); ax.text(c,yc+line/2,sub,ha='center',va='center',fontsize=FS_S)
    y=Hd+gap
    for k,row in enumerate(cells):
        h=rl[k]*line
        for c,cell in zip(cx,row):
            if cell is None: ax.text(c,y+h/2,'–',ha='center',va='center',fontsize=FS,color='#666'); continue
            nm,v=cell; ty=y+(h-(len(nm)+1)*line)/2+line/2
            for L in nm: ax.text(c,ty,L,ha='center',va='center',fontsize=FS); ty+=line
            ax.text(c,ty,v,ha='center',va='center',fontsize=FS_S,color='#444')
        y+=h+gap
    print('column widths (in):',[round(leftw,2)]+[round(w,2) for w in colw])
    save(fig,name)
columns_fitted('table1_taxonomy',T1)
f=lambda irr,lo,hi,p: f'{irr} [{lo}, {hi}], p = {p}'
T2=[('Morphology','ref. Rectangular',[('Asymmetrical/Irregular',f('1.115','0.942','1.319','0.205')),('Composite/Sprawling',f('1.258','0.932','1.698','0.134')),('L-Shape',f('0.899','0.738','1.095','0.290'))]),
    ('Model Cohort','ref. Nano Banana Pro',[('GPT-Image-1',f('1.152','0.952','1.394','0.145')),('GPT-Image-2',f('0.983','0.806','1.199','0.869')),('Nano Banana 2',f('0.994','0.808','1.224','0.958'))]),
    ('Layout Size','per-unit effects',[('Bedroom Count',f('1.031','0.930','1.144','0.562')),('Floor Area (m²)',f('0.999','0.997','1.001','0.349')),('Prompted Room Count (log)',f('1.164','0.869','1.560','0.307'))])]
table('table2_negbin',('Predictors','group / ref. level'),('Estimates','IRR [95% CI], p'),T2)
