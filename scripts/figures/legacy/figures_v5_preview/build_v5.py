"""Four main-figure previews. Read-only scientific inputs; no Word or ED access.

AREA_WEIGHTS = {
 'Figure1': {'a': 0.12, 'b': 0.14, 'c': 0.14, 'd': 0.51, 'e': 0.09},
 # d is the decisive 16-context matched-null evidence; UMAP is context only.
 'Figure2': {'a': 0.18, 'b': 0.17, 'd': 0.34, 'c': 0.09, 'e': 0.09, 'g': 0.09, 'f': 0.04},
 # Nine biological identities receive >3x the plot area of either four-layer audit.
 'Figure3': {'a': 0.34, 'b': 0.10, 'c': 0.10, 'e': 0.26, 'f': 0.15, 'd': 0.05},
 # Real cell map and the existing 16335-gene transfer evaluation are two heroes.
 'Figure4': {'a': 0.19, 'b': 0.19, 'c': 0.30, 'd': 0.09, 'e': 0.14, 'f': 0.09},
 # The compact framework yields space to actual nested-target observations.
}
Weights describe intended reading emphasis, not a claim of exact ink coverage.
All reported statistics are read from the existing source tables, never refit.
"""
from pathlib import Path
import os,sys,json,hashlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch,Ellipse
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

ROOT=Path(__file__).resolve().parent
PROJECT=Path(__file__).resolve().parents[4];DATA=PROJECT/'data';QA=ROOT/'qa';QA.mkdir(exist_ok=True)
sys.path.insert(0,str(PROJECT/'scripts'/'utilities'))
from audit_panel_alignment import require_matplotlib_panel_alignment
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],
 'font.size':6.5,'axes.labelsize':6.7,'axes.titlesize':7.2,'xtick.labelsize':6,'ytick.labelsize':6,
 'text.color':'#22252A','axes.labelcolor':'#22252A','axes.edgecolor':'#3D4248',
 'xtick.color':'#3D4248','ytick.color':'#3D4248','axes.linewidth':.65,
 'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,
 'legend.frameon':False,'legend.fontsize':6,'lines.linewidth':.65,
 'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white'})
BLUE='#4B6F8A';VIOLET='#776A8A';OCHRE='#B99155';OLIVE='#738267';RUST='#B36F5D'
INK='#22252A';AXIS='#3D4248';MID='#A6ADB5';LIGHT='#D9DDE2';PALE='#EEF0F2'
COND=['Mild.Drought','Moderate.Drought','Mild.Salinity','Moderate.Salinity']
CTEXT=['Mild drought','Moderate drought','Mild salinity','Moderate salinity']
CSHORT=['Mild\ndrought','Moderate\ndrought','Mild\nsalinity','Moderate\nsalinity']
CC=[BLUE,VIOLET,OCHRE,OLIVE]
INPUTS={}

def read(rel):
 p=DATA/rel;INPUTS[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return pd.read_csv(p,sep='\t')

F=[read('paper_figures_v2/'+s) for s in ['fig1_candidate_stability_and_problem_definition.tsv','fig2_future_time.tsv','fig3_missing_state_and_cross_experiment.tsv','fig4_reliability_stratification_and_boundary.tsv']]
B=read('publication_baseline_suite_v2_corrected/baseline_suite_long.tsv')

def sel(d,**kw):
 for key,val in kw.items():d=d[d[key].isin(val) if isinstance(val,list) else d[key].eq(val)]
 return d.copy()

def wide(d,section=None,idx='context'):
 if section:d=sel(d,section=section)
 return d.pivot_table(index=idx,columns='metric',values='value',aggfunc='first',sort=False)

figs=[];AUDITS={}
def canvas(name,h):
 f=plt.figure(figsize=(180/25.4,h/25.4),facecolor='white');figs.append((name,f));return f

def panel(f,r,letter,title='',off=False,pad=11):
 a=f.add_axes(r);a.set_title(title,loc='left',fontsize=7.2,pad=pad,weight='normal')
 a.annotate(letter,(0,1),xycoords='axes fraction',xytext=(-17,pad),textcoords='offset points',fontsize=8.8,weight='bold',va='baseline')
 a.tick_params(length=2.5,width=.6)
 if off:a.set(xlim=(0,1),ylim=(0,1));a.axis('off')
 return a

def say(a,x,y,s,ha='center',color=INK,size=6.2):
 return a.text(x,y,s,ha=ha,va='center',color=color,fontsize=size)

def arr(a,st,en,color=AXIS,lw=.65):
 a.add_patch(FancyArrowPatch(st,en,arrowstyle='-|>',mutation_scale=6,lw=lw,color=color,shrinkA=3,shrinkB=3))

def zero(a):a.axhline(0,c=LIGHT,ls='--',lw=.6,zorder=0)

def paired(a,d,metric,models):
 p=sel(d,metric=metric).pivot_table(index='context',columns='model',values='value',sort=False)
 for _,r in p.iterrows():a.plot([0,1],[r[models[0]],r[models[1]]],c=MID,lw=.55)
 a.scatter(np.zeros(len(p)),p[models[0]],s=17,c=BLUE,zorder=3)
 a.scatter(np.ones(len(p)),p[models[1]],s=17,c=MID,zorder=3)
 a.set_xticks([0,1],['History','Shared'],rotation=45,ha='right',rotation_mode='anchor')
 a.set(xlim=(-.3,1.3),ylim=(-.1,1));zero(a)

def qstrip(a,p,col,labels):
 assert (p[col]>0).all()
 a.scatter(range(len(p)),-np.log10(p[col]),c=[BLUE if q<.05 else RUST for q in p[col]],s=15,zorder=3)
 a.axhline(-np.log10(.05),c=LIGHT,ls='--',lw=.6)
 a.set(ylim=(0,4.5),ylabel='−log10(q)',xlim=(-.5,len(p)-.5))
 a.set_xticks(range(len(p)),labels);a.set_yticks([0,4]);a.tick_params(length=2)

# FIGURE 1 — majority area assigned to the empirical null, not context.
d=F[0];f=canvas('Figure1_v5',174)
a=panel(f,[.065,.675,.25,.28],'a','Rice leaf nuclei')
rice=pd.read_csv(ROOT/'rice_deposited_umap.tsv.gz',sep='\t');assert len(rice)==126139
graycolors={'Epidermis':'#B9C2C9','Meristem':'#7D8D99','Mesophyll':'#D5DADF','Vasculature':'#9EA9B2'}
for tissue,gg in rice.groupby('TissueSystem',sort=False):
 a.scatter(gg.UMAP1,gg.UMAP2,s=.18,c=graycolors[tissue],lw=0,rasterized=True)
centers=rice.groupby('TissueSystem')[['UMAP1','UMAP2']].median()
for tissue,r in centers.iterrows():
 a.annotate(tissue,(r.UMAP1,r.UMAP2),xytext=(0,0),textcoords='offset points',ha='center',fontsize=5.8,color=AXIS)
a.set(xlabel='UMAP 1',ylabel='UMAP 2');a.set_xticks([]);a.set_yticks([])
b=panel(f,[.44,.675,.20,.28],'b','Replicate concordance')
c=panel(f,[.765,.675,.20,.28],'c','Top100 recurrence')
for aa,sec,metric in [(b,'rice_full_vector','replicate_response_spearman'),(c,'rice_top100_random_expectation','observed_overlap')]:
 for i,ct in enumerate(COND):
  gg=sel(d,section=sec,metric=metric);v=gg[gg.context.str.startswith(ct)].value.to_numpy();assert len(v)==4
  aa.scatter(i+np.linspace(-.13,.13,4),v,s=15,c=CC[i],zorder=3)
  aa.plot([i-.23,i+.23],[np.median(v)]*2,c=AXIS,lw=.9)
 aa.set_xticks(range(4),CSHORT,rotation=35,ha='right',rotation_mode='anchor');aa.set_xlim(-.5,3.5)
b.set(ylim=(-.05,.30),ylabel='Replicate response ρ');zero(b)
c.set(ylim=(0,36),ylabel='Top100 overlap');c.axhline(sel(d,metric='random_expected_overlap').value.iloc[0],c=LIGHT,ls='--',lw=.6)
dd=panel(f,[.20,.105,.55,.43],'d','Matched empirical null',pad=21)
p=wide(d,'rice_top100_matched_null').reindex([ct+'|'+t for ct in COND for t in ['Epidermis','Meristem','Mesophyll','Vasculature']])
yy=np.arange(16)+np.repeat(np.arange(4)*.85,4)
dd.hlines(yy,0,p.null_q99,color=LIGHT,lw=3.5)
dd.scatter(p.null_mean,yy,c=MID,s=10,zorder=2)
dd.scatter(p.observed_overlap,yy,c=np.repeat(CC,4),s=20,zorder=3)
dd.set_yticks(yy,[x.split('|')[1] for x in p.index]);dd.invert_yaxis();dd.set(xlim=(0,48),xlabel='Top100 overlap')
dd.text(0,1.022,'16/16 BH q ≈ 1e-4',transform=dd.transAxes,fontsize=6)
for i,txt in enumerate(CTEXT):dd.text(47,yy[i*4:i*4+4].mean(),txt,ha='right',va='center',fontsize=6,color=CC[i])
e=panel(f,[.875,.19,.09,.25],'e','Enrichment')
p=wide(d,'rice_condition_summary').reindex(COND)
e.scatter(np.arange(4),p.median_matched_enrichment,s=20,c=CC)
e.set(xlim=(-.5,3.5),ylim=(8,46),ylabel='Median enrichment (fold)')
e.set_xticks(range(4),['MD','MoD','MS','MoS'],rotation=90,ha='right',rotation_mode='anchor')
f.legend(handles=[Line2D([],[],c=LIGHT,lw=3.5,label='Null: 0–99th percentile'),Line2D([],[],c=MID,marker='o',ls='',ms=3,label='Null mean')],loc='lower left',bbox_to_anchor=(.19,.018),ncol=2,handlelength=1.3,fontsize=6)
AUDITS[f]=([b,c],[[b,c]],[], 'The biological-context panel and the wide null hero have intentional unequal spans.')

# FIGURE 2 — nine biological identities make a large, readable data hero.
d=F[1];f=canvas('Figure2_v5',188)
a=panel(f,[.09,.765,.465,.19],'a','Eligible evidence',True)
for j,txt in enumerate(['State A','State B','State C']):say(a,-.035,.76-j*.28,txt,ha='right',size=6)
for x,txt in [(.09,'0 h'),(.44,'3 h'),(.90,'8 h target')]:say(a,x,.99,txt,size=6)
for x,y,target in [(.09,.76,False),(.44,.76,False),(.90,.76,True),(.44,.48,False),(.44,.20,False)]:
 a.add_patch(Ellipse((x,y),.032,.062,fc='white' if target else BLUE,ec=AXIS if target else BLUE,lw=.7,ls='--' if target else '-'))
arr(a,(.12,.76),(.41,.76),LIGHT)
arr(a,(.48,.76),(.87,.76),BLUE,1.1)
# Thin lines indicate alternative evidence without a fan of arrowheads.
a.plot([.47,.87],[.48,.72],c=MID,lw=.45);a.plot([.47,.87],[.20,.72],c=MID,lw=.45)
say(a,.68,.89,'Same-state history',color=BLUE,size=6);say(a,.72,.18,'Cross-state evidence',color=AXIS,size=5.8)
b=panel(f,[.715,.745,.25,.21],'b','Salt temporal baselines')
for i,ds in enumerate(['GSE153103_salt_6h','GSE153103_salt_12h']):
 gg=sel(B,dataset=ds);ms=[next(m for m in gg.model if m.startswith(pre)) for pre in ['latest','history_mean','linear_last2']]
 ms.append(next(m for m in gg.model if m.startswith('linear_') and not m.startswith('linear_last2')))
 b.plot(range(4),gg.set_index('model').loc[ms].spearman_full_vector,'o-',c=[BLUE,VIOLET][i],ms=3,label=['6 h','12 h'][i])
b.set_xticks(range(4),['Latest','Mean','Last-two\ntrend','Baseline\ntrend'],rotation=35,ha='right',rotation_mode='anchor')
b.set(ylim=(-.1,1),ylabel='Full-response ρ');zero(b)
b.legend(ncol=2,loc='upper right',fontsize=6,handlelength=1.1,columnspacing=.9)
b.text(.02,.16,'Top100: 36 vs. 30 (6 h)',transform=b.transAxes,fontsize=5.8)
dd=panel(f,[.26,.29,.32,.335],'d','Lateral-root cell types',pad=21)
gg=sel(d,dataset='lateral_root_8h_to_20h',section='frozen_reveal_primary_results')
p=sel(gg,metric='full_rho').pivot_table(index='context',columns='model',values='value',sort=False)
names={'Ambiguous Stele Cells':'Ambiguous stele','Columella/ Root Cap':'Columella / root cap','Epidermis':'Epidermis','Lateral Root Endodermis':'Lateral-root endodermis','Mature Pericycle':'Mature pericycle','Phloem':'Phloem','Phloem Pole Pericycle':'Phloem-pole pericycle','Xylem':'Xylem','Xylem Pole Pericycle':'Xylem-pole pericycle'}
for i,(_,r) in enumerate(p.iterrows()):
 dd.plot([r.same_state_latest,r.shared_other_states],[i,i],c=MID,lw=.6)
 dd.scatter(r.same_state_latest,i,c=BLUE,s=20,zorder=3);dd.scatter(r.shared_other_states,i,c=MID,s=16,zorder=3)
dd.set_yticks(range(9),[names[x] for x in p.index]);dd.invert_yaxis();dd.set(xlim=(-.1,1),xlabel='Full-response ρ at 20 h')
dd.text(0,1.024,'9/9 wins; exact P = 0.00195',transform=dd.transAxes,fontsize=6)
dd.text(.01,-.18,'● Same-state 8 h',transform=dd.transAxes,fontsize=6,color=BLUE)
dd.text(.64,-.18,'● Shared',transform=dd.transAxes,fontsize=6,color=AXIS)
c=panel(f,[.705,.345,.105,.25],'c','1→3 h\nroot layers')
paired(c,sel(d,dataset='GSE46205_salt_3h_corrected',section='post_reveal_preprocessing_corrected_results'),'spearman_full_vector',['same_state_history','shared_other3']);c.set_ylabel('Full-response ρ')
e=panel(f,[.86,.345,.105,.25],'e','3→8 h\nroot layers')
paired(e,sel(d,dataset='GSE46205_salt_8h_corrected',section='post_reveal_preprocessing_corrected_results'),'spearman_full_vector',['same_state_history','shared_other3'])
f.text(.68,.265,'GSE46205 · corrected post-reveal',fontsize=5.8,color=AXIS)
g=panel(f,[.10,.08,.58,.055],'g','Matched-null significance',pad=18)
p=wide(d,'matched_empirical_null',idx=['dataset','context']);qstrip(g,p,'q_bh_all_future_time_tests_corrected',['S6','S12','Su']+['L'+str(i) for i in range(1,10)]+['R'+str(i) for i in range(1,5)])
g.text(0,1.05,'13/16 significant; lateral root 6/9',transform=g.transAxes,fontsize=5.8)
ff=panel(f,[.845,.08,.12,.055],'f','Sulfur depletion',pad=18)
p=wide(sel(d,dataset='sulfur_3h_to_12h'),'matched_empirical_null')
ff.hlines([0],0,p.null_q99,color=LIGHT,lw=3.5);ff.scatter(p.null_mean,[0],c=MID,s=10);ff.scatter(p.observed_overlap,[0],c=BLUE,s=20,zorder=3)
ff.set(xlim=(0,30),ylim=(-1,1),xlabel='Top100 overlap');ff.set_yticks([]);ff.set_xticks([0,30]);ff.text(0,1.06,'ρ = 0.555',transform=ff.transAxes,fontsize=5.8)
AUDITS[f]=([c,e],[[c,e]],[], 'Shallow design, nine-identity hero and narrow validation strips intentionally differ in span.')

# FIGURE 3 — two real-data heroes: deposited UMAP and frozen cross-study vectors.
d=F[2];f=canvas('Figure3_v5',200)
xy=pd.read_csv(str(PROJECT/'source_data'/'figure3'/'plot_inputs'/'heatshock_All_reducedDimS.tsv'),sep='\t',index_col=0)
meta=read('external_heatshock/GSE121619/GSE121619_Control_Heatshock_pData.tsv.gz').set_index('cell')
emb=xy.join(meta);assert len(emb)==2085 and emb.louvain_component.notna().all()
src=emb.louvain_component.isin([5,6,8]);tgt=emb.louvain_component.isin([2,3,10]);sealed=tgt&emb.treatment_id.eq('HeatShock');assert sealed.sum()==830
a=panel(f,[.10,.54,.46,.41],'a','Held-out root states')
a.scatter(emb.iloc[:,0],emb.iloc[:,1],c=LIGHT,s=2.4,lw=0,rasterized=True)
a.scatter(emb.loc[src].iloc[:,0],emb.loc[src].iloc[:,1],c=BLUE,s=5,lw=0,rasterized=True)
a.scatter(emb.loc[sealed].iloc[:,0],emb.loc[sealed].iloc[:,1],s=7,fc='white',ec='#68717A',lw=.3,rasterized=True)
for num in [5,6,8,2,3,10]:
 z=emb[emb.louvain_component.eq(num)].iloc[:,:2].median()
 a.text(z.iloc[0],z.iloc[1],str(num),fontsize=6.5,color=INK,ha='center',va='center',bbox=dict(fc='white',ec='none',pad=.7))
a.set(xlabel='Deposited UMAP 1',ylabel='Deposited UMAP 2');a.set_xticks([]);a.set_yticks([])
a.text(0,-.135,'Sources 5/6/8 · targets 2/3/10: 830 HS cells sealed',transform=a.transAxes,fontsize=5.8)
b=panel(f,[.76,.785,.205,.165],'b','Candidate recovery')
p=wide(sel(d,dataset='heatshock_missing_state'),'frozen_reveal_primary_results')
for j,(_,r) in enumerate(p.iterrows()):
 b.hlines(j,0,r.top100_overlap,color=LIGHT,lw=1);b.scatter(r.top100_overlap,j,s=21,c=BLUE,zorder=3)
 b.text(r.top100_overlap+4,j,f'ρ = {r.full_rho:.3f}',fontsize=5.8,va='center')
b.set_yticks(range(3),['Cluster 2','Cluster 3','Cluster 10']);b.invert_yaxis();b.set(xlim=(0,114),xlabel='Top100 overlap');b.margins(y=.25)
c=panel(f,[.76,.535,.205,.14],'c','Shared-response baselines',pad=21)
mods=['signed_mean','signed_median','rank1_common_axis','strongest_source_pretarget','pooled_donor'];labs=['Signed mean','State median','Rank-1 axis','Strongest source','Pooled donors']
for i,m in enumerate(mods):
 v=sel(B,dataset='GSE121619_heatshock',model=m).spearman_full_vector
 c.scatter(v,np.full(len(v),i),s=10,c=MID);c.scatter(np.median(v),i,s=22,c=BLUE,zorder=3)
c.set_yticks(range(5),labs);c.invert_yaxis();c.set(xlim=(.3,.72),xlabel='Full-response ρ')
c.text(0,1.035,'Large points: target median',transform=c.transAxes,fontsize=5.8,color=AXIS)
e=panel(f,[.10,.17,.46,.24],'e','Independent heat-shock transfer')
z=pd.read_csv(ROOT/'crossstudy_plot_values.tsv.gz',sep='\t');assert len(z)==16335 and z.top100_overlap.sum()==64
e.scatter(z.shared_prediction,z.target_response,s=2,c=MID,alpha=.35,lw=0,rasterized=True)
hit=z[z.top100_overlap]
e.scatter(hit.shared_prediction,hit.target_response,s=8,c=BLUE,lw=0,rasterized=True)
zero(e);e.axvline(0,c=LIGHT,ls='--',lw=.6,zorder=0)
e.set(xlabel='Frozen source heat-response score',ylabel='Seedling response (Δ log2[CPM + 1])')
e.text(.025,.96,'ρ = 0.415\nTop100 overlap = 64\nDirection among overlap = 100%',transform=e.transAxes,ha='left',va='top',fontsize=6,linespacing=1.4)
ff=panel(f,[.76,.17,.205,.24],'f','Perturbation specificity')
p=wide(d,'cross_experiment_wrong_perturbation_controls')
for j,(_,r) in enumerate(p.iterrows()):
 ff.plot([r.wrong_program_rho,r.true_heat_rho],[j,j],c=MID,lw=.6);ff.scatter(r.wrong_program_rho,j,c=MID,s=13);ff.scatter(r.true_heat_rho,j,c=BLUE,s=19,zorder=3)
ff.set_yticks(range(7),['Salt 6 h','Salt 3 h','Pep1','SCOOP10','SCREW1','flg22','Chitin']);ff.invert_yaxis();ff.set(xlim=(-.15,.9),xlabel='Paired-universe ρ')
dd=panel(f,[.10,.055,.64,.035],'d','Matched empirical null')
p=wide(d,'heatshock_matched_null')
dd.vlines(range(4),0,p.null_q99,color=LIGHT,lw=4.5);dd.scatter(range(4),p.null_mean,c=MID,s=8);dd.scatter(range(4),p.observed_overlap,c=BLUE,s=18,zorder=3)
dd.set_xticks(range(4),['Cluster 2','Cluster 3','Cluster 10','Seedling']);dd.set(ylim=(0,100),xlim=(-.5,3.5),ylabel='Top100');dd.set_yticks([0,100]);f.text(.81,.061,'4/4 q ≈ 1e-4',fontsize=6)
AUDITS[f]=([a,e],[],[[a,e]], 'The real cell map and gene-wise transfer share a left column; support and boundary have unequal spans.')

# FIGURE 4 — compact type-based framework, then the large nested-target result.
d=F[3];f=canvas('Figure4_v5',180)
S=wide(d,'stimulus_level_source_vs_target');order=['Pep1','SCOOP10','SCREW1','flg22','Chitin'];high=order[:3]
T=wide(sel(d,section='immune_frozen_reveal',context=['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']),idx=['dataset','context'])
a=panel(f,[.065,.70,.40,.25],'a','Evidence eligibility',True)
say(a,0,.90,'What is missing?',ha='left',size=6.5)
for yy,label,res in [(.70,'Future time','Recent same-state history'),(.47,'Unseen state / study','Shared cross-state response')]:
 say(a,0,yy,label,ha='left',size=6);arr(a,(.40,yy),(.51,yy),BLUE);say(a,.55,yy,res,ha='left',color=BLUE,size=6)
arr(a,(.75,.39),(.75,.29),AXIS);say(a,.75,.22,'Source capability',size=6.5)
arr(a,(.68,.16),(.43,.04));arr(a,(.81,.16),(.95,.04))
say(a,.33,-.04,'Fine output',color=BLUE,size=6);say(a,.92,-.04,'Coarse / abstain',color=RUST,size=6)

def labeled(aa,x,y,xlab,ylab,offsets):
 for s in order:
  r=S.loc[s];aa.scatter(r[x],r[y],s=23,c=BLUE if s in high else RUST,marker='o' if s in high else 's',zorder=3)
  dx,dy=offsets[s];aa.annotate(s,(r[x],r[y]),xytext=(dx,dy),textcoords='offset points',fontsize=5.8,ha='left' if dx>=0 else 'right',va='bottom' if dy>=0 else 'top')
 aa.margins(x=.30,y=.28);aa.set(xlabel=xlab,ylabel=ylab)

b=panel(f,[.64,.70,.325,.25],'b','Frozen source stratification')
ob={'Pep1':(5,5),'SCOOP10':(-5,-8),'SCREW1':(5,5),'flg22':(5,3),'Chitin':(5,-6)}
labeled(b,'median_pairwise_rho','effective_state_fraction','Median pairwise source ρ','Effective-state fraction',ob)
c=panel(f,[.11,.315,.52,.255],'c','Nested target performance',pad=21)
for i,s in enumerate(order):
 v=T.loc[s].full_rho;assert len(v)==4
 c.scatter(i+np.linspace(-.13,.13,4),v,c=MID,s=15,zorder=2)
 c.scatter(i,v.mean(),c=BLUE if s in high else RUST,s=30,zorder=3)
c.set_xticks(range(5),order);c.set(ylim=(-.1,1),ylabel='Target full-response ρ',xlim=(-.5,4.5));zero(c)
c.text(0,1.035,'n = 5 stimuli; 4 targets/stimulus',transform=c.transAxes,fontsize=6)
di=panel(f,[.80,.465,.165,.105],'d','Group contrast')
dj=f.add_axes([.80,.30,.165,.095]);dj.tick_params(length=2.5,width=.6)
H=wide(d,'hierarchical_stimulus_inference')
for aa,key,label,lim,ticks in [(di,'full_rho','Δ full-response ρ',(-.05,.45),[0,.2,.4]),(dj,'top100_overlap','Δ Top100 hits',(-5,50),[0,20,40])]:
 r=H.loc[key];aa.errorbar(r.difference_allow_minus_reject,0,xerr=[[r.difference_allow_minus_reject-r.bootstrap_ci025],[r.bootstrap_ci975-r.difference_allow_minus_reject]],fmt='o',c=BLUE,ms=3,capsize=2,lw=.7)
 aa.axvline(0,c=LIGHT,ls='--',lw=.6);aa.set(ylim=(-1,1),xlim=lim,xlabel=label);aa.set_yticks([]);aa.set_xticks(ticks)
di.text(0,.90,'95% CI',transform=di.transAxes,fontsize=5.8)
f.text(.755,.225,'Exact permutation floor P = 0.10',fontsize=5.8)
e=panel(f,[.69,.065,.275,.105],'e','Geometry–amplitude confounding')
oe={'Pep1':(5,4),'SCOOP10':(-6,7),'SCREW1':(5,-6),'flg22':(5,4),'Chitin':(5,-6)}
labeled(e,'shared_rms','rank1_energy','Shared-response RMS','Rank-1 energy',oe)
e.set_ylim(.15,.90)
ff=panel(f,[.10,.075,.40,.055],'f','Matched-null significance',pad=18)
p=wide(d,'immune_matched_null');qstrip(ff,p,'empirical_q_bh_20_targets',[str(i) if i in [1,5,10,15,20] else '' for i in range(1,21)])
ff.text(0,1.07,'19/20 significant',transform=ff.transAxes,fontsize=5.8)
ff.set_xlabel('Held-out target',fontsize=6)
AUDITS[f]=([di,dj],[],[[di,dj]], 'Large nested-target hero, compact framework, boundary scatter and significance strip have intentionally unequal sizes.')

# Final 180 mm render, vector output, 300 dpi previews, and alignment evidence.
hierarchy={}
for name,ff in figs:
 if os.environ.get('V5_EXPORT_ONLY') and name not in os.environ['V5_EXPORT_ONLY'].split(','):continue
 ff.canvas.draw();axes,rows,cols,reason=AUDITS[ff];ids={aa:str(i) for i,aa in enumerate(axes)}
 require_matplotlib_panel_alignment(ff,axes=axes,panel_ids=ids,row_groups=[[ids[a] for a in g] for g in rows],column_groups=[[ids[a] for a in g] for g in cols],json_out=QA/(name+'.alignment.json'),strict=True)
 hierarchy[name]={'intentional_unequal_span_reason':reason,'rectangles':[list(a.get_position().bounds) for a in ff.axes]}
 ff.savefig(ROOT/(name+'.pdf'));ff.savefig(ROOT/(name+'.svg'));ff.savefig(ROOT/(name+'.png'),dpi=300)
 print('EXPORTED',name,flush=True)
with PdfPages(ROOT/'Figures1-4_v5_preview.pdf') as pdf:
 for _,ff in figs:pdf.savefig(ff)
(QA/'layout_hierarchy.json').write_text(json.dumps(hierarchy,indent=2),encoding='utf-8')
for path,digest in INPUTS.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
(QA/'source_input_sha256.json').write_text(json.dumps(INPUTS,ensure_ascii=False,indent=2),encoding='utf-8')
print('Four v5 main-figure previews complete; scientific inputs unchanged.',flush=True)
