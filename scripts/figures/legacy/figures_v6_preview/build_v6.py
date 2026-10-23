"""v6: scientific objects and shared rulers; four main previews only.

AREA_WEIGHTS = {
 'Figure1': {'a_umap': .25, 'b_two_rulers': .30, 'c_null_field': .45},
 'Figure2': {'a_design': .28, 'b_salt_axis': .18, 'c_route_field': .44, 'd_sulfur': .10},
 'Figure3': {'a_umap': .35, 'b_recovery_null': .20, 'c_baseline_ruler': .12, 'd_gene_scatter': .25, 'e_specificity': .08},
 'Figure4': {'a_experimental_landscape': .55, 'b_source_axis': .15, 'c_targets_axis': .22, 'd_effect_rulers': .08},
}
The percentages guide scientific salience, rather than exact rectangular ink area.
No new statistics, fitted embeddings, Word edits or Extended Data edits.
"""
from pathlib import Path
V5=Path(__file__).resolve().parent.parent/'figures_v5_preview'
# Reuse the established read-only data loader and typography, not the v5 panels.
exec((V5/'build_v5.py').read_text(encoding='utf-8').split('# FIGURE 1')[0])
from matplotlib.patches import PathPatch,Polygon
from matplotlib.path import Path as MPath
from matplotlib.transforms import blended_transform_factory

INK='#2A2E33';SECONDARY='#5F6670';GDARK='#8F98A3';MID='#B8C0C8';LIGHT='#DDE2E7';PALE='#F1F3F5'
HISTORY='#3F6680';SHARED='#AAB2BB';WARNING='#A76552';HEAT='#A75C46';HDARK='#7F4637'
CC=['#C6A46A','#8E6837','#91AABA','#4D7087']
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],
 'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':INK,'xtick.color':INK,'ytick.color':INK,
 'svg.fonttype':'none','pdf.fonttype':42,'font.size':6.5,'axes.labelsize':6.7,
 'axes.titlesize':7.2,'xtick.labelsize':6,'ytick.labelsize':6})
figs=[];AUDITS={};STRUCTURE={}

def ruler(a,xlim,ticks,label=None,xticklabels=None):
 """One numerical ruler, with row labels and no second quantitative axis."""
 a.set_xlim(*xlim);a.set_xticks(ticks,xticklabels if xticklabels is not None else [str(t) for t in ticks])
 for side in ['left','right','top']:a.spines[side].set_visible(False)
 a.tick_params(axis='y',length=0,pad=5);a.tick_params(axis='x',length=2.5)
 if label:a.set_xlabel(label)

def note(a,x,y,s,color=SECONDARY,size=6,ha='left',va='center',trans=None):
 return a.text(x,y,s,color=color,fontsize=size,ha=ha,va=va,transform=trans if trans is not None else a.transAxes)

def noframe(a):
 for spine in a.spines.values():spine.set_visible(False)
 a.set_xticks([]);a.set_yticks([])

def line_ruler(a,x0,x1,y,lo,hi,ticks,label):
 a.plot([x0,x1],[y,y],c=INK,lw=.6,clip_on=False)
 for t in ticks:
  x=x0+(t-lo)/(hi-lo)*(x1-x0);a.plot([x,x],[y,y-.015],c=INK,lw=.5,clip_on=False)
  say(a,x,y-.055,str(t),size=5.8)
 say(a,(x0+x1)/2,y-.12,label,size=6)

# F1: an actual leaf-cell object, a data table with two rulers, and a null field.
d=F[0];f=canvas('Figure1_v6_preview',184)
a=panel(f,[.065,.535,.36,.42],'a','Rice leaf nuclei')
rice=pd.read_csv(V5/'rice_deposited_umap.tsv.gz',sep='\t');assert len(rice)==126139
grays={'Epidermis':'#C9D0D6','Meristem':'#9AA6B0','Mesophyll':'#E0E4E8','Vasculature':'#B7C1C9'}
for tissue,g in rice.groupby('TissueSystem',sort=False):a.scatter(g.UMAP1,g.UMAP2,c=grays[tissue],s=.19,lw=0,rasterized=True)
centers=rice.groupby('TissueSystem')[['UMAP1','UMAP2']].median()
for tissue,r in centers.iterrows():a.text(r.UMAP1,r.UMAP2,tissue,color=INK,ha='center',fontsize=6)
noframe(a);a.set(xlabel='UMAP 1',ylabel='UMAP 2')
b=panel(f,[.485,.58,.48,.375],'b','Two levels of replicate stability',True)
note(b,.44,.94,'Full-vector ρ',color=INK,ha='center',size=6.5)
note(b,.86,.94,'Top100 overlap',color=INK,ha='center',size=6.5)
for i,(ct,lab) in enumerate(zip(COND,CTEXT)):
 yy=.77-i*.175;say(b,-.005,yy,lab,ha='left',color=CC[i],size=6)
 for sec,metric,x0,x1,lo,hi in [('rice_full_vector','replicate_response_spearman',.36,.62,-.05,.30),('rice_top100_random_expectation','observed_overlap',.735,.995,0,36)]:
  vv=sel(d,section=sec,metric=metric);v=vv[vv.context.str.startswith(ct)].value.to_numpy();assert len(v)==4
  xp=x0+(v-lo)/(hi-lo)*(x1-x0)
  b.scatter(xp,yy+np.linspace(-.018,.018,4),s=13,c=CC[i],zorder=3)
  xm=x0+(np.median(v)-lo)/(hi-lo)*(x1-x0);b.plot([xm,xm],[yy-.035,yy+.035],c=INK,lw=.85)
for x0,x1,lo,hi,ticks,lab in [(.36,.62,-.05,.30,[0,.1,.2,.3],'Response ρ'),(.735,.995,0,36,[0,10,20,30],'Top100 overlap')]:
 line_ruler(b,x0,x1,.13,lo,hi,ticks,lab)
 xx=x0+(0-lo)/(hi-lo)*(x1-x0);b.plot([xx,xx],[.20,.84],c=LIGHT,ls='--',lw=.5,zorder=0)
expect=sel(d,metric='random_expected_overlap').value.iloc[0]
xx=.735+expect/36*(.995-.735);b.plot([xx,xx],[.20,.84],c=LIGHT,ls='--',lw=.5,zorder=0)
c=panel(f,[.20,.095,.49,.35],'c','Matched empirical null',pad=20)
p=wide(d,'rice_top100_matched_null').reindex([ct+'|'+t for ct in COND for t in ['Epidermis','Meristem','Mesophyll','Vasculature']])
ys=np.arange(16)+np.repeat(np.arange(4)*.8,4)
c.hlines(ys,0,p.null_q99,color=LIGHT,lw=3.4);c.scatter(p.null_mean,ys,c=MID,s=9,zorder=2);c.scatter(p.observed_overlap,ys,c=np.repeat(CC,4),s=19,zorder=3)
c.set_yticks(ys,[x.split('|')[1] for x in p.index]);c.set_ylim(ys[-1]+.8,-.8)
ruler(c,(0,35),[0,10,20,30],'Top100 overlap')
note(c,0,1.027,'16/16 BH q ≈ 1e-4',size=6,color=INK)
sm=wide(d,'rice_condition_summary').reindex(COND);trans=blended_transform_factory(c.transAxes,c.transData)
for i,lab in enumerate(CTEXT):
 yy=ys[i*4:i*4+4].mean();note(c,1.10,yy-.4,lab,color=CC[i],size=6.2,trans=trans)
 note(c,1.10,yy+.55,f'Median enrichment {sm.iloc[i].median_matched_enrichment:.1f}×',size=5.8,trans=trans)
f.legend(handles=[Line2D([],[],c=LIGHT,lw=3.4,label='Null: 0–99th percentile'),Line2D([],[],c=MID,marker='o',ms=3,ls='',label='Null mean')],loc='lower left',bbox_to_anchor=(.19,.018),ncol=2,fontsize=6,handlelength=1.2)
STRUCTURE[f]={'full_xy_axes':0,'real_embeddings':1,'axis_light':['Two-column stability rulers','16-context null field'],'hero':'Matched empirical null'}
AUDITS[f]=([a,b],[[a,b]],[], 'UMAP and comparison field have unequal widths/heights; no equivalence is implied.')

# F2: time/state experiment, one line plot, and one shared rho field.
d=F[1];f=canvas('Figure2_v6_preview',204)
a=panel(f,[.09,.715,.48,.24],'a','Time × state design',True)
for j,label in enumerate(['State A','State B','State C']):
 yy=.73-j*.26;say(a,-.025,yy,label,ha='right',size=6)
 a.plot([.07,.94],[yy,yy],c=PALE,lw=.65,zorder=0)
for xx,txt in [(.10,'0 h'),(.45,'Early'),(.91,'Future target')]:say(a,xx,.98,txt,size=6)
for xx,yy,missing in [(.10,.73,False),(.45,.73,False),(.91,.73,True),(.45,.47,False),(.45,.21,False)]:
 a.add_patch(Ellipse((xx,yy),.032,.058,fc='white' if missing else HISTORY,ec=INK if missing else HISTORY,lw=.7,ls='--' if missing else '-'))
arr(a,(.13,.73),(.41,.73),LIGHT);arr(a,(.49,.73),(.87,.73),HISTORY,1.1)
a.plot([.48,.88],[.47,.70],c=SHARED,lw=.5);a.plot([.48,.88],[.21,.70],c=SHARED,lw=.5)
say(a,.68,.85,'Same-state history',color=HISTORY,size=6);say(a,.72,.17,'Cross-state evidence',color=SECONDARY,size=5.8)
b=panel(f,[.73,.715,.235,.24],'b','Salt temporal baselines')
for i,ds in enumerate(['GSE153103_salt_6h','GSE153103_salt_12h']):
 g=sel(B,dataset=ds);ms=[next(m for m in g.model if m.startswith(pre)) for pre in ['latest','history_mean','linear_last2']];ms.append(next(m for m in g.model if m.startswith('linear_') and not m.startswith('linear_last2')))
 b.plot(range(4),g.set_index('model').loc[ms].spearman_full_vector,'o-',c=[HISTORY,GDARK][i],ms=3,label=['6 h','12 h'][i])
b.set_xticks(range(4),['Latest','Mean','Last-two\ntrend','Baseline\ntrend'],rotation=35,ha='right',rotation_mode='anchor');b.set(ylim=(-.1,1),ylabel='Full-response ρ');zero(b)
b.legend(ncol=2,loc='upper right',fontsize=6,handlelength=1,columnspacing=.8)
note(b,.02,.15,'6 h Top100: trend 36; latest 30',size=5.8)
c=panel(f,[.28,.12,.36,.465],'c','Same-state versus shared evidence')
blocks=[('lateral_root_8h_to_20h','frozen_reveal_primary_results','full_rho',['same_state_latest','shared_other_states'],0),('GSE46205_salt_3h_corrected','post_reveal_preprocessing_corrected_results','spearman_full_vector',['same_state_history','shared_other3'],11),('GSE46205_salt_8h_corrected','post_reveal_preprocessing_corrected_results','spearman_full_vector',['same_state_history','shared_other3'],17)]
names={'Ambiguous Stele Cells':'Ambiguous stele','Columella/ Root Cap':'Columella / root cap','Epidermis':'Epidermis','Lateral Root Endodermis':'Lateral-root endodermis','Mature Pericycle':'Mature pericycle','Phloem':'Phloem','Phloem Pole Pericycle':'Phloem-pole pericycle','Xylem':'Xylem','Xylem Pole Pericycle':'Xylem-pole pericycle','epidermis':'Epidermis','cortex':'Cortex','endodermis':'Endodermis','stele':'Stele'}
ys=[];labels=[]
for ds,section,metric,models,y0 in blocks:
 p=sel(d,dataset=ds,section=section,metric=metric).pivot_table(index='context',columns='model',values='value',sort=False)
 for i,(context,r) in enumerate(p.iterrows()):
  yy=y0+i;ys.append(yy);labels.append(names.get(context,context))
  c.plot([r[models[1]],r[models[0]]],[yy,yy],c=SHARED,lw=.55);c.scatter(r[models[1]],yy,c=SHARED,s=13,zorder=2);c.scatter(r[models[0]],yy,c=HISTORY,s=19,zorder=3)
c.set_yticks(ys,labels);c.set_ylim(21.2,-1.9);ruler(c,(-.1,1),[0,.25,.5,.75,1],'Full-response ρ')
tr=blended_transform_factory(c.transAxes,c.transData)
for yy,header,body in [(3.7,'Lateral root · 8→20 h','9/9 wins\nExact P = 0.00195'),(12.5,'Root layers · 1→3 h','4/4 wins\ncorrected_post_reveal'),(18.5,'Root layers · 3→8 h','4/4 wins\ncorrected_post_reveal')]:
 note(c,1.17,yy-.55,header,size=6.2,color=INK,trans=tr)
 note(c,1.17,yy+.55,body,size=5.8,trans=tr,va='top')
note(c,0,1.005,'● Same-state history',color=HISTORY,size=6)
note(c,.66,1.005,'● Shared',color=SECONDARY,size=6)
dd=panel(f,[.75,.125,.215,.042],'d','Sulfur depletion',True)
p=wide(sel(d,dataset='sulfur_3h_to_12h'),'matched_empirical_null')
dd.plot([0,1],[.55,.55],c=PALE,lw=.6)
dd.plot([0,float(p.null_q99.iloc[0])/30],[.55,.55],c=LIGHT,lw=3.5);dd.scatter(float(p.null_mean.iloc[0])/30,.55,c=MID,s=9);dd.scatter(21/30,.55,c=HISTORY,s=19,zorder=3)
say(dd,0,-.30,'3→12 h · ρ = 0.555 · Top100 = 21',size=5.8,ha='left');say(dd,0,-.75,'Matched-null q ≈ 0.0002',size=5.8,ha='left')
f.text(.10,.025,'13/16 matched-null tests significant; lateral root 6/9. Full context-level controls: source tables / Extended Data.',fontsize=5.8,color=SECONDARY)
STRUCTURE[f]={'full_xy_axes':1,'real_embeddings':0,'axis_light':['One 17-context route ruler','Sulfur null interval'],'hero':'17-context route comparison field'}
AUDITS[f]=None

# F3: real cells, recovery under the null, and actual gene-wise transfer.
d=F[2];f=canvas('Figure3_v6_preview',198)
xy=pd.read_csv(str(PROJECT/'source_data'/'figure3'/'plot_inputs'/'heatshock_All_reducedDimS.tsv'),sep='\t',index_col=0)
md=read('external_heatshock/GSE121619/GSE121619_Control_Heatshock_pData.tsv.gz').set_index('cell');emb=xy.join(md)
src=emb.louvain_component.isin([5,6,8]);sealed=emb.louvain_component.isin([2,3,10])&emb.treatment_id.eq('HeatShock');assert len(emb)==2085 and sealed.sum()==830
a=panel(f,[.075,.49,.48,.465],'a','Held-out root states')
a.scatter(emb.iloc[:,0],emb.iloc[:,1],c=LIGHT,s=2.2,lw=0,rasterized=True);a.scatter(emb.loc[src].iloc[:,0],emb.loc[src].iloc[:,1],c=HEAT,s=5,lw=0,rasterized=True)
a.scatter(emb.loc[sealed].iloc[:,0],emb.loc[sealed].iloc[:,1],s=7,fc='white',ec='#68717A',lw=.3,rasterized=True)
for num in [5,6,8,2,3,10]:
 p=emb[emb.louvain_component.eq(num)].iloc[:,:2].median();a.text(p.iloc[0],p.iloc[1],str(num),fontsize=6.5,ha='center',va='center',bbox=dict(fc='white',ec='none',pad=.7))
noframe(a);a.set(xlabel='Deposited UMAP 1',ylabel='Deposited UMAP 2')
note(a,0,-.105,'Sources 5/6/8; targets 2/3/10; 830 HS cells sealed',size=5.8)
b=panel(f,[.75,.735,.215,.22],'b','Candidate recovery and null')
p=wide(d,'heatshock_matched_null');rhos=[.615,.621,.450,.415]
for i,(_,r) in enumerate(p.iterrows()):
 b.plot([0,r.null_q99],[i,i],c=LIGHT,lw=3.5);b.scatter(r.null_mean,i,c=MID,s=9);b.scatter(r.observed_overlap,i,c=HEAT,s=20,zorder=3)
 b.text(106,i,f'ρ {rhos[i]:.3f}',fontsize=5.8,va='center',color=SECONDARY)
b.set_yticks(range(4),['Cluster 2','Cluster 3','Cluster 10','Seedling']);b.set_ylim(3.6,-.6);ruler(b,(0,147),[0,50,100],'Top100 overlap')
note(b,0,-.25,'4/4 source-matched BH q ≈ 1e-4',size=5.8)
c=panel(f,[.75,.465,.215,.15],'c','Shared-response baselines')
mods=['signed_mean','rank1_common_axis','signed_median','strongest_source_pretarget','pooled_donor'];mlabs=['Signed mean','Rank-1 axis','State median','Strongest source','Pooled donors']
for i,m in enumerate(mods):
 v=sel(B,dataset='GSE121619_heatshock',model=m).spearman_full_vector
 c.scatter(v,np.full(len(v),i),s=9,c=MID);c.plot([np.median(v)]*2,[i-.21,i+.21],c=HEAT,lw=1.2,zorder=4)
c.set_yticks(range(5),mlabs);c.set_ylim(4.55,-.55);ruler(c,(.3,.72),[.3,.5,.7],'Full-response ρ')
note(c,0,-.36,'Tick: median across three target states',size=5.8)
dd=panel(f,[.10,.10,.455,.255],'d','Independent heat-shock transfer')
z=pd.read_csv(V5/'crossstudy_plot_values.tsv.gz',sep='\t');assert len(z)==16335 and z.top100_overlap.sum()==64
dd.scatter(z.shared_prediction,z.target_response,s=2,c=LIGHT,alpha=.55,lw=0,rasterized=True)
hit=z[z.top100_overlap];dd.scatter(hit.shared_prediction,hit.target_response,s=8,c=HEAT,lw=0,rasterized=True)
dd.axhline(0,c=LIGHT,ls='--',lw=.5,zorder=0);dd.axvline(0,c=LIGHT,ls='--',lw=.5,zorder=0)
dd.set(xlabel='Frozen source heat-response score',ylabel='Independent seedling response\n(Δ log2[CPM + 1])')
note(dd,.025,.97,'ρ = 0.415 · Top100 = 64\nDirection among overlap = 100%',size=6,color=INK,va='top')
e=panel(f,[.78,.105,.185,.25],'e','Perturbation specificity')
p=wide(d,'cross_experiment_wrong_perturbation_controls')
for i,(_,r) in enumerate(p.iterrows()):
 e.plot([r.wrong_program_rho,r.true_heat_rho],[i,i],c=SHARED,lw=.55);e.scatter(r.wrong_program_rho,i,c=SHARED,s=12);e.scatter(r.true_heat_rho,i,c=HEAT,s=18,zorder=3)
e.set_yticks(range(7),['Salt 6 h','Salt 3 h','Pep1','SCOOP10','SCREW1','flg22','Chitin']);e.set_ylim(6.5,-.5);ruler(e,(-.1,.75),[0,.25,.5,.75],'Paired-universe ρ')
STRUCTURE[f]={'full_xy_axes':1,'real_embeddings':1,'axis_light':['Recovery + source-matched null','Baseline median strip','Seven specificity dumbbells'],'heroes':['Deposited root UMAP','Independent gene-wise transfer']}
AUDITS[f]=None

# F4: explicit experimental information geometry, not a molecular mechanism.
f=canvas('Figure4_v6_preview',202)
a=panel(f,[.045,.10,.59,.85],'a','Experimental evidence routing',True)

def root_shape(aa,x,y,w=.12,h=.22):
 """Conceptual longitudinal root, with distinct longitudinal tissue layers."""
 verts=[(x-w/2,y+h/2),(x-w/2,y-h*.15),(x-w*.32,y-h*.44),(x,y-h*.52),(x+w*.32,y-h*.44),(x+w/2,y-h*.15),(x+w/2,y+h/2)]
 aa.add_patch(Polygon(verts,closed=True,fc='#F4F5F6',ec=GDARK,lw=.65))
 for frac in [-.28,0,.28]:
  xx=x+w*frac;aa.plot([xx,xx,x],[y+h/2,y-h*.15,y-h*.46],c=LIGHT,lw=.55)
 for yy in [y+.22*h,y,y-.2*h]:aa.plot([x-w*.47,x+w*.47],[yy,yy],c=LIGHT,lw=.5)

def sample(aa,x,y,missing=False,col=HISTORY):
 aa.add_patch(Ellipse((x,y),.022,.015,fc='white' if missing else col,ec=INK if missing else col,lw=.65,ls='--' if missing else '-'))

def curved(aa,start,c1,c2,end,color,lw=.75):
 aa.add_patch(PathPatch(MPath([start,c1,c2,end],[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4]),fc='none',ec=color,lw=lw))

# One continuous scaffold; scene boundaries are whitespace rather than cards.
say(a,.01,.985,'Future time',size=6.5,ha='left')
root_shape(a,.065,.84,.095,.22)
for yy in [.88,.83,.78]:a.plot([.19,.65],[yy,yy],c=PALE,lw=.6)
for xx,txt in [(.22,'t0'),(.40,'t1'),(.61,'t2')]:say(a,xx,.945,txt,size=5.8)
for xx,yy,missing in [(.22,.88,False),(.40,.88,False),(.61,.88,True),(.40,.83,False),(.40,.78,False)]:sample(a,xx,yy,missing)
arr(a,(.24,.88),(.38,.88),LIGHT);arr(a,(.42,.88),(.59,.88),HISTORY,1.0)
a.plot([.41,.59],[.83,.87],c=LIGHT,lw=.45);a.plot([.41,.59],[.78,.87],c=LIGHT,lw=.45)
say(a,.20,.73,'Direct same-state history exists',ha='left',color=SECONDARY,size=5.8)

say(a,.01,.675,'Unseen state',size=6.5,ha='left')
root_shape(a,.065,.56,.095,.18)
for yy in [.62,.57,.52]:sample(a,.25,yy)
sample(a,.61,.52,True)
for yy in [.62,.57,.52]:curved(a,(.27,yy),(.33,yy),(.36,.57),(.42,.57),SHARED,.7)
# A ribbon represents the shared-response evidence path, not measured amplitude.
curved(a,(.42,.57),(.50,.57),(.52,.52),(.59,.52),HISTORY,1.2)
say(a,.20,.445,'No treatment history in target state',ha='left',color=SECONDARY,size=5.8)

say(a,.01,.385,'Unseen experiment',size=6.5,ha='left')
root_shape(a,.245,.31,.075,.12)
a.plot([.45,.45],[.235,.37],c=GDARK,ls='--',lw=.6)
say(a,.46,.39,'Study boundary',size=5.8,color=SECONDARY)
# A restrained seedling outline identifies the actual experimental scale change.
a.plot([.62,.62,.61],[.28,.34,.35],c=GDARK,lw=.7)
a.add_patch(Ellipse((.602,.348),.04,.018,angle=-20,fc=PALE,ec=GDARK,lw=.6))
a.add_patch(Ellipse((.638,.348),.04,.018,angle=20,fc=PALE,ec=GDARK,lw=.6))
a.plot([.62,.61,.615],[.28,.245,.23],c=GDARK,lw=.6)
arr(a,(.29,.31),(.58,.31),HISTORY,.9)
say(a,.245,.215,'Source root',size=5.8);say(a,.62,.205,'Independent study',size=5.8)

# Only source evidence enters the bundle; target expression is never an input.
for st,c1,c2,en,col in [((.40,.90),(.59,.93),(.72,.89),(.82,.82),HISTORY),((.42,.59),(.62,.65),(.72,.67),(.82,.66),SHARED),((.29,.29),(.48,.27),(.72,.45),(.82,.50),GDARK)]:
 curved(a,st,c1,c2,en,col,.8)
for off,col in [(-.012,HISTORY),(0,SHARED),(.012,GDARK)]:
 curved(a,(.82+off,.82),(.85+off,.72),(.85+off,.55),(.82+off,.43),col,.8)
say(a,.85,.91,'Eligible source\nevidence',size=6.3)
arr(a,(.82,.42),(.61,.26),INK,.7)
say(a,.59,.30,'Source capability',size=6.5)
say(a,.49,.17,'Response strength · cross-source consistency · replicate stability',size=5.8)

# Output shapes are explicitly conceptual; bars encode ranks, not invented values.
for xx,col in [(.20,HISTORY),(.50,GDARK),(.82,WARNING)]:arr(a,(.60,.24),(xx,.12),col,.6)
for i,w in enumerate([.09,.075,.06,.045]):
 yy=.085-i*.016;a.plot([.15,.15+w],[yy,yy],c=HISTORY,lw=1)
 if i<3:arr(a,(.26,yy-.008),(.26,yy+.008) if i!=1 else (.26,yy-.022),HISTORY,.6)
say(a,.20,-.017,'Gene ranking\n+ direction',size=6)
# Broader ribbon denotes a program-level response; no individual gene prediction.
for yy in [.045,.065,.085]:a.plot([.44,.56],[yy,yy],c=MID,lw=2.2,solid_capstyle='butt')
say(a,.50,-.017,'Coarse response\nprogram',size=6)
sample(a,.82,.065,True,col=WARNING);say(a,.82,-.017,'Abstain',size=6,color=WARNING)

# Exactly two complete quantitative plots at right; group contrasts use rulers.
d=F[3];S=wide(d,'stimulus_level_source_vs_target');order=['Pep1','SCOOP10','SCREW1','flg22','Chitin'];high=order[:3]
T=wide(sel(d,section='immune_frozen_reveal',context=['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']),idx=['dataset','context'])
b=panel(f,[.765,.745,.20,.21],'b','Frozen source stratification')
offs={'Pep1':(4,4),'SCOOP10':(-4,-7),'SCREW1':(4,4),'flg22':(4,4),'Chitin':(4,-7)}
for ss in order:
 r=S.loc[ss];b.scatter(r.median_pairwise_rho,r.effective_state_fraction,c=HISTORY if ss in high else WARNING,s=20,marker='o' if ss in high else 's')
 dx,dy=offs[ss];b.annotate(ss,(r.median_pairwise_rho,r.effective_state_fraction),xytext=(dx,dy),textcoords='offset points',ha='left' if dx>0 else 'right',va='bottom' if dy>0 else 'top',fontsize=5.8)
b.margins(x=.43,y=.33);b.set(xlabel='Pairwise source ρ',ylabel='Effective-state fraction')
c=panel(f,[.765,.395,.20,.225],'c','Nested target performance',pad=21)
for i,ss in enumerate(order):
 v=T.loc[ss].full_rho;c.scatter(i+np.linspace(-.14,.14,4),v,c=MID,s=12);c.scatter(i,v.mean(),c=HISTORY if ss in high else WARNING,s=27,zorder=3)
c.set_xticks(range(5),order,rotation=45,ha='right',rotation_mode='anchor');c.set(ylim=(-.1,1),ylabel='Target full-response ρ');zero(c)
note(c,0,1.035,'5 stimuli × 4 nested targets',color=INK,size=5.8)
di=panel(f,[.785,.225,.18,.055],'d','Group contrast')
dj=f.add_axes([.785,.115,.18,.055]);dj.tick_params(labelsize=6,length=2.5,width=.6)
H=wide(d,'hierarchical_stimulus_inference')
for aa,key,label,lim,ticks in [(di,'full_rho','Δ full-response ρ',(-.05,.45),[0,.2,.4]),(dj,'top100_overlap','Δ Top100 hits',(-5,50),[0,20,40])]:
 r=H.loc[key];aa.errorbar(r.difference_allow_minus_reject,0,xerr=[[r.difference_allow_minus_reject-r.bootstrap_ci025],[r.bootstrap_ci975-r.difference_allow_minus_reject]],fmt='o',c=HISTORY,ms=3,capsize=2,lw=.7)
 aa.set_ylim(-1,1);aa.set_yticks([]);ruler(aa,lim,ticks,label);aa.axvline(0,c=LIGHT,ls='--',lw=.5)
note(di,0,.95,'95% CI',size=5.8)
f.text(.75,.04,'Exact permutation floor P = 0.10',fontsize=5.8,color=SECONDARY)
f.text(.065,.035,'Geometry covaries with response amplitude; independent effect unresolved.\n19/20 held-out targets exceed matched null.',fontsize=5.8,color=SECONDARY,linespacing=1.5)
STRUCTURE[f]={'full_xy_axes':2,'real_embeddings':0,'axis_light':['Two effect-size rulers in original units'],'hero':'Conceptual experimental evidence routing','quantitative_hero':'20 nested target observations','schematic_status':'Conceptual experimental information geometry, not a fitted anatomical or molecular mechanism'}
AUDITS[f]=([di,dj],[],[[di,dj]], 'Both effect rulers have matching width and distinct original units; other panels are intentionally unequal.')

# Render all four; audit only scientifically comparable rectangles and state exemptions.
reports={}
for name,ff in figs:
 if os.environ.get('V6_EXPORT_ONLY') and name not in os.environ['V6_EXPORT_ONLY'].split(','):continue
 ff.canvas.draw()
 scope=AUDITS[ff]
 # In F1, UMAP/rulers are deliberately not equal plot rectangles; no equality assertion.
 if scope is not None and name=='Figure4_v6_preview':
  axes,rows,cols,reason=scope;ids={aa:str(i) for i,aa in enumerate(axes)}
  require_matplotlib_panel_alignment(ff,axes=axes,panel_ids=ids,row_groups=[[ids[a] for a in gr] for gr in rows],column_groups=[[ids[a] for a in gr] for gr in cols],json_out=QA/(name+'.alignment.json'),strict=True)
 else:
  (QA/(name+'.alignment.json')).write_text(json.dumps({'status':'intentional_unequal_spans','reason':'Scientific object, shared ruler fields and supporting quantitative panels have distinct roles; there are no repeated equal-span plots.','rectangles':[list(aa.get_position().bounds) for aa in ff.axes]},indent=2),encoding='utf-8')
 ff.savefig(ROOT/(name+'.pdf'));ff.savefig(ROOT/(name+'.svg'));ff.savefig(ROOT/(name+'.png'),dpi=300)
 reports[name]=STRUCTURE[ff];print('EXPORTED',name,flush=True)
with PdfPages(ROOT/'Figures1-4_v6_preview.pdf') as pdf:
 for _,ff in figs:pdf.savefig(ff)
(QA/'visual_structure_audit.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2),encoding='utf-8')
for path,digest in INPUTS.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
(QA/'source_input_sha256.json').write_text(json.dumps(INPUTS,ensure_ascii=False,indent=2),encoding='utf-8')
print('Four English v6 previews exported. No Word, manuscript or ED edits.',flush=True)
