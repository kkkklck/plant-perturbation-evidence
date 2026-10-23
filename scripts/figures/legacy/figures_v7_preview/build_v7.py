"""Compact scientific narrative, main figures only; old artifacts remain intact.

AREA_WEIGHTS = {
 'Figure1': {'real_umap':.33,'two_level_field':.29,'null_evidence':.38},
 'Figure2': {'design':.16,'salt_axis':.22,'route_atlas':.56,'sulfur_inset':.06},
 'Figure3': {'umap':.32,'recovery_null':.16,'baseline':.10,'gene_transfer':.30,'specificity':.12},
 'Figure4': {'context_to_capability_map':.58,'source_stratification':.19,'nested_targets_with_CI':.23},
}
Reading paths: F1 biological source → two levels of stability → null evidence;
F2 temporal design → baseline → common 17-context route atlas;
F3 held-out cell map → recovery/null → simple baselines → transfer → specificity;
F4 three observed contexts → four source-only channels → capability gate → resolution,
with the five-stimulus calibration immediately adjacent.
"""
from pathlib import Path
V6=Path(__file__).resolve().parent.parent/'figures_v6_preview'
exec((V6/'build_v6.py').read_text(encoding='utf-8').split('# F1')[0])
NAVY='#496A82';DUST='#7896AA';OCHRE='#B58A4A';DOCHRE='#8E6C35';BRICK='#A45F49'
INK='#2E3135';SECONDARY='#5C636B';MID='#A9B0B8';LIGHT='#DCE1E5';PALE='#F2F4F5'
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],
 'svg.fonttype':'none','pdf.fonttype':42,'font.size':7,'axes.titlesize':7.5,'axes.labelsize':7,
 'xtick.labelsize':6.5,'ytick.labelsize':6.5,'text.color':INK,'axes.edgecolor':INK,'axes.labelcolor':INK})
figs=[];STRUCTURE={};SCOPES={}

def p7(f,r,letter,title,off=False,pad=9):
 a=f.add_axes(r);a.set_title(title,loc='left',fontsize=7.5,pad=pad,weight='normal')
 a.annotate(letter,(0,1),xycoords='axes fraction',xytext=(-16,pad),textcoords='offset points',fontsize=9,weight='bold',va='baseline')
 a.tick_params(length=2.5,width=.6)
 if off:a.set(xlim=(0,1),ylim=(0,1));a.axis('off')
 return a

def t7(a,x,y,s,size=6.5,color=INK,ha='left',va='center',tr=None):
 return a.text(x,y,s,fontsize=size,color=color,ha=ha,va=va,transform=tr if tr is not None else a.transAxes)

def curve7(a,st,c1,c2,en,color=SECONDARY,lw=.8):
 a.add_patch(PathPatch(MPath([st,c1,c2,en],[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4]),fc='none',ec=color,lw=lw))

def f2():
 d=F[1];f=canvas('Figure2_v7_preview',150)
 a=p7(f,[.075,.715,.43,.23],'a','Future-time design',True)
 for i,label in enumerate(['State A','State B','State C']):
  yy=.72-i*.25;t7(a,-.035,yy,label,size=6.5,ha='right');a.plot([.07,.94],[yy,yy],c=PALE,lw=.6)
 for xx,txt in [(.10,'0 h'),(.44,'Early'),(.91,'Future')]:t7(a,xx,.99,txt,ha='center')
 for xx,yy,target in [(.10,.72,False),(.44,.72,False),(.91,.72,True),(.44,.47,False),(.44,.22,False)]:
  a.add_patch(Ellipse((xx,yy),.029,.058,fc='white' if target else NAVY,ec=INK if target else NAVY,lw=.7,ls='--' if target else '-'))
 arr(a,(.13,.72),(.41,.72),LIGHT,.7);arr(a,(.48,.72),(.88,.72),NAVY,1.1)
 a.plot([.47,.88],[.47,.69],c=MID,lw=.55);a.plot([.47,.88],[.22,.69],c=MID,lw=.55)
 t7(a,.64,.85,'Same-state history',color=NAVY,ha='center');t7(a,.73,.15,'Cross-state',color=SECONDARY,ha='center')
 b=p7(f,[.69,.715,.27,.23],'b','Salt temporal baselines')
 for i,ds in enumerate(['GSE153103_salt_6h','GSE153103_salt_12h']):
  gg=sel(B,dataset=ds);ms=[next(m for m in gg.model if m.startswith(pre)) for pre in ['latest','history_mean','linear_last2']];ms.append(next(m for m in gg.model if m.startswith('linear_') and not m.startswith('linear_last2')))
  b.plot(range(4),gg.set_index('model').loc[ms].spearman_full_vector,'o-',c=[NAVY,MID][i],ms=3.5,lw=.9,label=['6 h','12 h'][i])
 b.set_xticks(range(4),['Latest','Mean','Last 2','+ baseline'],rotation=30,ha='right',rotation_mode='anchor');b.set(ylim=(-.1,1),ylabel='Full-response ρ');zero(b)
 b.legend(ncol=2,loc='upper right',fontsize=6.5,handlelength=1,columnspacing=.8)
 t7(b,.025,.15,'6 h Top100: trend 36; latest 30',size=6.5,color=SECONDARY)
 c=p7(f,[.26,.13,.375,.425],'c','Route comparison atlas',pad=18)
 blocks=[('lateral_root_8h_to_20h','frozen_reveal_primary_results','full_rho',['same_state_latest','shared_other_states'],0),('GSE46205_salt_3h_corrected','post_reveal_preprocessing_corrected_results','spearman_full_vector',['same_state_history','shared_other3'],10.5),('GSE46205_salt_8h_corrected','post_reveal_preprocessing_corrected_results','spearman_full_vector',['same_state_history','shared_other3'],16)]
 names={'Ambiguous Stele Cells':'Ambiguous stele','Columella/ Root Cap':'Columella / root cap','Epidermis':'Epidermis','Lateral Root Endodermis':'Lateral-root endodermis','Mature Pericycle':'Mature pericycle','Phloem':'Phloem','Phloem Pole Pericycle':'Phloem-pole pericycle','Xylem':'Xylem','Xylem Pole Pericycle':'Xylem-pole pericycle','epidermis':'Epidermis','cortex':'Cortex','endodermis':'Endodermis','stele':'Stele'}
 yy=[];labs=[]
 for ds,sec,metric,models,y0 in blocks:
  p=sel(d,dataset=ds,section=sec,metric=metric).pivot_table(index='context',columns='model',values='value',sort=False)
  for i,(ct,r) in enumerate(p.iterrows()):
   y=y0+i;yy.append(y);labs.append(names[ct]);c.plot([r[models[1]],r[models[0]]],[y,y],c=MID,lw=.9);c.scatter(r[models[1]],y,c=MID,s=13,zorder=2);c.scatter(r[models[0]],y,c=NAVY,s=22,zorder=3)
 c.set_yticks(yy,labs);c.set_ylim(19.6,-1.5);ruler(c,(-.1,1),[0,.5,1],'Full-response ρ')
 t7(c,0,1.025,'● Same-state history',color=NAVY);t7(c,.69,1.025,'● Shared',color=SECONDARY)
 tr=blended_transform_factory(c.transAxes,c.transData)
 for yy,header,body in [(3.2,'Lateral root | 8→20 h','9/9 wins; exact P = 0.00195'),(11.2,'Root layers | 1→3 h','4/4 wins; corrected_post_reveal'),(17.0,'Root layers | 3→8 h','4/4 wins; corrected_post_reveal')]:
  t7(c,1.14,yy,header,size=7,tr=tr)
  t7(c,1.14,yy+1.05,body,size=6.5,tr=tr,color=SECONDARY)
 # Sulfur is an attached evidence item, without a separate letter or axes.
 inset=f.add_axes([.70,.56,.26,.065]);inset.set(xlim=(0,30),ylim=(0,1));inset.axis('off')
 p=wide(sel(d,dataset='sulfur_3h_to_12h'),'matched_empirical_null')
 inset.plot([0,float(p.null_q99.iloc[0])],[.51,.51],c=LIGHT,lw=3.5);inset.scatter(float(p.null_mean.iloc[0]),.51,c=MID,s=10);inset.scatter(21,.51,c=NAVY,s=22,zorder=3)
 t7(inset,0,1.12,'Sulfur 3→12 h',size=7)
 t7(inset,0,-.12,'ρ = 0.555 · Top100 = 21 · q ≈ 0.0002',size=6.5,color=SECONDARY)
 f.text(.11,.028,'13/16 matched-null tests significant; lateral root 6/9.',fontsize=6.5,color=SECONDARY)
 STRUCTURE[f]={'full_xy_axes':1,'shared_rulers':1,'contexts_in_atlas':17,'insets':['Sulfur, no independent panel'],'hero':'17-context route comparison atlas','traditional_panels':3}
 SCOPES[f]=None

def f4():
 f=canvas('Figure4_v7_preview',150)
 a=p7(f,[.04,.12,.565,.82],'a','Experimental evidence routing',True)
 # Context 1: observed same-state trajectory; direct and alternative routes differ.
 t7(a,0,.975,'Future time',size=7.5)
 for xx,lab in [(.055,'t0'),(.235,'t1'),(.43,'target')]:t7(a,xx,.90,lab,ha='center')
 for xx,yy,target in [(.055,.83,False),(.235,.83,False),(.43,.83,True),(.235,.765,False)]:
  a.add_patch(Ellipse((xx,yy),.025,.019,fc='white' if target else NAVY,ec=INK if target else NAVY,lw=.7,ls='--' if target else '-'))
 arr(a,(.075,.83),(.215,.83),LIGHT,.7);arr(a,(.255,.83),(.41,.83),NAVY,1.1)
 a.plot([.255,.41],[.765,.82],c=MID,lw=.55)
 t7(a,0,.705,'Own observed history',size=6.5,color=SECONDARY)
 curve7(a,(.255,.845),(.40,.87),(.43,.69),(.50,.69),NAVY,1.0)
 # Context 2: the actual held-out root map, rather than a generic cluster icon.
 t7(a,0,.625,'Unseen state',size=7.5)
 xy=pd.read_csv(str(PROJECT/'source_data'/'figure3'/'plot_inputs'/'heatshock_All_reducedDimS.tsv'),sep='\t',index_col=0)
 md=pd.read_csv(DATA/'external_heatshock/GSE121619/GSE121619_Control_Heatshock_pData.tsv.gz',sep='\t').set_index('cell');em=xy.join(md)
 vals=em.iloc[:,:2].to_numpy();mins=vals.min(axis=0);span=np.ptp(vals,axis=0);norm=(vals-mins)/span
 xx=.015+norm[:,0]*.275;yy=.355+norm[:,1]*.215
 src=em.louvain_component.isin([5,6,8]);sealed=em.louvain_component.isin([2,3,10])&em.treatment_id.eq('HeatShock')
 a.scatter(xx,yy,c=LIGHT,s=.5,lw=0,rasterized=True);a.scatter(xx[src],yy[src],c=NAVY,s=1,lw=0,rasterized=True);a.scatter(xx[sealed],yy[sealed],fc='white',ec=SECONDARY,s=1,lw=.15,rasterized=True)
 t7(a,.01,.325,'Sources 5/6/8 → sealed 2/3/10',size=6.5,color=SECONDARY)
 curve7(a,(.30,.49),(.39,.49),(.43,.59),(.50,.59),MID,.9)
 # Context 3: actual source/target experiment conditions and the frozen boundary.
 t7(a,0,.265,'Unseen experiment',size=7.5)
 t7(a,.005,.19,'Root scRNA-seq\n38°C, 45 min',size=6.5)
 t7(a,.32,.19,'Whole seedlings\n37°C, 1 h',size=6.5)
 a.plot([.28,.28],[.08,.245],c=MID,lw=.65,ls='--')
 arr(a,(.12,.10),(.43,.10),NAVY,.9)
 t7(a,.27,.035,'Study boundary · frozen ranking',ha='center',size=6.5,color=SECONDARY)
 curve7(a,(.44,.115),(.56,.17),(.44,.35),(.50,.39),SECONDARY,.8)
 # Four source-only channels form the capability sieve; no measured values invented.
 t7(a,.66,.855,'Source-only gate',size=7,ha='center')
 criteria=['Response magnitude','Cross-state agreement','Replicate stability','Source coverage']
 for lab,y in zip(criteria,[.735,.635,.535,.435]):
  t7(a,.53,y,lab,size=6.5)
  a.plot([.52,.80],[y-.045,y-.045],c=MID,lw=.65)
 for x in [.805,.817,.829]:a.plot([x,x+.013],[.76,.34],c=MID,lw=.65)
 # Resolution classes are the outputs, rather than large illustrative icons.
 for y,end,lw,col in [(.69,.76,1.2,NAVY),(.55,.51,.8,SECONDARY),(.40,.28,.55,MID)]:
  curve7(a,(.835,y),(.89,y),(.90,end),(.96,end),col,lw)
 for y,head,body,col in [(.82,'High support','Ranked genes\n+ direction',NAVY),(.57,'Intermediate','Coarse program /\ncandidate class',SECONDARY),(.34,'Insufficient','Abstain',SECONDARY)]:
  t7(a,1.0,y,head,size=6.5,ha='center',color=col)
  t7(a,1.0,y-.10,body,size=6.5,ha='center',color=col)
 # Real five-stimulus calibration beside the conceptual map.
 d=F[3];S=wide(d,'stimulus_level_source_vs_target');order=['Pep1','SCOOP10','SCREW1','flg22','Chitin'];high=order[:3]
 T=wide(sel(d,section='immune_frozen_reveal',context=['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']),idx=['dataset','context'])
 b=p7(f,[.755,.69,.215,.25],'b','Frozen source stratification')
 offs={'Pep1':(4,-3),'SCOOP10':(5,6),'SCREW1':(4,-6),'flg22':(4,4),'Chitin':(4,-6)}
 for s in order:
  r=S.loc[s];b.scatter(r.median_pairwise_rho,r.effective_state_fraction,c=NAVY if s in high else BRICK,s=23,marker='o' if s in high else 's')
  dx,dy=offs[s];b.annotate(s,(r.median_pairwise_rho,r.effective_state_fraction),xytext=(dx,dy),textcoords='offset points',ha='left' if dx>0 else 'right',va='bottom' if dy>0 else 'top',fontsize=6.5)
 b.margins(x=.55,y=.38);b.set(xlabel='Pairwise source ρ',ylabel='Effective-state fraction')
 c=p7(f,[.72,.22,.165,.28],'c','Nested target performance',pad=18)
 for i,s in enumerate(order):
  v=T.loc[s].full_rho;c.scatter(i+np.linspace(-.14,.14,4),v,c=MID,s=13);c.scatter(i,v.mean(),c=NAVY if s in high else BRICK,s=26,zorder=3)
 c.set_xticks(range(5),order,rotation=90,ha='right',rotation_mode='anchor');c.set(ylim=(-.1,1),ylabel='Target full-response ρ');zero(c)
 t7(c,0,1.055,'5 stimuli × 4 targets',size=6.5)
 ci1=f.add_axes([.925,.405,.055,.075]);ci2=f.add_axes([.925,.255,.055,.075]);H=wide(d,'hierarchical_stimulus_inference')
 for aa,key,lim,ticks,label in [(ci1,'full_rho',(-.05,.45),[0,.4],'Δρ'),(ci2,'top100_overlap',(-5,50),[0,40],'ΔTop100')]:
  rr=H.loc[key];aa.errorbar(rr.difference_allow_minus_reject,0,xerr=[[rr.difference_allow_minus_reject-rr.bootstrap_ci025],[rr.bootstrap_ci975-rr.difference_allow_minus_reject]],fmt='o',c=NAVY,ms=3,capsize=2,lw=.8)
  aa.set(ylim=(-1,1));aa.set_yticks([]);ruler(aa,lim,ticks);aa.tick_params(labelsize=6.5);aa.axvline(0,c=LIGHT,lw=.55,ls='--')
  t7(aa,.5,1.30,label,size=6.5,ha='center')
 f.text(.905,.14,'95% CI\nP floor = 0.10',fontsize=6.5,color=SECONDARY,linespacing=1.4)
 f.text(.055,.035,'Geometry–amplitude covariation unresolved; 19/20 targets exceed matched null.',fontsize=6.5,color=SECONDARY)
 STRUCTURE[f]={'full_xy_axes':2,'shared_rulers':2,'CI_insets':2,'real_cell_map_in_context':2085,'hero':'Evidence-routing map with four observable source channels','traditional_panels':3,'conceptual_gate':'No numerical decision thresholds fitted; source criteria and output classes are conceptual.'}
 SCOPES[f]=([ci1,ci2],[],[[ci1,ci2]],'Two compact CI rulers preserve their original units and align beside the nested-target data.')

def f1():
 # 144 × 142 mm: the information determines the canvas, rather than full-page fill.
 d=F[0];f=plt.figure(figsize=(144/25.4,142/25.4),facecolor='white');figs.append(('Figure1_v7_preview',f))
 cc=[OCHRE,DOCHRE,DUST,NAVY]
 a=p7(f,[.065,.565,.34,.38],'a','Rice leaf nuclei')
 rice=pd.read_csv(V5/'rice_deposited_umap.tsv.gz',sep='\t');assert len(rice)==126139
 grays={'Epidermis':'#C6CFD6','Meristem':'#98A6B0','Mesophyll':'#E0E5E9','Vasculature':'#B4C0C8'}
 for tissue,g in rice.groupby('TissueSystem',sort=False):a.scatter(g.UMAP1,g.UMAP2,c=grays[tissue],s=.20,lw=0,rasterized=True)
 centers=rice.groupby('TissueSystem')[['UMAP1','UMAP2']].median()
 for tissue,r in centers.iterrows():a.text(r.UMAP1,r.UMAP2,tissue,ha='center',fontsize=6.5,color=INK)
 noframe(a);a.set(xlabel='UMAP 1',ylabel='UMAP 2')
 b=p7(f,[.47,.60,.49,.345],'b','Two-level stability',True)
 for x,lab in [(.49,'Full-vector ρ'),(.88,'Top100')]:t7(b,x,.94,lab,ha='center')
 for i,(ct,lab) in enumerate(zip(COND,CTEXT)):
  yy=.77-i*.18;t7(b,-.02,yy,lab,color=cc[i],size=6.5)
  for sec,metric,x0,x1,lo,hi in [('rice_full_vector','replicate_response_spearman',.37,.63,-.05,.30),('rice_top100_random_expectation','observed_overlap',.76,.995,0,36)]:
   g=sel(d,section=sec,metric=metric);v=g[g.context.str.startswith(ct)].value.to_numpy();assert len(v)==4
   xp=x0+(v-lo)/(hi-lo)*(x1-x0);b.scatter(xp,yy+np.linspace(-.018,.018,4),s=13,c=cc[i],zorder=3)
   xm=x0+(np.median(v)-lo)/(hi-lo)*(x1-x0);b.plot([xm,xm],[yy-.035,yy+.035],c=INK,lw=.9)
 for x0,x1,lo,hi,ticks in [(.37,.63,-.05,.30,[0,.15,.3]),(.76,.995,0,36,[0,15,30])]:
  b.plot([x0,x1],[.12,.12],c=INK,lw=.65)
  for t in ticks:
   xx=x0+(t-lo)/(hi-lo)*(x1-x0);b.plot([xx,xx],[.12,.10],c=INK,lw=.6);t7(b,xx,.025,str(t),ha='center')
 xx=.37+.05/.35*(.63-.37);b.plot([xx,xx],[.19,.83],c=LIGHT,lw=.55,ls='--')
 expected=sel(d,metric='random_expected_overlap').value.iloc[0];xx=.76+expected/36*(.995-.76);b.plot([xx,xx],[.19,.83],c=LIGHT,lw=.55,ls='--')
 c=p7(f,[.205,.105,.49,.365],'c','Matched empirical null',pad=18)
 p=wide(d,'rice_top100_matched_null').reindex([ct+'|'+t for ct in COND for t in ['Epidermis','Meristem','Mesophyll','Vasculature']])
 yy=np.arange(16)+np.repeat(np.arange(4)*.7,4)
 c.hlines(yy,0,p.null_q99,color=LIGHT,lw=3.2);c.scatter(p.null_mean,yy,c=MID,s=8,zorder=2);c.scatter(p.observed_overlap,yy,c=np.repeat(cc,4),s=20,zorder=3)
 c.set_yticks(yy,[x.split('|')[1] for x in p.index]);c.set_ylim(yy[-1]+.65,-.65);ruler(c,(0,35),[0,10,20,30],'Top100 overlap')
 t7(c,0,1.025,'16/16 BH q ≈ 1e-4',size=6.5)
 sm=wide(d,'rice_condition_summary').reindex(COND);tr=blended_transform_factory(c.transAxes,c.transData)
 for i,label in enumerate(['MD','MoD','MS','MoS']):t7(c,1.10,yy[i*4:i*4+4].mean(),f'{label}  {sm.iloc[i].median_matched_enrichment:.1f}×',size=7,color=cc[i],tr=tr)
 f.legend(handles=[Line2D([],[],c=LIGHT,lw=3.2,label='Null: 0–99th percentile'),Line2D([],[],c=MID,marker='o',ms=3,ls='',label='Null mean')],loc='lower left',bbox_to_anchor=(.18,.008),ncol=2,fontsize=6.5,handlelength=1.1,columnspacing=1.2)
 STRUCTURE[f]={'canvas_mm':[144,142],'full_xy_axes':0,'real_embeddings':1,'shared_rulers':3,'traditional_panels':3,'hero':'16-context matched-null evidence','height_reduction_vs_v6':'22.8%','biological_color_hues':2}
 SCOPES[f]=None

def f3():
 d=F[2];f=canvas('Figure3_v7_preview',160)
 a=p7(f,[.065,.53,.46,.415],'a','Held-out root states')
 xy=pd.read_csv(str(PROJECT/'source_data'/'figure3'/'plot_inputs'/'heatshock_All_reducedDimS.tsv'),sep='\t',index_col=0)
 md=read('external_heatshock/GSE121619/GSE121619_Control_Heatshock_pData.tsv.gz').set_index('cell');em=xy.join(md)
 src=em.louvain_component.isin([5,6,8]);sealed=em.louvain_component.isin([2,3,10])&em.treatment_id.eq('HeatShock');assert len(em)==2085 and sealed.sum()==830
 a.scatter(em.iloc[:,0],em.iloc[:,1],s=2.2,c=LIGHT,lw=0,rasterized=True)
 a.scatter(em.loc[src].iloc[:,0],em.loc[src].iloc[:,1],s=5,c=BRICK,lw=0,rasterized=True)
 a.scatter(em.loc[sealed].iloc[:,0],em.loc[sealed].iloc[:,1],s=7,fc='white',ec=SECONDARY,lw=.3,rasterized=True)
 for num in [5,6,8,2,3,10]:
  z=em[em.louvain_component.eq(num)].iloc[:,:2].median();a.text(z.iloc[0],z.iloc[1],str(num),fontsize=6.5,ha='center',va='center',bbox=dict(fc='white',ec='none',pad=.7))
 noframe(a);a.set(xlabel='Deposited UMAP 1',ylabel='Deposited UMAP 2')
 t7(a,.99,-.055,'Sources 5/6/8; targets 2/3/10\n830 HS cells sealed',ha='right',va='top',size=6.5,color=SECONDARY)
 b=p7(f,[.745,.755,.215,.19],'b','Candidate recovery and null',pad=18)
 p=wide(d,'heatshock_matched_null');rho=[.615,.621,.450,.415]
 for i,(_,r) in enumerate(p.iterrows()):
  b.plot([0,r.null_q99],[i,i],c=LIGHT,lw=3.5);b.scatter(r.null_mean,i,c=MID,s=9);b.scatter(r.observed_overlap,i,c=BRICK,s=22,zorder=3);b.text(105,i,f'ρ {rho[i]:.3f}',fontsize=6.5,va='center',color=SECONDARY)
 b.set_yticks(range(4),['Cluster 2','Cluster 3','Cluster 10','Seedling']);b.set_ylim(3.5,-.5);ruler(b,(0,149),[0,50,100],'Top100 overlap')
 t7(b,0,1.025,'4/4 BH q ≈ 1e-4',size=6.5)
 c=p7(f,[.745,.515,.215,.125],'c','Shared baselines',pad=18)
 mods=['signed_mean','rank1_common_axis','signed_median','strongest_source_pretarget','pooled_donor'];labs=['Signed mean','Rank-1 axis','State median','Strongest source','Pooled donors']
 for i,m in enumerate(mods):
  v=sel(B,dataset='GSE121619_heatshock',model=m).spearman_full_vector;c.scatter(v,np.full(len(v),i),s=9,c=MID);c.plot([np.median(v)]*2,[i-.20,i+.20],c=BRICK,lw=1.2,zorder=4)
 c.set_yticks(range(5),labs);c.set_ylim(4.5,-.5);ruler(c,(.30,.72),[.3,.5,.7],'Full-response ρ')
 t7(c,0,1.035,'Tick: target median',size=6.5,color=SECONDARY)
 dd=p7(f,[.095,.115,.47,.255],'d','Independent heat-shock transfer')
 z=pd.read_csv(V5/'crossstudy_plot_values.tsv.gz',sep='\t');assert len(z)==16335 and z.top100_overlap.sum()==64
 dd.scatter(z.shared_prediction,z.target_response,s=1.8,c=LIGHT,alpha=.6,lw=0,rasterized=True)
 hit=z[z.top100_overlap];dd.scatter(hit.shared_prediction,hit.target_response,s=8,c=BRICK,lw=0,rasterized=True)
 dd.axhline(0,c=LIGHT,ls='--',lw=.55,zorder=0)
 dd.set(xlabel='Frozen source heat-response score',ylabel='Seedling response\n(Δ log2[CPM + 1])')
 t7(dd,.065,.90,'ρ = 0.415\nTop100 overlap = 64\nOverlap direction = 100%',size=6.5,va='top')
 e=p7(f,[.76,.115,.20,.255],'e','Perturbation specificity')
 p=wide(d,'cross_experiment_wrong_perturbation_controls')
 for i,(_,r) in enumerate(p.iterrows()):
  e.plot([r.wrong_program_rho,r.true_heat_rho],[i,i],c=MID,lw=.9);e.scatter(r.wrong_program_rho,i,c=MID,s=13);e.scatter(r.true_heat_rho,i,c=BRICK,s=21,zorder=3)
 e.set_yticks(range(7),['Salt 6 h','Salt 3 h','Pep1','SCOOP10','SCREW1','flg22','Chitin']);e.set_ylim(6.5,-.5);ruler(e,(-.1,.75),[0,.25,.5,.75],'Paired-universe ρ')
 STRUCTURE[f]={'canvas_mm':[180,160],'full_xy_axes':1,'real_embeddings':1,'shared_rulers':3,'traditional_panels':5,'heroes':['Deposited root UMAP','16335-gene independent transfer'],'all_genes':16335,'highlighted_top100_overlap':64,'sealed_cells':830}
 SCOPES[f]=None

def export7():
 audit={}
 for name,f in figs:
  f.canvas.draw();scope=SCOPES.get(f)
  if scope:
   aa,rows,cols,reason=scope;ids={ax:str(i) for i,ax in enumerate(aa)}
   require_matplotlib_panel_alignment(f,axes=aa,panel_ids=ids,row_groups=[[ids[ax] for ax in g] for g in rows],column_groups=[[ids[ax] for ax in g] for g in cols],json_out=QA/(name+'.alignment.json'),strict=True)
  else:
   (QA/(name+'.alignment.json')).write_text(json.dumps({'status':'intentional_unequal_spans','reason':'The context/design object, dense shared-ruler field and embedded support have distinct sizes; no equal-span plots are claimed.','rectangles':[list(ax.get_position().bounds) for ax in f.axes]},indent=2),encoding='utf-8')
  f.savefig(ROOT/(name+'.pdf'));f.savefig(ROOT/(name+'.svg'));f.savefig(ROOT/(name+'.png'),dpi=300);audit[name]=STRUCTURE[f];print('EXPORTED',name,flush=True)
 if len(figs)==4:
  with PdfPages(ROOT/'Figures1-4_v7_preview.pdf') as pdf:
   for _,f in sorted(figs):pdf.savefig(f)
 (QA/'visual_structure_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
 for path,digest in INPUTS.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
 (QA/'source_input_sha256.json').write_text(json.dumps(INPUTS,ensure_ascii=False,indent=2),encoding='utf-8')

if __name__=='__main__':
 f2();f4()
 if os.environ.get('V7_STAGE')!='layout':f1();f3()
 export7()
