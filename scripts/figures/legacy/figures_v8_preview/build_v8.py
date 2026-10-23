"""v8 low-resolution PNG previews; internal PDFs are QA artifacts only.

F1 layout: hero = objectively selected 18×16 candidate membership matrix;
secondary = true rice map; companion = stability rulers; inset = matched null.
F2 layout: hero = 17×5 context/evidence matrix; secondary = temporal eligibility;
insets = salt baselines, sulfur; complete baseline/q audits remain outside main.
F3 layout: heroes = cell holdout map and 16335-gene transfer scatter;
gene identities = objective 18-gene two-response map and 6 scatter labels;
insets = recovery/null and specificity; full baseline suite remains in ED.
F4 layout: hero = real-data routing composite, centered on 5×4 source metrics;
secondary = 20 target observations; insets = original-unit CIs; boundary audits
are one caption line. Resolution classes are conceptual, not fitted thresholds.
"""
from pathlib import Path
V7=Path(__file__).resolve().parent.parent/'figures_v7_preview'
exec((V7/'build_v7.py').read_text(encoding='utf-8').split('def f2():')[0])
from matplotlib.colors import ListedColormap,LinearSegmentedColormap,Normalize,TwoSlopeNorm
from matplotlib.patches import Rectangle,Patch
from matplotlib.text import Text
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],
 'svg.fonttype':'none','pdf.fonttype':42,'font.size':6.5,'axes.titlesize':7.2,'axes.labelsize':6.5,
 'xtick.labelsize':6,'ytick.labelsize':6,'text.color':'#2F343B'})
NAVY='#466B86';OCHRE='#B18445';BRICK='#AD664F';DUST='#91AABA';DOCHRE='#8E6837'
INK='#2F343B';MID='#A9B1B9';LIGHT='#D8DDE2'
seq=LinearSegmentedColormap.from_list('response_rho',['#F5F6F7','#8AA5B8',NAVY])
div=LinearSegmentedColormap.from_list('difference',[BRICK,'#FFFFFF',NAVY])
signed=LinearSegmentedColormap.from_list('signed_response',['#68727D','#FFFFFF',BRICK])
grayseq=LinearSegmentedColormap.from_list('hits',['#F5F6F7','#5F6A76'])
figs=[];STRUCTURE={};INSETS={}
zheat=pd.read_csv(V5/'crossstudy_plot_values.tsv.gz',sep='\t')
chosen=pd.read_csv(ROOT/'crossstudy_gene_display_selection.tsv',sep='\t')
assert len(zheat)==16335 and zheat.top100_overlap.sum()==64 and len(chosen)==18
heatxy=pd.read_csv(str(PROJECT/'source_data'/'figure3'/'plot_inputs'/'heatshock_All_reducedDimS.tsv'),sep='\t',index_col=0)
heatmd=read('external_heatshock/GSE121619/GSE121619_Control_Heatshock_pData.tsv.gz').set_index('cell');em=heatxy.join(heatmd)
src=em.louvain_component.isin([5,6,8]);sealed=em.louvain_component.isin([2,3,10])&em.treatment_id.eq('HeatShock')

def p8(f,r,l,title,off=False,pad=8):
 a=p7(f,r,l,title,off,pad);a.title.set_fontsize(7.2);return a

def smallbar(f,rect,cmap,norm,label,ticks):
 ca=f.add_axes(rect);mpl.colorbar.ColorbarBase(ca,cmap=cmap,norm=norm,orientation='horizontal',ticks=ticks)
 ca.tick_params(labelsize=5.5,length=1.8,width=.5);ca.set_title(label,fontsize=5.8,pad=3);return ca

def rootmap(a,col=BRICK,labels=True,size=2):
 a.scatter(em.iloc[:,0],em.iloc[:,1],c=LIGHT,s=size,lw=0,rasterized=True)
 a.scatter(em.loc[src].iloc[:,0],em.loc[src].iloc[:,1],c=col,s=size*2,lw=0,rasterized=True)
 a.scatter(em.loc[sealed].iloc[:,0],em.loc[sealed].iloc[:,1],fc='white',ec=SECONDARY,s=size*2.7,lw=.25,rasterized=True)
 if labels:
  for n in [5,6,8,2,3,10]:
   q=em[em.louvain_component.eq(n)].iloc[:,:2].median();a.text(q.iloc[0],q.iloc[1],str(n),fontsize=6,ha='center',va='center',bbox=dict(fc='white',ec='none',pad=.5))
 noframe(a)

def genecloud(a,tiny=False):
 a.scatter(zheat.shared_prediction,zheat.target_response,s=.7 if tiny else 1.5,c=LIGHT,alpha=.6,lw=0,rasterized=True)
 h=zheat[zheat.top100_overlap];a.scatter(h.shared_prediction,h.target_response,s=2 if tiny else 7,c=BRICK,lw=0,rasterized=True)
 if tiny:noframe(a)

def figure1():
 f=canvas('Figure1_v8_preview',114);d=F[0];cc=[OCHRE,DOCHRE,DUST,NAVY]
 a=p8(f,[.04,.51,.20,.40],'a','Rice nuclei')
 rice=pd.read_csv(V5/'rice_deposited_umap.tsv.gz',sep='\t')
 grays={'Epidermis':'#C6CFD6','Meristem':'#98A6B0','Mesophyll':'#E0E5E9','Vasculature':'#B4C0C8'}
 for tissue,g in rice.groupby('TissueSystem',sort=False):a.scatter(g.UMAP1,g.UMAP2,c=grays[tissue],s=.15,lw=0,rasterized=True)
 for tissue,r in rice.groupby('TissueSystem')[['UMAP1','UMAP2']].median().iterrows():a.annotate(tissue,(r.UMAP1,r.UMAP2),xytext=(5 if tissue=='Meristem' else 0,0),textcoords='offset points',ha='center',fontsize=5.8)
 noframe(a);a.set(xlabel='UMAP 1',ylabel='UMAP 2')
 # 18 genes chosen by the documented recurrence/rank rule, never by name.
 vals=pd.read_csv(ROOT/'rice_recurrence_plot_values.tsv',sep='\t');selgenes=pd.read_csv(ROOT/'rice_gene_display_selection.tsv',sep='\t').gene.tolist()
 contexts=[co+'|'+ti for co in COND for ti in ['Epidermis','Meristem','Mesophyll','Vasculature']]
 mat=vals.pivot(index='gene',columns='context',values='membership_count').reindex(index=selgenes,columns=contexts)
 b=p8(f,[.40,.55,.55,.35],'b','Recurrent Top100 candidates',pad=25)
 b.imshow(mat,aspect='auto',interpolation='nearest',cmap=ListedColormap(['#FFFFFF','#CCD3D9',NAVY]),vmin=0,vmax=2)
 b.set_yticks(range(18),selgenes,fontsize=5.8);b.set_xticks(range(16),['E','M','P','V']*4,fontsize=5.8)
 b.tick_params(length=0,pad=2)
 for s in b.spines.values():s.set_visible(False)
 for i,lab in enumerate(['Mild drought','Moderate drought','Mild salinity','Moderate salinity']):
  b.add_patch(Rectangle((i*4-.5,-1.15),4,.32,fc=cc[i],ec='none',clip_on=False))
  b.text(i*4+1.5,-1.5,lab,fontsize=6,ha='center',va='bottom',color=cc[i],clip_on=False)
 b.set_ylim(17.5,-.5)
 f.legend(handles=[Patch(fc='white',ec=MID,lw=.5,label='Neither'),Patch(fc='#CCD3D9',label='One replicate'),Patch(fc=NAVY,label='Both replicates')],loc='lower left',bbox_to_anchor=(.40,.475),ncol=3,fontsize=5.8,handlelength=.9,columnspacing=1)
 c=p8(f,[.045,.12,.22,.285],'c','Two-level stability',True)
 t7(c,.49,.92,'ρ',size=6.2,ha='center');t7(c,.88,.92,'Top100',size=6.2,ha='center')
 for i,(ct,label) in enumerate(zip(COND,['MD','MoD','MS','MoS'])):
  yy=.75-i*.18;t7(c,0,yy,label,size=6.2,color=cc[i])
  for sec,metric,x0,x1,lo,hi in [('rice_full_vector','replicate_response_spearman',.36,.64,-.05,.30),('rice_top100_random_expectation','observed_overlap',.76,1.0,0,36)]:
   v=sel(d,section=sec,metric=metric);v=v[v.context.str.startswith(ct)].value.to_numpy();xx=x0+(v-lo)/(hi-lo)*(x1-x0)
   c.scatter(xx,yy+np.linspace(-.02,.02,4),s=11,c=cc[i],zorder=3);xm=x0+(np.median(v)-lo)/(hi-lo)*(x1-x0);c.plot([xm,xm],[yy-.04,yy+.04],c=INK,lw=.8)
 for lo,hi,x0,x1,ticks in [(-.05,.30,.36,.64,[0,.3]),(0,36,.76,1,[0,30])]:
  c.plot([x0,x1],[.12,.12],c=INK,lw=.6)
  for t in ticks:
   xx=x0+(t-lo)/(hi-lo)*(x1-x0);c.plot([xx,xx],[.12,.10],c=INK,lw=.5);t7(c,xx,.025,str(t),size=5.8,ha='center')
 dd=p8(f,[.40,.09,.40,.31],'d','Matched empirical null',pad=16)
 p=wide(d,'rice_top100_matched_null').reindex(contexts);yy=np.arange(16)+np.repeat(np.arange(4)*.35,4)
 dd.hlines(yy,0,p.null_q99,color=LIGHT,lw=2.8);dd.scatter(p.null_mean,yy,c=MID,s=7,zorder=2);dd.scatter(p.observed_overlap,yy,c=np.repeat(cc,4),s=12,zorder=3)
 dd.set_yticks(yy,[x.split('|')[1] for x in contexts],fontsize=5.8);dd.set_ylim(yy[-1]+.5,-.5);ruler(dd,(0,35),[0,10,20,30],'Top100 overlap');dd.tick_params(axis='x',labelsize=5.8)
 t7(dd,0,1.025,'16/16 BH q ≈ 1e-4',size=5.8)
 sm=wide(d,'rice_condition_summary').reindex(COND);tr=blended_transform_factory(dd.transAxes,dd.transData)
 for i,label in enumerate(['MD','MoD','MS','MoS']):t7(dd,1.12,yy[i*4:i*4+4].mean(),f'{label} {sm.iloc[i].median_matched_enrichment:.1f}×',size=6.2,color=cc[i],tr=tr)
 STRUCTURE[f]={'canvas_mm':[180,114],'hero':'18 genes × 16 true recurrence contexts','all_nuclei':126139,'selection_file':'rice_gene_display_selection.tsv','membership':'0 neither, 1 single replicate, 2 both; all original intersections reproduced exactly','secondary':'Full-vector and candidate stability; matched null','removed_to_ED':'Complete depth, full controls and detailed recurrence tables'}

def route_values():
 d=F[1];rows=[];names=[];ys=[];starts=[]
 specs=[('lateral_root_8h_to_20h','frozen_reveal_primary_results','full_rho','top100_overlap',['same_state_latest','shared_other_states'],0),('GSE46205_salt_3h_corrected','post_reveal_preprocessing_corrected_results','spearman_full_vector','top100_overlap',['same_state_history','shared_other3'],10.5),('GSE46205_salt_8h_corrected','post_reveal_preprocessing_corrected_results','spearman_full_vector','top100_overlap',['same_state_history','shared_other3'],16)]
 rename={'Ambiguous Stele Cells':'Ambiguous stele','Columella/ Root Cap':'Columella / root cap','Epidermis':'Epidermis','Lateral Root Endodermis':'Lateral-root endodermis','Mature Pericycle':'Mature pericycle','Phloem':'Phloem','Phloem Pole Pericycle':'Phloem-pole pericycle','Xylem':'Xylem','Xylem Pole Pericycle':'Xylem-pole pericycle'}
 for ds,sec,rho,hits,mods,y0 in specs:
  p=sel(d,dataset=ds,section=sec,metric=rho).pivot_table(index='context',columns='model',values='value',sort=False)
  h=sel(d,dataset=ds,section=sec,metric=hits).pivot_table(index='context',columns='model',values='value',sort=False)
  starts.append(y0)
  for i,(ct,r) in enumerate(p.iterrows()):
   rows.append([r[mods[0]],r[mods[1]],r[mods[0]]-r[mods[1]],h.loc[ct,mods[0]],h.loc[ct,mods[1]]]);names.append(rename.get(ct,ct.capitalize()));ys.append(y0+i)
 pd.DataFrame(rows,index=names,columns=['history_rho','shared_rho','delta_rho','history_top100','shared_top100']).to_csv(ROOT/'route_heatmap_plot_values.tsv',sep='\t')
 return np.array(rows),names,np.array(ys)

def figure2():
 f=canvas('Figure2_v8_preview',112)
 a=p8(f,[.05,.53,.18,.37],'a','Eligible evidence',True)
 for i,label in enumerate(['A','B','C']):t7(a,-.02,.72-i*.25,label,size=6,ha='right')
 for x,lab in [(.08,'0 h'),(.45,'Early'),(.93,'Target')]:t7(a,x,.96,lab,size=6,ha='center')
 for x,y,target in [(.08,.72,False),(.45,.72,False),(.93,.72,True),(.45,.47,False),(.45,.22,False)]:a.add_patch(Ellipse((x,y),.045,.06,fc='white' if target else NAVY,ec=INK if target else NAVY,lw=.7,ls='--' if target else '-'))
 arr(a,(.11,.72),(.41,.72),LIGHT,.6);arr(a,(.5,.72),(.88,.72),NAVY,1)
 a.plot([.49,.88],[.47,.68],c=MID,lw=.5);a.plot([.49,.88],[.22,.68],c=MID,lw=.5)
 t7(a,0,.02,'Own history exists',size=6.3,color=NAVY)
 values,names,yy=route_values();b=p8(f,[.47,.20,.28,.67],'b','',pad=31);f.text(.47,.965,'Context × evidence',fontsize=7.2)
 norms=[Normalize(-.1,1),Normalize(-.1,1),TwoSlopeNorm(vmin=-.5,vcenter=0,vmax=.5),Normalize(0,100),Normalize(0,100)];maps=[seq,seq,div,grayseq,grayseq]
 for i,row in enumerate(values):
  for j,val in enumerate(row):
   color=maps[j](norms[j](val));b.add_patch(Rectangle((j-.43,yy[i]-.43),.86,.86,fc=color,ec='none'))
   lum=.2126*color[0]+.7152*color[1]+.0722*color[2]
   b.text(j,yy[i],f'{val:.2f}' if j<3 else f'{val:.0f}',fontsize=5.8,ha='center',va='center',color='white' if lum<.55 else INK)
 b.set(xlim=(-.5,4.5),ylim=(19.5,-1.0));b.set_yticks(yy,names,fontsize=6)
 b.set_xticks(range(5),['History ρ','Shared ρ','Δρ','History hits','Shared hits'],rotation=50,ha='left',rotation_mode='anchor');b.xaxis.tick_top();b.tick_params(length=0,pad=3)
 for s in b.spines.values():s.set_visible(False)
 tr=blended_transform_factory(b.transAxes,b.transData)
 for y,lab in [(-.7,'Lateral root | 8→20 h'),(9.65,'Root layers | 1→3 h'),(15.15,'Root layers | 3→8 h')]:t7(b,-.85,y,lab,size=6.2,tr=tr,color=NAVY)
 # One compact real baseline inset; the rest of the baseline suite stays in ED4.
 c=p8(f,[.845,.60,.125,.27],'c','Salt baselines')
 for i,ds in enumerate(['GSE153103_salt_6h','GSE153103_salt_12h']):
  g=sel(B,dataset=ds);ms=[next(m for m in g.model if m.startswith(pre)) for pre in ['latest','history_mean','linear_last2']];ms.append(next(m for m in g.model if m.startswith('linear_') and not m.startswith('linear_last2')))
  c.plot(range(4),g.set_index('model').loc[ms].spearman_full_vector,'o-',c=[NAVY,MID][i],ms=2.4,lw=.75,label=['6 h','12 h'][i])
 c.set_xticks(range(4),['L','M','T','B']);c.set(ylim=(-.1,1),ylabel='Full-response ρ');c.tick_params(labelsize=5.8);zero(c)
 c.legend(loc='lower left',bbox_to_anchor=(-.05,-.47),ncol=2,fontsize=5.8,handlelength=.8,columnspacing=.5)
 ins=f.add_axes([.845,.30,.125,.06]);ins.set(xlim=(0,30),ylim=(0,1));ins.axis('off');p=wide(sel(F[1],dataset='sulfur_3h_to_12h'),'matched_empirical_null')
 ins.plot([0,float(p.null_q99.iloc[0])],[.5,.5],c=LIGHT,lw=3.2);ins.scatter(float(p.null_mean.iloc[0]),.5,c=MID,s=8);ins.scatter(21,.5,c=NAVY,s=18,zorder=3)
 t7(ins,0,1.6,'Sulfur 3→12 h',size=6.3)
 t7(ins,0,-.5,'ρ = 0.555\nTop100 = 21\nq ≈ 2e-4',size=6.2,color=SECONDARY,va='top')
 smallbar(f,[.47,.095,.075,.018],seq,Normalize(-.1,1),'Full ρ',[-.1,1]);smallbar(f,[.575,.095,.075,.018],div,Normalize(-.5,.5),'Δρ',[-.5,0,.5]);smallbar(f,[.68,.095,.075,.018],grayseq,Normalize(0,100),'Top100',[0,100])
 f.text(.045,.30,'History wins\n9/9 lateral root\nexact P = 0.00195\n4/4 + 4/4 root layers',fontsize=6.3,linespacing=1.5,color=SECONDARY)
 f.text(.29,.025,'GSE46205: corrected_post_reveal. Matched null: 13/16 significant; lateral root 6/9.',fontsize=6.1,color=SECONDARY)
 STRUCTURE[f]={'canvas_mm':[180,112],'hero':'17 contexts × 5 evidence measurements','values':'rho columns share scale; delta uses a centered scale; hits share a 0–100 scale','history_wins':17,'provenance':'Both root-layer blocks are corrected_post_reveal','insets':['Salt baselines','Sulfur null support'],'removed_to_ED':'Full baseline and absolute-persistence audit'}

def figure3():
 f=canvas('Figure3_v8_preview',124)
 a=p8(f,[.045,.56,.25,.36],'a','Held-out root states');rootmap(a,size=1.6)
 a.set(xlabel='Deposited UMAP 1',ylabel='Deposited UMAP 2')
 f.text(.045,.505,'Sources 5/6/8 → targets 2/3/10; 830 HS cells sealed',fontsize=6,color=SECONDARY)
 b=p8(f,[.40,.65,.215,.25],'b','Candidate recovery / null',pad=16);p=wide(F[2],'heatshock_matched_null')
 for i,(_,r) in enumerate(p.iterrows()):
  b.plot([0,r.null_q99],[i,i],c=LIGHT,lw=3.1);b.scatter(r.null_mean,i,c=MID,s=8);b.scatter(r.observed_overlap,i,c=BRICK,s=20,zorder=3);b.text(105,i,f'ρ {[.615,.621,.450,.415][i]:.3f}',fontsize=5.8,va='center')
 b.set_yticks(range(4),['Cluster 2','Cluster 3','Cluster 10','Seedling'],fontsize=6);b.set_ylim(3.5,-.5);ruler(b,(0,149),[0,50,100],'Top100 overlap');t7(b,0,1.03,'4/4 BH q ≈ 1e-4',size=5.8)
 # Objectively selected shared genes; separate response scales preserve source/target units.
 c=p8(f,[.815,.40,.135,.47],'c','',pad=31);f.text(.815,.96,'Shared genes',fontsize=7.2)
 for i,(_,r) in enumerate(chosen.iterrows()):
  for j,key,maximum in [(0,'shared_prediction',8),(1,'target_response',13)]:
   val=r[key];color=signed(Normalize(-maximum,maximum)(val));c.add_patch(Rectangle((j-.43,i-.43),.86,.86,fc=color,ec='none'))
  agree=np.sign(r.shared_prediction)==np.sign(r.target_response);c.text(2,i,'=' if agree else '≠',fontsize=6,ha='center',va='center',color=INK)
 c.set(xlim=(-.5,2.5),ylim=(17.5,-.5));c.set_yticks(range(18),chosen.gene,fontsize=6);c.set_xticks([0,1,2],['Root','Seedling','Sign'],rotation=50,ha='left',rotation_mode='anchor');c.xaxis.tick_top();c.tick_params(length=0,pad=2)
 for sp in c.spines.values():sp.set_visible(False)
 f.text(.815,.365,'Root: ±8; seedling: ±13\nSeparate response scales',fontsize=5.8,color=SECONDARY,linespacing=1.4)
 dd=p8(f,[.085,.13,.455,.315],'d','Independent heat-shock transfer')
 genecloud(dd);dd.axhline(0,c=LIGHT,ls='--',lw=.5,zorder=0)
 dd.set(xlabel='Frozen source heat-response score',ylabel='Seedling response\n(Δ log2[CPM + 1])')
 t7(dd,.055,.90,'ρ = 0.415\nTop100 overlap = 64\nOverlap direction = 100%',size=6.2,va='top')
 # Six labels selected by the same joint-rank rule, positioned in response order.
 dd.set_xlim(dd.get_xlim());dd.set_ylim(dd.get_ylim())
 labs=chosen.head(6).sort_values('target_response',ascending=False)
 for (_,r),y in zip(labs.iterrows(),np.linspace(.96,.55,6)):
  end=dd.transData.inverted().transform(dd.transAxes.transform((1.025,y)))
  dd.plot([r.shared_prediction,end[0]],[r.target_response,end[1]],c=MID,lw=.45,clip_on=False)
  t7(dd,1.04,y,r.gene,size=5.8,color=BRICK)
 e=p8(f,[.805,.115,.165,.205],'e','Specificity')
 p=wide(F[2],'cross_experiment_wrong_perturbation_controls')
 for i,(_,r) in enumerate(p.iterrows()):e.plot([r.wrong_program_rho,r.true_heat_rho],[i,i],c=MID,lw=.7);e.scatter(r.wrong_program_rho,i,c=MID,s=11);e.scatter(r.true_heat_rho,i,c=BRICK,s=17,zorder=3)
 e.set_yticks(range(7),['Salt 6 h','Salt 3 h','Pep1','SCOOP10','SCREW1','flg22','Chitin'],fontsize=6);e.set_ylim(6.5,-.5);ruler(e,(-.1,.75),[0,.5],'Paired-universe ρ')
 STRUCTURE[f]={'canvas_mm':[180,124],'heroes':['Real 2085-cell UMAP','16335-gene transfer'],'shared_gene_display':18,'scatter_labels':6,'selection':'Shared Top100 intersection, ascending sum of source and target absolute ranks; deterministic ties','response_heatmap':'Column-specific symmetric original-unit scales ±8 and ±13; sign agreement is measured','removed_to_ED':'Full baseline and ranking-depth audit'}

def figure4():
 f=canvas('Figure4_v8_preview',112)
 f.text(.035,.96,'a',fontsize=9,weight='bold');f.text(.065,.96,'Real-data evidence routing',fontsize=7.2)
 # Three true data thumbnails, rather than an invented 3×3 expression glyph.
 tp=pd.read_csv(ROOT/'temporal_thumbnail_plot_values.tsv',sep='\t')
 ta=f.add_axes([.095,.735,.13,.155]);tv=tp[['r1','r3','r6']].to_numpy();tm=float(np.ceil(np.abs(tv).max()));ta.imshow(tv,aspect='auto',cmap=div,vmin=-tm,vmax=tm)
 for i in range(3):
  for j in range(3):ta.text(j,i,f'{tv[i,j]:.1f}',ha='center',va='center',fontsize=5.5,color='white' if tv[i,j]/tm>.65 else INK)
 ta.set_yticks(range(3),tp.gene,fontsize=5.8);ta.set_xticks(range(3),['1 h','3 h','6 h'],fontsize=5.8);ta.tick_params(length=0,pad=2)
 for s in ta.spines.values():s.set_visible(False)
 f.text(.065,.915,'Future time',fontsize=6.5);f.text(.045,.69,'Recent same-state history',fontsize=6.2,color=NAVY)
 ua=f.add_axes([.055,.415,.17,.22]);rootmap(ua,col=NAVY,labels=False,size=.35)
 f.text(.065,.655,'Unseen state',fontsize=6.5);f.text(.045,.38,'Shared cross-state response',fontsize=6.2,color=SECONDARY)
 sa=f.add_axes([.07,.135,.155,.175]);genecloud(sa,tiny=True)
 f.text(.065,.33,'Unseen experiment',fontsize=6.5);f.text(.045,.105,'Frozen cross-study program',fontsize=6.2,color=SECONDARY)
 # Source capability is itself the measured 5×4 source-only evidence object.
 S=wide(F[3],'stimulus_level_source_vs_target');order=['Pep1','SCOOP10','SCREW1','flg22','Chitin'];high=order[:3]
 metrics=['rank1_energy','median_pairwise_rho','effective_state_fraction','shared_rms'];raw=S.reindex(order)[metrics];scaled=(raw-raw.min())/(raw.max()-raw.min())
 cap=f.add_axes([.365,.37,.18,.35]);cap.imshow(scaled,aspect='auto',vmin=0,vmax=1,cmap=seq)
 cap.set_yticks(range(5),order,fontsize=6);cap.set_xticks(range(4),['Rank-1','Pairwise ρ','Eff. fraction','RMS'],rotation=55,ha='left',rotation_mode='anchor',fontsize=6);cap.xaxis.tick_top();cap.tick_params(length=0,pad=3)
 for sp in cap.spines.values():sp.set_visible(False)
 for i in range(5):
  for j in range(4):cap.text(j,i,f'{raw.iloc[i,j]:.3f}' if j==3 else f'{raw.iloc[i,j]:.2f}',fontsize=5.8,ha='center',va='center',color='white' if scaled.iloc[i,j]>.65 else INK)
 f.text(.345,.845,'Source capability',fontsize=7.2)
 smallbar(f,[.365,.305,.18,.016],seq,Normalize(0,1),'Column-scaled; raw values shown',[0,1])
 # Minimal routes connect genuine objects; they are not learned weights.
 route=f.add_axes([0,0,1,1]);route.set_label('annotation routes');route.set(xlim=(0,1),ylim=(0,1));route.axis('off');route.patch.set_visible(False)
 curve7(route,(.235,.80),(.29,.81),(.30,.755),(.36,.745),NAVY,.75)
 curve7(route,(.235,.52),(.29,.53),(.28,.70),(.36,.745),MID,.65)
 curve7(route,(.235,.225),(.30,.22),(.30,.29),(.36,.295),MID,.65)
 # Actual frozen candidate list illustrates a gene-level output, separately from
 # the immune calibration. Resolution classes have no fitted cutoffs.
 f.text(.585,.96,'b',fontsize=9,weight='bold');f.text(.61,.96,'Output resolution',fontsize=7.2)
 f.text(.60,.89,'Gene ranking\n+ direction',fontsize=6.5,color=NAVY,linespacing=1.2)
 f.text(.60,.81,'Frozen heat example',fontsize=5.8,color=SECONDARY)
 for i,r in zheat.sort_values(['shared_prediction_abs_rank','gene']).head(6).iterrows():
  rank=int(r.shared_prediction_abs_rank);f.text(.60,.77-(rank-1)*.032,f'{rank}  {r.gene}  {"+" if r.shared_prediction>0 else "−"}',fontsize=6.1,color=NAVY)
 f.text(.60,.44,'Coarse response\nprogram',fontsize=6.5,color=SECONDARY,linespacing=1.3)
 f.text(.60,.22,'No gene-level call\n(abstain)',fontsize=6.5,color=SECONDARY,linespacing=1.3)
 arr(route,(.55,.64),(.59,.70),NAVY,.8);arr(route,(.55,.49),(.59,.46),MID,.65);arr(route,(.55,.39),(.59,.24),MID,.5)
 f.text(.32,.16,'Conceptual resolution; no fitted thresholds.',fontsize=6,color=SECONDARY)
 # The only primary quantitative axes: nested target performance.
 T=wide(sel(F[3],section='immune_frozen_reveal',context=['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']),idx=['dataset','context'])
 c=p8(f,[.795,.22,.175,.61],'c','Nested targets',pad=20)
 for i,s in enumerate(order):
  v=T.loc[s].full_rho;c.scatter(i+np.linspace(-.13,.13,4),v,c=MID,s=12);c.scatter(i,v.mean(),c=NAVY if s in high else BRICK,s=25,zorder=3)
 c.set_xticks(range(5),order,rotation=90,ha='right',rotation_mode='anchor');c.set(ylim=(-.1,1.15),ylabel='Target full-response ρ');c.set_yticks([0,.5,1]);zero(c)
 f.text(.795,.85,'5 stimuli × 4 targets',fontsize=6)
 H=wide(F[3],'hierarchical_stimulus_inference');ci_axes=[]
 for rect,key,lim,ticks,label in [([.805,.745,.06,.055],'full_rho',(-.05,.45),[0,.4],'Δρ'),([.905,.745,.06,.055],'top100_overlap',(-5,50),[0,40],'ΔTop100')]:
  aa=f.add_axes(rect);aa.set_label('CI inset '+key);ci_axes.append(aa);r=H.loc[key];aa.errorbar(r.difference_allow_minus_reject,0,xerr=[[r.difference_allow_minus_reject-r.bootstrap_ci025],[r.bootstrap_ci975-r.difference_allow_minus_reject]],fmt='o',c=NAVY,ms=2.5,capsize=1.8,lw=.65)
  aa.set_ylim(-1,1);aa.set_yticks([]);ruler(aa,lim,ticks);aa.tick_params(labelsize=5.5,length=1.5);aa.axvline(0,c=LIGHT,ls='--',lw=.5);t7(aa,.5,1.22,label,size=5.8,ha='center')
 INSETS[f]=ci_axes
 f.text(.81,.69,'95% CI; exact P floor 0.10',fontsize=5.8,color=SECONDARY)
 f.text(.045,.035,'Geometry and amplitude covary; 19/20 targets exceed source-matched nulls.',fontsize=6.1,color=SECONDARY)
 STRUCTURE[f]={'canvas_mm':[180,112],'hero':'Three real data thumbnails → actual 5-stimulus × 4-metric source capability object','source_metric_colors':'Within-column min-max display scaling; every raw metric printed','temporal_input':'3 genes × r1/r3/r6, objectively selected from pre-12h frozen source artifact','output_example':'Actual top 6 frozen heat genes; illustrative output, not a fit of immune stimulus thresholds','primary_quantitative_axes':1,'CI_insets':'Original units, existing hierarchical CIs','boundary':'Geometry/amplitude confounded; no independent geometry mechanism claimed.'}

def axes_overlap_check(f):
 findings=[]
 for i,a in enumerate(f.axes):
  ra=a.get_position()
  for j,b in enumerate(f.axes[i+1:],i+1):
   rb=b.get_position();w=min(ra.x1,rb.x1)-max(ra.x0,rb.x0);h=min(ra.y1,rb.y1)-max(ra.y0,rb.y0)
   if w*f.get_figwidth()*72>1 and h*f.get_figheight()*72>1:
    allow=('annotation routes' in [a.get_label(),b.get_label()] or a.get_label().startswith('CI inset') or b.get_label().startswith('CI inset'))
    findings.append({'axes':[i,j],'intentional':allow,'reason':'Transparent annotation canvas / CI inside target panel' if allow else 'Unexpected plot rectangle overlap'})
 return findings

def tightbbox_check(f):
 f.canvas.draw();ren=f.canvas.get_renderer();fb=f.bbox;bad=[]
 for obj in f.findobj(match=Text):
  if not obj.get_visible() or not obj.get_text().strip():continue
  try:b=obj.get_window_extent(ren)
  except Exception:continue
  if b.x0<fb.x0-1 or b.y0<fb.y0-1 or b.x1>fb.x1+1 or b.y1>fb.y1+1:bad.append({'text':obj.get_text(),'bbox':list(b.bounds)})
 return bad

if __name__=='__main__':
 figure1();figure2();figure3();figure4();report={}
 for name,f in figs:
  if os.environ.get('V8_ONLY') and name not in os.environ['V8_ONLY'].split(','):continue
  out=tightbbox_check(f)
  if f in INSETS:
   aa=INSETS[f];require_matplotlib_panel_alignment(f,axes=aa,panel_ids={a:str(i) for i,a in enumerate(aa)},row_groups=[['0','1']],column_groups=[],json_out=QA/(name+'.alignment.json'),strict=True)
  else:
   (QA/(name+'.alignment.json')).write_text(json.dumps({'status':'intentional_unequal_spans','reason':'Independent biological/data objects, shared-ruler companions and colorbars have different roles and unequal spans; no equal-size panel claim.','rectangles':[list(a.get_position().bounds) for a in f.axes]},indent=2),encoding='utf-8')
  overlaps=axes_overlap_check(f)
  report[name]={'outside_canvas':out,'axes_overlap_checks':overlaps,'design':STRUCTURE[f],'axes_rectangles':[list(a.get_position().bounds) for a in f.axes]}
  f.savefig(ROOT/(name+'.png'),dpi=220)
  # Internal vector PDF only for mandatory text/collision QA, not final production export.
  f.savefig(QA/(name+'.pdf'))
  print('PREVIEW',name,'outside_canvas=',len(out),flush=True)
 (QA/'renderer_checks.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
 for path,digest in INPUTS.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
