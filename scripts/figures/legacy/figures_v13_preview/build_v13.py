"""v13 architecture pass. Figure 2/3 reused byte-for-byte from v11; no analysis changes.
F1 claim: strong candidate recurrence is supported despite limited vector stability.
F4 claim: task structure determines eligibility; source capability informs resolution.
VISUAL_PATH: tissue/context -> stability -> same-context candidates -> matched null;
             task anchors -> same-stimulus source/target rows -> inference boundary.
No new inference, model, gene selection, decorative connector, Word or ED change.
"""
from pathlib import Path
V9=Path(__file__).resolve().parent.parent/'figures_v9_preview'
exec((V9/'build_v9.py').read_text(encoding='utf-8').split('def figure2():')[0])
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],'svg.fonttype':'none','pdf.fonttype':42,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7})
figs=[];STRUCTURE={};ALIGN={}
TISSUES=['Epidermis','Meristem','Mesophyll','Vasculature']
TCOL=[NAVY,DUST,DBRICK,BRICK]
for p in [V8/'rice_recurrence_plot_values.tsv',V8/'rice_gene_display_selection.tsv',V8/'temporal_thumbnail_plot_values.tsv',V5/'rice_deposited_umap.tsv.gz',V5/'crossstudy_plot_values.tsv.gz',Path(str(PROJECT/'source_data'/'figure3'/'plot_inputs'/'heatshock_All_reducedDimS.tsv'))]:
 INPUTS[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()

def title(a,x,y,letter,label):
 t=tt(a,x-.03,y,letter,size=9.5);t.set_weight('bold');tt(a,x,y,label,size=8.2)

def clean(a):
 for s in a.spines.values():s.set_visible(False)
 a.tick_params(length=0,pad=2)

def figure1():
 f=canvas('Figure1_v13_preview',100);a=scene(f);title(a,.05,.967,'a','Rice stress landscape');title(a,.43,.967,'b','Recurrent Top100 candidates')
 rice=pd.read_csv(V5/'rice_deposited_umap.tsv.gz',sep='\t');assert len(rice)==126139
 um=f.add_axes([.025,.61,.27,.27]);um.set_label('rice tissue map')
 grays=['#C6CFD6','#98A6B0','#E0E5E9','#B4C0C8']
 for ti,col in zip(TISSUES,grays):
  g=rice[rice.TissueSystem.eq(ti)];um.scatter(g.UMAP1,g.UMAP2,s=.2,c=col,lw=0,rasterized=True)
 for ti,r in rice.groupby('TissueSystem')[['UMAP1','UMAP2']].median().iterrows():
  um.annotate(ti,(r.UMAP1,r.UMAP2),xytext=(7 if ti=='Meristem' else 0,0),textcoords='offset points',ha='center',va='center',fontsize=6.5)
 noframe(um)
 tt(a,.025,.566,'4 tissue systems × 4 stress conditions',size=6.5)
 tt(a,.025,.535,'16 matched contexts',color=SECONDARY)
 # Same treatment rows; each point is the same tissue in both measurements.
 tt(a,.035,.382,'Replicate stability',size=7)
 tt(a,.225,.340,'Full-vector ρ',size=6.5,ha='center');tt(a,.332,.340,'Top100 hits',size=6.5,ha='center')
 stable=[]
 for i,(co,lab,col) in enumerate(zip(COND,['Mild drought','Moderate drought','Mild salinity','Moderate salinity'],TCOL)):
  y=.293-i*.043;tt(a,.035,y,lab,size=7,color=col)
  for sec,metric,x0,w,lo,hi in [('rice_full_vector','replicate_response_spearman',.195,.07,-.05,.30),('rice_top100_random_expectation','observed_overlap',.295,.07,0,36)]:
   s=sel(F[0],section=sec,metric=metric).set_index('context');vals=np.array([s.loc[co+'|'+ti,'value'] for ti in TISSUES]);xs=x0+(vals-lo)/(hi-lo)*w
   a.scatter(xs,y+np.linspace(-.008,.008,4),c=col,s=10,zorder=3);xm=x0+(np.median(vals)-lo)/(hi-lo)*w;a.plot([xm,xm],[y-.017,y+.017],c=INK,lw=.85)
   for ti,v in zip(TISSUES,vals):stable.append(dict(treatment=co,tissue=ti,metric=metric,value=v,row_y=y))
 for x0,w,lo,hi,ticks in [(.195,.07,-.05,.3,[0,.3]),(.295,.07,0,36,[0,30])]:
  a.plot([x0,x0+w],[.112,.112],c=SECONDARY,lw=.6)
  for v in ticks:
   x=x0+(v-lo)/(hi-lo)*w;a.plot([x,x],[.112,.105],c=SECONDARY,lw=.5);tt(a,x,.080,str(v),ha='center')
 tt(a,.025,.041,'Dots: tissues; tick: treatment median',color=SECONDARY)
 tt(a,.43,.040,'E: epidermis · M: meristem · P: mesophyll · V: vasculature',color=SECONDARY)
 # Fixed v8 gene list. Matrix and null have EXACTLY the same x centers and limits.
 vals=pd.read_csv(V8/'rice_recurrence_plot_values.tsv',sep='\t');genes=pd.read_csv(V8/'rice_gene_display_selection.tsv',sep='\t').gene.tolist()
 contexts=[co+'|'+ti for co in COND for ti in TISSUES]
 mat=vals.pivot(index='gene',columns='context',values='membership_count').reindex(index=genes,columns=contexts);assert mat.shape==(18,16) and not mat.isna().any().any()
 hm=f.add_axes([.43,.36,.54,.54]);hm.set_label('candidate context columns');hm.imshow(mat,aspect='auto',cmap=ListedColormap(['white',PALE,'#78838D']),vmin=0,vmax=2,interpolation='nearest');hm.set_yticks(range(18),genes,fontsize=6.5);hm.set_xticks(range(16),['E','M','P','V']*4,fontsize=7);clean(hm)
 for i,(lab,col) in enumerate(zip(['Mild drought','Moderate drought','Mild salinity','Moderate salinity'],TCOL)):
  xx=.43+.54*(i+.5)/4;tt(a,xx,.928,lab,ha='center',size=6.5,color=col);a.add_patch(Rectangle((.43+i*.135,.910),.135,.006,fc=col,ec='none'))
 for x,c,lab in [(.43,'white','Neither'),(.58,PALE,'One replicate'),(.77,'#78838D','Both replicates')]:
  a.add_patch(Rectangle((x,.297),.010,.012,fc=c,ec=MID,lw=.5));tt(a,x+.017,.303,lab)
 tt(a,.43,.270,'Matched null · 16/16 BH q ≈ 1e-4',size=7)
 nu=f.add_axes([.43,.130,.54,.085]);nu.set_label('matched null same context columns');p=wide(F[0],'rice_top100_matched_null').reindex(contexts)
 nu.vlines(range(16),0,p.null_q99,colors=PALE,lw=4);nu.scatter(range(16),p.null_mean,c=MID,s=8,zorder=3);nu.scatter(range(16),p.observed_overlap,c=np.repeat(TCOL,4),s=16,zorder=4)
 nu.set(xlim=(-.5,15.5),ylim=(0,36));nu.set_xticks([]);nu.set_yticks([0,30]);nu.tick_params(length=2,width=.5,labelsize=7);nu.spines['top'].set_visible(False);nu.spines['right'].set_visible(False);nu.spines['bottom'].set_visible(False)
 for k in [3.5,7.5,11.5]:
  hm.axvline(k,c='white',lw=.55);nu.axvline(k,c=PALE,lw=.55,zorder=0)
 hm.set_xlim(-.5,15.5)
 tt(a,.43,.105,'Top100 overlap · treatment median null enrichment below (×)',color=SECONDARY)
 sm=wide(F[0],'rice_condition_summary').reindex(COND)
 for i,col in enumerate(TCOL):tt(a,.43+.54*(i+.5)/4,.080,f'{sm.iloc[i].median_matched_enrichment:.1f}×',ha='center',color=col)
 tt(a,.43,.015,'Matched null controls baseline expression and instability',color=SECONDARY)
 tt(a,.43,.241,'Grey: null mean / 99th percentile; coloured: observed',color=SECONDARY)
 mat.to_csv(ROOT/'figure1_candidate_matrix.tsv',sep='\t');p.to_csv(ROOT/'figure1_null_values.tsv',sep='\t');pd.DataFrame(stable).to_csv(ROOT/'figure1_stability_values.tsv',sep='\t',index=False)
 pd.DataFrame({'column':range(16),'context':contexts,'matrix_x':range(16),'null_x':range(16)}).to_csv(ROOT/'figure1_context_mapping.tsv',sep='\t',index=False)
 ALIGN[f]=([hm,nu],[],[['matrix','null']],['matrix','null']);STRUCTURE[f]={'canvas_mm':[180,100],'contexts':16,'genes':18,'rice_nuclei':len(rice),'visual_path':'Rice tissues → same-treatment stability → candidate recurrence → same-column matched null','identity':'Exact matrix/null context column order; same treatment/tissue order for stability','letters':['a','b']}

def figure4():
 f=canvas('Figure4_v13_preview',108);a=scene(f);title(a,.05,.974,'a','Task anchors');title(a,.39,.974,'b','Source capability × transfer')
 # Three existing true-data thumbnails; proximity conveys task role, no arrows.
 tt(a,.035,.934,'Future time',size=7)
 tp=pd.read_csv(V8/'temporal_thumbnail_plot_values.tsv',sep='\t');tv=tp[['r1','r3','r6']].to_numpy();lim=float(np.ceil(np.abs(tv).max()))
 ta=f.add_axes([.105,.825,.21,.070]);ta.set_label('frozen temporal data');ta.imshow(tv,aspect='auto',cmap=div,vmin=-lim,vmax=lim)
 for i in range(3):
  for j in range(3):ta.text(j,i,f'{tv[i,j]:.1f}',ha='center',va='center',fontsize=6.5,color='white' if tv[i,j]/lim>.65 else INK)
 ta.set_yticks(range(3),tp.gene,fontsize=6.5);ta.set_xticks(range(3),['1 h','3 h','6 h'],fontsize=6.5);clean(ta)
 tt(a,.08,.777,'Recent same-state history',color=NAVY)
 tt(a,.035,.718,'Unseen state',size=7)
 um=f.add_axes([.025,.435,.25,.245]);um.set_label('held-out root thumbnail');rootmap(um,col=BRICK,labels=False,size=.45)
 tt(a,.035,.400,'Shared cross-state response',color=SECONDARY)
 tt(a,.055,.345,'Unseen experiment',size=7)
 sc=f.add_axes([.055,.185,.28,.13]);sc.set_label('cross-study gene thumbnail');genecloud(sc,tiny=True)
 tt(a,.055,.150,'Frozen cross-study program',color=SECONDARY)
 # A single stimulus row carries source-only raw metrics and all four targets.
 order=['Pep1','SCOOP10','SCREW1','flg22','Chitin'];ys=[.78,.685,.59,.455,.36]
 S=wide(F[3],'stimulus_level_source_vs_target').reindex(order);metrics=['rank1_energy','median_pairwise_rho','effective_state_fraction','shared_rms'];raw=S[metrics];scaled=(raw-raw.min())/(raw.max()-raw.min())
 T=wide(sel(F[3],section='immune_frozen_reveal',context=['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']),idx=['dataset','context'])
 for lo,hi,c in [(.545,.815,NAVY),(.315,.49,BRICK)]:a.plot([.325,.325],[lo,hi],c=c,lw=1.1,alpha=.35)
 capcmap=LinearSegmentedColormap.from_list('neutral capability',[PALE,'#78838D'])
 for x,lab in zip([.435,.505,.575,.645],['Rank-1\nenergy','Pairwise\nρ','Effective\nfraction','Shared\nRMS']):tt(a,x,.857,lab,size=7,ha='center')
 tt(a,.795,.857,'Target ρ',size=7,ha='center')
 tt(a,.39,.909,'Grey: four targets; coloured: stimulus mean',color=SECONDARY)
 exported=[]
 for i,(s,y) in enumerate(zip(order,ys)):
  col=NAVY if i<3 else BRICK;tt(a,.396,y,s,ha='right',size=7)
  for j,x in enumerate([.435,.505,.575,.645]):
   c=capcmap(scaled.iloc[i,j]);a.add_patch(Rectangle((x-.031,y-.028),.062,.056,fc=c,ec='none'));tt(a,x,y,f'{raw.iloc[i,j]:.3f}' if j==3 else f'{raw.iloc[i,j]:.2f}',ha='center',color='white' if scaled.iloc[i,j]>.65 else INK)
  v=T.loc[s].full_rho;off=np.linspace(-.012,.012,4);xs=.715+.16*v.to_numpy();a.scatter(xs,y+off,c=MID,s=13,zorder=3);a.scatter(.715+.16*v.mean(),y,c=col,s=26,zorder=4)
  for (context,val),oy in zip(v.items(),off):exported.append(dict(stimulus=s,context=context,target_rho=val,source_row_y=y,target_row_y=y,point_y=y+oy,stimulus_mean=v.mean()))
 tt(a,.75,.938,'Higher reliability',color=NAVY);tt(a,.75,.900,'Gene candidates + direction',color=SECONDARY)
 tt(a,.72,.278,'Lower reliability',color=BRICK);tt(a,.72,.251,'Reduce output resolution',color=SECONDARY)
 for val in [0,.5,1]:
  x=.715+.16*val;a.plot([x,x],[.810,.805],c=SECONDARY,lw=.5);tt(a,x,.829,str(val),ha='center')
 tt(a,.40,.296,'Column-scaled colour; raw values shown',color=SECONDARY)
 a.plot([.715,.875],[.810,.810],c=SECONDARY,lw=.5)
 tt(a,.86,.278,'No fitted thresholds',color=SECONDARY)
 # Original hierarchical inference in compact original-unit rulers.
 title(a,.39,.250,'c','Hierarchical contrast')
 H=wide(F[3],'hierarchical_stimulus_inference');ci=[]
 for rect,key,lims,ticks,label in [([.40,.115,.18,.035],'full_rho',(-.05,.45),[0,.4],'Δ full-response ρ'),([.72,.160,.21,.035],'top100_overlap',(-5,50),[0,40],'Δ Top100 hits')]:
  ax=f.add_axes(rect);ax.set_label('hierarchical '+key);ci.append(ax);r=H.loc[key]
  ax.errorbar(r.difference_allow_minus_reject,0,xerr=[[r.difference_allow_minus_reject-r.bootstrap_ci025],[r.bootstrap_ci975-r.difference_allow_minus_reject]],fmt='o',c=NAVY,ms=3,capsize=2,lw=.7)
  ax.set_ylim(-1,1);ax.set_yticks([]);ruler(ax,lims,ticks,'');ax.tick_params(labelsize=7,length=2,width=.5);ax.axvline(0,c=MID,lw=.5,ls='--')
  fmt='.3f' if key=='full_rho' else '.1f';ax.set_title(('Δρ ' if key=='full_rho' else 'ΔTop100 ')+format(r.difference_allow_minus_reject,'+'+fmt)+' ['+format(r.bootstrap_ci025,fmt)+', '+format(r.bootstrap_ci975,fmt)+']',fontsize=7,pad=3)
 tt(a,.39,.223,'95% bootstrap CI',color=SECONDARY)
 tt(a,.39,.197,'Highest independent unit: stimulus (n = 5)',color=SECONDARY)
 tt(a,.60,.080,'Exact one-sided stimulus-label permutation P = 0.10 (10 labelings)',color=SECONDARY)
 tt(a,.035,.047,'Geometry and amplitude covary; n = 5 does not identify an independent geometric effect.',color=SECONDARY)
 tt(a,.035,.021,'Matched null: 19/20 targets significant; lower reliability retains transferable information.',color=SECONDARY)
 raw.to_csv(ROOT/'figure4_source_raw_values.tsv',sep='\t');pd.DataFrame(exported).to_csv(ROOT/'figure4_target_row_mapping.tsv',sep='\t',index=False);H.to_csv(ROOT/'figure4_existing_inference.tsv',sep='\t');tp.to_csv(ROOT/'figure4_temporal_thumbnail.tsv',sep='\t',index=False)
 ALIGN[f]=None;STRUCTURE[f]={'canvas_mm':[180,108],'stimuli':5,'nested_targets':20,'highest_independent_unit':'stimulus n=5','visual_path':'Actual task anchors → same-stimulus source/target rows → output resolution → inferential boundary','normalization':'Within-column min-max for descriptive colour only; raw values printed','letters':['a','b','c'],'no_arrows':True}

def check_render(f):
 f.canvas.draw();ren=f.canvas.get_renderer();fb=f.bbox;outside=[]
 for t in f.findobj(match=Text):
  if not t.get_visible() or not t.get_text().strip():continue
  b=t.get_window_extent(ren)
  if b.x0<fb.x0-1 or b.y0<fb.y0-1 or b.x1>fb.x1+1 or b.y1>fb.y1+1:outside.append(t.get_text())
 return {'outside_canvas':outside,'axes_count':len(f.axes),'objects':[{'object':ax.get_label(),'bbox':list(ax.get_position().bounds)} for ax in f.axes]}

if __name__=='__main__':
 # Hash already approved figures to establish that they remain untouched.
 protected={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT.parent/'figures_v11_preview').glob('Figure[23]_v11_preview.*')}
 figure1();figure4();reports={}
 for name,f in figs:
  ch=check_render(f);assert not ch['outside_canvas'],ch;reports[name]={'design':STRUCTURE[f],**ch}
  if ALIGN[f] is not None:
   axes,rows,cols,ids=ALIGN[f];require_matplotlib_panel_alignment(f,axes=axes,panel_ids=dict(zip(axes,ids)),row_groups=rows,column_groups=cols,json_out=QA/(name+'.alignment.json'),strict=True)
  else:
   (QA/(name+'.alignment.json')).write_text(json.dumps({'status':'intentional floating insets','reason':'CI rulers have different original units and deliberately staggered positions/unequal widths; no equal-span comparison group. Stimulus identity rows remain numerically identical.','objects':ch['objects']},indent=2),encoding='utf-8')
  f.savefig(ROOT/(name+'.png'),dpi=600);f.savefig(ROOT/(name+'.pdf'));f.savefig(ROOT/(name+'.svg'))
  print(name,STRUCTURE[f],ch,flush=True)
 (QA/'renderer_checks.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2),encoding='utf-8')
 (ROOT/'input_hashes.json').write_text(json.dumps(INPUTS,ensure_ascii=False,indent=2),encoding='utf-8')
 for p,h in {**INPUTS,**protected}.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
 (QA/'protected_v11_hashes.json').write_text(json.dumps(protected,ensure_ascii=False,indent=2),encoding='utf-8')
