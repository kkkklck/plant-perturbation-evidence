"""v11 final polish, Figure 2/3 only. English PNG/PDF/SVG; no Word changes.

Claim F2: Direct same-state history is the robust default for future-time tasks.
Claim F3: The shared heat program recovers held-out states and transfers studies.
Archetype: context-linked atlas; identity-linked transfer scene.
VISUAL_PATH = actual times -> same context rows -> measured route differences;
              actual clusters -> recovery -> actual genes -> specificity.
No fabricated statistics or decorative cross-object connectors.
"""
from pathlib import Path
V9=Path(__file__).resolve().parent.parent/'figures_v9_preview'
exec((V9/'build_v9.py').read_text(encoding='utf-8').split('def figure2():')[0])
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],'svg.fonttype':'none','pdf.fonttype':42})
figs=[];STRUCTURE={};ANCHORS={};ALIGN={}
for plot_input in [V5/'crossstudy_plot_values.tsv.gz',V8/'crossstudy_gene_display_selection.tsv',Path(str(PROJECT/'source_data'/'figure3'/'plot_inputs'/'heatshock_All_reducedDimS.tsv'))]:
 INPUTS[str(plot_input)]=hashlib.sha256(plot_input.read_bytes()).hexdigest()

def timebone(a,x,w,y,times,edges,targets,arc=.034):
 xx={t:x+w*t/max(times) for t in times}
 a.plot([x,x+w],[y,y],c=MID,lw=.65)
 for t in times:
  a.plot(xx[t],y,'o',ms=2.8,mfc='white' if t in targets or t>max(targets) else MID,mec=NAVY if t in targets else MID,mew=.65)
  tt(a,xx[t],y+.030,str(t),ha='center')
 tt(a,x+w+.015,y+.030,'h')
 for u,v in edges:
  path(a,(xx[u],y+.003),(xx[u]+(xx[v]-xx[u])*.25,y+arc),(xx[v]-(xx[v]-xx[u])*.25,y+arc),(xx[v],y+.003),NAVY,.9)
 return xx

def figure2():
 f=canvas('Figure2_v11_preview',94.5);a=scene(f)
 # Same context groups span the temporal branch and measured evidence.
 for lo,hi in [(.459,.705),(.297,.408),(.135,.246)]:
  a.add_patch(Rectangle((.027,lo),.953,hi-lo,fc=INK,alpha=.045,ec='none',zorder=-3))
 tt(a,.02,.98,'a',size=9.5);tt(a,.05,.98,'Observed experimental times',size=8.2)
 tt(a,.05,.937,'Salt bulk · GSE153103',size=7.2)
 saltpos=timebone(a,.05,.47,.875,[0,1,3,6,12],[(3,6),(6,12)],{6,12},.045)
 sa=f.add_axes([.625,.867,.245,.092]);sa.set_label('salt timeline inset');sa.patch.set_visible(False)
 salt_out=[]
 for i,ds in enumerate(['GSE153103_salt_6h','GSE153103_salt_12h']):
  g=sel(B,dataset=ds);mods=[next(m for m in g.model if m.startswith(pre)) for pre in ['latest','history_mean','linear_last2']];mods.append(next(m for m in g.model if m.startswith('linear_') and not m.startswith('linear_last2')))
  v=g.set_index('model').loc[mods].spearman_full_vector.to_numpy()
  sa.plot(range(4),v,'o-',c=[NAVY,MID][i],ms=2.5,lw=.75)
  t=[6,12][i];tt(a,saltpos[t],.832,f'Latest {v[0]:.2f} / trend {v[2]:.2f}',ha='center',color=NAVY)
  salt_out.append(dict(dataset=ds,latest_rho=v[0],mean_rho=v[1],trend_rho=v[2],baseline_trend_rho=v[3]))
 sa.set(ylim=(0,1),xlim=(-.15,3.15));sa.set_yticks([0,1]);sa.set_xticks(range(4),['Latest','Mean','Trend','+baseline']);sa.tick_params(length=2,width=.6,labelsize=7)
 sa.spines['left'].set_visible(False)
 tt(a,.585,.947,'ρ',color=SECONDARY);tt(a,.765,.945,'6 h',color=NAVY);tt(a,.82,.945,'12 h',color=SECONDARY)
 # Each evidence row has a single cell-type label, shared with its temporal branch.
 vals,names,_,blocks=route_data();ys=np.r_[np.arange(9),np.arange(11,15),np.arange(17,21)];yy=.69-.027*ys
 tt(a,.452,.792,'b',size=9.5);tt(a,.48,.792,'Context × evidence',size=8.2)
 for x,lab in [(.50,'History ρ'),(.582,'Shared ρ'),(.664,'Δρ')]:tt(a,x,.746,lab,size=7,ha='center')
 tt(a,.826,.777,'Top100 hits',size=7.2,ha='center');tt(a,.765,.746,'History',ha='center');tt(a,.89,.746,'Shared',ha='center')
 tt(a,.455,.719,'9/9 full-response wins · exact P = 0.00195',size=7)
 tt(a,.455,.421,'4/4 full-response wins',size=7);tt(a,.455,.260,'4/4 full-response wins',size=7)
 tt(a,.20,.719,'Lateral root · GSE158761',size=7)
 lp=timebone(a,.035,.145,.735,[0,8,20],[(8,20)],{20})
 tt(a,.22,.421,'Root layers · GSE46205',size=7)
 tt(a,.22,.260,'3→8 h · root layers',size=7)
 rp1=timebone(a,.035,.145,.414,[0,1,3,8],[(1,3)],{3})
 rp2=timebone(a,.035,.145,.252,[0,1,3,8],[(3,8)],{8})
 groups=[(range(9),lp,8,20,.735,'GSE158761',20),(range(9,13),rp1,1,3,.414,'GSE46205',3),(range(13,17),rp2,3,8,.252,'GSE46205',8)]
 label_art=[];contexts=[]
 norms=[Normalize(-.1,1),Normalize(-.1,1),TwoSlopeNorm(vmin=-.5,vcenter=0,vmax=.5)]
 for i,(row,label,y) in enumerate(zip(vals,names,yy)):
  label_art.append(tt(a,.435,y,label,ha='right'))
  # The short guide joins the same context label to its measured evidence cells.
  a.plot([.439,.458],[y,y],c=MID,lw=.5,alpha=.32)
  for j,x in enumerate([.50,.582,.664]):
   col=[seq,seq,div][j](norms[j](row[j]));a.add_patch(Rectangle((x-.034,y-.013),.068,.026,fc=col,ec='none'))
   lum=.2126*col[0]+.7152*col[1]+.0722*col[2];tt(a,x,y,f'{row[j]:.2f}',ha='center',color='white' if lum<.55 else INK)
  for x,val,c in [(.744,row[3],NAVY),(.87,row[4],MID)]:
   a.plot([x,x+.083*val/100],[y,y],c=c,lw=3.1,solid_capstyle='butt');tt(a,.855 if x<.8 else .978,y,str(int(val)),ha='right')
 # The vertical forks belong to one experiment/time per group, not separate schematics.
 f.canvas.draw();renderer=f.canvas.get_renderer()
 for idx,xx,st,en,ytime,ds,target in groups:
  ids=list(idx);lo=yy[ids[-1]]
  for t in [st,en]:a.plot([xx[t],xx[t]],[ytime-.008,lo],c=PALE,lw=.55,zorder=0)
  for i in ids:
   y=yy[i];a.plot([xx[st],xx[en]],[y,y],c=NAVY,lw=.65,alpha=.7)
   a.plot(xx[st],y,'o',ms=1.7,c=NAVY);a.plot(xx[en],y,'o',ms=2,mfc='white',mec=NAVY,mew=.5)
   left=f.transFigure.inverted().transform((label_art[i].get_window_extent(renderer).x0,0))[0]
   a.plot([xx[en]+.004,left-.005],[y,y],c=MID,lw=.5,alpha=.32)
   contexts.append(dict(dataset=ds,target_h=target,context=names[i],row_y=y,timeline_row_y=y,evidence_row_y=y,history_h=st))
 # Only shared unit rulers; no independent Top100 axes per group.
 rulers=[]
 for rect,cmap,norm,title,ticks in [([.47,.074,.09,.014],seq,Normalize(-.1,1),'Full ρ',[-.1,1]),([.61,.074,.09,.014],div,Normalize(-.5,.5),'Δρ',[-.5,0,.5])]:
  cb=f.add_axes(rect);cb.set_label('shared color ruler');mpl.colorbar.ColorbarBase(cb,cmap=cmap,norm=norm,orientation='horizontal',ticks=ticks);cb.tick_params(length=1.8,width=.5,labelsize=7);cb.set_title(title,fontsize=7,pad=3);rulers.append(cb)
 for x in [.744,.87]:
  a.plot([x,x+.083],[.082,.082],c=SECONDARY,lw=.6)
  for v in [0,100]:
   pos=x+.083*v/100;a.plot([pos,pos],[.082,.077],c=SECONDARY,lw=.5);tt(a,pos,.055,str(v),ha='center')
 a.add_patch(Rectangle((.03,.059),.397,.04,fc=PALE,alpha=.5,ec='none',zorder=-1))
 tt(a,.039,.079,'Sulfur 3→12 h   ρ 0.555   Top100 21   q ≈ 2e-4',color=SECONDARY)
 tt(a,.035,.014,'GSE46205: corrected_post_reveal · Matched null: 13/16 significant; lateral root 6/9.',color=SECONDARY)
 pd.DataFrame(vals,index=names,columns=['history_rho','shared_rho','delta_rho','history_top100','shared_top100']).to_csv(ROOT/'figure2_context_values.tsv',sep='\t')
 pd.DataFrame(contexts).to_csv(ROOT/'figure2_context_row_mapping.tsv',sep='\t',index=False);pd.DataFrame(salt_out).to_csv(ROOT/'figure2_salt_outcomes.tsv',sep='\t',index=False)
 ALIGN[f]=rulers;STRUCTURE[f]={'canvas_mm':[180,94.5],'height_reduction_percent':10,'contexts':17,'visual_path':'Observed times → identity-shared row labels → history/shared evidence','removed':'Alternative-state fan; repeated row schematics; separate sulfur glyph','scientific_boundary':'Salt 6h/12h are separate frozen tasks; GSE46205 remains corrected_post_reveal.'}

def figure3():
 f=canvas('Figure3_v11_preview',102);um=f.add_axes([.035,.675,.35,.255]);um.set_label('real root UMAP');rootmap(um,col=BRICK,labels=False,size=1.5)
 meds={};ut={}
 for n in [5,6,8,2,3,10]:
  xy=em[em.louvain_component.eq(n)].iloc[:,:2].median().to_numpy();meds[n]=xy
  ut[n]=um.text(*xy,str(n),ha='center',va='center',fontsize=6.5,bbox=dict(fc='white',ec='none',pad=.6))
 a=scene(f);tt(a,.02,.978,'a',size=9.5);tt(a,.052,.978,'Held-out root states',size=8.2)
 tt(a,.05,.948,'Sources 5/6/8 → sealed targets 2/3/10',color=SECONDARY)
 tt(a,.52,.944,'Top100 recovery · source-matched null',size=7.2);tt(a,.52,.913,'4/4 BH q ≈ 1e-4',color=SECONDARY)
 f.canvas.draw();uv={n:f.transFigure.inverted().transform(um.transData.transform(meds[n])) for n in [2,3,10]}
 ordered=sorted(uv,key=lambda n:uv[n][1],reverse=True);rows=[(n,float(uv[n][1])) for n in ordered];rows.append(('seedling',min(y for _,y in rows)-.047))
 p=wide(F[2],'heatshock_matched_null');rho={2:.615,3:.621,10:.450,'seedling':.415};rt={}
 for n,y in rows:
  k='whole_seedling' if n=='seedling' else 'cluster_'+str(n);r=p.loc[k]
  rt[n]=tt(a,.48,y,'Seedling' if n=='seedling' else 'Cluster '+str(n),ha='right')
  x0=.52;w=.23;a.plot([x0,x0+w*r.null_q99/100],[y,y],c=PALE,lw=3.4);a.plot(x0+w*r.null_mean/100,y,'o',ms=2.4,c=MID)
  ox=x0+w*r.observed_overlap/100;a.plot(ox,y,'D' if n=='seedling' else 'o',ms=4,c=BRICK);tt(a,ox+.01,y,str(int(r.observed_overlap)),color=BRICK);tt(a,.79,y,f'ρ = {rho[n]:.3f}',color=SECONDARY)
 ry=.887;a.plot([.52,.75],[ry,ry],c=SECONDARY,lw=.6)
 for v in [0,50,100]:
  x=.52+.23*v/100;a.plot([x,x],[ry,ry-.006],c=SECONDARY,lw=.5);tt(a,x,ry-.024,str(v),ha='center')
 # Recovery leaders link label identity, never an arbitrary space or a numeric dot.
 f.canvas.draw();ren=f.canvas.get_renderer();cluster_map=[]
 for n in ordered:
  b1=ut[n].get_window_extent(ren);b2=rt[n].get_window_extent(ren);y=uv[n][1]
  x1=f.transFigure.inverted().transform((b1.x1+2,0))[0];x2=f.transFigure.inverted().transform((b2.x0-2,0))[0]
  f.add_artist(Line2D([x1,x2],[y,y],transform=f.transFigure,c=MID,lw=.55,alpha=.3,zorder=15))
  cluster_map.append(dict(cluster=n,centroid_x=meds[n][0],centroid_y=meds[n][1],label_figure_y=y,recovery_figure_y=y,guide_start_x=x1,guide_end_x=x2))
 pd.DataFrame(cluster_map).to_csv(ROOT/'figure3_cluster_mapping.tsv',sep='\t',index=False)
 # Five baseline methods share the same three target identities; individual points + median tick.
 ba=f.add_axes([.165,.547,.19,.112]);ba.set_label('embedded shared baseline ruler');ba.patch.set_visible(False)
 mods=['signed_mean','rank1_common_axis','signed_median','strongest_source_pretarget','pooled_donor'];labs=['Signed mean','Rank-1','Median','Strongest','Pooled']
 baseline_rows=[]
 for i,(m,lab) in enumerate(zip(mods,labs)):
  g=sel(B,dataset='GSE121619_heatshock',model=m);v=g.spearman_full_vector.to_numpy();assert len(v)==3
  ba.scatter(v,np.full(len(v),i),s=7,c=MID,zorder=2);ba.plot([np.median(v)]*2,[i-.25,i+.25],c=BRICK,lw=1.1,zorder=4)
  baseline_rows.extend(g.assign(display_method=lab).to_dict('records'))
 ba.set_yticks(range(5),labs);ba.set_ylim(4.25,-.25);ruler(ba,(-.1,.8),[0,.5],'');ba.tick_params(labelsize=6.5,length=1.8,width=.5)
 pd.DataFrame(baseline_rows).to_csv(ROOT/'figure3_baseline_values.tsv',sep='\t',index=False)
 tt(a,.405,.573,'Grey: null mean / 99th percentile; rust: observed',color=SECONDARY)
 # Merge the frozen source file with the established, separately derived target-response table.
 pred=read('external_heatshock_crossstudy/GSE290176/crossstudy_blind/frozen_crossstudy_prediction.tsv.gz')
 targets=zheat[['gene','target_response','target_abs_rank']].copy()
 universe=pred[['gene','shared_prediction','shared_prediction_abs_rank']].merge(targets,on='gene',how='inner',validate='one_to_one')
 q=universe.set_index('gene').reindex(zheat.gene);assert len(universe)==16335 and np.isfinite(q[['shared_prediction','target_response']]).all().all()
 np.testing.assert_allclose(q.shared_prediction,zheat.shared_prediction,rtol=0,atol=1e-12);np.testing.assert_allclose(q.target_response,zheat.target_response,rtol=0,atol=1e-12)
 universe['top100_overlap']=(universe.shared_prediction_abs_rank<=100)&(universe.target_abs_rank<=100);assert universe.top100_overlap.sum()==64
 universe['joint_rank']=universe.shared_prediction_abs_rank+universe.target_abs_rank
 six=universe[universe.top100_overlap].sort_values(['joint_rank','shared_prediction_abs_rank','gene']).head(6).sort_values('target_response',ascending=False)
 assert six.gene.nunique()==6 and set(six.gene)<=set(universe.gene)
 sc=f.add_axes([.075,.12,.49,.33]);sc.set_label('all-gene independent transfer');genecloud(sc);sc.patch.set_visible(False)
 sc.axhline(0,c=PALE,lw=.5,ls='--',zorder=0);sc.set(xlabel='Frozen root heat-response score',ylabel='Independent seedling response\n(Δ log2[CPM + 1])');sc.tick_params(length=2.4,width=.6,labelsize=7);sc.xaxis.label.set_size(7);sc.yaxis.label.set_size(7)
 sc.set_xlim(sc.get_xlim());sc.set_ylim(sc.get_ylim())
 tt(a,.02,.484,'b',size=9.5);tt(a,.052,.484,'Independent heat-shock transfer',size=8.2)
 tt(a,.095,.458,'ρ = 0.415 · Top100 overlap = 64 · overlap direction = 100%')
 tt(a,.700,.453,'Root',ha='center');tt(a,.750,.453,'Seedling',ha='center')
 label_rows=[];f.canvas.draw()
 for (_,r),y in zip(six.iterrows(),np.linspace(.427,.277,6)):
  x,y0=float(r.shared_prediction),float(r.target_response);loc=f.transFigure.inverted().transform(sc.transData.transform((x,y0)))
  np.testing.assert_allclose(sc.transData.inverted().transform(f.transFigure.transform(loc)),[x,y0],atol=1e-10)
  tt(a,.59,y,r.gene,color=BRICK)
  # Figure-level artists are ABOVE axes, so the entire leader remains visible inside the cloud.
  end=(.584,y);elbow=(.575,y)
  f.add_artist(Line2D([loc[0],elbow[0],end[0]],[loc[1],elbow[1],end[1]],transform=f.transFigure,c=MID,lw=.55,alpha=.65,zorder=20,solid_capstyle='butt'))
  sc.plot(x,y0,'o',ms=3,mfc='white',mec=DBRICK,mew=.65,zorder=8)
  for bx,val,maxval,c in [(.685,x,8,DBRICK),(.735,y0,13,BRICK)]:
   assert 0<=val<=maxval;a.plot([bx,bx+.031],[y,y],c=PALE,lw=.55);a.plot([bx,bx+.031*val/maxval],[y,y],c=c,lw=2.7,solid_capstyle='butt')
  label_rows.append(dict(gene=r.gene,source_score=x,target_response=y0,label_x=.59,label_y=y,point_figure_x=loc[0],point_figure_y=loc[1],source_abs_rank=int(r.shared_prediction_abs_rank),target_abs_rank=int(r.target_abs_rank),joint_rank=int(r.joint_rank)))
 pd.DataFrame(label_rows).to_csv(ROOT/'figure3_labelled_genes.tsv',sep='\t',index=False)
 tt(a,.59,.225,'Bar scales (original units)\nRoot 0–8; seedling 0–13',color=SECONDARY)
 sp=f.add_axes([.875,.14,.105,.29]);sp.set_label('attached perturbation specificity');sp.patch.set_visible(False)
 p=wide(F[2],'cross_experiment_wrong_perturbation_controls')
 for i,(_,r) in enumerate(p.iterrows()):
  sp.plot([r.wrong_program_rho,r.true_heat_rho],[i,i],c=MID,lw=.65);sp.scatter(r.wrong_program_rho,i,c=MID,s=10);sp.scatter(r.true_heat_rho,i,c=BRICK,s=14,zorder=3)
 sp.set_yticks(range(7),['Salt 6 h','Salt 3 h','Pep1','SCOOP10','SCREW1','flg22','Chitin']);sp.set_ylim(6.5,-.5);ruler(sp,(-.1,.75),[0,.5],'ρ');sp.tick_params(labelsize=6.5,length=1.7,width=.5)
 tt(a,.80,.478,'Specificity controls',size=7)
 STRUCTURE[f]={'canvas_mm':[180,102],'height_reduction_percent':100*8/110,'cells':2085,'sealed_cells':830,'genes':16335,'shared_top100':64,'labelled_genes':6,'identity_leaders':{'cluster':3,'gene':6},'gene_selection':'Top six by joint absolute rank; display descending target response','baseline_targets':3,'visual_path':'Cluster identity → recovery; frozen gene ID → independent gene response → specificity','removed':'14-gene side list, 8 unlabelled gene bars, hidden leader segments','display_boundary':'Deposited UMAP is post-reveal display, not prediction input; 100% direction applies only to 64 shared Top100 genes.'}

def check_render(f):
 f.canvas.draw();ren=f.canvas.get_renderer();fb=f.bbox;outside=[]
 for t in f.findobj(match=Text):
  if not t.get_visible() or not t.get_text().strip():continue
  b=t.get_window_extent(ren)
  if b.x0<fb.x0-1 or b.y0<fb.y0-1 or b.x1>fb.x1+1 or b.y1>fb.y1+1:outside.append({'text':t.get_text(),'bbox':list(b.bounds)})
 boxes=[{'index':i,'object':a.get_label(),'bbox':list(a.get_position().bounds)} for i,a in enumerate(f.axes)];over=[]
 for i,aa in enumerate(f.axes):
  for j,bb in enumerate(f.axes[i+1:],i+1):
   r,s=aa.get_position(),bb.get_position();w=min(r.x1,s.x1)-max(r.x0,s.x0);h=min(r.y1,s.y1)-max(r.y0,s.y0)
   if w*f.get_figwidth()*72>1 and h*f.get_figheight()*72>1:over.append({'axes':[i,j],'intentional':'linked annotation / data field' in [aa.get_label(),bb.get_label()]})
 return {'outside_canvas':outside,'axes_count':len(f.axes),'objects':boxes,'axes_overlaps':over}

if __name__=='__main__':
 figure2();figure3();reports={}
 for name,f in figs:
  ch=check_render(f);reports[name]={'design':STRUCTURE[f],**ch};assert not ch['outside_canvas'],ch['outside_canvas']
  assert all(o['intentional'] for o in ch['axes_overlaps']),ch['axes_overlaps']
  if f in ALIGN:
   axes=ALIGN[f];require_matplotlib_panel_alignment(f,axes=axes,panel_ids={ax:str(i) for i,ax in enumerate(axes)},row_groups=[['0','1']],column_groups=[],json_out=QA/(name+'.ruler-alignment.json'),strict=True)
  # Heterogeneous UMAP, scatter and embedded rulers are intentionally not equal panels.
  (QA/(name+'.alignment-exemptions.json')).write_text(json.dumps({'status':'intentional heterogeneous objects','reason':'Different data roles and physical spans; shared numeric rows checked in exported context and identity mappings.','objects':ch['objects']},indent=2),encoding='utf-8')
  f.savefig(ROOT/(name+'.png'),dpi=600);f.savefig(ROOT/(name+'.pdf'));f.savefig(ROOT/(name+'.svg'))
  print(name,'canvas_mm',STRUCTURE[f]['canvas_mm'],'axes/insets',len(f.axes),'clipping',ch['outside_canvas'],'overlaps',ch['axes_overlaps'],flush=True)
  for o in ch['objects']:print('OBJECT',o,flush=True)
 (QA/'renderer_checks.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2),encoding='utf-8')
 (ROOT/'input_hashes.json').write_text(json.dumps(INPUTS,ensure_ascii=False,indent=2),encoding='utf-8')
 for p,h in INPUTS.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
