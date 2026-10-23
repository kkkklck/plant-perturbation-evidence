"""Stage 1 only: Figure 2/3 linked previews. Figure 1/4 await user review.

VISUAL_PATH = {
 'Figure2': ['real experimental times','aligned context evidence','rho differences / candidate counts','embedded salt and sulfur support'],
 'Figure3': ['real source/target cells','cluster-linked recovery','independent gene response','linked gene bars and specificity'],
}
No new fitted statistics, embeddings, source tables, manuscript or Word changes.
"""
from pathlib import Path
V8=Path(__file__).resolve().parent.parent/'figures_v8_preview'
prefix=(V8/'build_v8.py').read_text(encoding='utf-8').split('def figure1():')[0]
prefix=prefix.replace("chosen=pd.read_csv(ROOT/'crossstudy_gene_display_selection.tsv'","chosen=pd.read_csv(V8/'crossstudy_gene_display_selection.tsv'")
exec(prefix)
NAVY='#4C718C';DUST='#88A9BE';BRICK='#B7674C';DBRICK='#98523B';INK='#30343B';MID='#C8CDD3';PALE='#E8EBEE'
seq=LinearSegmentedColormap.from_list('rho_v9',['#F5F6F7',DUST,NAVY])
div=LinearSegmentedColormap.from_list('difference_v9',[BRICK,'white',NAVY])
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],
 'font.size':6.5,'axes.labelsize':6.5,'axes.titlesize':7.5,'xtick.labelsize':6.5,'ytick.labelsize':6.5,
 'svg.fonttype':'none','pdf.fonttype':42,'text.color':INK,'axes.edgecolor':INK,'axes.labelcolor':INK})
figs=[];STRUCTURE={};ROLES={}
VISUAL_PATH={'Figure2':['real experimental times','aligned context evidence','rho differences / candidate counts','embedded salt and sulfur support'],'Figure3':['real source/target cells','cluster-linked recovery','independent gene response','linked gene bars and specificity']}

def scene(f):
 a=f.add_axes([0,0,1,1]);a.set(xlim=(0,1),ylim=(0,1));a.axis('off');a.patch.set_visible(False);a.set_label('linked annotation / data field');return a

def tt(a,x,y,s,size=6.5,color=INK,ha='left',va='center'):
 return a.text(x,y,s,fontsize=size,color=color,ha=ha,va=va,fontweight='bold' if s in ['a','b'] else 'normal')

def path(a,st,c1,c2,en,col=NAVY,lw=.65,alpha=1):
 p=PathPatch(MPath([st,c1,c2,en],[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4]),fc='none',ec=col,lw=lw,alpha=alpha);a.add_patch(p);return p

def timeline(a,start,width,y,times,edges,targets,label):
 xx={t:start+width*t/max(times) for t in times}
 tt(a,start,y+.065,label,size=7.2)
 a.plot([start,start+width],[y,y],c=PALE,lw=.65)
 for t in times:
  a.plot(xx[t],y,marker='o',ms=3.5,mec=NAVY if t in targets else MID,mfc='white' if t in targets else MID,mew=.7)
  tt(a,xx[t],y+.027,str(t),ha='center')
 tt(a,start+width+.022,y+.027,'h')
 for t0,t1 in edges:
  x0,x1=xx[t0],xx[t1];path(a,(x0,y+.008),(x0+(x1-x0)*.27,y+.065),(x0+(x1-x0)*.73,y+.065),(x1,y+.008),NAVY,.85)
 return xx

def route_data():
 rows=[];names=[];ys=[];blocks=[]
 specs=[('lateral_root_8h_to_20h','frozen_reveal_primary_results','full_rho',['same_state_latest','shared_other_states'],0),('GSE46205_salt_3h_corrected','post_reveal_preprocessing_corrected_results','spearman_full_vector',['same_state_history','shared_other3'],10.5),('GSE46205_salt_8h_corrected','post_reveal_preprocessing_corrected_results','spearman_full_vector',['same_state_history','shared_other3'],16)]
 rename={'Ambiguous Stele Cells':'Ambiguous stele','Columella/ Root Cap':'Columella / root cap','Epidermis':'Epidermis','Lateral Root Endodermis':'Lateral-root endodermis','Mature Pericycle':'Mature pericycle','Phloem':'Phloem','Phloem Pole Pericycle':'Phloem-pole pericycle','Xylem':'Xylem','Xylem Pole Pericycle':'Xylem-pole pericycle'}
 for ds,sec,metric,mods,y0 in specs:
  p=sel(F[1],dataset=ds,section=sec,metric=metric).pivot_table(index='context',columns='model',values='value',sort=False)
  h=sel(F[1],dataset=ds,section=sec,metric='top100_overlap').pivot_table(index='context',columns='model',values='value',sort=False)
  blocks.append((ds,y0,len(p)))
  for i,(ct,r) in enumerate(p.iterrows()):
   rows.append([r[mods[0]],r[mods[1]],r[mods[0]]-r[mods[1]],h.loc[ct,mods[0]],h.loc[ct,mods[1]]]);names.append(rename.get(ct,ct.capitalize()));ys.append(y0+i)
 return np.array(rows),names,np.array(ys),blocks

def figure2():
 f=canvas('Figure2_v9_preview',120);a=scene(f)
 tt(a,.025,.977,'a',size=9.5);tt(a,.055,.977,'Observed experimental times',size=7.5)
 # Different actual systems and proportional time spacing; 6h/12h are separate targets.
 timeline(a,.335,.39,.887,[0,1,3,6,12],[(3,6),(6,12)],{6,12},'Salt bulk · GSE153103')
 salt=f.add_axes([.815,.83,.17,.095]);salt.set_label('embedded salt baseline');salt.patch.set_visible(False)
 for i,ds in enumerate(['GSE153103_salt_6h','GSE153103_salt_12h']):
  g=sel(B,dataset=ds);ms=[next(m for m in g.model if m.startswith(pre)) for pre in ['latest','history_mean','linear_last2']];ms.append(next(m for m in g.model if m.startswith('linear_') and not m.startswith('linear_last2')))
  salt.plot(range(4),g.set_index('model').loc[ms].spearman_full_vector,'o-',ms=2.8,lw=.75,c=[NAVY,MID][i],label=['6 h','12 h'][i])
 salt.set_xticks(range(4),['Latest','Mean','Trend','+base']);salt.set(ylim=(0,1),xlim=(-.15,3.15));salt.set_yticks([0,1]);salt.tick_params(length=2,width=.6)
 salt.spines['left'].set_visible(False);salt.set_title('ρ',loc='left',pad=2,fontsize=6.5)
 salt.legend(ncol=2,loc='upper center',bbox_to_anchor=(.5,1.34),fontsize=6.5,handlelength=.7,columnspacing=.7,frameon=False)
 # Sulfur attaches to the real temporal-design scaffold rather than a separate panel.
 tt(a,.035,.885,'Sulfur 3→12 h',size=7.2)
 p=wide(sel(F[1],dataset='sulfur_3h_to_12h'),'matched_empirical_null')
 sx=.035;sw=.19;sy=.848
 a.plot([sx,sx+sw*float(p.null_q99.iloc[0])/30],[sy,sy],c=PALE,lw=3.5);a.plot(sx+sw*float(p.null_mean.iloc[0])/30,sy,'o',ms=2.5,c=MID);a.plot(sx+sw*21/30,sy,'o',ms=4,c=NAVY)
 tt(a,.035,.811,'ρ = 0.555 · Top100 = 21',color=SECONDARY)
 tt(a,.035,.785,'Matched-null q ≈ 2e-4',color=SECONDARY)
 # All 17 contexts share one row scaffold; no repeated context labels.
 vals,names,ys,blocks=route_data();global_y=.646-.0255*ys
 tt(a,.495,.773,'b',size=9.5);tt(a,.53,.773,'Context × evidence',size=7.5)
 for x,lab in [(.535,'History ρ'),(.604,'Shared ρ'),(.673,'Δρ')]:tt(a,x,.726,lab,ha='center')
 tt(a,.825,.754,'Top100 recurrence',size=7.2,ha='center')
 tt(a,.775,.726,'History',ha='center');tt(a,.90,.726,'Shared',ha='center')
 norms=[Normalize(-.1,1),Normalize(-.1,1),TwoSlopeNorm(vmin=-.5,vcenter=0,vmax=.5)];cmaps=[seq,seq,div]
 for i,(row,label,y) in enumerate(zip(vals,names,global_y)):
  tt(a,.493,y,label,ha='right')
  for j,x in enumerate([.535,.604,.673]):
   col=cmaps[j](norms[j](row[j]));a.add_patch(Rectangle((x-.029,y-.0105),.058,.021,fc=col,ec='none'))
   lum=.2126*col[0]+.7152*col[1]+.0722*col[2];tt(a,x,y,f'{row[j]:.2f}',ha='center',color='white' if lum<.55 else INK)
  for x,val,col in [(.744,row[3],NAVY),(.86,row[4],MID)]:
   a.plot([x,x+.077*val/100],[y,y],c=col,lw=2.6,solid_capstyle='butt');tt(a,x+.077*val/100+.007,y,str(int(val)))
 # Real state lanes terminate at the row labels and measured evidence.
 lx=timeline(a,.04,.235,.686,[0,8,20],[(8,20)],{20},'Lateral root · GSE158761')
 for y in global_y[:9]:
  a.plot([lx[8],lx[20]],[y,y],c=DUST,lw=.5,alpha=.65)
  a.plot(lx[8],y,'o',ms=1.8,c=NAVY);a.plot(lx[20],y,'o',ms=2,mfc='white',mec=NAVY,mew=.5)
 # The alternative-state fan is shown once, as a contrast to a direct temporal lane.
 path(a,(lx[8]-.032,global_y[1]),(lx[8]+.02,global_y[1]),(lx[20]-.03,global_y[0]),(lx[20]-.005,global_y[0]),MID,.5,.28)
 tt(a,.50,.682,'9/9 full-response wins · exact P = 0.00195',color=NAVY)
 tt(a,.035,.50,'Root layers',size=7.2)
 tt(a,.035,.475,'GSE46205',size=6.5,color=SECONDARY)
 for ytime,st,en,idx,header in [(.405,1,3,slice(9,13),'1→3 h · 4/4 full-response wins'),(.264,3,8,slice(13,17),'3→8 h · 4/4 full-response wins')]:
  xx={t:.035+.24*t/8 for t in [0,1,3,8]}
  a.plot([xx[0],xx[8]],[ytime,ytime],c=PALE,lw=.6)
  for t in [0,1,3,8]:
   a.plot(xx[t],ytime,'o',ms=2.7,mfc='white' if t in [3,8] else MID,mec=NAVY if t in [3,8] else MID,mew=.6);tt(a,xx[t],ytime+.023,str(t),ha='center')
  tt(a,xx[8]+.023,ytime+.023,'h')
  path(a,(xx[st],ytime+.005),(xx[st]+.02,ytime+.052),(xx[en]-.02,ytime+.052),(xx[en],ytime+.005),NAVY,.8)
  for y in global_y[idx]:
   a.plot([xx[st],xx[en]],[y,y],c=DUST,lw=.5,alpha=.65);a.plot(xx[st],y,'o',ms=1.8,c=NAVY);a.plot(xx[en],y,'o',ms=2,mfc='white',mec=NAVY,mew=.5)
  tt(a,.50,ytime+.001,header,color=NAVY)
 # Unit keys are shared across all context groups, not repeated per dataset.
 for rect,cmap,norm,title,ticks in [([.51,.074,.081,.015],seq,Normalize(-.1,1),'Full ρ',[-.1,1]),([.628,.074,.081,.015],div,Normalize(-.5,.5),'Δρ',[-.5,0,.5])]:
  cb=f.add_axes(rect);cb.set_label('shared color ruler');mpl.colorbar.ColorbarBase(cb,cmap=cmap,norm=norm,orientation='horizontal',ticks=ticks);cb.tick_params(length=1.8,labelsize=6.5,width=.5);cb.set_title(title,fontsize=6.5,pad=3)
 for x in [.744,.86]:
  a.plot([x,x+.077],[.083,.083],c=SECONDARY,lw=.6)
  for val in [0,100]:
   xx=x+.077*val/100;a.plot([xx,xx],[.083,.078],c=SECONDARY,lw=.5);tt(a,xx,.058,str(val),ha='center')
 tt(a,.035,.091,'GSE46205: corrected_post_reveal',color=SECONDARY)
 tt(a,.32,.018,'Matched null: 13/16 significant; lateral root 6/9.',color=SECONDARY)
 pd.DataFrame(vals,index=names,columns=['history_rho','shared_rho','delta_rho','history_top100','shared_top100']).to_csv(ROOT/'figure2_context_values.tsv',sep='\t')
 STRUCTURE[f]={'canvas_mm':[180,120],'visual_path':VISUAL_PATH['Figure2'],'actual_time_systems':3,'contexts':17,'full_quantitative_axes':1,'linked_objects':'Real time lanes share context rows with rho cells and Top100 bars. Salt/sulfur are embedded.','scientific_boundary':'6h and 12h salt targets belong to separate forecast tasks; 6h history becomes available only for the 12h task. Both GSE46205 blocks remain corrected_post_reveal.'}

def figure3():
 f=canvas('Figure3_v9_preview',124)
 um=f.add_axes([.045,.60,.365,.33]);um.set_label('real root UMAP')
 rootmap(um,col=BRICK,labels=False,size=1.4)
 medians={}
 for n in [5,6,8,2,3,10]:
  z=em[em.louvain_component.eq(n)].iloc[:,:2].median();medians[n]=z.to_numpy()
  um.text(z.iloc[0],z.iloc[1],str(n),fontsize=6.5,ha='center',va='center',bbox=dict(fc='white',ec='none',pad=.6))
 a=scene(f)
 tt(a,.02,.972,'a',size=9.5);tt(a,.052,.972,'Root states → target recovery',size=7.5)
 tt(a,.05,.558,'Sources 5/6/8 → targets 2/3/10; 830 HS cells sealed',color=SECONDARY)
 tt(a,.54,.928,'Top100 recovery against source-matched null',size=7.2)
 tt(a,.54,.894,'4/4 BH q ≈ 1e-4',color=SECONDARY)
 f.canvas.draw()
 # Target readout positions are linked to actual UMAP label heights, not a new grid.
 uv={n:f.transFigure.inverted().transform(um.transData.transform(medians[n])) for n in [2,3,10]}
 targetorder=sorted([2,3,10],key=lambda n:uv[n][1],reverse=True)
 rows=[(n,float(uv[n][1])) for n in targetorder];rows.append(('seedling',min(y for _,y in rows)-.064))
 p=wide(F[2],'heatshock_matched_null');rho={2:.615,3:.621,10:.450,'seedling':.415}
 keys=list(p.index)
 for n,y in rows:
  key=next(k for k in keys if (str(n) in str(k) if n!='seedling' else 'seedling' in str(k)))
  r=p.loc[key];name='Seedling' if n=='seedling' else 'Cluster '+str(n)
  tt(a,.495,y,name,ha='right')
  x0=.54;w=.19;a.plot([x0,x0+w*r.null_q99/100],[y,y],c=PALE,lw=3.4);a.plot(x0+w*r.null_mean/100,y,'o',ms=2.4,c=MID)
  ox=x0+w*r.observed_overlap/100;a.plot(ox,y,marker='D' if n=='seedling' else 'o',ms=4.2,c=BRICK);tt(a,ox+.009,y,str(int(r.observed_overlap)),color=BRICK)
  tt(a,.765,y,f'ρ = {rho[n]:.3f}',color=SECONDARY)
  if n!='seedling':
   ux,uy=uv[n];path(a,(ux+.013,uy),(ux+.045,uy),(.405,y),(.435,y),MID,.6,.28)
 # One compact ruler serves the complete readout.
 ry=.859;a.plot([.54,.73],[ry,ry],c=SECONDARY,lw=.6)
 for count in [0,50,100]:
  x=.54+.19*count/100;a.plot([x,x],[ry,ry-.006],c=SECONDARY,lw=.5);tt(a,x,ry-.025,str(count),ha='center')
 # The entire lower scene is one transfer object, including its gene annotations.
 sc=f.add_axes([.08,.12,.49,.365]);sc.set_label('all-gene cross-study transfer')
 genecloud(sc);sc.axhline(0,c=PALE,lw=.5,ls='--',zorder=0)
 sc.set(xlabel='Frozen root heat-response score',ylabel='Independent seedling response\n(Δ log2[CPM + 1])');sc.tick_params(labelsize=6.5,length=2.5,width=.6)
 sc.set_xlim(sc.get_xlim());sc.set_ylim(sc.get_ylim())
 tt(a,.02,.518,'b',size=9.5);tt(a,.052,.518,'Frozen program → independent whole seedlings',size=7.5)
 sc.text(.055,.90,'ρ = 0.415\nTop100 overlap = 64\nOverlap direction = 100%',transform=sc.transAxes,fontsize=6.5,va='top')
 display=chosen.head(14).copy();leaders=display.head(6).sort_values('target_response',ascending=False)
 order=pd.concat([leaders,display.iloc[6:]],ignore_index=True);assert len(order)==14 and order.gene.nunique()==14
 tt(a,.62,.535,'Shared genes',size=7.2);tt(a,.734,.512,'Root',ha='center');tt(a,.782,.512,'Seedling',ha='center')
 gy=np.linspace(.484,.118,14);f.canvas.draw()
 for i,((_,r),y) in enumerate(zip(order.iterrows(),gy)):
  tt(a,.62,y,r.gene,color=BRICK if i<6 else SECONDARY)
  # Separate original-unit scales, consistently retained from v8.
  for x,val,maxval,col in [(.72,r.shared_prediction,8,DBRICK),(.766,r.target_response,13,BRICK)]:
   assert val>=0 and val<=maxval
   a.plot([x,x+.031],[y,y],c=PALE,lw=.6);a.plot([x,x+.031*val/maxval],[y,y],c=col,lw=2.6,solid_capstyle='butt')
  if i<6:
   px,py=f.transFigure.inverted().transform(sc.transData.transform((r.shared_prediction,r.target_response)))
   path(a,(px,py),(px+.017,py),(.593,y),(.613,y),MID,.55,.55)
 order.to_csv(ROOT/'figure3_linked_gene_display.tsv',sep='\t',index=False)
 tt(a,.62,.073,'Bar scales: root 0–8; seedling 0–13',color=SECONDARY)
 # Specificity is a boundary strip attached to the same transfer scene, no new panel.
 spec=f.add_axes([.895,.14,.09,.33]);spec.set_label('attached specificity strip')
 p=wide(F[2],'cross_experiment_wrong_perturbation_controls')
 for i,(_,r) in enumerate(p.iterrows()):spec.plot([r.wrong_program_rho,r.true_heat_rho],[i,i],c=MID,lw=.7);spec.scatter(r.wrong_program_rho,i,c=MID,s=11);spec.scatter(r.true_heat_rho,i,c=BRICK,s=17,zorder=3)
 spec.set_yticks(range(7),['Salt 6 h','Salt 3 h','Pep1','SCOOP10','SCREW1','flg22','Chitin']);spec.set_ylim(6.5,-.5);ruler(spec,(-.1,.75),[0,.5],'ρ');spec.tick_params(labelsize=6.5,length=2,width=.5)
 spec.set_title('Specificity',loc='left',fontsize=7.2,pad=7)
 STRUCTURE[f]={'canvas_mm':[180,124],'visual_path':VISUAL_PATH['Figure3'],'cells':2085,'sealed_cells':830,'genes':16335,'top100_overlap':64,'linked_gene_count':14,'point_to_gene_leaders':6,'target_to_recovery_leaders':3,'recovery_order':'Actual UMAP target-label height; independent seedling follows below','selection':'Same joint-rank rule as v8: first 14 of 64 shared Top100 genes. First six displayed in descending target response to reduce leader crossings.','bar_scales':'Root 0–8; independent seedling 0–13, original column units.'}

def check_render(f):
 f.canvas.draw();renderer=f.canvas.get_renderer();fb=f.bbox;outside=[]
 for t in f.findobj(match=Text):
  if not t.get_visible() or not t.get_text().strip():continue
  b=t.get_window_extent(renderer)
  if b.x0<fb.x0-1 or b.y0<fb.y0-1 or b.x1>fb.x1+1 or b.y1>fb.y1+1:outside.append({'text':t.get_text(),'bbox':list(b.bounds)})
 rects=[{'index':i,'object':a.get_label(),'bbox':list(a.get_position().bounds)} for i,a in enumerate(f.axes)]
 overlaps=[]
 for i,aa in enumerate(f.axes):
  for j,bb in enumerate(f.axes[i+1:],i+1):
   r,s=aa.get_position(),bb.get_position();w=min(r.x1,s.x1)-max(r.x0,s.x0);h=min(r.y1,s.y1)-max(r.y0,s.y0)
   if w*f.get_figwidth()*72>1 and h*f.get_figheight()*72>1:
    intentional='linked annotation / data field' in [aa.get_label(),bb.get_label()]
    overlaps.append({'axes':[i,j],'intentional':intentional,'reason':'Transparent shared data/annotation scaffold' if intentional else 'Unexpected object overlap'})
 return {'outside_canvas':outside,'axes_count':len(f.axes),'objects':rects,'axes_overlaps':overlaps}

if __name__=='__main__':
 figure2();figure3();report={}
 for name,f in figs:
  checks=check_render(f);report[name]={'design':STRUCTURE[f],**checks}
  if name.startswith('Figure2'):
   # Comparable unit rulers alone share a span; all other objects deliberately differ.
   aa=f.axes[-2:]
   require_matplotlib_panel_alignment(f,axes=aa,panel_ids={ax:str(i) for i,ax in enumerate(aa)},row_groups=[['0','1']],column_groups=[],json_out=QA/(name+'.ruler-alignment.json'),strict=True)
  f.savefig(ROOT/(name+'.png'),dpi=600)
  # Internal vector QA only; final publication vectors await approval of direction.
  f.savefig(QA/(name+'.pdf'))
  (QA/(name+'.alignment.json')).write_text(json.dumps({'status':'intentional_shared_scaffold','reason':'Linked data fields are parts of a continuous scene, with no repeated equal-span subplot group. Numerical row alignment is shared by construction.','objects':checks['objects']},indent=2),encoding='utf-8')
  print(name,'canvas_mm',STRUCTURE[f]['canvas_mm'],'axes/insets',checks['axes_count'],'clipping',len(checks['outside_canvas']),'overlaps',checks['axes_overlaps'],flush=True)
  for o in checks['objects']:print('OBJECT',o,flush=True)
 (QA/'renderer_checks.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
 for p,h in INPUTS.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
