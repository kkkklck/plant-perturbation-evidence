"""Hero-first v15; read-only scientific inputs, Figure 1/4 only.
VISUAL_PATH F1: upper-left diagnostic / lower-left rice landscape -> raised
candidate mosaic -> same-context marginal null. F4: task satellites -> five
identity paths through source rulers -> true target means -> inference appendages.
Parallel-coordinate scaling is display only, never a new score or model.
"""
from pathlib import Path
V14=Path(__file__).resolve().parent.parent/'figures_v14_preview'
exec((V14/'build_v14.py').read_text(encoding='utf-8').split('\ndef figure4():')[0])
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],'svg.fonttype':'none','pdf.fonttype':42})
# Reuse the identical candidate/recurrence extraction; change only its geography.
code=code.replace("Figure1_v14_preview","Figure1_v15_preview")
swaps={
 "title(a,.52,.982,'b'":"title(a,.41,.978,'b'",
 "[.010,.420,.420,.480]":"[.008,.110,.315,.450]",
 "tt(a,.025,.931,'4 tissues × 4 stresses · 16 contexts',size=7)":"tt(a,.025,.590,'4 tissues × 4 stresses · 16 contexts',size=7)",
 "tt(a,.035,.382,'Replicate stability',size=7)":"tt(a,.035,.935,'Replicate stability',size=7)",
 "tt(a,.225,.340":"tt(a,.199,.903", "tt(a,.332,.340":"tt(a,.273,.903",
 "'replicate_response_spearman',.195,.07":"'replicate_response_spearman',.174,.05",
 "'observed_overlap',.295,.07":"'observed_overlap',.250,.047",
 "(.195,.07,-.05,.3,[0,.3])":"(.174,.05,-.05,.3,[0,.3])",
 "(.295,.07,0,36,[0,30])":"(.250,.047,0,36,[0,30])",
 "y=.293-i*.043":"y=.865-i*.030",
 "[.112,.112]":"[.738,.738]", "[.112,.105]":"[.738,.731]", "tt(a,x,.080,str(v)":"tt(a,x,.706,str(v)",
 "tt(a,.025,.041,'Dots: tissues; tick: treatment median'":"tt(a,.025,.059,'Dots: tissues; tick: treatment median'",
 "[.52,.275,.46,.640]":"[.41,.335,.57,.580]",
 "xx=.52+.46*(i+.5)/4":"xx=.41+.57*(i+.5)/4",
 "(.52+i*.115,.925),.115,.006":"(.41+i*.1425,.925),.1425,.006",
 "[(.52,'#FDFEFE','Neither'),(.65,'#CAD7E1','One replicate'),(.81,NAVY,'Both replicates')]":"[(.41,'#FDFEFE','Neither'),(.56,'#CAD7E1','One replicate'),(.76,NAVY,'Both replicates')]",
 "(x,.219)":"(x,.283)","tt(a,x+.017,.225":"tt(a,x+.017,.289",
 "tt(a,.52,.195":"tt(a,.41,.252",
 "[.52,.107,.46,.060]":"[.41,.160,.57,.055]",
 "tt(a,.52,.081":"tt(a,.41,.131",
 "tt(a,.52+.46*(i+.5)/4,.056":"tt(a,.41+.57*(i+.5)/4,.103",
 "tt(a,.52,.034":"tt(a,.41,.072", "tt(a,.52,.010":"tt(a,.41,.045",
}
for old,new in swaps.items():
 assert old in code,old
 code=code.replace(old,new)
exec(code)

def figure4():
 f=canvas('Figure4_v15_preview',108);a=scene(f)
 title(a,.05,.977,'a','What is missing?');title(a,.37,.977,'b','Shared evidence: source capability and transfer')
 tp=pd.read_csv(V8/'temporal_thumbnail_plot_values.tsv',sep='\t');tv=tp[['r1','r3','r6']].to_numpy();lim=float(np.ceil(np.abs(tv).max()))
 tt(a,.025,.925,'Future time',size=7)
 ta=f.add_axes([.073,.822,.130,.070]);ta.set_label('frozen temporal island');ta.imshow(tv,aspect='auto',cmap=div,vmin=-lim,vmax=lim)
 for i in range(3):
  for j in range(3):ta.text(j,i,f'{tv[i,j]:.1f}',fontsize=6.5,ha='center',va='center',color='white' if tv[i,j]/lim>.65 else INK)
 ta.set_yticks(range(3),tp.gene,fontsize=6.5);ta.set_xticks(range(3),['1 h','3 h','6 h'],fontsize=6.5);clean(ta)
 tt(a,.025,.776,'Recent same-state history',color=NAVY)
 tt(a,.025,.724,'Unseen state',size=7)
 um=f.add_axes([.015,.482,.180,.207]);um.set_label('root state island');rootmap(um,col=BRICK,labels=False,size=.4)
 tt(a,.015,.449,'Shared cross-state response',color=SECONDARY)
 tt(a,.025,.396,'Unseen experiment',size=7)
 sc=f.add_axes([.035,.225,.207,.131]);sc.set_label('independent gene island');genecloud(sc,tiny=True)
 tt(a,.025,.190,'Frozen cross-study program',color=SECONDARY)
 # One field, five exact stimulus paths; each source ruler has its OWN raw range.
 order=['Pep1','SCOOP10','SCREW1','flg22','Chitin'];colors=[NAVY,DUST,'#34576F',BRICK,DBRICK];markers=['o','s','^','D','p'];label_y=[.881,.785,.687,.485,.395]
 S=wide(F[3],'stimulus_level_source_vs_target').reindex(order);metrics=['rank1_energy','median_pairwise_rho','effective_state_fraction','shared_rms'];raw=S[metrics]
 mn=raw.min();mx=raw.max();scaled=(raw-mn)/(mx-mn);assert (mx>mn).all()
 T=wide(sel(F[3],section='immune_frozen_reveal',context=['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']),idx=['dataset','context'])
 xx=np.array([.410,.535,.660,.785,.925]);bottom=.380;height=.420
 for j,(x,lab) in enumerate(zip(xx,['Rank-1\nenergy','Pairwise\nρ','Effective\nfraction','Shared\nRMS','Target\nρ'])):
  a.plot([x,x],[bottom,bottom+height],c=MID,lw=.65,zorder=0)
  tt(a,x,.869,lab,ha='center',size=7)
  lo,hi=(mn.iloc[j],mx.iloc[j]) if j<4 else (0,1)
  fmt='.3f' if j==3 else '.2f'
  tt(a,x,.823,format(hi,fmt),ha='center',color=SECONDARY);tt(a,x,.356,format(lo,fmt),ha='center',color=SECONDARY)
  if j==4:
   a.plot([x-.003,x+.003],[bottom+height*.5]*2,c=MID,lw=.5);tt(a,x+.022,bottom+height*.5,'0.5',color=SECONDARY)
 source_map=[];target_map=[]
 for i,(stim,col,marker,ly) in enumerate(zip(order,colors,markers,label_y)):
  v=T.loc[stim].full_rho;mean=float(v.mean());yy=np.r_[bottom+height*scaled.loc[stim].to_numpy(),bottom+height*mean]
  a.plot(xx,yy,c=col,lw=.85,alpha=.85,zorder=2)
  a.scatter(xx[:4],yy[:4],c=col,marker=marker,s=13,zorder=4)
  a.scatter(xx[-1]+np.linspace(-.006,.006,4),bottom+height*v.to_numpy(),c=MID,s=10,zorder=3)
  a.scatter(xx[-1],yy[-1],c=col,marker=marker,s=28,zorder=5)
  # Identity legend carries the raw source tuple; no artificial label-to-node line.
  a.scatter(.245,ly,c=col,marker=marker,s=16);tt(a,.259,ly,stim,size=7,color=col)
  rr=raw.loc[stim];tt(a,.242,ly-.025,f'{rr.iloc[0]:.2f} / {rr.iloc[1]:.2f} / {rr.iloc[2]:.2f} / {rr.iloc[3]:.3f}',color=SECONDARY)
  for j,m in enumerate(metrics):source_map.append(dict(stimulus=stim,metric=m,raw_value=rr[m],raw_min=mn[m],raw_max=mx[m],display_0_1=scaled.loc[stim,m],point_x=xx[j],point_y=yy[j],marker=marker))
  for (ct,val),ox in zip(v.items(),np.linspace(-.006,.006,4)):target_map.append(dict(stimulus=stim,context=ct,target_rho=val,stimulus_mean=mean,point_x=xx[-1]+ox,point_y=bottom+height*val,mean_x=xx[-1],mean_y=yy[-1],path_endpoint_x=xx[-1],path_endpoint_y=yy[-1]))
 tt(a,.415,.929,'Higher reliability: gene candidates + direction supported',color=NAVY)
 tt(a,.300,.322,'Source axes: display 0–1; raw ranges shown at ruler ends',color=SECONDARY)
 tt(a,.300,.297,'Raw tuple order: rank-1 / pairwise ρ / effective fraction / RMS',color=SECONDARY)
 tt(a,.415,.267,'Lower reliability: reduce output resolution',color=BRICK)
 tt(a,.415,.242,'Grey: four target states; coloured endpoint: stimulus mean',color=SECONDARY)
 H=wide(F[3],'hierarchical_stimulus_inference');ci=[]
 for rect,key,lims,ticks in [([.42,.142,.17,.035],'full_rho',(-.05,.45),[0,.4]),([.74,.123,.19,.035],'top100_overlap',(-5,50),[0,40])]:
  ax=f.add_axes(rect);ax.set_label('attached inference '+key);r=H.loc[key];ci.append(ax)
  ax.errorbar(r.difference_allow_minus_reject,0,xerr=[[r.difference_allow_minus_reject-r.bootstrap_ci025],[r.bootstrap_ci975-r.difference_allow_minus_reject]],fmt='o',c=NAVY,ms=3,capsize=2,lw=.8)
  ax.set_ylim(-1,1);ax.set_yticks([]);ruler(ax,lims,ticks,'');ax.tick_params(labelsize=7,length=2,width=.5);ax.axvline(0,c=MID,lw=.5,ls='--')
  fmt='.3f' if key=='full_rho' else '.1f';ax.set_title(('Δρ ' if key=='full_rho' else 'ΔTop100 ')+format(r.difference_allow_minus_reject,'+'+fmt)+' ['+format(r.bootstrap_ci025,fmt)+', '+format(r.bootstrap_ci975,fmt)+']',fontsize=7,pad=3)
 tt(a,.300,.083,'95% bootstrap CI · highest independent unit: stimulus (n = 5)',color=SECONDARY)
 tt(a,.300,.060,'Exact one-sided permutation P = 0.10 · no fitted threshold',color=SECONDARY)
 tt(a,.025,.036,'Geometry and response amplitude covary; n = 5 does not identify an independent geometric effect.',color=SECONDARY)
 tt(a,.025,.012,'19/20 targets exceed source-matched nulls; lower reliability retains transferable information.',color=SECONDARY)
 raw.to_csv(ROOT/'figure4_source_raw_values.tsv',sep='\t');pd.DataFrame(source_map).to_csv(ROOT/'figure4_source_path_mapping.tsv',sep='\t',index=False);pd.DataFrame(target_map).to_csv(ROOT/'figure4_target_path_mapping.tsv',sep='\t',index=False);H.to_csv(ROOT/'figure4_existing_inference.tsv',sep='\t');tp.to_csv(ROOT/'figure4_temporal_thumbnail.tsv',sep='\t',index=False)
 ALIGN[f]=None;STRUCTURE[f]={'canvas_mm':[180,108],'stimuli':5,'nested_targets':20,'highest_independent_unit':'stimulus n=5','field_bounds':[.30,.35,.68,.54],'source_scaling':'Independent source-axis min-max for display only; raw values in identity tuples and min/max ruler labels','target_scaling':'Unscaled 0–1 rho; each path terminates exactly at its true stimulus mean','path_semantics':'Identity across measured descriptors/outcomes, not a temporal trajectory, mechanism, score or fitted model','letters':['a','b']}

def check_render(f):
 f.canvas.draw();ren=f.canvas.get_renderer();fb=f.bbox;outside=[]
 for t in f.findobj(match=Text):
  if not t.get_visible() or not t.get_text().strip():continue
  b=t.get_window_extent(ren)
  if b.x0<fb.x0-1 or b.y0<fb.y0-1 or b.x1>fb.x1+1 or b.y1>fb.y1+1:outside.append(t.get_text())
 return {'outside_canvas':outside,'axes_count':len(f.axes),'objects':[{'object':ax.get_label(),'bbox':list(ax.get_position().bounds)} for ax in f.axes]}

if __name__=='__main__':
 protected={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT.parent/'figures_v11_preview').glob('Figure[23]_v11_preview.*')}
 figure1();figure4();reports={}
 for name,f in figs:
  ch=check_render(f);assert not ch['outside_canvas'],ch;reports[name]={'design':STRUCTURE[f],**ch}
  if ALIGN[f] is not None:
   axes,rows,cols,ids=ALIGN[f];require_matplotlib_panel_alignment(f,axes=axes,panel_ids=dict(zip(axes,ids)),row_groups=rows,column_groups=cols,json_out=QA/(name+'.alignment.json'),strict=True)
  else:(QA/(name+'.alignment.json')).write_text(json.dumps({'status':'intentional continuous field and unequal appendages','reason':'Parallel-coordinate rulers are data objects on one shared canvas. Biological islands and original-unit CI appendages intentionally have different spans.','objects':ch['objects']},indent=2),encoding='utf-8')
  f.savefig(ROOT/(name+'.png'),dpi=600);f.savefig(ROOT/(name+'.pdf'));f.savefig(ROOT/(name+'.svg'));print(name,STRUCTURE[f],ch,flush=True)
 (QA/'renderer_checks.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2),encoding='utf-8');(ROOT/'input_hashes.json').write_text(json.dumps(INPUTS,ensure_ascii=False,indent=2),encoding='utf-8')
 for p,h in {**INPUTS,**protected}.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
 (QA/'protected_v11_hashes.json').write_text(json.dumps(protected,ensure_ascii=False,indent=2),encoding='utf-8')
