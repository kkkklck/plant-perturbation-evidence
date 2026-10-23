"""v14: rebuild the spatial skeleton, retaining all frozen scientific values.
VISUAL_PATH F1: organic tissue landscape -> recurrence mosaic -> statistical shadow.
VISUAL_PATH F4: surrounding task anchors -> staggered stimulus fingerprints with
attached target points -> interpretation -> compact inference footer.
No new analysis, no new selection, no scientific-source mutations, no Word work.
"""
from pathlib import Path
V13=Path(__file__).resolve().parent.parent/'figures_v13_preview'
base=(V13/'build_v13.py').read_text(encoding='utf-8')
exec(base.split('def figure1():')[0])
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],'svg.fonttype':'none','pdf.fonttype':42})
code=base.split('def figure1():')[1].split('def figure4():')[0]
code='def figure1():'+code
changes={
 "canvas('Figure1_v13_preview',100)":"canvas('Figure1_v14_preview',108)",
 "title(a,.43,.967,'b'":"title(a,.52,.982,'b'",
 '[.025,.61,.27,.27]':'[.010,.420,.420,.480]',
 "grays=['#C6CFD6','#98A6B0','#E0E5E9','#B4C0C8']":"grays=['#8EA6B8','#A89FB7','#A2BAB3','#B4BCC6']",
 "s=.2,c=col,lw=0":"s=.22,c=col,lw=0,alpha=.65",
 "tt(a,.025,.566,'4 tissue systems × 4 stress conditions',size=6.5)":"tt(a,.025,.931,'4 tissues × 4 stresses · 16 contexts',size=7)",
 " tt(a,.025,.535,'16 matched contexts',color=SECONDARY)\n":'',
 '[.43,.36,.54,.54]':'[.52,.275,.46,.640]',
 "ListedColormap(['white',PALE,'#78838D'])":"ListedColormap(['#FDFEFE','#CAD7E1',NAVY])",
 "xx=.43+.54*(i+.5)/4;tt(a,xx,.928":"xx=.52+.46*(i+.5)/4;tt(a,xx,.946",
 "(.43+i*.135,.910),.135,.006":"(.52+i*.115,.925),.115,.006",
 "[(.43,'white','Neither'),(.58,PALE,'One replicate'),(.77,'#78838D','Both replicates')]":"[(.52,'#FDFEFE','Neither'),(.65,'#CAD7E1','One replicate'),(.81,NAVY,'Both replicates')]",
 "(x,.297)":"(x,.219)","tt(a,x+.017,.303,lab)":"tt(a,x+.017,.225,lab)",
 "tt(a,.43,.270,'Matched null · 16/16 BH q ≈ 1e-4',size=7)":"tt(a,.52,.195,'Matched null · 16/16 BH q ≈ 1e-4',size=7)",
 '[.43,.130,.54,.085]':'[.52,.107,.46,.060]',
 'ylim=(0,36)':'ylim=(-2,40)',"nu.set_yticks([0,30])":"nu.set_yticks([0,20,40])",
 "tt(a,.43,.105,'Top100 overlap · treatment median null enrichment below (×)',color=SECONDARY)":"tt(a,.52,.081,'Top100 overlap · median null enrichment below (×)',color=SECONDARY)",
 "tt(a,.43+.54*(i+.5)/4,.080":"tt(a,.52+.46*(i+.5)/4,.056",
 " tt(a,.43,.241,'Grey: null mean / 99th percentile; coloured: observed',color=SECONDARY)\n":'',
 "tt(a,.43,.040,'E: epidermis · M: meristem · P: mesophyll · V: vasculature',color=SECONDARY)":"tt(a,.52,.034,'Grey: null mean / 99th percentile; coloured: observed',color=SECONDARY)",
 "tt(a,.43,.015,'Matched null controls baseline expression and instability',color=SECONDARY)":"tt(a,.52,.010,'Matched on baseline expression and instability',color=SECONDARY)",
 "'canvas_mm':[180,100]":"'canvas_mm':[180,108]"
}
for old,new in changes.items():
 assert old in code,old
 code=code.replace(old,new)
# Tissue-letter colours echo the actual embedding identities.
code=code.replace("hm.set_xlim(-.5,15.5)","hm.set_xlim(-.5,15.5)\n for k,t in enumerate(hm.get_xticklabels()):t.set_color(['#657F94','#837493','#668A7D','#737F90'][k%4])")
exec(code)

def figure4():
 f=canvas('Figure4_v14_preview',108);a=scene(f)
 title(a,.05,.970,'a','Task evidence');title(a,.40,.970,'b','Stimulus transfer field')
 # Actual task objects surround the transfer field at distinct scales/positions.
 tp=pd.read_csv(V8/'temporal_thumbnail_plot_values.tsv',sep='\t');tv=tp[['r1','r3','r6']].to_numpy();lim=float(np.ceil(np.abs(tv).max()))
 tt(a,.035,.918,'Future time',size=7)
 ta=f.add_axes([.100,.802,.230,.078]);ta.set_label('frozen temporal thumbnail');ta.imshow(tv,aspect='auto',cmap=div,vmin=-lim,vmax=lim)
 for i in range(3):
  for j in range(3):ta.text(j,i,f'{tv[i,j]:.1f}',fontsize=6.5,ha='center',va='center',color='white' if tv[i,j]/lim>.65 else INK)
 ta.set_yticks(range(3),tp.gene,fontsize=6.5);ta.set_xticks(range(3),['1 h','3 h','6 h'],fontsize=6.5);clean(ta)
 tt(a,.08,.757,'Recent same-state history',color=NAVY)
 tt(a,.025,.706,'Unseen state',size=7)
 um=f.add_axes([.01,.395,.32,.275]);um.set_label('root biological anchor');rootmap(um,col=BRICK,labels=False,size=.6)
 tt(a,.025,.360,'Shared cross-state response',color=SECONDARY)
 tt(a,.24,.292,'Unseen experiment',size=7)
 sc=f.add_axes([.245,.129,.285,.135]);sc.set_label('independent transfer anchor');genecloud(sc,tiny=True)
 tt(a,.24,.094,'Frozen cross-study program',color=SECONDARY)
 # Two offset lane clusters, not a rectangular five-row matrix.
 order=['Pep1','SCOOP10','SCREW1','flg22','Chitin'];ys=[.79,.70,.61,.435,.340];origins=[.40]*5
 S=wide(F[3],'stimulus_level_source_vs_target').reindex(order);metrics=['rank1_energy','median_pairwise_rho','effective_state_fraction','shared_rms'];raw=S[metrics];scaled=(raw-raw.min())/(raw.max()-raw.min())
 T=wide(sel(F[3],section='immune_frozen_reveal',context=['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']),idx=['dataset','context'])
 rows=[];fingerprints=[]
 # Both clusters retain the same metric identity/order; duplicate local labels
 # are necessary because the lower fingerprints have a different x origin.
 for ox,hy in [(.40,.853)]:
  for dx,lab in zip([.035,.105,.175,.245],['Rank-1','Pairwise ρ','Eff. fraction','RMS']):tt(a,ox+dx,hy,lab,ha='center',size=6.5)
 for i,(s,y,ox) in enumerate(zip(order,ys,origins)):
  col=NAVY if i<3 else BRICK
  tt(a,ox-.012,y,s,ha='right',size=7,color=col)
  for j,dx in enumerate([.035,.105,.175,.245]):
   x=ox+dx;left=x-.0275;length=.055*scaled.iloc[i,j]
   a.plot([left,left+.055],[y+.009,y+.009],c=PALE,lw=2.6,solid_capstyle='butt')
   if length>0:a.plot([left,left+length],[y+.009,y+.009],c=col,lw=2.6,solid_capstyle='butt')
   tt(a,x,y-.014,f'{raw.iloc[i,j]:.3f}' if j==3 else f'{raw.iloc[i,j]:.2f}',ha='center')
   fingerprints.append(dict(stimulus=s,metric=metrics[j],raw_value=raw.iloc[i,j],within_column_minmax=scaled.iloc[i,j],bar_start_x=left,bar_end_x=left+length))
  # The SAME stimulus row extends into its four real held-out outcomes.
  v=T.loc[s].full_rho;x0=.707;w=.180;off=np.linspace(-.010,.010,4)
  a.scatter(x0+w*v.to_numpy(),y+off,c=MID,s=12,zorder=3);a.scatter(x0+w*v.mean(),y,c=col,s=24,zorder=4)
  for (ct,val),oy in zip(v.items(),off):rows.append(dict(stimulus=s,context=ct,target_rho=val,stimulus_mean=v.mean(),source_row_y=y,target_row_y=y,point_y=y+oy,point_x=x0+w*val,mean_x=x0+w*v.mean(),lane_origin_x=ox))
 # A single quantitative ruler supplies a common rho length scale for all lanes.
 x0=.707;w=.180;ry=.893;a.plot([x0,x0+w],[ry,ry],c=SECONDARY,lw=.6)
 for v in [0,.5,1]:
  x=x0+w*v;a.plot([x,x],[ry,ry-.005],c=SECONDARY,lw=.5);tt(a,x,ry+.019,str(v),ha='center')
 tt(a,.90,.878,'Target ρ',size=7)
 tt(a,.41,.933,'Higher reliability: gene candidates + direction supported',color=NAVY)
 tt(a,.46,.548,'Grey: 4 targets; larger dot: stimulus mean',color=SECONDARY)
 tt(a,.55,.223,'Lower reliability: reduce output resolution',color=BRICK)
 tt(a,.40,.291,'Bars: within-column min–max; raw values shown',color=SECONDARY)
 # Inference floats within the lower-right space rather than a full-width footer panel.
 H=wide(F[3],'hierarchical_stimulus_inference');ci=[]
 for rect,key,lims,ticks in [([.58,.113,.17,.035],'full_rho',(-.05,.45),[0,.4]),([.79,.155,.18,.035],'top100_overlap',(-5,50),[0,40])]:
  ax=f.add_axes(rect);ax.set_label('original inference '+key);r=H.loc[key];ci.append(ax)
  ax.errorbar(r.difference_allow_minus_reject,0,xerr=[[r.difference_allow_minus_reject-r.bootstrap_ci025],[r.bootstrap_ci975-r.difference_allow_minus_reject]],fmt='o',c=NAVY,ms=3,capsize=2,lw=.8)
  ax.set_ylim(-1,1);ax.set_yticks([]);ruler(ax,lims,ticks,'');ax.tick_params(labelsize=7,length=2,width=.5);ax.axvline(0,c=MID,lw=.5,ls='--')
  fmt='.3f' if key=='full_rho' else '.1f';ax.set_title(('Δρ ' if key=='full_rho' else 'ΔTop100 ')+format(r.difference_allow_minus_reject,'+'+fmt)+' ['+format(r.bootstrap_ci025,fmt)+', '+format(r.bootstrap_ci975,fmt)+']',fontsize=7,pad=3)
 tt(a,.025,.066,'95% bootstrap CI · highest independent unit: stimulus (n = 5)',color=SECONDARY)
 tt(a,.56,.055,'Exact one-sided permutation P = 0.10 · no fitted thresholds',color=SECONDARY)
 tt(a,.025,.031,'Geometry and amplitude covary; n = 5 does not identify an independent geometric effect.',color=SECONDARY)
 tt(a,.025,.009,'19/20 targets exceed source-matched nulls; lower reliability retains transferable information.',color=SECONDARY)
 raw.to_csv(ROOT/'figure4_source_raw_values.tsv',sep='\t');pd.DataFrame(fingerprints).to_csv(ROOT/'figure4_fingerprint_mapping.tsv',sep='\t',index=False);pd.DataFrame(rows).to_csv(ROOT/'figure4_target_row_mapping.tsv',sep='\t',index=False);H.to_csv(ROOT/'figure4_existing_inference.tsv',sep='\t');tp.to_csv(ROOT/'figure4_temporal_thumbnail.tsv',sep='\t',index=False)
 ALIGN[f]=None;STRUCTURE[f]={'canvas_mm':[180,108],'stimuli':5,'nested_targets':20,'highest_independent_unit':'stimulus n=5','source_geometry':'Short-bar fingerprints with printed raw values, two unequal-spaced lane clusters, no matrix tiles.','target_geometry':'One global x=.707+.180*rho scale for all 20 observations and five means.','colour_semantics':'Hue is frozen reliability group; bar lengths use established within-column min-max scaling. No fitted threshold.','archetype':'Central fingerprint field, surrounding task anchors, floating inference','letters':['a','b'],'no_arrows':True}

# Preserve all established exporting, hashing, clipping and comparable-column QA.
exec(base.split('def check_render(f):')[1].split("if __name__=='__main__':")[0].join(['def check_render(f):','']))
if __name__=='__main__':
 protected={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT.parent/'figures_v11_preview').glob('Figure[23]_v11_preview.*')}
 figure1();figure4();reports={}
 for name,f in figs:
  ch=check_render(f);assert not ch['outside_canvas'],ch;reports[name]={'design':STRUCTURE[f],**ch}
  if ALIGN[f] is not None:
   axes,rows,cols,ids=ALIGN[f];require_matplotlib_panel_alignment(f,axes=axes,panel_ids=dict(zip(axes,ids)),row_groups=rows,column_groups=cols,json_out=QA/(name+'.alignment.json'),strict=True)
  else:(QA/(name+'.alignment.json')).write_text(json.dumps({'status':'intentional heterogeneous objects','reason':'Bar fingerprints, biological anchors and unequal floating original-unit CI rulers have different roles. No equal-span subplot group. All target positions use the same global rho transform.','objects':ch['objects']},indent=2),encoding='utf-8')
  f.savefig(ROOT/(name+'.png'),dpi=600);f.savefig(ROOT/(name+'.pdf'));f.savefig(ROOT/(name+'.svg'))
  print(name,STRUCTURE[f],ch,flush=True)
 (QA/'renderer_checks.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2),encoding='utf-8');(ROOT/'input_hashes.json').write_text(json.dumps(INPUTS,ensure_ascii=False,indent=2),encoding='utf-8')
 for p,h in {**INPUTS,**protected}.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
 (QA/'protected_v11_hashes.json').write_text(json.dumps(protected,ensure_ascii=False,indent=2),encoding='utf-8')
