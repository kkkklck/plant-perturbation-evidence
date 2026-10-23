"""Figure 4 only: compact stimulus glyph rows and quiet held-out endpoints.

Contract: task structure defines evidence eligibility; frozen source capability
supports the existing output-resolution grouping. Five stimuli are independent
units; twenty target states remain nested observations. No new inference.
Archetype: asymmetric mixed-modality figure; subordinate task references,
dominant stimulus rows, attached original-unit confidence intervals.
"""
from pathlib import Path
V15=Path(__file__).resolve().parent.parent/'figures_v15_preview'
exec((V15/'build_v15.py').read_text(encoding='utf-8').split("\nif __name__=='__main__':")[0])
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans'],'svg.fonttype':'none','pdf.fonttype':42})

def figure4():
    f=canvas('Figure4_v16_preview',108);a=scene(f)
    title(a,.050,.977,'a','What is missing?')
    title(a,.315,.977,'b','Source capability and held-out transfer')
    # Small, unchanged data objects reference the full earlier task evidence.
    tp=pd.read_csv(V8/'temporal_thumbnail_plot_values.tsv',sep='\t')
    tv=tp[['r1','r3','r6']].to_numpy();lim=float(np.ceil(np.abs(tv).max()))
    tt(a,.025,.918,'Future time',size=7)
    ta=f.add_axes([.074,.804,.125,.085]);ta.set_label('temporal task reference')
    ta.imshow(tv,aspect='auto',cmap=div,vmin=-lim,vmax=lim)
    for i in range(3):
        for j in range(3):ta.text(j,i,f'{tv[i,j]:.1f}',fontsize=6.5,ha='center',va='center',color='white' if tv[i,j]/lim>.65 else INK)
    ta.set_yticks(range(3),tp.gene,fontsize=6.5);ta.set_xticks(range(3),['1 h','3 h','6 h'],fontsize=6.5);clean(ta)
    tt(a,.025,.750,'Recent same-state history',color=NAVY)
    tt(a,.025,.705,'Unseen state',size=7)
    um=f.add_axes([.020,.530,.180,.130]);um.set_label('state task reference')
    rootmap(um,col=BRICK,labels=False,size=.32)
    tt(a,.025,.500,'Shared cross-state response',color=SECONDARY)
    tt(a,.025,.445,'Unseen experiment',size=7)
    sc=f.add_axes([.025,.270,.180,.120]);sc.set_label('frozen program task reference')
    genecloud(sc,tiny=True)
    tt(a,.025,.235,'Frozen cross-study program',color=SECONDARY)

    # Five identity rows. No connecting lines across descriptors or to targets.
    order=['Pep1','SCOOP10','SCREW1','flg22','Chitin']
    colors=[NAVY,DUST,'#34576F',BRICK,DBRICK];markers=['o','s','^','D','p']
    ys=[.780,.670,.560,.400,.290]
    metrics=['rank1_energy','median_pairwise_rho','effective_state_fraction','shared_rms']
    labels=['Rank-1\nenergy','Pairwise\nρ','Effective\nfraction','Shared\nRMS']
    xs=[.415,.503,.591,.679];gw=.044;target_left=.770;target_width=.205
    S=wide(F[3],'stimulus_level_source_vs_target').reindex(order);raw=S[metrics]
    # Each source glyph has a zero baseline; rho/fractions use 0–1.
    # RMS retains a separate original-unit ruler, 0–0.160. No joint score.
    maxima=[1.,1.,1.,.160]
    assert (raw.to_numpy()>=0).all() and (raw.to_numpy()<=np.asarray(maxima)+1e-12).all()
    T=wide(sel(F[3],section='immune_frozen_reveal',context=['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']),idx=['dataset','context'])
    tt(a,.282,.926,'Higher reliability: gene candidates + direction supported',color=NAVY)
    tt(a,.415,.896,'Source-only glyphs: 0–1; RMS: 0–0.160',color=SECONDARY)
    tt(a,.770,.943,'Grey: states; colour: mean',color=SECONDARY)
    tt(a,target_left+target_width/2,.913,'Held-out target ρ',ha='center',size=7)
    for x,label in zip(xs,labels):tt(a,x+gw/2,.858,label,ha='center',size=7)
    ry=.852
    a.plot([target_left,target_left+target_width],[ry,ry],c=SECONDARY,lw=.55)
    for value in [0,.5,1]:
        x=target_left+target_width*value;a.plot([x,x],[ry,ry-.007],c=SECONDARY,lw=.5)
        tt(a,x,.878,f'{value:g}',ha='center',color=SECONDARY)
    source_map=[];target_map=[];targets=[]
    for i,(stim,col,mk,y) in enumerate(zip(order,colors,markers,ys)):
        a.add_patch(Rectangle((.268,y-.046),.714,.092,fc='#F2F6F8' if i<3 else '#FBF4F1',ec='none',zorder=-5))
        a.scatter(.282,y,c=col,s=15,marker=mk,zorder=3)
        tt(a,.296,y,stim,size=7,color=col)
        for j,(x,m,maximum) in enumerate(zip(xs,metrics,maxima)):
            v=float(raw.loc[stim,m]);fraction=v/maximum;gy=y-.016
            a.plot([x,x+gw],[gy,gy],c='#DCE2E6',lw=.6,solid_capstyle='round',zorder=0)
            a.plot([x,x+gw*fraction],[gy,gy],c=col,lw=1.2,solid_capstyle='round',zorder=1)
            a.scatter(x+gw*fraction,gy,c=col,s=3,zorder=2)
            tt(a,x+gw/2,y+.014,format(v,'.3f' if j==3 else '.2f'),ha='center',color=SECONDARY)
            source_map.append(dict(stimulus=stim,metric=m,raw_value=v,raw_min=0.,raw_max=maximum,display_fraction=fraction,glyph_x=x,glyph_end_x=x+gw*fraction,glyph_y=gy,row_y=y))
        ax=f.add_axes([target_left,y-.027,target_width,.054]);ax.set_label(stim+' held-out targets')
        ax.patch.set_alpha(0);ax.set_xlim(0,1);ax.set_ylim(-.6,.6);ax.axis('off')
        v=T.loc[stim].full_rho;mean=float(v.mean());offsets=np.array([-.22,-.075,.075,.22])
        ax.scatter(v.to_numpy(),offsets,c='#B9C1C9',s=11,lw=0,zorder=2)
        ax.scatter(mean,0,c=col,s=28,marker=mk,lw=.4,edgecolor='white',zorder=3)
        targets.append(ax)
        for (context,val),oy in zip(v.items(),offsets):
            target_map.append(dict(stimulus=stim,context=context,target_rho=val,stimulus_mean=mean,row_y=y,point_x=target_left+target_width*val,point_y=y+.054*oy/1.2,mean_x=target_left+target_width*mean,mean_y=y))
    tt(a,.282,.473,'Lower reliability: reduce output resolution',color=BRICK)
    tt(a,.282,.213,'95% bootstrap CI · highest independent unit: stimulus (n = 5)',color=SECONDARY)
    # Original frozen hierarchical contrasts; no new panel or fitted threshold.
    H=wide(F[3],'hierarchical_stimulus_inference');cis=[]
    for rect,key,lims,ticks in [([.420,.120,.190,.036],'full_rho',(-.05,.45),[0,.4]),([.750,.120,.210,.036],'top100_overlap',(-5,50),[0,40])]:
        ax=f.add_axes(rect);ax.set_label('attached CI '+key);r=H.loc[key];cis.append(ax)
        ax.errorbar(r.difference_allow_minus_reject,0,xerr=[[r.difference_allow_minus_reject-r.bootstrap_ci025],[r.bootstrap_ci975-r.difference_allow_minus_reject]],fmt='o',c=NAVY,ms=3,capsize=2,lw=.8)
        ax.set_ylim(-1,1);ax.set_yticks([]);ruler(ax,lims,ticks,'');ax.tick_params(labelsize=7,length=2,width=.5);ax.axvline(0,c=MID,lw=.5,ls='--')
        fmt='.3f' if key=='full_rho' else '.1f'
        ax.set_title(('Δρ ' if key=='full_rho' else 'ΔTop100 ')+format(r.difference_allow_minus_reject,'+'+fmt)+' ['+format(r.bootstrap_ci025,fmt)+', '+format(r.bootstrap_ci975,fmt)+']',fontsize=7,pad=3)
    tt(a,.282,.060,'Exact one-sided stimulus-label permutation P = 0.10 · no fitted thresholds',color=SECONDARY)
    tt(a,.025,.036,'Geometry and response amplitude covary; n = 5 does not identify an independent geometric effect.',color=SECONDARY)
    tt(a,.025,.012,'19/20 targets exceed source-matched nulls; lower reliability retains transferable information.',color=SECONDARY)
    raw.to_csv(ROOT/'figure4_source_raw_values.tsv',sep='\t')
    pd.DataFrame(source_map).to_csv(ROOT/'figure4_source_glyph_mapping.tsv',sep='\t',index=False)
    pd.DataFrame(target_map).to_csv(ROOT/'figure4_target_row_mapping.tsv',sep='\t',index=False)
    H.to_csv(ROOT/'figure4_existing_inference.tsv',sep='\t');tp.to_csv(ROOT/'figure4_temporal_thumbnail.tsv',sep='\t',index=False)
    STRUCTURE[f]={'canvas_mm':[180,108],'hero':'Five stimulus rows with compact source glyphs and quiet target endpoints','stimuli':5,'targets':20,'source_ranges':dict(zip(metrics,maxima)),'target_range':[0,1],'source_lines_across_metrics':0,'highest_independent_unit':'stimulus n=5'}
    return f,targets,cis

def protected_files():
    paths=list((PROJECT/'figures_v15_preview').glob('Figure*_v15_preview.*'))
    paths+=list((PROJECT/'figures_v11_preview').glob('Figure[23]_v11_preview.*'))
    paths+=list(PROJECT.glob('*.docx'))
    return {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

if __name__=='__main__':
    protected=protected_files();f,targets,cis=figure4()
    check=check_render(f);assert not check['outside_canvas'],check
    ids=[ax.get_label() for ax in targets]
    require_matplotlib_panel_alignment(f,axes=targets,panel_ids=dict(zip(targets,ids)),column_groups=[ids[:3],ids[3:],[ids[0],ids[3]]],json_out=QA/'Figure4_v16_preview.alignment.json',strict=True)
    # Heterogeneous units and deliberately unequal task references are excluded
    # from comparable-target geometry. Same-height CI appendages share a row.
    require_matplotlib_panel_alignment(f,axes=cis,panel_ids=dict(zip(cis,['CI rho','CI hits'])),row_groups=[['CI rho','CI hits']],json_out=QA/'Figure4_v16_preview.CI_alignment.json',strict=True)
    f.savefig(ROOT/'Figure4_v16_preview.png',dpi=600)
    f.savefig(ROOT/'Figure4_v16_preview.pdf')
    f.savefig(ROOT/'Figure4_v16_preview.svg')
    (QA/'renderer_checks.json').write_text(json.dumps({'design':STRUCTURE[f],**check,'alignment_exclusions':'Task references have different data modalities and intentionally unequal extents; miniature source glyphs are marks on one scene canvas.'},indent=2),encoding='utf-8')
    (ROOT/'input_hashes.json').write_text(json.dumps(INPUTS,indent=2),encoding='utf-8')
    (QA/'protected_output_hashes.json').write_text(json.dumps(protected,indent=2),encoding='utf-8')
    for p,h in {**INPUTS,**protected}.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
    print('Rendered Figure 4 v16 only; scientific inputs and frozen outputs unchanged.',flush=True)
