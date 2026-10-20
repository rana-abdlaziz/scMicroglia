
ABBR=("AD, Alzheimer's disease; DEG, differentially expressed gene; HPC, hippocampus; log₂FC, log₂ fold change of inflammatory "
      "(MG2/MG8/MG10) microglial nuclei from AD donors versus homeostatic (MG0) nuclei from non-AD donors within each sex; "
      "P_adj, Bonferroni-adjusted P value from Seurat FindMarkers (Poisson generalized linear model with age at death as a latent variable; "
      "values below double precision are given as < 1 × 10⁻³⁰⁰); PFC, prefrontal cortex.")
def T1(doc):
    para(doc,"**TABLE 1** Study groups by region, sex and diagnosis in the ROSMAP microglia atlas.",size=10,after=6)
    t=TT[0]; hdr=["Region","Sex","Diagnosis","Samples (n)","Donors (n)ᵃ","Age at death, years\n(mean ± SD)","PMI, h\n(mean ± SD)","Early / late AD\nsamples (n)"]
    body=[]; prev=None
    for _,r in t.iterrows():
        if r['Region']!=prev:
            body.append(('GROUP', {'HPC':'Hippocampus (86 samples; 66 identified donors)','PFC':'Prefrontal cortex (469 samples; 423 donors)'}[r['Region']])); prev=r['Region']
        body.append([r['Region'],r['Sex'],r['Diagnosis'],r['Samples (n)'],r['Donors (n)'],r['Age at death, years (mean ± SD)'],r['PMI, h (mean ± SD)'],r['Early / late AD samples (n)']])
    add_table(doc,[hdr],body,[1.6,1.6,2.0,2.2,2.2,3.4,2.6,3.0],size=8,align=['left','left','left','right','right','right','right','right'])
    para(doc,"Note: Samples and donors per group; age at death and PMI are donor-level means ± SD; early/late AD are per-sample counts from the source metadata [12]. "
         "Kruskal–Wallis tests across the four sex-by-diagnosis groups: HPC, age P = 0.606, PMI P = 0.788; PFC, age P = 6.95 × 10⁻⁴, PMI P = 0.552; PFC age, women with AD versus men with AD, two-sided Mann–Whitney P = 0.0142. "
         "ᵃIn HPC, 20 of 86 samples lack a subject identifier and are counted as distinct donors (66 identified donors in total). "
         "Numbers of nuclei per group were not retained with the processed results. "
         "Abbreviations: AD, Alzheimer's disease; HPC, hippocampus; PFC, prefrontal cortex; PMI, post-mortem interval; ROSMAP, Religious Orders Study and Rush Memory and Aging Project; SD, standard deviation.",size=8,before=4)
def contrast_hdr(first):
    h1=first+["HPC, men","","HPC, women","","PFC, men","","PFC, women",""]
    h2=['']*len(first)+["log₂FC","P_adj"]*4
    return h1,h2
FG_SHORT={'Alarmins/complement/immune receptors & signalling':'Alarmins, complement & immune receptors',
          'Lipid / proteoglycan / glycan metabolism':'Lipid, proteoglycan & glycan metabolism','Proteostasis / translation / DNA repair / other':'Proteostasis, translation & other',
          'Cytoskeleton / small-GTPase regulation / motility':'Cytoskeleton, GTPases & motility','Solute/amino-acid/ion transport':'Solute & ion transport',
          'Other transcriptional/chromatin regulators':'Other transcriptional regulators'}
def T2(doc):
    para(doc,"**TABLE 2** The cross-region core: the ten most strongly induced and all eight repressed core genes in the four sex-by-region contrasts.",size=10,after=6)
    t=TT[1]; first=["Rank","Gene","Functional group"]; last=["Summed log₂FC","Effect stratum","State marker","Feedback regulator"]
    h1,h2=contrast_hdr(first); h1+=last; h2+=['']*4
    body=[]; prev=None
    for _,r in t.iterrows():
        if r['Block']!=prev:
            body.append(('GROUP',r['Block'].replace('(8)','(n = 8)').replace('Ten most strongly induced core genes','Ten most strongly induced core genes (of 101 increased)'))); prev=r['Block']
        body.append([r['Rank (summed log₂FC)'],r['Gene'],FG_SHORT.get(r['Functional group'],r['Functional group'])]+[r[a+' '+b+' '+c] for a in ['HPC','PFC'] for b in ['men','women'] for c in ['log₂FC','P_adj']]+[r['Summed log₂FC'],r['Effect stratum (both regions)'],r['State marker (list)'] or '–','yes' if r['Feedback regulator']=='yes' else '–'])
    W=[1.0,1.5,3.6]+[1.3,1.9]*4+[1.5,1.3,1.3,1.5]
    al=['right','left','left']+['right']*9+['left','left','left']
    add_table(doc,[h1,h2],body,W,size=7,italic_cols=(1,),align=al,merges=[(0,3,4),(0,5,6),(0,7,8),(0,9,10)])
    para(doc,"Note: Genes are ranked by the sum of their four log₂FC values (ranking only). All 109 core genes have P_adj < 0.05 and the same direction of change in all four contrasts "
         "(101 increased, 8 decreased); 38 are high-effect (|log₂FC| > 1.5 in both sexes and both regions) and 71 low-effect (|log₂FC| < 0.5). State marker, membership of the curated marker list used for state annotation (Supplementary Table S1). "
         "Feedback regulator, recognized negative-feedback regulator of myeloid inflammatory signaling. The complete table is Supplementary Table S3. Abbreviations: "+ABBR,size=8,before=4)
GROUPS3=[("Chemokines, cytokines and receptors",['CCL3','CCL3L1','CCL2','IL1B','CXCR4']),("Alarmins and immune receptors",['S100A8','S100A9','CD83','CD163','THBD']),
         ("Lipid-handling and proteoglycan genes",['SRGN','ACSL1','PLTP']),("AP-1, EGR and NR4A transcription factors",['FOSL2','JUN','EGR1','EGR3','NR4A1','NR4A3']),
         ("Other negative-feedback regulators",['BCL6','BCL3','SOCS3','DUSP1','ZFP36L1']),("High-effect repressed core gene",['PLCXD3'])]
def T3(doc):
    para(doc,"**TABLE 3** Representative core genes, including the co-induced negative-feedback regulators, in the four sex-by-region contrasts.",size=10,after=6)
    t=TT[2].set_index('Gene'); h1,h2=contrast_hdr(["Gene","Effect stratum"]); h1+=["Note"]; h2+=['']
    body=[]
    for gname,genes in GROUPS3:
        body.append(('GROUP',gname))
        for g in genes:
            r=t.loc[g]; body.append([g,r['Effect stratum']]+[r[a+' '+b+' '+c] for a in ['HPC','PFC'] for b in ['men','women'] for c in ['log₂FC','P_adj']]+[r['Note'] or '–'])
    assert sum(len(g) for _,g in GROUPS3)==25 and set(sum([g for _,g in GROUPS3],[]))==set(t.index)
    W=[1.6,1.4]+[1.4,2.0]*4+[5.2]
    al=['left','left']+['right']*8+['left']
    add_table(doc,[h1,h2],body,W,size=7.5,italic_cols=(0,),align=al,merges=[(0,2,3),(0,4,5),(0,6,7),(0,8,9)])
    para(doc,"Note: Core genes discussed in the text, grouped by function. All have P_adj < 0.05 and the same direction in all four contrasts. Effect stratum: high, |log₂FC| > 1.5 in both sexes and both regions; low, |log₂FC| < 0.5. "
         "Marker status refers to the curated state-annotation lists (Supplementary Table S1; IL1B, CCL3: MG8; CD83, EGR1, NR4A1, NR4A3: MG10; JUN: MG6), so part of the fold change of these genes is expected by design. "
         "Feedback regulator, recognized negative-feedback regulator of myeloid inflammatory signaling (NR4A1, NR4A3, BCL6, BCL3, SOCS3, DUSP1, ZFP36L1). Functional groups and all values: Supplementary Table S3. Abbreviations: "+ABBR,size=8,before=4)
