
def T4(doc):
    para(doc,"**TABLE 4** Sex-specific and sex-dimorphic DEGs in the hippocampus and prefrontal cortex.",size=10,after=6)
    t=TT[3]
    h1=["Gene","Men","","Women","","Category and log₂FC (men / women) in the other region"]; h2=["","log₂FC","P_adj","log₂FC","P_adj",""]
    order=[('HPC','male-specific','HPC, male-specific (n = 9)'),('HPC','sex-dimorphic','HPC, sex-dimorphic (n = 4)'),
           ('PFC','male-specific','PFC, male-specific (n = 5)'),('PFC','female-specific','PFC, female-specific (n = 2)'),('PFC','sex-dimorphic','PFC, sex-dimorphic (n = 1)')]
    body=[('GROUP','(A) Per-gene statistics')]
    for reg,cat,lab in order:
        sub=t[(t.Region==reg)&(t.Category==cat)]
        assert len(sub)==int(re.search(r'n = (\d+)',lab).group(1)), lab
        body.append(('GROUP',lab))
        for _,r in sub.iterrows():
            body.append([r['Gene'],r['Men log₂FC'],r['Men P_adj'],r['Women log₂FC'],r['Women P_adj'],r['Status in the other region']])
    b=TT[4].iloc[0]
    body.append(('GROUP','(B) HPC, female-specific (n = 45; all increased in women)'))
    W=[2.2,1.8,2.3,1.8,2.3,9.0]
    tb=add_table(doc,[h1,h2],body+[['']*6],W,size=7.5,italic_cols=(0,),align=['left','right','right','right','right','left'],merges=[(0,1,2),(0,3,4)])
    last=tb.rows[-1]; a=last.cells[0].merge(last.cells[5])
    genes=b['Gene symbols (ordered by P_adj in women)'].replace('↑','')
    glist=[g.strip() for g in genes.split(',')]; assert len(glist)==45
    write(a,"Ordered by P_adj in women: "+", ".join("*"+g+"*" for g in glist)+
          ". P_adj (women) "+b['Range of P_adj (women)']+"; same direction in men, "+b['Same direction in men (n)']+" (binomial P = 9.3 × 10⁻⁹); median (minimum) women:men |log₂FC| ratio, "+b[RATIO_COL]+"; sex-shared in PFC, "+b['Sex-shared in the other region (n)']+" (11 in the filtered and 6 in the unfiltered PFC list).",size=7.5)
    para(doc,"Note: Sex-specific, P_adj < 0.05 and |log₂FC| > 0.25 in one sex and nominal P > 0.5 in the other; sex-dimorphic, P_adj < 0.05 and |log₂FC| > 0.25 in both sexes with opposite directions. "
         "P_adj = 1 corresponds to nominal P > 0.5 in the non-significant sex. In HPC, 7 of the 9 male-specific genes were decreased in men and 8 of 9 changed in the same direction in women (binomial P = 0.039; median men:women |log₂FC| ratio 2.75). "
         "Last column: category and men / women log₂FC of the gene in the other region (filtered list unless stated); Absent, not in the filtered or unfiltered list of the other region. "
         "Full per-region tables: Supplementary Table S9; all values: Supplementary Table S2. Abbreviations: "+ABBR,size=8,before=4)
