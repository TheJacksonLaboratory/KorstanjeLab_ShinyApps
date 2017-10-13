# Copying required data for shiny app
# A good practice for typing your passowrd a 100 times without making mistakes,
# otherwise IT will lock yo ufor for 30 mins!


ssh cadillac

# DOQTL objects and genoprob data
cd /hpcdata/ytakemon/Col4a5xDO/GBRS_reconstruction/reconstruct/best.compiled.genoprob

# Move qtl permuations
list=(Alb10.F.perm Alb15.F.perm Alb6.F.perm GFR.F.perm Alb10.Hi.perm Alb15.Hi.perm Alb6.Hi.perm GFR.Hi.perm Alb10.M.perm Alb15.M.perm Alb6.M.perm GFR.M.perm Alb10.perm Alb15.perm Alb6.perm GFR.sex.int.perm Alb10.sex.int.perm Alb15.sex.int.perm Alb6.sex.int.perm)
for i in ${list[*]}
do
scp -r qtl.perm/${i}/1000.perms.txt ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/qtl.perm/${i}/
done

# Move gwas permutations
list=(ACR10.perm ACR15.perm ACR6.perm Alb10.perm Alb15.perm Alb6.perm gfr.c2.perm)
for i in ${list[*]}
do
scp -r gwas.perm/${i}/all.perm.txt ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/gwas.perm/${i}/
done

# Ensemble data
scp -r EnsemblID_GRCm38.p4.Rdata ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/
# GM_snps
scp -r GM_snps.Rdata ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/
# Sex covar
scp -r sex.covar.Rdata ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/
# Genoprobs
scp -r genoprobs/best.genoprobs.192.Rdata ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/genoprobs/
# Gwas R objects
scp -r gwas ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/
# Kinship r obj
scp -r kinship ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/
# qtl perm objects
scp -r perm1000 ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/
# qtl plots
scp -r plot ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/
# move qtl objects
scp -r qtl ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/
# move Gene allele
scp -r Gene_allele.Rdata ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/
# move RNA seq data directory
scp -r RNA_seq_Rdata ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO/best.compiled.genoprob/



# Other data and resouces used
cd /hpcdata/ytakemon/Col4a5xDO/
scp -r Phenotype ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO
scp -r resources ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO
scp -r Sample_list ytakemon@helix.jax.org:/projects/ytakemon/Col4a5xDO

# Need to copy into shinyapp directly
cp /projects/ytakemon/Col4a5xDO/best.compiled.genoprob/plot/RNA_qtl/Complete_eQTL_plots/* /home/ytakemon/ShinyApps/Col4a5xDO/eQTL/www/
