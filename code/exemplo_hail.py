# pip install hail
import hail as hl

arquivo = '/home/matheus/Documents/gnomeAD-tb/data/gnomad.exomes.v4.1.sites.chr18.vcf.bgz'

hl.init()


# Ler VCF
mt = hl.import_vcf(arquivo, force_bgz=True, reference_genome='GRCh38')

# Salvar no formato de matrix do hail
mt.write('meu_arquivo_hail_matrix')

# Printar infos
print(mt.count())
print(mt.summarize())

