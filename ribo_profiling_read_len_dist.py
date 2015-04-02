import pysam, numpy, sys
import matplotlib
matplotlib.use('Agg') # do this before importing pyplot
from matplotlib import pyplot
pyplot.ioff()

ribo_profile_trna = pysam.Samfile(sys.argv[1],'r')
ribo_profile_rrna = pysam.Samfile(sys.argv[2],'r')
ribo_profile_other_rna = pysam.Samfile(sys.argv[3],'r')
ribo_profile_no_rna = pysam.Samfile(sys.argv[4],'r')
ribo_profile_all = pysam.Samfile(sys.argv[5],'r')

iter_trna = ribo_profile_trna.fetch()
iter_rrna = ribo_profile_rrna.fetch()
iter_other_rna = ribo_profile_other_rna.fetch()
iter_no_rna = ribo_profile_no_rna.fetch()
iter_all = ribo_profile_all.fetch()

def fix_name(name):
    return name.split("_")[0]

lens_no_rna=[]
for read in iter_no_rna:
    if read.is_unmapped==False:
        lens_no_rna.append(read.alen)

lens_trna=[]
for read in iter_trna:
    if read.is_unmapped==False:
        lens_trna.append(read.alen)

lens_rrna=[]
for read in iter_rrna:
    if read.is_unmapped==False:
        lens_rrna.append(read.alen)

lens_other_rna=[]
for read in iter_other_rna:
    if read.is_unmapped==False:
        lens_other_rna.append(read.alen)

lens_all=[]
for read in iter_all:
    if read.is_unmapped==False:
        lens_all.append(read.alen)
     
ribo_profile_trna.close()
ribo_profile_rrna.close()
ribo_profile_other_rna.close()
ribo_profile_no_rna.close()
ribo_profile_all.close()

hist_all = numpy.histogram(lens_all, bins=70,range=(1,71))
hist_no_rna = numpy.histogram(lens_no_rna, bins=70,range=(1,71))
hist_trna = numpy.histogram(lens_trna, bins=70, range=(1,71))
hist_rrna = numpy.histogram(lens_rrna, bins=70, range=(1,71))
hist_other_rna = numpy.histogram(lens_other_rna, bins=70, range=(1,71))

pyplot.figure(figsize=(11,8.5))

pyplot.subplot(411)
pyplot.title('%s: Summary'%fix_name(sys.argv[4]))
#pyplot.bar(hist_all[1][1:],hist_all[0],color='blue',label='All')
pyplot.bar(hist_no_rna[1][1:],hist_no_rna[0],color='red',label='mRNA',linewidth=0)
pyplot.bar(hist_trna[1][1:],hist_trna[0],color='green',label='tRNA', bottom=hist_no_rna[0],linewidth=0)
pyplot.bar(hist_rrna[1][1:],hist_rrna[0],color='yellow',label='rRNA', bottom=hist_no_rna[0]+hist_trna[0],linewidth=0)
pyplot.bar(hist_other_rna[1][1:],hist_other_rna[0],color='purple',label='other RNA', bottom=hist_no_rna[0]+hist_trna[0]+hist_rrna[0],linewidth=0)
pyplot.xlabel('Read lengths')
pyplot.ylabel('Counts')
pyplot.legend()

pyplot.subplot(412)
pyplot.title('%s: tRNA'%fix_name(sys.argv[1]))
pyplot.bar(hist_all[1][1:],hist_all[0],color='blue',label='All',linewidth=0)
pyplot.bar(hist_trna[1][1:],hist_trna[0],color='green',label='tRNA',linewidth=0)
pyplot.xlabel('Read lengths')
pyplot.ylabel('Counts')
pyplot.legend()

pyplot.subplot(413)
pyplot.title('%s: rRNA'%fix_name(sys.argv[2]))
pyplot.bar(hist_all[1][1:],hist_all[0],color='blue',label='All',linewidth=0)
pyplot.bar(hist_rrna[1][1:],hist_rrna[0],color='yellow',label='rRNA',linewidth=0)
pyplot.xlabel('Read lengths')
pyplot.ylabel('Counts')
pyplot.legend()

pyplot.subplot(414)
pyplot.title('%s: other-RNA'%fix_name(sys.argv[3]))
pyplot.bar(hist_all[1][1:],hist_all[0],color='blue',label='All',linewidth=0)
pyplot.bar(hist_other_rna[1][1:],hist_other_rna[0],color='purple',label='other RNA',linewidth=0)
pyplot.xlabel('Read lengths')
pyplot.ylabel('Counts')
pyplot.legend()

pyplot.tight_layout()
pyplot.savefig('%s_plots.pdf'%(sys.argv[1][:-24]))

