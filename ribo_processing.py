import pysam as ps
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
from matplotlib import pyplot as plt
from scipy import stats, optimize

if __name__=='__main__':
    bamFile = ''
    center_weighting_from_sam(bamFile, min_read_length=20,max_read_length=42,reference='NC_000913')

def center_weighting_from_sam(bamFile,min_read_length=20,max_read_length=42,reference='NC_000913'):
    """Input: Bamfile of ribosome profile reads mapped to the reference genome
    Output1: Write .gff file of read densities at every nucleotide position on + and - strands in current directory
    Output2: Returns dataframe of read densities at every nucleotide position on + and - strands 
    Function: Center weighting is achieved by taking off 1/2 the min_read_length at each end of each read and weighting the 
    remaining nucleotides in the center equally. This is summed up for all reads at one position. Minimum and
    maximum read length allowed can be set (default of 20/42)"""
    
    bamfile = ps.Samfile( bamFile )
    base_name = bamFile[:-4]
    
    ### NEED TO GENERALIZE THIS TO ACCOUNT FOR MULTIPLE CHROMOSOMES...SEE ALI'S SCRIPT
    #Checks for reference in bamfile
    if not reference in bamfile.references:
            print 'Reference ID is invalid. These are your options:\n'
            for ref in bamfile.references: print ref+'\n'
            return
    for i,ref in enumerate(bamfile.references):
        if ref == reference:
            ref_id = i
    
    center_weight_plus=np.zeros(int(bamfile.lengths[0]))  #to make sure that all positions in the genome are accounted for, even those without reads
    center_weight_minus=np.zeros(int(bamfile.lengths[0]))
    for read in bamfile:
        if read.is_unmapped == True:
            continue
        if not read.rname == ref_id:
            continue
        if min_read_length > read.qlen or max_read_length < read.qlen:
            continue
        center_left = read.pos + min_read_length/2
        center_right = int(read.aend) - min_read_length/2
        center_length = center_right - center_left + 1
	#center weighting is done above by removing the min read length from each end and giving each position a weight inversely proportional to the resulting read length
        if read.is_reverse==False:
            center_weight_plus[center_left:center_right+1]+= (1./center_length)
        else:
            center_weight_minus[center_left:center_right+1]+= (1./center_length)
        #center weights for each of the reads are summed up at every nucleotide position for the plus and minus strands separately
    gff = open('%s_center_weighted.gff'%base_name,'w')
    for item in center_weight_plus:
        if center_weight_plus[item]!=0:
            gff.write('%s\tcenter_weighted\t%s\t%i\t%i\t%f\t+\t.\t.\n'%(reference,base_name,item,item,center_weight_plus[item]))
    for item in center_weight_minus:
        if center_weight_minus[item]!=0:
            gff.write('%s\tcenter_weighted\t%s\t%i\t%i\t%f\t-\t.\t.\n'%(reference,base_name,item,item,center_weight_minus[item]))
    gff.close()
    
    df= pd.DataFrame({'plus':center_weight_plus,'minus':center_weight_minus}, columns = ["plus", "minus"], index = range(1, len(center_weight_plus) + 1))
    df.to_csv('%s_center_weighted_df.csv'%base_name)
    
    
    
    return df
    
    
def geneFrame(genbank_file):
    """Input: Genbank file
    Output1: Dataframe containing all genes, start/stop locations and strand (plus/minus), 
    as well as gene function, product and amino acid sequence
    Output2: Dataframe where every nucleotide is a row in the sequence, used for... 
    Function: Parses genbank file to dataframe format for easier downstream utilization"""
    from Bio import SeqIO
    
    infile = SeqIO.read(genbank_file,'gb')
    genes = []
    name = []
    product = []
    func = []
    strand = []
    start = []
    stop = []
    aaseq = []
    cds_seq = []
    
    genome_seq_df = pd.DataFrame({'sequence':list(infile.seq.tostring())},index=range(1,len(infile.seq.tostring())+1))
    for feature in infile.features:
        if feature.type == 'CDS' and 'product' in feature.qualifiers:  #Only cares for coding sequences which are not pseudogenes
            genes.append(feature.qualifiers['locus_tag'][0])
            try: name.append(feature.qualifiers['gene'][0])
            except: name.append('')
            product.append(feature.qualifiers['product'][0])
            cds_seq.append(feature.location.extract(infile.seq).tostring())
            if 'function' in feature.qualifiers:                       #not all genes have known functions
                func.append(feature.qualifiers['function'][0])
            else:
                func.append("N/A")
            aaseq.append(feature.qualifiers['translation'][0])
            if feature.strand == 1:
                strand.append("plus") 
                start.append(feature.location.start.real+1)  
                stop.append(feature.location.end.real)
            elif feature.strand == -1:
                strand.append("minus") 
                start.append(feature.location.start.real+1)
                stop.append(feature.location.end.real)
    df = pd.DataFrame({"gene": genes, "name": name, "product": product, "function": func, "strand": strand, "start": start, "stop": stop, "cds_seq":cds_seq,"aaseq": aaseq},
                          columns = ["gene", "name", "function", "product", "strand", "start", "stop", "cds_seq","aaseq"])
    df = df.set_index("gene")
    return df, genome_seq_df  #FOR HAYTHEM: I'm not sure what the genome_seq_df is for

def countReads(rawreaddensities):
    """Input: Dataframe of center weighted read densities
    Output: Float count
    Function: counts the total number of reads, used for normalisation"""
    count = rawreaddensities.sum().sum()
    return count


def RPM_normed_gene_expression(genes, readdensities, totalreads, aa_ends_excluded=5):
    """Input1(genes): Dataframe produced by geneFrame()
    Input2(readdensities): Dataframe of center weighted read densities produced by center_weighting_from_sam()
    Input3(totalreads): Float from countReads() 
    Input4(aa_ends_excluded): Integer number of codons to exclude from each end to account for effects of translation initiation and termination. Default = 5
    Output: Dataframe of RPKM normalized read count for each gene
    Function: This script will calculate 'RPKM' for all genes in a given dataset"""
    expression = []
    ind_list=[]
    for index, row in genes.iterrows():
        start = row['start'] + aa_ends_excluded*3+1 #exclude first and last 5 codons to remove effects of translation initiation and termination
        stop = row['stop'] - aa_ends_excluded*3-1
        strand = row['strand']
        length = row['stop'] - row['start'] + 1
        genesum = readdensities[strand].ix[start: stop].sum() * 1000000000 / length / totalreads #sum up all read densities along gene
        expression.append(genesum)
        ind_list.append(index)
    gene_expression_df = pd.DataFrame({'gene_expression':expression},index=ind_list)
    return gene_expression_df

def absolute_synthesis_rate(gene_df, gene_expression_df, protein_mass_per_cell):  #HAYTHEM: Need help with this one, not sure where to find protein_mass_per_cell
    """Input1(gene_df): Dataframe from geneFrame()
    Input2(gene_expression_df): Dataframe of normalized gene read count from RPM_normed_gene_expression()
    Input3(protein_mass_per_cell): 
    Output: Dataframe of absolute synthesis rate for each gene
    Function: Calculates the absolute protein synthesis rate for each gene by multiplying its relative translation rate to the input protein mass per cell
    """
    from Bio.SeqUtils.ProtParam import ProtParamData, ProteinAnalysis
    #sum the total protein mass*readdensity
    total_protein_mass=0

    for index,row in gene_df.ix[:,['aaseq']].iterrows():
        if 'U' in row['aaseq']: #Checks for selenocystine AA and removes from analysis
            continue
        total_protein_mass += ProteinAnalysis(row['aaseq']).molecular_weight()*float(gene_expression_df.ix[index,'gene_expression'])
    #finds total protein mass of the cell based on translation rates and protein molecular weight of each gene

    #print index, ProteinAnalysis(row['aaseq']).molecular_weight()
    absolute_rate = {}
    for index,row in gene_expression_df.iterrows():
        absolute_rate[index]={'absolute_rate':float(row['gene_expression'])*protein_mass_per_cell*1e-15/total_protein_mass*6.0221413e23}
    #finds absolute synthesis rate of each protein as a fraction of the total protein mass
    absolute_rate_df = pd.DataFrame.from_dict(absolute_rate, orient='index')
    
    return absolute_rate_df

def meta_gene(genes, readdensities,distance=999):
    """
    Input1(genes): Dataframe from geneFrame()
    Input2(readdensities): Dataframe of center weighted read densities from center_weighting_from_SAM
    Input3(distance): Integer, how far along each gene to read. Default = 999
    Output1: Dataframe of averaged read densities of all front ends of genes(from start-30nt to start+distance)
    Output2: Dataframe of averaged read densities of all back ends of genes(from end-distance to end+30nt)
    """
    meta_gene_start_df=pd.DataFrame(index=genes.index,columns=range(-10,distance/3+1))
    meta_gene_stop_df=pd.DataFrame(index=genes.index,columns=range(-10,distance/3+1))
    
    for index, row in genes.iterrows():
        start = row['start']
        stop = row['stop']
        strand = row['strand']
        length = row['stop'] - row['start'] + 1
        if length < distance:
            continue #only count genes which are longer than the read distance or this will skew the averaging
        
        regionsum = float(readdensities[strand].ix[start-30: stop+30].sum())
        
        if strand=='plus':
            for n,i in enumerate(range(start-30,start+distance,3)):
                meta_gene_start_df[n-10][index]=readdensities[strand].ix[i:i+2].sum()/(regionsum/length)
            for n,i in enumerate(range(stop+30,stop-distance,-3)):
                meta_gene_stop_df[n-10][index]=readdensities[strand].ix[i-2:i].sum()/(regionsum/length)
            
        if strand=='minus':
            for n,i in enumerate(range(stop+30,stop-distance,-3)):
                meta_gene_start_df[n-10][index]=readdensities[strand].ix[i-2:i].sum()/(regionsum/length)
            for n,i in enumerate(range(start-30,start+distance,3)):
                meta_gene_stop_df[n-10][index]=readdensities[strand].ix[i:i+2].sum()/(regionsum/length)
        
    meta_gene_start_df = meta_gene_start_df.fillna(0)
    meta_gene_stop_df = meta_gene_stop_df.fillna(0)
    
    return meta_gene_start_df, meta_gene_stop_df
    
def half_gene_densities(genes, readdensities,min_read_count=1,aa_ends_excluded=5):
    """compares the read densities along the first half of the gene and the second half of the gene"""
    first = []
    second = []
    for index, row in genes.iterrows():
        start = row['start'] + aa_ends_excluded*3+1 #exclude first N codons to remove effects of translation initiation
        stop = row['stop'] - + aa_ends_excluded*3-1 #exclude last N codons to remove effects of translation termination
        strand = row['strand']
        length = stop - start
        firsthalf = 0.
        secondhalf = 0.
        if length > 60:   #Only count the gene if it is more than 20 codons
            firsthalf = readdensities[strand][start: start + length / 2].sum() #sum up first and second halves separately and append to the list
            secondhalf = readdensities[strand][start + length / 2: stop].sum()   #we don't count the first and last 5 codons
            if firsthalf > min_read_count and secondhalf > min_read_count:
                first.append(firsthalf / length * 3)
                second.append(secondhalf / length * 3)
    result = [first, second]
    #print stats.pearsonr(log(first),log(second))
    return result

def gene_length_dropoff_function(genes, readdensities,min_read_count=1,min_length=150,aa_ends_excluded=5,p0=[.446,6e-3,1]):
    """tries to check for dropoff of read densities along gene by comparing read densities in 
    50 codon windows against first 50 codons. This function takes a gene dataframe and a read
    density dataframe"""
    window_counter = {1:[]}     #all genes we look at will have first window
    window_averages = []
    for index, row in genes.iterrows():
        start = row['start'] + aa_ends_excluded*3+1   #exclude first and last 5 codons
        stop = row['stop'] - aa_ends_excluded*3-1
        strand = row['strand']
        length = stop - start
        if readdensities[strand][start:stop].sum() > min_read_count and length > min_length: #only take those with more than 128 reads
            window = 1          #Starting from second window
            #plus strand
            if strand == 'plus':
                firstwindow = readdensities[strand][start:start + min_length].sum()
                if firstwindow == 0.:
                    continue
                window_counter[window].append(1.)
                window += 1
                while start + min_length * (window) < stop:
                    if window in window_counter:
                        #append the read density for each window relative to the read density of the first window
                        window_counter[window].append(readdensities[strand][start + min_length * (window - 1): start + min_length * (window)].sum() / firstwindow)
                    else:
                        window_counter[window] = [readdensities[strand][start + min_length * (window - 1): start + min_length * (window)].sum() / firstwindow]
                    window += 1
                    if readdensities[strand][start + min_length * (window - 1): start + min_length * (window)].sum() / firstwindow > 49000.459866486657:
                        print index
            
            #minus strand
            else:
                firstwindow = readdensities[strand][stop - min_length:stop].sum()
                if firstwindow == 0.:
                    continue
                window_counter[window].append(1.)
                window += 1
                while stop - min_length * (window) > stop:
                    if window in window_counter:
                        #append the read density for each window relative to the read density of the first window
                        window_counter[window].append(readdensities[strand][stop - min_length * (window + 1): stop - min_length * (window)].sum() / firstwindow)
                    else:
                        window_counter[window] = [readdensities[strand][stop - min_length * (window - 1): stop - min_length * (window)].sum() / firstwindow]
                    window += 1
                    
    for i in range(1, len(window_counter) + 1):
        window_averages.append(np.median(window_counter[i]))
    
    plotx = range(25,875,50)
    plot_x = np.array(plotx)
    popt, pcov = optimize.curve_fit(model_func,plot_x,window_averages[0:len(plotx)],p0=p0)
    A, K, C = popt[0], popt[1], popt[2]
    residuals = window_averages[0:len(plotx)] - model_func(plot_x,A,K,C)
    sqr_res = sum(residuals**2)
    print 'The calculated function parameters are: \nA=%f\nK=%f\nC=%f\nSum of residuals squared=%f'%(A,K,C,sqr_res)
    plt.scatter(plotx,window_averages[0:len(plotx)])
    plt.plot(plot_x,model_func(plot_x,A,K,C))
    return A,K,C


    
def model_func(t, A, K, C):
    """Function: Exponential dropoff fitting function"""
    return A * np.exp(-K * t) + C


def exp_dropoff_correction(genes, readdensities,min_read_count=1,min_length=150,aa_ends_excluded=5,p0=[.446,6e-3,1]):
    
    A,K,C = gene_length_dropoff_function(genes, readdensities,min_read_count=min_read_count,\
                                        min_length=min_length,aa_ends_excluded=aa_ends_excluded,p0=p0)
    i = range(1,10000)
    multiplier = []
    for nt in i:
        j = model_func(float(nt)/3, A, K, C)
        multiplier.append(j)
    
    dropoff_corrected_read_densities = readdensities.copy()
    dropoff_corrected_read_densities['plus_corrected'] = 0.000
    dropoff_corrected_read_densities['minus_corrected'] = 0.000    #make a copy of the read densities and set the initial values for correction columns to 0.
    
    for index, row in genes.iterrows():
        #print index
        start = row['start'] + aa_ends_excluded*3+1   #exclude first and last 5 codons
        stop = row['stop'] - aa_ends_excluded*3-1
        strand = row['strand']
        length = stop - start
        
        if strand == 'plus':
            for i in range(length + 1):
                dropoff_corrected_read_densities['plus_corrected'][start + i] = dropoff_corrected_read_densities['plus'][start + i] / multiplier[i]
        
        else:
            for i in range(length + 1):
                dropoff_corrected_read_densities['minus_corrected'][stop - i] = dropoff_corrected_read_densities['minus'][stop - i] / multiplier[i]
                
    del dropoff_corrected_read_densities['plus']
    del dropoff_corrected_read_densities['minus']
    dropoff_corrected_read_densities.columns = ['plus','minus']
    return dropoff_corrected_read_densities
