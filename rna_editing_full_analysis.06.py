#!venv/bin/python3
# coding: utf-8

#check presence of virtual env
try:
    ###libraries
    import os
    import sys
    import gzip
    import subprocess
    import shutil
    import numpy as np
    import pysam
    import collections
    import math
    import time
    import string
    import multiprocessing
    from multiprocessing import Process,Queue
    import colorsys
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import datetime
    import random
    import gc
    import argparse
    from argparse import RawTextHelpFormatter
except Exception as e:
    print("Error when loading libraries.")
    print(e)
    sys.exit()

###user arguments
helper=['UMI extraction','Adapter clipping','Mapping Ref 1','Read mapping 1','Alignment sorting','Indexing 1','Deduplication','PCR template_reconstruction','SNP ident.','focus analysis']

parser=argparse.ArgumentParser(description='input table (-f): DNA/RNA \\t sample ID \\t fwd/rev \\t path_to_file, no header\n\nlist for omitting steps:\n%s'%('\n'.join([f'{e+1} - {h}' for e,h in enumerate(helper)])),formatter_class=RawTextHelpFormatter)
parser.add_argument('--version', action='version', version='6.0')
parser.add_argument('-f',dest='infile',default='',help='list of sequencing read input files [""]')
parser.add_argument('-r',dest='ref',default='',help='mapping ref fasta [""]')
parser.add_argument('-a',dest='annotation',default='',help='path to reference annotation file [""]')
parser.add_argument('-n',dest='cores',default=1,help='number of CPUs [1]')
parser.add_argument('-o',dest='outname',default='',help='output filename prefix [""]')
parser.add_argument('--cov',dest='cov',default=5,help='min coverage [5]')
parser.add_argument('--direction',dest='direction',default=1,help='min reads per orientation for support [1]')
parser.add_argument('--nomap',dest='no_mapping',default=False,action='store_true',help='skip mapping step b/c data exists already [False]')
parser.add_argument('-x',dest='phix',default='',help='PhiX fasta sequence [""]')
parser.add_argument('-l',dest='log',default='pipeline.log',help='log file ["pipeline.log"]')
parser.add_argument('-d',dest='devel',default=False,action='store_true',help='developer mode [OFF]')
parser.add_argument('-s',dest='skip',default='1'*len(helper),help=F'omit analysis steps (binary) ["{"1"*len(helper)}"]')
parser.add_argument('--att',dest='attrition',default=True,action='store_false',help='create read attrition plot [True]')
parser.add_argument('--no_umi',dest='no_umi',default=False,action='store_true',help='omit umis and use regular deduplication')
parser.add_argument('--term',dest='terminate',default=0,type=int,help='terminate analysis prematurely after this processing step [0=none]')
parser.add_argument('--askip',dest='auto_skip',default=False,action='store_true',help='PARTIALLY IMPLEMENTED: skip existing files except latest one')
options=parser.parse_args()

###defaults
metrics={'delete_files':[],'programs':{},'site_QC_criteria':{},'references':{},'log':''}
if not os.path.exists('config.txt'):
    print('Please create the required config file "config.txt" in the same directory as this script.')
    print('Alternatively, you can run the "Installer.sh" script to create the config file.')
    sys.exit('Terminating.')
else:
    with open('config.txt','r') as infile:
        for line in infile:
            if not line.startswith('#') and len(line)>1:
                line=line.strip('\r\n').strip('\n').split(' = ')
                metrics['programs'][line[0]]=line[1][1:-1]
metrics['references']['cutadapt']={'fwd':['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC','TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC'],'rev':['AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC','GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT']}
#metrics setup
metrics['site_QC_criteria']['min_site_qs']=30
metrics['site_QC_criteria']['max_base_count_per_site']=2
metrics['site_QC_criteria']['min_alt_freq']=0.01# fraction, not percent
metrics['site_QC_criteria']['max_homopol_len']=4
metrics['site_QC_criteria']['min_dist_to_paralogs']=50
metrics['site_QC_criteria']['min_dist_from_read_end']=4#bp
#additional QC metrics
metrics['site_QC_criteria']['DNA_max_cov_diff_SD']=2.0
metrics['site_QC_criteria']['DNA_cov_win_size']=100
metrics['skip_steps']=helper

###functions
def create_attrition_plot(attrition,outfile,width=6,height=4):
    sampleIDs=[]
    for seq_typ in ['DNA','RNA']:
        sampleIDs+=[seq_typ+'_'+x for x in sorted(list(attrition[seq_typ].keys()))]
    colors=[[.2,.2,.2]]+[[x for x in colorsys.hsv_to_rgb((i+1)*1.0/(len(sampleIDs)+1),.95,.95)] for i in range(len(sampleIDs)-1)]
    colors={sampleIDs[i]:colors[i] for i in range(len(sampleIDs))}

    # print(attrition)

    fig,ax1=plt.subplots(figsize=(width,height), dpi=300)
    maxy=0
    for seq_typ in ['DNA','RNA']:
        for sID in attrition[seq_typ]:
            if seq_typ=='RNA': linestyle=':'
            else: linestyle='-'
            plt.plot(list(range(len(attrition['stages']))),attrition[seq_typ][sID],color=colors[seq_typ+'_'+sID],lw=1,linestyle=linestyle,label=seq_typ+'_'+sID)
            if max(attrition[seq_typ][sID])>maxy: maxy=max(attrition[seq_typ][sID])
    ax1.set_xlabel('data processing stage')
    ax1.set_ylabel('read counts')
    ax1.set_ylim(0,maxy*1.5)
    ax1.set_xticks(list(range(len(attrition['stages']))))
    ax1.set_xticklabels(attrition['stages'],ha='center',size='xx-small',rotation=90)
    plt.legend(fontsize='x-small',loc='upper center', ncol=4,framealpha=0.75, scatterpoints=1)
    plt.tight_layout()

    if not outfile.endswith('.svg'): outfile+='.svg'
    plt.savefig(outfile,format='svg')
    plt.close()
    print(F'Figure file {outfile} has been created.')

    with open(outfile.rsplit('.',1)[0]+'.txt','w') as outstream:
        outstream.write('\t'.join(['seq_typ','sID']+attrition['stages'])+'\n')
        for seq_typ in ['DNA','RNA']:
            for sID in attrition[seq_typ]:
                outstream.write('\t'.join([seq_typ,sID]+[str(x) for x in attrition[seq_typ][sID]])+'\n')
    print(F"Figure data file {outfile.rsplit('.',1)[0]+'.txt'} has been written.")

    return(outfile)

def analyze_coverages(arr_cov,annotation):
    results={}
    length=arr_cov.size

    gene_cov=[]

    pos_arr=np.zeros(shape=(length),dtype='u1')

    missing=[]

    for start,end,strand,description in annotation:
        avg_cov=sum(arr_cov[start:end+1])*1.0/(end-start)
        gene_cov.append(avg_cov)
        pos_arr[start:end+1]=1
        if avg_cov==0: missing.append((start,end,strand,description))

    intronic=average([arr_cov[e] for e in range(arr_cov.size) if pos_arr[e]==0])
    gene_cov=average(gene_cov)
    overall_cov=average(arr_cov)

    return(overall_cov,gene_cov,intronic,missing)

def read_input_table(datei):#returns {DNA/RNA: list with ID,direction,file path}
    results={'DNA':[],'RNA':[]}
    matching_types={}
    with open(datei,'r') as infile:
        for line in infile:
            if not line.startswith('#') and len(line)>1:
                if '\t' in line:
                    line=line.strip('\r\n').strip('\n').split('\t')
                else:
                    line=line.strip('\r\n').strip('\n').split()
                if len(line)!=4:
                    sys.exit(f'Invalid number of columns, expecting 4 but seeing {len(line)}.\nLine:\n{line}\n Aborting.')
                else:
                    if not line[0].upper() in ['DNA','RNA']:
                        sys.exit('Missing DNA/RNA information in column 1. Aborting.')
                    else:
                        if not os.path.exists(line[3]):
                            sys.exit('Incorrect file location for: %s. Aborting.'%(' > '.join(line)))
                        else:
                            results[line[0].upper()].append(line[1:])
                            results[line[0].upper()][-1][1]=results[line[0].upper()][-1][1].lower()
                            if not results[line[0].upper()][-1][1] in ['fwd','rev']:
                                print('Incorrect column information (direction) for %s. Only "fwd" or "rev" are allowed. Omitting file.')
                                results[line[0].upper()].pop(-1)
                            else:
                                if not line[1] in matching_types: matching_types[line[1]]=[False,False]#DNA,RNA
                                if line[0].upper()=='DNA': matching_types[line[1]][0]=True
                                elif line[0].upper()=='RNA': matching_types[line[1]][1]=True
    for sampleID,liste in list(matching_types.items()):
        if liste[0]*liste[1]==True: matching_types.pop(sampleID)

    if len(matching_types)>0:
        print('The following samples are missing matching DNA and RNA data and will be removed from analysis: %s'%(', '.join(matching_types.keys())))
        for typ in results.keys():
            results[typ]=[l for l in results[typ] if not l[1] in matching_types]

    if len(results['DNA'])==0:
        sys.exit('No input sample remain with matching DNA and RNA datasets.')

    return(results)

def read_fasta(datei):
    if datei.endswith('gz'):
        infile=gzip.open(datei,'rt')
    else:
        infile=open(datei,'r')

    seq=''
    title=''
    for line in infile:
        line=line.strip()
        if line.startswith('>') and title=='': title=line[1:]
        elif line.startswith('>') and title!='': break
        else: seq+=line

    infile.close()
    return(seq,title)

def parse_parameters(options,metrics):
    if options.devel:
        options.infile=''
        options.ref=''
        options.annotation=''
        metrics['CPU_count']=50
        options.outname='dev_test'

    try:
        metrics['CPU_count']=int(options.cores)
        if metrics['CPU_count']<1:
            print('Invalid selection of number of CPUs. Setting to 1.')
            metrics['CPU_count']=1
    except:
        print('Invalid selection of number of CPUs. Setting to 1.')
        metrics['CPU_count']=1

    metrics['skip']=options.skip
    if len(metrics['skip'].replace('1','').replace('0',''))>0:
        sys.exit('Incorrect format for analysis step omission. Please use binary string format only.')
    elif len(metrics['skip'])<len(metrics['skip_steps']):
        while len(metrics['skip'])<len(metrics['skip_steps']): metrics['skip']+='1'
        # sys.exit(f"Please provide the correct number of steps for omission [{len(metrics['skip_steps'])}]. Aborting.")

    if options.phix=='':
        sys.exit('Missing PhiX fasta file. Aborting.')
    elif not os.path.exists(options.phix):
        sys.exit('PhiX fasta file "%s" does not exist. Aborting.'%options.phix)
    else:
        metrics['PhiX']=options.phix

    if options.infile=='':
        sys.exit('Missing list of read files. Aborting.')
    elif not os.path.exists(options.infile):
        sys.exit('Input list file "%s" does not exist. Aborting.'%options.infile)
    else:
        metrics['input_files']=read_input_table(options.infile)

    if options.ref=='':
        sys.exit('Missing reference sequence fasta file. Aborting.')
    elif not os.path.exists(options.ref):
        sys.exit('Reference sequence file "%s" does not exist. Aborting.'%options.ref)
    else:
        metrics['ref_file']=options.ref
        metrics['ref_seq'],metrics['ref_ID']=read_fasta(metrics['ref_file'])
        metrics['ref_ID']=metrics['ref_ID'].split(' ')[0]

    if options.annotation=='':
        sys.exit('Missing reference annotation file. Aborting.')
    elif not os.path.exists(options.annotation):
        sys.exit('Reference annotation file "%s" does not exist. Aborting.'%options.annotation)
    else:
        metrics['ref_anno_file']=options.annotation

    if options.outname=='':
        metrics['prefix']=''
    else:
        metrics['prefix']=options.outname
        if not metrics['prefix'][-1]=='.': metrics['prefix']+='.'

    try:
        metrics['site_QC_criteria']['min_read_count_per_direction']=int(options.direction)
    except:
        sys.exit('Unable to parse direction flag with value "%s". Please use integers only. Aborting'%options.direction)

    try:
        metrics['site_QC_criteria']['min_cov']=int(options.cov)
    except:
        sys.exit('Unable to parse coverage flag with value "%s". Please use integers only. Aborting'%options.cov)

    metrics['skip_mapping']=options.no_mapping

    if options.log=='':
        metrics['log']='pipeline.log'
    else:
        metrics['log']=options.log

    # print(F"Number of cores: {metrics['CPU_count']}")

    return metrics

def calc_SD(liste):
  if len(liste)>1:
    avg=average(liste)
    tmp=sum((i-avg)**2 for i in liste)
    if tmp!=0: return (tmp/(len(liste)-1))**.5
    else: return 0.0
  else:
    return 'n/a'

def average(liste):#returns average value
  if len(liste)>0:
    return sum(liste)*1.0/len(liste)
  else:
    return 'n/a'

def load_coverage_pysam(bam_file, chrom=None, pad=True):
    """
    Calculate coverage statistics for a given BAM file using pysam library.

    Parameters:
        bam_file (str): Path to the input BAM file.
        pad (bool, optional): If True, include padding for positions with zero coverage. 
                              If False, exclude padding. Default is True.
        chrom (str, optional): The reference ID (name) for which coverage should be calculated.
                                If None (default), coverage will be calculated for all references.

    Returns:
        dict: A dictionary containing coverage statistics.
              The keys are reference names and the values are lists of integers representing coverage at each position.
    """
    # Dictionary to store coverage for each reference
    coverage_dict = {}

    # Open the BAM file for reading
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Get a list of reference names
        ref_names = bam.references

        # If chrom is provided, calculate coverage for only that reference
        if chrom is not None:
            if chrom not in ref_names:
                raise ValueError("Reference ID not found in the BAM file.")
            ref_names = [chrom]

        # Iterate through each reference
        for ref_name in ref_names:
            # Get the length of the reference
            ref_length = bam.get_reference_length(ref_name)
            # Initialize coverage list for this reference
            coverage_list = [0] * ref_length

            # Iterate through pileup data for this reference
            for pileup_column in bam.pileup(ref_name, stepper="samtools", pad=True, ignore_overlaps=True):
                # Get the 0-based position in the reference
                pos = pileup_column.reference_pos

                # Calculate the coverage at this position (considering gaps)
                coverage = pileup_column.nsegments
                coverage_list[pos] = coverage

            # Store the coverage list for this reference in the dictionary
            coverage_dict[ref_name] = coverage_list

    if chrom is not None:
        return coverage_dict[chrom]
    else:
        return coverage_dict

def qualify_mapping_results(bamfile,ref_ID,typ,cores):
    bamfile_pointer=pysam.AlignmentFile(bamfile,'rb',threads=cores)
    ref_index=bamfile_pointer.references.index(ref_ID)
    ref_len=bamfile_pointer.get_reference_length(ref_ID)
    #arr contains all read starts and ends
    arr_ends=np.zeros(shape=(ref_len),dtype='u4')
    count=0
    avg_read_length=0
    for read in bamfile_pointer.fetch():

        if not read.is_duplicate and not read.is_unmapped and read.reference_id==ref_index:
            count+=1
            arr_ends[read.reference_start]+=1
            arr_ends[read.reference_end-1]+=1
            if count%100000==0: print('Currently processed read count: {:,}'.format(count))
            avg_read_length+=read.query_length

    avg_read_length=int(round(avg_read_length*1.0/count))

    #array with coverage
    arr_cov=load_coverage_pysam(bamfile, chrom=ref_ID, pad=True)
    arr_cov=np.array(arr_cov,dtype='u4')


    critical_ranges_and_sites=ident_outlier(arr_ends,arr_cov,avg_read_length)

    print('%i potential critial structural aberrations identified.'%len(critical_ranges_and_sites))
    # print(critical_ranges_and_sites)

    return(arr_cov,critical_ranges_and_sites)

def ident_outlier(arr_ends,arr_cov,read_length,SD_factor=2.0):
    ctr=True
    threshold=max(arr_ends)+1
    excluded=0
    while ctr:
        ends=[v for v in arr_ends if v>0.0 and v<threshold]
        hist=collections.Counter(ends)
        a,b=linear_fit(hist.items())
        a,b=linear_fit([(k,math.log(v,2)) for k,v in hist.items()])
        threshold_new=math.floor(-b*1.0/a)
        # threshold_new=-b*1.0/a
        excluded_new=len([v for v in ends if v>=threshold_new])
        if excluded_new>0:
            threshold=threshold_new
        elif threshold_new<=2: ctr=False
        else: ctr=False

    SD=[arr_ends[p]*1.0/arr_cov[p] for p in range(len(arr_cov)) if arr_cov[p]>2]
    avg_ratio=average(SD)
    SD=calc_SD(SD)
    ratio_threshold=avg_ratio+SD_factor*SD

    critical_sites=[]
    for pos in range(len(arr_ends)):
        if arr_ends[pos]>threshold:
            if arr_ends[pos]*1.0/arr_cov[pos]>ratio_threshold:
                critical_sites.append((pos,arr_ends[pos]*100.0/arr_cov[pos]))

    filtered=[]
    for pos in range(len(critical_sites)-1):
        if critical_sites[pos+1][0]-critical_sites[pos][0]>read_length: filtered.append((critical_sites[pos][0],critical_sites[pos+1][0],critical_sites[pos][1],critical_sites[pos+1][1]))

    filtered+=critical_sites
    return filtered#returns list with potentially problematic regions, followed by all individual sites

def linear_fit(dataset_pairs):#returns f(x)=ax+b
  mean_x=sum([x for x,y in dataset_pairs])*1.0/len(dataset_pairs)
  mean_y=sum([y for x,y in dataset_pairs])*1.0/len(dataset_pairs)
  sum_ss_xy=sum((x-mean_x)*(y-mean_y) for x,y in dataset_pairs)
  sum_ss_xx=sum((x-mean_x)**2 for x,y in dataset_pairs)
  if sum_ss_xx!=0: a=sum_ss_xy/sum_ss_xx
  else: a=0
  b=mean_y-a*mean_x
  return [a,b]

def mask_genome(metrics):#1D array w/ ref len, 1==masked
    ref_len=len(metrics['ref_seq'])
    mask=np.zeros(shape=(ref_len),dtype='u1')

    #homopol exclusion
    batch_size=int(ref_len/metrics['CPU_count'])
    q=Queue()
    thread_list=[]
    for i in range(metrics['CPU_count']-1):
        thread_list.append((i*batch_size,Process(target=homopol_worker,args=(metrics['ref_seq'][i*batch_size:(i+1)*batch_size],i*batch_size,metrics['site_QC_criteria']['max_homopol_len'],q))))
        thread_list[-1][1].deamon=True
        thread_list[-1][1].start()
    thread_list.append(((metrics['CPU_count']-1)*batch_size,Process(target=homopol_worker,args=(metrics['ref_seq'][(metrics['CPU_count']-1)*batch_size:],(metrics['CPU_count']-1)*batch_size,metrics['site_QC_criteria']['max_homopol_len'],q))))
    thread_list[-1][1].deamon=True
    thread_list[-1][1].start()

    while len(thread_list)>0:
        if not q.empty():
            ID,res=q.get()
            for pos in res: mask[pos]=1
            thread_list.pop([t[0] for t in thread_list].index(ID))

    return(mask)

def homopol_worker(seq,offset,threshold,q):
    pos_list=[]
    old_base=''
    length=1
    start=0
    while seq[start]==seq[0]: start+=1
    for pos in range(start,len(seq)):
        if seq[pos]==old_base: length+=1
        else:
            if length>threshold:
                for j in range(length):
                    pos_list.append(offset+j-(pos-length))
            length=1
            old_base=seq[pos]

    q.put((offset,pos_list))
    return(True)

def perform_SNP_ident(bamfileIDs,ref_ID,metrics):
    snp_list=[]
    excl_snp_list=[]
    q=Queue()
    thread_list=[]
    batch_size=int(len(metrics['ref_seq'])/metrics['CPU_count'])

    for i in range(metrics['CPU_count']-1):
        thread_list.append((i,Process(target=snp_worker2,args=(i,bamfileIDs,ref_ID,i*batch_size,(i+1)*batch_size,metrics,q))))
        thread_list[-1][1].deamon=True
        thread_list[-1][1].start()
    thread_list.append((metrics['CPU_count']-1,Process(target=snp_worker2,args=(metrics['CPU_count']-1,bamfileIDs,ref_ID,(metrics['CPU_count']-1)*batch_size,len(metrics['ref_seq']),metrics,q))))
    thread_list[-1][1].deamon=True
    thread_list[-1][1].start()

    while len(thread_list)>0:
        if not q.empty():
            ID,res,excl_res=q.get()
            snp_list+=res
            excl_snp_list+=excl_res
            pro_index=[t[0] for t in thread_list].index(ID)
            thread_list[pro_index][1].join()
            thread_list.pop(pro_index)
            # print(ID,ID*batch_size,(ID+1)*batch_size,res[0],res[-1])

    return(snp_list,excl_snp_list)

def snp_worker2(ID,bamfiles,ref_ID,start,end,metrics,q):#pileup implementation
    arr=[]
    for bamfile in bamfiles:
        bamfile_pointer=pysam.Samfile(bamfile,'rb')
        total_sites=end+1-start
        column_steps=int(round(total_sites*.05))
        column_count=0
        for pileupcolumn in bamfile_pointer.pileup(ref_ID,start,end+1,truncate=True,ignore_overlaps=False,ignore_orphans=False):
            column_count+=1

            counts={b:[0,0] for b in ['A','C','G','T','N','*']}
            templates={}

            if ID==0 and column_count%column_steps==0: print(f'Process 0 - sites processed: {column_count:,} ({column_count*100/total_sites:.1f}%)')
            for pileupread in pileupcolumn.pileups:
                try:
                    read_value=1.0/pileupread.get_tag('NH')
                except:
                    read_value=1.0
                # if not pileupread.is_del and not pileupread.is_refskip:
                merge_tag=False
                merge_rev=-1
                double_read=False
                if '_FR' in pileupread.alignment.query_name:
                    merge_tag=bin(int(pileupread.alignment.query_name.split('_FR')[-1],16))[3:]
                    if merge_tag[-1]=='0':
                        merge_rev=pileupread.alignment.query_length-1
                        while merge_tag[merge_rev-1]=='0': merge_rev-=1
                    double_read=bool(int(merge_tag[pileupread.query_position_or_next]))
                if pileupread.query_position!=None:
                    base=pileupread.alignment.query_sequence[pileupread.query_position]
                    query_pos=pileupread.query_position
                else:
                    base='*'
                    query_pos=pileupread.query_position_or_next
                reverse=pileupread.alignment.is_reverse
                if not double_read and merge_tag and query_pos>=merge_rev: reverse=True
                within_range=True
                if abs(pileupread.alignment.query_alignment_start-query_pos)<metrics['site_QC_criteria']['min_dist_from_read_end'] or abs(query_pos-pileupread.alignment.query_alignment_end)<metrics['site_QC_criteria']['min_dist_from_read_end']: within_range=False

                if within_range:
                    if double_read or not reverse: counts[base][0]+=1
                    if double_read or reverse: counts[base][1]+=1

                    templates[pileupread.alignment.query_name]=(base,read_value)

            arr.append((pileupcolumn.pos,sum([sum(l) for l in counts.values()]),sum([v[1] for v in templates.values()]),counts,templates.values()))

    results=[]
    excluded_results=[]

    for pos,total,true_total,counts,templates_list in arr:
        excluded=True
        if total>=metrics['site_QC_criteria']['min_cov']:
            res=[]
            for b in ['A','C','G','T','N','*']:
                if min(counts[b])>=metrics['site_QC_criteria']['min_read_count_per_direction']:
                    if sum(counts[b])/total>=metrics['site_QC_criteria']['min_alt_freq']:
                        res.append((sum([v[1] for v in templates_list if v[0]==b])*1.0/true_total,b))
                        excluded=False
            if len(res)>0:
                results.append((pos,res))
        if excluded:
            excl_res=(pos,[(sum([v[1] for v in templates_list if v[0]==b])*1.0/true_total,b) for b in ['A','C','G','T','N','*'] if sum([v[1] for v in templates_list if v[0]==b])>0])
            if len(excl_res[1])>0:
                excluded_results.append(excl_res)

    q.put((ID,results,excluded_results))
    return(True)

def snp_worker(ID,bamfileIDs,ref_ID,start,end,metrics,q):#doesn't deal w/ insertions right now
    arr=[{} for i in range(end-start)]
    true_count=[{} for i in range(end-start)]
    for bamfile in bamfileIDs:
        bamfile_pointer=pysam.AlignmentFile(bamfile,'rb')
        read_count=0
        expected_read_count=bamfile_pointer.count(ref_ID,start,end)
        for read in bamfile_pointer.fetch(ref_ID,start,end):
            if not read.is_duplicate and not read.is_unmapped:# and not read.is_secondary:
                read_count+=1
                if ID==0 and read_count%10000==0: print('Process 0 - current read count: %s (%.1f%%)'%('{:,}'.format(read_count),read_count*100.0/expected_read_count))
                pos_pairs=read.get_aligned_pairs()
                merge_tag=False
                merge_rev=-1
                if '_FR' in read.query_name:
                    merge_tag=bin(int(read.query_name.split('_FR')[-1],16))[3:]
                    if merge_tag[-1]=='0':
                        merge_rev=read.query_length-1
                        while merge_tag[merge_rev-1]=='0': merge_rev-=1
                while len(pos_pairs)>0:
                    qpos,rpos=pos_pairs.pop(0)
                    if (rpos!=None and start<=rpos<end and metrics['site_QC_criteria']['min_dist_from_read_end']<=rpos-read.reference_start and metrics['site_QC_criteria']['min_dist_from_read_end']<=read.reference_end-rpos) or (qpos!=None and start<=read.reference_start+qpos<end and metrics['site_QC_criteria']['min_dist_from_read_end']<=read.query_alignment_start-qpos and qpos<read.query_alignment_end-metrics['site_QC_criteria']['min_dist_from_read_end']):
                        try:
                            read_value=1.0/read.get_tag('NH')
                        except:
                            read_value=1.0
                        if qpos==None:
                            if merge_tag and ((len(merge_tag)>(rpos-read.reference_start) and merge_tag[rpos-read.reference_start]=='1') or (len(merge_tag)<=(rpos-read.reference_start) and merge_tag[-1]=="1")):
                                for both_dir in [True,False]:
                                    if not (both_dir,'*') in arr[rpos-start]: arr[rpos-start][(both_dir,'*')]=read_value
                                    else: arr[rpos-start][(both_dir,'*')]+=read_value
                                if not '*' in true_count[rpos-start]: true_count[rpos-start]['*']=read_value
                                else: true_count[rpos-start]['*']+=read_value

                            elif merge_tag and (rpos-read.reference_start)>=merge_rev:
                                if not (True,'*') in arr[rpos-start]: arr[rpos-start][(True,'*')]=read_value
                                else: arr[rpos-start][(True,'*')]+=read_value
                                if not '*' in true_count[rpos-start]: true_count[rpos-start]['*']=read_value
                                else: true_count[rpos-start]['*']+=read_value

                            else:
                                if not (read.is_reverse,'*') in arr[rpos-start]: arr[rpos-start][(read.is_reverse,'*')]=read_value
                                else: arr[rpos-start][(read.is_reverse,'*')]+=read_value
                                if not '*' in true_count[rpos-start]: true_count[rpos-start]['*']=read_value
                                else: true_count[rpos-start]['*']+=read_value
                        else:
                            if rpos!=None and start<=rpos<end+1:
                                q_base=read.query_sequence[qpos]
                                q_qs=[read.query_qualities[qpos]]
                                while len(pos_pairs)>0 and pos_pairs[0][1]==None:
                                    qpos2,rpos2=pos_pairs.pop(0)
                                    q_base+=read.query_sequence[qpos2]
                                    q_qs.append(read.query_qualities[qpos2])
                                if len(pos_pairs)==0:#to avoid having skipped read ends attached to the last base
                                    q_base=q_base[0]
                                    q_qs=q_qs[:1]
                                if min(q_qs)>=metrics['site_QC_criteria']['min_site_qs']:
                                    if merge_tag and '1' in merge_tag[qpos:qpos+len(q_base)] and math.ceil(merge_tag[qpos:qpos+len(q_base)].count('1')/2)>=len(q_base):
                                        for both_dir in [True,False]:
                                            if not (both_dir,q_base) in arr[rpos-start]: arr[rpos-start][(both_dir,q_base)]=read_value
                                            else: arr[rpos-start][(both_dir,q_base)]+=read_value
                                    elif merge_tag and qpos>=merge_rev:
                                        if not (True,q_base) in arr[rpos-start]: arr[rpos-start][(True,q_base)]=read_value
                                        else: arr[rpos-start][(True,q_base)]+=read_value
                                    else:
                                        if not (read.is_reverse,q_base) in arr[rpos-start]: arr[rpos-start][(read.is_reverse,q_base)]=read_value
                                        else: arr[rpos-start][(read.is_reverse,q_base)]+=read_value
                                    if not q_base in true_count[rpos-start]: true_count[rpos-start][q_base]=read_value
                                    else: true_count[rpos-start][q_base]+=read_value

    results=[]
    excluded_results=[]
    for pos in range(len(arr)):
        total=sum(arr[pos].values())
        true_total=sum(true_count[pos].values())
        excluded=True
        if total>=metrics['site_QC_criteria']['min_cov']:
            bases=set([base[1] for base in arr[pos].keys()])
            #filter for min read directions and freq
            res=[]
            for b in bases:
                if not (False,b) in arr[pos]: arr[pos][(False,b)]=0
                if not (True,b) in arr[pos]: arr[pos][(True,b)]=0
                if arr[pos][(False,b)]>=metrics['site_QC_criteria']['min_read_count_per_direction']:
                    if arr[pos][(True,b)]>=metrics['site_QC_criteria']['min_read_count_per_direction']:
                        if (arr[pos][(False,b)]+arr[pos][(True,b)])*1.0/total>=metrics['site_QC_criteria']['min_alt_freq']:
                            res.append((true_count[pos][b]*1.0/true_total,b))
                            excluded=False
            if len(res)>0:
                results.append((pos+start,res))
        else:
            if excluded:
                bases=set([base[1] for base in arr[pos].keys()])
                excl_res=(pos+start,[(true_count[pos][b]*1.0/true_total,b) for b in bases if true_count[pos][b]>0])
                if len(excl_res[1])>0:
                    excluded_results.append(excl_res)

    q.put((ID,results,excluded_results))
    return(True)

def getMedian(values):#calculation of median
  values = sorted(values)
  count = len(values)
  if count % 2 == 1:
    return values[int(round((count+1)/2-1))]
  else:
    return (float(values[int(round(count/2-1))] + values[int(round(count/2))])) / 2

def read_annotation(datei):#corrects coordinates to 0-based
    annotation=[]
    if datei.endswith('.gz'):
        infile=gzip.open(datei,'rt')
    else:
        infile=open(datei,'r')
    for line in infile:
        if not line.startswith('#'):
            line=line.strip('\r\n').strip('\n').split('\t')
            if len(line)>=9:
                annotation.append((int(line[3])-1,int(line[4])-1,line[6],[v for v in line[8].split(';') if v.split('=')[0].lower() in ['name','description','_description','note','symbol','_symbol']]))
    infile.close()
    return(annotation)

def compare_datasets(SNP_results,mask,annotation,critical_sites,metrics):#list w/ pos,DNA bases and freqs, RNA bases and freqs, critical, annotation
    results=[]
    #converting DNA to dict for lookup
    DNA=dict(SNP_results['DNA'][0])
    DNA_excl=dict(SNP_results['DNA'][1])
    count=0
    for pos,snp_list in SNP_results['RNA'][0]:
        count+=1
        if count%100000==0: print('Processed {:,} sites'.format(count))
        if len(snp_list)<=metrics['site_QC_criteria']['max_base_count_per_site']:
            if mask[pos]==0:
                if pos in DNA:
                    dna_bases=set([x[1] for x in DNA[pos]])
                    snp_list=[(x[1],x[0]) for x in snp_list if not x[1] in dna_bases]
                    dna_list=[(y[1],y[0]) for y in DNA[pos]]
                    critical='none'
                elif pos in DNA_excl:
                    dna_bases=set([x[1] for x in DNA_excl[pos]])
                    snp_list=[(x[1],x[0]) for x in snp_list if not x[1] in dna_bases]
                    dna_list=[(y[1],y[0]) for y in DNA_excl[pos]]
                    critical='insufficient DNA data'
                    # print(F'Insufficient DNA data for {pos}')
                else:
                    snp_list=[(x[1],x[0]) for x in snp_list]
                    dna_list='missing DNA info'
                    critical='none'
                if len(snp_list)>0:
                    results.append([pos,dna_list,snp_list])

                    #critical sites
                    # critical='none'
                    for c in critical_sites:
                        if len(c)==4:
                            if c[0]<=pos<=c[1]:
                                critical=c
                        elif len(c)==2:
                            if c[0]==pos:
                                critical=c
                    results[-1].append(critical)

                    #annotation
                    gene=False
                    for start,end,strand,description in annotation:
                        if start<=pos<=end and not gene:
                            results[-1].append((start,end,strand,description))
                            gene=True
                        elif start<=pos<=end and gene:
                            results.append(results[-1])
                            results[-1][-1]=(start,end,strand,description)
                    if not gene:
                        results[-1].append('none')
    return(results)

def write_results(phase,replicateID,results,metrics):
    mode='wt'
    if os.path.exists('%sRNA_editing_results.txt.gz'%metrics['prefix']): mode='at'
    outfile=gzip.open('%sRNA_editing_results.txt.gz'%metrics['prefix'],mode)
    if mode=='wt': outfile.write('\t'.join(['phase','replicateID','pos','DNA_info','RNA_info','critical_site','annotation'])+'\n')
    for liste in results[phase][replicateID]['filtered']:
        outfile.write('\t'.join([phase,replicateID]+[str(l) for l in liste])+'\n')
    outfile.close()

    mode='wt'
    if os.path.exists('%sRNA_editing_supplement.txt.gz'%metrics['prefix']): mode='at'
    outfile=gzip.open('%sRNA_editing_supplement.txt.gz'%metrics['prefix'],mode)
    if mode=='wt': outfile.write('\t'.join(['phase','replicateID','dataset_type','metric','result'])+'\n')
    for m_name,m_ID in [('average_coverage','avg_cov'),('average_gene_coverage','gene_cov'),('intergenic_coverage','intergenic_cov')]:
        for typ in results[phase][replicateID][m_ID].keys():
            outfile.write('\t'.join([phase,replicateID,typ,m_name,str(results[phase][replicateID][m_ID][typ])])+'\n')
    for typ in results[phase][replicateID]['missing_genes'].keys():
        for liste in results[phase][replicateID]['missing_genes'][typ]:
            outfile.write('\t'.join([phase,replicateID,typ,'uncovered_genes',str(liste)])+'\n')
    outfile.close()

    return('%sRNA_editing_results.txt.gz'%metrics['prefix'],'%sRNA_editing_supplement.txt.gz'%metrics['prefix'])

def reverse_complement(string):
  table=str.maketrans('ACGTNWSMKRYBDHV.acgtnwsmkrybdhv','TGCANWSKMYRVHDBNtgcanwskmyrvhdb')
  string=string.translate(table)[::-1]
  return(string)

def focus_analysis(in_name,prefix):
    print(f'reading result file: {in_name}')
    dataset={}
    with gzip.open(in_name,'rt') as infile:
        head=False
        for line in infile:
            if head:
                line=line.strip('\n').split('\t')
                if not line[0] in dataset: dataset[line[0]]={}
                if not line[1] in dataset[line[0]]: dataset[line[0]][line[1]]={}
                for rel_list_pos in [3,4]:
                    if not line[rel_list_pos].startswith('missing '):
                        line[rel_list_pos]=[(x.split(', ')[0][1:-1],float(x.split(', ')[1])) for x in line[rel_list_pos][2:-2].split('), (')]
                if line[-1]!='none':
                    line[-1]=line[-1][1:-1].split(', ')
                    line[-1][0]=int(line[-1][0])
                    line[-1][1]=int(line[-1][1])
                    line[-1][2]=line[-1][2][1]
                    try:
                        line[-1][3]=dict([x[1:-1].split('=') for x in ','.join(line[-1][3:])[1:-1].split(',')])
                    except:
                        line[-1][3]={}
                    line[-1]=line[-1][:4]
                if not int(line[2]) in dataset[line[0]][line[1]]:

                    dataset[line[0]][line[1]][int(line[2])]=line[3:]
                    dataset[line[0]][line[1]][int(line[2])][-1]=[dataset[line[0]][line[1]][int(line[2])][-1]]
                else:
                    dataset[line[0]][line[1]][int(line[2])][-1].append(line[-1])
            else:
                head=True

    # filter
    print('filtering against missing DNA info')
    filtered={}
    for phase in dataset:
        filtered[phase]={}
        for rep in sorted(dataset[phase].keys()):
            filtered[phase][rep]={}
            count=0
            for pos in dataset[phase][rep].keys():
                if dataset[phase][rep][pos][0]=='missing DNA info': count+=1
                else:
                    filtered[phase][rep][pos]=dataset[phase][rep][pos]

    # transition counts:
    print('counting replicate appearances')
    transition_counts={}
    for phase in filtered:
        transition_counts[phase]={}
        for rep in filtered[phase]:
            for pos in filtered[phase][rep]:
                for b1 in filtered[phase][rep][pos][0]:
                    b1=b1[0]
                    for b2 in filtered[phase][rep][pos][1]:
                        b2=b2[0]
                        if filtered[phase][rep][pos][3][0][2]=='-':
                            b1=reverse_complement(b1)
                            b2=reverse_complement(b2)
                        if not (b1,b2) in transition_counts[phase]: transition_counts[phase][(b1,b2)]={}
                        if pos in transition_counts[phase][(b1,b2)]: transition_counts[phase][(b1,b2)][pos].append(rep)
                        else: transition_counts[phase][(b1,b2)][pos]=[rep]

    # for phase in transition_counts:
    #     for pair in transition_counts[phase]:
    #         count=0
    #         for pos in transition_counts[phase][pair]:
    #             if len(transition_counts[phase][pair][pos])>0: count+=1

    pos_results={}
    for phase in transition_counts:
        for direction in [('A','G')]:#,('C','T')]:
            if not '->'.join(direction) in pos_results: pos_results['->'.join(direction)]={}
            if direction in transition_counts[phase]:
                for pos in transition_counts[phase][direction]:
                    for rep in transition_counts[phase][direction][pos]:
                        if not pos in pos_results['->'.join(direction)]: pos_results['->'.join(direction)][pos]={phase:[rep]}
                        elif not phase in pos_results['->'.join(direction)][pos]: pos_results['->'.join(direction)][pos][phase]=[rep]
                        else: pos_results['->'.join(direction)][pos][phase].append(rep)

    for direction in pos_results.keys():
        count=[]
        for pos in list(pos_results[direction].keys()):
            ctr=True
            for phase in pos_results[direction][pos]:
                if len(pos_results[direction][pos][phase])>1: ctr=False
            if ctr:
                pos_results[direction].pop(pos)
                count.append(pos)
        print('excluded due to not sharing among reps ',direction,len(count))#,sorted(count))

    pos_reporting_xtr={}
    with open('%sRNA_edit_focus.results.txt'%prefix,'w') as outfile:
        header=['edit_direction','shared_phases','pos','gene']
        column_order=[]
        for phase in filtered:
            for rep in filtered[phase]:
                column_order.append((phase,rep))
                header.append(F"{phase}_{rep}")
        outfile.write('\t'.join(header)+'\n')

        # outfile.write('edit_direction\tshared_phases\tpos\tgene\t%s\n'%('\t'.join(['\t'.join(['%s_%s'%(phase,rep) for rep in ['01','02','03']]) for phase in ['na']])))
        for direction in pos_results.keys():
            if not direction in pos_reporting_xtr: pos_reporting_xtr[direction]=set()
            for length in [3,2,1]:
                for pos in sorted(pos_results[direction].keys()):
                    if len(list(pos_results[direction][pos].keys()))==length:
                        gene=''
                        complement='+'
                        nucs=[]
                        for phase,rep in column_order:
                            if phase in filtered:
                                if rep in filtered[phase]:
                                    if pos in filtered[phase][rep]:
                                        if filtered[phase][rep][pos][3][0][2]=='-':
                                            temp=['%s(%.3f)'%(reverse_complement(b),f) for b,f in filtered[phase][rep][pos][1]]
                                            complement='-'
                                        else:
                                            temp=['%s(%.3f)'%(b,f) for b,f in filtered[phase][rep][pos][1]]
                                        nucs.append(';'.join(temp))
                                        if gene=='':
                                            gene=str(filtered[phase][rep][pos][3])
                                    else:
                                        nucs.append('n/a')
                                else:
                                    nucs.append('n/a')
                            else:
                                nucs.append('n/a')
                        outfile.write('%s\tna\t%i\t%s\t%s\n'%(direction,pos,gene,'\t'.join(nucs)))
                        pos_reporting_xtr[direction].add((pos,complement))
    return(list(pos_reporting_xtr.items()))

def xtr_site(bamfile,pos,metrics,orient):
    counts=[[0,0],[0,0],[0,0],[0,0],[0,0]]
    base_dict={'A':0,'C':1,'G':2,'T':3,'*':4}
    if orient=='-': base_dict={'A':3,'C':2,'G':1,'T':0,'*':4}
    bamfile_pointer=pysam.AlignmentFile(bamfile,'rb')
    results={}
    for read in bamfile_pointer.fetch(metrics['ref_ID'],pos,pos+1):
        if not read.is_duplicate and not read.is_unmapped:
            pos_pairs=read.get_aligned_pairs()
            index=[p[1] for p in pos_pairs].index(pos)
            index2=index+1
            while read.query_alignment_end>index2 and pos_pairs[index2][1]==None: index2+=1
            key=None
            if (index-read.query_alignment_start)>=metrics['site_QC_criteria']['min_dist_from_read_end']:
                if (read.query_alignment_end-index2-1)>=metrics['site_QC_criteria']['min_dist_from_read_end']:
                    pos_pairs=pos_pairs[index:index2]
                    rseq=''.join([read.query_sequence[p[0]] if p[0]!=None else '*' for p in pos_pairs])
                    rSQ=[read.query_qualities[p[0]] if p[0]!=None else 40 for p in pos_pairs]
                    if min(rSQ)>=metrics['site_QC_criteria']['min_site_qs']:
                        key=(read.is_reverse,rseq)
            if key!=None:
                if key in results: results[key]+=1
                else: results[key]=1
    bamfile_pointer.close()

    for key,val in results.items():
        if key[1] in base_dict:
            if key[0]:#reverse
                counts[base_dict[key[1]]][1]+=val
            else:
                counts[base_dict[key[1]]][0]+=val

    return(counts)

def comprehensive_site_reporting(sites,metrics):

    data_order={'na':set([x[0] for x in metrics['input_files']['DNA']])}

    header=['direction','pos','orientation']
    for h in ['DNA','RNA']:
        for p in ['na']:
            for s in sorted(data_order[p]):
                for n in ['A','C','G','T','*']:
                    for r in ['fwd','rev']:
                        header.append('_'.join([h,p,s,n,r]))

    data_order=[[(phase,x) for x in sorted(data_order[phase])] for phase in ['na']]
    data_order=[x for y in data_order for x in y]
    file_order=['%s%s_DNA_sort.bam'%(metrics['mapping_file_base_path'],b) for a,b in data_order]+['%s%s_RNA_sort.bam'%(metrics['mapping_file_base_path'],b) for a,b in data_order]

    with open('%scomprehensive_results.txt'%(metrics['prefix']),'w') as outfile:
        outfile.write('\t'.join(header)+'\n')
        for direction,liste in sites:
            for pos,orient in sorted(liste):
                line=[direction,str(pos),orient]
                for datei in file_order:
                    # print('Extracting pos %i for %s in %s'%(pos,direction,datei.split('/')[-1]))
                    line+=[str(x) for y in xtr_site(datei,pos,metrics,orient) for x in y]
                outfile.write('\t'.join(line)+'\n')
    return(True)

def orchestrate_analysis(metrics):
    print('Creating baseline folder tree')
    if not os.path.exists('analysis_temp/'): os.makedirs('analysis_temp')
    path='analysis_temp/'
    if not metrics['skip_mapping']:
        print('Performing read mapping')
        mapping_alternative(metrics,path)

    #masking genome based on low complexity (homopol)## and paralogs
    print('Masking genome')
    mask=mask_genome(metrics)

    print('Reading annotation file')
    annotation=read_annotation(metrics['ref_anno_file'])

    #create sampleID based data structure for processing
    print('Switching to sample ID based data analysis structure')

    # mapping_file_base_path='analysis_temp/07_final_sort/'
    mapping_file_base_path='analysis_temp/05_template_reconstruction/'


    metrics['mapping_file_base_path']=mapping_file_base_path

    analysis_results={}

    # print('Preparing pooled DNA datasets for all phases across replicates')
    # DNA_files_for_pooling={}
    # for sampleID in sorted(set([x[0] for x in metrics['input_files']['DNA']])):
    #     phase=sampleID[3:5]
    #     replicateID=sampleID[-2:]
    #     if not phase in DNA_files_for_pooling:
    #         DNA_files_for_pooling[phase]=[]
    #     DNA_files_for_pooling[phase].append(sampleID)
    #
    #
    # combinded_DNA={}
    # for phase in DNA_files_for_pooling:
    #     combinded_DNA[phase]=perform_SNP_ident(['%s%s_DNA_sort.bam'%(mapping_file_base_path,f) for f in DNA_files_for_pooling[phase]],metrics['ref_ID'],metrics)
    #
    starttime=time.time()
    counter=0
    counter_total=len(set([x[0] for x in metrics['input_files']['DNA']]))

    outfile_name=''
    outfile_supplement=''


    if os.path.exists('%sRNA_editing_results.txt.gz'%metrics['prefix']):
        text='NOT '
        if metrics['skip'][8]=='1':
            text=''
            answer=input(F"Removing old files {metrics['prefix']}RNA_editing_results.txt.gz and {metrics['prefix']}RNA_editing_supplement.txt.gz? [y/n] ({text}recommended):")
            if answer[0].lower()=='y':

                os.remove('%sRNA_editing_results.txt.gz'%metrics['prefix'])
                os.remove('%sRNA_editing_supplement.txt.gz'%metrics['prefix'])

    if metrics['skip'][8]=='1':

        for sampleID in sorted(set([x[0] for x in metrics['input_files']['DNA']])):
            phase='na'
            replicateID=sampleID
            if not phase in analysis_results: analysis_results[phase]={}
            if not replicateID in analysis_results[phase]:
                analysis_results[phase][replicateID]={'avg_cov':{},'gene_cov':{},'intergenic_cov':{},'missing_genes':{},'critical_sites':[],'filtered':[]}

            #looking for mapping artefacts
            for typ in ['RNA','DNA']:
                print('Looking for structural abberations in %s data for %s'%(typ,sampleID))
                arr_cov,critical_ranges_and_sites=qualify_mapping_results('%s%s_%s_sort.bam'%(mapping_file_base_path,sampleID,typ),metrics['ref_ID'],typ,metrics['CPU_count'])
                if typ=='DNA': analysis_results[phase][replicateID]['critical_sites']=critical_ranges_and_sites
                overall_cov,gene_cov,intergenic,missing=analyze_coverages(arr_cov,annotation)
                analysis_results[phase][replicateID]['avg_cov'][typ]=overall_cov
                analysis_results[phase][replicateID]['gene_cov'][typ]=gene_cov
                analysis_results[phase][replicateID]['intergenic_cov'][typ]=intergenic
                analysis_results[phase][replicateID]['missing_genes'][typ]=missing

            #identifying annotated sites w/o coverage = blind sites

            #identifying SNPs
            SNP_results={}
            for typ in ['RNA','DNA']:
                print(F'Identifying SNPs in {sampleID} {typ}')
                SNP_results[typ]=perform_SNP_ident(['%s%s_%s_sort.bam'%(mapping_file_base_path,sampleID,typ)],metrics['ref_ID'],metrics)

            # SNP_results['DNA']=combinded_DNA[phase]

            print(F'Comparing RNA to DNA results and applying additional filter for {sampleID}')
            analysis_results[phase][replicateID]['filtered']=compare_datasets(SNP_results,mask,annotation,analysis_results[phase][replicateID]['critical_sites'],metrics)

            print(F'writing results to file for {sampleID}')
            outfile_name,outfile_supplement=write_results(phase,replicateID,analysis_results,metrics)

            counter+=1
            print(f"est. remaining runtime for ind. sample analysis {((time.time()-starttime)/counter)*(counter_total-counter)/60:,.2f} min")


    if metrics['skip'][9]=='1':
        outfile_name='%sRNA_editing_results.txt.gz'%(metrics['prefix'])
        sites=focus_analysis(outfile_name,metrics['prefix'])

        comprehensive_site_reporting(sites,metrics)

        create_coverage_plots(sites,metrics,mapping_file_base_path)

    #cleaning up
    print('Cleaning up',flush=True)
    for datei in os.listdir('./'):
        if datei.startswith('pipeline.log'):
            os.remove(datei)

    #packaging
    print('packing up result files',flush=True)
    now=datetime.datetime.now()
    outfolder=F"{datetime.datetime.strftime(now,'%y%m%d%H%M')}_results_{metrics['prefix'][:-1]}"
    os.makedirs(outfolder)
    outfolder+='/'
    for suffix in ['RNA_edit_focus.results.txt','RNA_editing_results.txt.gz','RNA_editing_supplement.txt.gz','comprehensive_results.txt']:
        shutil.move(F"{metrics['prefix']}{suffix}",F"{outfolder}")
    shutil.move('analysis_temp/cov_plots',F"{outfolder}cov_plots")
    if os.path.exists(f"analysis_temp/{metrics['prefix']}attrition.svg"):
        shutil.move(f"analysis_temp/{metrics['prefix']}attrition.svg",f"{outfolder}")

    #report run parameters
    print('Reporting run parameters',flush=True)
    with open(F"{outfolder}{metrics['prefix']}analysis_run_parameters.txt",'w') as outfile:
        outfile.write("Program options:\n")
        for attr in ['infile','ref','annotation','cores','outname','cov','direction','no_mapping','phix','skip','no_umi']:
            outfile.write(F"{attr}\t{getattr(options,attr)}\n")
        outfile.write('\nQC filter criteria:\n')
        for attr in ['min_site_qs','max_base_count_per_site','min_alt_freq','max_homopol_len','min_dist_to_paralogs','min_dist_from_read_end','DNA_max_cov_diff_SD','DNA_cov_win_size']:
            outfile.write(F"{attr}\t{metrics['site_QC_criteria'][attr]}\n")


    return(True)

def worker_cov_plot(target_bin,source,pos,ID,q):
    pos_range=30

    class read_obj:
        def __init__(self,start,end,reverse,secondary,trim_list):
            self.start=start
            self.end=end
            self.reverse=reverse
            self.secondary=secondary
            self.trim=trim_list

    def convert_read_to_obj(read,pos,pos_range):
        start=read.reference_start
        if start<pos-pos_range:
            pre_start=start=pos-pos_range
        else:
            if read.query_alignment_start>0:
                pre_start=start-read.query_alignment_start
                if pre_start<pos-pos_range: pre_start=pos-pos_range
            else:
                pre_start=start

        end=read.reference_end
        if end>pos+pos_range:
            post_end=end=pos+pos_range
        else:
            if read.query_alignment_end<read.query_length:
                post_end=end+(read.query_length-read.query_alignment_end)
                if post_end>pos+pos_range: post_end=pos+pos_range
            else:
                post_end=end

        reverse=read.is_reverse
        secondary=read.is_secondary
        trim_list=[]
        if pre_start<start: trim_list.append((pre_start,start))
        if post_end>end: trim_list.append((end,post_end))
        return(read_obj(start,end,reverse,secondary,trim_list))

    def create_plot(read_list,position,pos_range,ref_ID,infile,target_bin,outfile=''):

        width,height=[6,4]
        fig,ax1=plt.subplots(figsize=(width,height), dpi=300)

        colors={(False,False):(.3,.3,.8),(False,True):(.6,.6,.8),(True,False):(.8,.2,.2),(True,True):(.8,.5,.5),'trim':'gray'}#(reverse,secondary)
        color_codes=[((False,False),'fwd primary'),((True,False),'rev primary'),('trim','trim')]

        read_dist=int(round(pos_range/15))
        temp=[]
        for e,obj in enumerate(read_list):
            if len(obj.trim)==0:
                temp.append((obj.start,obj.end,e,obj))
            else:
                temp.append((min([obj.start]+[t[0] for t in obj.trim]),max([obj.end]+[t[1] for t in obj.trim]),e,obj))
        read_list=sorted(temp)

        height_table=[-read_dist-1]
        for s,e,_,obj in read_list:
            ctr=True
            for pos in range(len(height_table)):
                if height_table[pos]<s:
                    ax1.broken_barh([(obj.start,obj.end-obj.start)],(pos,1),facecolors=colors[(obj.reverse,obj.secondary)],edgecolor=(1,1,1),lw=.1)
                    height_table[pos]=e
                    ctr=False
                    h=pos
                    break
            if ctr:
                height_table.append(e)
                ax1.broken_barh([(s,e-s)],(len(height_table)-1,1),facecolors=colors[(obj.reverse,obj.secondary)],edgecolor=(1,1,1),lw=.1)
                h=len(height_table)-1
            if obj.trim!=[]:
                for ts,te in obj.trim:
                    ax1.broken_barh([(ts,te-ts)],(h,1),facecolor=colors['trim'],edgecolor=(1,1,1))

        ax1.axvline(x=position,color=(1,1,1),lw=1.5)
        ax1.axvline(x=position,color=(.8,.5,0),lw=1.0)

        ax1.set_xlim(position-pos_range,position+pos_range)
        ax1.set_ylim(0,len(height_table)*1.1)
        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize('x-small')

        legend_elements=[]
        for key,name in color_codes:
            legend_elements.append(Patch(facecolor=colors[key],edgecolor=(1,1,1),label=name))
        plt.legend(handles=legend_elements,fontsize='x-small',loc='upper center', ncol=len(color_codes),framealpha=1, scatterpoints=1)

        ax1.set_xticks([position-pos_range,position,position+pos_range])
        ax1.set_xticklabels([F"-{pos_range}bp",F"{position:,}",F"+{pos_range}bp"],ha='center',fontsize='x-small')

        plt.title(F"{infile.split('/')[-1]}: ref ID: {ref_ID}",fontsize='small')

        ax1.set_xlabel('Reference position',fontsize='small')
        ax1.set_ylabel('Read coverage',fontsize='small')
        plt.tight_layout()

        if outfile=='': outfile=target_bin+infile.split('/')[-1].rsplit('.',1)[0]+F".{position}.png"
        plt.savefig(outfile,format='png')
        print(F"Coverage plot file {outfile} written.",flush=True)
        plt.close()

        return(True)

    try:
        read_list=[]
        bamfile_pointer=pysam.Samfile(source,'rb',check_sq=False)
        ref_len,ref_ID=sorted([(bamfile_pointer.get_reference_length(r),r) for r in bamfile_pointer.references],reverse=True)[0]
        for read in bamfile_pointer.fetch(ref_ID,pos-pos_range,pos+pos_range):
            if not read.is_unmapped:
                if not read.is_secondary:
                    read_list.append(convert_read_to_obj(read,pos,pos_range))
        if len(read_list)>200:
            read_list=random.sample(read_list,200)

        if len(read_list)>0:
            create_plot(read_list,pos,pos_range,ref_ID,source,target_bin)
    except:
        print(F"Coverage plot file {outfile} failed!",flush=True)

    q.put(ID)
    return(True)

def create_coverage_plots(sites,metrics,mapping_file_base_path):
    pos_list=set()
    for direction,liste in sites:
        pos_list|=set([pos for pos,orient in liste])
    pos_list=sorted(list(pos_list))

    to_do_list=[]
    for mapping_location,file_filter,bin in [('analysis_temp/04_deduplication/','clean.bam','non_merged'),(mapping_file_base_path,'.bam','merged')]:
        for datei in os.listdir(mapping_location):
            if datei.endswith('.bam') and file_filter in datei:
                for pos in pos_list:
                    to_do_list.append(('analysis_temp/cov_plots/'+bin,mapping_location+datei,pos))

    for bin in ['merged','non_merged']:
        if not os.path.exists('analysis_temp/cov_plots/'+bin): os.makedirs('analysis_temp/cov_plots/'+bin)

    thread_list=[]
    q=Queue()
    ID_count=0
    # report_length=False
    start_time=None
    while len(thread_list)>0 or len(to_do_list)>0:
        while len(thread_list)<metrics['CPU_count'] and len(to_do_list)>0:
            target_bin,source,pos=to_do_list.pop()
            thread_list.append((ID_count,Process(target=worker_cov_plot,args=(target_bin,source,pos,ID_count,q))))
            thread_list[-1][1].deamon=True
            thread_list[-1][1].start()
            ID_count+=1
        # if not report_length and len(to_do_list)==0:
        #     print("\nThe plotting todo list is empty.\n")
        #     report_length=True
        if not q.empty():
            start_time=time.time()
            tID=q.get()
            index=[i for i,p in thread_list].index(tID)
            # thread_list[index][1].join()
            thread_list.pop(index)
        if start_time!=None and time.time()-start_time>60:
            print('Waiting time is up. Starting to kill subpocesses.')
            for i,p in thread_list:
                p.terminate()
            thread_list=[]
    return(True)

def reading_bamfile(infile,SE,merging_q,writing_q,process_ctr,cores):
    buffer={}
    bamfile=pysam.AlignmentFile(infile,'rb')
    count=0
    for read in bamfile.fetch(until_eof=True):
        count+=1
        if count%100000==0: print(F'Parsed {count:,} reads.')
        if (not read.is_unmapped) and (not read.mate_is_unmapped) and (not read.is_secondary):
            if not read.query_name in buffer:
                buffer[read.query_name]=read_obj(read)
            else:
                if read.reference_start<=buffer[read.query_name].reference_start<read.reference_end or read.reference_start<buffer[read.query_name].reference_end<=read.reference_end or (read.reference_start==buffer[read.query_name].reference_start and buffer[read.query_name].reference_end==read.reference_end):#overlap
                    #making sure, that the first read is in fwd orientation
                    if buffer[read.query_name].is_reverse:
                        merging_q.put((read_obj(read),buffer.pop(read.query_name)))
                    else:
                        merging_q.put((buffer.pop(read.query_name),read_obj(read)))
                else:
                    writing_q.put(read_obj(read))
                    writing_q.put(buffer.pop(read.query_name))
        elif (not read.is_unmapped) and (read.mate_is_unmapped or SE) and (not read.is_secondary):
            writing_q.put(read_obj(read))
    for read in buffer.values():#to empty buffer
        writing_q.put(read)
    process_ctr.put('reading')
    for i in range(cores):
        merging_q.put('done')#one for each merge process to terminate them
    writing_q.put('done')
    bamfile.close()
    return(True)

def write_file(filename_prefix,writing_q,countdown,process_ctr):
    count=0
    outfile_single=gzip.open(filename_prefix+'.single.fq.gz','wt')
    outfile_fwd=gzip.open(filename_prefix+'.fwd.fq.gz','wt')
    outfile_rev=gzip.open(filename_prefix+'.rev.fq.gz','wt')
    buffer={}
    while countdown>0:
        if not writing_q.empty():
            read_obj=writing_q.get()
            if read_obj=='done':
                countdown-=1
            else:
                if read_obj.merged:
                    outfile_single.write('@%s\n%s\n+\n%s\n'%(read_obj.query_name,read_obj.query_sequence,''.join([chr(qs+33) for qs in read_obj.query_qualities])))
                else:
                    if read_obj.query_name in buffer:
                        if read_obj.read1:
                            outfile_fwd.write('@%s\n%s\n+\n%s\n'%(read_obj.query_name,read_obj.query_sequence,''.join([chr(qs+33) for qs in read_obj.query_qualities])))
                            outfile_rev.write('@%s\n%s\n+\n%s\n'%(read_obj.query_name,reverse_complement(buffer[read_obj.query_name].query_sequence),''.join([chr(qs+33) for qs in buffer[read_obj.query_name].query_qualities[::-1]])))
                            del buffer[read_obj.query_name]

                        else:
                            outfile_fwd.write('@%s\n%s\n+\n%s\n'%(read_obj.query_name,reverse_complement(buffer[read_obj.query_name].query_sequence),''.join([chr(qs+33) for qs in buffer[read_obj.query_name].query_qualities[::-1]])))
                            outfile_rev.write('@%s\n%s\n+\n%s\n'%(read_obj.query_name,read_obj.query_sequence,''.join([chr(qs+33) for qs in read_obj.query_qualities])))
                            del buffer[read_obj.query_name]
                    else:
                        buffer[read_obj.query_name]=read_obj
                count+=1
                if count%100000==0: print(F'{count:,} reads written')
    outfile_single.close()
    outfile_fwd.close()
    outfile_rev.close()
    print(F'Finished writing of {count:,} reads for file prefix "{filename_prefix}".')
    process_ctr.put('writing')
    return(True)

def merge(pair,qs_str=False):
    l1=[([None,0] if tup[0]==None else [pair[0].query_sequence[tup[0]],pair[0].query_qualities[tup[0]]] ,tup[1]) for tup in pair[0].pairs]
    l2=[([None,0] if tup[0]==None else [pair[1].query_sequence[tup[0]],pair[1].query_qualities[tup[0]]] ,tup[1]) for tup in pair[1].pairs]

    p1=0
    while l1[p1][1]==None: p1+=1
    l1=l1[p1:]


    min_ref_pos=l1[0][1]

    p2=0
    while l2[p2][1]==None or l2[p2][1]<min_ref_pos: p2+=1
    l2=l2[p2:]
    p2=len(l2)-1
    while l2[p2][1]==None: p2-=1
    l2=l2[:p2+1]

    max_ref_pos=l2[-1][1]

    p1=len(l1)-1
    while l1[p1][1]==None or l1[p1][1]>max_ref_pos: p1-=1
    l1=l1[:p1+1]

    liste=[]
    t1=l1.pop(0)
    t2=l2.pop(0)
    ctr1=0
    ctr2=0
    while len(l1)>0 and len(l2)>0:
        if t1[1]==t2[1]:
            ctr1=1
            ctr2=1
            liste.append([t1[0],t2[0],ctr1*ctr2])
            t1=l1.pop(0)
            t2=l2.pop(0)
        elif t1[1]==None or (t2[1]!=None and t1[1]<t2[1]):
            ctr1=1
            liste.append([t1[0],[None,0],ctr1*ctr2])
            t1=l1.pop(0)
        elif t2[1]==None or (t1[1]!=None and t2[1]<t1[1]):
            ctr2=1
            liste.append([[None,0],t2[0],ctr1*ctr2])
            t2=l2.pop(0)

    if t1[1]==t2[1]:
        ctr1=1
        ctr2=1
        liste.append([t1[0],t2[0],ctr1*ctr2])
    elif t1[1]==None or (t2[1]!=None and t1[1]<t2[1]):
        ctr1=1
        liste.append([t1[0],[None,0],ctr1*ctr2])
        ctr1=0
        liste.append([[None,0],t2[0],ctr1*ctr2])
    elif t2[1]==None or (t1[1]!=None and t2[1]<t1[1]):
        ctr2=1
        liste.append([[None,0],t2[0],ctr1*ctr2])
        ctr2=0
        liste.append([t1[0],[None,0],ctr1*ctr2])

    liste+=[[[None,0],t[0],0] for t in l2]#should there be anything left in l1, it goes beyond the 5prim start of the reverse read and must be an artifact

    #filling in gap-quality scores
    for rank in [0,1]:
        start=0
        while liste[start][rank][0]==None: start+=1
        end=len(liste)-1
        while liste[end][rank][0]==None: end-=1
        pos=start
        while pos<end:#to modify blocks if necessary
            if liste[pos][rank][0]==None:
                pos2=pos
                while liste[pos2][rank][0]==None: pos2+=1
                avg_qs=int(round(average([liste[pos-1][rank][1],liste[pos2][rank][1]])))
                for pos3 in range(pos,pos2+1):
                    liste[pos3][rank][1]=avg_qs
                pos=pos2+1
            else:
                pos+=1

    #creating consensus seq and new quality scores
    consensus=[[],[],[],[]]#seq,QS,overlap_info,ref_pos
    for t1,t2,m in liste:
        if m==1:#merged
            # if t1[0]!=None and t2[0]!=None:#otherwise, this template has a deletion
            consensus[2].append('1')
            if t1[0]==t2[0]:
                if t1[0]!=None:#ignoring alignment gaps supported by both reads
                    consensus[0].append(t1[0])
                    new_qs=t1[1]+t2[1]
                    if new_qs>60: new_qs=60
                    consensus[1].append(new_qs)
            else:
                if t1[1]==t2[1]:
                    consensus[0].append('N')
                    consensus[1].append(1)
                else:
                    temp=sorted([(t1[1],t1[0]),(t2[1],t2[0])])
                    if temp[1][1]!=None:#ignoring alignment gaps supported byhigher QS
                        consensus[0].append(temp[1][1])
                        consensus[1].append(temp[1][0]-temp[0][0])

        else:
            consensus[2].append('0')
            if not (t1[0]==None and t2[0]==None):
                if t1[0]==None:#ignoring alignment gaps in non-merged read fragments
                    consensus[0].append(t2[0])
                    consensus[1].append(t2[1])
                elif t2[0]==None:
                    consensus[0].append(t1[0])
                    consensus[1].append(t1[1])

    merged=merge_obj()
    merged.query_name=pair[0].query_name+'_FR'+hex(int('1'+''.join(consensus[2]),2))#to encode overlapping region
    merged.reference_name=pair[0].reference_name
    merged.query_sequence=''.join(consensus[0])
    if qs_str:
        merged.query_qualities=''.join([chr(qs+33) for qs in consensus[1]])
    else:
        merged.query_qualities=consensus[1]
    # merged.reference_start=ref_start
    # merged.reference_end=ref_end
    merged.merged=True

    if len(merged.query_sequence)!=len(merged.query_qualities):
        print(f'error QS length of {len(merged.query_qualities):,}')
        sys.exit()
    # if max(consensus[1])>60:
    #     sys.exit('QS error')

    return(merged)

def merge_read_worker(merging_q,writing_q,ID):
    waiting=True
    count=0
    while waiting:
        if not merging_q.empty():
            pair=merging_q.get()
            count+=2
            if pair=='done':
                waiting=False
                writing_q.put('done')
                # print('Merging worker #%i finished. (%i reads)'%(ID,count*2))
            else:
                writing_q.put(merge(pair))
    return(True)

class merge_obj:
    def __init__(self):
        self.query_name=''
        self.query_sequence=''
        self.reference_name=''
        self.reference_start=0
        self.reference_end=0
        self.query_qualities=[]
        self.pairs=[]

class read_obj:
    def __init__(self,read,merged=False):
        self.query_name=read.query_name
        self.read1=read.is_read1
        self.query_sequence=read.query_sequence
        self.reference_name=read.reference_name
        self.reference_start=read.reference_start
        self.reference_end=read.reference_end
        self.query_qualities=read.query_qualities
        self.merged=merged
        self.is_reverse=read.is_reverse
        self.pairs=read.get_aligned_pairs()
        self.cigartuples=read.cigartuples
        self.mapping_quality=read.mapping_quality

def merging_manager2(bamfile,cores,filename_prefix):

    def read_bamfile2(infile,cores):
        liste=[]
        bamfile=pysam.AlignmentFile(infile,'rb',threads=cores)
        count=0
        for read in bamfile.fetch(until_eof=True):
            count+=1
            if count%100000==0: print(F'Parsed {count:,} reads from {infile.split("/")[-1]}.')
            if (not read.is_unmapped) and (not read.is_secondary):
                liste.append((read.query_name,read_obj(read)))
        bamfile.close()
        dataset={n:[] for n,o in liste}
        for n,o in liste:
            dataset[n].append(o)
        return(dataset)

    def merge_read_worker2(ID,read_pairs,q):
        results=[]
        step=int(len(read_pairs)/10)
        # count=0
        for pair in read_pairs:
            results.append(merge(pair,qs_str=True))
            # count+=1
            # if ID==0 and count%step==0: print(F"{count/step*10}% of eligible read merged.",flush=True)
        q.put((ID,results))
        return(True)

    def merge_read_worker3(ID,read_pairs,q):
        try:
            step=int(len(read_pairs)/10)
            count=0
            outname=F"mtemp_{ID}.fq.gz"
            with gzip.open(outname,'wt') as outfile:
                for pair in read_pairs:
                    mread=merge(pair,qs_str=True)
                    outfile.write('@%s\n%s\n+\n%s\n'%(mread.query_name,mread.query_sequence,mread.query_qualities))
                    count+=1
                    if ID==0 and count%step==0: print(F"{count/step*10}% of eligible read merged.",flush=True)
            q.put((ID,outname))
        except:
            print(F"\n\nERROR in merge read worker3 #{ID}!\n\n",flush=True)
        return(True)

    if not os.path.exists(bamfile):
        sys.exit('Input file "%s" does not exist. Aborting.'%bamfile)

    start_time_check=time.time()

    #reading all in list as read obj
    # init dict with read names
    #add read objects
    dataset=read_bamfile2(bamfile,cores)

    IDs=dataset.keys()

    gzip_threads=8
    if cores<gzip_threads: gzip_threads=cores

    outfile_single=gzip.open(filename_prefix+'.single.fq.gz','wt',compresslevel=4)
    outfile_fwd=gzip.open(filename_prefix+'.fwd.fq.gz','wt',compresslevel=4)
    outfile_rev=gzip.open(filename_prefix+'.rev.fq.gz','wt',compresslevel=4)

    merge_IDs=[]#(ID,correct order)
    #write single end reads to file directly
    written_count=0
    for ID in IDs:
        if len(dataset[ID])==1:
            outfile_single.write('@%s\n%s\n+\n%s\n'%(dataset[ID][0].query_name,dataset[ID][0].query_sequence,''.join([chr(qs+33) for qs in dataset[ID][0].query_qualities])))
            written_count+=1
        elif dataset[ID][0].reference_start<=dataset[ID][1].reference_start<dataset[ID][0].reference_end or dataset[ID][0].reference_start<dataset[ID][1].reference_end<=dataset[ID][0].reference_end:
            #making sure, that the first read is in fwd orientation
            if dataset[ID][0].is_reverse:
                merge_IDs.append((ID,False))
            else:
                merge_IDs.append((ID,True))
        else: #handle paired reads - if not overlap, write to file
            written_count+=2
            if dataset[ID][0].is_reverse:
                outfile_fwd.write('@%s\n%s\n+\n%s\n'%(dataset[ID][1].query_name,dataset[ID][1].query_sequence,''.join([chr(qs+33) for qs in dataset[ID][1].query_qualities])))
                outfile_rev.write('@%s\n%s\n+\n%s\n'%(dataset[ID][0].query_name,reverse_complement(dataset[ID][0].query_sequence),''.join([chr(qs+33) for qs in dataset[ID][0].query_qualities][::-1])))
            else:
                outfile_fwd.write('@%s\n%s\n+\n%s\n'%(dataset[ID][0].query_name,dataset[ID][0].query_sequence,''.join([chr(qs+33) for qs in dataset[ID][0].query_qualities])))
                outfile_rev.write('@%s\n%s\n+\n%s\n'%(dataset[ID][1].query_name,reverse_complement(dataset[ID][1].query_sequence),''.join([chr(qs+33) for qs in dataset[ID][1].query_qualities][::-1])))

    m,s=divmod(time.time()-start_time_check,60)
    h,m=divmod(m,60)
    print(F"{written_count:,} non-merged reads directly written to file (runtime {int(h)}:{int(m)}:{int(s)}).",flush=True)

    # handle paired reads
    # split remaining pairs in equal batches and distribute to merge workers

    batch_size=int(math.floor(len(merge_IDs)/cores))
    temp_cores=cores
    if batch_size<1:
        batch_size=1
        temp_cores=len(merge_IDs)

    q=Queue()

    thread_list=[]
    batch_start=0
    temp_files=[]
    while batch_start<len(merge_IDs):
        thread_list.append((batch_start,Process(target=merge_read_worker3,args=(batch_start,[dataset[ID] if ctr else [dataset[ID][1],dataset[ID][0]] for ID,ctr in merge_IDs[batch_start:batch_start+batch_size]],q))))
        thread_list[-1][1].deamon=True
        thread_list[-1][1].start()
        batch_start+=batch_size

    del merge_IDs
    del dataset

    while len(thread_list)>0:# have one writing function
        if not q.empty():
            # ID,results=q.get()
            ID,temp_f=q.get()
            temp_files.append(temp_f)
            # for mread in results:
            #     written_count+=1
            #     outfile_single.write('@%s\n%s\n+\n%s\n'%(mread.query_name,mread.query_sequence,''.join([chr(qs+33) for qs in mread.query_qualities])))
            #     if written_count%100000==0: print(F"{written_count:,} total reads written to file (thread queue length: {len(thread_list)})",flush=True)
            index=[t[0] for t in thread_list].index(ID)
            print(F"Merging worker ID {ID}, queue length: {len(thread_list)}",flush=True)
            thread_list[index][1].terminate()
            thread_list[index][1].join()
            thread_list.pop(index)

    outfile_single.close()
    outfile_fwd.close()
    outfile_rev.close()

    for datei in temp_files:
        command=F"gzip -4 -c {datei} >> {filename_prefix}.single.fq.gz"
        print(command,flush=True)
        proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while proc.poll() is None:#similar to wait until program terminates
            l=proc.stdout.readline().strip()
            if len(l)>0: print(l,flush=True)
        l=proc.stdout.read()
        proc.wait()
        if len(l)>0: print(l,flush=True)
        os.remove(datei)
        print(F"Deleted temp file {datei}",flush=True)

    #return alle 3 file names
    return([f'{filename_prefix}.single.fq.gz',f'{filename_prefix}.fwd.fq.gz',f'{filename_prefix}.rev.fq.gz'])

def manager(bamfile,cores,outfile_prefix,fake=False,SE=False):
    if not os.path.exists(bamfile):
        sys.exit('Input file "%s" does not exist. Aborting.'%bamfile)

    if not fake:
        merging_q=Queue()
        writing_q=Queue()
        process_ctr=Queue()

        reading_process=Process(target=reading_bamfile,args=(bamfile,SE,merging_q,writing_q,process_ctr,cores))
        reading_process.deamon=True
        reading_process.start()

        writing_process=Process(target=write_file,args=(outfile_prefix,writing_q,1+cores,process_ctr))
        writing_process.deamon=True
        writing_process.start()

        thread_list=[]
        for i in range(cores):
            thread_list.append(Process(target=merge_read_worker,args=(merging_q,writing_q,i+1)))
            thread_list[-1].deamon=True
            thread_list[-1].start()

        waiting=True
        read_ctr=True
        while waiting:
            if not process_ctr.empty():
                pID=process_ctr.get()
                if pID=='reading':
                    print('Finished reading.')
                    read_ctr=False
                elif pID=='writing':
                    print('Closing writing stream.')
                    waiting=False

    return([f'{outfile_prefix}.single.fq.gz',f'{outfile_prefix}.fwd.fq.gz',f'{outfile_prefix}.rev.fq.gz'])

def read_artefact_removal(in_filename,out_filename,cores):
    seen=set()
    with gzip.open(out_filename,'wt') as outfile:
        print(F"Processing {in_filename} for artefact removal",flush=True)
        with gzip.open(in_filename,'rt') as infile:
            title=''
            seq=''
            qs=''
            plus=False
            for line in infile:
                line=line.strip()
                if line.startswith('@') and title=='': title=line[1:]
                elif line.startswith('@') and title!='' and len(qs)==len(seq):
                    if len(seq)==len(qs):
                        if not title in seen:
                            seen.add(title)
                            outfile.write(F'@{title}\n{seq}\n+\n{qs}\n')
                    else:
                        sys.exit('error, inconsistent lengths')
                    title=line[1:]
                    seq=''
                    qs=''
                    plus=False
                elif line=='+': plus=True
                elif plus: qs+=line
                else: seq+=line
            if len(seq)==len(qs):
                if not title in seen:
                    seen.add(title)
                    outfile.write(F'@{title}\n{seq}\n+\n{qs}\n')
    return(True)

def UMI_xtr_worker(command,ID,q):
    proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while proc.poll() is None:#similar to wait until program terminates
        l=proc.stdout.readline().strip()
        if len(l)>0: print(l,flush=True)
    l=proc.stdout.read()
    proc.wait()
    if len(l)>0: print(l,flush=True)
    q.put(ID)
    return(True)

def ident_file_type(datei):
    datei=datei.lower()
    if datei.endswith('fa') or datei.endswith('fas') or datei.endswith('fasta'): return('fa')
    elif datei.endswith('fa.gz') or datei.endswith('fas.gz') or datei.endswith('fasta.gz'): return('fagz')
    elif datei.endswith('fq') or datei.endswith('fastq') or datei.endswith('fasta'): return('fq')
    elif datei.endswith('fq.gz') or datei.endswith('fastq.gz') or datei.endswith('fastq.gz'): return('fqgz')
    elif datei.endswith('.bam'): return('bam')
    else: return('unknown')

def read_file(datei,ID,q):
    # stepping=1000000
    seqs_count=0
    read_length=[]
    file_type=ident_file_type(datei)

    if file_type=='fa' or file_type=='fagz':
        seq=''
        title=''
        if file_type=='fa': infile=open(datei,'r')
        else: infile=gzip.open(datei,'rt')
        for line in infile:
            line=line.strip('\r\n').strip('\n')
            if line.startswith('>') and len(seq)>0:
                seqs_count+=1
                # read_length.append(len(seq))
                seq=''
                title=line[1:]
                # if seqs_count%stepping==0:
                #         print(f"{seqs_count:,} seqs read for {datei.split('/')[-1]}")
            elif line.startswith('>') and len(seq)==0: title=line[1:]
            else: seq+=line
        infile.close()
        seqs_count+=1
        # read_length.append(len(seq))

    elif file_type=='fq' or file_type=='fqgz':
        seq=''
        qs=''
        title=''
        plus=False
        if file_type=='fq':
          infile=open(datei,'r')
        else:
          infile=gzip.open(datei,'rt')
        for line in infile:
            line=line.strip('\r\n').strip('\n')
            if line.startswith('@') and len(seq)>0 and len(seq)==len(qs):
                seqs_count+=1
                # read_length.append(len(seq))
                seq=''
                qs=''
                title=line[1:]
                plus=False
                # if seqs_count%stepping==0:
                #         print(f"{seqs_count:,} seqs read for {datei.split('/')[-1]}")
            elif line.startswith('@') and len(seq)==0 and len(seq)==len(qs):
                title=line[1:]
            elif line.startswith('+') and len(line)==1: plus=True
            elif len(seq)>0 and plus: qs+=line
            else: seq+=line.upper()
        seqs_count+=1
        # read_length.append(len(seq))
        infile.close()
    elif file_type=='bam':
        try:
            bamfile=pysam.Samfile(datei,'rb')
            for read in bamfile.fetch(until_eof=True):
                if (not read.is_unmapped) and (not read.is_secondary):
                    seqs_count+=1
                    # read_length.append(len(read.query_sequence))
                # if seqs_count%stepping==0:
                #         print(f"{seqs_count:,} seqs read for {datei.split('/')[-1]}")
        except:
            print(f'Error in file {datei}',flush=True)
    else:
        print(f'Error! Unknown file type for {datei}',flush=True)
    q.put((seqs_count,ID))
    return True

def get_read_counts(track_files,attrition,cores,stage,filter=''):
    job_list=[]
    for seq_typ in track_files:
        for sampleID in track_files[seq_typ]:
            attrition[seq_typ][sampleID].append(0)
            if type(track_files[seq_typ][sampleID])==type([]):
                for datei in track_files[seq_typ][sampleID]:
                    if datei!=None and filter in datei:
                        job_list.append(((seq_typ,sampleID),datei))
            elif type(track_files[seq_typ][sampleID])==type(''):
                if track_files[seq_typ][sampleID]!=None and filter in track_files[seq_typ][sampleID]:
                    job_list.append(((seq_typ,sampleID),track_files[seq_typ][sampleID]))
            else:
                sys.exit('Error! Unknown file tracking structure in function "get_read_counts". Aborting.')

    q=Queue()
    thread_list=[]
    while len(job_list)>0 or len(thread_list)>0:
        while len(job_list)>0 and len(thread_list)<cores:
            ID,datei=job_list.pop()
            thread_list.append((ID,Process(target=read_file,args=(datei,ID,q))))
            thread_list[-1][1].deamon=True
            thread_list[-1][1].start()
        if not q.empty():
            count,ID=q.get()
            attrition[ID[0]][ID[1]][-1]+=count
            index=[x[0] for x in thread_list].index(ID)
            thread_list[index][1].join()
            thread_list.pop(index)

    attrition['stages'].append(stage)
    return(attrition)

def bam_to_fastq(infile,filename_prefix,fake):
    if not os.path.exists(infile):
        sys.exit('Input file "%s" does not exist. Aborting.'%bamfile)
    if not fake:
        outfile_single=gzip.open(filename_prefix+'.single.fq.gz','wt')
        outfile_fwd=gzip.open(filename_prefix+'.fwd.fq.gz','wt')
        outfile_rev=gzip.open(filename_prefix+'.rev.fq.gz','wt')
        bamfile=pysam.AlignmentFile(infile,'rb')
        count=0
        for read in bamfile.fetch(until_eof=True):
            count+=1
            if count%100000==0: print(F'Parsed {count:,} reads.')
            if (not read.is_unmapped) and (not read.is_secondary):
                outfile_single.write('@%s\n%s\n+\n%s\n'%(read.query_name,read.query_sequence,''.join([chr(qs+33) for qs in read.query_qualities])))
        outfile_single.close()
        outfile_fwd.close()
        outfile_rev.close()
        print(F'Finished writing of {count:,} reads for file prefix "{filename_prefix}".')
    return([f'{filename_prefix}.single.fq.gz',f'{filename_prefix}.fwd.fq.gz',f'{filename_prefix}.rev.fq.gz'])

def convert_read(read):
    obj=pysam.AlignedSegment()
    obj.query_sequence=read.query_sequence
    obj.query_name=read.query_name
    obj.query_qualities=read.query_qualities
    obj.reference_start=read.reference_start
    obj.is_read1=read.read1
    obj.is_reverse=read.is_reverse
    obj.is_unmapped=False
    obj.is_secondary=False
    obj.cigartuples=read.cigartuples
    obj.reference_id=0
    obj.mapping_quality=read.mapping_quality
    return(obj)

def merge_full(pair):
    obj=pysam.AlignedSegment()

    sequences=[]
    ref_range=[pair[0].reference_start,pair[1].reference_end]
    for ppos in range(2):
        sequences.append([])
        seq_temp=list(pair[ppos].query_sequence)
        qs_temp=list(pair[ppos].query_qualities)
        ref_pos=pair[ppos].reference_start

        for cig in pair[ppos].cigartuples:
            for subpos in range(cig[1]):
                if cig[0]==0 or cig[0]==7 or cig[0]==8:#match, seq match, seq mismatch
                    sequences[ppos].append([seq_temp.pop(0),qs_temp.pop(0),ref_pos,True])#base,QS,ref pos,read covering pos
                    ref_pos+=1
                elif cig[0]==4: #soft clip
                    seq_temp.pop(0)
                    qs_temp.pop(0)
                elif cig[0]==5 or cig[0]==6: #hard clip, silent padding
                    pass
                elif cig[0]==1:#insertion
                    sequences[ppos].append([seq_temp.pop(0),qs_temp.pop(0),None,True])
                elif cig[0]==2 or cig[0]==3:#deletion, skip ref
                    sequences[ppos].append([None,0,ref_pos,True])
                    ref_pos+=1

        #trimming accoring to outermost mapping coordinates
        left=0
        while sequences[ppos][left][2]==None or sequences[ppos][left][2]<ref_range[0]: left+=1
        right=len(sequences[ppos])-1
        while sequences[ppos][right][2]==None or sequences[ppos][right][2]>ref_range[1]:
            right-=1
        sequences[ppos]=sequences[ppos][left:right+1]

    # print(sequences[0])
    # print(sequences[1])

    aligned=[[],[]]

    #aligning sequences
    s1=''
    s2=''
    while len(sequences[0])>0 or len(sequences[1])>0:
        if s1=='' and len(sequences[0])>0: s1=sequences[0].pop(0)
        if s2=='' and len(sequences[1])>0: s2=sequences[1].pop(0)
        if s1!='' and s2!='' and s1[2]==s2[2]:
            aligned[0].append(s1)
            aligned[1].append(s2)
            s1=''
            s2=''
        elif (s1!='' and s1[2]==None) or s2=='':
            aligned[0].append(s1)
            aligned[1].append([None,0,None,False])
            s1=''
        elif (s2!='' and s2[2]==None) or s1=='':
            aligned[0].append([None,0,None,False])
            aligned[1].append(s2)
            s2=''
        elif s1[2]<s2[2]:
            aligned[0].append(s1)
            aligned[1].append([None,0,None,False])
            s1=''
        elif s1[2]>s2[2]:
            aligned[0].append([None,0,None,False])
            aligned[1].append(s2)
            s2=''
        else:
            print(s1,s2)
            sys.exit('Error while aligning read pairs to each other.')

    #read continuity
    index=len(aligned[0])-1
    while not aligned[0][index][3]: index-=1
    for i in range(0,index):
        aligned[0][i][3]=True

    index=0
    while not aligned[1][index][3]: index+=1
    for i in range(index+1,len(aligned[1])):
        aligned[1][i][3]=True

    consensus_seq=[]
    consensus_qs=[]
    consensus_merge_binary=[]

    consensus_cigar=[]
    cigar_status=None
    cigar_length=0
    ref_start=ref_range[0]
    ref_pos=ref_start

    alignment_center=int(math.ceil(len(aligned[0])/2))

    for pos in range(len(aligned[0])):
        #Seq and QS
        ref_liste=[aligned[i][pos] for i in range(2) if aligned[i][pos][3]]

        if len(ref_liste)==1:
            consensus_seq.append(ref_liste[0][0])
            consensus_qs.append(ref_liste[0][1])
        else:
            if ref_liste[0][0]==ref_liste[1][0]:#match
                consensus_seq.append(ref_liste[0][0])
                consensus_qs.append(ref_liste[0][1]+ref_liste[1][1])
            elif ref_liste[0][0]!=None and ref_liste[1][0]!=None:#mismatch, handling through QS
                if ref_liste[0][1]>ref_liste[1][1]:
                    consensus_seq.append(ref_liste[0][0])
                    consensus_qs.append(ref_liste[0][1]-ref_liste[1][1])
                elif ref_liste[0][1]<ref_liste[1][1]:
                    consensus_seq.append(ref_liste[1][0])
                    consensus_qs.append(-ref_liste[0][1]+ref_liste[1][1])
                elif ref_liste[0][1]==ref_liste[1][1]:
                    consensus_seq.append('N')
                    consensus_qs.append(1)
            else:#InDels
                if pos<=alignment_center:#trusting fwd read more, because closer to read start
                    consensus_seq.append(ref_liste[0][0])
                    consensus_qs.append(ref_liste[0][1])
                else:
                    consensus_seq.append(ref_liste[1][0])
                    consensus_qs.append(ref_liste[1][1])


        '''
        if aligned[0][pos][0]==aligned[1][pos][0]:
            consensus_seq.append(aligned[0][pos][0])
            consensus_qs.append(aligned[0][pos][1]+aligned[1][pos][1])
        else:
            if aligned[0][pos][0]!=None and aligned[1][pos][0]!=None:
                if aligned[0][pos][1]>aligned[1][pos][1]:
                    consensus_seq.append(aligned[0][pos][0])
                    consensus_qs.append(aligned[0][pos][1]-aligned[1][pos][1])
                elif aligned[0][pos][1]<aligned[1][pos][1]:
                    consensus_seq.append(aligned[1][pos][0])
                    consensus_qs.append(-aligned[0][pos][1]+aligned[1][pos][1])
                elif aligned[0][pos][1]==aligned[1][pos][1]:
                    consensus_seq.append('N')
                    consensus_qs.append(1)
            else:
                if aligned[0][pos][2]!=None: ref_pos_dist=aligned[0][pos][2]
                else: ref_pos_dist=aligned[1][pos][2]

                if aligned[0][pos][3] and not aligned[1][pos][3]:
                    consensus_seq.append(aligned[0][pos][0])
                    consensus_qs.append(aligned[0][pos][1])
                elif aligned[1][pos][3] and not aligned[0][pos][3]:
                    consensus_seq.append(aligned[1][pos][0])
                    consensus_qs.append(aligned[1][pos][1])
                elif ref_pos_dist!=None and abs(ref_pos_dist-ref_range[0])<=abs(ref_range[1]-ref_pos_dist):
                    consensus_seq.append(aligned[0][pos][0])
                    if aligned[0][pos][1]>0:
                        consensus_qs.append(aligned[0][pos][1])
                    else:
                        consensus_qs.append(aligned[1][pos][1])
                elif ref_pos_dist!=None and abs(ref_pos_dist-ref_range[0])>abs(ref_range[1]-ref_pos_dist):
                    consensus_seq.append(aligned[1][pos][0])
                    if aligned[1][pos][1]>0:
                        consensus_qs.append(aligned[1][pos][1])
                    else:
                        consensus_qs.append(aligned[0][pos][1])
                else:
                    sys.exit('Error')
        '''

        consensus_merge_binary.append(aligned[0][pos][3]*aligned[1][pos][3])

        #cigar
        ref_liste=[aligned[i][pos] for i in range(2) if aligned[i][pos][3]]
        if len(ref_liste)==1:
            if ref_liste[0][2]==None:#insertion
                temp_cig=1
                temp_cig_len=1
            elif ref_liste[0][0]==None:#deletion
                temp_cig=2
                temp_cig_len=1
            else:
                temp_cig=0
                temp_cig_len=1
        else:
            if consensus_seq[-1]==None:#deletion
                temp_cig=2
                temp_cig_len=1
            elif ref_liste[0][2]==None and ref_liste[1][2]==None:#insertion
                temp_cig=1
                temp_cig_len=1
            else:
                temp_cig=0
                temp_cig_len=1


        # #cigar
        # if consensus_seq[-1]!=None and (aligned[0][pos][2]==None or aligned[1][pos][2]==None):#insertion
        #     temp_cig=1
        #     temp_cig_len=1
        # # elif aligned[0][pos][2]>ref_pos+1:#deletion
        # #     temp_cig=2
        # #     temp_cig_len=aligned[0][pos][2]+1-ref_pos
        # elif (aligned[0][pos][2]!=None or aligned[1][pos][2]!=None) and consensus_seq[-1]==None:#deletion
        #     temp_cig=2
        #     temp_cig_len=1
        # else:
        #     temp_cig=0
        #     temp_cig_len=1
        #
        # if aligned[0][pos][2]!=None:
        #     ref_pos=aligned[0][pos][2]

        if cigar_status==None:
            cigar_status=temp_cig
            cigar_length=temp_cig_len
        else:
            if cigar_status==temp_cig:
                cigar_length+=temp_cig_len
            else:
                consensus_cigar+=[cigar_status]*cigar_length
                cigar_status=temp_cig
                cigar_length=temp_cig_len
    consensus_cigar+=[cigar_status]*cigar_length

    consensus_qs=[consensus_qs[pos] for pos in range(len(consensus_seq)) if consensus_seq[pos]!=None]
    consensus_merge_binary=[consensus_merge_binary[pos] for pos in range(len(consensus_seq)) if consensus_seq[pos]!=None]
    # do not purge cigar, or you will loose within read deletion information: consensus_cigar=[consensus_cigar[pos] for pos in range(len(consensus_seq)) if consensus_seq[pos]!=None]
    consensus_seq=[consensus_seq[pos] for pos in range(len(consensus_seq)) if consensus_seq[pos]!=None]

    temp=consensus_cigar


    consensus_merge_binary=['1' if c else '0' for c in consensus_merge_binary]
    consensus_cigar_temp=[]
    for cig in consensus_cigar:
        if len(consensus_cigar_temp)==0 or consensus_cigar_temp[-1][0]!=cig:
            consensus_cigar_temp.append([cig,1])
        else:
            consensus_cigar_temp[-1][1]+=1
    consensus_cigar=consensus_cigar_temp

    #softclipping
    for pos in [0,-1]:
        if consensus_cigar[pos][0]==1:  consensus_cigar[pos][0]=4

    consensus_qs=[qs if qs<=60 else 60 for qs in consensus_qs]#capping QS

    # if pair[0].query_name=='A00417:34:HG3HFDRXX:2:1260:14570:14857_GGTCGTTG':
    #     for pos in range(len(aligned[0])):
    #         print(pos,aligned[0][pos],aligned[1][pos],consensus_seq[pos],temp[pos])

    # if not len(consensus_seq)==sum([c[1] for c in consensus_cigar if c[0] in [0,1,4,7,8]]):
    #     if len(consensus_seq)<100:
    #         print('error')
    #         print(consensus_cigar,f"sum {sum([c[1] for c in consensus_cigar if c[0] in [0,1,4,7,8]])}")
    #         print(len(consensus_seq))
    #
    #         for pos in range(len(aligned[0])):
    #             print(pos,aligned[0][pos],aligned[1][pos],consensus_seq[pos],temp[pos])
    #         print(consensus_cigar)
    #         sys.exit()


    obj.query_sequence=''.join(consensus_seq)
    obj.query_name=pair[0].query_name+'_FR'+hex(int('1'+''.join(consensus_merge_binary),2))#to encode overlapping region
    obj.query_qualities=consensus_qs
    obj.reference_start=ref_start
    obj.cigartuples=consensus_cigar
    obj.mapping_quality=int(round((pair[0].mapping_quality+pair[1].mapping_quality)/2))

    obj.is_read1=True
    obj.is_reverse=False
    obj.is_unmapped=False
    obj.is_secondary=False
    obj.reference_id=0
    return(obj)

def template_reconstruction(infile,outfile,sID,temp_path):
    #read bam file in packages and filter for reference
    if not os.path.exists(infile):
        sys.exit('Input file "%s" does not exist. Aborting.'%infile)

    def merge_worker(pairs,pID,process_q,temp_path):
        bamfile_out=pysam.AlignmentFile(F"temp_{pID}a.bam",'wb',reference_names=[metrics['ref_ID']],reference_lengths=[len(metrics['ref_seq'])])
        for pair in pairs:
            bamfile_out.write(merge_full(pair))
        bamfile_out.close()

        command=F"{metrics['programs']['path_to_samtools']} sort -O bam -T {temp_path} -o temp_{pID}.bam temp_{pID}a.bam"
        proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while proc.poll() is None:#similar to wait until program terminates
            l=proc.stdout.readline().strip()
            if len(l)>0: print(l,flush=True)
        l=proc.stdout.read()
        if len(l)>0: print(l,flush=True)

        os.remove(F"temp_{pID}a.bam")

        process_q.put(pID)
        return(True)

    start_time_check=time.time()

    subcores=int(metrics['CPU_count']/2)
    if subcores<1: subcores=1

    merge_pairs=[]
    bamfile_in=pysam.AlignmentFile(infile,'rb',threads=subcores)
    bamfile_out=pysam.AlignmentFile('temp_0.bam','wb',reference_names=[metrics['ref_ID']],reference_lengths=[len(metrics['ref_seq'])],threads=subcores)

    #multiprocessing variables
    process_q=Queue()
    batch_size=100000
    thread_list=[]
    worker_ID=1
    join_list=['temp_0.bam']

    #reading_parameters
    r_count=0
    w_count=0
    buffer={}

    for read in bamfile_in.fetch(contig=metrics['ref_ID']):
        r_count+=1
        if r_count%200000==0: print(F'Parsed {r_count:,} reads from {infile.split("/")[-1]}. {w_count:,} reads directly written.',flush=True)

        if (not read.is_unmapped) and (not read.is_secondary):
            if read.mate_is_unmapped:#singleton, nothing to merge
                w_count+=1
                #####
                # REACTIVATE
                #####
                bamfile_out.write(read)
            else:
                if not read.query_name in buffer:
                    buffer[read.query_name]=read_obj(read)
                else:
                    if read.reference_start<=buffer[read.query_name].reference_start<read.reference_end or read.reference_start<buffer[read.query_name].reference_end<=read.reference_end or (read.reference_start==buffer[read.query_name].reference_start and buffer[read.query_name].reference_end==read.reference_end):#overlap
                        #making sure, that the first read is in fwd orientation
                        if buffer[read.query_name].is_reverse:
                            merge_pairs.append((read_obj(read),buffer.pop(read.query_name)))
                            #####
                            # DELETE LINES BELOW
                            #####
                            # bamfile_out.write(convert_read(merge_pairs[-1][0]))
                            # bamfile_out.write(convert_read(merge_pairs[-1][1]))
                        else:
                            merge_pairs.append((buffer.pop(read.query_name),read_obj(read)))
                            #####
                            # DELETE LINES BELOW
                            #####
                            # bamfile_out.write(convert_read(merge_pairs[-1][0]))
                            # bamfile_out.write(convert_read(merge_pairs[-1][1]))
                    else:
                        #####
                        # REACTIVATE
                        #####
                        bamfile_out.write(read)
                        bamfile_out.write(convert_read(buffer.pop(read.query_name)))
                        w_count+=2

        if len(merge_pairs)>=batch_size and len(thread_list)<metrics['CPU_count']:
            thread_list.append((worker_ID,Process(target=merge_worker,args=(merge_pairs[:batch_size],worker_ID,process_q,temp_path))))
            thread_list[-1][1].deamon=True
            thread_list[-1][1].start()
            merge_pairs=merge_pairs[batch_size:]
            worker_ID+=1

        if not process_q.empty():
            pID=process_q.get()
            join_list.append(F"temp_{pID}.bam")
            pindex=[t[0] for t in thread_list].index(pID)
            thread_list[pindex][1].join(1)
            thread_list[pindex][1].terminate()
            # thread_list[pindex][1].close()
            thread_list.pop(pindex)
            gc.collect()


        # if r_count>10000:
        #     print("\n\nREMOVE LINE BELOW\n\n",flush=True)
        #     break

    #####
    # REACTIVATE
    #####
    for read in buffer.values():#to empty buffer
        w_count+=1
        bamfile_out.write(convert_read(read))

    print(F'Parsed {r_count:,} reads from {infile.split("/")[-1]}. {w_count:,} reads directly written.',flush=True)

    del buffer
    gc.collect()

    bamfile_in.close()
    bamfile_out.close()

    while len(merge_pairs)>0 or len(thread_list)>0:
        while len(thread_list)<metrics['CPU_count'] and len(merge_pairs)>0:
            thread_list.append((worker_ID,Process(target=merge_worker,args=(merge_pairs[:batch_size],worker_ID,process_q,temp_path))))
            thread_list[-1][1].deamon=True
            thread_list[-1][1].start()
            merge_pairs=merge_pairs[batch_size:]
            worker_ID+=1

        if not process_q.empty():
            pID=process_q.get()
            join_list.append(F"temp_{pID}.bam")
            pindex=[t[0] for t in thread_list].index(pID)
            thread_list[pindex][1].join(1)
            thread_list[pindex][1].terminate()
            thread_list.pop(pindex)
            gc.collect()
            print(F"Workers: {len(thread_list)}; remaining read pairs to merge: {len(merge_pairs):,}",flush=True)

    process_q.close()

    del merge_pairs
    gc.collect()

    #join files and remove temp files
    command=F"{metrics['programs']['path_to_samtools']} merge -f -@ {metrics['CPU_count']} - {' '.join(join_list)} | {metrics['programs']['path_to_samtools']} sort -@ {metrics['CPU_count']} -O bam -T {temp_path} -o {outfile} - && {metrics['programs']['path_to_samtools']} index {outfile}"
    print(command)
    proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while proc.poll() is None:#similar to wait until program terminates
        l=proc.stdout.readline().strip()
        if len(l)>0: print(l,flush=True)
    l=proc.stdout.read()
    if len(l)>0: print(l,flush=True)

    #remove temp files
    for datei in join_list:
        os.remove(datei)

    m,s=divmod(time.time()-start_time_check,60)
    h,m=divmod(m,60)
    print(F"Runtime {int(h):02d}:{int(m):02d}:{int(s):02d}.",flush=True)

    # sys.exit('checkpoint')

    return(outfile)

def mapping_alternative(metrics,path,mapper='bwa'):#replicating mapping approach
    if not path[-1]=='/': path+='/'

    track_files={'DNA':{},'RNA':{}}
    track_PE_protocol={'DNA':{},'RNA':{}}
    attrition={'DNA':{},'RNA':{},'stages':[]}

    #filling in sample structure
    for seq_typ in ['DNA','RNA']:
        for ID,direction,location in metrics['input_files'][seq_typ]:
            if not ID in track_files[seq_typ]:
                track_files[seq_typ][ID]=[None,None]
                attrition[seq_typ][ID]=[]
                track_PE_protocol[seq_typ][ID]=[False,False]
            if direction=='fwd':
                track_files[seq_typ][ID][0]=location
                track_PE_protocol[seq_typ][ID][0]=True
            elif direction=='rev':
                track_files[seq_typ][ID][1]=location
                track_PE_protocol[seq_typ][ID][1]=True
            else:
                sys.exit('Unknown read direction typ')

    #keeping track of PE vs SE seq protocol
    for seq_typ in track_PE_protocol:
        for ID in track_PE_protocol[seq_typ]:
            track_PE_protocol[seq_typ][ID]=track_PE_protocol[seq_typ][ID][0]*track_PE_protocol[seq_typ][ID][1]

    #UMI extraction
    if not options.no_umi:
        stage='UMI extraction'
        print(stage)
        if not os.path.exists(path+'00_umi_xtr'): os.makedirs(path+'00_umi_xtr')
        seq_typ='RNA'
        # starttime=time.time()
        # counter=0
        # counter+=1
        # print(f"est. remaining runtime {round(((time.time()-starttime)/counter)*(len(track_files[seq_typ])-counter)/60):.2f,} min")
        thread_list=[]
        commands=[]
        q=Queue()
        for sampleID in track_files[seq_typ]:
            print(f"starting UMI extraction for {sampleID}")
            if track_files[seq_typ][sampleID][0]!=None:
                new_file_fwd=path+'00_umi_xtr/'+track_files[seq_typ][sampleID][0].split('/')[-1]
            else:
                new_file_fwd=None
            if track_files[seq_typ][sampleID][1]!=None:
                new_file_rev=path+'00_umi_xtr/'+track_files[seq_typ][sampleID][1].split('/')[-1]
            else:
                new_file_rev=None

            if track_files[seq_typ][sampleID][0]!=None:
                if track_files[seq_typ][sampleID][1]!=None:
                    command=F"{metrics['programs']['path_to_umi_tools']} extract --bc-pattern=NNNN --bc-pattern2=NNNN --stdin={track_files[seq_typ][sampleID][0]} --read2-in={track_files[seq_typ][sampleID][1]} --stdout={new_file_fwd} --read2-out={new_file_rev} --log={metrics['log']} >> {metrics['log']}"
                    print(command)
                else:
                    command=F"{metrics['programs']['path_to_umi_tools']} extract --bc-pattern=NNNN --stdin={track_files[seq_typ][sampleID][0]} --stdout={new_file_fwd} --log={metrics['log']} >> {metrics['log']}"
            elif track_files[seq_typ][sampleID][1]!=None:
                command=F"{metrics['programs']['path_to_umi_tools']} extract --bc-pattern=NNNN --stdin={track_files[seq_typ][sampleID][1]} --stdout={new_file_rev} --log={metrics['log']} >> {metrics['log']}"
            commands.append(command)
            track_files[seq_typ][sampleID]=[new_file_fwd,new_file_rev]

        #running tool in parallel
        counter=0
        if metrics['skip'][0]=='1':
            while len(commands)>0 or len(thread_list)>0:
                if len(commands)>0 and len(thread_list)<metrics['CPU_count']:
                    thread_list.append((counter,Process(target=UMI_xtr_worker,args=(commands.pop(),counter,q))))
                    thread_list[-1][1].deamon=True
                    thread_list[-1][1].start()
                    counter+=1
                if not q.empty():
                    ID=q.get()
                    index=[x[0] for x in thread_list].index(ID)
                    thread_list[index][1].join()
                    thread_list.pop(index)
        print(f'Finished {stage}')

        if options.attrition:
            print(f'collecting read counts for stage {stage}')
            attrition=get_read_counts(track_files,attrition,metrics['CPU_count'],stage)
            if options.attrition:
                create_attrition_plot(attrition,path+'%sattrition.svg'%metrics['prefix'])

    #cutadapt clipping
    stage='Adapter clipping'
    print(stage)
    if not os.path.exists(path+'01_cutadapt'): os.makedirs(path+'01_cutadapt')
    for seq_typ in ['DNA','RNA']:
        for sampleID in track_files[seq_typ]:
            print(f"processing file {sampleID}")
            if track_files[seq_typ][sampleID][0]!=None:
                new_file_fwd=path+'01_cutadapt/'+track_files[seq_typ][sampleID][0].split('/')[-1]
            else:
                new_file_fwd=None
            if track_files[seq_typ][sampleID][1]!=None:
                new_file_rev=path+'01_cutadapt/'+track_files[seq_typ][sampleID][1].split('/')[-1]
            else:
                new_file_rev=None

            if new_file_fwd!=None:
                if new_file_rev!=None:
                    command=f"{metrics['programs']['path_to_cutadapt']} {' '.join(['-a %s'%s for s in metrics['references']['cutadapt']['fwd']])} {' '.join(['-A %s'%s for s in metrics['references']['cutadapt']['rev']])} -j {metrics['CPU_count']} -q 10 -m 22 --overlap=3 -o {new_file_fwd} -p {new_file_rev} {track_files[seq_typ][sampleID][0]} {track_files[seq_typ][sampleID][1]} >> {metrics['log']}"
                else:
                    command=f"{metrics['programs']['path_to_cutadapt']} {' '.join(['-a %s'%s for s in metrics['references']['cutadapt']['fwd']])} -j {metrics['CPU_count']} -q 10 -m 22 --overlap=3 -o {new_file_fwd} {track_files[seq_typ][sampleID][0]} >> {metrics['log']}"
            elif new_file_rev!=None:
                command=f"{metrics['programs']['path_to_cutadapt']} {' '.join(['-a %s'%s for s in metrics['references']['cutadapt']['rev']])} -j {metrics['CPU_count']} -q 10 -m 22 --overlap=3 -o {new_file_rev} {track_files[seq_typ][sampleID][1]} >> {metrics['log']}"

            if metrics['skip'][1]=='1':
                print(command)
                proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                while proc.poll() is None:#similar to wait until program terminates
                    l=proc.stdout.readline().strip()
                    if len(l)>0: print(l,flush=True)
                l=proc.stdout.read()
                if len(l)>0: print(l,flush=True)

            track_files[seq_typ][sampleID]=[new_file_fwd,new_file_rev]

    if options.attrition:
        print(f'collecting read counts for stage {stage}')
        attrition=get_read_counts(track_files,attrition,metrics['CPU_count'],stage)
        if options.attrition:
            create_attrition_plot(attrition,path+'%sattrition.svg'%metrics['prefix'])

    if metrics['skip'][2]=='1':
        if not os.path.exists(path+'mapping_ref'): os.makedirs(path+'mapping_ref')
        #star mapping index
        print('Creating mapping reference')
        #combining ref fasta with phix fasta
        with open(path+'mapping_ref/mapping_ref.fas','w') as outfile:
            for datei in [metrics['ref_file'],metrics['PhiX']]:
                seq,ID=read_fasta(datei)
                outfile.write(F">{ID}\n{seq}\n")
        # command=F"cat {metrics['ref_file']} {metrics['PhiX']} > {path+'mapping_ref/mapping_ref.fas'}"
        # proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # while proc.poll() is None:#similar to wait until program terminates
        #     l=proc.stdout.readline().strip()
        #     if len(l)>0: print(l,flush=True)
        # l=proc.stdout.read()
        # if len(l)>0: print(l,flush=True)

        if mapper=='star':
            command=F"{metrics['programs']['path_to_star']} --runThreadN {metrics['CPU_count']} --runMode genomeGenerate  --genomeSAindexNbases {min(14,int(round(math.log(len(metrics['ref_seq']),2)/2-1)))} --genomeDir {path+'mapping_ref'} --genomeFastaFiles {path+'mapping_ref/mapping_ref.fas'}"
        elif mapper=='bwa':
            command=F"{metrics['programs']['path_to_bwa']} index {path+'mapping_ref/mapping_ref.fas'}"
        proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while proc.poll() is None:#similar to wait until program terminates
            l=proc.stdout.readline().strip()
            if len(l)>0: print(l,flush=True)
        l=proc.stdout.read()
        if len(l)>0: print(l,flush=True)

    #read mapping
    stage='read mapping'
    print(stage)
    if not os.path.exists(path+'02_mappings'): os.makedirs(path+'02_mappings')
    for seq_typ in ['DNA','RNA']:
        for sampleID in track_files[seq_typ]:
            print(f"processing sample {sampleID}")
            if mapper=='star':
                if track_files[seq_typ][sampleID][0]!=None:
                    if track_files[seq_typ][sampleID][1]!=None:
                        command=F"{metrics['programs']['path_to_star']} --runThreadN {metrics['CPU_count']} --genomeDir {path+'mapping_ref'} --outSAMtype BAM Unsorted --readFilesIn {track_files[seq_typ][sampleID][0]} {track_files[seq_typ][sampleID][1]} --readFilesCommand zcat --outSAMattributes All --outFilterType Normal --outMultimapperOrder Random --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFileNamePrefix {path+'02_mappings/'+sampleID+'_'+seq_typ+'_'}"
                    else:
                        command=F"{metrics['programs']['path_to_star']} --runThreadN {metrics['CPU_count']} --genomeDir {path+'mapping_ref'} --outSAMtype BAM Unsorted --readFilesIn {track_files[seq_typ][sampleID][0]} --readFilesCommand zcat --outSAMattributes All --outFilterType Normal --outMultimapperOrder Random --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFileNamePrefix {path+'02_mappings/'+sampleID+'_'+seq_typ+'_'}"
                elif track_files[seq_typ][sampleID][1]!=None:
                    command=F"{metrics['programs']['path_to_star']} --runThreadN {metrics['CPU_count']} --genomeDir {path+'mapping_ref'} --outSAMtype BAM Unsorted --readFilesIn {track_files[seq_typ][sampleID][1]} --readFilesCommand zcat --outSAMattributes All --outFilterType Normal --outMultimapperOrder Random --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFileNamePrefix {path+'02_mappings/'+sampleID+'_'+seq_typ+'_'}"

            elif mapper=='bwa':
                if track_files[seq_typ][sampleID][0]!=None:
                    if track_files[seq_typ][sampleID][1]!=None:
                        command=F'{metrics["programs"]["path_to_bwa"]} mem -t {metrics["CPU_count"]} -R "@RG\\tID:{sampleID}\\tLB:lib1\\tPL:illumina\\tSM:None" {path+"mapping_ref/mapping_ref.fas"} {track_files[seq_typ][sampleID][0]} {track_files[seq_typ][sampleID][1]} | {metrics["programs"]["path_to_samtools"]} view -b -o {path}02_mappings/{sampleID}_{seq_typ}_Aligned.out.bam --threads {metrics["CPU_count"]} -'
                    else:
                        command=F'{metrics["programs"]["path_to_bwa"]} mem -t {metrics["CPU_count"]} -R "@RG\\tID:{sampleID}\\tLB:lib1\\tPL:illumina\\tSM:None" {path+"mapping_ref/mapping_ref.fas"} {track_files[seq_typ][sampleID][0]} | {metrics["programs"]["path_to_samtools"]} view -b -o {path}02_mappings/{sampleID}_{seq_typ}_Aligned.out.bam --threads {metrics["CPU_count"]} -'
                elif track_files[seq_typ][sampleID][1]!=None:
                    command=F'{metrics["programs"]["path_to_bwa"]} mem -t {metrics["CPU_count"]} -R "@RG\\tID:{sampleID}\\tLB:lib1\\tPL:illumina\\tSM:None" {path+"mapping_ref/mapping_ref.fas"} {track_files[seq_typ][sampleID][1]} | {metrics["programs"]["path_to_samtools"]} view -b -o {path}02_mappings/{sampleID}_{seq_typ}_Aligned.out.bam --threads {metrics["CPU_count"]} -'

            if metrics['skip'][3]=='1':
                proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                while proc.poll() is None:#similar to wait until program terminates
                    l=proc.stdout.readline().strip()
                    if len(l)>0: print(l,flush=True)
                l=proc.stdout.read()
                if len(l)>0: print(l,flush=True)
            track_files[seq_typ][sampleID]=f"{path}02_mappings/{sampleID}_{seq_typ}_Aligned.out.bam"

    if options.attrition:
        print(f'collecting read counts for stage {stage}')
        attrition=get_read_counts(track_files,attrition,metrics['CPU_count'],stage,filter="_Aligned.out.bam")
        if options.attrition:
            create_attrition_plot(attrition,path+'%sattrition.svg'%metrics['prefix'])

    #sorting
    stage='Sorting alignments'
    print(stage)
    if not os.path.exists(path+'03_sorted'): os.makedirs(path+'03_sorted')
    for seq_typ in ['DNA','RNA']:
        for sampleID in track_files[seq_typ]:
            print(f"processing sample {sampleID}")
            new_file=f"{path}03_sorted/{sampleID}_{seq_typ}_sort.bam"
            command=F"{metrics['programs']['path_to_samtools']} sort -@ {metrics['CPU_count']} -O bam -T {path+'03_sorted/tmp'} -o {new_file} {track_files[seq_typ][sampleID]} >> {metrics['log']}"
            if metrics['skip'][4]=='1':
                proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                while proc.poll() is None:#similar to wait until program terminates
                    l=proc.stdout.readline().strip()
                    if len(l)>0: print(l,flush=True)
                l=proc.stdout.read()
                if len(l)>0: print(l,flush=True)
            track_files[seq_typ][sampleID]=new_file

    if options.attrition:
        print(f'collecting read counts for stage {stage}')
        attrition=get_read_counts(track_files,attrition,metrics['CPU_count'],stage)
        if options.attrition:
            create_attrition_plot(attrition,path+'%sattrition.svg'%metrics['prefix'])

    #indexing
    print('bam file indexing')
    for seq_typ in ['DNA','RNA']:
        for sampleID in track_files[seq_typ]:
            command=F"{metrics['programs']['path_to_samtools']} index {track_files[seq_typ][sampleID]} >> {metrics['log']}"
            if metrics['skip'][5]=='1':
                proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                while proc.poll() is None:#similar to wait until program terminates
                    l=proc.stdout.readline().strip()
                    if len(l)>0: print(l,flush=True)
                l=proc.stdout.read()
                if len(l)>0: print(l,flush=True)

    #deduplication (RNA: UMI based, DNA: samtools based)
    stage='Read deduplication'
    print(stage)
    if not os.path.exists(path+'04_deduplication'): os.makedirs(path+'04_deduplication')
    commands_list=[]
    for seq_typ in ['DNA','RNA']:
        if seq_typ=='DNA' or options.no_umi:
            print(F'Deduplication of {seq_typ} data using samtools engine')
            for sampleID in track_files[seq_typ]:
                print(f"processing sample {sampleID}")
                new_file=f"{path}04_deduplication/{sampleID}_{seq_typ}_dedup.bam"
                command=f"{metrics['programs']['path_to_samtools']} rmdup --reference {path+'mapping_ref/mapping_ref.fas'} {track_files[seq_typ][sampleID]} {new_file} >>{metrics['log']} && {metrics['programs']['path_to_samtools']} index {new_file} >>{metrics['log']}"
                commands_list.append(command)
                track_files[seq_typ][sampleID]=new_file
                # print(command)
        else:
            print(F'Deduplication of {seq_typ} data using UMIs')
            for sampleID in track_files[seq_typ]:
                print(f"processing sample {sampleID}")
                new_file=f"{path}04_deduplication/{sampleID}_{seq_typ}_dedup.bam"
                command=F"{metrics['programs']['path_to_umi_tools']} dedup --edit-distance-threshold=0 --paired --stdin={track_files[seq_typ][sampleID]} --output-stats={metrics['log']} --log={metrics['log']} | {metrics['programs']['path_to_samtools']} sort -@ {metrics['CPU_count']} -O bam -T {path+'04_deduplication/tmp'} -o {new_file} >> {metrics['log']} && {metrics['programs']['path_to_samtools']} index {new_file} >> {metrics['log']}"
                commands_list.append(command)
                track_files[seq_typ][sampleID]=new_file

    if metrics['skip'][6]=='1':
        for command in commands_list:
            proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            while proc.poll() is None:#similar to wait until program terminates
                l=proc.stdout.readline().strip()
                if len(l)>0: print(l,flush=True)
            l=proc.stdout.read()
            if len(l)>0: print(l,flush=True)

    if options.attrition:
        print(f'collecting read counts for stage {stage}')
        attrition=get_read_counts(track_files,attrition,metrics['CPU_count'],stage)
        if options.attrition:
            create_attrition_plot(attrition,path+'%sattrition.svg'%metrics['prefix'])

    stage='PCR_template_reconstruction'
    print(stage,flush=True)

    if not os.path.exists(path+'05_template_reconstruction'): os.makedirs(path+'05_template_reconstruction')


    for seq_typ in ['DNA','RNA']:
        for sampleID in track_files[seq_typ]:
            print(f"processing sample {sampleID} for {seq_typ}")
            new_file=f"{path}05_template_reconstruction/{sampleID}_{seq_typ}_sort.bam"
            if metrics['skip'][7]=='1':
                template_reconstruction(track_files[seq_typ][sampleID],new_file,sampleID,path+'05_template_reconstruction/')
            track_files[seq_typ][sampleID]=new_file

    if options.attrition:
        print(f'collecting read counts for stage {stage}')
        attrition=get_read_counts(track_files,attrition,metrics['CPU_count'],stage)
        if options.attrition:
            create_attrition_plot(attrition,path+'%sattrition.svg'%metrics['prefix'])

    # sys.exit('checkpoint_stage')


    # #demultiplexing references
    # stage='reference demultiplexing'
    # print(stage)
    # for seq_typ in ['DNA','RNA']:
    #     for sampleID in track_files[seq_typ]:
    #         print(f"processing sample {sampleID}")
    #         new_file=f"{path}04_deduplication/{sampleID}_{seq_typ}_dedup_clean.bam"
    #         command=F"{metrics['programs']['path_to_samtools']} view -@ {metrics['CPU_count']} -b {track_files[seq_typ][sampleID]} {metrics['ref_ID']} | {metrics['programs']['path_to_samtools']} sort -@ {metrics['CPU_count']} > {new_file} 2>> {metrics['log']} && {metrics['programs']['path_to_samtools']} index {new_file} 2>> {metrics['log']}"
    #         if metrics['skip'][7]=='1':
    #             proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #             while proc.poll() is None:#similar to wait until program terminates
    #                 l=proc.stdout.readline().strip()
    #                 if len(l)>0: print(l,flush=True)
    #             l=proc.stdout.read()
    #             if len(l)>0: print(l,flush=True)
    #         track_files[seq_typ][sampleID]=new_file
    #
    # if options.terminate==8: sys.exit('Terminating on user request!')
    #
    # if options.attrition:
    #     print(f'collecting read counts for stage {stage}')
    #     attrition=get_read_counts(track_files,attrition,metrics['CPU_count'],stage,filter='_dedup_clean.bam')
    #     if options.attrition:
    #         create_attrition_plot(attrition,path+'%sattrition.svg'%metrics['prefix'])
    #
    # #read merging
    # stage='reads merging'
    # print(stage)
    # start_time_check=time.time()
    # if not os.path.exists(path+'05_merging'): os.makedirs(path+'05_merging')
    #
    # skip_existing_files=[e[1] for e in sorted([(os.path.getmtime(f"{path}05_merging/{d}"),d) for d in os.listdir(f"{path}05_merging/") if d.endswith('single.fq.gz')])[:-1]]
    #
    # for seq_typ in ['RNA','DNA']:
    #     for sampleID in track_files[seq_typ]:
    #         print(f"processing sample {sampleID} for {seq_typ}")
    #         new_file_prefix=f"{path}05_merging/{sampleID}_{seq_typ}"
    #         fake=True
    #         outfiles=[f'{new_file_prefix}.single.fq.gz',f'{new_file_prefix}.fwd.fq.gz',f'{new_file_prefix}.rev.fq.gz']
    #         if metrics['skip'][8]=='1':
    #             fake=False
    #             SE=not track_PE_protocol[seq_typ][sampleID]
    #             if not SE:
    #                 if (options.auto_skip and not F"{sampleID}_{seq_typ}.single.fq.gz" in skip_existing_files) or (not options.auto_skip):
    #                     # outfiles=merging_manager2(track_files[seq_typ][sampleID],metrics['CPU_count'],new_file_prefix)
    #                     outfiles=manager(track_files[seq_typ][sampleID],metrics['CPU_count'],new_file_prefix,fake=fake,SE=SE)
    #             else:
    #                 if (options.auto_skip and not F"{sampleID}_{seq_typ}.single.fq.gz" in skip_existing_files) or (not options.auto_skip):
    #                     outfiles=bam_to_fastq(track_files[seq_typ][sampleID],new_file_prefix,fake=fake)
    #         track_files[seq_typ][sampleID]=outfiles
    #
    # m,s=divmod(time.time()-start_time_check,60)
    # h,m=divmod(m,60)
    # print(F"Total runtime {int(h)}:{int(m)}:{int(s)}",flush=True)
    # # print('previous runtime (old,1 core): 3:28:50')
    # # print('previous runtime (old, 80 core): 2:25:5')
    # # print('previous runtime (new, 80 cores): 1:36:44')
    # # sys.exit('checkpoint')
    #
    #
    # if options.attrition:
    #     print(f'collecting read counts for stage {stage}')
    #     attrition=get_read_counts(track_files,attrition,metrics['CPU_count'],stage)
    #     if options.attrition:
    #         create_attrition_plot(attrition,path+'%sattrition.svg'%metrics['prefix'])
    #
    # #removing read duplicate artefacts
    # stage='removing mapping duplication'
    # print(stage)
    # if not os.path.exists(path+'05a_read_dup_artefacts'): os.makedirs(path+'05a_read_dup_artefacts')
    #
    # skip_existing_files=[e[1] for e in sorted([(os.path.getmtime(f"{path}05a_read_dup_artefacts/{d}"),d) for d in os.listdir(f"{path}05a_read_dup_artefacts/") if d.endswith('.fq.gz')])[:-1]]
    #
    # for seq_typ in ['DNA','RNA']:
    #     for sampleID in track_files[seq_typ]:
    #         print(f"processing sample {sampleID}")
    #         for file_pos in range(3):
    #             new_file=path+'05a_read_dup_artefacts/'+track_files[seq_typ][sampleID][file_pos].split('/')[-1]
    #             if metrics['skip'][9]=='1' and ((options.auto_skip and not new_file.split('/')[-1] in skip_existing_files) or (not options.auto_skip)):
    #                 read_artefact_removal(track_files[seq_typ][sampleID][file_pos],new_file,metrics['CPU_count'])
    #             track_files[seq_typ][sampleID][file_pos]=new_file
    #
    # if options.attrition:
    #     print(f'collecting read counts for stage {stage}')
    #     attrition=get_read_counts(track_files,attrition,metrics['CPU_count'],stage)
    #     if options.attrition:
    #         create_attrition_plot(attrition,path+'%sattrition.svg'%metrics['prefix'])
    #
    # #creating single reference mapping reference
    # print('Creating PURE mapping reference')
    # if not os.path.exists(path+'mapping_ref_pure'): os.makedirs(path+'mapping_ref_pure')
    # if mapper=='star':
    #     command=[F"{metrics['programs']['path_to_star']} --runThreadN {metrics['CPU_count']} --runMode genomeGenerate  --genomeSAindexNbases {min(14,int(round(math.log(len(metrics['ref_seq']),2)/2-1)))} --genomeDir {path+'mapping_ref_pure'} --genomeFastaFiles {metrics['ref_file']}"]
    # elif mapper=='bwa':
    #     command=[F"cp {metrics['ref_file']} {path}mapping_ref_pure/mapping_ref_pure.fas",F"{metrics['programs']['path_to_bwa']} index {path}mapping_ref_pure/mapping_ref_pure.fas"]
    # if metrics['skip'][10]=='1':
    #     for cmd in command:
    #         proc=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #         while proc.poll() is None:#similar to wait until program terminates
    #             l=proc.stdout.readline().strip()
    #             if len(l)>0: print(l,flush=True)
    #         l=proc.stdout.read()
    #         if len(l)>0: print(l,flush=True)
    #
    # #remapping of merged reads
    # stage='final mapping'
    # print(stage)
    # if not os.path.exists(path+'06_2nd_mapping'): os.makedirs(path+'06_2nd_mapping')
    # for seq_typ in ['DNA','RNA']:
    #     for sampleID in track_files[seq_typ]:
    #         print(f"processing sample {sampleID}")
    #         temp_files=[]
    #         for dt in ['SE','PE']:
    #             new_file=f"{path}06_2nd_mapping/{sampleID}_{seq_typ}_{dt}_"
    #
    #             if mapper=='star':
    #                 if dt=='SE':
    #                     command=F"{metrics['programs']['path_to_star']} --runThreadN {metrics['CPU_count']} --genomeDir {path+'mapping_ref_pure'} --outSAMtype BAM Unsorted --readFilesIn {track_files[seq_typ][sampleID][0]} --readFilesCommand zcat --outSAMattributes All --outFilterType Normal --outMultimapperOrder Random --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFileNamePrefix {new_file}"
    #                 else:
    #                     command=F"{metrics['programs']['path_to_star']} --runThreadN {metrics['CPU_count']} --genomeDir {path+'mapping_ref_pure'} --outSAMtype BAM Unsorted --readFilesIn {track_files[seq_typ][sampleID][1]} {track_files[seq_typ][sampleID][2]} --readFilesCommand zcat --outSAMattributes All --outFilterType Normal --outMultimapperOrder Random --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFileNamePrefix {new_file}"
    #             elif mapper=='bwa':
    #                 if dt=='SE':
    #                     command=F"{metrics['programs']['path_to_bwa']} mem -t {metrics['CPU_count']} -R \"@RG\\tID:{sampleID}\\tLB:lib1\\tPL:illumina\\tSM:None\" {path}mapping_ref_pure/mapping_ref_pure.fas {track_files[seq_typ][sampleID][0]} | {metrics['programs']['path_to_samtools']} view -b -o {path}06_2nd_mapping/{sampleID}_{seq_typ}_{dt}_Aligned.out.bam --threads {metrics['CPU_count']} -"
    #                 else:
    #                     command=F"{metrics['programs']['path_to_bwa']} mem -t {metrics['CPU_count']} -R \"@RG\\tID:{sampleID}\\tLB:lib1\\tPL:illumina\\tSM:None\" {path}mapping_ref_pure/mapping_ref_pure.fas {track_files[seq_typ][sampleID][1]} {track_files[seq_typ][sampleID][2]} | {metrics['programs']['path_to_samtools']} view -b -o {path}06_2nd_mapping/{sampleID}_{seq_typ}_{dt}_Aligned.out.bam --threads {metrics['CPU_count']} -"
    #
    #             if metrics['skip'][11]=='1':
    #                 proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #                 while proc.poll() is None:#similar to wait until program terminates
    #                     l=proc.stdout.readline().strip()
    #                     if len(l)>0: print(l,flush=True)
    #                 l=proc.stdout.read()
    #                 if len(l)>0: print(l,flush=True)
    #             temp_files.append(f"{path}06_2nd_mapping/{sampleID}_{seq_typ}_{dt}_Aligned.out.bam")
    #
    #         #merging bams
    #         new_file=f"{path}06_2nd_mapping/{sampleID}_{seq_typ}_Aligned.out.bam"
    #         command=F"{metrics['programs']['path_to_samtools']} merge -f {new_file} {' '.join(temp_files)}"
    #         if metrics['skip'][11]=='1':
    #             proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #             while proc.poll() is None:#similar to wait until program terminates
    #                 l=proc.stdout.readline().strip()
    #                 if len(l)>0: print(l,flush=True)
    #             print(proc.stdout.read(),flush=True)
    #         track_files[seq_typ][sampleID]=new_file
    #
    # if options.attrition:
    #     print(f'collecting read counts for stage {stage}')
    #     attrition=get_read_counts(track_files,attrition,metrics['CPU_count'],stage)
    #     if options.attrition:
    #         create_attrition_plot(attrition,path+'%sattrition.svg'%metrics['prefix'])
    #
    # #final sorting and indexing
    # print('Final sorting of merged paired end reads, followed by indexing')
    # if not os.path.exists(path+'07_final_sort'): os.makedirs(path+'07_final_sort')
    #
    # skip_existing_files=[e[1] for e in sorted([(os.path.getmtime(f"{path}07_final_sort/{d}"),d) for d in os.listdir(f"{path}07_final_sort/") if d.endswith('.bam')])[:-1]]
    #
    # for seq_typ in ['DNA','RNA']:
    #     for sampleID in track_files[seq_typ]:
    #         print(f"processing sample {sampleID}")
    #         new_file=f"{path}07_final_sort/{sampleID}_{seq_typ}_sort.bam"
    #         #sorting
    #         command=F"{metrics['programs']['path_to_samtools']} sort -@ {metrics['CPU_count']} -O bam -T {path+'07_final_sort/tmp'} -o {new_file} {track_files[seq_typ][sampleID]} >> {metrics['log']} && {metrics['programs']['path_to_samtools']} index {new_file} >>{metrics['log']}"
    #         if metrics['skip'][12]=='1':
    #             if (options.auto_skip and not new_file.split('/')[-1] in skip_existing_files) or (not options.auto_skip):
    #                 proc=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #                 while proc.poll() is None:#similar to wait until program terminates
    #                     l=proc.stdout.readline().strip()
    #                     if len(l)>0: print(l,flush=True)
    #                 l=proc.stdout.read()
    #                 if len(l)>0: print(l,flush=True)
    #         track_files[seq_typ][sampleID]=new_file

    # if options.attrition:
    #     create_attrition_plot(attrition,path+'%sattrition.svg'%metrics['prefix'])

    return(True)

###main
print('parsing program parameters')
metrics=parse_parameters(options,metrics)


print('Initiating analysis')
orchestrate_analysis(metrics)


print('\nAnalysis finished successfully! Have a great day.\n')
