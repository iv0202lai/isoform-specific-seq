# Find specific region for a transcript

""" In all of the following, the list of intervals must be sorted and non-overlapping. 
We also assume that x is in tp(start, end) iff start <= x and x < end."""

import os
import time
import vcf
from cogent.db.ensembl import Genome

# Functions for interval calculation
class tp():
   def __repr__(self):
       return '(%d,%d)' % (self.start, self.end)
   def __init__(self,start,end): 
       self.start=start
       self.end=end

def rev(n):
    if n == 'A': n='T'
    elif n == 'T': n='A'
    elif n == 'C': n='G'
    elif n == 'G': n='C'
    return n

def flatten(list_of_tps):
  #Convert a list of intervals to a list of endpoints
  return reduce(lambda ls, ival: ls + [ival.start, ival.end],
                list_of_tps,
                [])

def unflatten(list_of_endpoints):
   #Convert a list of endpoints, with an optional terminating sentinel,
   #into a list of intervals
   return [tp(list_of_endpoints[i], list_of_endpoints[i + 1])
          for i in range(0, len(list_of_endpoints) - 1, 2)]

def merge(a_tps, b_tps, op):
  #Merge two lists of intervals according to the boolean function operator
  a_endpoints = flatten(a_tps)
  b_endpoints = flatten(b_tps)

  sentinel = max(a_endpoints[-1], b_endpoints[-1]) + 1
  a_endpoints += [sentinel]
  b_endpoints += [sentinel]

  a_index = 0
  b_index = 0

  res = []

  scan = min(a_endpoints[0], b_endpoints[0])
  while scan < sentinel:
    in_a = not ((scan < a_endpoints[a_index]) ^ (a_index % 2))
    in_b = not ((scan < b_endpoints[b_index]) ^ (b_index % 2))
    in_res = op(in_a, in_b)

    if in_res ^ (len(res) % 2): res += [scan]
    if scan == a_endpoints[a_index]: a_index += 1
    if scan == b_endpoints[b_index]: b_index += 1
    scan = min(a_endpoints[a_index], b_endpoints[b_index])

  return unflatten(res)

def interval_diff(a, b):
  return merge(a, b, lambda in_a, in_b: in_a and not in_b)

def interval_union(a, b):
  return merge(a, b, lambda in_a, in_b: in_a or in_b)

def interval_intersect(a, b):
  return merge(a, b, lambda in_a, in_b: in_a and in_b)

# Functions for finding sequence
def specific_region(transID):
    # Return gene object, transcript object, and isoform-specific regions of a transID
    # specific_regions contain tuples of (interval,sequence)
      
    start = time.time()
    print "\nGetting data of transcript",transID,"from Ensembl"
    target = mouse.getTranscriptByStableId(StableId=transID)
    
    if target is None:
        # If no transcript found regarding to the transID
        end = time.time()
        print "Transcript ID", transID, "not found!!!"
        print "\nExecution time:", end - start, "seconds."
        return None,None,[]
    else:
        print "...Getting gene data"
        gene = target.Gene
        
        specific_region = _specific_region(gene,target)
        end = time.time()
        
        # Automatically add SNP label if reference vcf exists
        if os.path.isfile('CAST_snps.vcf.gz') == False:
            print "CAST_snps.vcf.gz doesn't exist. Print out results without SNP labels."
        elif os.path.isfile('CAST_snps.vcf.gz.tbi') == False:
            print "CAST_snps.vcf.gz.tbi doesn't exist. Print out results without SNP labels."
        else:
            specific_region = label_snp(target,specific_region)
        
        print "\nExecution time:", end - start, "seconds."
        return gene,target,specific_region
        
def _specific_region(gene,target):
    
    #Get range union of other isoforms of the same gene
    print "...Listing transcripts of gene",gene.Symbol,"(" + gene.StableId + ")"
    transcripts = list(gene.Transcripts)
    
    if len(transcripts) == 1:
        print "..."+target.StableId, "has no other isoform in the gene."
        range_seqs = ranges_seq(gene,exon_union(target))
        return range_seqs
    else:
        print "..."+str(len(transcripts)).strip(),"transcripts found in",gene.Symbol,\
        "(" + gene.StableId + ")"
    
    #Get specific ranges in the target transcript
    transcripts.remove(target)
    diff = interval_diff(exon_union(target),transcript_union(transcripts))
    range_seqs = ranges_seq(gene,diff)
    
    if diff == []:
        print "\nNo specific region found LOL"
    else:
        if len(diff) == 1:
            print "\nAn isoform-specific region has been found!"
        else:
            print "\n"+str(len(diff)).strip(),"isoform-specific regions found!"
    
    return range_seqs

def exon_union(transcript):
    #Return ranges for a transcript object in order
    ranges = []
    print "...Contatenating exons in transcript",transcript.StableId
    for exon in transcript.Exons:
        location = exon.Location
        range = tp(int(location.EnsemblStart),int(location.EnsemblEnd)+1)
        ranges.append(range)
    if str(transcript.Location).split(':')[-1] == '-1':
        ranges.reverse() # Note that exons come in decreasing order if on reverse strand!
    return ranges
    
def transcript_union(transcripts):
    #Convert list of transcript objects into range union
    ranges = []
    for transcript in transcripts:
        if ranges == []:
            ranges = exon_union(transcript)
        else: 
            ranges = interval_union(ranges,exon_union(transcript))
    return ranges

def ranges_seq(gene,ranges):
    #Return tuples of (range,sequence) from (gene, ranges list)
    print "...Generating sequences from regions found"
    
    geneseq = gene.Seq
    range_seqs = []
    
    if str(gene.Location).split(':')[-1] == '-1':
        # for genes on reverse strand
        end = gene.Location.EnsemblEnd
        for range in ranges:
            interval = [end - range.end + 1, end - range.start + 1]
            seq = geneseq[interval[0]:interval[1]]
            range_seqs.append((range,seq))
    else:
        # for genes on forward strand
        start = gene.Location.EnsemblStart
        for range in ranges:
            interval = [range.start - start, range.end - start]
            seq = geneseq[interval[0]:interval[1]]
            range_seqs.append((range,seq))
    
    return range_seqs

def label_snp(transcript,range_seqs):
    range_seq_snp = []
    chrom = str(transcript.Location).split(":")[2]
    strand = str(transcript.Location).split(':')[-1]
    
    if strand == '-1':
        # for genes on reverse strand
        for index,(range,seq) in enumerate(range_seqs):
            print "\nCreating SNP data for Region",index+1,\
            "("+str(range.start)+"-"+str(range.end-1)+")"
            seq = _label_snp_rev(chrom,range,seq)
            range_seq_snp.append((range,seq))
    else:
        # for genes on forward strand
        for index,(range,seq) in enumerate(range_seqs):
            print "\nCreating SNP data for Region",index+1,\
            "("+str(range.start)+"-"+str(range.end-1)+")"
            seq = _label_snp_for(chrom,range,seq)
            range_seq_snp.append((range,seq))
    return range_seq_snp
        
def _label_snp_rev(chrom,range,seq):
    vcf_reader = vcf.Reader(filename='CAST_snps.vcf.gz')
    snps = []
    
    for snp in vcf_reader.fetch(chrom, range.start, range.end):
        if snp.samples[0]['GT'] != '1/1': continue
        if len(snp.ALT) > 1: continue
        if seq[range.end - snp.POS - 1] != rev(snp.REF):
            print "Data not coincident:",snp
            continue
        print "...",snp
        label = '[' + rev(snp.REF) + '/' + rev(str(snp.ALT[0])) + ']'
        snps.append((snp.POS,label))
    
    if snps == []:
       print "... No SNPs found Orz"
    else:       
       for pos,label in snps:
           pos = range.end - pos - 1
           seq = str(seq[:pos]) + label + str(seq[pos+1:])
        
    return seq

def _label_snp_for(chrom,range,seq):
    vcf_reader = vcf.Reader(filename='CAST_snps.vcf.gz')
    snps = []
    
    for snp in vcf_reader.fetch(chrom, range.start, range.end):
        if snp.samples[0]['GT'] != '1/1': continue
        if len(snp.ALT) > 1: continue
        if seq[snp.POS - range.start] != snp.REF:
            print "Data not coincident:",snp
            continue
        print "...",snp
        label = '[' + snp.REF + '/' + str(snp.ALT[0]) + ']'
        snps.append((snp.POS,label))
        snps.reverse()
     
    if snps == []:
        print "... No SNPs found Orz"
    else:          
        for pos,label in snps:
            pos = pos - range.start
            seq = str(seq[:pos]) + label + str(seq[pos+1:])
        
    return seq

def print_seq(range_seqs):
    for index,(range,seq) in enumerate(range_seqs):
        print "\nRegion",index+1,
        print "   Position:",str(range.start)+"-"+str(range.end-1),
        print "   Length:", range.end - range.start,
        snp_count = seq.count("/")
        if snp_count > 0:
            print "   SNP count:", snp_count
        else:
            print "   no SNP"
        print str(seq)

def seq_output(gene,transcript,range_seqs):
    # Export result file in ./isoform_specific_seq
    dir = os.path.join(os.path.dirname(__file__), 'isoform_specific_seq')
    if not os.path.exists(dir):
        os.mkdir(dir)
    filename = gene.Symbol+'_'+transcript.StableId+'.txt'
    output = open(os.path.join(dir,filename), 'w')
    
    # Write in basic data of the transcript
    output.write("Transcript ID: "+transcript.StableId)
    output.write("\nGene: "+gene.Symbol+'('+gene.StableId+')')
    output.write("\nChromosome: "+str(transcript.Location).split(':')[2])
    strand = str(transcript.Location).split(':')[-1]
    if strand == '1':
        output.write("       Strand: Forward")
    else:
        output.write("       Strand: Reverse")
    
    # Write in regions found 
    if range_seqs == []:
        output.write("\nNo isoform-specific region found.")
    for index,(range,seq) in enumerate(range_seqs):
        header = "\nRegion " + str(index+1) \
        + "   Position: " + str(range.start) + "-" + str(range.end-1)\
        + "   Length: " + str(range.end - range.start)
        snp_count = seq.count("/")
        if snp_count > 0:
            header = header + "   SNP count: " + str(snp_count)
        else:
            header = header + "   no SNP" 
        output.write("\n"+header)
        output.write("\n"+str(seq))
        output.close

#Setting ensembl parameters 
release = 81
species = 'Mus musculus'
mouse = Genome(Species=species,Release=release)

input = raw_input('\nEnter ensembl mouse transcript ID or txt file: ')

if input.strip().split('.')[-1] == 'txt':
    try:
        fh = open(input.strip(),'r')
        for transID in fh:
            transID = transID.strip().split('.')[0]
            print "\n"+"============ Start to manipulate",transID,"============"
            gene,transcript,range_seqs = specific_region(transID)
            if transcript is None: continue
            seq_output(gene,transcript,range_seqs)
        fh.close()
    except IOError:
        print "File", input, "not found!!"
    
else:
    while True:
        to_output = raw_input('Export result to an txt file? (Y/N)  ').lower()
        if to_output == 'y' or to_output == 'n': break
    transID = input.split('.')[0]
    gene,transcript,range_seqs = specific_region(transID)
    print_seq(range_seqs)
    if to_output == 'y' and transcript is not None:
        seq_output(gene,transcript,range_seqs)
    
        
        