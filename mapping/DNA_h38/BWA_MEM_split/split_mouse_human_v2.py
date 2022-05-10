#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 14:10:41 2017

@author: Hilmar Berger

Hilmar Berger 
(C) Max Planck Institute for Infection Biology, 2018
"""

import os
import pysam
import string
import copy

def species(contig_name):
    return ("mm" if contig_name[0:3]=="mm_" else "hs")
        
def total_mapping_len(read):
    return(sum([e[1] for e in read.cigartuples if e[0]=="M"]))

def split_reads(input_file, human_output, mouse_output, cross_output, min_length, min_temp_length):
    """
    split human and mouse reads
    """
    total_read_cnt = 0
    excluded_read_cnt = 0
    hs_cnt = 0
    mm_cnt = 0
    
    ifile = pysam.AlignmentFile(input_file,"rb")

    ref_dict = {}
    for n in xrange(0,ifile.nreferences):
        ref_dict[n] = ifile.getrname(n)

    #print ref_dict

    header_orig = ifile.header
    header_hs = copy.deepcopy(header_orig)
    header_mm = copy.deepcopy(header_orig)
    
    new_sq = []
    for contig in header_hs["SQ"]:
        if contig["SN"][0:3]=="mm_":
            continue
        new_sq.append(contig)
    header_hs["SQ"] = new_sq

    new_sq = []
    for contig in header_mm["SQ"]:
        cname = contig["SN"]
        if cname[0:3]!="mm_":
            continue
        contig["SN"] = cname.replace("mm_","")
        new_sq.append(contig)        
    header_mm["SQ"] = new_sq
    
    hs_file = pysam.AlignmentFile(human_output, "wb", header=header_hs)
    mm_file = pysam.AlignmentFile(mouse_output, "wb", header=header_mm)
    cross_file = pysam.AlignmentFile(cross_output, "wb", header=header_orig)

    ref_dict_hs = {}
    ref_dict_hs_rev = {}
    for n in xrange(0,hs_file.nreferences):
        ref_dict_hs[n] = hs_file.getrname(n)
        ref_dict_hs_rev[hs_file.getrname(n)] = n

    ref_dict_mm = {}
    ref_dict_mm_rev = {}
    for n in xrange(0,mm_file.nreferences):
        ref_dict_mm[n] = mm_file.getrname(n)
        ref_dict_mm_rev[mm_file.getrname(n)] = n
    
#    print mm_file.header
#
#    print "====================== HUMAN ====================="
#    for e in header_hs["SQ"]:
#        print "%s\t%d" % (str(e),ref_dict_hs_rev[e["SN"]])
#    print "====================== MOUSE ====================="
#    for e in header_mm["SQ"]:
#        print "%s\t%d" % (str(e),ref_dict_mm_rev[e["SN"]])
    
    
    for rr in ifile.fetch(until_eof=True):
        total_read_cnt += 1
        # both reads unmapped
        if rr.is_unmapped and rr.mate_is_unmapped:
            cross_file.write(rr)
            excluded_read_cnt += 1
            continue
        
        # Treattment of mapped reads: 
        # - both ends mapped: 
        #    - both mm -> mm
        #    - both hs -> hs
        #    - mixed hs/mm -> cross mapping
        # - single end mapped: 
        #    - hs or mm, according to to mapped contig
        # Additional, unimplemented idea: use XA tag if present: count mismatches between XA contigs and follow rule above
        
        # one might expect that for unmapped mates the next contig is set to "*", however, that is not the case
        # instead, next contig is set to the same contig as the mapped read of the pair (first=last). 
        # We therefore trust that all next contigs are pointing to a valid contig
        # one might expect that for unmapped mates the next contig is set to "*", however, that is not the case
        contig_orig = ref_dict[rr.reference_id]
        next_contig_orig = ref_dict[rr.next_reference_id]
        
        if not rr.mate_is_unmapped and not rr.is_unmapped: # both mates mapped
            if species(contig_orig) == species(next_contig_orig):
                if species(contig_orig) == "hs":
                    rr.reference_id = ref_dict_hs_rev[contig_orig]
                    rr.next_reference_id = ref_dict_hs_rev[next_contig_orig]
                    hs_file.write(rr)
                    hs_cnt += 1
                else:
                    rr.reference_id = ref_dict_mm_rev[contig_orig.replace("mm_","")]
                    rr.next_reference_id = ref_dict_mm_rev[next_contig_orig.replace("mm_","")]
                    mm_file.write(rr)
                    mm_cnt += 1
            else:
                cross_file.write(rr)
                excluded_read_cnt += 1
        elif not rr.is_unmapped: # this read mapped but mate unmapped
            if species(contig_orig)=="hs":
                rr.reference_id = ref_dict_hs_rev[contig_orig]
                rr.next_reference_id = ref_dict_hs_rev[next_contig_orig]
                hs_file.write(rr)
                hs_cnt += 1
            else:
                rr.reference_id = ref_dict_mm_rev[contig_orig.replace("mm_","")]
                rr.next_reference_id = ref_dict_mm_rev[next_contig_orig.replace("mm_","")]
                mm_file.write(rr)
                mm_cnt += 1
        else: # this read unmapped but mate mapped
            if species(next_contig_orig)=="hs":
                rr.reference_id = ref_dict_hs_rev[contig_orig]
                rr.next_reference_id = ref_dict_hs_rev[next_contig_orig]
                hs_file.write(rr)
                hs_cnt += 1
            else:
                rr.reference_id = ref_dict_mm_rev[contig_orig.replace("mm_","")]
                rr.next_reference_id = ref_dict_mm_rev[next_contig_orig.replace("mm_","")]
                mm_file.write(rr)
                mm_cnt += 1
            
        
        #tags_orig = rr.tags
        
    hs_file.close()
    mm_file.close()
    cross_file.close()

    print "Total reads: %d\tHuman: %d (Prop. %g)\tMouse: %d\tExcluded: %d" % (total_read_cnt, hs_cnt, float(hs_cnt)/total_read_cnt, mm_cnt, excluded_read_cnt)

####################### Main ################################
# This is the interface to the command line
#############################################################    
if __name__ == '__main__':
    import argparse
    from operator import itemgetter
    
    description_text = """Split reads mapped from human and mouse genomes into different files.
    Assumes that all mouse contigs are marked with mm_ and all other contigs are not. 
    Fragments with both mates in either Hs or Mm are stored in corresponding files. 
    Fragments with only a single mapped mate will stored according to the mapped mate. 
    Fragments with one mate mapping to Hs and the other to Mm as well as all unmapped reads will be stored in the BAM for excluded reads. 
    """
    
    parser = argparse.ArgumentParser(description=description_text, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_file', type=str, help='Input BAM file')

    parser.add_argument('-o','--output-folder', action='store', type=str, default=None,
                       help='Output file folder. Default: same as input files.')

    args = parser.parse_args()

    ###########################################################################
    # Main program                                                            #
    ###########################################################################
    inp_file = args.input_file
    
    bn = os.path.splitext(os.path.basename(inp_file))[0]
    if args.output_folder is None:
        dn = os.path.dirname(inp_file)
    else:
        dn = args.output_folder
    
    hs_filename = os.path.join(dn, bn + "_hs.bam")
    mm_filename = os.path.join(dn, bn + "_mm.bam")
    excluded_filename = os.path.join(dn, bn + "_excluded.bam")
    
    split_reads(inp_file, hs_filename, mm_filename, excluded_filename, 1e6, 1e6)
    
#runfile('/data_genome2/projects/CINOCA/Cervix_Clones/src/split_mouse_human.py',args='"/data_genome2/projects/CINOCA/Cervix_Clones/mapping/DNA/BWA_MEM_with_mm38/test.bam"')       