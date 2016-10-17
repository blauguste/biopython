from Bio import SeqIO, SeqFeature
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
import sys

def make_IGR_feature(start_locus, end_locus):
    igr_feature_location = FeatureLocation(start_locus, end_locus, strand=None)
    igr_feature_type = "misc_RNA"
    igr_feature = SeqFeature(igr_feature_location, type=igr_feature_type, strand=None)
    return igr_feature

#The genbank file in which gene regions are defined and 
#to which new IGRs will be appended
def get_ds_igrs(genbank_path, outpath, min_igr_length=50):

    base_record = SeqIO.read(open(genbank_path), "genbank")
    gene_list = []
    ds_igrs = []
    igr_count = 0

    #Loop over the genome file, add each of the genes to a list
    for feature in base_record.features:
        if feature.type == "gene":
            left_pos = feature.location.start
            right_pos = feature.location.end
            gene_name = feature.qualifiers["locus_tag"][0]
            ref_strand = feature.strand
            gene_list.append((gene_name, left_pos, right_pos, ref_strand))

    #Sort list of genes by start position
    gene_list.sort(key=lambda tup: tup[1])

    base_gene_list_position = 0
    comparison_gene_list_position = 1

    #Make the IGR between "beginning" of genome and the gene with the smallest start position
    if gene_list[0][1] > 0:
        first_gene_left = gene_list[0][1]
        first_gene_name = gene_list[0][0]
        intergene_seq = base_record.seq[:first_gene_left]
        igr_feature = make_IGR_feature(0, first_gene_left)
        igr_feature.qualifiers["locus_tag"] = "start_IGR_" + first_gene_name
        base_record.features.append(igr_feature)
        igr_count = 1

    #Compare the right boundary of the "first" gene with the left boundary of the subsequent gene
    for i, gene_record in enumerate(gene_list[comparison_gene_list_position:]):
        base_gene_right = gene_list[base_gene_list_position][2]
        comparison_gene_left = gene_list[comparison_gene_list_position][1]
        comparison_gene_right = gene_list[comparison_gene_list_position][2]
        base_gene_name = gene_list[base_gene_list_position][0]
        comparison_gene_name = gene_list[comparison_gene_list_position][0]
        #If there's a gap of desired size between the two genes...
        if comparison_gene_left - base_gene_right >= min_igr_length:
            #Create new features and add them to the genbank file
            intergene_seq = base_record.seq[base_gene_right:comparison_gene_left]
            igr_feature = make_IGR_feature(base_gene_right, comparison_gene_left)
            igr_feature.qualifiers["locus_tag"] = base_gene_name + "_IGR_" + comparison_gene_name
            base_record.features.append(igr_feature)
            igr_count += 1
            #Reset the base gene and comparison position (move to next gene)
            base_gene_list_position = comparison_gene_list_position
            comparison_gene_list_position = base_gene_list_position + 1
        #If the right boundary of the comparison gene is greater than that
        #of the base gene, make the comparison gene the new base gene
        elif comparison_gene_right > base_gene_right:
            base_gene_list_position = comparison_gene_list_position
            comparison_gene_list_position = base_gene_list_position + 1
        #If none of these conditions are met, repeat with the next comparison gene
        else:
            comparison_gene_list_position += 1

    #Sort the gene list by end position
    gene_list.sort(key=lambda tup: tup[2])

    #Make the IGR between the "last" gene and the "end" of the genome
    if gene_list[-1][1] < len(base_record.seq):
        last_gene_right = gene_list[-1][2]
        last_gene_name = gene_list[-1][0]
        intergene_seq = base_record.seq[last_gene_right:]
        igr_feature = make_IGR_feature(last_gene_right, len(base_record.seq))
        igr_feature.qualifiers["locus_tag"] = last_gene_name + "_IGR_end"
        base_record.features.append(igr_feature)
        igr_count += 1

    #Write new IGR features to the genbank
    SeqIO.write(base_record, open(outpath,"w"), "genbank")

    print(str(igr_count) + " intergenic regions identified and appended to genbank file.")

if __name__ == '__main__':
    if len(sys.argv) == 3:
         get_ds_igrs(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
         get_ds_igrs(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    else:
         print("Usage: get_ds_igrs.py gb_file_in gb_file_out [intergenic_length]")
         sys.exit(0)