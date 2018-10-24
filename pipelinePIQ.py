import pandas as pd
import CGATCore.IOTools as IOTools
import numpy as np
from CGATCore.Pipeline import cluster_runnable
import matplotlib as mpl
mpl.use('Agg') # So that matplotlib.pyplot doesn't try to find the X-server and give an error
import matplotlib.pyplot as plt
import networkx as nx
import tempfile
import os


# Gets a calls file with row numbers (calls_processed_file) and extracts the
# significant calls (by row number) using sign_calls_ids_processed_file
# into output_file
# 
# Inputs:
#     -calls_processed_file: File containing all calls for the binding sites 
#        (whether significant or not). First row is the row number. Format:
#        1       chr1    42640   42655   PB01051Arid3a2  583     +
#        2       chr1    73772   73787   PB01051Arid3a2  587     +
#        3       chr1    144505  144520  PB01051Arid3a2  587     +
#
#     -sign_calls_ids_processed_file: File containing the row numbers of the significant matches
#        and the purity.
#        Format: 
#        10      0.744213649851632
#        22      0.742762221167537
#        27      0.735926305015353
#
#     -output_file: Only the significant binding sites without the row numbers (substitutes purity for score)
#        NOTE: Score is purity x 1000 but has less precision
#        chr1    1380372 1380389 PB01061Arid5a2  0.713     +
#        chr1    2062819 2062836 PB01061Arid5a2.RC       0.738     -
#        chr1    2392460 2392477 PB01061Arid5a2.RC       0.739     -
def filterSignCalls(calls_processed_file,
                    sign_calls_ids_processed_file,
                    output_file):

    bed_table = pd.read_table(calls_processed_file,
                              compression='gzip',
                              header=None,
                              delimiter="\t")
    
    id_table = pd.read_table(sign_calls_ids_processed_file,
                             compression='gzip',
                             header=None,
                             delimiter="\t")
 
    # For each table specify the column name fo the columns
    bed_table_columns = list(bed_table.columns)
     
    bed_table_columns[0] = "hit_id"
    bed_table_columns[1] = "chr"
    bed_table_columns[2] = "start"
    bed_table_columns[3] = "end"
    bed_table_columns[4] = "TF"
    bed_table_columns[5] = "score"
    bed_table_columns[6] = "strand"
    
     
    bed_table.columns = bed_table_columns
     
     
    id_table_columns = list(id_table.columns)
     
    id_table_columns[0] = "hit_id"
    id_table_columns[1] = "purity"
     
    id_table.columns = id_table_columns
     
    result = pd.merge(id_table, bed_table, how="left", on=["hit_id"])
    
    # Drop the hit_id and score columns
    result = result.drop(['hit_id','score'], 1) 
    
    # Reorder
    result = result[['chr', 'start', 'end', 'TF', 'purity', 'strand']]
     
    # Output the table
    result.to_csv(output_file,
                 compression="gzip",
                 sep="\t",
                 header=None,
                 mode="w",
                 index_label=None,
                 index=False,
                 line_terminator="\n")
              
    

# Based on PipelineChromHMM.getUniquePatientsIDsFromCellMarkFileTable
# Gets unique peaks file from the bam_to_peaks_table
#
# Inputs:
#     -bam_to_peaks_table: Generated before with structure:
#        CCDN1_merged.bam        /mnt/fastdata/mbp15ja/atac_seq_balanced_consensus_peaks_MM_subgrouping_all_quality_MM_ND_samples_by_14_07_2018_CCND1/common_peaks.dir/all_merged_narrow_broad_peaks.gz
#        HD-merged_new.bam       /mnt/fastdata/mbp15ja/atac_seq_balanced_consensus_peaks_MM_subgrouping_all_quality_MM_ND_samples_by_14_07_2018_HD/common_peaks.dir/all_merged_narrow_broad_peaks.gz
#
#     -bam_file: bam file to look for to return the corresponding peaks file
#
# Outputs: Corresponding peaks file
#
# Exception:
#     -If no corresponding bam file is found
def getPeaksFileFromBam(bam_to_peaks_table, bam_file):

    peaks_file = ""

    with IOTools.open_file(bam_to_peaks_table, "r") as reader:

        for line in reader:

            line = line.rstrip('\n')

            bam = line.split("\t")[0]

            if(bam == bam_file):
               peaks_file = line.split("\t")[1]
               break

    reader.close()

    if(peaks_file == ""):
        raise Exception("No corresponding bam file found")

    return peaks_file





# For each gene, merge the regions (gene body and 2500bp upstream of the TSS).
# Any space between these two regions for each gene is considered to be part of the region too.
# Both input lists are checked to have the same genes.
# Inputs:
#     -gene_body_table: A table with the gene bodies of a set of genes. Format:
#        chr1    11868   14409   ENSG00000223972 0       +
#        chr1    14403   29570   ENSG00000227232 0       -
#        chr1    17368   17436   ENSG00000278267 0       -
#
#     -tss_upstream_table: A table with the gene TSS of a set of genes. Format:
#        chr1    11868   14409   ENSG00000223972 0       +
#        chr1    14403   29570   ENSG00000227232 0       -
#        chr1    17368   17436   ENSG00000278267 0       -
#
#     -output_file: The output table. Format:
#        chr1    11868   14409   ENSG00000223972 0       +
#        chr1    14403   29570   ENSG00000227232 0       -
#        chr1    17368   17436   ENSG00000278267 0       -
#
# Exception:
#     -If the gene_body_table and tss_upstream_table contain different genes.
def merge_TSS_gene_body_per_gene(gene_body_table,
                                 tss_upstream_table,
                                 output_file):

    gene_body_df = pd.read_table(gene_body_table,
                             header=None,
                             delimiter="\t")


    tss_upstream_df = pd.read_table(tss_upstream_table,
                             header=None,
                             delimiter="\t")


    gene_body_df.columns = ["chr_gene_body", "start_gene_body", "end_gene_body", "gene", "score_gene_body", "strand_gene_body"]

    tss_upstream_df.columns = ["chr_tss_upstream", "start_tss_upstream", "end_tss_upstream", "gene", "score_tss_upstream", "strand_tss_upstream"]

    # Merge on the gene id
    # We are going to check if there are any genes not occuring in both
    result_table = gene_body_df.merge(tss_upstream_df, on="gene", how="outer")

    # Report if some gene is not in the list
    if (result_table.isnull().values.any()):
        raise Exception("There are genes which are not in both lists")

    # Now for each row:
    # if the strand is +, take the start of the tss_upstream as the new start and end of the gene_body as the new end
    # if the strand is -, take the end of the tss_upstream as the new end and the start of the gene body as the new start
    for index, row in result_table.iterrows():
        if row["strand_gene_body"] == "+":
            result_table.at[index, "new_start"] = row["start_tss_upstream"]
            result_table.at[index, "new_end"] = row["end_gene_body"]

        elif row["strand_gene_body"] == "-":
            result_table.at[index, "new_start"] = row["start_gene_body"]
            result_table.at[index, "new_end"] = row["end_tss_upstream"]

    # Subselect only the rows with the chr, new coordinates, gene, score and strand
    # Output only the coordinate columns
    output_df = result_table.loc[:,["chr_tss_upstream", "new_start", "new_end","gene","score_gene_body","strand_gene_body"]]

    # Output
    output_df.to_csv(output_file,
                 float_format="%.f", # Integers may have been converted to floats
                 sep="\t",
                 header=False,
                 compression="gzip",
                 mode="w",
                 index_label=None,
                 index=False,
                 line_terminator="\n")


# Starts with a file containing location of the binding site, Transcription factor, purity score, nearest gene (TSS)
# and distance, calculates the interaction score (as per Rendeiro et. al 2016) for each TF (file) - gene and outputs a table
# with the interaction score between the TF represented by the file name and each gene.
#
# Inputs:
#   -tf_gene_purity_distance: location of the binding site, Transcription factor, purity score, nearest gene (TSS)
#       and distance. Format:
#       chrX    100636257       100636262       MA03661RGM1     0.8057527592550009      +       ENSG00000000003 3730
#       chr20   50948128        50948133        MA03661RGM1     0.7328857020819399      +       ENSG00000000419 10423
#       chr20   50948382        50948387        MA03661RGM1     0.7972367381429591      +       ENSG00000000419 10169
#
#   -output_file: interaction score between the TF represented by the file name and each gene
#       gene    interaction_score
#       ENSG00000000971 0.5421240520604342
#       ENSG00000002726 0.15164771529245755
@cluster_runnable
def calculateInteractionScore(tf_gene_purity_distance, output_file):

    tf_gene_purity_distance_df = pd.read_table(tf_gene_purity_distance,
                              compression='gzip',
                              header=None,
                              delimiter="\t")

    # Rename columns
    tf_gene_purity_distance_df.columns = ["chr", "start", "end", "tf", "purity", "strand", "gene", "distance"]

    # Get unique genes
    unique_genes = tf_gene_purity_distance_df.gene.unique()

    # Get the number of unique genes
    num_uniq_genes = len(unique_genes)

    # Create empty dataframe with the number of unique genes as rows
    output_cols = ["gene", "interaction_score"]

    output_df = pd.DataFrame(np.nan, index=range(0,num_uniq_genes), columns=output_cols)

    # To know which row to update
    output_row_counter = 0

    # Create an interaction score column empty
    tf_gene_purity_distance_df["interaction_score"] = np.nan

    # Go through each unique gene
    for unique_gene in unique_genes:

        # Get the interaction score for each gene
        tf_gene_purity_distance_df.loc[tf_gene_purity_distance_df['gene'] == unique_gene, "interaction_score"] = 2 * (tf_gene_purity_distance_df["purity"] - 0.5) * (10 ** (-tf_gene_purity_distance_df["distance"] / 100000))

        # Get the sum of the scores
        sum_interaction_score = tf_gene_purity_distance_df.loc[tf_gene_purity_distance_df['gene'] == unique_gene, "interaction_score"].sum()

        # Update the gene
        output_df.iloc[output_row_counter, 0] = unique_gene

        # Update the sum
        output_df.iloc[output_row_counter, 1] = sum_interaction_score

        # Update the position counter
        output_row_counter = output_row_counter + 1

    # Output the file with the header
    output_df.to_csv(output_file,
                         sep="\t",
                         header=True,
                         compression="gzip",
                         mode="w",
                         index_label=None,
                         index=False,
                         line_terminator="\n")


# Copied from pipelineAtacseq.histogramFromOneDValues
#
# Requirements: The values are equal or higher to 0
# Computes a histogram taking the values from values_file and
# using the specified number_bins. Stores the output on the outfile.
#
# Inputs:
#     -values_file: File containing one value per line.
#     -number_bins: Number of bins to produce for the histogram.
#     -outfile: Outfile containing the histogram
#
# Outputs:
#     -A histogram with the characteristics specified
@cluster_runnable
def histogramFromOneDValues(values_file, number_bins, outfile):
    # We are going to do two passes:
    # 1) Get the maximum score of all the series (for the bins).
    # 2) Create the numpy array of instances to populate the histogram

    # We are going to read the values_file in chunks at a time, this indicates the size of each chunk
    chunk_size = 10000000

    # Pass 1
    # Go through the values_file, get:
    # -the maximum score
    maximum_score = float(0)

    # Read a chunk of lines in the values_file
    for chunk in pd.read_table(values_file,
                               names=['ends_distance'],
                               header=None,
                               chunksize=chunk_size,
                               delimiter='\t'):

        # Update maximum distance if it has been surpassed
        max_score_chunk = chunk['ends_distance'].max()

        if (max_score_chunk > maximum_score):
            maximum_score = max_score_chunk

    # Now we are ready for pass 2

    # Get the bin edges for the histogram based on the max
    # min is always 0.0
    bin_edges = np.linspace(0.0, maximum_score, number_bins + 1)

    # These are going to be the actual bins, where in each
    # position we put the count
    # np.ndarray filled with np.uint32 zeros, CHANGED TO int64
    total_hist_data = np.zeros(number_bins, np.int64)

    # Read a chunk of lines in the values_file
    for chunk in pd.read_table(values_file,
                               names=['ends_distance'],
                               header=None,
                               chunksize=chunk_size,
                               delimiter='\t'):

        # Iterates rows
        for row in chunk.itertuples():
            # Each row will contain score and the counts for that score
            # Bin the data for one score
            subtotal, edges = np.histogram(row[1], bins=bin_edges)

            # accumulate bin counts over chunks
            total_hist_data += subtotal

    # Turn interactive plotting off, we don't want to show it
    plt.ioff()

    # Plot the histogram
    plt.hist(bin_edges[:-1], bins=bin_edges, weights=total_hist_data)

    # Save the file, remove whitespace around the image
    plt.savefig(outfile)

    # Clear the plot for the next one
    plt.clf()




# Gets the header from infile, adds "tf" and separates by tab in the first column and outputs it to outfile
#
# Inputs:
#   -infile: The infile with format:
#       gene    interaction_score
#       ENSG00000002726 0.8585400707523888
#       ENSG00000002746 0.5130992616557325
#       ENSG00000003402 0.15516727449470588
#
#   -outfile: The outfile where the new header is outputted
#       tf  gene    interaction_score
def addTFToHeader(infile, outfile):

    first_line = ""

    with IOTools.open_file(infile, "r") as reader:

            for line in reader:

                first_line = line

                break


    reader.close()

    first_line = "tf\t" + first_line

    with IOTools.open_file(outfile, "w") as writer:

        writer.write(first_line)

    writer.close()


# Copied from AuxiliaryPrograms/Table_operations/table_operations
# Gets a list of all the fields in the first line after applying the separator
# Inputs:
#     -infile: The file to parse
#     -separator: The separator of the fields.
#
# Outputs:
#     -list of fields.
def getFieldsInFirstLine(infile, separator="\t"):
    # Get the number of fields and the header and the text of the header

    header_fields = []

    # First determine the number of fields
    with IOTools.open_file(infile, "r") as reader:
        for line in reader:
            # Remove the new line character from the line (avoids taking this into account downstream)
            line = line.rstrip('\n')

            header_fields = line.split(separator)

            break

    reader.close()

    return header_fields



# Copied from AuxiliaryPrograms/Genomic_regions/genomic_regions
# Starts with "infile", reads the infile by chunks of predefined size and extracts the top highest "max_num_edges" number
# of rows from the column "column_name" in the full file. Outputs the results into outfile
#
# Inputs:
#   -infile: A file WITH HEADER containing the "column_name" column. For example, format:
#       gene    tf      interaction_score       gene_desc
#       CFH     AGL3    0.476420431210367       complement factor H
#       STPG1   AGL3    0.319511602570206       sperm tail PG-rich repeat containing 1
#
#   -outfile: The resulting outfile WITH HEADER containing the "column_name" column. For example, format:
#       gene    tf      interaction_score       gene_desc
#       CFH     AGL3    0.476420431210367       complement factor H
#       STPG1   AGL3    0.319511602570206       sperm tail PG-rich repeat containing 1
#
#   -max_num_edges: Top number of rows to extract from infile (largest values in "column_name")
#
#   -column_name: The column name in infile to use as filtering criteria.
#
# Outputs:
#   -The filtered outfile
#
# Exception:
#   -If the column is not found
def getTopInteractionScoreMotifs(infile, outfile, max_num_edges, column_name):

    # Read the first line and determine if the column_name exists
    list_fields = getFieldsInFirstLine(infile, separator="\t")

    # Check that the column name exists
    if not column_name in list_fields:

        raise Exception("The column "+column_name+" is not found in the table.")

    # We are going to read the values_file in chunks at a time, this indicates the size of each chunk
    chunk_size = 1000000

    # First chunk
    first_chunk = True

    # Accumulate the top max_num_edges (interaction_score) for each dataframe in here
    # Then we sort and get the top interaction_score to get the final list
    top_row_chunks = None

    # Read a chunk of lines in the values_file, by default "NA" values are taken as na
    for chunk in pd.read_table(infile,
                               compression='gzip',
                               header=0,
                               delimiter="\t",
                               chunksize=chunk_size):

        if first_chunk:

            # Sort the values of the chunk by interaction_score
            # Get the top max_num_edges and add them to the top
            top_row_chunks = (chunk.sort_values(column_name,
                                                axis=0,
                                                ascending=False,
                                                inplace=False,
                                                kind='quicksort',
                                                na_position='last')).iloc[0:max_num_edges]

            first_chunk = False

        # If there are other chunks read, concatenate them by rows (note that the row index of each chunk comes from the original file),
        # so indexes can't overlap
        else:

            top_row_chunks = pd.concat([top_row_chunks,
                                        (chunk.sort_values(column_name,
                                                           axis=0,
                                                           ascending=False,
                                                           inplace=False,
                                                           kind='quicksort',
                                                           na_position='last')).iloc[0:max_num_edges]],
                                       axis=0)

    # Once everything is red, sort to get the final list

    # Sort by interaction score (inplace)
    top_row_chunks.sort_values(column_name,
                               axis=0,
                               ascending=False,
                               inplace=True,
                               kind='quicksort',
                               na_position='last')

    # Reindex to get index going from 0 -> num rows
    top_row_chunks = top_row_chunks.reset_index(drop=True)

    # Get the top max_num_edges hits
    top_row_chunks = top_row_chunks.iloc[0:max_num_edges]

    # Output the table with the header
    top_row_chunks.to_csv(outfile,
                  compression="gzip",
                  sep="\t",
                  header=True,
                  na_rep="NA",
                  mode="w",
                  index_label=None,
                  index=False,
                  line_terminator="\n")




def preventNodeTextOverlapping(position_dictionary,
                               top_bottom_nodes_correct,
                               distance_difference):

    # To avoid collisions in node texts, since the nodes will be in circles get the highest and the lowest nodes
    sorted_by_value_ascending = sorted(position_dictionary.items(), key=lambda kv: kv[1][1])

    highest_nodes = [sorted_by_value_ascending[0:top_bottom_nodes_correct]][0]

    lowest_nodes = [sorted_by_value_ascending[-(top_bottom_nodes_correct):]][0]



    # Get highest and lowest heights from the highest nodes and from the lowest nodes

    highest_nodes_heights = []

    for node in highest_nodes:
        # Get the node positions
        node_name, node_position = node

        highest_nodes_heights.append(node_position[1])

    max_highest_node_height = max(highest_nodes_heights)

    min_highest_node_height = min(highest_nodes_heights)

    lowest_nodes_heights = []

    for node in lowest_nodes:
        # Get the node positions
        node_name, node_position = node

        lowest_nodes_heights.append(node_position[1])

    max_lowest_node_height = max(lowest_nodes_heights)

    min_lowest_node_height = min(lowest_nodes_heights)





    initial_nodes_highest = {}

    # Order them left to right
    for node in highest_nodes:
        # Get the node positions
        node_name, node_position = node

        # Get the subset of keys from the initial dict
        initial_nodes_highest[node_name] = position_dictionary[node_name]

    sorted_by_left_right_highest = sorted(initial_nodes_highest.items(), key=lambda kv: kv[1][0])

    highest_nodes_left_right = [sorted_by_left_right_highest][0]


    initial_nodes_lowest = {}

    for node in lowest_nodes:
        # Get the node positions
        node_name, node_position = node

        # Get the subset of keys from the initial dict
        initial_nodes_lowest[node_name] = position_dictionary[node_name]

    sorted_by_left_right_lowest = sorted(initial_nodes_lowest.items(), key=lambda kv: kv[1][0])

    lowest_nodes_left_right = [sorted_by_left_right_lowest][0]


    node_counter = 1

    for node in lowest_nodes_left_right:

        # Get the node positions
        node_name, node_position = node

        # Put one high and one low
        if node_counter % 2 == 0:
            node_height = max_lowest_node_height - distance_difference
        else:
            node_height = max_lowest_node_height

        # Update node height
        old_node_pos = position_dictionary[node_name]

        new_node_pos = (old_node_pos[0], node_height)

        position_dictionary[node_name] = new_node_pos

        node_counter += 1





    node_counter = 1

    for node in highest_nodes_left_right:
        # Get the node positions
        node_name, node_position = node

        # Put one high and one low
        if node_counter % 2 == 0:
            node_height = max_highest_node_height
        else:
            node_height = max_highest_node_height + distance_difference

        # Update node height
        old_node_pos = position_dictionary[node_name]

        new_node_pos = (old_node_pos[0], node_height)

        position_dictionary[node_name] = new_node_pos

        node_counter += 1


    return position_dictionary





@cluster_runnable
def generateNetworkPlots(edges_file, tf_name, outfile, max_edges_to_show):

    top_row_chunks = pd.read_table(edges_file,
                                    compression='gzip',
                                    header=0,
                                    delimiter="\t")

    # If the table is empty don't do anything
    if top_row_chunks.shape[0] != 0:

        # Sort by interaction_score and get the max_edges_to_show
        top_row_chunks = top_row_chunks.sort_values("interaction_score",
                          axis=0,
                          ascending=False,
                          inplace=False,
                          kind='quicksort',
                          na_position='last').iloc[0:max_edges_to_show]


        # Get the old columns
        old_cols = top_row_chunks.columns

        new_cols = []

        # Replace "interaction_score" for "Weight"
        for col in old_cols:

            if col == "interaction_score":
                new_cols.append("Weight")

            else:
                new_cols.append(col)

        # Assign new columns
        top_row_chunks.columns = new_cols

        # Reduce the precision of the weight to 2 decimal places
        top_row_chunks.Weight = top_row_chunks.Weight.round(2)

        # Based on https://networkx.github.io/documentation/stable/auto_examples/drawing/plot_directed.html#sphx-glr-auto-examples-drawing-plot-directed-py
        # Author: Rodrigo Dorantes-Gilardi (rodgdor@gmail.com)
        directed_graph = nx.from_pandas_edgelist(top_row_chunks,
                                                 source='tf', # Source node column name
                                                 target='gene', # Target node column name
                                                 edge_attr="Weight", # Attributes to use for the edges
                                                 create_using=nx.DiGraph())


        # Increase space between nodes
        pos = nx.nx_agraph.graphviz_layout(directed_graph, prog='circo')
        #pos = nx.spring_layout(directed_graph, k=0.5, iterations=600)

        # If there are a lot of edges to show
        if max_edges_to_show >= 60:

            # Create a space between the top nodes and bottom nodes in the circle
            # so that the labels don't overlap
            pos = preventNodeTextOverlapping(pos,
                                             top_bottom_nodes_correct=15,
                                             distance_difference=50)





        # Get the weights on each edge
        labels = nx.get_edge_attributes(directed_graph, 'Weight')

        # Get each label weight
        edge_colors = list(labels.values())


        # node_sizes = [3 + 10 * i for i in range(len(directed_graph))]
        number_of_edges = directed_graph.number_of_edges()

        nodes = nx.draw_networkx_nodes(directed_graph,
                                       pos,
                                       node_color='red',
                                       node_size=100,
                                       alpha=0.5)  # node_size=node_sizes
        edges = nx.draw_networkx_edges(directed_graph, pos, arrowstyle='->',  # node_size=node_sizes,
                                       arrowsize=10, edge_color=edge_colors,
                                       edge_cmap=plt.cm.Wistia, width=1) # edge_cmap=plt.cm.Blues


        # Draw the labels on the nodes
        nx.draw_networkx_labels(directed_graph, pos, font_size=2)

        # Add the arrow labels
        #nx.draw_networkx_edge_labels(directed_graph, pos, labels)

        #nx.draw_networkx(directed_graph, pos, with_labels=True, node_size=100, font_size=6)

        # Set alpha value for each edge based on weight
        # Alpha will go from relative transparent in the least weight edge to totally opaque (1.0) in the highest weight edge

        # # Get the index of the original list when sorted by weight
        # # The index at position 0 from this list is the original edge list index for the lowest weight
        # weight_sorted_old_index = [b[0] for b in sorted(enumerate(edge_colors), key=lambda i: i[1])]
        #
        # # Get alphas to go from relative transparent 0.3 (pos 0) to 1.0 (pos number_of_edges) based on the number of edges
        # edge_alphas = [(0.5+ ((1-0.5)/(number_of_edges-1)*i)) for i in range(number_of_edges)]
        #
        # # Even alpha from 0 - 1.0
        # # edge_alphas = [(5 + i) / (number_of_edges + 4) for i in range(number_of_edges)]
        #
        # edge_alphas_pos = 0
        #
        # # Assign to the edges list (original) the corresponding alpha
        # for old_index in weight_sorted_old_index:
        #     edges[old_index].set_alpha(edge_alphas[edge_alphas_pos])
        #
        #     edge_alphas_pos += 1

        pc = mpl.collections.PatchCollection(edges, cmap=plt.cm.Wistia) #cmap=plt.cm.Blues
        pc.set_array(edge_colors)

        # Turn interactive plotting off, we don't want to show it
        plt.ioff()

        plt.colorbar(pc)

        ax = plt.gca()
        ax.set_axis_off()
        # plt.show()

        plt.savefig(outfile, dpi=300)

        # Clear the plot for the next one
        plt.clf()






# Copied from Motif_tools.motif_db_formatting_py3
#
# Inputs:
#     -infile: A jaspar motif format file. Such as:
#        /shared/sudlab1/General/apps/bio/PIQ_human/pwms/jasparfix.txt
#
# >MA0001.1;AGL3
# A  [ 0  3 79 40 66 48 65 11 65  0 ]
# C  [94 75  4  3  1  2  5  2  3  3 ]
# G  [ 1  0  3  4  1  0  5  3 28 88 ]
# T  [ 2 19 11 50 29 47 22 81  1  6 ]
# >MA0002.1;RUNX1
# A  [10 12  4  1  2  2  0  0  0  8 13 ]
# C  [ 2  2  7  1  0  8  0  0  1  2  2 ]
# G  [ 3  1  1  0 23  0 26 26  0  0  4 ]
# T  [11 11 14 24  1 16  0  0 25 16  7 ]
#
#    -outfile: An outfile with a header and format
#       file_position   gene
#
#    -field_start_db: For the database file specified, the field number (separated by ;) where the gene name begins.
#       Examples:
#       -HOCOMOCO, field start = 0 /shared/sudlab1/General/apps/bio/PIQ_human/pwms/HOCOMOCOv11_core_pwm_HUMAN_mono_jaspar_format.txt
#       >AHR_HUMAN.H11MO.0.B
#       41      11      22      3       1       3       0       0       43
#       18      12      44      1       150     1       3       0       67
#       56      35      21      146     1       149     1       154     16
#       39      96      67      4       2       1       150     0       28
#
#       Gets AHR_HUMAN.H11MO.0.B as gene name
#
#       -JASPAR, field_start = 1 /shared/sudlab1/General/apps/bio/PIQ_human/pwms/jasparfix.txt
#       >MA0001.1;AGL3
#       A  [ 0  3 79 40 66 48 65 11 65  0 ]
#       C  [94 75  4  3  1  2  5  2  3  3 ]
#       G  [ 1  0  3  4  1  0  5  3 28 88 ]
#       T  [ 2 19 11 50 29 47 22 81  1  6 ]
#
#       Gets AGL3 as gene name
#
#
# Outputs:
#     A table of motif pos (starting from 1 at the top) - gene
def getMotifGeneTablePIQJasparfix(infile, outfile, field_start_db):

    pos_motifs = 1

    # The output list of dictionaries ({pos, gene})
    dict_list = []

    with IOTools.open_file(infile, "r") as reader:

        for line in reader:

            # New motif
            if line.startswith(">"):

                # Remove the new line character from the line (avoids taking this into account downstream)
                line = line.rstrip('\n')

                # Remove the leading ">"
                line = line.lstrip(">")

                # New motif
                # Some motifs have the id and then the gene name with multiple ";" separated fields
                # Eg. MF0002.1;bZIP;CREB/G-box-like;subclass
                gene = line.split(";")[field_start_db:]

                gene = "_".join(gene)

                dict_motif_gene = {'pos': pos_motifs,
                                   'gene': gene}

                dict_list.append(dict_motif_gene)

                # Update the position
                pos_motifs += 1


    reader.close()

    # Output to the file as a tab separated table
    with IOTools.open_file(outfile, "w") as writer:

        # Output a header
        writer.write("file_position\tgene\n")

        for dict_motif_gene in dict_list:

            writer.write(str(dict_motif_gene['pos'])+"\t"+dict_motif_gene['gene']+"\n")

    writer.close()




# Copied from AuxiliaryPrograms/Genomic_regions/genomic_regions
# The infile at the moment must include header. It can be included the combination
# no header + column ids. See how header/no header is done in delete_rows_containing_character_any_cell_line
#
# Accepts compressed inputs, outputs compressed outputs if .gz file specified.
# Default mode (allowed=True): Goes row by row and checks in the column
# names provided (column_names) how many columns specified contain any of the
# strings specified in allowed_strings. Then it maintains the row if
# there are at least min_row_ocurrences or filters it out otherwise.
# If only rows containing any of the allowed_words in the specified columns are to be allowed:
# allowed = True, min_row_ocurrences = len(column_names)
#
# (allowed=False): In this case, the allowed_strings becomes not allowed strings,
# the rest strings become allowed and columns containing strings other than allowed_strings
# are marked as positive.
# Then it maintains the row if there are at least min_row_ocurrences positive occurrences
# or filters it out otherwise.
# If only rows containing not containing any of the allowed_words in the specified columns are to be allowed:
# allowed = False, min_row_ocurrences = len(column_names)
#
# Eg.
# row_no    chr     start   end     cell1   cell2   cell3
# 0         chr1    200     400       .       E1      E2
# 1         chr1    600     800      E1       .       E1
# 2         chr1    800    1000       .       E2      .
#
# column_names = ['cell1', 'cell3']
# allowed_words = ['E1', 'E2']
# allowed=True
# min_row_ocurrences = 2
# row 1: Eliminated because it only contains one specified columns (cell3) with either 'E1' or 'E2'
# row 2: Maintained because it contains two specified columns: cell1 and cell3 being 'E1' or 'E2'.
# row 3: Eliminated because it contains zero specified columns with either 'E1' or 'E2'.
#
# column_names = ['cell1', 'cell3']
# allowed_words = ['E1', 'E2']
# allowed=False (allowed_words become prohibited)
# min_row_ocurrences = 2
# row 1: Eliminated because it only contains one specified columns (cell1) without containing either 'E1' or 'E2'
# row 2: Eliminated because it contains zero specified columns without containing either 'E1' or 'E2'
# row 3: Maintained because it contains two specified columns (cell1 and cell3) with either 'E1' or 'E2'.
#
# Outputs the filtered file and a file with the filtered out rows.
# Does the process by chunks to save on memory.
#
#
# Inputs:
#     -infile: genomic regions and counts for each cell line. Must include header.
#        Format (including header):
#             chrom   start   end     cell_line1    cell_line2    cell_line3
#             chr1    156123600       156123800       .       E2       E3
#             chr1    160714000       160714200       E1    E1    E3
#
#     -allowed_strings: list of Strings (have to be Strings). Only these strings are allowed,
#         if any of these strings is not found (exact match as the
#         full string) in one of the cell lines, the row will be excluded.
#     -column_names: List of names of the columns in which to filter
#     -min_row_ocurrences: Each row must have at least this number of specified columns as positive.
#     -outfile: the name of the outfile
#     -outfile_filtered_out: the file to stored the filtered out rows
#     -allowed: by default (True) the behaviour is to allow only the strings provided
#             by setting this to False we disallow the strings provided in allowed_strings
#
# Outputs: writes the file in the specified outfile.
#
# Raises exception:
#     -if any of the columns contains a data type not int, float or string. Eg dates (for all the rows in the chunk).
#     -if any column name given can't be converted to column id (the column with the name doesn't exist).
@cluster_runnable
def delete_rows_not_containing_allowed_string_number_specified_cell_lines(infile,
                                                                          allowed_strings,
                                                                          column_names,
                                                                          min_row_ocurrences,
                                                                          outfile,
                                                                          outfile_filtered_out,
                                                                          allowed=True):
    # For debugging
    debug = False

    # For debugging
    rows_seen = 0

    # We are going to read the matrix file in chunks at a time, this indicates the size of each chunk
    chunk_size = 100000

    # To create the outfile we need to know when to append headers (first chunk)
    first_chunk = True

    # Read the infile by chunks specifying the header
    for chunk in pd.read_table(infile,
                               header=0,
                               chunksize=chunk_size,
                               delimiter='\t'):

        if debug:
            print
            "Rows seen: " + str(rows_seen) + " To: " + str(rows_seen + chunk_size - 1) + "\n"

            rows_seen += chunk_size

        # Extract all the rows of the selected columns, use loc with row/column names
        matrix_values_remaining = chunk.loc[:, column_names]

        # Go through all the cell line columns to check they are correct
        # and convert the columns
        for column_name in column_names:

            column_id = 0

            # Change the column_names to column_ids to check if the column exists
            try:
                column_id = chunk.columns.get_loc(column_name)
            except KeyError:
                raise Exception("The column " + column_name + " doesn't exist")

            # Convert each column to string to be compared
            # If all the numbers in the column are integers
            # to avoid being converted to floats before string
            # (https://stackoverflow.com/questions/22276503/how-to-i-change-data-type-of-pandas-data-frame-to-string-with-a-defined-format)
            if "int" in str(chunk.dtypes[column_id]):

                # Do in multiple steps to avoid http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
                # loc with column names
                matrix_values_remaining.loc[:, column_name] = (matrix_values_remaining.loc[:, column_name]).astype(int,
                                                                                                                   copy=False)

                matrix_values_remaining.loc[:, column_name] = (matrix_values_remaining.loc[:, column_name]).astype(
                    'str', copy=False)


            # If they are floats
            elif ("float" in str(chunk.dtypes[column_id]) or "object" in str(chunk.dtypes[column_id])):

                # If it's float, convert to string
                if "float" in str(chunk.dtypes[column_id]):
                    matrix_values_remaining.loc[:, column_name] = (matrix_values_remaining.loc[:, column_name]).astype(
                        'str', copy=False)

            else:

                # If its not int, float or object, raise an exception
                raise Exception("Type of column for id: " + str(column_id) + " not integer, float or string")

                # After the conversion filter rows

        # If the list of strings is allowed
        if allowed:
            # First we get for each cell in the table with all rows and columns whether it
            # contains one of the values in the allowed_strings
            # Element = True if it isn't in the allowed_strings
            table_match_value = matrix_values_remaining.isin(allowed_strings)

        # If the list of strings is not allowed
        else:
            # First we get for each cell in the table with all rows and columns whether it
            # doesn't contain one of the values in the allowed_strings
            # Element = True if it isn't in the allowed_strings
            table_match_value = ~matrix_values_remaining.isin(allowed_strings)

        # Then we see for every row if there are at least the minimum number of columns
        # positive (True) as required
        # table_match_value.sum(axis=1) -> Get the sum of matches (True) per row
        # table_match_value[table_match_value.sum(axis=1)>=min_row_ocurrences] -> get rows where the sum of matches is at least the threshold
        remaining_rows_names = (table_match_value[table_match_value.sum(axis=1) >= min_row_ocurrences]).index.tolist()

        # Get starting index names (not ids) for the whole dataframe
        starting_rows_names = chunk.index.tolist()

        # Get the difference: Starting - remaining as a set
        # (row names of the filtered out rows to be outputted)
        filtered_out_rows_names = list(set(starting_rows_names) - set(remaining_rows_names))

        # Recreate the output dataframe and the filtered out output dataframe
        # Since we have row names of the index, we use loc
        filtered_out_rows = chunk.loc[filtered_out_rows_names]

        remaining_rows = chunk.loc[remaining_rows_names]

        # Once done all filterings in the chunk, output it
        # Append headers
        if (first_chunk):

            # For the first chunk, output the header if it contains it
            remaining_rows.to_csv(outfile,
                                  sep="\t",
                                  header=True,
                                  compression='gzip',
                                  mode="w",
                                  index_label=None,
                                  index=False,
                                  line_terminator='\n')

            # For the first chunk, output the header if it contains it
            filtered_out_rows.to_csv(outfile_filtered_out,
                                     sep="\t",
                                     header=True,
                                     compression='gzip',
                                     mode="w",
                                     index_label=None,
                                     index=False,
                                     line_terminator='\n')

            first_chunk = False


        # Don't append headers
        else:

            remaining_rows.to_csv(outfile,
                                  sep="\t",
                                  header=False,
                                  compression='gzip',
                                  mode="a",
                                  index_label=None,
                                  index=False,
                                  line_terminator='\n')

            filtered_out_rows.to_csv(outfile_filtered_out,
                                     sep="\t",
                                     header=False,
                                     compression='gzip',
                                     mode="a",
                                     index_label=None,
                                     index=False,
                                     line_terminator='\n')



# Copied from AuxiliaryPrograms/Genomic_regions/genomic_regions
# Gets the unique values from a table in "infile" and the column specified by "column_name"
# If the file is too big, it reads it by chunks. It makes sure the values are unique among chunks by concatenating
# the unique values together on each pass and getting the unique values again.
#
# Inputs:
#   -infile: A file pointing to a table with a header.
#   -column_name: The column name to look for unique values in the table.
#
# Outputs:
#   -A list of unique values.
#
# Exception:
#   -If the column_name is not found in the table
def get_unique_values_in_colname(infile, column_name):

    # Read the first line and determine if the column_name exists
    list_fields = getFieldsInFirstLine(infile, separator="\t")

    # Check that the column name exists
    if not column_name in list_fields:
        raise Exception("The column " + column_name + " is not found in the table.")

    # We are going to read the values_file in chunks at a time, this indicates the size of each chunk
    chunk_size = 1000000

    # First chunk
    first_chunk = True

    # Accumulate the unique data for each dataframe in here
    uniq_value_chunks = None

    # Read a chunk of lines in the values_file, by default "NA" values are taken as na
    for chunk in pd.read_table(infile,
                               compression='gzip',
                               header=0,
                               delimiter="\t",
                               chunksize=chunk_size):

        # Extract all the rows of the selecte columns
        # It's one dimension so it becomes a series
        series_values_remaining = chunk.loc[:, column_name]

        # Remove the duplicates
        series_values_remaining = np.unique(series_values_remaining)

        if first_chunk:

            # Put the remaining values from the first chunk
            uniq_value_chunks = series_values_remaining

            first_chunk = False

        # Concatenate the values, the indexes come from the same starting file so indexes can't overlap
        else:

            uniq_value_chunks = np.concatenate([uniq_value_chunks,
                           series_values_remaining])

            # Drop the duplicates
            uniq_value_chunks = np.unique(uniq_value_chunks)

    # Once everything is red, turn into list

    # Return the list of values
    return uniq_value_chunks.tolist()


# NOTE: The relationship can be many to many (returned through multiple dictionaries of 1 element of each list each)
# NOTE: Empty strings in any list are ignored
# Compares all the elements from list_elements_contained with all the elements from list_elements_containing
# Gets all the pairs of elements where the element from list_elements_contained is completely contained in the element of
# list_elements_containing.
# For each pair, outputs the dictionary with the details.
#
#
# Inputs:
#   -list_elements_contained: A list of strings. Each of these elements is tested to be contained within each element
#       from list_elements_containing. This means that the string from this list is exactly the same as one from the other
#       or contained within one of the other.
#
#   -list_elements_containing: A list of strings. The strings here must fully contain a string from list_elements_contained.
#
#
# Outputs:
#   -A list of dictionaries of format:
#       {'element_contained':''
#        'element_containing':'' }
#
def get_elements_contained_in_list(list_elements_contained, list_elements_containing):

    output_list = []


    # Go through each list
    for element_contained in list_elements_contained:

        # Ignore empty strings
        if element_contained == "":

            continue

        for element_containing in list_elements_containing:

            # Ignore empty strings
            if element_containing == "":
                continue

            # If the element contained
            if element_contained in element_containing:

                dict = {}
                dict["element_contained"] = element_contained

                dict["element_containing"] = element_containing

                output_list.append(dict)

    return output_list



# Gets the files header_table_contained and header_table_containing. It gets the columns field_contained and field_containing
# respectively. Then sees for each element in header_table_contained if its contained in each element in header_table_containing.
# Outputting a file of contained relationships.
#
# Inputs:
#   -header_table_contained: A table with header containing a column with strings to check if they are contained.
#   -field_contained: The column to get the elements from.
#   -header_table_containing: A table with header containing a column with strings to check if they are containing the strings from header_table_contained.
#   -field_containing: The column to get the elements from.
#   -output_file: A file with header containing the successful contained relationships.
#
# Outputs:
#   -The table with header containing the successful relationships
#
# Exception:
#   -If any of the table doesn't have their corresponding column name
def generate_conversion_table_contained(header_table_contained,
                              field_contained,
                              header_table_containing,
                              field_containing,
                              output_file):


    # Check that the fields specified are present

    # Read the first line and determine if the column_name exists
    list_fields_contained = getFieldsInFirstLine(header_table_contained, separator="\t")

    # Check that the column name exists
    if not field_contained in list_fields_contained:
        raise Exception("The column " + field_contained + " is not found in the table "+header_table_contained)


    list_fields_containing = getFieldsInFirstLine(header_table_containing, separator="\t")

    # Check that the column name exists
    if not field_containing in list_fields_containing:
        raise Exception("The column " + field_containing + " is not found in the table " + header_table_containing)



    # Read the file
    contained_table = pd.read_table(header_table_contained,
                                           compression='gzip',
                                           header=0,
                                           delimiter="\t")


    # Get all the unique
    unique_contained = (contained_table.loc[:,field_contained]).unique().tolist()




    # Read the file
    containing_table = pd.read_table(header_table_containing,
                                    compression='gzip',
                                    header=0,
                                    delimiter="\t")

    # Get all the unique
    unique_containing = (containing_table.loc[:, field_containing]).unique().tolist()


    output_dict_list = get_elements_contained_in_list(list_elements_contained=unique_contained,
                                                              list_elements_containing=unique_containing)

    # Write the header
    header_output = "contained\tcontaining\n"

    with IOTools.open_file(output_file, "w") as writer:

        writer.write(header_output)

        # Write each dictionary
        for output_dict in output_dict_list:

            writer.write(output_dict['element_contained']+"\t"+output_dict['element_containing']+"\n")

    writer.close()






# Gets the files header_table_contained and header_table_containing. It gets the columns field_contained and field_containing
# respectively. Then sees for each element in header_table_contained if its contained in each element in header_table_containing.
# Outputting a file of contained relationships.
#
# Inputs:
#   -header_table_contained: A table with header containing a column with strings to check if they are contained.
#   -field_contained: The column to get the elements from.
#   -header_table_containing: A table with header containing a column with strings to check if they are containing the strings from header_table_contained.
#   -field_containing: The column to get the elements from.
#   -output_file: A file with header containing the successful contained relationships.
#
# Outputs:
#   -The table with header containing the successful relationships
#
# Exception:
#   -If any of the table doesn't have their corresponding column name
def generate_table_same_ids(header_table1,
                              field_table1,
                              header_table2,
                              field_table2,
                              output_file):


    # Check that the fields specified are present

    # Read the first line and determine if the column_name exists
    list_fields_1 = getFieldsInFirstLine(header_table1, separator="\t")

    # Check that the column name exists
    if not field_table1 in list_fields_1:
        raise Exception("The column " + field_table1 + " is not found in the table "+header_table1)


    list_fields_2 = getFieldsInFirstLine(header_table2, separator="\t")

    # Check that the column name exists
    if not field_table2 in list_fields_2:
        raise Exception("The column " + field_table2 + " is not found in the table " + header_table2)



    # Read the file
    table1 = pd.read_table(header_table1,
                           compression='gzip',
                           header=0,
                           delimiter="\t")


    # Get all the unique
    unique_1 = (table1.loc[:,field_table1]).unique().tolist()




    # Read the file
    table2 = pd.read_table(header_table2,
                           compression='gzip',
                           header=0,
                           delimiter="\t")

    # Get all the unique
    unique_2 = (table2.loc[:, field_table2]).unique().tolist()


    common_ids = list(set(unique_1).intersection(set(unique_2)))

    # Write the header
    header_output = "same_ids\n"

    with IOTools.open_file(output_file, "w") as writer:

        writer.write(header_output)

        # Write each dictionary
        for common_id in common_ids:

            writer.write(common_id+"\n")

    writer.close()




# Starts by getting all motifs in motif_gene_conv, goes through each motif tf and gets any genes symbols in tfs_contained_in_TSS
# and TSS_contained_in_tfs. Additionally it also gets any alternative TSS gene symbols from tss_alt_symbols and stores all the unique
# alternative names in tf_conversions_table with the corresponding TF (with the same name as in the motif database).
#
# Inputs:
#   -tfs_contained_in_TSS: tfs contained in TSS names produced automatically and then manually curated.
#       Format:
#           contained       containing
#           RUNX1   RUNX1T1
#           NR2F1   NR2F1-AS1
#
#   -TSS_contained_in_tfs: TSS contained in tfs: TSS contained in tfs names produced automatically and then manually curated.
#       Format:
#           contained       containing
#           RUNX1   RUNX1T1
#           NR2F1   NR2F1-AS1
#
#   -tss_alt_symbols: Primary gene symbol with alternative symbols in the line
#       Format: First column: main gene symbol. Following columns: alternative names
#       OR4F16
#       SAMD11  MRS
#       NOC2L   PPP1R112        NET15   NIR     NET7
#
#   -motif_gene_conv: Contains all the TF gene names from the motif database.
#       Format:
#       file_position   gene
#       1       AGL3
#       2       RUNX1
#
# Outputs:
#   -tf_conversions_table: Writes this table with the TF name (as in the motif database) in the first column and all the possible
#       alternative symbols in the next columns.
#       Format:
#       HAT5
#       T       SAVA    TBXT    T       TFT
#       br_Z1
def combine_tf_to_gene_symbol_table(tfs_contained_in_TSS,
                                    TSS_contained_in_tfs,
                                    tss_alt_symbols,
                                    motif_gene_conv,
                                    tf_conversions_table):

    # Read the motif gene conversions
    motif_gene_conv_table = pd.read_table(motif_gene_conv,
                                               compression='gzip',
                                               header=0,
                                               delimiter="\t")

    # Read the tfs contained in TSS
    tfs_contained_in_TSS_table = pd.read_table(tfs_contained_in_TSS,
                                               compression='gzip',
                                               header=0,
                                               delimiter="\t")

    # Read the TSS contained in tfs
    TSS_contained_in_tfs_table = pd.read_table(TSS_contained_in_tfs,
                                               compression='gzip',
                                               header=0,
                                               delimiter="\t")


    # Create a list with each element being a list of symbols
    alt_symbols_list = []

    with IOTools.open_file(tss_alt_symbols, "r") as reader:

        for line in reader:

            line = line.rstrip('\n')

            elements = line.split("\t")

            # If there isn't an empty line
            if not elements == [""]:

                alt_symbols_list.append(elements)

    reader.close()

    # Get all the tf involved
    unique_tfs = motif_gene_conv_table.gene.unique().tolist()

    with IOTools.open_file(tf_conversions_table, "w") as writer:


        # Go through each tf
        for unique_tf in unique_tfs:

            # First symbol (reference symbol to be outputted)
            ref_symbol = None


            # Get the rows of each table corresponding to that TF
            tfs_contained_in_TSS_table_tf = tfs_contained_in_TSS_table.loc[tfs_contained_in_TSS_table["contained"] == unique_tf]

            # By default create an empty list of symbols returned
            TSS1 = []

            # If rows are returned update the list of symbols
            if tfs_contained_in_TSS_table_tf.shape[0] != 0:

                # Get the corresponding TSS
                TSS1 = tfs_contained_in_TSS_table_tf.containing.unique().tolist()

                # Assign the first element
                ref_symbol = TSS1[0]




            # Get the rows of each table corresponding to that TF
            TSS_contained_in_tfs_table_tf = TSS_contained_in_tfs_table.loc[TSS_contained_in_tfs_table["containing"] == unique_tf]

            # By default create an empty list of symbols returned
            TSS2 = []

            # If rows are returned update the list of symbols
            if TSS_contained_in_tfs_table_tf.shape[0] != 0:

                # Get the corresponding TSS
                TSS2 = TSS_contained_in_tfs_table_tf.contained.unique().tolist()

                # Assign the first element
                if ref_symbol == None:
                    ref_symbol = TSS2[0]



            # Get the alternative gene symbols
            gene_symbols = []

            for alt_symbols in alt_symbols_list:

                if unique_tf in alt_symbols:

                    # For every instance found, concatenate the lists
                    gene_symbols = gene_symbols + alt_symbols

                    # Assign the first element
                    if ref_symbol == None:
                        ref_symbol = gene_symbols[0]


            # Concatenate all the gene symbols
            tf_conversions = TSS1 + TSS2 + gene_symbols

            # Get only the unique gene symbols
            tf_conversions = list(set(tf_conversions))

            # In the first column write the tf
            writer.write(unique_tf)

            # If there are other symbols write them tab separated
            if(len(tf_conversions) != 0):

                # Remove from the list the ref_symbol
                tf_conversions.remove(ref_symbol)

                # Write it first
                writer.write("\t" +ref_symbol)

                # If there are elements remaining
                if (len(tf_conversions) != 0):

                    writer.write("\t"+"\t".join(tf_conversions))

            # End of line
            writer.write("\n")




    writer.close()












def get_relevant_TFs_for_top_TF_gene_edges_per_TF(file_all_edges,
                                                  top_edges_tf_df,
                                                  tf,
                                                  tf_gene_symbols_table,
                                                  outfile,
                                                  tmp_dir):

    # Temp table with gene names
    tfs_as_genes_temp_file = (tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name

    # Temp table without gene names
    filtered_out_temp_file = (tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name

    # Create a df containing the top edges and if we find anything we will add
    # the rows with tf as TF and including also tf as Gene
    output_df = top_edges_tf_df

    # Create a column called "linking name" which contains the TF as in the "TF" column to link
    # the related nodes
    output_df.loc[:, 'linking_name'] = tf



    # TF and gene replacement symbol
    gene_replacing_symbol = None

    # First determine if there is a replacement symbol
    # Match the TF to get the gene symbol
    for tf_gene in tf_gene_symbols_table:

        # The first column contains the TF
        if tf_gene[0] == tf:

            # If it exists and it is the first coincidence in the table with a reference symbol,
            # use the reference gene symbol found to replace the TF name top hits and the any Gene found related
            if len(tf_gene) > 1:

                gene_replacing_symbol = tf_gene[1]

                break


    if gene_replacing_symbol != None:

        # Start by replacing the top interactions whose TF is in the list
        output_df.loc[:, 'tf'] = gene_replacing_symbol


    # Match the TF to get the gene symbol
    for tf_gene in tf_gene_symbols_table:

        # The first column contains the TF
        if(tf_gene[0] == tf):


            # Now we get all the genes where the symbol is in the list of genes overlapping tfs
            delete_rows_not_containing_allowed_string_number_specified_cell_lines(file_all_edges,
                                                                                  allowed_strings=tf_gene,
                                                                                  column_names=["gene"],
                                                                                  min_row_ocurrences=1,
                                                                                  outfile=tfs_as_genes_temp_file,
                                                                                  outfile_filtered_out=filtered_out_temp_file,
                                                                                  allowed=True)



            # Read the relevant edges
            relevant_edges = pd.read_table(tfs_as_genes_temp_file,
                                                   compression='gzip',
                                                   header=0,
                                                   delimiter="\t")

            # If there are rows to add
            if relevant_edges.shape[0] != 0:

                # Create a column called "linking name" which contains the TF as in the "TF" column to link
                # the related nodes
                relevant_edges.loc[:, 'linking_name'] = tf

                if gene_replacing_symbol != None:

                    # Replace the relevant edges genes with the symbol
                    relevant_edges.loc[:, "gene"] = gene_replacing_symbol


                # Concatenate (attach rows) of the top edges (with which we begin) corresponding to the TF with the relevant
                # (TF as gene)
                output_df = pd.concat([output_df,
                            relevant_edges],
                            axis=0)


            break

    # Delete the temporal files created
    os.unlink(tfs_as_genes_temp_file)

    os.unlink(filtered_out_temp_file)


    # Whether or not we found additional edges, output the file containing the top + any found
    # Output the table with the header
    output_df.to_csv(outfile,
                          compression="gzip",
                          sep="\t",
                          header=True,
                          na_rep="NA",
                          mode="w",
                          index_label=None,
                          index=False,
                          line_terminator="\n")





# 1) The maximum number of edges are extracted from the table. An edge is a "tf" - "gene" row
# 2) For each "tf" column present in the list in 1) all edges where this TF appears as "gene" are
# extracted and combined with the list in 1)
# 3) The combined list is removed from duplicates
# Inputs:
#   -infile: Format:
#       gene    tf      interaction_score       gene_desc
#       CFH     AGL3    0.476420431210367       complement factor H
#       STPG1   AGL3    0.319511602570206       sperm tail PG-rich repeat containing 1
#       HS3ST1  AGL3    0.219486570405203       heparan sulfate-glucosamine 3-sulfotransferase 1
#
#
@cluster_runnable
def get_relevant_TFs_for_top_TF_gene_edges(interaction_scores_file,
                                       tf_gene_symbols_conv_file,
                                       outdir,
                                       max_num_edges,
                                       tmp_dir):

    # Create a list of dictionaries, one for each TF where we specify the TF and the file containing all the
    # relevant edges for that TF
    tf_relevant_edges_dict_list = []

    # Temp table with gene names
    temp_top_edges = (tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name

    # The first step is to get the top edges (tf - gene) rows by "interaction_score"
    getTopInteractionScoreMotifs(interaction_scores_file, temp_top_edges, max_num_edges, "interaction_score")

    # Read the top interaction_score edges
    top_interaction_scores = pd.read_table(temp_top_edges,
                              compression='gzip',
                              header=0,
                              delimiter="\t")

    # Delete the temp file
    os.unlink(temp_top_edges)

    # If we don't have any top interaction scores, return the empty dictionary
    if top_interaction_scores.shape[0] == 0:

        return tf_relevant_edges_dict_list


    # Get all the tf involved
    unique_tfs = top_interaction_scores.tf.unique().tolist()

    # Parse the file of TF - gene symbols to search
    # We create a list with a list of gene symbols where the first element is the TF motif name
    alt_symbols_list = []

    with IOTools.open_file(tf_gene_symbols_conv_file, "r") as reader:

        for line in reader:

            line = line.rstrip('\n')

            elements = line.split("\t")

            # If there isn't an empty line
            if not elements == [""]:
                alt_symbols_list.append(elements)

    reader.close()




    for unique_tf in unique_tfs:

        # Create a file to store the relevant edges table
        edges_file = os.path.join(outdir, (unique_tf+".tsv.gz"))

        # Subset the edges from the top interaction score corresponding to that TF
        top_interaction_scores_tf_df = top_interaction_scores.loc[top_interaction_scores["tf"] == unique_tf]

        get_relevant_TFs_for_top_TF_gene_edges_per_TF(file_all_edges=interaction_scores_file,
                                                      top_edges_tf_df=top_interaction_scores_tf_df,
                                                      tf=unique_tf,
                                                      tf_gene_symbols_table=alt_symbols_list,
                                                      outfile=edges_file,
                                                      tmp_dir=tmp_dir)

        tf_relevant_edges_dict = {}

        tf_relevant_edges_dict["tf"] = unique_tf

        tf_relevant_edges_dict["file"] = edges_file

        tf_relevant_edges_dict_list.append(tf_relevant_edges_dict)


    return tf_relevant_edges_dict_list





def generate_top_rows_histogram(dataframe_file,
                                sample,
                                num_top_hits,
                                outfile,
                                x_lab,
                                y_lab,
                                title,
                                tmp_dir):

    # Temp table with top numbers
    temp_top_rows = (tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name

    # Get the top rows in the sample by number
    getTopInteractionScoreMotifs(infile = dataframe_file,
                                 outfile = temp_top_rows,
                                 max_num_edges = num_top_hits,
                                 column_name = sample)

    # Read the top interaction_score edges
    top_rows = pd.read_table(temp_top_rows,
                             compression='gzip',
                             header=0,
                             delimiter="\t")

    # Delete temp file
    os.unlink(temp_top_rows)

    # Sort them
    top_rows = top_rows.sort_values(sample,
                                    axis=0,
                                    ascending=False,
                                    inplace=False,
                                    kind='quicksort',
                                    na_position='last').iloc[0:num_top_hits]

    # Turn interactive plotting off, we don't want to show it
    plt.ioff()


    ax = top_rows.plot(x = "TF", y = sample, kind='bar', title=title, figsize=(15, 10), legend=True, fontsize=12)
    ax.set_xlabel(x_lab, fontsize=12)
    ax.set_ylabel(y_lab, fontsize=12)



    # Save the file, remove whitespace around the image
    plt.savefig(outfile)

    # Clear the plot for the next one
    plt.clf()



