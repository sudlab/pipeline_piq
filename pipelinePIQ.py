import pandas as pd

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
 
    # For each table specify the column name for the first column
    bed_table_columns = list(bed_table.columns)
     
    bed_table_columns[0] = "hit_id"
     
    bed_table.columns = bed_table_columns
     
     
    id_table_columns = list(id_table.columns)
     
    id_table_columns[0] = "hit_id"
     
    id_table.columns = id_table_columns
     
    result = pd.merge(id_table, bed_table, how="left", on=["hit_id"])
    
    # Drop the hit_id column
    result = result.drop('hit_id', 1) 
     
    # Output the table appending
    result.to_csv(output_file,
                 compression="gzip",
                 sep="\t",
                 header=None,
                 mode="w",
                 index_label=None,
                 index=False,
                 line_terminator="\n")
              
    
    