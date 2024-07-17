#### CREATE A SCRIPT FOR OBTAINING THE ANNOTATIONS OF REGIONS IDENTIFIED FOR NUCLEOTIDE DIVERSITY AND FST ####

import sys
import pandas as pd
import re

file_for_analyse = sys.argv[1]
gtf_genome = sys.argv[2]

#Function for load the file as a dataframe for pandas:

def load_gtf_file(gtf_genome):
    # Create an empty dictionary to store the columns
    gtf_df = {"CHROM": [], "FEATURE": [], "START": [], "END": [], "STRAND": [], "FRAME": [], "INFO": []}
    
    # Open the GTF file
    with open(gtf_genome, 'r') as gtf_file:
        # Iterate through each line in the file
        for line in gtf_file:
            line = line.strip()
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            # Split the line by tabs to extract column values
            columns = line.split("\t")
            
            # Extract and append values to the respective columns in the dictionary
            #chr = columns[0]
            gtf_df["CHROM"].append(columns[0])
            #feature = columns[2]
            gtf_df["FEATURE"].append(columns[2])
            #start_pos = columns[3]
            gtf_df["START"].append(int(columns[3]))
            #end_pos = columns[4]
            gtf_df["END"].append(int(columns[4]))
            #strand = columns[6]
            gtf_df["STRAND"].append(columns[6])
            #frame = columns[7]
            gtf_df["FRAME"].append(columns[7])
            #information = columns[8]
            gtf_df["INFO"].append(columns[8])
    
    # Convert the dictionary to a DataFrame
    gtf_df = pd.DataFrame(gtf_df)
    
    return gtf_df

 #Create a function for extract and add important information: gene_id, transcript_id, annotation

def create_information_gtf(gtf_df):
    gene_ids = []
    transcript_ids = []
    annotations = []
    protein_ids = []
    #Iterate over each information data
    for info_data in gtf_df["INFO"]:
        # Established a empty value as first
        gene_id, transcript_id, annotation, prot_id = None, None, None, None
        #Divide the information parts into elements
        info_parts = info_data.split(";")
        for element in info_parts:
            if "gene_id" in element:
                gene_id = element.split('"')[1]
            elif "transcript_id" in element and not "orig_transcript_id" in element:
                transcript_id = element.split('"')[1]
            elif "note" in element:
                annotation = element.split('"')[1]
            elif "protein_id" in element and not "orig_protein_id" in element:
                prot_id = element.split('"')[1]
        #Add the information to the list for each values
        gene_ids.append(gene_id)
        transcript_ids.append(transcript_id)
        annotations.append(annotation)
        protein_ids.append(prot_id)

    #Create the columns from the list
    gtf_df["GENE_ID"] = gene_ids
    gtf_df["TRANSCRIPT_ID"] = transcript_ids
    gtf_df["NOTE"] = annotations
    gtf_df["PROTEIN_ID"] = protein_ids
    #Add NAs to the columns that are empty
    gtf_df.fillna(value=pd.NA, inplace=True)
    del gtf_df["INFO"]

    return gtf_df


#Function for extracting transcript_ids associated with the coordinates

def extracting_transcripts_ids(file_df, gtf_df):
    # Short the dataframes by chromosome and position for using the merge_asof
    file_df = file_df.sort_values(by=['CHROM', 'BIN_START'])
    gtf_df = gtf_df.sort_values(by=['CHROM', 'START'])

    # Convert the chromosome columns to a string
    file_df['CHROM'] = file_df['CHROM'].astype(str)
    gtf_df['CHROM'] = gtf_df['CHROM'].astype(str)

    # Create a index column for maintaining the original order
    file_df['index'] = file_df.index

    # Create an empty list for each row for the column TRANSCRIPT_ID
    file_df['TRANSCRIPT_ID'] = [None for _ in range(len(file_df))]

    # Use merge_asof for obtaining coincidences based on intervals
    #For each chromosome in the dataframe, subset a dataframe for the file (by chromosome) and another for the gtf file (for the entries with that chromosome and for transcripts)
    for chrom in file_df['CHROM'].unique():
        file_chr_df = file_df[file_df['CHROM'] == chrom]
        gtf_chr_df = gtf_df[(gtf_df['CHROM'] == chrom) & (gtf_df['FEATURE'] == 'transcript')]

        # Create merge dataframe between the two subsets created, doing the merge focusing on the intervals
        merged_df = pd.merge_asof(file_chr_df, gtf_chr_df, left_on='BIN_START', right_on='START', direction='forward', suffixes=('', '_gtf'))

        # Once the merge is done, go row by row extracting the transcripts within each interval, and add it to the list of the original dataframe
        for idx, row in merged_df.iterrows():
            transcripts = gtf_chr_df[(gtf_chr_df['START'] <= row['BIN_END']) & (gtf_chr_df['END'] >= row['BIN_START'])]['TRANSCRIPT_ID'].tolist()
            if len(transcripts) == 0:
                file_df.at[row['index'], 'TRANSCRIPT_ID'] = None
            else:
                file_df.at[row['index'], 'TRANSCRIPT_ID'] = ",".join(transcripts)

    # After cathing all the transcripts within the intervals, eliminate the index column
    file_df = file_df.drop(columns=['index'])

    return file_df

# Create a function for identifying transcripts that are classify only with one significance

def obtain_the_transcripts_with_an_unique_significance(file_df):
    # Create a new df with all the transcripts separate by rows:
    extracted_df = file_df[["TRANSCRIPT_ID", "Significance"]].copy()
    extracted_df["TRANSCRIPT_ID"] = extracted_df["TRANSCRIPT_ID"].str.split(",")
    extracted_df = extracted_df.explode("TRANSCRIPT_ID")
    extracted_df = extracted_df.dropna(subset="TRANSCRIPT_ID")
    #print(extracted_df.head(50))
    print()
    #Group bu ynique values of significances
    transcripts_number_of_significances = extracted_df.groupby("TRANSCRIPT_ID")["Significance"].nunique()
    #print(transcripts_number_of_significances.head(100))
    #Filter for obtaining the transcripts with more that one significance and the transcript with a unique significance
    transcripts_with_multiple_significance = transcripts_number_of_significances[transcripts_number_of_significances > 1].reset_index()
    transcript_with_one_significance = transcripts_number_of_significances[transcripts_number_of_significances == 1].reset_index()
    print("Number of transcripts with multiple significance:", len(transcripts_with_multiple_significance), sep ="\t")
    print("Number of transcripts with one significance:", len(transcript_with_one_significance), sep = "\t")
    #Obtain a dataframe with the transcripts with a unique significance
    transcript_with_one_significance.reset_index(drop=True, inplace=True)
    #transcripts_ids_multiple_significance = pd.DataFrame(transcripts_with_multiple_significance)
    print()
    #print(transcript_with_one_significance.head(30))
    return transcript_with_one_significance

# Create a function for obtain a dataframe with the transcripts that have the selected significance value
def obtain_transcripts_selected_for_significance(new_file_df, filter_value, transcripts_with_unique_signif):
    #Obtain a subset based on the value given for filtering by significance
    selected_df = new_file_df[new_file_df["Significance"] == filter_value].copy()
    #Obtain the values separated by "," to list
    selected_df["TRANSCRIPT_ID"] = selected_df["TRANSCRIPT_ID"].str.split(",")
    #Explode this column in order to have as many rows as transcripts have each row
    new_df = selected_df.explode("TRANSCRIPT_ID")
    #Eliminate intervals that do not have transcripts
    new_df = new_df.dropna(subset="TRANSCRIPT_ID")
    #print(new_df.head(50))
    # Calculate the number of unique transcripts
    number_unique_transcripts = new_df["TRANSCRIPT_ID"].unique()
    print()
    print("The number of unique transcripts from the file filtered is: ", str(len(number_unique_transcripts)), sep="")
    #Do a left_merge in order to retain only the transcripts which have been annotated with an unique significance
    df_transcripts_selected = pd.merge(new_df, transcripts_with_unique_signif, on= "TRANSCRIPT_ID", how="left")
    #Due to is a left join, keep only the rows that do not have NAs
    df_transcripts_selected = df_transcripts_selected.dropna()
    df_transcripts_selected = df_transcripts_selected.drop_duplicates(subset="TRANSCRIPT_ID")
    number_unique_transcripts_filtered = df_transcripts_selected["TRANSCRIPT_ID"].unique()
    print()
    print("The number of transcripts filtered that only have a significance value is: ", str(len(number_unique_transcripts_filtered)),sep="")
    print()
    return df_transcripts_selected
    

#Function for adding columns with notes, protein_ids and Prudu
def adding_extra_information(df_transcripts_selected, gtf_df):
    # Obtain a subset with the notes, the protein_ids
    subset_notes = gtf_df[["TRANSCRIPT_ID", "NOTE"]].dropna()
    subset_notes = subset_notes.drop_duplicates(subset="TRANSCRIPT_ID")
    subset_proteins = gtf_df[["TRANSCRIPT_ID", "PROTEIN_ID"]].dropna(subset="PROTEIN_ID")
    subset_proteins = subset_proteins.drop_duplicates(subset="PROTEIN_ID")
    # Once we get that, do a join based on the transcript id
    subset_notes_proteins = pd.merge(subset_proteins, subset_notes, on="TRANSCRIPT_ID", how="left")
    #print(subset_notes_proteins.head(100))
    # Once we had that, do the merge
    merged_df = pd.merge(df_transcripts_selected, subset_notes_proteins, on="TRANSCRIPT_ID", how="left")
    #Check the number of rows of the previous df and the merge
    n_rows_previous = df_transcripts_selected.shape[0]
    n_rows_merged = merged_df.shape[0]
    print("Number of rows from the previous dataframe: ", str(n_rows_previous), sep ="")
    print("Number of rows from the merged dataframe:", str(n_rows_merged), sep="")
    # Create thr Prudu id, using the function created
    merged_df = create_prudu_id_column_for_the_dataset(merged_df)
    return merged_df


#For having a enrichment analysis of GO terms, we need to create a new variable which contains the PRUDU id (contained in TRANSCRIPT_ID)
def create_prudu_id_column_for_the_dataset(df):
    # Create the new variable and asign the value
    df["PRUDU"] = df["TRANSCRIPT_ID"].apply(lambda x: x.split("|")[2] if isinstance(x,str) else None)
    return df


# We create a function to upload the fastas in a dataframe
def create_dataframe_fasta(protein_fasta):
    #We create an empty dataframe, with two columns, one for the header and the second for the adn sequence
    fasta_df = {"Id" : [], "Seq" : []}
    with open(protein_fasta, "r") as fasta_file:
        id = None
        seq = ""
        for line in fasta_file:
            line = line.strip()
            #If it is a header
            if line.startswith(">"):
                #If we find another header, we save the previous one and its sequence
                if id is not None:
                    fasta_df['Id'].append(id)
                    fasta_df['Seq'].append(seq)
                #After we establish the new data, with a new id
                id = line[1:].split(" ")[0]
                seq = ""
            #If it is not a header, its the sequence of the header
            else:
                seq += line
        #We add the last fasta sequence
        if id is not None:
            fasta_df['Id'].append(id)
            fasta_df['Seq'].append(seq)
    return pd.DataFrame(fasta_df)

#Create a function for obtaining the sequence of the selected proteins
def obtain_sequences_from_selected_proteins(final_df, fasta_df):
    #First obtain a list with all the proteins selected
    proteins_selected = final_df["PROTEIN_ID"].dropna()
    proteins_selected = proteins_selected.drop_duplicates()
    list_proteins_selected = proteins_selected.tolist()
    #Create an output file 
    salida_file = input("Enter the file name to save the fasta (with extension, e.g., '.fasta): ")
    #Now write the fasta file
    print()
    print("The fasta is being written")
    print()
    with open(salida_file, "w") as salida:
        for protein in list_proteins_selected:
            row = fasta_df.loc[fasta_df["Id"] == protein]
            #print(row)
            sequence = row["Seq"].values[0]
            salida.write(">" + protein + "\n")
            salida.write(sequence + "\n")



##### MAIN PROGRAM #####

#load the csv file:
file_df = pd.read_csv(file_for_analyse)
print(file_df.columns)
gtf_df = load_gtf_file(gtf_genome)
gtf_df = create_information_gtf(gtf_df)
print()
new_file_df = extracting_transcripts_ids(file_df, gtf_df)

print(new_file_df.head(10))

#print(new_file_df.head(30))

print()
transcripts_with_unique_signif = obtain_the_transcripts_with_an_unique_significance(new_file_df)


#Once it is done, ask the user if he wants to filter by some value of a column
# Once all the files are loaded and processed, create a loop for the program
cont = 0

while cont != 4:
    print("==================================================")
    print()
    #print each column of the dataframe
    
    # Choose the desired option
    print("Options: ")
    print()
    print("\t" + "1. Obtain the note of each transcript filtering by Significance")
    print("\t" + "2. End the program")
    print()

    cont = int(input("Introduce the desired option: "))
    
    if cont == 1:
        # Filter by column and value
        values_column_selected = new_file_df["Significance"].unique()
        print("These are the values for filtering in this column", values_column_selected, sep = "\n")
        filter_value = input("Introduce the value that do you want to keep: ")
        # Obtain the transcripts that have that value for significance
        df_transcripts_selected = obtain_transcripts_selected_for_significance(new_file_df, filter_value, transcripts_with_unique_signif)
        #print(df_transcripts_selected.head(50))
        #Now add the note, the protein id, and the Prudu id from the gft file associated with the transcript
        final_df = adding_extra_information(df_transcripts_selected, gtf_df)
        print(final_df.head(50))
        #Now save the dataframe
        rute_and_name_file = input("Write the route and the name file for saving this information: ")
        final_df.to_csv(rute_and_name_file, index=False)
        #Ask the user if he wants to obtain also a fasta file with all the proteins found
        ask_for_fasta = str(input("Do you want to obtain a fasta ptotein file with the protein ids found? Write yes or not: ")).upper()
        #If the user wants to obtain the fasta file, use the function to load the fasta file and obtain the file with the proteins_ids and sequences
        if ask_for_fasta == "YES":
            protein_fasta = str(input("Write the path and the file of the fasta protein file: "))
            fasta_df = create_dataframe_fasta(protein_fasta)
            obtain_sequences_from_selected_proteins(final_df, fasta_df)

    if cont == 2:
        cont = 4
        break





