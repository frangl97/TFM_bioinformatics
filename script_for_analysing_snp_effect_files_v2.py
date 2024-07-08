####### SCRIPT FOR OBTAINING THE ANNOTATIONS OF THE SNP EFFECT ##################

import pandas as pd
#import matplotlib.pyplot as plt

import sys

import time

#Load the file

snp_effect_file = sys.argv[1]
gft_genome = sys.argv[2]
name_samples_file = sys.argv[3]

# Function for load the sample names into a list

def create_list_samples(name_samples_file):
    with open(name_samples_file, 'r') as samples:
        list_samples = []
        for line in samples:
            line = line.split()
            #Add the sample to the list
            list_samples.append(line)
        return list_samples


# this function will recognize the headers of the VCF file
def select_columns_and_samples(snp_effect_file):
    list_headers = []
    #for each line of the file, we only want to know the name of the columns of the vcf file
    with open(snp_effect_file, 'r') as vcf:
        for line in vcf:
            line = line.strip()
            #Select the line with the headers "#"
            if line[0] == "#" and line[1] != "#":
                columns = line.split("\t")
                # For each element of the headers, each one sould be load in a list
                for element_header in columns:
                    # If it is the first element, the "#" have to be eliminated
                    if element_header[0] == "#":
                        element_header = element_header[1:]
                    #We add the header into a list
                    list_headers.append(element_header)
                # Stop reading the file if we had the headers
                break
    return list_headers


# Function for creating the dataframe of the snp effect file with pandas
def create_database(snp_effect_file, list_headers):

    # Open the vcf file and add to a list all the lines that are not headers or information (contains '#')
    data = []
    with open(snp_effect_file, 'r') as vcf:
        for line in vcf:
            line = line.strip()
            if not line.startswith("#"):
                data.append(line.strip().split("\t"))
    vcf_eff_df = pd.DataFrame(data, columns=list_headers)
    return vcf_eff_df



# Create a function in order to know the type of variant identificated for each one


def identifying_variants(vcf_eff_df):
    # Split the INFO column once and store the result to avoid repetitive splitting
    info_split = vcf_eff_df["INFO"].str.split("|")
    ref_info = vcf_eff_df["REF"]
    alt_info = vcf_eff_df["ALT"]
    # Extract information and assign to respective columns
    vcf_eff_df["ANNOTATION"] = info_split.str[1]
    vcf_eff_df["ANNOTATION_IMPACT"] = info_split.str[2]
    vcf_eff_df["GENE_NAME"] = info_split.str[3]
    vcf_eff_df["GENE_ID"] = info_split.str[4]
    possible_transcript = info_split.str[5]

    # Extract transcript data and concatenate
    data_transcript = info_split.str[6].astype(str) + "|" + info_split.str[7].astype(str) + "|" + info_split.str[8].astype(str)
    # Assign transcript data to TRANSCRIPT_ID column where possible_transcript is "transcript"
    vcf_eff_df.loc[possible_transcript == "transcript", "TRANSCRIPT_ID"] = data_transcript

    # Established if the variant is a SNP or an INDEL
    vcf_eff_df["TYPE_VARIANT"] = vcf_eff_df.apply(lambda row: 'INDEL' if len(row['REF']) != len(row['ALT']) else 'SNP', axis=1)

    #Eliminate the INFO column once all the needed information has been extracted
    del vcf_eff_df["INFO"]
    return vcf_eff_df

# Load the information of the gtf file

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


# Create a function for extract and add important information: gene_id, transcript_id, annotation

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


def show_information_about_samples(list_samples):
    print()
    print("Number of samples analysed:", len(list_samples))
    print()
    print("Name of samples:")
    print()
    for element in list_samples:
        print(element)
    print()


# This function allows the identification of different types of variants identificated, as well as the different number of them by chromosome, and by sample
def calculate_type_annotation_variants(vcf_eff_df):
    # First we select the columns of the IDs and the annotations
    selected_information = vcf_eff_df[["ID", "ANNOTATION"]]
    # Group the data by its annotation feature
    effect_annotation_counts = selected_information.groupby("ANNOTATION").count()
    effect_annotation_counts = effect_annotation_counts.rename(columns={"ID" : "Counts"})
    #print(selected_information)
    print(effect_annotation_counts)
    # If the user wants to write the information into a csv file
    print()
    save_or_not = str(input("Do you want to save this information? Indicate yes or not ")).upper()
    #Save the information
    if save_or_not == "YES":
        rute_and_name_file = input("Write the route and the name file for saving this information: ")
        effect_annotation_counts.to_csv(rute_and_name_file, index=True)

# This function allows the counting of the different annotations and types of variants in the snpEffect file
def calculate_annotation_impact(vcf_eff_df):
    # First we select the columns of the IDs and the annotations
    selected_information = vcf_eff_df[["ID", "ANNOTATION_IMPACT"]]
    
    # Group the data by type of variant and annotation, then count the occurrences
    annotation_impact_counts = selected_information.groupby("ANNOTATION_IMPACT").size().reset_index(name="COUNT")
    annotation_impact_counts = annotation_impact_counts.rename(columns={"ID" : "Counts"})
    print(annotation_impact_counts)
    # If the user wants to write the information into a csv file
    save_or_not = input("Do you want to save this information? Indicate 'yes' or 'no': ").strip().lower()
    if save_or_not == "yes":
        file_name = input("Enter the file name to save the information (with extension, e.g., 'output.csv'): ")
        annotation_impact_counts.to_csv(file_name)
    
    return annotation_impact_counts
    

# This function allows the counting of the different annotations in the snpEffect file
def calculate_type_variants(vcf_eff_df):
    type_variants = selected_information = vcf_eff_df[["ID", "TYPE_VARIANT"]]
    # Group the data by its annotation feature
    type_variants_counts = type_variants.groupby("TYPE_VARIANT").count()
    type_variants_counts = type_variants_counts.rename(columns={"ID" : "Counts"})
    print()
    print(type_variants_counts)
    print()

# This function allows the counting of the type of variants found in the snpEffect classified by chromosomes
def calculate_type_variants_by_chromosome(vcf_eff_df):
    type_variants = selected_information = vcf_eff_df[["ID", "CHROM", "TYPE_VARIANT"]]
    type_variants_per_chr_counts = type_variants.groupby(["CHROM", "TYPE_VARIANT"]).size().reset_index(name="COUNT")
    
    # Pivot the table to have types of annotation as columns
    type_variants_per_chr_counts_pivot = type_variants_per_chr_counts.pivot(index="CHROM", columns="TYPE_VARIANT", values="COUNT")
    print(type_variants_per_chr_counts_pivot)
    # If the user wants to write the information into a csv file
    save_or_not = input("Do you want to save this information? Indicate 'yes' or 'no': ").strip().lower()
    if save_or_not == "yes":
        file_name = input("Enter the file name to save the information (with extension, e.g., 'output.csv'): ")
        type_variants_per_chr_counts_pivot.to_csv(file_name)


# Create a function for combining the two dataframes
def combine_dataframes_and_eliminate_info_columns(vcf_eff_df, gtf_new):
    #First we need to obtain a datraframe only with feature = transcript, and eliminate he protein id column
    selected_gtf_transcripts = gtf_new[gtf_new["FEATURE"] == "transcript"]
    del selected_gtf_transcripts["PROTEIN_ID"]
    #print("transcript_database")
    #print(selected_gtf_transcripts)
    #Secondly, select all the cds in order to keep the protein ids associated with each transcript
    selected_gtf_proteins = gtf_new[["TRANSCRIPT_ID", "PROTEIN_ID"]]
    #Eliminate all the transcripts that do not have a protein id and eliminate duplicates
    selected_gtf_proteins = selected_gtf_proteins.drop_duplicates(subset=["PROTEIN_ID"])
    selected_gtf_proteins = selected_gtf_proteins.dropna(subset=["PROTEIN_ID"])
    #print("protein_database")
    #print(selected_gtf_proteins)

    # Combine these information in order to have a dataframe with all the transcripts, their annotations if the have and the protein id if the have left merge on "TRANSCRIPT_ID"
    gtf_transcripts = pd.merge(selected_gtf_transcripts, selected_gtf_proteins, on= "TRANSCRIPT_ID", how="left")
    #print("merge")
    #print(gtf_transcripts)

    # Obtend the name of the columns for this new dataframe gtf_transcripts in order to view if there is other columns repeated apart from transcript_id with the vcf_eff_df
    gtf_transcripts_columns = gtf_transcripts.columns
    vcf_eff_df_columns = vcf_eff_df.columns

    list_repet_gtf_transcripts = []
    for col in gtf_transcripts_columns:
        if col in vcf_eff_df_columns:
            list_repet_gtf_transcripts.append(col)
    
    #now eliminate all the repeat columns with the exception of TRANSCRIPT_ID in order to the merge
    for element in list_repet_gtf_transcripts:
        if element != "TRANSCRIPT_ID":
            del gtf_transcripts[element]

    #Combine the vcf_effect_df with the gft_transcripts

    combine_data_frame = pd.merge(vcf_eff_df, gtf_transcripts, on= "TRANSCRIPT_ID", how = "left")

    # Add the Prudu_id column
    combine_data_frame = create_prudu_id_column_for_the_dataset(combine_data_frame)

    return combine_data_frame

# Create a function for provide the user all the possible information that could be viewed

def create_dataframe_for_viewing_information(combine_data_frame, list_samples):
    #Firstly, show the name of the combine dataframe and the others, and its headers, as well as the samples names
    columns_combine_data_frame = combine_data_frame.columns
    print("------------------------------------------------------------------")
    print("Combine_data_frame")
    print( "\t" + "Columns:" )
    print()
    for column in columns_combine_data_frame:
        print("\t" + column)
    print("------------------------------------------------------------------")
    print("Samples:")
    print(list_samples)
    print()
    columns = str(input("Write the name of the columns that do you want to observe separated with ',', without blank spaces: "))
    sample_names = str(input("Write the name of the samples do you want to observe separated with ',', without blank spaces (write none if you do not want to observe any sample, or all if you want to observe all): "))
    #Now we test if the name of the columns and the samples are correct
    # Create a variable for not execute the selection if the columns are incorrect
    not_work = 0
    selected_columns = []
    columns = columns.split(",")
    sample_names = sample_names.split(",")
    for col in columns:
        if col not in columns_combine_data_frame:
            print(col + " is not contained in the dataframe")
            break
        else: selected_columns.append(col)
    # if the name of the columns are not correct break
        
    for sample in sample_names:
        if sample == "none":
            break
        if sample == "all":
            for element in list_samples:
                selected_columns.append(element)
            break
        if sample not in list_samples:
            print(sample + " is not contained in the dataframe")
            break
        else: 
            selected_columns.append(sample)
    data_selected = combine_data_frame[selected_columns]
    print(data_selected)
    return data_selected
    

#Create a function for eliminatinf NAs from all the dataset or only from a selected column from the selected dataset

def eliminate_NAs(information_data_frame):
    columns_df = information_data_frame.columns
    print()
    print("Columns from our selected dataframe")
    print()
    for element in columns_df:
        print("\t" + element)
    print("If you want to eliminate all the NAs from the data write 'all'", 
    "If you only want to eliminate NAs from one or more columns write the columns with commas and without blank spaces", 
    "If you do not want to eliminate any NAs write 'none'", sep = "\n")
    eliminate_option = str(input("Choose your option: ")).upper()
    # If the user wants to eliminate all the NAs, 
    if eliminate_option == "ALL":
        information_data_frame = information_data_frame.dropna() 
    #If the user do not want to eliminate the NAs, we continue
    if eliminate_option == "NONE":
        print("The NAs have not been removed")
     #If the user want to eliminate data that have NAs from some columns
    else:
        #First of all we assure that the columns indicated exist in our dataframe
        list_columns_for_filtering = eliminate_option.split(",")
        cont_exist = 0
        for element in list_columns_for_filtering:
            if element not in columns_df:
                print("The column " + element + " does not exist in the selected dataframe")
                cont_exist += 1
        #If all the columns exists, eliminate the NAs
        if cont_exist == 0:
            information_data_frame = information_data_frame.dropna(subset=list_columns_for_filtering)
    return information_data_frame

# Create a function for group by some columns the selected dataframe

def groupby_selected_dataframe(information_data_frame):
    #Fist we provide again the columns to the user
    columns_combine_data_frame = information_data_frame.columns
    print("Columns from the dataframe selected")
    print()
    for element in columns_combine_data_frame:
        print("\t" + element)
    print()
    #Indicate one mandatory column for grouping
    mand_column = input(str("Select the 'ID' for the SNP  or the 'TRANSCRIPT_ID' as a reference: ")).upper()
    if mand_column != "ID" and mand_column != "TRANSCRIPT_ID":
        print("The mandatory column is not correct, try again")
    else:
        #Now ask for the column or columns they want to filter
        selected_for = str(input("Write one or two that you want to explore separated by ',' without black spaces: "))
        #create a selected_dataframe
        selected_for = selected_for.split(",")
        cont_exist = 0
        for element in selected_for:
            if element not in columns_combine_data_frame:
                print("The column " + element + " does not exist in the selected dataframe")
                cont_exist += 1
        #If all the columns exists, eliminate the NAs
        if cont_exist == 0:
            selected = selected_for
            selected.append(mand_column)
            select_for_viewing_df = information_data_frame[selected]
            #Now do the grouping
            data_grouping = select_for_viewing_df.groupby(selected_for).size().reset_index(name="COUNT")
            # Ask what data do you want in columns and what in rows
            column_variable = str(input("Write the column that do you want to be in the column: "))
            #Check again
            if column_variable in selected_for:
                index_variable = [valor for valor in selected_for if valor != column_variable]
                data_grouping_pivot = data_grouping.pivot(index=index_variable, columns=column_variable, values="COUNT")
                print(data_grouping_pivot)
                save_or_not = input("Do you want to save this information? Indicate 'yes' or 'no': ").strip().lower()
                if save_or_not == "yes":
                    file_name = input("Enter the file name to save the information (with extension, e.g., 'output.csv'): ")
                    data_grouping_pivot.to_csv(file_name)

            else:
                print("Something is wrong. Try again")
                
        # If there are columns that do not exist, give a error message
        else: print("Some names of the columns are wrong, try again")

#Function for change the code of the samples 
def convert_genotypes(vcf_eff_df,list_samples):
    for col in vcf_eff_df.columns:
        if col in list_samples:
            vcf_eff_df[col] = vcf_eff_df.apply(lambda row: "".join([row["REF"] if g == "0" else (row["ALT"] if g == "1" else ".") for g in row[col].split("/")]), axis=1)
    return vcf_eff_df

# Function for adding the prudu id to the dataset

def create_prudu_id_column_for_the_dataset(combine_data_frame):
    # Create the new variable and asign the value
    combine_data_frame["Prudu"] = combine_data_frame["TRANSCRIPT_ID"].apply(lambda x: x.split("|")[2] if isinstance(x,str) else None)
    return combine_data_frame

################# MAIN PROGRAM ################################3


def main(snp_effect_file, gft_genome, list_samples):
    # First of all, we upload the neccessary files and give information
    #Create the sampple list
    # Obtend the list of headers
    list_headers = select_columns_and_samples(snp_effect_file)
    print()
    print("The vcf file is loaded")
    print()
    vcf_eff_df = create_database(snp_effect_file, list_headers)
    #Change the code of the genotypes for the SNPs
    print("The process of converting the genotypes from 0 and 1 to alleles will take some minutes")
    vcf_eff_df= convert_genotypes(vcf_eff_df, list_samples)
    #Add a column identifying if the variant is a SNP or an INDEL
    vcf_eff_df = identifying_variants(vcf_eff_df)
    #vcf_eff_df.to_csv("/home/frangl97/Escritorio/Data_frame.csv", index=False)
    print()
    print("The data frame of the variants and SNPEffect is created")
    #print(vcf_eff_df)
    print()
    gtf = load_gtf_file(gft_genome)
    print("The gft file is being loaded")
    gtf_new = create_information_gtf(gtf)
    print()
    #print(gtf_new)
    print("The dataframe about the annotations of the gtf file is created")

    # Iniciate the bucle
    cont = 0
    while cont != 5:
        print("=======================================================================")
        print()
        print("Choose the desired option:")
        print()
        print("\t" + "1. Show the number of samples and their names")
        print("\t" + "2. Statistics")
        print("\t" + "3. Combine and extract information from the dataframes")
        print("\t" + "4. Save the snpEffect dataframe as a csv file")
        print("\t" + "5. Switch off the program")
        print()
        action = int(input("Introduce the number of the desire option: "))
        print()
        print("=======================================================================")

        if action == 1:
            # Execute the function that shows the information of the samples
            show_information_about_samples(list_samples)
        if action == 2:
            action_2 = 0
            while action_2 != 4:
                print("Choose the statistics or information do you want to know:")
                print()
                print("\t" + "1. Number of SNPs and Indels Identificated")
                print("\t" + "2. Statistics of snpEffect results")
                print("\t" + "3. Go back")
                print()
                action_2 = int(input("Choose the option: "))
                print()
                # If the option 1 is chosen, show a message indicating the number and types of variants identificated
                if action_2 == 1:
                    print()
                    print("Choose the analysis that you want: ")
                    print("\t" + "1. Number of SNPs and Indels identificated in all samples")
                    print("\t" + "2. Number of SNPs and Indels identificated in all samples by chromosome")
                    print()
                    action_2_1 = int(input("Choose the analysis that you want: "))
                    if action_2_1 == 1:
                        calculate_type_variants(vcf_eff_df)
                    if action_2_1 == 2:
                        calculate_type_variants_by_chromosome(vcf_eff_df)
                # If the option 2 is chosen, show a message indicating the type of statistics and if they want to know these statistics from all the samples,
                # from all the samples by chromosome, for one sample, or for a selection of samples individually
                if action_2 == 2:
                    print()
                    print("Choose the analysis that you want: ")
                    print("\t" + "1. Effects of the annotations for all the variants")
                    print("\t" + "2. Annotation impact for all the variants") 
                    print()
                    action_2_2 = int(input("Choose the analysis that you want: "))
                    # For obtaing information about the effects identificated by SNPEffect within all the variants 
                    if action_2_2 == 1:
                        calculate_type_annotation_variants(vcf_eff_df)
                    # For obtaining information about the effects identificated by SNPEffect per type of variants
                    if action_2_2 == 2:
                        calculate_annotation_impact(vcf_eff_df)
                
                if action_2 == 3:
                    break
                
                
                
        # If the option 3 is chosen, a message about the different options will appear, such as combind the information from gtf dataframe with the SNPEffect information, 
                #change the genotype of each sample to the combination of the REF/ALT alleles, 

        if action == 3:
            #First of all, the two dataframes are combined, and the "INFO" column is eliminated
            print()
            print("The dataframes will be combined")
            combine_data_frame = combine_dataframes_and_eliminate_info_columns(vcf_eff_df, gtf_new)
            information_data_frame = create_dataframe_for_viewing_information(combine_data_frame, list_samples)
            back_process = 0
            while back_process == 0:
                print()
                print("\t" + "1. Save the entire combine dataframe")
                print("\t" + "2. Select columns from the combine dataframe")
                print("\t" + "3. Eliminate NAs from all or a desire column from the entire combine dataframe")
                print("\t" + "4. Go back")
                print()
                action_3 = int(input("Choose the option do you want: "))

                #If the user wants to save the dataframe, we require him two things, the name of the file and the route
                if action_3 == 1:
                    file_name = input("Enter the file name to save the information (with extension, e.g., 'output.csv'): ")
                    information_data_frame.to_csv(file_name)
                #If the user wants to explore more the dataframe, some tools are provide, such as group by, compare samples if there are samples selected... also the option for eliminate NAs
                # is provided
                if action_3 == 2:
                    #First we ask for eliminate NAs
                    information_data_frame = eliminate_NAs(information_data_frame)
                    print(information_data_frame)
                    #After, we ask the user what we wants to do
                    print()
                    print("\t" + "1. Group by columns")
                    print("\t" + "2. Filter by value of a column")
                    print("\t" + "3. Go back")
                    print()
                    action_3_2 = int(input("Choose one option:"))
                    # If the user wants to group by one or more columns
                    if action_3_2 == 1:
                        groupby_selected_dataframe(information_data_frame)
                        # If the user wants to filter by value of a column
                    if action_3_2 == 2:
                        print()
                        print(information_data_frame.columns)
                        print()
                        variable = str(input("Write the column that do you want to filter: "))
                        #Print the unique values for the desire column
                        print(information_data_frame[variable].unique())
                        value_for_filtering = str(input("Write the value for filtering: "))
                        data_frame_filtered = information_data_frame[(information_data_frame[variable] == value_for_filtering)]
                        # Save the file
                        file_name = input("Enter the file name to save the information (with extension, e.g., 'output.csv'): ")
                        data_frame_filtered.to_csv(file_name)
                    else:
                        break                    


                #If the user wants only to eliminate NAs
                if action_3 == 3:
                    information_data_frame = eliminate_NAs(information_data_frame)
                    print(information_data_frame)
                    save_or_not = input("Do you want to save this information? Indicate 'yes' or 'no': ").strip().lower()
                    if save_or_not == "yes":
                        file_name = input("Enter the file name to save the information (with extension, e.g., 'output.csv'): ")
                        information_data_frame.to_csv(file_name)

                if action_3 == 4:
                    back_process = 1




        # If the option for saving the file as a csv is indicated, do that
        if action == 4:
            file_name = input("Enter the file name to save the information (with extension, e.g., 'output.csv'): ")
            vcf_eff_df.to_csv(file_name)


                        

                # It must be asked which will act as reference, the group or the sample
        if action == 5:
            print("The program will end")
            cont = 5
            






#Create the list of samples:
list_samples = []
with open(name_samples_file, 'r') as file_names:
    for line in file_names:
        line = line.strip()
        list_samples.append(line)

print(list_samples)

main(snp_effect_file, gft_genome, list_samples)
