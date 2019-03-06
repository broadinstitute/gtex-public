"""Subsets gene expression data by tissue.

This script allows the user to subset GTEx gene expression data by tissue.  It supports subsetting either gene read
counts or gene TPM.

Prerequisites:

    The user must create the following directories and download the following files into those directories
    from the https://gtexportal.org:

        annotations/
            GTEx_v7_Annotations_SampleAttributesDS.txt
        rna_seq/
            GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct
            GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct

Example:

    The user must specify the tissue to subset.  By default, the script will subset gene TPM values, but gene_reads
    can also be specified.

    For example, the following creates a tab-delimited file containing gene TPM values for the tissue
    Adipose_Subcutaneous:

        $ python gene_expression_by_tissue.py -t Adipose_Subcutaneous


    The following creates a tab-delimited file containing gene read counts for the tissue Whole_Blood:

        $ python gene_expression_by_tissue.py -t Whole_Blood -e gene_reads

Output:

    This script creates an output file in a directory named 'output'.  The output filename will be of the format

        <expression_type_id>_<tissue_id>.txt

    where:

        expression_type_id is one of: 'gene_reads' or 'gene_tpm"
        tissue_id is a lower-case version of a GTEx tissue ID.

    For example, the command:

        $ python gene_expression_by_tissue.py -t Adipose_Subcutaneous

    will result in an output file named:

        output/gene_tpm_adipose_subcutaneous.txt
"""

import pandas as pd
import argparse
import os
import errno

__author__ = 'Jared Nedzel'

tissue_hash = {
    'Adipose_Subcutaneous': 'Adipose - Subcutaneous',
    'Adipose_Visceral_Omentum': 'Adipose - Visceral (Omentum)',
    'Adrenal_Gland': 'Adrenal Gland',
    'Artery_Aorta': 'Artery - Aorta',
    'Artery_Coronary': 'Artery - Coronary',
    'Artery_Tibial': 'Artery - Tibial',
    'Bladder': 'Bladder',
    'Brain_Amygdala': 'Brain - Amygdala',
    'Brain_Anterior_cingulate_cortex_BA24': 'Brain - Anterior cingulate cortex (BA24)',
    'Brain_Caudate_basal_ganglia': 'Brain - Caudate (basal ganglia)',
    'Brain_Cerebellar_Hemisphere': 'Brain - Cerebellar Hemisphere',
    'Brain_Cerebellum': 'Brain - Cerebellum',
    'Brain_Cortex': 'Brain - Cortex',
    'Brain_Frontal_Cortex_BA9': 'Brain - Frontal Cortex (BA9)',
    'Brain_Hippocampus': 'Brain - Hippocampus',
    'Brain_Hypothalamus': 'Brain - Hypothalamus',
    'Brain_Nucleus_accumbens_basal ganglia': 'Brain - Nucleus accumbens (basal ganglia)',
    'Brain_Putamen_basal_ganglia': 'Brain - Putamen (basal ganglia)',
    'Brain_Spinal_cord_cervical_c-1)': 'Brain - Spinal cord (cervical c-1)',
    'Brain_Substantia_nigra': 'Brain - Substantia nigra',
    'Breast_Mammary_Tissue': 'Breast - Mammary Tissue',
    'Cells_EBV-transformed_lymphocytes': 'Cells - EBV-transformed lymphocytes',
    'Cells_Transformed_fibroblasts': 'Cells - Transformed fibroblasts',
    'Cervix_Ectocervix': 'Cervix - Ectocervix',
    'Cervix_Endocervix': 'Cervix - Endocervix',
    'Colon_Sigmoid': 'Colon - Sigmoid',
    'Colon_Transverse': 'Colon - Transverse',
    'Esophagus_Gastroesophageal_Junction': 'Esophagus - Gastroesophageal Junction',
    'Esophagus_Mucosa': 'Esophagus - Mucosa',
    'Esophagus_Muscularis': 'Esophagus - Muscularis',
    'Fallopian_Tube': 'Fallopian Tube',
    'Heart_Atrial_Appendage': 'Heart - Atrial Appendage',
    'Heart_Left_Ventricle': 'Heart - Left Ventricle',
    'Kidney_Cortex': 'Kidney - Cortex',
    'Liver': 'Liver',
    'Lung': 'Lung',
    'Minor_Salivary_Gland': 'Minor Salivary Gland',
    'Muscle_Skeletal': 'Muscle - Skeletal',
    'Nerve_Tibial': 'Nerve - Tibial',
    'Ovary': 'Ovary',
    'Pancreas': 'Pancreas',
    'Pituitary': 'Pituitary',
    'Prostate': 'Prostate',
    'Skin_Not_Sun_Exposed_Suprapubic': 'Skin - Not Sun Exposed (Suprapubic)',
    'Skin_Sun_Exposed_Lower_leg': 'Skin - Sun Exposed (Lower leg)',
    'Small_Intestine_Terminal_Ileum': 'Small Intestine - Terminal Ileum',
    'Spleen': 'Spleen',
    'Stomach': 'Stomach',
    'Testis': 'Testis',
    'Thyroid': 'Thyroid',
    'Uterus': 'Uterus',
    'Vagina': 'Vagina',
    'Whole_Blood': 'Whole Blood'
}

RNA_SEQ_DIR = 'rna_seq'
ANNOTATIONS_DIR = 'annotations'
OUTPUT_DIR = 'output'


expression_type_hash = {
    'gene_reads': 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct',
    'gene_tpm': 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct'
}


def sample_ids_by_tissue(tissue_name):
    """
    Reads the sample annotation file, filters for the given tissue, and for the rna_seq samples in the
    analysis freeze, and returns the sample IDs that match.
    :param tissue_name: A GTEx tissue name
    :return: A list of sample IDs containing the filtered sample annotations.
    """
    # read sample attributes
    sample_annotations_file = ANNOTATIONS_DIR + os.sep + 'GTEx_v7_Annotations_SampleAttributesDS.txt'
    print('Starting to read sample annotations file: ' + sample_annotations_file)
    samples_df = pd.read_csv(sample_annotations_file, sep='\t')

    # filter samples to just tissue
    tissue_filtered_samples_df = samples_df.loc[samples_df['SMTSD'] == tissue_name]

    # filter to just the samples used in the analysis freeze
    rna_seq_df = tissue_filtered_samples_df.loc[tissue_filtered_samples_df['SMAFRZE'] == 'RNASEQ']

    # get the  sample ids
    sample_ids = rna_seq_df['SAMPID']

    print('Finished reading and filtering sample IDs')
    return sample_ids


def expression_by_tissue(tissue_id, expression_type_id):
    """
    Subsets the expression type by the given tissue.
    :param tissue_id: The GTEx tissue_id by which to subset the expression data.
    :param expression_type_id: The type of expression data to subset (either gene_reads or gene_tpm)
    :return: Filtered Pandas dataframe containing the expression data of the desired type.
    """
    print('Subsetting ' + expression_type_id + ' by tissue: ' + tissue_id)

    tissue_name = tissue_hash[tissue_id]
    # sample_ids = sample_annotations_by_tissue(tissue_name)
    sample_ids = sample_ids_by_tissue(tissue_name)

    expression_df = expression_by_samples(sample_ids, expression_type_id)
    print('Subset complete')
    return expression_df


def expression_by_samples(sample_ids, expression_type_id):
    """
    Subsets the expression type by the given tissue.
    :param sample_ids: A Pandas dataframe containing the sample annotations, filtered by the GTEx tissue_id
    of interest.
    :param expression_type_id: The type of expression data to subset (either gene_reads or gene_tpm).
    :return: Filtered Pandas dataframe containing the expression data of the desired type.
    """

    expression_filename = RNA_SEQ_DIR + os.sep + expression_type_hash[expression_type_id]
    # read gene reads file
    print('Starting to read ' + expression_type_id + ' file: ' + expression_filename)
    desired_columns = ['Name', 'Description']
    desired_columns.extend(sample_ids)
    expression_df = pd.read_csv(expression_filename, sep='\t', skiprows=2, usecols=desired_columns)
    print('Finished reading ' + expression_type_id + ' file')
    return expression_df


def write_subset_file(expression_type_id, tissue_id, expression_df):
    """
    Writes the expression data in expression_df to a tab-delimited file.  The file will be written to a subdirectory
    named 'output', which will be created if it doesn't exist.  The filename will be of the form:

        <expression_type_id>_<tissue_id>.txt

    where:

        expression_type_id is one of: 'gene_reads' or 'gene_tpm"
        tissue_id is a lower-case version of a GTEx tissue ID.

    :param expression_type_id: Either gene_reads or gene_tpm
    :param tissue_id: A valid GTEx tissue ID.
    :param expression_df: A Pandas dataframe containing the expression data to write to a file
    :return: None
    """

    output_filename = OUTPUT_DIR + os.sep + expression_type_id + '_' + tissue_id.lower() + '.txt'

    print('Starting to write output file: ' + output_filename)
    create_output_dir()
    expression_df.to_csv(output_filename, sep='\t', index=False)
    print('Finished writing output file: ' + output_filename)


def create_output_dir():
    """
    Creates the directory named 'output' if it doesn't already exist.
    :return: None
    """
    try:
        os.makedirs(OUTPUT_DIR)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def parse_args():
    """
    Parses the command-line arguments.
    :return: A tuple containing two strings: the tissue_id (a legal GTEx tissue ID) and expression_type_id
    (either gene_reads or gene_tpm).
    """
    parser = argparse.ArgumentParser('Gene expression by tissue.')
    tissue_ids = tissue_hash.keys()
    expression_type_ids = expression_type_hash.keys()
    parser.add_argument('-t', '--tissue', help='Specify a tissue to subset gene expression data',
                        choices=tissue_ids, required=True)
    parser.add_argument('-e', '--expression_type', help='Specify the type of expression data to subset',
                        choices=expression_type_ids, default='gene_tpm')
    args = parser.parse_args()

    tissue_id = args.tissue
    expression_type_id = args.expression_type

    return tissue_id, expression_type_id


if __name__ == "__main__":
    tissue_identifier, expression_type_identifier = parse_args()
    expression_dataframe = expression_by_tissue(tissue_identifier, expression_type_identifier)
    write_subset_file(expression_type_identifier, tissue_identifier, expression_dataframe)
