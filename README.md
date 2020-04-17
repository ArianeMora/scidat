# Sci-dat: Download Annotate TCGA

A package developed to enable the download an annotation of TCGA data from `https://portal.gdc.cancer.gov/`

## Install

```
pip install scidat
```

## Use
### API
The API combines the functions in Download and Annotation. It removes some of the ability to set specific directories etc but makes it easier to perform the functions.

See example notebook for how we get the following from the TCGA site:
```
    1. manifest_file
    2. gdc_client
    3. clinical_file
    4. sample_file
```

```
api = API(manifest_file, gdc_client, clinical_file, sample_file, requires_lst=None, clin_cols=None,
                 max_cnt=100, sciutil=None, split_manifest_dir='.', download_dir='.', meta_dir='.', sep='_')

```
Step 1. Download manifest data
```
# Downloads every file using default parameters in the manifest file
api.download_data_from_manifest()
# This will also unzip and copy the files all into one directory
```
Step 2. Annotation 
```
# Builds the annotation information
api.build_annotation()
```
Step 3. Download mutation data
```
# Downloads all the mutation data for all the cases in the clinical_file
api.download_mutation_data()
```
Step 4. Generate RNAseq dataframe
```
# Generates the RNA dataframe from the downloaded folder
api.build_rna_df()
```
Step 5. Get cases that have any mutations or specific mutations
```
# Returns a list of cases that have mutations (either in any gene if gene_list = None or in specific genes)
list_of_cases = api.get_cases_with_mutations(gene_list=None, id_type='symbol')

# Get genes with a small deletion
filter_col = 'ssm.consequence.0.transcript.gene.symbol'
genes = api.get_mutation_values_on_filter(filter_col, ['Small deletion'], 'ssm.mutation_subtype')

# Get genes with a specifc genomic change: ssm.genomic_dna_change
filter_col = 'case_id'
cases =  api.get_mutation_values_on_filter(filter_col, ['chr13:g.45340134A>G'], 'ssm.genomic_dna_change')

```
Step 6. Get cases with specific metadata information

Metadata list:
```
submitter_id
project_id
age_at_index
gender
race
vital_status
tumor_stage
normal_samples
tumor_samples
case_files
tumor_stage_num
example: {'gender': ['female'], 'tumor_stage_num': [1, 2]}
```
Method can be `any` i.e. it satisfies any of the conditions, or `all`, a case has to satisfy all the conditions in the meta_dict

```
# Returns cases that have the chosen metadata information e.g. gender, race, tumour_stage_num
cases_list = api.get_cases_with_meta(meta: dict, method="all")
```
Step 7. Get genes with mutations
```
# Returns a list of genes with mutations for specific cases
list_of_genes = api.get_genes_with_mutations(case_ids=None, id_type='symbol')
```
Step 8. Get values from the dataframe
```
# Returns the values, columns, dataframe of a subset of the RNAseq dataframe
values, columns, dataframe = get_values_from_df(df: pd.DataFrame, gene_id_column: str, case_ids=None, gene_ids=None,
                           column_name_includes=None, column_name_method="all")

```

### Download

```
# Downloads data using a manifest file
download = Download(manifest_file, split_manifest_dir, download_dir, gdc_client, max_cnt=100)
download.download()
```

```
# Downloads data from API to complement data from manifest file
# example datatype = mutation (this is the only one implemented for now)
download.download_data_using_api(case_ids: list, data_type: str)
```

### Annotate

** Generate annotation using clinical information from TCGA **
```
annotator = Annotate(output_dir: str, clinical_file: str, sample_file: str, manifest_file: str, file_types: list,
                 sep='_', clin_cols=None)
# Generate the annotate dataframe
annotator.build_annotation()

# Save the dataframe to a csv file
annotator.save_annotation(output_directory: str, filename: str)

# Save the clinical information to a csv file
annotator.save_annotated_clinical_df(output_directory: str, filename: str)

```

** Download mutation data for the cases of interest **
Note we first need to download the data using the `download_data_using_api` from above.
```
annotator.build_mutation_df(mutation_dir)

# Get that dataframe
mutation_df = annotator.get_mutation_df()

# Save the mutation dataframe to a csv
annotator.save_mutation_df(output_directory: str, filename: str)

```



