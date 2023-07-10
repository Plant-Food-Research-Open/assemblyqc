# Configuring the Pipeline for a Complete Run

This document explains how to correctly configure the pipeline for a complete run. For a minimal example which documents basic configuration along with pipeline execution, refer to [Quick Start: A Minimal Example](./minimal_example.md).

- [Configuring the Pipeline for a Complete Run](#configuring-the-pipeline-for-a-complete-run)
  - [Complete Example Configuration File](#complete-example-configuration-file)
  - [ASSEMBLATHON STATS](#assemblathon-stats)
  - [NCBI FCS Adaptor](#ncbi-fcs-adaptor)
  - [NCBI FCS GX](#ncbi-fcs-gx)
  - [BUSCO](#busco)
  - [TIDK](#tidk)
  - [LAI](#lai)
  - [KRAKEN2](#kraken2)
  - [HIC](#hic)
  - [SYNTENY](#synteny)

## Complete Example Configuration File

This document explains the pipeline configuration using an example configuration file packaged with the pipeline. Refer to 'conf/test_full.config'. This configuration is an expansion of the 'conf/test_minimal.config' covered in [Quick Start: A Minimal Example](./minimal_example.md).

## ASSEMBLATHON STATS

There is only one configurable parameter for this module: `n_limit`. This is the number of 'N's for the unknown gap size. This number is used to split the scaffolds into contigs to compute contig-related stats. NCBI's recommendation for unknown gap size is 100 <https://www.ncbi.nlm.nih.gov/genbank/>.

> ⚙️ From conf/test_full.config

```groovy
assemblathon_stats {
    n_limit = 100
}
```

## NCBI FCS Adaptor

This module has only one parameter: `empire`. The permissible values are: `euk` for Eukaryotes and `prok` for Prokaryotes.

> ⚙️ From conf/test_full.config

```groovy
ncbi_fcs_adaptor {
    empire = 'euk'
}
```

## NCBI FCS GX

Following parameters must be configured:

- `tax_id`: The taxonomy ID for all the target assemblies listed in the `target_assemblies` parameter. A taxonomy ID can be obtained by searching a *Genus species* at <https://www.ncbi.nlm.nih.gov/taxonomy>. A single ID for all assemblies implies that the pipeline is designed to be used for checking one or more assemblies of the same *species* in one run.
- `db_manifest_url`: This URL specifies the database version used by the pipeline.
- `db_path`: This is the path to the database files stored on a directory accessible to the pipeline. The data placed inside this directory should match with the `db_manifest_url`. Otherwise, the pipeline fails with an error. Before running the pipeline, the user must ensure that the database is correctly downloaded and placed in a directory accessible to the pipeline. Setup instructions are available at <https://github.com/ncbi/fcs/wiki/FCS-GX>. The database directory should contain following files:

```bash
all.assemblies.tsv
all.blast_div.tsv.gz
all.gxi
all.gxs
all.manifest
all.meta.jsonl
all.README.txt
all.seq_info.tsv.gz
all.taxa.tsv
```

> ⚙️ From conf/test_full.config

```groovy
ncbi_fcs_gx {
    tax_id          = "35717"
    db_manifest_url = "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/r2023-01-24/all.manifest"
    db_path         = "/workspace/ComparativeDataSources/NCBI/FCS/GX/r2023-01-24"
}
```

## BUSCO

Following parameters must be configured:

- `mode`: geno or genome, for genome assemblies (DNA), tran or transcriptome, for transcriptome assemblies (DNA); and prot or proteins, for annotated gene sets (protein).
- `lineage_datasets`: A list of BUSCO lineages. Any number of lineages can be specified. Each target assembly is assessed against each of the listed lineage. To select a lineage, refer to <https://busco.ezlab.org/list_of_lineages.html>
- `download_path`: A directory where the BUSCO can download and cache its databases. BUSCO manages download and validation of the databases itself, therefore, the user does not have to manually setup these databases.

> ⚙️ From conf/test_full.config

```groovy
busco {
    mode              = "geno"
    lineage_datasets  = ["fungi_odb10", "hypocreales_odb10"]
    download_path     = "/workspace/ComparativeDataSources/BUSCO/assembly_qc"
}
```

## TIDK

Following parameter must be configured:

- `repeat_seq`: The telomere search sequence. To select an appropriate sequence, see <http://telomerase.asu.edu/sequences_telomere.html>. Commonly used sequences are: TTTAGGG (Plant), TTAGGG (Fungus, Vertebrates), TTAGG (Insect).

The following parameters are optional:

- `filter_by_size`: Set this flag to 1 to filter out assembly sequences smaller than the size specified by the next parameter (default: 0).
- `filter_size_bp`: Minimum size of the assembly sequence processed by TIDK (default: 1000000 (1Mbp)).

> ⚙️ From conf/test_full.config

```groovy
tidk {
    repeat_seq = "TTAGGG"
}
```

In the example configuration above, the `filter_by_size` and `filter_size_bp` are not set. The pipeline will pick up their default values from 'nextflow.config' file.

## LAI

Following parameters must be configured:

- `mode`: "" for Standard, "-q" for Quick and "-qq" for Very Quick. The Very Quick mode produces Raw LAI which is only valid for intra-species comparisons.
- `pass_list`: A list of lists which specifies the `*pass.list` file for each assembly. This file is generated by EDTA or LTR_retriever software. If this parameter is `[]`, the pipeline first runs EDTA and itself generates the `*.pass.list` file for each assembly. If `*.pass.list` file for each assembly is available, specify it in the `pass_list` parameter along with the tag for the corresponding target assembly. Default value for this parameter is `[]`. Here is an example:

```groovy
target_assemblies   = [
    ["hap1", "/input/test_data/default/test_data1.fasta.gz"],
    ["hap2", "/input/test_data/default/test_data2.fasta"]
]

lai {
    pass_list       = [
        ["hap1", "/input/test_data/default/test_data1.pass.list"],
        ["hap2", "/input/test_data/default/test_data2.pass.list"]
    ]
}
```

Notice that the tags (hap1, hap2) in target_assemblies and pass_list are matched.

- `out_file`: Similar to `*.pass.list`, an `*.out` file is also required for LAI calculation. The specification rules for `out_file` are exactly same as those for `pass_list`. Default value for this parameter is `[]`. Here is an example:

```groovy
target_assemblies   = [
    ["hap1", "/input/test_data/default/test_data1.fasta.gz"],
    ["hap2", "/input/test_data/default/test_data2.fasta"]
]

lai {
    pass_list       = [
        ["hap1", "/input/test_data/default/test_data1.pass.list"],
        ["hap2", "/input/test_data/default/test_data2.pass.list"]
    ]

    out_file        = [
        ["hap1", "/input/test_data/default/test_data1.out"],
        ["hap2", "/input/test_data/default/test_data2.out"]
    ]
}
```

> ⚠ **WARNING**: EDTA modifies the sequence names in a fasta file if they are not short (<=13 characters) and simple (i.e, letters, numbers, and underscore; no comments allowed in the ID field). When specifying the `*pass.list` and `*.out` files, ensure that these files have the same sequence IDs as those in the fasta files specified by the `target_assemblies` parameter.

- `monoploid_seqs`: A list of lists which specifies the `-mono` parameter-file for LAI when processing a polyploid assembly. The `-mono` parameter-file is a single column text file listing IDs of the monoploid sequences for a polyploid assembly. If this parameter is not needed, it can be set to `[]`. If only some of the assemblies listed in `target_assemblies` are polyploid, the `-mono` parameter-file can be specified only for those assemblies. Similar to the `pass_list` parameter, an assembly is identified by its tag. Here are the contents of an example `-mono` parameter-file:

```TSV
CP031385.1
CP031386.1
CP031387.1
CP031388.1
CP031389.1
CP031390.1
CP031391.1
```

- `edta::is_sensitive`: "--sensitive" parameter for the EDTA software. Set to 1 to turn on Sensitive (very slow) and 0 to turn off Sensitive. Default is 0. For Tair10, the LAI score with Sensitive turned off is 18.06 and with Sensitive turned on is 18.09.

> ⚙️ From conf/test_full.config

```groovy
lai {
    mode                = "" // Standard
    
    pass_list           = []
    out_file            = []

    monoploid_seqs      = [
        ["FI1", "./docs/test_files/FI1.monoploid.seqs.txt"]
    ]

    edta {
        is_sensitive    = 0
    }
}
```

Notice that the default values are used for `pass_list` and `out_file`. This means that the pipeline will first run EDTA to perform the repeat annotation and then calculate the LAI score. Moreover, the pipeline will only consider the sequences listed in `monoploid_seqs` when calculating LAI.

## KRAKEN2

Following parameters must be configured:

- `db_url`: The URL for the KRAKEN2 database. To select a DB, see <https://benlangmead.github.io/aws-indexes/k2>. By default, the pipeline uses PlusPFP database.
- `download_path`: The directory where the pipeline can download and cache the data. The pipeline can download the database itself. However, care must be taken when switching databases. If the `db_url` is changed, a new `download_path` should be provided, otherwise, the pipeline will fail with error saying that the downloaded database is corrupted.

> ⚙️ From conf/test_full.config

```groovy
kraken2 {
    db_url        = "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20230314.tar.gz"
    download_path = "/workspace/ComparativeDataSources/kraken2db/k2_pluspfp_20230314"
}
```

## HIC

Following parameter must be configured:

- `paired_reads`: A relative or absolute path to paired reads in fastq.gz format, or a SRA ID. The format for file path is `*R{1,2}*.(fasta|fq).gz`. An example is '/input/genomic/fungal/Neonectria/Genome/20190506_CAGRF19591_CGYCF_HiC/PG_PETUNIA_HiC_CGYCF_CACTCA_L001_R{1,2}.fastq.gz'.

> ⚙️ From conf/test_full.config

```groovy
hic {
    paired_reads = "SRR8238190"
}
```

## SYNTENY

Following parameters must be configured:

- `assembly_seq_list`: This is a list of lists which specifies the `*.seq.list` file for each assembly declared in the `target_assemblies` parameter. A `*.seq.list` file is a two column tab-delimited txt file listing fasta sequence ids (first column) and labels for the synteny plots (second column). An example file is shown below:

```TSV
CP031385.1  FI1_1
CP031386.1  FI1_2
CP031387.1  FI1_3
CP031388.1  FI1_4
CP031389.1  FI1_5
CP031390.1  FI1_6
CP031391.1  FI1_7
```

This parameter is specified as a list of lists. Each sub-list has two elements. The first element identifies the assembly `target_assemblies` using the assembly tag. The second element is the path to the `*.seq.list` file. Here is an example:

```groovy
target_assemblies = [
    ["hap1", "/workspace/assembly_qc/test_data/default/test_data1.fasta.gz"],
    ["hap2", "/workspace/assembly_qc/test_data/default/test_data2.fasta"]
]

assembly_seq_list = [
    ["hap1", "/workspace/assembly_qc/test_data/default/test_data1.seq.list"],
    ["hap2", "/workspace/assembly_qc/test_data/default/test_data2.seq.list"]
]
```

- `xref_assemblies`: A list of lists which specifies reference assemblies against which the synteny should be performed. This parameter can be set to an empty list `[]` if reference assemblies are not available. To specify a reference assembly, three items must be declared. First a unique tag for the reference assembly, second a fasta file (fasta, fasta.gz) for the assembly, and, third, a `*.seq.list` file. Here is an example:

```groovy
xref_assemblies = [
    [
        "GenomeA",
        "/workspace/assembly_qc/test_data/default/test_data3.fasta",
        "/workspace/assembly_qc/test_data/default/test_data3.seq.list"
    ]
]
```

The following parameters are optional:

- `between_target_asm`: Set it to 1 to create syntenic plots between each pair of target_assemblies. Default is 1. This parameter is useful if multiple assemblies are specified by the `target_assemblies` parameter and the user needs control over whether syntenic plots are created between each pair of assemblies or not.
- `many_to_many_align`: Set it to 1 to include alignment blocks with many-to-many mappings or set to 0 to only include 1-to-1 mappings. Default is 0. See the documentation of `dnadiff` for further details: <https://github.com/mummer4/mummer/blob/master/docs/dnadiff.README>
- `max_gap`: Alignments within this distance are bundled together. Default: 1000000 (1 Mbp).
- `min_bundle_size`: After bundling, any bundle smaller than this size is filtered out. Default: 1000 (1 Kbp).
- `plot_1_vs_all`: Set it to 1 to create a separate synteny plot for each contig of the target assembly versus all contigs of the reference assembly. Set it to 0 to create a single plot for each target assembly against each reference assembly. This joint plot is also created when `plot_1_vs_all` is set to 1. Default: 0.
- `color_by_contig`: Set it to 1 to color the synteny plot by contig. Set it to 0 to color the synteny plot by the number of links in a bundle. Default: 1.

> ⚙️ From conf/test_full.config

```groovy
synteny {

    assembly_seq_list   = [
        ["FI1", "./docs/test_files/FI1.seq.list"]
    ]

    xref_assemblies     = [
        [
            "TT_2021a",
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/021/950/295/GCA_021950295.1_ASM2195029v1/GCA_021950295.1_ASM2195029v1_genomic.fna.gz",
            "./docs/test_files/TT_2021a.seq.list"
        ],
    ]
}
```
