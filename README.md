# Assembly QC Report Generator

Welcome to the Assembly QC report generator. This software is a Nextflow pipeline that can be used to perform BUSCO searches on fasta data and will generate an easy-to-read html report. This is a preliminary release, more capabilities will be added in the future.

<br>

# Installation

1. Copy the Github repository URL and run the following in your target folder:

```
$ git clone https://github.com/PlantandFoodResearch/assembly_qc.git
```

2. Navigate into the project

```
$ cd assembly_qc/
```

3. Create input_data folder

```
$ mkdir input_data
```

4. Create the Python virtual environment (venv):

```
$ python3 -m venv venv
```

5. Activate the Python virtual environment:

```
$ source venv/bin/activate
```

6. Download the Python dependencies:

```
$ pip install -r requirements.txt
```

7. Deactivate the virtual environment:

```
$ deactivate
```

<br>

# Getting sample data

In order to retrieve dummy data to test the pipeline with, run the following:

```
$ cp /output/genomic/fairGenomes/Fungus/Neonectria/ditissima/sex_na/1x/assembly_rs324p/v1/Nd324_canupilon_all.sorted.renamed.fasta \
./input_data/test_data.fasta
```

The data will now be in the input_data folder. It will be named test_data.fasta.

<br>

# Running the Pipeline

1. Load the required Nextflow module:

```
$ ml nextflow/22.10.4
```

2. Run the pipeline:

```
$ nextflow main.nf
```

The test data will take around 15 minutes to run. When the pipeline has finished running you will see the output of "Complete!" in the terminal along with the standard BUSCO output.

You will now see a results folder which will contain a file named 'report.html' and can be opened in your browser.

---

Note: If you are using your own data, please place it into the input_data folder. You will need to update the 'inputFilePath' value in the nextflow.config file to match the path to your data. The nextflow.config file also contains the parameters for the BUSCO search which can be changed to suit the user.

---

After running the pipeline, if you wish to clean up the logs and work folder, you can run the following:

```
$ ./cleanNXF.sh
```

<br>

# Final notes

This tool is designed to make your life easier. If you have any suggestions for improvements please feel free to contact me to discuss!
