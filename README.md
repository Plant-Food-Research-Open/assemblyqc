# Assembly QC Report Generator

Welcome to the Assembly QC report generator. This software is a Nextflow pipeline that can be used to perform BUSCO searches on fasta data and will generate an easy-to-read html report. This is a preliminary release, more capabilities will be added in the future.

<br>

# Installation

1. Copy the Github repository URL and run the following in your project folder:

```
$ git clone https://github.com/PlantandFoodResearch/assembly_qc.git
```

2. Create the Python virtual environment (venv) in the same directory:

```
$ python3 -m venv venv
```

3. Activate the Python virtual environment:

```
$ source venv/bin/activate
```

4. Download the Python dependencies:

```
$ pip install -r requirements.txt
```

5. Deactivate the virtual environment:

```
$ deactivate
```

<br>

# Getting sample data

In order to retrieve a sample of data to test the pipeline, SSH into PowerPlant and run the following bash script into your terminal:

```
$ seqkit sample -p 0.05 /workspace/hraczw/github/Karaka-genomics/09.hifiasm_phasing/karaka_phasing.hic.hap1.p_ctg.fa > {path_to_your_project}/input_data/test_sample.fa
```

You may adjust the size of the sample by changing the number from '0.05' to your desired size. For example, running the following will give you a sample of 10%:

```
$ seqkit sample -p 0.10 ...
```

<br>

# Usage

1. Load the required Nextflow module:

```
$ ml nextflow/22.10.4
```

2. Run the pipeline:

```
$ nextflow main.nf
```

The test data should take about 4-5 minutes to run. When the pipeline has finished running you will see the output of "Complete!" in the terminal along with the standard BUSCO output.

You will now see a results folder which will contain a file named 'report.html' and can be opened in Live Server.

---

Note: If you are using the sampled test data as shown in the previous section, it will be found in the input_data folder. If you are using your own data, please place it into the input_data folder. You will need to update the inputFilePath value in the nextflow.config file to match the path to the sample data. The nextflow.config file also contains the parameters for the BUSCO search which can be changed to suit the user.

---

After running the pipeline, if you wish to clean up the logs and work folder, you can run the following:

```
$ ./cleanNXF.sh
```
