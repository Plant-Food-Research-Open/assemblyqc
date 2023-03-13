# Change Log

## Version 0.7 (14-Mar-2023)

1. Added Kraken2 and NCBI FCS Adaptor tools.
2. Added Assemblathon stats.
3. Added `Genometools gt stat` statistics for gff3 files.
4. Simplified pipeline flow chart.
5. Simplified conda environment.
6. TIDK process now uses a container instead of conda.

## Version 0.6.1 (8-Mar-2023)

1. Included results_dict and dependencies dict (without html formatting) to json.
2. Removed completed items in readme.
3. Fixed json dump repeating image url.

## Version 0.6 (17-Feb-2023)

1. Added LAI.
2. Now sorting sequences by size before feeding to TIDK.
3. Added skip switches for all the tools.
4. Added configuration annotations.
5. Optimised resource allocation.

## Version 0.5.1

1. Changed report parsers to allow alphanumeric ([a-zA-Z0-9_]) characters in the haplotype names.

## Version 0.5

1. Added TIDK

## Version 0.4

1. Added ability run BUSCO for multiple augustus species simultaneously
2. Formatted tabs into a drop down list for ease of navigation
3. Summary page has been added
4. BUSCO plots are now rendered on the summary page
5. Styling has been changed for better user experience

## Version 0.3

1. Added ability to run BUSCO for multiple haplotypes simultaneously
2. Updated README for new functionality
3. Adjusted styling for easier comparisons between reports
4. Incorporated conda instead of python venv

## Version 0.2

1. Added ability to run BUSCO for multiple lineages simultaneously
2. Removed intermediary outputDir
3. Standardised naming conventions across the tool
4. Updated README for new functionality
5. Change report.html layout to tab view
