# J1994---NCT04117087
R code for analysis of bulk TCRseq data for the colorectal cohort of NCT04117087. Instructions for running code are below.

A. Data structure
1. TCR sequencing files should be placed in /Data_Bulk within individual subfolders according to patient ID. File names should match the provided metadata.xlsx file
2. Code output will be in /Output folder with individualized patient output in subfolders according to patient ID.
3. Output of peripheral and tumor TCR-seq analysis will need to be placed in /Input/Bulk_Putative/. File names should match the provided metadata.xlsx file. Additional list of Published mKRAS TCRs are also provided. Additional instructions are provided below.

B. Peripheral TCRseq analysis
1. 'BulkTCR Peripheral Analysis' is a R script with code for annotating mKRAS-specific and treatment-enriched TCRs.
2. Code will need to be rerun for each patient based on patient ID.
3. After running this code, please place the file 'xx_Bulk_TCR_Summary' into /Input/Bulk_Putative/ for each patient
4. Also place the file 'xx_Bulk_Treatment_Enriched_Summary' into /Input/Bulk_Putative/Treatment-enriched for each patient

C. Merging mKRAS TCRs
1. 'BulkTCR Public TCR' is a R script with code for merging mKRAS TCRs identified across all patients as well as identifying Public mKRAS TCRs. This code requires that B.3 is completed
2. First time running the code, do not run section E which requires Public TCRs identified from the tumor

D. Tumor TCRseq analysis
1. 'BulkTCR Tumor Analysis' is a R script with code for determining mKRAS or treatment-enriched TCR tumor infiltration as well as their aggregate frequencies. This code requires that B.1-2 is completed
2. Code will need to be rerun for each patient based on patient ID.
3. After running this code, please place the file 'xx_Tumor_New_Public_TCR' into /Input/Bulk_Putative/Tumor for each patient
4. Also place the file 'xx_Tumor_TCR_Summary' into /Input/Bulk_Putative/Tumor-infiltration for each patient
5. Now, can run 'BulkTCR Public TCR' to update the full Public mKRAS TCRs list

E. Peripheral & tumor TCR frequency correlation
1. 'BulkTCR Spearmans' is a R script with code for determining whether mKRAS or treatment-enriched TCR frequencies in the tumor correlates with their frequencies in the periphery. This code requires that C.2 and D.3-5 is completed
