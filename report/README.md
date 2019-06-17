# RNAseq Report
Silvia Morini  

The directory _report_ should contain:

1) **/references.bibtex** : I suggest that this file is uploaded to the repository, and every time someone adds a reference,this changing is uploaded in the repo - by doing so, we save ourselves time since we are not going to have to go and look for ourselves again.

2) **./logo.png, ./corp-styles.css** : Provides a function "DEXSeqDataSetFromFeatureCounts", to load the output of featureCounts as a dexSeq dataset (dxd) object.

3) **./summary.tsv** : this file is downloadable (https://portal.qbic.uni-tuebingen.de/portal/web/qbic/user-d -> Projects). Until now, it contains information about (see example in parenthesis):
    - PI affiliation (WG Ullrich)
    - PI address
    - Contact person
    - Contact affiliation and address (could also be removed: they are the same as PI affiliation and address)
    - PM
    - Manager affiliation (QBiC): could also be removed, the PM will always be from QBiC
    - Manager address (Auf der Morgenstelle … ): could also be removed, since the address will not often change, it can be written fixed in the template.

It would be great to also have directly parseable information regarding:
    - **PM and PI email addresses**
    - **university/institution/department**
    - **Description**: this could be directly parsed from the description section in the qPortal: if the description is not written well enough and/or we need to add something, we are free to. NOTE: which is the description that gets parsed? The one at the very top, under the project name, or "Detailed description" at the bottom?
 
4) **./results/pipeline_info/software_versions.csv** : this file, that is now output by the dev version, and will be officially present in the next release of the nextflow/rna-seq pipeline, is parsed in order to retrieve the information about the versions of the software used by the pipeline.

5) same for the folder **./results/MultiQC/multiqc_plots/svg(or /png or /pdf)** : in this way, we don’t have to manually extrapolate pictures from the multiqc report. present in the dev version, will be there from the next release (1.4) on.

6) **.results/MultiQC/multiqc_data/multiqc_general_stats.txt** : pase this file for the initial statistics, in order to avoid having to do a screenshot of the initial table in the report. Column currently displayed (can be easily personalised): "Sample Name", "% Assigned", "Assigned", "% Aligned", "% Trimmed", “% Dups", "% GC", "Seqs" + there will be a "Secondary Name" columns containing the customer's sample names, from the metadata sheet.

7) **./DESeq2/*** : all of the subfolders, either for figures/tables needed in the report or for those that end up in the report.zip folder in the end.

8) **./results/MultiQC/multiqc_report.html**

9) **./results/fastqc_reports.zip**: do not forget to zip this folder.
