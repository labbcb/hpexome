## hpexome: an automated workflow for processing whole-exome sequencing data

Whole-exome sequencing (WES) consists in the capture and sequence of all exons of protein coding genes in the human genome. It has been used in clinical applications like disease diagnosis. [Genome Analysis Toolkit (GATK)](https://software.broadinstitute.org/gatk/) is a programming framework for development of NGS data analysis tools. GATK provides tools for depth of coverage calculation analysis and single nucleotide polymorphism calling. These tools were used by large-scale sequencing projects such as 1000 Genomes Project and The Cancer Genome Atlas.

[GATK-Queue](http://gatkforums.broadinstitute.org/gatk/discussion/1306/overview-of-queue) is a command-line scripting framework for defining multi-stage genomic analysis pipelines combined with an execution manager that runs those pipelines from end-to-end. It has been used to implement workflows for WES data analysis. Our workflow for processing WES data was based on good practices described by authors of GATK software (Van der Auwera et al. 2002). It takes one or more BAM files, reference genome file, known SNP and indels databases, and outputs one or more VCF files.

Our workflow for processing WES data was based on good practices described by authors of GATK software. It takes one or more BAM files, reference genome file, known SNP and indels databases, and outputs one or more VCF files.

![Workflow for processing WES data](hpexome_workflow.png)

We used GATK-Queue workflow management to implement the workflow. It consists in a QScript file written in Scala programming language. The input parameters from every walker are set to default values but user can change theses values through the R interface.

## References

__From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline__ Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M, 2013 _CURRENT PROTOCOLS IN BIOINFORMATICS 43:11.10.1-11.10.33_. [Article](http://dx.doi.org/10.1002/0471250953.bi1110s43)
