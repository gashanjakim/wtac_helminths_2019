---
layout: paper
title: "Using prototyping to choose a bioinformatics workflow management system"
year: "2021"
shortref: "Jackson <i>et al.</i> <i>PLoS Comput. Biol</i> 2021"
nickname: prototyping-workflows
journal: "PLoS Computational Biology"
volume: 17
issue: 
pages: e1008622
authors: "Jackson M, Wallace EWJ, Kavoussanakis K"
image: /assets/images/papers/default-paper.svg
redirect_from: 
fulltext: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008622
pdflink: https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1008622&type=printable
github: https://github.com/riboviz/workflows
biorxiv: https://www.biorxiv.org/content/10.1101/2020.08.04.236208v1.full
pmid: 
pmcid: 
f1000: 
doi: "10.1371/journal.pcbi.1008622"
dryad_doi:
figshare_doi: 
altmetric_id: 87231880
category: paper
# Note: 'published' is a Jekyll keyword and does not refer to whether the paper is published, but rather to whether this Markdown should be part of the rendered site.
published: true
preprint: true
embargo: false	
peerreview: false
review: true
tags: [riboviz, bioinformatics, workflows, education]
---
{% include JB/setup %}

# Abstract 

Workflow management systems represent, manage, and execute multistep computational analyses and offer many benefits to bioinformaticians. They provide a common language for describing analysis workflows, contributing to reproducibility and to building libraries of reusable components. They can support both incremental build and re-entrancy—the ability to selectively re-execute parts of a workflow in the presence of additional inputs or changes in configuration and to resume execution from where a workflow previously stopped. Many workflow management systems enhance portability by supporting the use of containers, high-performance computing (HPC) systems, and clouds. Most importantly, workflow management systems allow bioinformaticians to delegate how their workflows are run to the workflow management system and its developers. This frees the bioinformaticians to focus on what these workflows should do, on their data analyses, and on their science.

RiboViz is a package to extract biological insight from ribosome profiling data to help advance understanding of protein synthesis. At the heart of RiboViz is an analysis workflow, implemented in a Python script. To conform to best practices for scientific computing which recommend the use of build tools to automate workflows and to reuse code instead of rewriting it, the authors reimplemented this workflow within a workflow management system. To select a workflow management system, a rapid survey of available systems was undertaken, and candidates were shortlisted: Snakemake, cwltool, Toil, and Nextflow. Each candidate was evaluated by quickly prototyping a subset of the RiboViz workflow, and Nextflow was chosen. The selection process took 10 person-days, a small cost for the assurance that Nextflow satisfied the authors’ requirements. The use of prototyping can offer a low-cost way of making a more informed selection of software to use within projects, rather than relying solely upon reviews and recommendations by others.

# Author summary 

Data analysis involves many steps, as data are wrangled, processed, and analysed using a succession of unrelated software packages. Running the right steps, in the right order, and putting the right outputs in the right places, is a major source of frustration. Workflow management systems require that each data analysis step be “wrapped” in a structured way, describing its inputs, parameters, and outputs. By writing these wrappers, the scientist can focus on the meaning of each step, and how they fit together, which is the interesting part. The system uses these wrappers to decide what steps to run and how to run these and takes charge of running the steps, including reporting on errors. This makes it much easier to repeatedly run the analysis and to run it transparently upon different computers. To select a workflow management system, we surveyed available tools and chose 4 in which we developed prototype implementations to evaluate their suitability for our project. We conclude that many similar multistep data analysis workflows can be rewritten in a workflow management system, and we advocate prototyping as a low-cost (both time and effort) way of making an informed selection of software for use within a research project.

