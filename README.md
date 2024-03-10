# Rakeiora Primary Care / pharmacogenomic

[Rakeiora](http://rakeiora.ac.nz)

---
Primary Care public workflow aka pharmacogenomic

---

NeSI does not claim ownership or authorship of this workflow.
Authorship credits will follow later.

Based on Ben Halliday's work.

Adjusted to be used within Rakeiora by Matt Pestle.

The defined input for this workflow in Rakeiora should be variants,
multi-selected, merged globally into a file named whatever the value
of inputvcffile in the config.yaml nominates.

The bedfile used when running the workflow should match that
of the gene specified in the ALLELES variable. The values of this
variable should be uniform in the gene portion and the
appropriate/corresponding files must be present in the pharmvar data.

This bedfile can be determined at
[Rakeiora test - gene bedefile] (https://test.rakeiora.ac.nz/app/sandbox/gene.xsql)
Depending on your data sources, you may need to tick the "drop chr" option
(The UCSC reference require the chr, while the Ensemble ones don't).

You can peruse this data from your jupyter hub in the Rakeiora
test area -
[Rakeiora test] (https://test.rakeiora.ac.nz)

Rakeiora will retrieve the bedfile specified regions for the
sample selected and merge them using bcftools into a file
named as you specify.

DRUGS should be in the "drug(s)" column of pharmgkb file
specified in the config.yaml
(also perusable from the /shared area in the jupyter hub
of Rakeiora test area)

The GetDrugList Rule currently just generates random data.
It will be replaced with calls to ALEA to get real data
once all the necessary security and setup steps are done.
