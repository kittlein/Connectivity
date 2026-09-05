# Supplementary Material



## Finding where environment impacts recruitment: contributions of Population Genetics to Fisheries Biology



Marcelo J. Kittlein

This document provides supplementary material associated with the manuscript *Finding where environment impacts recruitment: contributions of Population Genetics to Fisheries Biology*.
## Genomic data

Genome-wide reference allele-frequency data for Atlantic herring
(*Clupea harengus*) were obtained from the Dryad Digital Repository.

The allele-frequency file `60.Neff.freq.gz` was downloaded from the dataset:

Pettersson, M. & Andersson, L. (2020). Reference allele frequencies for
populations pools of Atlantic Herring (*Clupea harengus*). Dryad.

DOI: [10.5061/dryad.pnvx0k6kr](https://doi.org/10.5061/dryad.pnvx0k6kr)

The dataset provides genome-wide allele frequencies and genomic positions
for geographically distributed population samples. These data were used as
the input for the genomic filtering and pair-wise $F_{ST}$ analyses described
below.

## Genomic data filtering and pair-wise FST estimation

Genome-wide allele-frequency data were filtered prior to the estimation of genetic differentiation among Atlantic herring samples. The analysis retains SNPs located on the 26 main chromosomes and excludes non-target population samples. To focus on broad-scale background differentiation, SNPs were filtered according to their allele-frequency variation among populations. Retained loci had a minor-allele-frequency proxy greater than 0.01, an among-population standard deviation in allele frequency between 0.02 and 0.08, and allele-frequency information for at least 90% of the populations.

To reduce the contribution of densely clustered genomic variants, the filtered dataset was further thinned by retaining one eligible SNP per 1-kb genomic window. Pair-wise genetic differentiation was then calculated among the 53 retained Atlantic herring samples.

### Pair-wise FST estimation

For each pair of populations, multilocus \(F_{ST}\) was estimated as the ratio of the summed difference between total and within-population expected heterozygosity to the summed total expected heterozygosity across loci:

$$
    F_{ST} = \frac{\sum(H_T-H_S)}{\sum H_T}.
$$

where $H_S$ is the mean expected heterozygosity within the two populations and $H_T$ is the expected heterozygosity calculated from their mean allele frequency. Only loci with allele-frequency information available for both populations were included in each pair-wise comparison.

### Output

The analysis generates:

- a complete pair-wise $F_{ST}$ matrix;
- a long-format table containing one observation per population pair;
- the number of SNPs contributing to each pair-wise estimate;
- a summary of the SNP-filtering procedure; and
- an RDS file containing the main intermediate and final analysis objects.

### Script

The complete genomic filtering and pair-wise $F_{ST}$ estimation workflow is available in:

[Filter genomic data and estimate pair-wise FST.r](Filter_genomic_data_and_estimate_pair-wise_FST.r)

## Isolation by sea distance

The relationship between geographic separation and genetic differentiation
was evaluated using the shortest distance by sea between population samples.
Pair-wise $F_{ST}$ estimates were obtained directly from the genomic
filtering and $F_{ST}$ estimation procedure described above, whereas marine
distances were obtained from the previously estimated shortest paths
connecting each pair of sampling localities.

Genetic differentiation was linearized as:

$$
\frac{F_{ST}}{1-F_{ST}}
$$

and related to the corresponding shortest distance by sea, expressed in
kilometers. Pair-wise genetic and geographic distance matrices were
constructed using the same population ordering, and isolation by distance
was evaluated with a Pearson Mantel test using 9,999 permutations.

The analysis generates a table containing the shortest sea distance and
genetic differentiation for each population pair, a table with the Mantel
statistic and permutation-based $P$ value, and a graphical representation
of the relationship between shortest sea distance and linearized $F_{ST}$.

### Script

The complete isolation-by-distance analysis is available in:

[Evaluate isolation by sea distance.r](Evaluate_isolation_by_sea_distance.r)

## Iterative Random Forest resistance-surface analysis

To identify oceanographic conditions potentially associated with genetic
differentiation, an iterative resistance-surface procedure was implemented
using the linearized pair-wise genetic differentiation estimates,
$F_{ST}/(1-F_{ST})$, obtained from the genomic analysis described above.

Oceanographic predictor layers were first assembled into a common raster
stack. Sampling localities falling outside valid raster cells were moved to
the nearest valid cell to ensure that all populations could be connected
through the marine domain.

### Initial shortest paths

The procedure was initialized with a spatially uniform resistance surface
whose value corresponded to the mean linearized $F_{ST}$ across all
population pairs. A transition matrix was constructed from this surface and
the shortest marine path between every pair of sampling localities was
calculated.

For each shortest path, the mean value of every oceanographic predictor
intersected by the route was extracted. These path-level environmental
values were then used as predictors of pair-wise linearized $F_{ST}$ in a
Random Forest model. The `nodesize` and `mtry` parameters were selected using
`tune.rfsrc`, and the fitted model was subsequently used to predict
linearized genetic differentiation across all valid raster cells.

### Iterative updating of the resistance surface

The spatial predictions from the Random Forest were used as the resistance
surface for the next iteration. Shortest paths were then recalculated on
this updated surface, oceanographic conditions were extracted along the new
routes, and a new Random Forest was fitted.

This sequence of:

1. calculating shortest paths;
2. extracting environmental conditions along each path;
3. fitting a Random Forest to linearized $F_{ST}$; and
4. updating the resistance surface from the Random Forest predictions

was repeated for 30 iterations.

For each iteration, the resistance surface, shortest paths, environmental
predictor values, in-sample $R^2$, and out-of-bag $R^2$ were retained. The
final predicted resistance surface was also exported as a GeoTIFF file.

### Script

The complete iterative shortest-path and Random Forest procedure is
available in:

[Iterative random forest resistance surface.r](Iterative_random_forest_resistance_surface.r)

## Stock-specific FST extremes

To identify spatial locations representing the lowest and highest predicted genetic differentiation within each herring management unit, the final predicted FST surface was intersected separately with the nine ICES stock polygons.

For each stock polygon, the raster was cropped and masked to the polygon extent, and all cells with valid predicted FST values were evaluated. The cell with the minimum predicted FST and the cell with the maximum predicted FST were then identified independently. Their geographic coordinates and corresponding raster-cell identifiers were retained for downstream extraction of stock-specific oceanographic conditions.

This procedure provides a reproducible, stock-specific spatial rule for selecting contrasting locations along the predicted differentiation surface, rather than relying on a common set of external reference points for all management units.

The implementation is available in:

[`extract_stock_fst_extremes.r`](extract_stock_fst_extremes.r)

The script generates tables containing, for each stock:

- stock identifier
- minimum or maximum predicted FST
- raster cell identifier
- longitude
- latitude

The selected coordinates can subsequently be used to extract temporal environmental information from oceanographic products for recruitment modelling.

## Annual oceanographic extremes at stock-specific FST cells

Annual oceanographic conditions were extracted at the two stock-specific
locations defined by the minimum and maximum predicted FST values within each
ICES management polygon.

Monthly oceanographic data for 1993–2025 were obtained from NetCDF files stored
in `ocean_monthly.7z`. The variables included chlorophyll-a concentration,
sea-surface temperature, sea-surface salinity, meridional current velocity, and
zonal current velocity.

For each stock, variable, year, and FST location, the annual minimum and maximum
values were calculated from the 12 monthly observations. This produced 20
environmental predictors per stock and year:

- five oceanographic variables
- two annual summaries: minimum and maximum
- two spatial locations: minimum-FST and maximum-FST cells

Predictor names encode these three components. For example,
`min_uo_maxFST` denotes the annual minimum zonal current velocity at the cell
with the maximum predicted FST within that stock polygon.

The NetCDF files are extracted from the compressed archive one at a time to a
temporary directory, processed, and removed before the next variable is read.
This avoids storing all uncompressed oceanographic files simultaneously.

The resulting Excel workbook contains one worksheet for each of the nine herring
stocks, with one row per year and the 20 annual oceanographic predictors. It
also includes supporting worksheets containing the extraction coordinates,
incomplete-year diagnostics, and the long-format annual data.

The implementation is available in:

[`extract_annual_oceanographic_extremes_at_fst_cells.r`](extract_annual_oceanographic_extremes_at_fst_cells.r)

## Monthly oceanographic predictors and recruitment modelling

Monthly chlorophyll-a, zonal and meridional surface-current velocity,
sea-surface salinity, and sea-surface temperature were extracted at the
minimum- or maximum-\(F_{ST}\) location selected for each herring stock.

For each stock, recruitment was related to environmental conditions using
moving monthly windows across a range of temporal lags. Recruitment age was
incorporated into the lag structure, and chlorophyll-a was allowed to occur
later than the physical oceanographic variables through an additional delay
parameter.

Candidate predictors were screened with Boruta, and model performance was evaluated using both
fitted and out-of-bag \(R^2\). The best temporal configuration was selected
according to out-of-bag performance.

Accumulated local effects (ALE) were calculated for the predictors retained in
the selected model. Spatial specificity was additionally evaluated by fitting
the selected model configuration at up to 1000 randomly sampled marine raster
cells within the corresponding stock polygon. This comparison keeps predictor
selection, temporal configuration, Random Forest hyperparameters, and
recruitment years fixed.

For each stock, the script saves an `RDS` object containing model settings,
selected predictors, fitted and OOB predictions, model-search results, ALE
objects, spatial-randomization results, and information required to reproduce
the associated figures without repeating the complete analysis.

Associated script:
[`monthly_ocean_and_recruitmen.r`](monthly_ocean_and_recruitment.r)
