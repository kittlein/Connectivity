# Supplementary Material

## Finding where environment impacts recruitment: contributions of Population Genetics to Fisheries Biology

Marcelo J. Kittlein

This document provides supplementary material associated with the manuscript *Finding where environment impacts recruitment: contributions of Population Genetics to Fisheries Biology*.

## Oceanographic covariates

Oceanographic covariates were derived from monthly surface fields obtained from the Copernicus Marine Environment Monitoring Service (CMEMS: `https://data.marine.copernicus.eu/viewer`). Physical variables, including sea surface temperature, sea surface salinity, and surface current velocity, were downloaded from the Global Ocean Physics Reanalysis product `GLOBAL_MULTIYEAR_PHY_001_030`, corresponding to the GLORYS12V1 global ocean reanalysis. This product provides eddy-resolving ocean physical fields at approximately 1/12 degree horizontal resolution and covers the period 1993-2023. Chlorophyll-*a* concentration was obtained from CMEMS ocean-colour data products. Because chlorophyll-*a* was provided at a different native spatial resolution, it was resampled to the spatial resolution of the physical variables before subsequent spatial processing and extraction of the environmental predictors used in the models.

## Neutral SNP filtering for connectivity analysis

SNP genotypes were obtained from the `herring` dataset distributed in the `FinePop2` package, which provides the data as a GENEPOP-formatted text object. Before estimating connectivity, loci previously reported in the literature as putative outliers were removed in order to represent gene flow using only neutral markers. This filtering step was intended to minimize the influence of loci potentially affected by divergent selection, because such loci can reflect adaptive population structure rather than contemporary or historical demographic exchange.

The filtering script parsed the GENEPOP header, standardized locus names by removing optional `Cha_` prefixes and parenthetical annotations, and matched them against a literature-based list of Atlantic herring outlier loci. Specifically, loci reported by Limborg et al. (2012) as global outliers in their Table 2 were excluded from the SNP panel before estimating connectivity. Using this procedure, 16 loci were removed from the original marker set, yielding a neutral dataset of 265 loci written to a GENEPOP-formatted file for estimating asymmetric connectivity.

Directional relative migration was estimated from the neutral GENEPOP file with the function `divMigrate()` in the `diveRsity` package. The analysis was run with 1000 bootstrap replicates and without filtering weak connections prior to estimation. The resulting matrix of relative migration values (`dRelMig`) was used as the connectivity matrix for subsequent visualization and spatial interpretation. Before plotting, diagonal values were set to zero, and edge colours were scaled according to connection strength to emphasize relative differences in inferred gene flow among sampled localities.

Associated script: [`filter_neutral_loci_and_estimate_connectivity.R`](filter_neutral_loci_and_estimate_connectivity.R).

## Connectivity surface modelling and polygon selection

To extend the locality-based connectivity estimates to the broader study region, a raster stack of seasonal environmental predictors was used as the spatial basis for modelling connectivity. This stack contained gridded oceanographic covariates for all marine cells in the study area, and only cells with complete environmental information were retained for model fitting and spatial prediction. Sampling localities were imported as spatial points and projected to the same coordinate reference system as the raster layers. The neutral-marker connectivity matrix derived from the relative migration analysis was used as the empirical representation of pairwise genetic connectivity among the 18 sampled sites.

As an initial spatial approximation, the mean connectivity value across all pairwise estimates was assigned to the raster domain to build a continuous transition surface. Using the `gdistance` package, this surface was converted into a transition object and corrected for geographic distance. Least-cost paths were then calculated between each locality and all remaining localities, thereby generating a spatial network of potential connectivity routes across the seascape. Mean values of the environmental covariates were extracted along these paths, and the extracted environmental summaries were linked to the observed pairwise connectivity values.

Connectivity was modelled as a function of the environmental covariates using Random Forest regression implemented with `randomForestSRC`. Model tuning was performed with `tune.rfsrc()` to optimize `nodesize` and `mtry`, and the fitted model was then used to predict connectivity for all raster cells with complete environmental information. This prediction step yielded a continuous regional surface of expected connectivity. Model performance was evaluated using both the apparent coefficient of determination ($R^2$) and the out-of-bag coefficient of determination ($R^2_{OOB}$).

The procedure was then iterated five additional times. In each iteration, the predicted connectivity surface from the previous step was used to recalculate least-cost paths, new environmental summaries were extracted along the updated paths, and the Random Forest model was refitted. This iterative approach allowed the connectivity surface and the inferred pathways to become progressively more consistent with one another. Among the fitted versions, the model iteration showing the highest out-of-bag $R^2$ was retained, and its predicted connectivity surface was used for final interpretation.

After selecting the best-performing iteration, the resulting connectivity surface was visualized together with the inferred least-cost paths and coastlines. Four focal polygons (`pol1b.kml` to `pol4b.kml`) were then overlaid on the regional connectivity surface. These polygons represented sectors of comparatively high predicted connectivity and were used as the spatial units for extracting seasonal environmental summaries. The extracted polygon-specific environmental time series constituted the covariates subsequently used in the Random Forest models fitted to the recruitment series of the different herring stocks.

Associated script: [`model_connectivity_surface.R`](model_connectivity_surface.R).

## Seasonal extraction of polygon-based environmental covariates for recruitment models

Monthly environmental values were extracted for each focal polygon from NetCDF (`.nc`) files spanning 1993-2023 and distributed in the data directory of this GitHub repository. These files contain the environmental-variable time series derived from Copernicus Marine products: `chla_monthly_1993-2023-herring.nc` for chlorophyll-*a*, `sst_sss_monthly_1993-2023-herring.nc` for sea surface temperature and sea surface salinity, and `vo_uo_monthly_1993-2023-herring.nc` for the zonal and meridional surface-current components. The NetCDF files were downloaded from Copernicus Marine for the spatial extent of the study region and contain monthly environmental time series covering 1993–2023.

Monthly environmental values were extracted for each focal polygon from three NetCDF files spanning 1993-2023: one containing chlorophyll-*a*, one containing sea surface temperature and sea surface salinity, and one containing the two surface current components. The four polygon layers were read from the KML files and standardized to a common geographic reference system before extraction. For each environmental variable, the corresponding raster time series was loaded as a multilayer object, and when necessary longitude coordinates were rotated from a 0-360 representation to a -180 to 180 domain to ensure spatial consistency with the polygon boundaries.

For every monthly raster layer, mean values were extracted within each polygon using area-weighted averaging so that partial cell overlap contributed proportionally to the polygon summary. The resulting monthly records were then assigned to year and month, and aggregated into two seasonal windows defined a priori for the recruitment analyses: spring (April-June) and fall (October-December). For each year, polygon, season, and environmental variable, the mean across the corresponding three monthly values was calculated. This procedure produced a polygon-specific seasonal time series for chlorophyll-*a*, sea surface temperature, sea surface salinity, and the zonal and meridional current components. The final output was organized in a wide table with one row per year and one column per variable-polygon-season combination, and this table was used as the predictor matrix for fitting the recruitment time-series models.

Associated script: [`extract_seasonal_polygon_covariates.R`](extract_seasonal_polygon_covariates.R).

## Recruitment time-series modelling

Annual recruitment series were modelled separately for each herring stock using Random Forest regression. For each stock, a modelling table was assembled by combining recruitment observations for 1993-2022, catch data lagged by one year (1992-2021), and the polygon-based seasonal environmental covariates described above. Environmental predictor names were standardized by abbreviating polygon identifiers (for example, `pol1` to `p1`) and seasonal labels (`spring` to `S` and `fall` to `F`) in order to simplify model handling and comparison among predictors.

For the main stock-specific fits, Random Forest models were estimated with the `rfsrc()` function in the `randomForestSRC` package. Annual recruitment was used as the response variable, whereas predictors included lagged catch together with all seasonal environmental covariates extracted for the four polygons. For each stock, `tune.rfsrc()` was first used to select `nodesize` and `mtry`, and those optimized values were then used in the final forest fitted with 2000 trees. Model performance was summarized with an apparent coefficient of determination (`R^2`) calculated from predicted versus observed recruitment values.

To evaluate the robustness of predictor effects, a repeated analysis of permutation importance was also performed for each stock. Using the stock-specific tuning parameters selected above, 200 replicate forests were fitted and permutation importance was computed for every predictor in every replicate. Importance values were summarized across replicates by their mean, standard deviation, median, 2.5% quantile, 97.5% quantile, and the probability of taking positive values. Mean positive importance values were additionally rescaled to relative percentages within each stock, thereby providing a standardized measure of the contribution of lagged catch and environmental covariates to recruitment variability.

Associated script: [`fit_recruitment_time_series.R`](fit_recruitment_time_series.R).

## Code annex

Curated versions of the core analysis scripts are provided as separate files in this GitHub repository. The scripts were simplified to retain only the steps required for the main analyses and to remove auxiliary plotting and post-processing routines on which the tables and figures presented in the main manuscript are derived.

- Neutral-locus filtering and directional connectivity estimation: [`filter_neutral_loci_and_estimate_connectivity.R`](filter_neutral_loci_and_estimate_connectivity.R)
- Iterative connectivity-surface modelling: [`model_connectivity_surface.R`](model_connectivity_surface.R)
- Seasonal extraction of polygon-based environmental covariates: [`extract_seasonal_polygon_covariates.R`](extract_seasonal_polygon_covariates.R)
- Recruitment time-series modelling: [`fit_recruitment_time_series.R`](fit_recruitment_time_series.R)

