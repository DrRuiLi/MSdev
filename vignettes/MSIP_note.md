---

title: "MSIP_note"
output: 
  html_document:
    toc: true
    toc_depth: 3
    toc_float: true
vignette: >
  %\VignetteIndexEntry{MSIP_note}
  %\VignetteEngine{knitr::rmarkdown}

##   %\VignetteEncoding{UTF-8}

# FUNCTION

## MSIP_get_isotopomer_data.targeted

### *Prompt*

```text
Write a function 'MSIP_get_isotopomer_data'
Input: 
    1. MSdev object
    2. mode, could be "untargeted" (default) or "targeted"
Output: 
    1. MSdev object with obj@advancedAna$MSIP$isotopomer_data
Untargeted Workflow: 
    1. Check if there are obj@advancedAna$MSIP$compound_table, if not, message and stop
    2. Check the format of the compound_table
    3. For cycle for each compound in compound_table:
        a. Calculate the mz of all possible isotopologue, and their adduct (M+H and M-H)
        b. Check if Spectra (in obj@spectra) exist for the isotopologue by matching mz and rt
        c. Separate the spectra by isotopologue and sample.source
        d. Construct cfmdata using get_CFM_data_from_smiles
        e. Construct msipcoredata using get_MSIPCoreData
        f. Solve by MSIPCore_solve
        g. Construct the list of MSIPcoredata as a matrix, row as isotopologue, such as ([13]C1, [13]C2…) and column as sample.source
    4. Package the result as a list of MSIPcoredata matrix
    5. Store the list in obj@advancedAna$MSIP$isotopomer_data
```

### *Review*

```text
Update1: 
    1. If creat tempdir for cfmd cache, save it into msdev obj
    2. iso_count_max default set to 30
    3. sp.ms2$sample.source  should be determined in  MSdev_extract_Spectra(), not match in this function
    4. Use sp.pol$isolationWindowTargetMz instead of precursorMZ, rtime instead of retentionTime
```

```text
Update2:
    1. In step 2, the iso_form are not integrated in all_matched
    2. In step 3 and the following step, do not use iso_counts, just iso_forms
```

```text
Update3:
    1. There no all_matched$iso_mz
```

```text
Update4:
    1. Do not use iso_labels (M0, M1, M2…), just use iso_form as label([13]C0, [13]C1, [13]C2… )
```

```text
Update5:
    1. You miss the polarity of the spectra, thus all spectra are input into pos or neg adduct 
    2. When message solved… also include adduct, and if there are no spectra, skip if there are no sp, and not message
```

```text
Update6: 
    1. The cycle for "Process each isotopologue and sample combination" should be:
        a. Iso_form
        b. Sample_source
        c. Polarity, After solve pos and neg, MSIPCore_merge() them
```

```text
Update7: 
    1. The sp.ms2 <- onDiskData_retrieve(object@spectra$MS2_Spectra) should be placed outside the loop, not repeatly I/O
```

```text
Update8: 
    1. The default cache dir should be get_dir_expand_from_onedrive("/Code/R/data/MSDB/CompoundDB/CFM_predicted_kegg.compdb_cfmd") when CompoundDB_path miss
    2.  the above procedure should be placed under a ifelse(mode == "targeted")

```

```text
# Thu Apr 23 15:10:36 2026 ------------------------------
update MSIP data structure msdev.aa@advancedAna$MSIP$isotopomer_data:
1. remove 'metadata', these param: ppm, rt_tol, thresh, et al. should be stored in the MSIPCore's 'solve' slot
2. remove 'CompoundInfo' 
3. the new structure should be :
  - msdev
    - advancedAna
      - MSIP
        - compound_table: a data.frame of compound info
        - isotopomer_data: a list of compound
          - MSIPCoreData_matrix: named as compound_id
          - MSIPCoreData_matrix: named as compound_id
          ......
          - MSIPCoreData_matrix: named as compound_id
```

```text
# Thu Apr 23 15:22:57 2026 ------------------------------
when process isotopomer_data, the 'iso_count_max' of 'MSIPCoreData' should be determined by the 'iso_form' and 'iso_ele', such as iso_form == "[13]C1', the iso_count_max should be 1
```

```text
# Thu Apr 23 16:00:21 2026 ------------------------------
1. in 'get_MSIP_isotopomer_data.targeted', when '# Split spectra by isotopologue form and polarity and sample.source', the step '# Get spectra indices for this isotopologue' only select 'iso_mz_val' of the first one (usually be positive) and thus miss the negative spectra, fix it
2. in 'heatmap_MSIPFragmentMap', the column name of M0, M1... should be update to M+0, M+1 ....
```

```text
# Thu Apr 23 16:04:35 2026 ------------------------------

1. the function 'MSIPCore_merge' are not updated for the current class and functions, update it to adapt current data structure, the suffix is no need with '_0' and '_1' already added in MSIPATOMMAP
```

```text
# Thu Apr 23 16:14:46 2026 ------------------------------

1. in 'get_MSIP_isotopomer_data.targeted',  when extract spectra, first match mz before match rt, if there are spectras with mz matched, but rt miss, message: "n sp targeted to xxx with rt range: xx - xx, selected m sp with rt range: xx - xx"
2. the step of 'MSIPCore_solve' should move after 'MSIPCore_merge'
```

```text
# Thu Apr 23 17:30:58 2026 ------------------------------

1. 'heatmap_MSIPFragmentMap':
  - do not cluster clomun, arrange with isotopomer set
2. 'plot_MSIPCore_spectra_consistency_hm': 
  - remove black border in rect_gp 
  - arrange row by rt, and replace the 'Log int' row_anno to a point plot to indicate the log10(Spectra$totIonCurrent)
2. 'plot_Spectra_CE', 'combineSpectra_groupby_ce' before process data and plot
```

```text
# Thu Apr 23 18:18:14 2026 ------------------------------

1. 'plot_MSIPCore_spectra_consistency_hm': 
  - add row annotation to show polarity
  - the row should be splited by both CE and polarity
2. there are codes to show the fragment coverage in "Script/Agent.R", but it was writen for the previous data structure, update it to adapt new data structure. write a function 'plot_MSIPCore_fragment_coverage', it input a MSIPCoreData and MSIPAtomMap, output the coverage figure
```

```text
# Thu Apr 23 18:47:17 2026 ------------------------------

1. 'plot_MSIPCore_spectra_consistency_hm': 
  - the Log10 TIC panel's x axis should be reversed: from right to left, min to max
```

#### Thu Apr 23 19:11:18 2026 ------------------------------

```text

1. add a param 'sp_top' in 'MSIP_get_isotopomer_data', to filter spectra including in MSIPCore
  - this need a helper function 'Spectra_filter_TIC':
    - input sp, topN (default 10), split_var (default: c("sample.source","CE","polarity"))
    - split sp by split_var
    - select topN TIC sp in every splited sp
    - return selected sp
  - use 'Spectra_filter_TIC' to filter sp_top sp and input to MSIPCore
2. in 'get_MSIP_isotopomer_data.targeted', 'Spectra_filter_TIC' should be placed after spliting spectra for every sample, polarity, isotopologue, and before 'get_MSIPCoreData', not applied to total spectra
```

#### Thu Apr 23 19:48:51 2026 ------------------------------

```text

the 'plot_MSIPCore_fragment_coverage' run error because there are gap between input msipcore and msipatommap: the msipcore is merged from both pos and neg, but msipatommap not, so we should make them align.
1. write function 'MSIPAtomMap_merge', this function merge MSIPAtomMap like 'MSIPCore_merge'
2. write function 'get_MSIPAtomMap_cached'
  - input compound_id, temp_dir (default same as 'get_MSIPAtomMap_from_smiles' )
  - retrieve cached MSIPAtomMap
  - MSIPAtomMap_merge
  - retrn merged MSIPAtomMap
3. update plot_MSIPCore_fragment_coverage to adapt both merged and un-merged MSIPAtomMap and MSIPCoreData, if anyone input is unmerged, message
```

#### Fri Apr 24 11:26:15 2026 ------------------------------

```text
in 'plot_MSIPCore_fragment_coverage', the color scale should be fixed:
ggplot2::scale_fill_manual(values = c("grey", "#FE9D71", "#CE4736", "#B80C27")) should be replaced with a named vec, from "1e0" to "1e6", if the value more than 1e6, colored with '1e6'
```

#### Tue Apr 28 14:19:35 2026 ------------------------------

```text
write a function 'MSIP_xcms_processing.targeted':
  - this function performance xcms analysis like 'MSdev_xcmsProcessing':
    - filter samples according to sample.info
    - find peaks 
    - rt adjust
    - group peaks 
    - fill features
  - the difference to 'MSdev_xcmsProcessing' include:
    - find peaks with param 'CentWaveParam(roiList )'
    - the 'roiList' should be calculated by the 'msdev@advancedAna$MSIP$compound_table':
      - simulate all possible isotopologue of each compound
      - calculate isotopologues' mz
      - construct a roiList according to the mz and rt
  - the return should be similar to 'MSdev_xcmsProcessing', stored feature's defination  
```

#### Tue Apr 28 16:33:52 2026 ------------------------------

```text
1. the part of roilist construction should be warpped in another function 'get_xcms_roi_list':
  - input a matrix with column: mz and rt, return a centwave aceeptable roilist
  - in given ppm and rt tolerance, calc mzmin, mzmax, rtmin, rtmax
  - match rtmin and rtmax to scan number
  - the scan number in different sample are different, pick the union range for each 
2. in 'MSIP_xcms_processing.targeted', just call the 'get_xcms_roi_list'
```

#### Wed Apr 29 12:39:36 2026 ------------------------------

```
I update the MSIP data structure, please read the part of '# At-a-glance relationship map' in 'vignettes/MSIP_Workflow.Rmd' to known the new data structure, these class and function need to update:
  - update class 'MSIPIsotopologueData'
  - write a new function 'get_MSIPIsotopologueData'
  - .Deprecated the old function 'MSIP_get_isotopologues_data'
before you go, you should read the old version of 'MSIPIsotopologueData', the new version should containing the same information but with updated data structure. you can refer to the new class and function about 'MSIPCoreData', some of information now move to MSIPCoreData or compound_table, these data should not included:
  - MSIP_get_isotopologues_data
  - get_MSIP_Isotopomer_SE
  - MSIP_get_isotopologues_data
```

#### Wed Apr 29 13:34:03 2026 ------------------------------

```
1. the msdev@advancedAna$MSIP$isotopologues_table should be removed, use rowData(MSIPIsotopologueData) instead
2. fix bug: Error in get_MSIPIsotopologueData(msdev.aa) : 
  isotopologues_table missing required columns: iso_count 
```

#### Wed Apr 29 13:45:16 2026 ------------------------------

```text
1. the 'get_MSIPIsotopologueData' should return a list of MSIPIsotopologueData, 'MSIP_get_isotopologues_data' return a MSdev object
2. update the doc "vignettes/MSIP_Workflow.Rmd" according to the updated data structure and functions
```

#### Thu Apr 30 20:04:14 2026 ------------------------------

```text

update 'MSIP_xcms_processing.targeted' to adpat 'MSIP_get_isotopologues_data':
1. add the 'MSIP_get_isotopologues_data' required var in xcms fdf after xcms processing according to the pre-constructed roilist
2. sometimes the M+0 feature is missing, manually add this feature, with rt using the corresponding isotopologues' rt and mz using the the theoritical mz +- 5 ppm, their featurevalue fill with 0
```

### *To Do*

## MSIP_get_isotopologues_data.untargeted

### *Prompt*

*Wed Apr 15 12:43:07 2026* Write a function `MSIP_get_isotopologues_data.untargeted` in new file

- input: MSdev object
- output: MSdev object with isotopomer data
- to do:
  1. check wheather there are isotopologue table

## get_MSIP_Isotopomer_SE

### *Prompt*

Fri Apr 24 13:19:06 2026 ------------------------------

```text
write a function 'get_MSIP_Isotopomer_SE' to extract quantification data from 'msdev@advancedAna$MSIP$isotopomer_data':
  - the return should be a SummarizedExperiment, column as sample, row as isotopomer
    - colData include: 
      - sample.source: from colname of msipcorematrix, also as column name
      - group: from sample.info
    - rowData include: 
      - isotopomer_id: paste0(compound_id,iso_form), also asrow name
      - compound_id
      - isotopomer_set: in which msip.core@Solve$MSIPIsotopomerMap@solve$isotopomer.set
      - isotopomer_total: count of isotopomer in the isotopomer_set
      - isotopologue_form: such as '[13]C1', '[13]C2'
  - for each samples' isotopomer data, extract from 'msip.core@Solve$MSIPIsotopomerMap@isotopomer.probability'
```

### *Review*

#### Fri Apr 24 13:37:09 2026 ------------------------------

```text
update 'get_MSIP_Isotopomer_SE':
1. in the result se's rowdata:
  - the 'isotopomer_id' should be like 'HMDB0000123_101', not 'HMDB0000123[13]C101'
  - add 'isotopologue_id', the value should be like 'HMDB0000123_[13]C1'
  - remove 'isotopomer_total', 'isotopomer_set'
  - add 'isotopomer_average_mix', calculated by:
    - for each sample, each isotopomer have a value 'mix', represent by the count of isotopomers in its isotopomer.set
    - for each isotopomer, calc the average of 'mix' across samples
```

#### Fri Apr 24 13:49:27 2026 ------------------------------

```text
in the result se's rowdata:
  - add 'compound_name'
  - add 'label.isotopomer', use compound name + isotopomer form, such as 'Glutamate_10000'
  - add 'label.isotopologue', use compound name + isotopologue form, such as 'Glutamate_M+1', 'Glutamate_M+2'
```

#### Fri Apr 24 14:06:30 2026 ------------------------------

```text
in the result se's rowdata:
  - the 'label.isotopomer' is fixed '001', i want it as the isotopomer's form like 'isotopomer_id'
  - the 'label.isotopologue' just extract number, "[13]C1" become 'M+131', but it should be "M+1" 
```

#### Fri Apr 24 20:20:50 2026 ------------------------------

```text
in the result se's rowdata:
  - I fix the problem where 'label.isotopomer' is fixed '001', the ifelse input length 1 vec and than output the first ele of meta$iso_names
  - the 'iso_name' should also be included in the rowdata, rename it to 'isotopome_form' and export to return
```

#### Thu Apr 30 19:49:28 2026 ------------------------------

```text
1. there are no 'object@advancedAna$MSIP$isotopologues_matrix$ratio_to_seed', do not rely on this data anymore
2. there should be 6 assay:
  - intensity.positive and intensity.negative, extract from featureValue (use 'get_xcms_quantify_MSIP')
  - ratio.positive and ratio.negative, from 'get_xcms_iso_fraction', read this function and make sure it called properly
  - purity.positive and purity.negative, from 'get_xcms_feature_purity_matrix'
```

```text

1. when '.get_iso_map_from_fdf', if there are no var need for isotopologue match, please message
2. the purity matrix is still all NA
```

#### Wed May 6 19:55:49 2026 ------------------------------

```text
'MSIP_xcms_processing.targeted' should update:
1. when fill the feautre of M+0, it's feature_id should be combination of 'FT' and number, such as 'FTF0001'
2. the iso_seed should be recorded for each isotopologue feature, indicate, what's the M+0 feature
3. the column in the compound_table should be append to the xcms.fdf, such as kegg.id, smiles (if exist).
```

### *To Do*

## get_MSIP_Isotopologue_SE

### *Review*

#### Thu May 7 16:05:25 2026 ------------------------------

```
write a function 'get_MSIP_Isotopologue_SE' to extract isotopologue data and construct a SummarizedExperiment object
1. read 'MSIP_xcms_processing.targeted' to understand the isotopologue annotation in xcms.fdf
2. read 'get_MSIP_Isotopomer_SE' to understand the data structure, the isotopologue se is similar to isotopomer se
```

```{text date="Fri May  8 14:23:24 2026"}
  the 'get_MSIP_Isotopologue_SE' return se with features as row, not compound's isotopologue, please merge the positve and negative ratio assay into one, name as 'ratio' and set as the default assay
```

```{text # Fri May  8 13:50:26 2026 ------------------------------}
1. keep the raito.positive and ratio.negative assay
2. fill the NA with 0
3. the rowData should include:
  - polarity: 0, 1, or 0;1 (both positive and negative)
  - similarity: the similarity of positive and negative ratio vector, using cosine similarity. all the isotopologue of the compound should be counted.
```

```{text date="Fri May  8 14:23:24 2026"}
1. the 'get_MSIP_Isotopologue_SE' should return a list of MSIPIsotopologueData, 'MSIP_get_isotopologues_data' return a MSdev object
2. update the doc "vignettes/MSIP_Workflow.Rmd" according to the updated data structure and functions

```

### *To Do*

## MSdev_find_isotope_label

### *Review*

```{text date="Fri May  8 14:23:24 2026"}

  rename 'MSdev_find_isotope_label' to 'MSIP_find_traced_isotopologue' and update:
  1. read the 'MSIP_get_isotopologues_data.targeted' to understand the labeled(traced) isotopologue identification
  2. be clear the isotopologue and traced isotopologue are different, the isotopologue is the original isotopologue (from naturally occurring isotopes), the traced isotopologue is the labeled(traced) isotopologue (from labeled isotopes tracer)
  3. there are two method to find the traced isotopologue:
    - method 1: compare the labeled and un-labeled sample ( the old version of 'MSdev_find_isotope_label')
    - method 2: compare the ratio of isotopologue to the theoretical ratio ( calculated by the MSCC::chemform_isotopes_pattern_enviPat)
  4. the actuall update should be in 'xcms_get_feature_isotope_label', rename it to "xcms_get_feature_traced_isotopologue", add a param 'method' to choose the method, default is 'method 1', implement the method 2
  5. make sure the output keep the same as the old version of 'MSdev_find_isotope_label'

```

```{text date="Fri May  8 14:23:24 2026"}
1. use "untraced_compare" and "natural_based" to name the two methods
2. when call 'untraced_compare', check the untraced sample and message use cli 
3. if any polarity is missing, just message use cli, not stop
```

### *To Do*

# Feature

## MSIPIsotopologueData

### plot_MSIPIsotopologueData_ratio


write a function 'plot_MSIPIsotopologueData_ratio' to plot the ratio of MSIPIsotopologueData
1. the input should be a MSIPIsotopologueData
2. the output should be a ggplot object
3. the ratio should be plot as a circlized bar plot
4. the bar should be colored by the isotopologue

Review:
1. plot should be split by sample.source, plot for each sample.source and then patchwork them together
2. the radius should be the ratio
3. the color should be the isotopologue
4. the legend should be the isotopologue
5. the title should be the compound name

Update:
1. add rt, intensity (recorded in se rowdata )  and formula (recorded in compound_table) after compound name in the title of the plot
  - the intensity should be shown as scientific notation, such as 1.23e+05
2. the legend remove the compound name, only show the isotopologue. such as 'M+1', 'M+2', 'M+3', not 'Glutamate_M+1', 'Glutamate_M+2', 'Glutamate_M+3'
3. place the compound name in title, others in subtitle
4. calculate the average of the ratio for each isotopologue, group the average ratio < min_ratio (add this arg to the function) into one group, and plot as a bar, they are 'other', colored with 'grey' and should be placed in the bottom of legend


### MSIP_report_isotopologue_ratio


write a function 'MSIP_report_isotopologue_ratio' to report the isotopologue ratio of MSIPIsotopologueData
1. the input should be a MSIP-msdev object
2. the output should be a report file, default as msdev@projectInfo$projectDir/report/MSIP_report_isotopologue_ratio.pdf
3. the report should include the isotopologue ratio for each sample.source
4. use a for loop to plot the isotopologue ratio for each isotopologue and then save the plot as a pdf file

Review:
1. add message_with_time to show the progress 
2. use export_graph2pdf to save the plot, not pdf()


### MSIP_get_isotopologues_data

1. remvoe '.fix_isotopologue_se_rownames' in 'MSIP_get_isotopologues_data', this step should be done in 'get_MSIPIsotopologueData'
2. refactor the code to make it more readable and maintainable:
  - do not use '.get_pol_mats' to wrap the three matrix calculation, directly call the function
  - show message for the three matrix calculation
3. the feature info should be named as: 
  - rt, average of positive and negative rt
  - rt.positive
  - rt.negative
  - intensity, average of positive and negative intensity
  - intensity.positive
  - intensity.negative
4. add a param 'purity' to control the purity calculation, default is TRUE, if FALSE, the purity matrix will be all NA
5. add rowdata for isotopologuedata:
  - rtmin, rtmax, rtsd: collect all isotopologue's rt, mz, and calculate the min, max, and sd
  - feature_id.positive and feature_id.negative: the feature_id of positive and negative 
6. add param 'ratio.aggregate' to control the ratio aggregation:
  - 'mean', mean of positive and negative ratio
  - 'wmean', weighted mean, the weight is the intensity of the polarity
  - 'max_intensity', default,  ratio using the ratio of the polarity with the highest intensity
7. the default assay of the isotopologue se should be 'ratio'
8. for every ratio assay, the value should be normalized to 1 across isotopologue, the value should be the ratio of the isotopologue to the sum of all isotopologue of the compound in the same sample.source





### MSIP_xcms_processing.targeted
1. Refactor the part of roiList construction, warpped in function and called in 'MSIP_xcms_processing.targeted': 
  - get_MSIP_xcms_roi_list
    - get_xcms_roi_list_from_compound_table
2. fill the default value for ion_mode, default is 0 (negative), 'get_MSIP_xcms_roi_list' do not need xcms.xcms as input, just use the msdev@xcmsData$ion_mode
3. refactor the part of '.annotate_fdf_with_iso', warpped in function 'MSIP_annotate_with_iso_grid' and called in 'MSIP_xcms_processing.targeted'
  - input: msdev object, iso_grid, iso_ele, mz_ppm, rt_tol, max_iso ...
  - output: msdev object with xcms's annotated fdf
4. remove param 'max_iso', when calculating the isotopologue, use the max element count of each compound
5. store the 'iso_grid' in msdev@advancedAna$MSIP$temp$iso_grid
6. remove the part of 'MSIP_annotate_with_iso_grid' from 'MSIP_xcms_processing.targeted', this step should be done in 'MSIP_annotate_with_iso_grid'
7. use featureValues(value = "maxo" ), not 'get_xcms_quantify_MSIP' to extract the intensity of the isotopologue features
8. the xcms chromPeaks sometimes containg too wide peaks (such as mz range of 0 ~ 200), add a function 'filter_xcms_chromPeaks_mz_width' to filter out the wide peaks (use ppm, default 20) before groupPeaks
9. it seems 'filter_xcms_chromPeaks_mz_width' remove the necessary peaks, doc it, and write another function 'fix_xcms_chromPeaks_mz_width' to replace the filter function in 'MSIP_xcms_processing.targeted':
  - this function find the peaks with mz width > ppm, and calculate the correct mzmin and mzmax, and then replace the peaks in the chromPeaks matrix
10. expose the param 'param' of 'findChromPeaks' to 'MSIP_xcms_processing.targeted':
  - default is NULL, if NULL, use get_MSdev_param to get the param and message the used param
  - if input a centwave param, message the used param and integrate the roilist to the param
11. refactor 'MSIP_xcms_processing.targeted', it now accept two types of input 'param':
  - a 'CentWaveParam' object, the current implementation
  - a 'CentWavePredIsoParam' object, the new implementation:
    - a iso_grid2 should be constructed by removing the M+n isotopologue, only keep the M+0 isotopologue from iso_grid, then inputted to 'findChromPeaks'
    - the 'maxIso' of 'CentWavePredIsoParam' should be the max value in the iso_grid
    - iso_grid2 only used for a temporary purpose, the original iso_grid still stored in the msdev object
12. move the M+0 inject process to 'MSIP_xcms_processing.targeted', after 'xcmsProcessingMS1', warpped in function 'xcms_inject_iso_seed'：
  - input: xcms object, iso_grid, iso_ele, mz_ppm
  - output: xcms object with M+0 feature added to the xcms.fdf
  - call this function after 'xcmsProcessingMS1'
13. add another fillChromPeaks step after 'xcms_inject_iso_seed' in  'MSIP_xcms_processing.targeted'


### MSIP_annotate_with_iso_grid
1. remove param 'iso_grid' and 'ion_mode'
2. read the msdev@advancedAna$MSIP$temp$iso_grid, if not exist, use the msdev@advancedAna$MSIP$compound_table to construct the iso_grid
3. loop for polarity, construct the fdf for each polarity
4. when annotate the iso_grid:
  - loop for each compound:
    - match compound's mz and rt to the iso_grid
    - group matched features by rt, with rt_tol.isotopologue (default 5)
    - for each group, calculate the sum of intensity of the features
    - select the feature with the highest intensity as the isotopologues of the compound
    - record the selected feature's rt.center, rt.sd
  - integrate all compounds' isotopologues features:
    - filter out conflicting features, which assgined to different compounds' isotopologues
    - message "XXX features are assgined to conflicting isotopologues: compound_id1, compound_id2, ..."
    - calcuate the rt diff of conflicted feature to the rt.center of the isotopologue, the feature with the smallest rt diff should be selected as the isotopologue of the compound
5. when inject the M+0 feature:
  - do not add peaks (chromPeaks) to the xcms object, just add the feature to the xcms.fdf (featureDefinitions and featureValues)
6. I checked the result, some of features are assigned to the wrong isotopologue: 
  - M+0 rt = 840, M+1 rt = 824, M + 3 rt = 857, they are absolutely wrong groupped (rt diff >  rt_tol.isotopologue)
7. the scoring for select feature as the isotopologue should be performed after all match for each polarity:
  - collect all matched features, both neg and pos
  - scoring use a block match:
    - all features (POS AND NEG) with mz matched should be grouped by rt, with rt_tol.isotopologue (default 5)
    - for each group, calculate the sum of intensity of the features
      - if M+n exist both positive and negative, the intensity should be the sum of positive and negative intensity    - the features group with the highest intensity should be selected as the isotopologues of the compound
    - the rda should update according to the selected features:
      - rt, rt.center, rt.sd, intensity, feature_id.positive, feature_id.negative, polarity, rtmin, rtmax
8. these var should be updated in the rda:
  - rt.min, rt.max (rename from rtmin, rtmax), rt.sd, rt.center: calculate by all of the selected isotopologue features, not by the pos and neg features of one isotopologue
  - remove rtsd, rt.sd is enough
9. rename param 'rt_tol' to 'rt_tol.reference', default is 60
10. add rda:
  - mz.positive and mz.negative: the mz of the positive and negative features
11. refactor the core loop of 'MSIP_annotate_with_iso_grid': (discard)
  - determine the rt.center for each compound:
    - loop for each compound:
      - confirm rt as the previous step 5:
        - match iso_grid and fdf
        - Cluster hits by RT
        - Pick the best RT cluster for this compound
        - Record rt_center and rt_sd for that compound’s chosen cluster
  - construct a fdf from iso_grid:
    - all isotopologue features should be added to the fdf, not only the detected features
    - mzmed as the theoretical mz, mzmin and mzmax as the theoretical mz +- 5 ppm
    - rtmed as the rt.center, rtmin/rtmax as the average of all detected features' rtmin/rtmax
  - replace the xcms.fdf with the constructed fdf
  - before fillChromPeaks, search chromPeaks by mzr/rtr against iso_grid to set peakidx; drop features with empty peakidx
  - fillChromPeaks
  - annotate the filled xcms.xcms object with iso_grid, use rt.center as the reference


### MSIP_export_isotopologue_acquisition_list
1. read the old function 'MSIP_export_isotopologue_acquisition_list' to understand the data structure of the acquisition list
2. read the 'MSIP_xcms_processing.targeted' and 'MSIP_get_isotopologues_data' to understand the isotopologue annotation in xcms.fdf
3. refactor the function to export the acquisition list for each isotopologue, the acquisition list should include the precursor m/z, collision energy, and retention time window
4. the acquisition list should be exported to the projectDir/acquisition_list/isotopologue_acquisition_list.csv








