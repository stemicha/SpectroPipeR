# SpectroPipeR 0.4.6

* feature: add directLFQ to HTML reporting

# SpectroPipeR 0.4.5

* bugfix: directLFQ column selection: use SpectroPipeR report columns now

# SpectroPipeR 0.4.4

* feature: integrate directLFQ for label-free quantification, using the directLFQ python package [MannLabs/directlfq](https://github.com/MannLabs/directlfq) 
  (Ammar, C., Schessner, J. P., Willems, S., Michaelis, A. C. & Mann, M. 
  "Accurate label-free quantification by directLFQ to compare unlimited numbers of proteomes", 2023, doi:10.1101/2023.02.17.528962), 
  adapted via the code of [RdirectLFQ](https://github.com/mildpiggy/Rdirectlfq) package (author: Zhenhuan Feng) with modifications.
  
# SpectroPipeR 0.4.3

* feature: add automatic missing value imputation by halve all values within the lowest 1% quantile distribution of ion quantities to create a new dataset used for imputing missing values.

# SpectroPipeR 0.4.2

* FIX: for duplicated raw file name in path

# SpectroPipeR 0.4.1

* add Spectronaut colors to XIC (MS1 and MS2) plot (optional: Spectronaut_colors = TRUE)

# SpectroPipeR 0.4.0

* add graphical user interface (GUI) function
* add identified and not-identified ion intensity plot
* updated HTML reports
* detailed description of condition-wise filtering & removing of oxidized methionine ions
* updated manual

# SpectroPipeR 0.3.0

* Initial CRAN submission.
