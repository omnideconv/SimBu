pkg.globals <- new.env()
pkg.globals$sfaira_streamline_feature_version <- "104"

pkg.globals$epic_scaling <- c(
  0.4016,
  1.4196,
  1.4196,
  0.1300,
  0.4396,
  0.3952,
  0.3952,
  0.3952,
  0.3952,
  0.3952,
  0.4000,
  0.4000
)
names(pkg.globals$epic_scaling) <- c(
  "B cells",
  "Macrophages",
  "Monocytes",
  "Neutrophils",
  "NK cells",
  "T cells",
  "T cells CD4",
  "T cells CD8",
  "T helper cells",
  "T regulatory cells",
  "otherCells",
  "default"
)

pkg.globals$abis_scaling <- c(
  22824.322779,
  21456.905980,
  17827.199342,
  14262.064545,
  10660.952413,
  24162.648934,
  15682.087676,
  13365.625992,
  28128.898864,
  9546.742675,
  25008.840273,
  8610.360645,
  22344.663722,
  20837.569855,
  57322.183941,
  28488.194770,
  325800.989149
)
names(pkg.globals$abis_scaling) <- c(
  "Monocytes",
  "NK cells",
  "T cells CD8 central memory",
  "T cells CD4",
  "T cells CD8",
  "B Naive",
  "T cells CD4 central memory",
  "MAIT",
  "T gd Vd2",
  "Neutrophils",
  "T gd non-Vd2",
  "Basophils LD",
  "Monocytes NC+I",
  "B cells",
  "Dendritic cells",
  "Dendritic cells plasmacytoid",
  "Plasma cells"
)

pkg.globals$qtsq_scaling <- c(
  65.66148,
  138.11520,
  119.35447,
  130.65455,
  27.73634,
  117.71584,
  63.87200,
  70.25659,
  72.55110,
  140.76091,
  68.89323
)
names(pkg.globals$qtsq_scaling) <- c(
  "B cells",
  "Macrophages",
  "MacrophagesM2",
  "Monocytes",
  "Neutrophils",
  "NK cells",
  "T cells CD4",
  "T cells CD8",
  "T regulatory cells",
  "Dendritic cells",
  "T cells"
)

