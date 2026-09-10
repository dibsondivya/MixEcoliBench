# set correct working directory. it should contain
#  - a data folder 
#  - a results folder

set_conf = 0.68

# import libraries
required_packages <- c( "dplyr","readr",  "stringr", "tidyverse",  "ggplot2", "tibble", "knitr", "kableExtra", "rstudioapi" )
missing_packages <- required_packages[ !required_packages %in% rownames(installed.packages()) ] 
if(length(missing_packages) > 0) { install.packages( missing_packages, repos = "https://cloud.r-project.org" ) }
invisible( lapply( required_packages, library, character.only = TRUE ) )
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ifelse(!dir.exists(file.path('figures')), dir.create(file.path('figures')), "Directory Exists")
ifelse(!dir.exists(file.path('tables')), dir.create(file.path('tables')), "Directory Exists")

# gene table ----------------------------------------------------------------------------
arg_genes <- c("blaTEM","sul1","sul2","sul3","tetA","tetB","tetD")
vg_genes  <- c("eae","stx1","stx2")
amrfinder <- read_tsv("results/assembly-amrfinderplus.tsv", show_col_types = FALSE)
rgi <- read_tsv("results/assembly-rgi.tsv", show_col_types = FALSE)
abricate <- read_tsv("results/abricate.tsv", show_col_types = FALSE)
vf <- read_tsv("results/virulencefinder.tsv", show_col_types = FALSE)
arg_srst2 <- read_tsv("results/arg-srst2.tsv", show_col_types = FALSE)
vg_srst2 <- read_tsv("results/virulence-srst2.tsv", show_col_types = FALSE)

# solve differences in gene naming because subtypes
collapse_vg <- function(df, tool) {
  if (tool == "abricate") {df %>% transmute(name, eae = eae, stx1 = pmax(stx1a, stx1b), stx2 = stx2)
  } else if (tool == "virulencefinder") {df %>% transmute(name, eae = eae, stx1 = stx1, stx2 = stx2)
  } else if (tool == "srst2") {df %>% transmute(name, eae = eae, stx1 = pmax(stx1a, stx1b), stx2 = pmax(stx2a, stx2b))
  }
}
abricate_vg <- collapse_vg(abricate, "abricate")
vf_vg       <- collapse_vg(vf, "virulencefinder")
srst2_vg    <- collapse_vg(vg_srst2, "srst2")

# gene recovery rate
calc_grr <- function(df, genes) {
  df %>% rowwise() %>% mutate(GRR = sum(c_across(all_of(genes))) / length(genes)) %>% ungroup() %>% select(name, GRR)
}


arg_tbl <- calc_grr(amrfinder, arg_genes) %>%
  rename(AMRFinderPlus = GRR) %>%
  left_join(calc_grr(rgi, arg_genes) %>% rename(RGI = GRR),by = "name")

arg_tbl

arg_agrr <- c(mean(arg_tbl$AMRFinderPlus), mean(arg_tbl$RGI))

vg_tbl <- calc_grr(abricate_vg, vg_genes) %>% rename(ABRicate = GRR) %>%
  left_join(calc_grr(vf_vg, vg_genes) %>% rename(VirulenceFinder = GRR), by = "name")

vg_tbl

vg_agrr <- c(mean(vg_tbl$ABRicate), mean(vg_tbl$VirulenceFinder))

arg_srst2_tbl <- calc_grr(arg_srst2, arg_genes)
vg_srst2_tbl  <- calc_grr(srst2_vg, vg_genes)

srst2_tbl <- tibble(Input = paste0(sprintf("ARG%02d", 1:10)," | ",sprintf("VG%02d", 1:10)),
                    `SRST2_ARG` = arg_srst2_tbl$GRR,
                    `SRST2_VG` = vg_srst2_tbl$GRR
)

#srst2_tbl

srst2_agrr <- c(mean(srst2_tbl$SRST2_ARG), mean(srst2_tbl$SRST2_VG))

# make arg table
table3_arg <- arg_tbl %>%mutate(AMRFinderPlus = sprintf("%.2f", AMRFinderPlus), RGI = sprintf("%.2f", RGI))
table3_arg <- bind_rows(table3_arg, tibble(name = "AGRR",
                                           AMRFinderPlus = sprintf("%.2f%%", 100 * mean(arg_tbl$AMRFinderPlus)),
                                           RGI = sprintf("%.2f%%", 100 * mean(arg_tbl$RGI))
  )
)

table3_vg <- vg_tbl %>% mutate(ABRicate = sprintf("%.2f", ABRicate), VirulenceFinder = sprintf("%.2f", VirulenceFinder))

table3_vg <- bind_rows(table3_vg, tibble(name = "AGRR",
                                         ABRicate = sprintf("%.2f%%", 100 * mean(vg_tbl$ABRicate)),
                                         VirulenceFinder = sprintf("%.2f%%", 100 * mean(vg_tbl$VirulenceFinder))
  )
)

table3_combined <- bind_rows(tibble(Input = "\\multicolumn{3}{c}{\\textbf{ARG}}", A = "", B = ""),
                             tibble(Input = "Input", A = "\\textbf{AMRFinderPlus}", B = "\\textbf{RGI}"),
                             rename(table3_arg, Input = name, A = AMRFinderPlus, B = RGI),
                             tibble(Input = "", A = "", B = ""),
                             tibble(Input = "\\multicolumn{3}{c}{\\textbf{VG}}", A = "", B = ""),
                             tibble(Input = "Input", A = "\\textbf{ABRicate}", B = "\\textbf{VirulenceFinder}"),
                             rename(table3_vg, Input = name, A = ABRicate, B = VirulenceFinder))

table3_tex <- capture.output(
  kbl(
    table3_combined,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    col.names = c("", "", ""),
    align = "lcc"
  )
)

table3_tex_full <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\caption{Gene recovery rates per dataset for antibiotic resistance gene (ARG) and virulence gene (VG) detection across the assembly methods of AMRFinderPlus, RGI, ABRicate and VirulenceFinder.}",
  "\\label{tab:Table3}",
  table3_tex,
  "\\end{table}"
)

writeLines(table3_tex_full, "tables/table3_gene_recovery_assembly.tex")

table4_tbl <- tibble(Input = paste0( sprintf("ARG%02d", 1:10), " $\\mid$ ", sprintf("VG%02d", 1:10)),
  ARG = sprintf("%.2f", arg_srst2_tbl$GRR),
  VG = sprintf("%.2f", vg_srst2_tbl$GRR)
)

table4_tbl <- bind_rows(
  table4_tbl,
  tibble(Input = "AGRR", ARG = sprintf("%.2f%%", 100 * mean(arg_srst2_tbl$GRR)), VG = sprintf("%.2f%%", 100 * mean(vg_srst2_tbl$GRR)))
)

table4_tex <- capture.output(
  kbl(
    table4_tbl,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = "lcc",
    col.names = c(
      "Input",
      "SRST2 (ARG)",
      "SRST2 (VG)"
    )
  )
)

table4_tex_full <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\caption{Gene recovery rates for antibiotic resistance gene (ARG) and virulence-associated gene (VG) detection using read-based SRST2.}",
  "\\label{tab:Table4}",
  table4_tex,
  "\\end{table}"
)

writeLines(table4_tex_full, "tables/table4_gene_recovery_srst2.tex")


# data description  ----------------------------------------------------------------------------
fixed <- read_tsv("data/st_serotype_mixed_range_fixed_manifest.tsv",show_col_types = FALSE)
rare <- read_tsv("data/st_serotype_mixed_range_rare_manifest.tsv", show_col_types = FALSE)

# make st, serotype description table -----------------------------------------------------------------------------------
# main descriptors of interest: dataset, sample #, ST, serotype
fixed_plot <- fixed %>%
  transmute(
    Dataset = case_when(
      dominant_percentage == 100 ~ "Single-Strain (100%)",
      dominant_percentage == 90  ~ "Dominant-Strain (90%)",
      dominant_percentage == 80  ~ "Dominant-Strain (80%)",
      dominant_percentage == 70  ~ "Dominant-Strain (70%)",
      dominant_percentage == 60  ~ "Dominant-Strain (60%)"
    ),
    simulation,
    ST = paste0("ST", dominant_ST),
    ST_Proportion = realised_ST_abundance,
    Serotype = dominant_Serotype,
    Serotype_Proportion = realised_Serotype_abundance
  )

rare_plot <- rare %>%
  transmute(
    Dataset = "Multi-Strain (Log-normal)",
    simulation,
    ST = paste0("ST", dominant_ST),
    ST_Proportion = dominant_ST_abundance,
    Serotype = dominant_Serotype,
    Serotype_Proportion = dominant_Serotype_abundance
  )

combined_st_sero_data <- bind_rows(fixed_plot, rare_plot) %>%
  mutate(
    Dataset = factor(
      Dataset,
      levels = c(
        "Single-Strain (100%)",
        "Dominant-Strain (90%)",
        "Dominant-Strain (80%)",
        "Dominant-Strain (70%)",
        "Dominant-Strain (60%)",
        "Multi-Strain (Log-normal)"
      )
    )
  )

# produce the st, serotype description table
combined_st_sero_data

# make mlst results table -----------------------------------------------------------------------------------
sample_to_name <- function(x) {
  sim <- as.integer(str_extract(x, "(?<=sim_)\\d+"))
  case_when(
    str_detect(x, "range_rare") ~ sprintf("RrST%02d", sim),
    str_detect(x, "range_60")   ~ sprintf("Mx60ST%02d", sim),
    str_detect(x, "range_70")   ~ sprintf("Mx70ST%02d", sim),
    str_detect(x, "range_80")   ~ sprintf("Mx80ST%02d", sim),
    str_detect(x, "range_90")   ~ sprintf("Mx90ST%02d", sim),
    str_detect(x, "range_100")  ~ sprintf("AbST%02d", sim),
    TRUE ~ NA_character_
  )
}

mlst <- read_tsv("results/ecoli_torsten_mlst.tsv", show_col_types = FALSE) %>%
  mutate(name = sample_to_name(FILE), mlst = ST) %>% select(name, mlst)
mlst_cge <- read_tsv("results/compiled_mlst-cge_mlst.tsv", show_col_types = FALSE) %>%
  mutate(name = sample_to_name(filename), `mlst-cge` = sequence_type) %>% select(name, `mlst-cge`)
srst2 <- read_tsv("results/compiled_srst2_mlst.tsv", show_col_types = FALSE) %>%
  mutate(name = sample_to_name(sample), srst2 = ST) %>% select(name, srst2)

mlst_truth <- combined_st_sero_data %>%
  mutate(
    name = case_when(
      Dataset == "Multi-Strain (Log-normal)" ~ sprintf("RrST%02d", simulation),
      Dataset == "Dominant-Strain (60%)"     ~ sprintf("Mx60ST%02d", simulation),
      Dataset == "Dominant-Strain (70%)"     ~ sprintf("Mx70ST%02d", simulation),
      Dataset == "Dominant-Strain (80%)"     ~ sprintf("Mx80ST%02d", simulation),
      Dataset == "Dominant-Strain (90%)"     ~ sprintf("Mx90ST%02d", simulation),
      Dataset == "Single-Strain (100%)"      ~ sprintf("AbST%02d", simulation)
    ),
    truth_dominant = str_remove(ST, "^ST"),
    truth_abundance = ST_Proportion
  ) %>%
  select(name, truth_dominant, truth_abundance)

mlst_summary <- mlst_truth %>%
  left_join(mlst, by = "name") %>%
  left_join(mlst_cge, by = "name") %>%
  left_join(srst2, by = "name") %>%
  mutate(metamlst = "-") %>%
  select(name, mlst, `mlst-cge`, truth_dominant, truth_abundance, metamlst, srst2) %>%
  mutate(
    group = case_when(
      str_detect(name, "^Rr")   ~ 1,
      str_detect(name, "^Mx60") ~ 2,
      str_detect(name, "^Mx70") ~ 3,
      str_detect(name, "^Mx80") ~ 4,
      str_detect(name, "^Mx90") ~ 5,
      str_detect(name, "^Ab")   ~ 6
    )
  ) %>%
  arrange(group, name) %>%
  select(-group)

# produce mlst summary
mlst_summary

# look at mlst predictions for rare dataset  -----------------------------------------------------------------------------------
mlst_predictions_tbl <- mlst_summary %>% filter(str_detect(name, "^Rr")) %>%
  mutate( Input = sprintf("ST%02d", seq_len(n()))) %>%
  transmute(Input, mlst, `mlst-cge`, metaMLST = metamlst, SRST2 = srst2, `Dominant ST` = truth_dominant)

mlst_predictions_tex <- capture.output(
  kable(
    mlst_predictions_tbl,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = "lccccc",
    caption = paste(
      "MLST predictions returned by assembly-based and read-based",
      "tools for simulated metagenomes."
    ),
    label = "tab:MLST_predictions"
  )
  
)


is_assigned <- function(x) {
  !is.na(x) & !str_detect(x, "NF") & !(x %in% c("-", "Unknown","Not Found"))
}

rare_tool_summary <- mlst_summary %>% filter(str_detect(name, "^Rr")) %>%
  pivot_longer(cols = c(mlst, `mlst-cge`, metamlst, srst2), names_to = "tool", values_to = "call") %>%
  mutate(assigned = is_assigned(call),
         correct = map2_lgl(call, truth_dominant, is_correct_st)) %>%
  group_by(tool) %>%
  summarise(datasets = n(),
            assignments = sum(assigned),
            unassigned = sum(!assigned),
            correct_assignments = sum(correct),
            accuracy = 100 * mean(correct),
            .groups = "drop"
  )

rare_summary_tbl <- tibble(
  Metric = c("Datasets analysed", "ST assignments", "Unassigned ST calls", "Correct ST assignments", "ST assignment accuracy"),
  mlst = c(
    rare_tool_summary$datasets[rare_tool_summary$tool=="mlst"],
    rare_tool_summary$assignments[rare_tool_summary$tool=="mlst"],
    rare_tool_summary$unassigned[rare_tool_summary$tool=="mlst"],
    rare_tool_summary$correct_assignments[rare_tool_summary$tool=="mlst"],
    paste0( round(rare_tool_summary$accuracy[rare_tool_summary$tool=="mlst"],1),"\\%")
  ),
  `mlst-cge` = c(
    rare_tool_summary$datasets[rare_tool_summary$tool=="mlst-cge"],
    rare_tool_summary$assignments[rare_tool_summary$tool=="mlst-cge"],
    rare_tool_summary$unassigned[rare_tool_summary$tool=="mlst-cge"],
    rare_tool_summary$correct_assignments[rare_tool_summary$tool=="mlst-cge"],
    paste0(round(rare_tool_summary$accuracy[rare_tool_summary$tool=="mlst-cge"],1),"\\%")
  ),
  
  metaMLST = c(
    rare_tool_summary$datasets[rare_tool_summary$tool=="metamlst"],
    rare_tool_summary$assignments[rare_tool_summary$tool=="metamlst"],
    rare_tool_summary$unassigned[rare_tool_summary$tool=="metamlst"],
    rare_tool_summary$correct_assignments[rare_tool_summary$tool=="metamlst"],
    paste0(round(rare_tool_summary$accuracy[rare_tool_summary$tool=="metamlst"],1),"\\%")
  ),
  
  SRST2 = c(
    rare_tool_summary$datasets[rare_tool_summary$tool=="srst2"],
    rare_tool_summary$assignments[rare_tool_summary$tool=="srst2"],
    rare_tool_summary$unassigned[rare_tool_summary$tool=="srst2"],
    rare_tool_summary$correct_assignments[rare_tool_summary$tool=="srst2"],
    paste0(round(rare_tool_summary$accuracy[rare_tool_summary$tool=="srst2"],1),"\\%")
  )
  
)


rare_summary_tex <- capture.output(
  kbl(
    rare_summary_tbl,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = "lcccc"
  ) %>%
    add_header_above(c(" " = 1, "Assembly-based" = 2, "Read-based" = 2))
)

rare_summary_tex_full <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\caption{Summary performance of MLST prediction tools across the ten multi-strain (log-normal) metagenomes.}",
  "\\label{tab:MLST_summary_rare}",
  rare_summary_tex,
  "\\end{table}"
)


writeLines(mlst_predictions_tex, "tables/table5_mlst_predictions.tex")
writeLines(rare_summary_tex_full, "tables/supplementary_table2_rare_mlst_summary.tex")


# make serotype results table -----------------------------------------------------------------------------------
ectyper_assembly <- read_tsv("results/assembly_ectyper_serotype.tsv", show_col_types = FALSE) %>%
  mutate(name = sample_to_name(name), ectyper = serotype) %>% select(name, ectyper)
serotypefinder <- read_tsv("results/assembly_serotypefinder_serotype.tsv", show_col_types = FALSE) %>%
  mutate(name = sample_to_name(sample), serotypefinder = serotype) %>% select(name, serotypefinder)
ectyper_reads <- read_tsv("results/ectyper_combined_serotypes.tsv", show_col_types = FALSE) %>%
  mutate(name = sample_to_name(sample), ectyper_read = serotype) %>% select(name, ectyper_read)
srst2_serotype <- read_tsv("results/compiled_srst2_serotype.tsv", show_col_types = FALSE) %>%
  mutate(name = sample_to_name(sample), srst2 = serotype) %>% select(name, srst2)

truth_serotype <- combined_st_sero_data %>%
  mutate(
    name = case_when(
      Dataset == "Single-Strain (100%)" ~ sprintf("AbST%02d", simulation),
      Dataset == "Multi-Strain (Log-normal)" ~ sprintf("RrST%02d", simulation),
      Dataset == "Dominant-Strain (60%)" ~ sprintf("Mx60ST%02d", simulation),
      Dataset == "Dominant-Strain (70%)" ~ sprintf("Mx70ST%02d", simulation),
      Dataset == "Dominant-Strain (80%)" ~ sprintf("Mx80ST%02d", simulation),
      Dataset == "Dominant-Strain (90%)" ~ sprintf("Mx90ST%02d", simulation)
    ),
    truth_dominant = Serotype,
    truth_abundance = Serotype_Proportion
  ) %>%
  select(
    name,
    truth_dominant,
    truth_abundance
  )

serotype_summary <- truth_serotype %>%
  left_join(ectyper_assembly, by = "name") %>%
  left_join(serotypefinder, by = "name") %>%
  left_join(ectyper_reads, by = "name") %>%
  left_join(srst2_serotype, by = "name") %>%
  select(name, ectyper, serotypefinder, truth_dominant, truth_abundance, ectyper_read, srst2) %>%
  mutate(
    group = case_when(
      str_detect(name, "^Rr")   ~ 1,
      str_detect(name, "^Mx60") ~ 2,
      str_detect(name, "^Mx70") ~ 3,
      str_detect(name, "^Mx80") ~ 4,
      str_detect(name, "^Mx90") ~ 5,
      str_detect(name, "^Ab")   ~ 6
    )
  ) %>%
  arrange(group, name) %>%
  select(-group)

# produce serotype summary
serotype_summary
View(serotype_summary)

# make serotype antigen results table -----------------------------------------------------------------------------------
split_O <- function(x) ifelse(x == "-", "-", str_split_fixed(x, ":", 2)[,1])
split_H <- function(x) ifelse(x == "-", "-", str_split_fixed(x, ":", 2)[,2])

rr_perf <- serotype_summary %>%
  filter(str_detect(name, "^Rr")) %>%
  mutate( truth_O = split_O(truth_dominant), truth_H = split_H(truth_dominant)) %>%
  pivot_longer(cols = c(ectyper, serotypefinder, ectyper_read, srst2), names_to = "tool", values_to = "call") %>%
  mutate(call_O = split_O(call), 
         call_H = split_H(call), 
         O_correct = map2_lgl(call_O, truth_O, ~ .y %in% str_split(.x, "/")[[1]]), 
         H_correct = map2_lgl(call_H, truth_H, ~ .y %in% str_split(.x, "/")[[1]])
  )

rr_summary <- rr_perf %>% group_by(tool) %>%
  summarise(
    O_k = sum(O_correct),
    H_k = sum(H_correct),
    Overall_k = O_k + H_k,
    Complete_k = sum(O_correct & H_correct),
    .groups = "drop"
  ) %>%
  mutate(
    O = sprintf("%d/10 (%.1f%%)", O_k, 100 * O_k / 10),
    H = sprintf("%d/10 (%.1f%%)", H_k, 100 * H_k / 10),
    Overall = sprintf("%d/20 (%.1f%%)", Overall_k, 100 * Overall_k / 20),
    Complete = sprintf("%d/10 (%.1f%%)", Complete_k, 100 * Complete_k / 10)
  ) %>%
  select(tool, O, H, Overall, Complete) %>%
  mutate(
    tool = recode(
      tool,
      ectyper = "ECTyper",
      serotypefinder = "SerotypeFinder",
      ectyper_read = "ECTyper (reads)",
      srst2 = "SRST2"
    )
  )
# produce antigen results summary for serotypes
rr_summary

# write latex code
rr_kable <- rr_summary %>%
  mutate(tool = factor(tool, levels = c("ECTyper","SerotypeFinder","ECTyper (reads)","SRST2"))) %>%
  arrange(tool) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    col.names = c(
      "Tool",
      "O-antigen",
      "H-antigen",
      "Overall antigen recovery",
      "Complete serotype recovery"
    )
  )

rr_latex <- paste(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\caption{Performance of serotype prediction tools. O-antigen recovery and H-antigen recovery were calculated as the proportion of expected antigens recovered. Complete serotype recovery required recovery of both the expected O and H antigens within a sample.}",
  "\\label{tab:serotype_summary}",
  rr_kable,
  "\\end{table}",
  sep = "\n"
)

writeLines(rr_latex,"tables/table6_table_overall_serotype_rare.tex")


# make figures -----------------------------------------------------------------------------------
wilson_ci <- function(k, n, conf = set_conf) {
  z <- qnorm(1 - (1 - conf) / 2)
  phat <- k / n
  denom <- 1 + z^2 / n
  centre <- (phat + z^2 / (2 * n)) / denom
  margin <- (z / denom) * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2))
  tibble(ymin = pmax(0, centre - margin), ymax = pmin(1, centre + margin))
}

make_accuracy_plot <- function(data, tools, truth_col, abundance_col, score_fun, tool_labels, assembly_tools, ylab, mode = c("simple", "serotype")) {
  mode <- match.arg(mode)
  perf_data <- data %>%
    filter(!str_detect(name, "^Rr")) %>%
    pivot_longer(cols = all_of(tools), names_to = "tool", values_to = "call") %>%
    mutate(abundance_pct = .data[[abundance_col]] * 100, 
           category = if_else(tool %in% assembly_tools, "Assembly-based", "Read-based"),
           tool = recode(tool, !!!tool_labels),
           tool = factor(tool, levels = unname(tool_labels))
    )
  
  if (mode == "simple") {
    perf_data <- perf_data %>% mutate(correct = map2_lgl(call, .data[[truth_col]], score_fun))
    summary_data <- perf_data %>% group_by(tool, category, abundance_pct) %>% 
      summarise(n = n(), 
                k = sum(correct),
                accuracy = mean(correct),
                .groups = "drop") %>%
      rowwise() %>%
      mutate(ci = list(wilson_ci(k, n))) %>%
      unnest(ci)
  }
  
  if (mode == "serotype") {
    split_O <- function(x) ifelse(x == "-", "-", str_split_fixed(x, ":", 2)[,1])
    split_H <- function(x) ifelse(x == "-", "-", str_split_fixed(x, ":", 2)[,2])
    
    perf_data <- perf_data %>% mutate( truth_O = split_O(.data[[truth_col]]),
                                       truth_H = split_H(.data[[truth_col]]),
                                       call_O = split_O(call),
                                       call_H = split_H(call),
                                       O_correct = map2_lgl(call_O, truth_O, ~ .y %in% str_split(.x, "/")[[1]]),
                                       H_correct = map2_lgl(call_H, truth_H, ~ .y %in% str_split(.x, "/")[[1]]))
    by_group <- perf_data %>% group_by(tool, category, abundance_pct)
    summary_data <- bind_rows(by_group %>% summarise(k=sum(O_correct),n=n(),.groups="drop") %>% mutate(metric="O-antigen"),
                              by_group %>% summarise(k=sum(H_correct),n=n(),.groups="drop") %>% mutate(metric="H-antigen"),
                              by_group %>% summarise(k=sum(O_correct)+sum(H_correct), n=2*n(), .groups="drop") %>%
                                mutate(metric="Overall antigen recovery"),
                              by_group %>% summarise( k=sum(O_correct & H_correct), n=n(), .groups="drop") %>%
                                mutate(metric="Complete recovery")) %>%
      mutate(accuracy = k/n, 
             metric = factor(metric, levels = c("O-antigen", "H-antigen", "Overall antigen recovery", "Complete recovery"))) %>%
      rowwise() %>%
      mutate(ci = list(wilson_ci(k, n))) %>%
      unnest(ci)
  }
  
  if (mode == "simple") {
    stopifnot(!is.null(score_fun))
    p <- ggplot(summary_data, aes(x = abundance_pct, y = accuracy, colour = category, group = tool)) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 2, linewidth = 0.4, alpha = 0.7) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 2, shape = 21, fill = "white", stroke = 0.9) +
      scale_colour_manual(values = c("Assembly-based" = "#1b6ca8", "Read-based" = "#c2571a")) +
      scale_x_continuous(breaks = c(60,70,80,90,100), limits = c(58,102)) +
      scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1)) +
      labs(x = "Dominant-strain abundance (%)", y = ylab) +
      theme_bw() +
      theme(legend.position = "bottom") +
      facet_wrap(~tool, nrow = 1)
    
  } else {
    p <- ggplot(summary_data, aes(x = abundance_pct, y = accuracy, colour = category, group = tool)) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 2, linewidth = 0.4, alpha = 0.7) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 2, shape = 21, fill = "white", stroke = 0.9) +
      scale_colour_manual(values = c("Assembly-based" = "#1b6ca8", "Read-based" = "#c2571a")) +
      scale_x_continuous(breaks = c(60,70,80,90,100), limits = c(58,102)) +
      scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1)) +
      labs(x = "Dominant-strain abundance (%)", y = ylab) +
      theme_bw() +
      theme(legend.position = "bottom") +
      facet_grid(rows = vars(metric), cols = vars(tool))
  }
  
  return(list(perf_data = perf_data, summary_data = summary_data, plot = p))
}


is_correct_st <- function(call, truth) {
  if (is.na(call))
    return(FALSE)
  if (call %in% c("-", "Unknown", "NF", "NF?", "NF*?"))
    return(FALSE)
  truth %in% str_extract_all(call, "\\d+")[[1]]
}

mlst_results <- make_accuracy_plot(
  data = mlst_summary,
  tools = c("mlst", "mlst-cge", "metamlst", "srst2"),
  truth_col = "truth_dominant",
  abundance_col = "truth_abundance",
  score_fun = is_correct_st,
  assembly_tools = c("mlst", "mlst-cge"),
  tool_labels = c("mlst" = "mlst", "mlst-cge" = "mlst-cge", "metamlst" = "metaMLST", "srst2" = "SRST2"),
  ylab = "Sequence type prediction accuracy",
  mode = 'simple'
)

mlst_results$plot
ggsave("figures/fig1_mixed_range_mlst_performance.png", mlst_results$plot, width = 180, height = 100, units = "mm", dpi = 300)

serotype_results <- make_accuracy_plot(
  data = serotype_summary,
  tools = c("ectyper", "serotypefinder", "ectyper_read", "srst2"),
  truth_col = "truth_dominant",
  abundance_col = "truth_abundance",
  score_fun = NULL,
  assembly_tools = c("ectyper", "serotypefinder"),
  tool_labels = c(
    "ectyper"="ECTyper (assembly)",
    "serotypefinder"="SerotypeFinder",
    "ectyper_read"="ECTyper (reads)",
    "srst2"="SRST2"
  ),
  ylab = "Prediction accuracy",
  mode = "serotype"
)

serotype_results$plot
ggsave("figures/fig2_mixed_range_serotype_performance.png", serotype_results$plot, width = 180, height = 190, units = "mm", dpi = 300)

# make supplementary data description figures -----------------------------------------------------------------------------------
p_st <- ggplot(combined_st_sero_data, aes(x = factor(simulation), y = ST_Proportion, fill = ST)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = scales::percent(ST_Proportion, accuracy = 1), vjust = ifelse(ST_Proportion == 1, 1.2, -0.3)), size = 3) +
  facet_wrap(~Dataset, ncol = 1) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Simulated metagenome", y = "Proportion of dominant ST", fill = "Dominant ST") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())

p_sero <- ggplot(combined_st_sero_data, aes(x = factor(simulation), y = Serotype_Proportion, fill = Serotype)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = scales::percent(Serotype_Proportion, accuracy = 1), vjust = ifelse(Serotype_Proportion == 1, 1.2, -0.3)), size = 3) +
  facet_wrap(~Dataset, ncol = 1) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Simulated metagenome", y = "Proportion of dominant serotype", fill = "Dominant serotype") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave("figures/supplementary_fig1_dominant_st_description.png", p_st, width = 290, height = 300, units = "mm", dpi = 300)
ggsave("figures/supplementary_fig2_dominant_serotype_description.png", p_sero, width = 290, height = 300, units = "mm", dpi = 300)


# make supplementary serotype tables -----------------------------------------------------------------------------------
perf_data_tbl <- serotype_summary %>%
  mutate(truth_O = split_O(truth_dominant), truth_H = split_H(truth_dominant)) %>%
  pivot_longer(cols = c(ectyper, serotypefinder, ectyper_read,srst2), names_to = "tool", values_to = "call") %>%
  mutate(
    call_O = split_O(call),
    call_H = split_H(call),
    O_correct = map2_lgl(call_O, truth_O, ~ .y %in% str_split(.x, "/")[[1]]),
    H_correct = map2_lgl(call_H, truth_H, ~ .y %in% str_split(.x, "/")[[1]]),
    Dataset = case_when(
      str_detect(name, "^Ab")   ~ "Single-strain (100\\%)",
      str_detect(name, "^Mx90") ~ "Dominant-strain (90\\%)",
      str_detect(name, "^Mx80") ~ "Dominant-strain (80\\%)",
      str_detect(name, "^Mx70") ~ "Dominant-strain (70\\%)",
      str_detect(name, "^Mx60") ~ "Dominant-strain (60\\%)",
      str_detect(name, "^Rr")   ~ "Multi-strain"
    ),
    tool = recode(tool, ectyper = "ECTyper", serotypefinder = "SerotypeFinder", ectyper_read = "ECTyper_read", srst2 = "SRST2"
    )
  )


make_metric_table <- function(data, metric = c("overall", "O", "H","complete")) {
  metric <- match.arg(metric)
  dat <- data
  summary_tbl <- switch(metric,
    O = dat %>% group_by(Dataset, tool) %>% summarise(value = mean(O_correct), .groups = "drop"),
    H = dat %>% group_by(Dataset, tool) %>% summarise(value = mean(H_correct), .groups = "drop"),
    complete = dat %>% group_by(Dataset, tool) %>% summarise(value = mean(O_correct & H_correct), .groups = "drop"),
    overall = dat %>% group_by(Dataset, tool) %>% summarise(value =(sum(O_correct) + sum(H_correct)) / (2 * n()), .groups = "drop")
  )
  
  summary_tbl %>% mutate(value = sprintf("%.1f\\%%", 100 * value)) %>%
    pivot_wider(names_from = tool, values_from = value) %>%
    select(Dataset, ECTyper, SerotypeFinder, ECTyper_read, SRST2) %>%
    mutate(
      Dataset = factor(
        Dataset,
        levels = c(
          "Single-strain (100\\%)",
          "Dominant-strain (90\\%)",
          "Dominant-strain (80\\%)",
          "Dominant-strain (70\\%)",
          "Dominant-strain (60\\%)",
          "Multi-strain"
        )
      )
    ) %>%
    arrange(Dataset)
}

overall_tbl <- make_metric_table(perf_data_tbl, "overall")
o_tbl <- make_metric_table(perf_data_tbl, "O")
h_tbl <- make_metric_table(perf_data_tbl, "H")
complete_tbl <- make_metric_table(perf_data_tbl, "complete")

overall_kable  <- paste(capture.output(kable(overall_tbl, escape = FALSE, format="latex", booktabs=TRUE)), collapse="\n")
o_kable        <- paste(capture.output(kable(o_tbl, escape = FALSE, format="latex", booktabs=TRUE)), collapse="\n")
h_kable        <- paste(capture.output(kable(h_tbl, escape = FALSE, format="latex", booktabs=TRUE)), collapse="\n")
complete_kable <- paste(capture.output(kable(complete_tbl, escape = FALSE, format="latex", booktabs=TRUE)), collapse="\n")

overall_tex <- paste(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\caption{Overall antigen recovery rates across dataset compositions.}",
  "\\label{tab:overall_accuracy}",
  overall_kable,
  "\\end{table}",
  sep = "\n"
)

o_latex <- paste(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\caption{Overall O-antigen recovery rates across dataset compositions.}",
  "\\label{tab:O_accuracy}",
  o_kable,
  "\\end{table}",
  sep = "\n"
)

h_latex <- paste(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\caption{Overall H-antigen recovery rates across dataset compositions.}",
  "\\label{tab:H_accuracy}",
  h_kable,
  "\\end{table}",
  sep = "\n"
)

complete_latex <- paste(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\caption{Overall complete antigen recovery rates across dataset compositions.}",
  "\\label{tab:Complete_accuracy}",
  complete_kable,
  "\\end{table}",
  sep = "\n"
)

writeLines(overall_tex,  "tables/supplementary_table3_serotype_overall.tex")
writeLines(o_latex,        "tables/supplementary_table4_serotype_O.tex")
writeLines(h_latex,        "tables/supplementary_table5_serotype_H.tex")
writeLines(complete_latex, "tables/supplementary_table6_serotype_complete.tex")
