# Multi-omics integration analysis

.joint_req <- function(p){m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];if(length(m))stop("Missing required R packages: ",paste(m,collapse=", "))}
.joint_out <- function(p){if(dir.exists(p)&&length(list.files(p,all.files=TRUE,no..=TRUE)))stop("Output directory is not empty: ",p);dir.create(p,recursive=TRUE,showWarnings=FALSE);p}
.joint_tmp <- function(p){if(!file.exists(p))stop("Input file not found: ",p);e<-tools::file_ext(p);z<-tempfile(fileext=if(nzchar(e))paste0(".",e)else"");if(!file.copy(p,z,overwrite=TRUE))stop("Failed to read input file: ",p);z}
.joint_csv <- function(p){z<-.joint_tmp(p);on.exit(unlink(z));read.csv(z,check.names=FALSE,stringsAsFactors=FALSE)}
.joint_tsv <- function(p){z<-.joint_tmp(p);on.exit(unlink(z));read.delim(z,check.names=FALSE,stringsAsFactors=FALSE)}
.joint_table <- function(p,sheet="Data"){
  ext<-tolower(tools::file_ext(p))
  if(ext%in%c("xlsx","xls")){z<-.joint_tmp(p);on.exit(unlink(z));as.data.frame(readxl::read_excel(z,sheet=sheet),check.names=FALSE)}
  else if(ext=="csv"){.joint_csv(p)}
  else{.joint_tsv(p)}
}
.joint_pair <- function(x,s){write.csv(x,paste0(s,".csv"),row.names=FALSE,fileEncoding="UTF-8");openxlsx::write.xlsx(x,paste0(s,".xlsx"),overwrite=TRUE)}
.joint_save <- function(p,s,w,h){ggplot2::ggsave(paste0(s,".pdf"),p,width=w,height=h,device=grDevices::cairo_pdf);tmp<-tempfile(fileext=".tiff");on.exit(unlink(tmp),add=TRUE);ggplot2::ggsave(tmp,p,width=w,height=h,dpi=600,compression="lzw");if(!file.copy(tmp,paste0(s,".tiff"),overwrite=TRUE))stop("Failed to write TIFF file: ",s)}
.joint_kegg <- function(gene_file,name_file){g<-.joint_tsv(gene_file);n<-.joint_tsv(name_file);if(!all(c("term","gene")%in%names(g))||!all(c("term","term_name")%in%names(n)))stop("KEGG mapping tables are missing required columns.");g$term<-sub("^ko","map",g$term);g$gene<-as.character(g$gene);n$term<-sub("^ko","map",n$term);list(g=unique(g[,c("term","gene")]),n=unique(n[,c("term","term_name")]))}
.joint_enrich <- function(genes,universe,k){
  universe<-intersect(unique(as.character(universe)),unique(as.character(k$g$gene)));fg_all<-unique(as.character(genes));fg<-intersect(fg_all,universe);map<-merge(k$g,k$n,by.x="term",by.y="term",all.x=TRUE);map<-map[map$gene%in%universe,];sets<-split(map$gene,map$term);sets<-lapply(sets,unique);sets<-sets[lengths(sets)>=10 & lengths(sets)<=500];N<-length(universe);n<-length(fg);r<-do.call(rbind,lapply(names(sets),function(id){bg<-sets[[id]];hit<-intersect(bg,fg);M<-length(bg);z<-length(hit);data.frame(pathway_id=id,pathway_name=map$term_name[match(id,map$term)],Bg_Count=M,Count=z,Pvalue=phyper(z-1,M,N-M,n,lower.tail=FALSE),GeneRatio=z/n,BgRatio=M/N,FoldEnrichment=(z/n)/(M/N),RichFactor=z/M,hit_geneIDs=paste(sort(hit),collapse="/"),Foreground_Total=length(fg_all),Foreground_Annotated=n,Background_Annotated_Tested=N,Min_GS_Size=10,Max_GS_Size=500,stringsAsFactors=FALSE)}));r$FDR<-p.adjust(r$Pvalue,"BH");r[order(r$Pvalue,-r$FoldEnrichment),]
}

# Stage 01 - Fig.5
stage_01_fig5 <- function(genes_file, metabolites_file, output_dir,
                          metabolite_id_file = NULL, common_samples = NULL,
                          gene_sheet = "Data", metabolite_sheet = "Data",
                          id_sheet = "Data", save_plot = TRUE) {
  .joint_req(c("readxl", "dplyr", "tidyr", "purrr", "tibble", "igraph",
               "ggraph", "ggplot2", "ggrepel", "scales"))
  suppressPackageStartupMessages({
    library(readxl); library(dplyr); library(tidyr); library(purrr); library(tibble)
    library(igraph); library(ggraph); library(ggplot2); library(ggrepel); library(scales)
  })

  out_dir <- .joint_out(output_dir)

  gene_color <- "#1F5C50"
  metabolite_color <- "#6A3D8F"
  positive_color <- "#E31A1C"
  negative_color <- "#0055D4"

  genes_raw <- .joint_table(genes_file, sheet = gene_sheet)
  metabolites_raw <- .joint_table(metabolites_file, sheet = metabolite_sheet)

  if (!"gene_name" %in% names(genes_raw)) stop("Gene table must contain a 'gene_name' column.")
  if (!"Metabolite" %in% names(metabolites_raw)) stop("Metabolite table must contain a 'Metabolite' column.")

  if (is.null(common_samples)) {
    gene_cand <- setdiff(names(genes_raw), "gene_name")
    gene_num <- gene_cand[vapply(genes_raw[gene_cand], is.numeric, logical(1))]
    met_cand <- setdiff(names(metabolites_raw), c("Metabolite", "Compound_ID"))
    met_num <- met_cand[vapply(metabolites_raw[met_cand], is.numeric, logical(1))]
    common_samples <- intersect(gene_num, met_num)
  } else {
    common_samples <- as.character(common_samples)
    if (!all(common_samples %in% names(genes_raw)) || !all(common_samples %in% names(metabolites_raw))) {
      stop("Some specified common_samples are missing from the gene or metabolite table.")
    }
  }
  if (length(common_samples) < 3L) {
    stop("At least 3 shared numeric sample columns are required for Pearson correlation analysis.")
  }

  if (!is.null(metabolite_id_file)) {
    raw_metabolite_ids <- .joint_table(metabolite_id_file, sheet = id_sheet) %>%
      transmute(Metabolite = trimws(as.character(Metabolite)),
                Compound_ID = as.character(Compound_ID))
  } else if ("Compound_ID" %in% names(metabolites_raw)) {
    raw_metabolite_ids <- metabolites_raw %>%
      transmute(Metabolite = trimws(as.character(Metabolite)),
                Compound_ID = as.character(Compound_ID)) %>%
      distinct()
  } else {
    raw_metabolite_ids <- tibble(Metabolite = character(), Compound_ID = character())
  }

  standard_abbrev <- c(
    "S-Lactoylglutathione" = "SLG",
    "Lysopc 14:0" = "LPC 14:0", "L-Glutathione (reduced)" = "GSH",
    "L-Ergothioneine" = "EGT",
    "Homotaurine" = "3-APS",
    "L-Glutathione oxidized" = "GSSG", "Leukotriene C4" = "LTC4"
  )
  abbrev_sources <- c(
    "SLG" = "https://pubmed.ncbi.nlm.nih.gov/24333633/",
    "GSH" = "https://pubmed.ncbi.nlm.nih.gov/7795391/",
    "gamma-GluCys" = "https://pubmed.ncbi.nlm.nih.gov/7795391/",
    "EGT" = "https://pubmed.ncbi.nlm.nih.gov/40220171/",
    "3-APS" = "https://pubchem.ncbi.nlm.nih.gov/compound/Tramiprosate",
    "GSSG" = "https://pubmed.ncbi.nlm.nih.gov/35204068/",
    "LTC4" = "https://www.lipidmaps.org/databases/shorthand/overview",
    "LPC 14:0" = "https://www.lipidmaps.org/databases/shorthand/overview"
  )
  metabolite_labels <- tibble(Metabolite = trimws(as.character(metabolites_raw$Metabolite))) %>%
    left_join(raw_metabolite_ids, by = "Metabolite") %>%
    mutate(display_name = case_when(
      grepl("Glutamylcysteine$", Metabolite) ~ "gamma-GluCys",
      Metabolite %in% names(standard_abbrev) ~ unname(standard_abbrev[Metabolite]),
      TRUE ~ Metabolite
    ),
    label_basis = case_when(
      grepl("Glutamylcysteine$", Metabolite) ~ "literature/database abbreviation",
      Metabolite %in% names(standard_abbrev) ~ "literature/database abbreviation",
      grepl("^(LPC |LPC O-|LPS |PC |PC O-|Cer |Hex2Cer |SPB )", Metabolite) ~
        "LIPID MAPS shorthand retained from input",
      TRUE ~ "full name retained; no established abbreviation verified"
    ),
    source = case_when(
      display_name %in% names(abbrev_sources) ~ unname(abbrev_sources[display_name]),
      grepl("^(LPC |LPC O-|LPS |PC |PC O-|Cer |Hex2Cer |SPB )", Metabolite) ~
        "https://www.lipidmaps.org/databases/shorthand/overview",
      TRUE ~ ""
    )) %>%
    select(display_name, Metabolite, Compound_ID, label_basis, source)
  write.csv(metabolite_labels, file.path(out_dir, "metabolite_abbreviation_mapping.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  metabolite_labels <- read.csv(file.path(out_dir, "metabolite_abbreviation_mapping.csv"),
                                stringsAsFactors = FALSE, check.names = FALSE)

  prepare_matrix <- function(df, id_col, sample_cols) {
    ids <- trimws(as.character(df[[id_col]]))
    if (anyNA(ids) || any(ids == "") || anyDuplicated(ids)) {
      stop("Empty or duplicate identifiers detected in column: ", id_col)
    }
    x <- df %>% select(all_of(sample_cols)) %>% mutate(across(everything(), as.numeric))
    if (anyNA(x)) stop("Missing or non-numeric values detected in shared sample columns for: ", id_col)
    if (any(as.matrix(x) < 0)) stop("Negative abundance values cannot be log2(x+1)-transformed for: ", id_col)
    log_mat <- log2(as.matrix(x) + 1)
    rownames(log_mat) <- ids
    row_sd <- apply(log_mat, 1, sd)
    if (any(row_sd == 0)) stop("Zero-variance features detected for: ", id_col)
    z_mat <- t(scale(t(log_mat)))
    list(log = log_mat, z = z_mat)
  }

  gene <- prepare_matrix(genes_raw, "gene_name", common_samples)
  metabolite <- prepare_matrix(metabolites_raw, "Metabolite", common_samples)

  edge_tests <- function(mat1, mat2 = NULL, type, r_cutoff, p_cutoff, within = FALSE) {
    if (within) {
      pairs <- t(combn(rownames(mat1), 2))
      result <- as_tibble(pairs, .name_repair = ~c("from", "to")) %>%
        mutate(test = map2(from, to, ~cor.test(mat1[.x, ], mat1[.y, ], method = "pearson")))
    } else {
      result <- expand_grid(from = rownames(mat1), to = rownames(mat2)) %>%
        mutate(test = map2(from, to, ~cor.test(mat1[.x, ], mat2[.y, ], method = "pearson")))
    }
    result %>%
      mutate(r = map_dbl(test, ~unname(.x$estimate)),
             p = map_dbl(test, "p.value"),
             q_BH = p.adjust(p, method = "BH"),
             type = type,
             selected = abs(r) >= r_cutoff & p < p_cutoff,
             sign = if_else(r >= 0, "Positive", "Negative"),
             weight = abs(r)) %>%
      select(-test)
  }

  gene_gene_all <- edge_tests(gene$log, type = "gene_gene",
                              r_cutoff = 0.92, p_cutoff = 0.005, within = TRUE)
  metabolite_all <- edge_tests(metabolite$log, type = "metabolite_metabolite",
                               r_cutoff = 0.92, p_cutoff = 0.005, within = TRUE)
  cross_all <- edge_tests(gene$log, metabolite$log, type = "gene_metabolite",
                          r_cutoff = 0.97, p_cutoff = 0.005, within = FALSE)
  all_tests <- bind_rows(gene_gene_all, metabolite_all, cross_all)
  edges <- all_tests %>% filter(selected) %>% arrange(type, p, from, to)

  write.csv(cbind(gene_name = rownames(gene$log),
                  setNames(as.data.frame(gene$log), paste0("log2_", common_samples)),
                  setNames(as.data.frame(gene$z), paste0("z_", common_samples))),
            file.path(out_dir, "DEGs_processed_log2_zscore.csv"), row.names = FALSE,
            fileEncoding = "UTF-8")
  write.csv(cbind(Metabolite = rownames(metabolite$log),
                  setNames(as.data.frame(metabolite$log), paste0("log2_", common_samples)),
                  setNames(as.data.frame(metabolite$z), paste0("z_", common_samples))),
            file.path(out_dir, "DEMs_processed_log2_zscore.csv"), row.names = FALSE,
            fileEncoding = "UTF-8")
  write.csv(all_tests, file.path(out_dir, "all_correlation_tests_BH.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(edges, file.path(out_dir, "network_selected_edges.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  nodes <- tibble(
    name = c(rownames(gene$log), rownames(metabolite$log)),
    layer = c(rep("Gene", nrow(gene$log)), rep("Metabolite", nrow(metabolite$log)))
  ) %>%
    left_join(metabolite_labels %>% select(name = Metabolite, metabolite_display = display_name),
              by = "name") %>%
    mutate(display_name = if_else(layer == "Gene", name, metabolite_display)) %>%
    select(-metabolite_display)
  degree_table <- tibble(name = c(edges$from, edges$to)) %>% count(name, name = "degree")
  nodes <- nodes %>% left_join(degree_table, by = "name") %>% mutate(degree = replace_na(degree, 0L))

  layout_layer <- function(df, layer_name, center_y, radius_x, radius_y) {
    n_nodes <- nrow(df)
    angle_slots <- pi / 2 - (seq_len(n_nodes) - 1) * 2 * pi / n_nodes
    facing_risk <- if (layer_name == "Gene") -sin(angle_slots) else sin(angle_slots)
    slot_order <- order(facing_risk, decreasing = TRUE)
    node_order <- order(nchar(df$display_name), -df$degree, df$name)
    assigned_theta <- numeric(n_nodes)
    assigned_slot <- integer(n_nodes)
    for (k in seq_len(n_nodes)) {
      assigned_theta[node_order[k]] <- angle_slots[slot_order[k]]
      assigned_slot[node_order[k]] <- slot_order[k]
    }
    tier_values <- if (layer_name == "Gene") {
      c(0.22, 0.42, 0.62, 0.82, 1.02)
    } else {
      c(0.38, 0.78, 1.18, 1.58, 1.98)
    }
    top_gene_tiers <- c(0.35, 1.25, 2.15, 3.05, 3.95)
    df %>%
      mutate(theta = assigned_theta,
             center_y = center_y,
             radius_x = radius_x,
             radius_y = radius_y,
             tier_offset = if_else(
               layer_name == "Gene" & sin(theta) > 0.45,
               top_gene_tiers[(assigned_slot - 1) %% length(top_gene_tiers) + 1],
               tier_values[(assigned_slot - 1) %% length(tier_values) + 1]
             ),
             x = radius_x * cos(theta),
             y = center_y + radius_y * sin(theta),
             label_x = x + tier_offset * cos(theta),
             label_y = y + tier_offset * sin(theta),
             label_angle = 0,
             label_hjust = case_when(
               abs(sin(theta)) > 0.72 ~ 0.5,
               cos(theta) < 0 ~ 1,
               TRUE ~ 0
             ))
  }
  nodes <- bind_rows(
    layout_layer(filter(nodes, layer == "Gene"), "Gene", 12.5, 18.0, 8.8),
    layout_layer(filter(nodes, layer == "Metabolite"), "Metabolite", -13.5, 23.0, 11.0)
  ) %>%
    mutate(display_plot = mapply(
      function(x, current_layer) paste(strwrap(x, width = ifelse(current_layer == "Gene", 99, 11)),
                                       collapse = "\n"),
      display_name, layer, USE.NAMES = FALSE
    ))
  write.csv(nodes, file.path(out_dir, "network_nodes.csv"), row.names = FALSE,
            fileEncoding = "UTF-8")

  graph <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  p <- ggraph(graph, layout = "manual", x = x, y = y) +
    geom_edge_link(aes(edge_colour = sign, edge_width = weight, edge_alpha = type)) +
    geom_node_point(aes(fill = layer, size = degree, shape = layer),
                    colour = "white", stroke = 0.95) +
    geom_text(data = nodes,
              aes(x = label_x, y = label_y, label = display_plot, angle = label_angle,
                  hjust = label_hjust),
              inherit.aes = FALSE, vjust = 0.5, size = 5.85,
              fontface = "bold", colour = "#111111", check_overlap = FALSE) +
    scale_fill_manual(values = c(Gene = gene_color, Metabolite = metabolite_color)) +
    scale_shape_manual(values = c(Gene = 22, Metabolite = 21), guide = "none") +
    scale_edge_colour_manual(values = c(Positive = positive_color, Negative = negative_color)) +
    scale_edge_width(range = c(0.90, 2.35), guide = "none") +
    scale_edge_alpha_manual(values = c(gene_gene = 0.46,
                                       metabolite_metabolite = 0.46,
                                       gene_metabolite = 0.50), guide = "none") +
    scale_size_continuous(range = c(8.0, 16.0), guide = "none") +
    guides(
      fill = guide_legend(override.aes = list(size = 10)),
      edge_colour = guide_legend(override.aes = list(edge_width = 3, edge_alpha = 1))
    ) +
    coord_cartesian(xlim = c(-38.0, 38.0), ylim = c(-32.5, 27.5), clip = "off") +
    labs(title = "Two-layer gene-metabolite correlation network",
         subtitle = "Pearson correlation after log2(x+1) transformation across paired samples",
         fill = "Node layer", edge_colour = "Interaction") +
    theme_void(base_family = "Arial") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "#666666"),
          legend.position = "right",
          legend.title = element_text(size = 28, face = "bold"),
          legend.text = element_text(size = 26),
          legend.key.height = grid::unit(1.8, "cm"),
          legend.key.width = grid::unit(2.6, "cm"),
          plot.margin = margin(18, 35, 18, 35))

  if (isTRUE(save_plot)) {
    .joint_save(p, file.path(out_dir, "Fig5"), 12, 10)
  }
  invisible(list(plot = p, nodes = nodes, edges = edges, all_tests = all_tests))
}

# Stage 02 - Fig.6A
stage_02_fig6a <- function(c53_csv,all_gene_results_csv,term2gene_file,term2name_file,metabolite_kegg_csv,output_dir){
  .joint_req(c("openxlsx","ggplot2"));out<-.joint_out(output_dir);c53<-.joint_csv(c53_csv);all<-.joint_csv(all_gene_results_csv);met<-.joint_csv(metabolite_kegg_csv);if(nrow(c53)!=53L)stop("Expected 53 overlapping candidate genes in input.")
  k<-.joint_kegg(term2gene_file,term2name_file);ge<-.joint_enrich(c53$geneID,all$geneID,k);detail<-merge(k$g[k$g$term%in%ge$pathway_id&k$g$gene%in%as.character(c53$geneID),],c53[,c("geneID","gene_name","log2FC","Direction")],by.x="gene",by.y="geneID");detail<-merge(detail,k$n,by="term",all.x=TRUE);names(detail)[names(detail)=="term"]<-"pathway_id";names(detail)[names(detail)=="term_name"]<-"pathway_name"
  .joint_pair(ge,file.path(out,"C53_KEGG_enrichment"));.joint_pair(detail,file.path(out,"C53_KEGG_gene_hits"))
  need<-c("pathway_id","pathway_name","Hits","Pvalue");if(!all(need%in%names(met)))stop("Metabolite KEGG table does not match the Stage 05 output format.")
  gx<-data.frame(pathway_id=ge$pathway_id,pathway_name=ge$pathway_name,Omics="Gene",Count=ge$Count,Pvalue=ge$Pvalue,FDR=ge$FDR,FoldEnrichment=ge$FoldEnrichment,stringsAsFactors=FALSE);mx<-data.frame(pathway_id=met$pathway_id,pathway_name=met$pathway_name,Omics="Metabolite",Count=met$Hits,Pvalue=met$Pvalue,FDR=met$FDR,FoldEnrichment=met$EnrichmentFactor,stringsAsFactors=FALSE);mx$pathway_name<-k$n$term_name[match(mx$pathway_id,k$n$term)];joint<-rbind(gx[gx$Count>0&gx$Pvalue<.05,],mx[mx$Count>0&mx$Pvalue<.05,]);joint$Shared<-joint$pathway_id%in%intersect(gx$pathway_id[gx$Pvalue<.05],mx$pathway_id[mx$Pvalue<.05]);joint$minus_log10_p<--log10(pmax(joint$Pvalue,.Machine$double.xmin))
  if(!any(joint$pathway_id=="map04216"&joint$Omics=="Gene")||!any(joint$pathway_id=="map04216"&joint$Omics=="Metabolite"))stop("Ferroptosis (map04216) was not detected in both omics layers for Fig.6A.")
  .joint_pair(joint,file.path(out,"Fig6A_joint_KEGG_source_data"));ord<-aggregate(FoldEnrichment~pathway_name,joint,max);ord<-ord[order(ord$FoldEnrichment,-xtfrm(ord$pathway_name)),];joint$pathway_name<-factor(joint$pathway_name,levels=ord$pathway_name);joint$Omics<-factor(joint$Omics,levels=c("Gene","Metabolite"));joint$fdr<-joint$FDR<.05;ferro<-which(levels(joint$pathway_name)=="Ferroptosis");p<-ggplot2::ggplot(joint,ggplot2::aes(FoldEnrichment,pathway_name))+ggplot2::annotate("rect",xmin=-Inf,xmax=Inf,ymin=ferro-.46,ymax=ferro+.46,fill="#D6A84B",alpha=.12)+ggplot2::geom_point(ggplot2::aes(shape=Omics,size=Count,fill=minus_log10_p,colour=fdr),stroke=1,alpha=.96)+ggplot2::scale_shape_manual(values=c(Gene=21,Metabolite=24),name="Omics layer")+ggplot2::scale_size_continuous(range=c(3.3,8),name="Hit count")+ggplot2::scale_fill_gradientn(colours=c("#1F5A6B","#D6A84B","#C65D3B"),name=expression(-log[10](P)))+ggplot2::scale_colour_manual(values=c(`FALSE`="white",`TRUE`="#263238"),guide="none")+ggplot2::labs(title="Cross-omics KEGG convergence",x="Fold enrichment",y=NULL)+ggplot2::theme_bw(base_size=9,base_family="Arial")+ggplot2::theme(panel.grid.major.x=ggplot2::element_line(colour="#D9DEDF",linetype="dashed"),panel.grid.major.y=ggplot2::element_line(colour="#E7EBEC"),panel.grid.minor=ggplot2::element_blank(),legend.position="right");.joint_save(p,file.path(out,"Fig6A"),8.2,7.2);invisible(list(gene_enrichment=ge,gene_detail=detail,joint=joint))
}

# Stage 03 - Fig.6B
stage_03_fig6b <- function(c53_ferroptosis_gene_hits_csv,dem51_ferroptosis_metabolite_hits_csv,output_dir){
  .joint_req(c("openxlsx","ggplot2"));out<-.joint_out(output_dir);g<-.joint_csv(c53_ferroptosis_gene_hits_csv);m<-.joint_csv(dem51_ferroptosis_metabolite_hits_csv)
  g<-unique(g[g$pathway_id=="map04216",c("gene","gene_name","log2FC"),drop=FALSE]);m<-unique(m[m$pathway_id=="map04216",c("Compound_ID","Metabolite","logFC"),drop=FALSE]);g$gene<-as.character(g$gene);g$Label<-ifelse(g$gene=="122880792","acsl1-like",ifelse(is.na(g$gene_name)|g$gene_name=="",g$gene,g$gene_name));m$Compound_ID<-trimws(as.character(m$Compound_ID));m$Label<-m$Metabolite;m$Label[grepl("Glutathione \\(reduced\\)",m$Metabolite,ignore.case=TRUE)]<-"GSH";m$Label[m$Compound_ID=="Com_2995_pos"|grepl("GluCys|Glutamylcysteine",m$Metabolite,ignore.case=TRUE)]<-"\u03b3-GluCys";m$Label[grepl("Glutathione oxidized",m$Metabolite,ignore.case=TRUE)]<-"GSSG"
  if(nrow(g)!=3L||nrow(m)!=3L)stop(sprintf("Ferroptosis anchor count mismatch: genes=%d (expected 3), metabolites=%d (expected 3); genes=%s; metabolites=%s",nrow(g),nrow(m),paste(g$gene,collapse=","),paste(m$Compound_ID,collapse=",")))
  tab<-rbind(data.frame(ID=g$gene,Label=g$Label,Omics="Gene",log2FC=g$log2FC),data.frame(ID=m$Compound_ID,Label=m$Label,Omics="Metabolite",log2FC=m$logFC));tab$Direction<-ifelse(tab$log2FC>=0,"Higher in F","Lower in F");tab$Omics<-factor(tab$Omics,levels=c("Gene","Metabolite"));tab<-do.call(rbind,lapply(split(tab,tab$Omics),function(z){z<-z[order(z$log2FC),];z$Label<-factor(z$Label,levels=z$Label);z}));.joint_pair(tab,file.path(out,"Fig6B_ferroptosis_factors"));p<-ggplot2::ggplot(tab,ggplot2::aes(log2FC,Label,colour=Direction))+ggplot2::geom_vline(xintercept=0,colour="#8A9499",linewidth=.45)+ggplot2::geom_segment(ggplot2::aes(x=0,xend=log2FC,yend=Label),linewidth=1.15)+ggplot2::geom_point(ggplot2::aes(fill=Direction,shape=Omics),size=5,colour="white",stroke=.65)+ggplot2::facet_grid(Omics~.,scales="free_y",space="free_y")+ggplot2::scale_colour_manual(values=c("Higher in F"="#C65D3B","Lower in F"="#1F5A6B"),guide="none")+ggplot2::scale_fill_manual(values=c("Higher in F"="#C65D3B","Lower in F"="#1F5A6B"),name=NULL)+ggplot2::scale_shape_manual(values=c(Gene=21,Metabolite=24),guide="none")+ggplot2::guides(fill=ggplot2::guide_legend(override.aes=list(shape=21,size=4,colour="white")))+ggplot2::labs(title="Six ferroptosis anchors",x=expression(log[2]*" fold change"),y=NULL)+ggplot2::theme_classic(base_size=10,base_family="Arial")+ggplot2::theme(strip.background=ggplot2::element_rect(fill="#F2F5F4",colour=NA),strip.text.y.right=ggplot2::element_text(angle=-90,face="bold"),legend.position=c(.77,.89));.joint_save(p,file.path(out,"Fig6B"),7.2,6.2);invisible(tab)
}
