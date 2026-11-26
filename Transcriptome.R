# Transcriptome analysis

.tx_require <- function(pkgs) { m <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if(length(m)) stop("Missing required R packages: ", paste(m, collapse=", ")) }
.tx_out <- function(path) { if(dir.exists(path) && length(list.files(path, all.files=TRUE, no..=TRUE))) stop("Output directory is not empty: ", path); dir.create(path, recursive=TRUE, showWarnings=FALSE); path }
.tx_tmp <- function(path) { if(!file.exists(path)) stop("Input file not found: ", path); e<-tools::file_ext(path); z<-tempfile(fileext=if(nzchar(e)) paste0(".",e) else ""); if(!file.copy(path,z,overwrite=TRUE)) stop("Failed to read input file: ",path); z }
.tx_read_delim <- function(path) { z<-.tx_tmp(path); on.exit(unlink(z)); read.delim(z, check.names=FALSE, stringsAsFactors=FALSE) }
.tx_read_csv <- function(path) { z<-.tx_tmp(path); on.exit(unlink(z)); read.csv(z, check.names=FALSE, stringsAsFactors=FALSE) }
.tx_read_rds <- function(path) { z<-.tx_tmp(path); on.exit(unlink(z)); readRDS(z) }
.tx_pair <- function(x, stem) { write.csv(x,paste0(stem,".csv"),row.names=FALSE,fileEncoding="UTF-8"); openxlsx::write.xlsx(x,paste0(stem,".xlsx"),overwrite=TRUE) }
.tx_device <- function(path,type,w,h){if(type=="pdf") pdf(path,w,h,useDingbats=FALSE) else tiff(path,w,h,units="in",res=600,compression="lzw")}
.tx_save_plot <- function(p, stem, w, h){ggplot2::ggsave(paste0(stem,".pdf"),p,width=w,height=h,device=grDevices::cairo_pdf);tmp<-tempfile(fileext=".tiff");on.exit(unlink(tmp),add=TRUE);ggplot2::ggsave(tmp,p,width=w,height=h,dpi=600,compression="lzw");if(!file.copy(tmp,paste0(stem,".tiff"),overwrite=TRUE))stop("Failed to write TIFF file: ",stem)}

.tx_pal <- c(teal="#1F5A6B",rust="#C65D3B",gold="#D6A84B",ink="#263238",
             grey="#8A9499",pale="#F2F5F4",grid="#D9DEDF",aqua="#73AEB8")
.tx_theme <- function(base_size=10){
  ggplot2::theme_classic(base_size=base_size,base_family="Arial")+
    ggplot2::theme(text=ggplot2::element_text(colour=.tx_pal[["ink"]]),
      axis.line=ggplot2::element_line(linewidth=.35,colour=.tx_pal[["ink"]]),
      plot.title=ggplot2::element_text(face="bold",hjust=.5),
      plot.subtitle=ggplot2::element_text(colour=.tx_pal[["grey"]],hjust=.5))
}

.tx_resolve_samples <- function(cnt, sample_info=NULL) {
  if(!is.null(sample_info)) {
    if(is.character(sample_info) && length(sample_info)==1L) {
      info <- if(tolower(tools::file_ext(sample_info))=="csv") .tx_read_csv(sample_info) else .tx_read_delim(sample_info)
    } else if(is.data.frame(sample_info)) {
      info <- sample_info
    } else {
      stop("sample_info must be a data.frame or a path to a CSV/TSV file with 'sample' and 'group' columns.")
    }
    if(!all(c("sample","group") %in% names(info))) stop("sample_info must contain 'sample' and 'group' columns.")
    samples <- as.character(info$sample); groups <- as.character(info$group)
    if(anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) stop("sample_info contains missing, empty, or duplicated sample IDs.")
    if(!all(samples %in% names(cnt))) stop("Missing required sample columns in read count table: ", paste(setdiff(samples, names(cnt)), collapse=", "))
  } else {
    cand <- setdiff(names(cnt), "geneID")
    samples <- cand[vapply(cnt[cand], is.numeric, logical(1))]
    if(!length(samples)) stop("No numeric sample columns found in read count table.")
    groups <- ifelse(grepl("^F", samples), "F", ifelse(grepl("^L", samples), "L", NA_character_))
    if(anyNA(groups)) stop("Unable to infer group labels from sample names; use 'F'/'L' prefixes or provide sample_info.")
  }
  if(!identical(sort(unique(groups)), c("F","L"))) stop("Expected exactly two sample groups ('F' and 'L'); found: ", paste(sort(unique(groups)), collapse=", "))
  list(samples=samples, groups=groups)
}

.tx_annotation <- function(path) {
  a <- .tx_read_delim(path)
  if(!all(c("gene_id","gene_name") %in% names(a))) stop("Annotation file is missing gene_id/gene_name columns.")
  a$gene_id <- as.character(a$gene_id); a
}

.tx_kegg_inputs <- function(term2gene_file, term2name_file) {
  t2g <- .tx_read_delim(term2gene_file); t2n <- .tx_read_delim(term2name_file)
  if(!all(c("term","gene") %in% names(t2g)) || !all(c("term","term_name") %in% names(t2n))) stop("KEGG mapping tables are missing required columns.")
  t2g$term <- sub("^ko","map",t2g$term); t2g$gene <- as.character(t2g$gene)
  t2n$term <- sub("^ko","map",t2n$term)
  list(t2g=unique(t2g[,c("term","gene")]),t2n=unique(t2n[,c("term","term_name")]))
}

.tx_enrich <- function(genes, universe, t2g, t2n, min_size=10, max_size=500) {
  universe<-intersect(unique(as.character(universe)),unique(as.character(t2g$gene)))
  fg_all<-unique(as.character(genes)); genes<-intersect(fg_all,universe)
  map<-merge(unique(t2g[,c("term","gene")]),unique(t2n[,c("term","term_name")]),by="term",all.x=TRUE)
  map<-map[map$gene%in%universe,,drop=FALSE]; sets<-split(map$gene,map$term)
  sets<-lapply(sets,unique); sets<-sets[lengths(sets)>=min_size & lengths(sets)<=max_size]
  if(!length(sets)) return(data.frame())
  N<-length(universe); n<-length(genes)
  r<-do.call(rbind,lapply(names(sets),function(id){bg<-sets[[id]];hit<-intersect(bg,genes);M<-length(bg);k<-length(hit);data.frame(pathway_id=id,pathway_name=map$term_name[match(id,map$term)],Bg_Count=M,Count=k,Pvalue=phyper(k-1,M,N-M,n,lower.tail=FALSE),GeneRatio=k/n,BgRatio=M/N,FoldEnrichment=(k/n)/(M/N),RichFactor=k/M,hit_geneIDs=paste(sort(hit),collapse="/"),Foreground_Total=length(fg_all),Foreground_Annotated=n,Background_Annotated_Tested=N,Min_GS_Size=min_size,Max_GS_Size=max_size,stringsAsFactors=FALSE)}))
  r$FDR<-p.adjust(r$Pvalue,"BH");r[order(r$Pvalue,-r$FoldEnrichment),]
}

.tx_bubble <- function(tab, title, top_n=20) {
  if(!nrow(tab)) stop("No enrichment results available for plotting.")
  x <- tab[order(tab$pvalue),,drop=FALSE]; x <- head(x,top_n); x$pathway_name <- factor(x$pathway_name,levels=rev(x$pathway_name)); x$minus_log10_p <- -log10(pmax(x$pvalue,.Machine$double.xmin)); x$GeneRatio <- vapply(strsplit(as.character(x$GeneRatio),"/",fixed=TRUE),function(z)as.numeric(z[1])/as.numeric(z[2]),numeric(1))
  ggplot2::ggplot(x,ggplot2::aes(x=GeneRatio,y=pathway_name,size=Count,color=minus_log10_p))+ggplot2::geom_point()+ggplot2::scale_color_gradient(low="#3D9F3C",high="#D73027")+ggplot2::labs(x="Gene ratio",y=NULL,color="-log10(P)",size="Count",title=title)+ggplot2::theme_bw(base_size=11)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank())
}

# Stage 01 - Fig.2A
stage_01_fig2a <- function(read_count_file, annotation_file, output_dir, sample_info=NULL) {
  .tx_require(c("DESeq2","openxlsx","ggplot2","ggrepel")); out<-.tx_out(output_dir)
  cnt <- .tx_read_delim(read_count_file)
  if(!"geneID" %in% names(cnt)) stop("Read count table is missing the 'geneID' column.")
  if(anyDuplicated(cnt$geneID)) stop("Duplicate geneID entries detected.")
  meta <- .tx_resolve_samples(cnt, sample_info=sample_info)
  samples <- meta$samples
  mat<-round(as.matrix(cnt[,samples,drop=FALSE])); storage.mode(mat)<-"integer"; rownames(mat)<-as.character(cnt$geneID)
  coldata<-data.frame(group=factor(meta$groups,levels=c("L","F")),row.names=samples)
  dds<-DESeq2::DESeqDataSetFromMatrix(mat,coldata,~group); dds<-DESeq2::DESeq(dds,quiet=TRUE)
  rr<-as.data.frame(DESeq2::results(dds,contrast=c("group","F","L"))); rr$geneID<-rownames(rr)
  ann<-.tx_annotation(annotation_file); rr<-merge(rr,ann[,c("gene_id","gene_name")],by.x="geneID",by.y="gene_id",all.x=TRUE,sort=FALSE)
  rr$Direction<-ifelse(!is.na(rr$padj)&rr$padj<.05&rr$log2FoldChange>=1,"Up_in_F",ifelse(!is.na(rr$padj)&rr$padj<.05&rr$log2FoldChange<=-1,"Up_in_L","Not_significant"))
  rr<-rr[!is.na(rr$pvalue),,drop=FALSE]
  rr<-rr[,c("geneID","gene_name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","Direction")]
  names(rr)[names(rr)=="log2FoldChange"]<-"log2FC"
  deg<-rr[rr$Direction!="Not_significant",,drop=FALSE]
  if(nrow(rr)!=22513L || nrow(deg)!=79L || sum(deg$Direction=="Up_in_F")!=43L || sum(deg$Direction=="Up_in_L")!=36L) stop(sprintf("DESeq2 validation failed: tested=%d, DEG=%d, Up_in_F=%d, Up_in_L=%d.",nrow(rr),nrow(deg),sum(deg$Direction=="Up_in_F"),sum(deg$Direction=="Up_in_L")))
  .tx_pair(rr,file.path(out,"all_DESeq2_results")); .tx_pair(deg,file.path(out,"DEG79"))
  plot_rr<-rr;plot_rr$neglog10P<--log10(pmax(plot_rr$pvalue,.Machine$double.xmin))
  lab<-do.call(rbind,lapply(split(deg,deg$Direction),function(z)head(z[order(z$padj,z$pvalue,-abs(z$log2FC),z$geneID),],3)))
  lab$label<-ifelse(is.na(lab$gene_name)|lab$gene_name=="",lab$geneID,lab$gene_name);lab$neglog10P<--log10(pmax(lab$pvalue,.Machine$double.xmin))
  p<-ggplot2::ggplot(plot_rr,ggplot2::aes(log2FC,neglog10P))+
    ggplot2::geom_vline(xintercept=c(-1,1),linetype="dashed",linewidth=.3,colour=.tx_pal[["grey"]])+
    ggplot2::geom_hline(yintercept=-log10(.05),linetype="dashed",linewidth=.3,colour=.tx_pal[["grey"]])+
    ggplot2::geom_point(data=plot_rr[plot_rr$Direction=="Not_significant",],colour=.tx_pal[["grid"]],size=.55,alpha=.72)+
    ggplot2::geom_point(data=plot_rr[plot_rr$Direction!="Not_significant",],ggplot2::aes(colour=Direction),size=1.2,alpha=.9)+
    ggrepel::geom_text_repel(data=lab,ggplot2::aes(label=label),family="Arial",size=2.3,seed=20260716,max.overlaps=Inf,box.padding=.42,point.padding=.22,force=3,min.segment.length=0,segment.size=.22,colour=.tx_pal[["ink"]])+
    ggplot2::scale_colour_manual(values=c(Up_in_L=.tx_pal[["teal"]],Up_in_F=.tx_pal[["rust"]]),labels=c(Up_in_L="Down",Up_in_F="Up"),name=NULL)+
    ggplot2::labs(title="79-DEG landscape",subtitle="43 increased, 36 decreased\nAdjusted P < 0.05 and |log2FC| >= 1",x=expression(log[2]*" fold change (F versus L)"),y=expression(-log[10](P)))+.tx_theme()+ggplot2::theme(legend.position="bottom")
  .tx_save_plot(p,file.path(out,"Fig2A"),7.5,6.0); invisible(list(all=rr,deg79=deg))
}

# Stage 02 - 53-gene intersection
stage_02_c53 <- function(deg79_csv, wgcna_rds, output_dir) {
  .tx_require(c("openxlsx","WGCNA")); out<-.tx_out(output_dir); deg<-.tx_read_csv(deg79_csv); b<-.tx_read_rds(wgcna_rds)
  if(!all(c("geneID","Direction")%in%names(deg)) || nrow(deg)!=79L) stop("Input is not the 79-DEG table from Stage 01.")
  if(is.null(b$gene_ids)||is.null(b$module_colors)||is.null(b$module_mapping)) stop("Input must be the unified WGCNA Stage 01 RDS bundle.")
  mm<-b$module_mapping; if(mm$DisplayModule[match("brown",mm$ModuleColor)]!="MEnine") stop("Inconsistent mapping between brown module and MEnine.")
  brown<-as.character(b$gene_ids[as.character(b$module_colors)=="brown"]); c53<-deg[as.character(deg$geneID)%in%brown,,drop=FALSE]
  if(nrow(c53)!=53L) stop("Expected 53 overlapping genes between brown/MEnine and 79 DEGs, got: ",nrow(c53))
  .tx_pair(c53,file.path(out,"C53_intersection_genes")); invisible(c53)
}

.tx_go_map <- function(annotation_file) {
  .tx_require(c("AnnotationDbi","GO.db"));a<-.tx_annotation(annotation_file); if(!"GO_pairs"%in%names(a)) stop("Annotation file is missing the GO_pairs column.")
  rows<-strsplit(ifelse(is.na(a$GO_pairs),"",a$GO_pairs),";",fixed=TRUE);lens<-lengths(rows);x<-data.frame(gene=rep(a$gene_id,lens),pair=trimws(unlist(rows)),stringsAsFactors=FALSE);x<-x[grepl("GO:[0-9]+",x$pair),];x$GO_ID<-sub(".*?(GO:[0-9]+).*","\\1",x$pair);x$pair_name<-sub("^GO:[0-9]+\\|","",x$pair);x<-x[!duplicated(x[,c("GO_ID","gene")]),]
  meta<-AnnotationDbi::select(GO.db::GO.db,keys=unique(x$GO_ID),columns=c("TERM","ONTOLOGY"),keytype="GOID");meta<-meta[!duplicated(meta$GOID),]
  x$GO_Name<-meta$TERM[match(x$GO_ID,meta$GOID)];x$GO_Name[is.na(x$GO_Name)]<-x$pair_name[is.na(x$GO_Name)];code<-meta$ONTOLOGY[match(x$GO_ID,meta$GOID)];x$Ontology<-c(BP="Biological process",CC="Cellular component",MF="Molecular function")[code];x[,c("GO_ID","GO_Name","Ontology","gene")]
}

# Stage 03 - Fig.2B
stage_03_fig2b <- function(deg79_csv, all_results_csv, annotation_file, output_dir) {
  .tx_require(c("openxlsx","ggplot2","AnnotationDbi","GO.db"));out<-.tx_out(output_dir);deg<-.tx_read_csv(deg79_csv);all<-.tx_read_csv(all_results_csv);g<-.tx_go_map(annotation_file)
  universe<-intersect(as.character(all$geneID),unique(g$gene));fg<-intersect(as.character(deg$geneID),universe);gm<-g[g$gene%in%universe & !is.na(g$Ontology),];sets<-split(gm$gene,gm$GO_ID);sets<-lapply(sets,unique);sets<-sets[lengths(sets)>=10 & lengths(sets)<=500];N<-length(universe);n<-length(fg)
  r<-do.call(rbind,lapply(names(sets),function(id){bg<-sets[[id]];hit<-intersect(bg,fg);M<-length(bg);k<-length(hit);q<-gm[match(id,gm$GO_ID),];data.frame(GO_ID=id,GO_Name=q$GO_Name,Ontology=q$Ontology,Bg_Count=M,Count=k,Pvalue=phyper(k-1,M,N-M,n,lower.tail=FALSE),GeneRatio=k/n,BgRatio=M/N,FoldEnrichment=(k/n)/(M/N),RichFactor=k/M,hit_geneIDs=paste(sort(hit),collapse="/"),Foreground_Total=nrow(deg),Foreground_Annotated=n,Background_Annotated_Tested=N,Min_GS_Size=10,Max_GS_Size=500,stringsAsFactors=FALSE)}));r$FDR<-p.adjust(r$Pvalue,"BH");r<-r[order(r$Pvalue,-r$FoldEnrichment),]
  detail<-merge(gm[gm$GO_ID%in%r$GO_ID & gm$gene%in%fg,],deg[,c("geneID","gene_name","log2FC","Direction")],by.x="gene",by.y="geneID")
  .tx_pair(r,file.path(out,"Fig2B_GO_enrichment"));.tx_pair(detail,file.path(out,"Fig2B_GO_gene_hits"))
  nom<-r[r$Pvalue<.05 & !is.na(r$Ontology),];bp0<-nom[nom$Ontology=="Biological process",];bp<-head(bp0[order(bp0$Pvalue),],8);cc<-nom[nom$Ontology=="Cellular component",];mf0<-nom[nom$Ontology=="Molecular function",];mf<-head(mf0[order(mf0$Pvalue),],8);plot<-rbind(bp,cc,mf);plot<-plot[order(match(plot$Ontology,c("Biological process","Cellular component","Molecular function")),plot$FoldEnrichment),];plot$term_label<-vapply(plot$GO_Name,function(z)paste(strwrap(z,39),collapse="\n"),character(1));plot$term_key<-factor(paste(plot$term_label,plot$GO_ID),levels=unique(paste(plot$term_label,plot$GO_ID)))
  p<-ggplot2::ggplot(plot,ggplot2::aes(FoldEnrichment,term_key))+ggplot2::geom_point(ggplot2::aes(size=Count,colour=-log10(Pvalue)),alpha=.9)+ggplot2::facet_grid(Ontology~.,scales="free_y",space="free_y",labeller=ggplot2::labeller(Ontology=c("Biological process"="BP","Cellular component"="CC","Molecular function"="MF")))+ggplot2::scale_y_discrete(labels=function(z)sub(" GO:[0-9]+$","",z))+ggplot2::scale_colour_gradient(low=.tx_pal[["teal"]],high=.tx_pal[["rust"]],name=expression(-log[10](P)))+ggplot2::scale_size_continuous(range=c(2.2,6),name="Genes")+ggplot2::labs(title="GO enrichment",x="Fold enrichment",y=NULL)+ggplot2::theme_bw(base_size=8,base_family="Arial")+ggplot2::theme(axis.text.y=ggplot2::element_text(size=6.7),panel.grid.major.x=ggplot2::element_line(colour=.tx_pal[["grid"]],linetype="dashed"),panel.grid.major.y=ggplot2::element_blank(),strip.background=ggplot2::element_rect(fill=.tx_pal[["pale"]]),strip.text.y.right=ggplot2::element_text(angle=0,face="bold"),legend.position="bottom")
  .tx_save_plot(p,file.path(out,"Fig2B"),8.2,6.8);invisible(r)
}

.tx_parse_hierarchy <- function(path) {
  z<-.tx_tmp(path); on.exit(unlink(z)); lines<-readLines(z,warn=FALSE); A<-B<-NA_character_; out<-list()
  for(s in lines){if(grepl("^A",s)){A<-trimws(sub("^A","",s));B<-NA_character_} else if(grepl("^B",s)){B<-trimws(sub("^B","",s))} else if(grepl("^C",s)){m<-regexec("^C\\s+([0-9]{5})\\s+(.+)$",s);q<-regmatches(s,m)[[1]];if(length(q)==3) out[[length(out)+1]]<-data.frame(pathway_id=paste0("map",q[2]),pathway_name=q[3],Category=A,Subcategory=B)}}
  do.call(rbind,out)
}

# Stage 04 - Fig.2C
stage_04_fig2c <- function(deg79_csv, term2gene_file, term2name_file, hierarchy_file, output_dir) {
  .tx_require(c("openxlsx","ggplot2")); out<-.tx_out(output_dir); deg<-.tx_read_csv(deg79_csv); k<-.tx_kegg_inputs(term2gene_file,term2name_file); h<-.tx_parse_hierarchy(hierarchy_file)
  hits<-unique(k$t2g[k$t2g$gene%in%as.character(deg$geneID),]); tab<-aggregate(gene~term,hits,function(x)length(unique(x))); names(tab)<-c("pathway_id","Count"); tab<-merge(tab,h,by="pathway_id"); tab<-tab[tab$Count>=2,,drop=FALSE]
  class_order<-c("Cellular Processes","Environmental Information Processing","Human Diseases","Metabolism","Organismal Systems");tab<-tab[tab$Category%in%class_order & tab$pathway_id!="map01100",];bg_sets<-split(k$t2g$gene,k$t2g$term);tab$Bg_Count<-lengths(bg_sets[tab$pathway_id]);N<-length(unique(k$t2g$gene));n<-length(intersect(as.character(deg$geneID),unique(k$t2g$gene)));tab$Pvalue<-mapply(function(M,z)phyper(z-1,M,N-M,n,lower.tail=FALSE),tab$Bg_Count,tab$Count);tab$Category<-factor(tab$Category,levels=class_order);tab<-tab[order(tab$Category,tab$Count,tab$Pvalue),];tab$pathway_key<-factor(paste(tab$Category,tab$pathway_name,sep="___"),levels=unique(paste(tab$Category,tab$pathway_name,sep="___")))
  .tx_pair(tab,file.path(out,"Fig2C_KEGG_classification"));cols<-c("Cellular Processes"=.tx_pal[["teal"]],"Environmental Information Processing"=.tx_pal[["rust"]],"Human Diseases"=.tx_pal[["grey"]],"Metabolism"=.tx_pal[["aqua"]],"Organismal Systems"=.tx_pal[["gold"]]);p<-ggplot2::ggplot(tab,ggplot2::aes(Count,pathway_key,fill=Category))+ggplot2::geom_col(width=.72)+ggplot2::geom_text(ggplot2::aes(label=Count),hjust=-.2,size=2.7)+ggplot2::facet_grid(Category~.,scales="free_y",space="free_y")+ggplot2::scale_y_discrete(labels=function(z)sub("^.*___","",z))+ggplot2::scale_fill_manual(values=cols,guide="none")+ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(0,.15)))+ggplot2::labs(title="KEGG classification of 79 differentially expressed genes",x="Number of DEGs",y=NULL)+ggplot2::theme_bw(base_size=8,base_family="Arial")+ggplot2::theme(panel.grid.major.x=ggplot2::element_line(colour=.tx_pal[["grid"]],linetype="dashed"),panel.grid.major.y=ggplot2::element_blank(),strip.background=ggplot2::element_rect(fill=.tx_pal[["pale"]]),strip.text.y.right=ggplot2::element_text(angle=0,face="bold"));.tx_save_plot(p,file.path(out,"Fig2C"),9.5,8);invisible(tab)
}

# Stage 05 - Fig.2D
stage_05_fig2d <- function(deg79_csv, all_results_csv, term2gene_file, term2name_file, output_dir) {
  .tx_require(c("openxlsx","ggplot2"));out<-.tx_out(output_dir);deg<-.tx_read_csv(deg79_csv);all<-.tx_read_csv(all_results_csv);k<-.tx_kegg_inputs(term2gene_file,term2name_file)
  r<-.tx_enrich(deg$geneID,all$geneID,k$t2g,k$t2n,10,500);detail<-merge(k$t2g[k$t2g$term%in%r$pathway_id & k$t2g$gene%in%deg$geneID,],deg[,c("geneID","gene_name","log2FC","Direction")],by.x="gene",by.y="geneID");detail<-merge(detail,k$t2n,by="term",all.x=TRUE)
  .tx_pair(r,file.path(out,"Fig2D_KEGG_enrichment"));.tx_pair(detail,file.path(out,"Fig2D_KEGG_gene_hits"));plot<-head(r[r$Count>0 & r$Pvalue<.05,],12);plot$pathway_name<-factor(plot$pathway_name,levels=plot$pathway_name[order(plot$Pvalue,decreasing=TRUE)]);p<-ggplot2::ggplot(plot,ggplot2::aes(FoldEnrichment,pathway_name))+ggplot2::geom_point(ggplot2::aes(size=Count,fill=-log10(Pvalue)),shape=21,colour="white",stroke=.65)+ggplot2::scale_fill_gradientn(colours=c(.tx_pal[["teal"]],.tx_pal[["gold"]],.tx_pal[["rust"]]),name=expression(-log[10](P)))+ggplot2::scale_size_continuous(range=c(3,7),name="Genes")+ggplot2::labs(title="79-DEG KEGG enrichment",x="Fold enrichment",y=NULL)+ggplot2::theme_bw(base_size=8,base_family="Arial")+ggplot2::theme(panel.grid.major.x=ggplot2::element_line(colour=.tx_pal[["grid"]],linetype="dashed"),panel.grid.major.y=ggplot2::element_line(colour="#E7EBEC"),panel.grid.minor=ggplot2::element_blank(),legend.position="right");.tx_save_plot(p,file.path(out,"Fig2D"),7.4,7.2);invisible(r)
}
