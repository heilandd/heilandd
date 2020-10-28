#' @title joinWithGeneSetsFromExpMatrix
#' @author Dieter Henrik Heiland
#' @description Add on SPATA
#' @inherit 
#' @return 
#' @examples 
#'#Compare Neftel with ne classification 
#'### Plot four states 
#'getGeneSets(object, index="Nef")
#'gene_sets=c("Neftel_OPClike", "Neftel_NPC_Comb", "Neftel_AClike","Neftel_Mes_Comb", "Neftel_G2.M")
#'#Get df with states (mean)
#'
#'Exp_Matrix=normalized_counts(cds)
#'
#'Cluster <- 
#'  cds@colData %>% 
#'  as.data.frame() %>%
#'  dplyr::select(barcodes, Age, Patients, Histo) %>%
#'  joinWithGeneSetsFromExpMatrix(object, coords_df = . ,Exp_Matrix=Exp_Matrix, gene_sets=gene_sets, method_gs = "zscore" ,smooth = F)
#'Cluster$x=Cluster$y=1
#'Cluster[,5:9]=scale(Cluster[,5:9])
#'
#'
#'p=plotFourStates(Cluster, states=gene_sets[1:4], color_to=gene_sets[5], pt_alpha = 0.9)
#'p+scale_color_viridis_c("Cycling Cells",option="inferno", limits = c(1.2, 3), oob = scales::squish)
#'plotFourStates(Cluster, states=gene_sets[1:4], color_to="Age", pt_alpha = 0.9)
#' 
#' @export
#' 





joinWithGeneSetsFromExpMatrix <- function(object,
                                          Exp_Matrix,
                                          coords_df,
                                          gene_sets,
                                          method_gs = "mean",
                                          smooth = FALSE,
                                          smooth_span = 0.02,
                                          normalize = F,
                                          verbose = TRUE){
  
  # 1. Control --------------------------------------------------------------
  
  # lazy check
  
  #check_object(object)
  #check_coords_df(coords_df)
  #check_smooth(df = coords_df, smooth = smooth, smooth_span = smooth_span)
  #check_method(method_gs = method_gs)
  
  
  # adjusting check
  gene_sets <- SPATA::check_gene_sets(object, gene_sets = gene_sets)
  
  # -----
  
  # 2. Extract gene set data and join with coords_df ------------------------
  
  rna_assay <- base::as.matrix(Exp_Matrix)
  gene_set_df <- object@used_genesets
  joined_df <- coords_df
  
  for(i in base::seq_along(gene_sets)){
    
    # get gene names of specified gene sets
    genes <-
      gene_set_df %>%
      dplyr::filter(ont %in% gene_sets[i]) %>%
      dplyr::filter(gene %in% base::rownames(rna_assay)) %>%
      dplyr::pull(gene)
    
    # apply specified method to handle gene sets
    if(method_gs == "mean"){
      
      geneset_vls <-
        base::colMeans(rna_assay[genes, ]) %>%
        base::as.data.frame() %>%
        magrittr::set_colnames(value = gene_sets[i]) %>%
        tibble::rownames_to_column(var = "barcodes")
      
      if(verbose){
        
        base::message(stringr::str_c(
          "Calculating expression score for gene set ",
          "(",i, "/", base::length(gene_sets), ")", "  '",
          gene_sets[i],
          "' according to method: '",
          method_gs,
          "'.",
          sep = ""))
        
      }
      
      
    } else if(method_gs %in% c("gsva", "ssgsea", "zscore", "plage")) {
      
      if(verbose){
        
        base::message(stringr::str_c(
          "Calculating expression score for gene set ",
          "(",i, "/", base::length(gene_sets), ")", "  '",
          gene_sets[i],
          "' according to method: '",
          method_gs,
          "'. This might take a few moments.",
          sep = ""))
        
      }
      
      geneset_vls <-
        GSVA::gsva(expr = rna_assay[genes,],
                   gset.idx.list = gene_set_df,
                   mx.diff = 1,
                   parallel.sz = 2,
                   method = method_gs,
                   verbose = FALSE) %>%
        t() %>%
        as.data.frame() %>%
        magrittr::set_colnames(value = gene_sets[i]) %>%
        tibble::rownames_to_column(var = "barcodes")
      
    } else {
      
      stop("Please enter valid method to handle gene sets.")
      
    }
    
    # gradually add gene_set columns to joined_df
    joined_df <-
      dplyr::left_join(x = joined_df, y = geneset_vls, by = "barcodes")
    
  }
  
  # -----
  
  
  # 3. Smooth and normalize if specified ------------------------------------
  
  if(base::isTRUE(smooth)){
    
    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_smooth,
                      coords_df = joined_df,
                      verbose = verbose,
                      smooth_span = smooth_span,
                      aspect = "gene set",
                      subset = gene_sets)
    
  }
  
  if(base::isTRUE(normalize)){
    
    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_normalize_imap,
                      aspect = "Gene set",
                      verbose = verbose,
                      subset = gene_sets)
    
  }
  
  # -----
  
  base::return(joined_df)
  
}