#' @title Transform SPATA object into a python list (json-file)
#' @author Dieter Henrik Heiland
#' @description Test
#' @param object SPATA object 
#' @param of_sample Name of the sample of the Spata Object
#' @param tissue_pos File of the possition in the tissue
#' @inherit verbose params
#' @return No return but saves a json file to workdir
#' @examples 
#' 
#' #Get Spata Object
#' library(SPATA)
#'input_paths = c("~/Desktop/SpatialTranscriptomics/Visium/Visium/275_T")
#'sample_names = c("275_T")
#'
#'spata_obj <- initiateSpataObject_10X(input_paths=input_paths, 
#'                                     sample_names=sample_names,
#'                                     object_name="Obj_275_T",
#'                                     output_path = getwd())
#'                                     
#'tissue_pos=read.csv(paste0(input_paths,"/outs/spatial/tissue_positions_list.csv"), header=F)
#'names(tissue_pos)=c("barcodes", "tissue", "x_rank", "y_rank", "x", "y")
#'tissue_pos$barcodes=paste0(tissue_pos$barcodes, "_", sample_names) 
#'
#'#Tranform data to list 
#'Transform_gene_array(object=spata_obj, of_sample=sample_names, tissue_pos)
#'
#' 
#' @export
#' 

Transform_gene_array=function(object, of_sample, tissue_pos){
  
  dim=base::dim(SPATA::exprMtr(object, of_sample))
  genes= SPATA::getGenes(object)
  
  base::message("... get joint genes ...")
  
  step1 <-
    SPATA::joinWithGenes(object = spata_obj,
                  spata_df = SPATA::getCoordinates(spata_obj),
                  genes = SPATA::getGenes(object), verbose = F) %>% 
    dplyr::left_join(., tissue_pos, by="barcodes")
  
  base::message("... create matrix ...")
  
  
  message("... compute array...take some time...")
  gene_array_list=pbapply::pblapply( 1:dim[1], function(i){
    gene_inp=genes[i]
    step2<-step1 %>% dplyr::select(x_rank,y_rank, !!sym(gene_inp))
    #new Changed
    mat1=tidyr::pivot_wider(step2, names_from = y_rank, values_from = !!sym(gene_inp)) %>% as.data.frame()
    base::rownames(mat1)=mat1$x_rank ; mat1$x_rank = NULL ; mat1=base::as.matrix(mat1)
    mat1=mat1[base::order(as.numeric(rownames(mat1))), base::order(as.numeric(colnames(mat1)))]
    base::return(mat1)
  })
  
  mase::names(gene_array_list)=genes
  
  base::message("... write output...")
  jsonlite::write_json(gene_array_list, paste0(of_sample,".json"))
  
}





