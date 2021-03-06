% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/joinWithGeneSetsFromMatrix.R
\name{joinWithGeneSetsFromExpMatrix}
\alias{joinWithGeneSetsFromExpMatrix}
\title{joinWithGeneSetsFromExpMatrix}
\usage{
joinWithGeneSetsFromExpMatrix(object, Exp_Matrix, coords_df, gene_sets,
  method_gs = "mean", smooth = FALSE, smooth_span = 0.02,
  normalize = F, verbose = TRUE)
}
\description{
Add on SPATA
}
\examples{
#Compare Neftel with ne classification 
### Plot four states 
getGeneSets(object, index="Nef")
gene_sets=c("Neftel_OPClike", "Neftel_NPC_Comb", "Neftel_AClike","Neftel_Mes_Comb", "Neftel_G2.M")
#Get df with states (mean)

Exp_Matrix=normalized_counts(cds)

Cluster <- 
 cds@colData \%>\% 
 as.data.frame() \%>\%
 dplyr::select(barcodes, Age, Patients, Histo) \%>\%
 joinWithGeneSetsFromExpMatrix(object, coords_df = . ,Exp_Matrix=Exp_Matrix, gene_sets=gene_sets, method_gs = "zscore" ,smooth = F)
Cluster$x=Cluster$y=1
Cluster[,5:9]=scale(Cluster[,5:9])


p=plotFourStates(Cluster, states=gene_sets[1:4], color_to=gene_sets[5], pt_alpha = 0.9)
p+scale_color_viridis_c("Cycling Cells",option="inferno", limits = c(1.2, 3), oob = scales::squish)
plotFourStates(Cluster, states=gene_sets[1:4], color_to="Age", pt_alpha = 0.9)

}
\author{
Dieter Henrik Heiland
}
