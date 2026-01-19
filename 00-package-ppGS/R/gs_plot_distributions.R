#' Create a plot for color trait densitydistributioons.
#'
#' @param filename of table with values of HCL trait.
#' @import ggplot2 
#' @export
gs_plot_distributions <- function (HCLFilename) {
	datos = read.csv (HCLFilename)
	datos = datos [,-1] # Remove "Register" column
	datos = reshape::melt (datos[,]);
	colnames(datos) = c("Trait", "Value");
	datos = datos [!is.na(datos$Value),]

	ggplot (data = datos, aes (x=Value)) +
			geom_density() +
			theme(text = element_text(size = 7)) +
			facet_wrap (~Trait, scale="free", ncol=6)

  	ggsave ("gs-plot-distributions-colorTraits.pdf", width=10)
}
#graficarDistribuciones (datos)
