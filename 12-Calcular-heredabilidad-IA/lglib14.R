
#HOME   <<- Sys.getenv ("MULTIGWAS_HOME")
#.libPaths (paste0(HOME, "/opt/Rlibs"))

args = commandArgs(trailingOnly = TRUE)
options (width=300)
#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
view <- function (data, n=5,m=6, suffix="") {
	filename = deparse (substitute (data))
	name = paste (deparse (substitute (data)),":  ")
	if (is.null (dim (data))) {
		dimensions = paste (length (data))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (length (data) < 6) n = length(data)
		print (data[1:n])
	}else {
		dimensions = paste0 (unlist (dim (data)),sep=c(" x ",""))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (n==0 | nrow (data) < 5) n = nrow(data)
		if (m==0 | ncol (data) < 6) m = ncol(data)
		print (data[1:n,1:m])
	}
	write.csv (data, paste0("x-", filename, suffix, ".csv"))
}
viewx <- function (...){
	view (...)
	quit()
}

#-------------------------------------------------------------
# Add label to filename and new extension (optional)
#-------------------------------------------------------------
addLabel <- function (filename, label, newExt=NULL)  {
	nameext = strsplit (filename, split="[.]")
	name    = nameext [[1]][1] 
	if (is.null (newExt))
		ext     = nameext [[1]][2] 
	else
		ext     = newExt
	newName = paste0 (nameext [[1]][1], "-", label, ".", ext )
	return (newName)
}
hd=view


#----------------------------------------------------------
# Utility for create dir, if it exists, it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", dirname (newDir), basename (newDir))
			if (dir.exists (oldDir) == T) checkOldDir (oldDir)
			file.rename (newDir, oldDir)
		}
	}
	checkOldDir (newDir)
	dir.create (sprintf (newDir))
}
#--------------------------------------------------------------
# Convert a dataframe to list of vectors
#--------------------------------------------------------------
dataFrameToList <- function (dataFrame) {
	dataList = split(dataFrame,seq(nrow(dataFrame)))
	dataList = setNames(dataList, rownames(dataFrame))
	return (dataList)
}

hst <- function (n=NULL) {
	filename = "/tmp/hst.txt"

	savehistory (filename)
	#system ("cat /tmp/hst.txt|xclip -sel cli")
	system ("tail -2 /tmp/hst.txt|head -1|xclip -sel cli")
}

