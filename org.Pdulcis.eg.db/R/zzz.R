datacache <- new.env(hash=TRUE, parent=emptyenv())

org.Pdulcis.eg <- function() showQCData("org.Pdulcis.eg", datacache)
org.Pdulcis.eg_dbconn <- function() dbconn(datacache)
org.Pdulcis.eg_dbfile <- function() dbfile(datacache)
org.Pdulcis.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.Pdulcis.eg_dbInfo <- function() dbInfo(datacache)

org.Pdulcis.egORGANISM <- "Prunus dulcis"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.Pdulcis.eg.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    db <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, db, envir=ns)
    namespaceExport(ns, dbNewname)
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.Pdulcis.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.Pdulcis.eg_dbconn())
}

