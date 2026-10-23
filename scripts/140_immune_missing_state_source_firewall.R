args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('usage: Rscript 140_immune_missing_state_source_firewall.R <DE.rds> <outdir>')
infile <- args[[1]]
outdir <- args[[2]]

stim_map <- list(
  flg22   = list(top='DE.flg22',  suffix='flg22'),
  Pep1    = list(top='DE.Pep1',   suffix='Pep1'),
  Chitin  = list(top='DE.Chitin', suffix='Chitin'),
  SCOOP10 = list(top='DE.SCOOP10',suffix='SCOOP10'),
  SCREW1  = list(top='DE.SCREW1', suffix='SCREW1')
)

source_prefix <- c(
  'Rootcap_Columella',
  'Meristem_Epidermis','Meristem_Endodermis','Meristem_Stele',
  'Elongation_Epidermis','Elongation_Cortex','Elongation_Stele',
  'Maturation_Epidermis','Maturation_Cortex','Maturation_Endodermis'
)
target_prefix <- c('Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele')

x <- readRDS(infile)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
manifest <- list(source_states=source_prefix, sealed_target_states=target_prefix, stimuli=names(stim_map), outputs=list())

write_one <- function(tab, path) {
  req <- c('avg_logFC')
  if (!all(req %in% colnames(tab))) stop(paste('missing avg_logFC in', path))
  z <- data.frame(gene=rownames(tab), tab, check.names=FALSE, stringsAsFactors=FALSE)
  con <- gzfile(path, open='wt')
  on.exit(close(con), add=TRUE)
  write.table(z, con, sep='\t', row.names=FALSE, quote=FALSE)
}

for (stim in names(stim_map)) {
  spec <- stim_map[[stim]]
  if (!(spec$top %in% names(x))) stop(paste('missing top list', spec$top))
  d <- x[[spec$top]]
  expected_all <- c(paste0(source_prefix, '_', spec$suffix), paste0(target_prefix, '_', spec$suffix))
  if (!all(expected_all %in% names(d))) {
    miss <- setdiff(expected_all, names(d))
    stop(paste('missing expected states for', stim, paste(miss, collapse=',')))
  }
  sd <- file.path(outdir, stim)
  dir.create(sd, recursive=TRUE, showWarnings=FALSE)
  rows <- list()
  for (p in source_prefix) {
    nm <- paste0(p, '_', spec$suffix)
    tab <- d[[nm]]
    op <- file.path(sd, paste0(p, '.tsv.gz'))
    write_one(tab, op)
    rows[[p]] <- nrow(tab)
  }
  manifest$outputs[[stim]] <- rows
  cat('SOURCE_ONLY', stim, 'states=', length(source_prefix), 'rows=', paste(unlist(rows), collapse=','), '\n')
}

saveRDS(manifest, file.path(outdir, 'SOURCE_FIREWALL_MANIFEST.rds'))
cat('IMMUNE_SOURCE_FIREWALL_SUCCESS\n')
cat('TARGET_TABLES_EXPORTED=0\n')
