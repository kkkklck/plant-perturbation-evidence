args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('usage: Rscript 142_immune_missing_state_target_firewall.R <DE.rds> <outdir>')
infile <- args[[1]]
outdir <- args[[2]]

stim_map <- list(
  flg22   = list(top='DE.flg22',  suffix='flg22'),
  Pep1    = list(top='DE.Pep1',   suffix='Pep1'),
  Chitin  = list(top='DE.Chitin', suffix='Chitin'),
  SCOOP10 = list(top='DE.SCOOP10',suffix='SCOOP10'),
  SCREW1  = list(top='DE.SCREW1', suffix='SCREW1')
)

target_prefix <- c('Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele')
source_prefix <- c(
  'Rootcap_Columella',
  'Meristem_Epidermis','Meristem_Endodermis','Meristem_Stele',
  'Elongation_Epidermis','Elongation_Cortex','Elongation_Stele',
  'Maturation_Epidermis','Maturation_Cortex','Maturation_Endodermis'
)

x <- readRDS(infile)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
manifest <- list(target_states=target_prefix, source_states_not_exported=source_prefix, stimuli=names(stim_map), outputs=list())

write_one <- function(tab, path) {
  if (!('avg_logFC' %in% colnames(tab))) stop(paste('missing avg_logFC in', path))
  z <- data.frame(gene=rownames(tab), tab, check.names=FALSE, stringsAsFactors=FALSE)
  con <- gzfile(path, open='wt')
  on.exit(close(con), add=TRUE)
  write.table(z, con, sep='\t', row.names=FALSE, quote=FALSE)
}

for (stim in names(stim_map)) {
  spec <- stim_map[[stim]]
  if (!(spec$top %in% names(x))) stop(paste('missing top list', spec$top))
  d <- x[[spec$top]]
  expected <- paste0(target_prefix, '_', spec$suffix)
  if (!all(expected %in% names(d))) {
    miss <- setdiff(expected, names(d))
    stop(paste('missing target states for', stim, paste(miss, collapse=',')))
  }
  td <- file.path(outdir, stim)
  dir.create(td, recursive=TRUE, showWarnings=FALSE)
  rows <- list()
  for (p in target_prefix) {
    nm <- paste0(p, '_', spec$suffix)
    tab <- d[[nm]]
    op <- file.path(td, paste0(p, '.tsv.gz'))
    write_one(tab, op)
    rows[[p]] <- nrow(tab)
  }
  manifest$outputs[[stim]] <- rows
  cat('TARGET_ONLY', stim, 'states=', length(target_prefix), 'rows=', paste(unlist(rows), collapse=','), '\n')
}

saveRDS(manifest, file.path(outdir, 'TARGET_FIREWALL_MANIFEST.rds'))
cat('IMMUNE_TARGET_FIREWALL_SUCCESS\n')
cat('SOURCE_TABLES_EXPORTED=0\n')
