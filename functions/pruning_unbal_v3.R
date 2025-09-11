pruning_unbal_v3 <- function(y, firmid, id, id_old, firmid_old, controls, prov_indicator = rep(1, length(y))) {
  n_of_bad_workers <- 1

  while (n_of_bad_workers >= 1) {
    # find movers
J <- max(firmid, na.rm = TRUE)

# Directly calculate move without multiple assignments
gcs <- c(NA, head(id, -1))
lagfirmid <- c(NA, head(firmid, -1))
move <- ifelse(is.na(firmid) | is.na(lagfirmid) | id != gcs, 0, firmid != lagfirmid)
move <- as.integer(ave(move, id, FUN = sum, na.rm = TRUE)) > 0
move <- move[!is.na(id)]

# Timer start
tic()

# Keep original id_movers before normalization
id_movers_orig <- id[move]

# Normalizing id_movers and firmid_movers
id_movers <- as.integer(factor(id_movers_orig))
firmid_movers <- firmid[move]
n = max(id_movers)
m = max(firmid_movers, na.rm = TRUE)
# Direct sparse matrix creation for bipartite graph
# Create the zero blocks without specifying 'i', 'j', or 'x'
      zero_block_n <- sparseMatrix(
        i = n,      # Row indices for non-zero entries
        j = n,  # Column indices for non-zero entries
        x = 0,  # Values at the non-zero entries (all 1s)
        dims = c(n, n)  # Dimensions of the matrix
      )
      zero_block_J <- sparseMatrix(
        i = m,      # Row indices for non-zero entries
        j = m,  # Column indices for non-zero entries
        x = 0,  # Values at the non-zero entries (all 1s)
        dims = c(m, m)  # Dimensions of the matrix
      )



B <- sparseMatrix(i = id_movers, j = firmid_movers, x = 1, dims = c(n, m))

# Construct the block matrix G using rbind and cbind to bind the matrices by rows and columns
      G <- rbind2(cbind2(zero_block_n, B), cbind2(t(B), zero_block_J))


# Find articulation points
g <- graph_from_adjacency_matrix(G, mode = 'undirected', diag = FALSE)
bad_workers <- articulation_points(g)[articulation_points(g) <= n]

   # Translate back to original worker IDs
sel <- id_movers %in% bad_workers
bad_workers <- unique(id_movers_orig[sel])
n_of_bad_workers <- length(bad_workers)
cat("Number of bad workers", n_of_bad_workers, "\n")


# Remove bad workers and select largest connected set
sel <- !id %in% bad_workers
flush.console()
# Define a function to apply the selection based on data type
apply_selection <- function(data, sel) {
  if(is.matrix(data) || is.array(data)) {
    return(data[sel, , drop = FALSE])
  } else {
    return(data[sel])
  }
}

# Update all data vectors/matrices in one go
data_list <- list(y, firmid, id, firmid_old, id_old, controls, prov_indicator)
data_list <- lapply(data_list, apply_selection, sel = sel)

# Reset IDs
data_list[[2]] <- match(data_list[[2]], unique(data_list[[2]]))  # firmid
data_list[[3]] <- match(data_list[[3]], unique(data_list[[3]]))  # id

flush.console()
# Largest connected set once removed bad workers
A <- build_adj(data_list[[3]], data_list[[2]])
g <- graph_from_adjacency_matrix(A, mode = 'undirected', diag=FALSE)
firm_list <- which(components(g)$membership == which.max(components(g)$csize))
flush.console()
# Update all data vectors/matrices in one go
sel <- data_list[[2]] %in% firm_list
data_list <- lapply(data_list, apply_selection, sel = sel)

# Reset IDs one last time
data_list[[2]] <- match(data_list[[2]], unique(data_list[[2]]))  # firmid
data_list[[3]] <- match(data_list[[3]], unique(data_list[[3]]))  # id
flush.console()
# Extract back the data from data_list
y <- data_list[[1]]
firmid <- data_list[[2]]
id <- data_list[[3]]
firmid_old <- data_list[[4]]
id_old <- data_list[[5]]
controls <- data_list[[6]]
prov_indicator <- data_list[[7]]
gc()
cat("End of one loop to find largest leave-one-out connected set. Will continue if #bad workers = 0\n")
flush.console()
toc()
  }

  return(list(y = y, firmid = firmid, id = id, id_old = id_old, firmid_old = firmid_old, controls = controls, prov_indicator = prov_indicator))
}