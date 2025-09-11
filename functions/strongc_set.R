strongc_set <- function(y, id, firmid, controls, min_degree=1 ) {
  
  # Assume that y, id, firmid, and controls are available and controls is a matrix
  DT <- data.table(y = y, id = id, firmid)
  DT_controls <- as.data.table(controls)
  
  # Lagfirmid
  DT[, lagfirmid := shift(firmid), by = id]
  DT[is.na(id) | shift(is.na(id), type = "lead", fill = TRUE), lagfirmid := NA]
  
  
  # Find connected set
  y <- DT$y
  id <- DT$id
  firmid <- DT$firmid
  lagfirmid <- DT$lagfirmid
  controls <- as.matrix(DT_controls)
  
  
  # Save old IDs
  firmid_old <- firmid
  id_old <- id
  
  # Create a selection vector for non-NA lagfirmid
  sel <- !is.na(lagfirmid)
  
  # Relabel firms
  firms <- unique(c(firmid, lagfirmid[sel]))
  ids <- unique(id)
  
  
  firmid <- match(firmid, firms)
  lagfirmid[sel] <- match(lagfirmid[sel], firms)
  id <- match(id, ids)
  
  
  # Create adjacency matrix: sparseMatrix is already efficient
  A <- Matrix::sparseMatrix(i = lagfirmid[sel], j = firmid[sel], x = 1)
  
  
  # Make the matrix square by extending the smaller dimension
  max_dim <- max(dim(A))
  A <- Matrix::sparseMatrix(i = c(lagfirmid[sel], max_dim), j = c(firmid[sel], max_dim), x = c(rep(1, length(lagfirmid[sel])), 0))
  
  cat("Dimensions of A after adjustment:", dim(A), "\n")
  
  
  # Remove disconnected vertices: using sparse matrix operation
  degrees <- Matrix::rowSums(A) + Matrix::colSums(A)
  connected_firms <- which(degrees > 1)
  A <- A[connected_firms, connected_firms]  # subsetting sparse matrix is efficient
  
  # Create igraph object and find connected components
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag=FALSE, add.colnames=NA, add.rownames=NA)
  d_g <- degree(g)
  cat("Mean Degree of g is:", mean(d_g), "\n")
  # Filtering based on largest connected component's membership
  firmlst <- connected_firms[d_g >= min_degree]
  
  
  # Select only observations in the largest connected set
  sel <- firmid %in% firmlst
  # Directly create data.tables within the function
  DT_result <- data.table(
    y = y[sel],
    id = id[sel],
    firmid = firmid[sel],
    id_old = id_old[sel],
    firmid_old = firmid_old[sel]
  )
  
  DT_controls_result <- as.data.table(controls[sel, , drop = FALSE])
  
  # Return a list containing both data.tables
  return(list(DT = DT_result, DT_controls = DT_controls_result))
  
}
