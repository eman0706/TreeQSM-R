# This file is part of the R port of TREEQSM.
# Original TREEQSM: Raumonen et al. (2013), distributed under the GNU General Public License, version 3 or (at your option) any later version.
# R port: Barca, Greco, Tartarino (2026); developed with an AI assistant (see README).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Create input parameters for TreeQSM
#'
#' @param name Name of the tree/model
#' @param point_cloud Point cloud (optional, used for auto-parameter determination)
#' @param PatchDiam1 Patch diameter for first cover (default 0.05)
#' @param PatchDiam2Min Minimum patch diameter for second cover (default 0.02)
#' @param PatchDiam2Max Maximum patch diameter for second cover (default 0.05)
#' @param BallRad1 Ball radius for first cover (default 0.05)
#' @param BallRad2 Ball radius for second cover (default 0.05)
#' @param nmin1 Minimum points in BallRad1 balls (default 3)
#' @param nmin2 Minimum points in BallRad2 balls (default 1)
#' @param OnlyTree Whether point cloud contains only tree (default TRUE)
#' @param Tria Whether to create triangulation (default FALSE)
#' @param Dist Whether to compute distances (default TRUE)
#' @param MinCylRad Minimum cylinder radius (default 0.005)
#' @param ParentCor Whether to apply parent radius correction (default TRUE)
#' @param TaperCor Whether to apply taper correction (default TRUE)
#' @param GrowthVolCor Whether to apply growth volume correction (default FALSE)
#' @param GrowthVolFac Growth volume factor (default 0.3)
#' @param tree Tree index (default 1)
#' @param model Model index (default 1)
#' @param savemat Whether to save as .mat file (default FALSE)
#' @param savetxt Whether to save as .txt files (default FALSE)
#' @param plot Plot level: 0=none, 1=models, 2=distributions (default 0)
#' @param disp Display level: 0=name, 1=parameters, 2=times (default 1)
#' @return List with input parameters
#' @export
create_input <- function(name = "Tree", 
                        point_cloud = NULL,
                        PatchDiam1 = 0.05,
                        PatchDiam2Min = 0.02,
                        PatchDiam2Max = 0.05,
                        BallRad1 = 0.05,
                        BallRad2 = 0.05,
                        nmin1 = 3,
                        nmin2 = 1,
                        OnlyTree = TRUE,
                        Tria = FALSE,
                        Dist = TRUE,
                        MinCylRad = 0.005,
                        ParentCor = TRUE,
                        TaperCor = TRUE,
                        GrowthVolCor = FALSE,
                        GrowthVolFac = 0.3,
                        tree = 1,
                        model = 1,
                        savemat = FALSE,
                        savetxt = FALSE,
                        plot = 0,
                        disp = 1) {
  
  # Auto-determine parameters from point cloud if provided
  if (!is.null(point_cloud) && is.matrix(point_cloud)) {
    P <- point_cloud[, 1:3]
    
    # Estimate appropriate patch diameters based on point cloud spread
    heights <- range(P[, 3])
    height_range <- heights[2] - heights[1]
    
    # Determine patch diameters based on tree size
    PatchDiam1 <- define_input(P, "PatchDiam1", height_range)
    PatchDiam2Min <- define_input(P, "PatchDiam2Min", height_range)
    PatchDiam2Max <- define_input(P, "PatchDiam2Max", height_range)
    BallRad1 <- define_input(P, "BallRad1", height_range)
    BallRad2 <- define_input(P, "BallRad2", height_range)
  }
  
  inputs <- list(
    name = name,
    PatchDiam1 = PatchDiam1,
    PatchDiam2Min = PatchDiam2Min,
    PatchDiam2Max = PatchDiam2Max,
    BallRad1 = BallRad1,
    BallRad2 = BallRad2,
    nmin1 = nmin1,
    nmin2 = nmin2,
    OnlyTree = OnlyTree,
    Tria = Tria,
    Dist = Dist,
    MinCylRad = MinCylRad,
    ParentCor = ParentCor,
    TaperCor = TaperCor,
    GrowthVolCor = GrowthVolCor,
    GrowthVolFac = GrowthVolFac,
    tree = tree,
    model = model,
    savemat = savemat,
    savetxt = savetxt,
    plot = plot,
    disp = disp
  )
  
  return(inputs)
}

#' Define input parameter based on point cloud
define_input <- function(P, param_name, height_range) {
  # Simple heuristic based on tree height
  if (height_range < 2) {
    # Small tree
    switch(param_name,
           "PatchDiam1" = 0.03,
           "PatchDiam2Min" = 0.015,
           "PatchDiam2Max" = 0.03,
           "BallRad1" = 0.03,
           "BallRad2" = 0.03,
           0.05)
  } else if (height_range < 10) {
    # Medium tree
    switch(param_name,
           "PatchDiam1" = 0.05,
           "PatchDiam2Min" = 0.02,
           "PatchDiam2Max" = 0.05,
           "BallRad1" = 0.05,
           "BallRad2" = 0.05,
           0.05)
  } else {
    # Large tree
    switch(param_name,
           "PatchDiam1" = 0.1,
           "PatchDiam2Min" = 0.03,
           "PatchDiam2Max" = 0.1,
           "BallRad1" = 0.1,
           "BallRad2" = 0.1,
           0.1)
  }
}

#' Load point cloud from file
#'
#' @param filename Path to file (.txt, .csv, .ply, .las)
#' @param format Format: "auto", "txt", "csv", "ply", "las"
#' @return Point cloud matrix (n x 3)
#' @export
load_point_cloud <- function(filename, format = "auto") {
  ext <- tools::file_ext(filename)
  
  if (format == "auto") {
    format <- tolower(ext)
  }
  
  switch(format,
         "txt" = {
           P <- read.table(filename, header = FALSE)
           as.matrix(P[, 1:3])
         },
         "csv" = {
           P <- read.csv(filename, header = FALSE)
           as.matrix(P[, 1:3])
         },
         "ply" = {
           # Basic PLY reader (vertices only)
           read_ply_vertices(filename)
         },
         "las" = {
           # Requires lidR package
           if (requireNamespace("lidR", quietly = TRUE)) {
             las <- lidR::readLAS(filename)
             return(as.matrix(las@data[, c("X", "Y", "Z")]))
           } else {
             stop("lidR package required for .las files")
           }
         },
         stop("Unsupported format: ", format))
}

#' Basic PLY vertex reader
read_ply_vertices <- function(filename) {
  lines <- readLines(filename)
  
  # Find vertex count
  vertex_line <- grep("^element vertex", lines)
  if (length(vertex_line) == 0) {
    stop("No vertex count found in PLY file")
  }
  n_verts <- as.integer(gsub("element vertex ", "", lines[vertex_line]))
  
  # Find header end
  header_end <- grep("^end_header", lines)
  
  # Parse vertices
  vertices <- matrix(0, n_verts, 3)
  for (i in 1:n_verts) {
    idx <- header_end + i
    parts <- strsplit(lines[idx], "\\s+")[[1]]
    vertices[i, ] <- as.numeric(parts[1:3])
  }
  
  return(vertices)
}
