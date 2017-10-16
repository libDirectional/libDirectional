#ifndef BINGHAM_HYPERSPHERE_H
#define BINGHAM_HYPERSPHERE_H

#include "tetramesh.h"
#include "octetramesh.h"


/** Tools for creating a finite element representation of a hypersphere **/


/*
 * Fast NN-radius and K-NN searches:
 *  - Add a cell-based data structure based on the tetramesh tessellation,
 *    where for each cell we store an ordered list of the min distances to all the other cells.
 *  - When we create a new point set on the hypersphere, add the points to a data structure
 *    containing an array of point lists per cell.
 *  - On lookup (NN-radius query), lookup the query point's cell with a kdtree NN lookup, then iterate through the
 *    points in cells with min distance < radius.
 */


typedef struct {
  int n;                   // number of cells
  int d;                   // hypersphere dimension + 1
  union {                  // cell mesh
    tetramesh_t *tetramesh;
    // trimesh_t *trimesh;
  };
  double **centroids;              // tetrahedron centroids
  double *volumes;                 // tetrahedron volumes
  double *radii;                   // radii of tetrahedron circumspheres
  int **proximity_queues;          // min-dist-ordered cell lists
  double **proximity_queue_dists;  // distances of cells in proximity queues
  int proximity_queue_size;        // length of proximity queues (i.e. # cells)
} hypersphere_tessellation_t;


typedef struct {
  hypersphere_tessellation_t *tessellation;
  double **points;                           // list of points to store
  int num_points;
  int **cell_members;                        // list of point indices in each cell
  int *cell_counts;                          // number of points in each cell
} hypersphere_pointset_t;


void hypersphere_init();
hypersphere_tessellation_t *tessellate_S3(int n);
hypersphere_pointset_t *new_hypersphere_pointset(hypersphere_tessellation_t *tessellation, double **points, int n);

#endif
