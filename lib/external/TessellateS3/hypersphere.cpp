#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "util.h"
#include "hypersphere.h"

// cached tessellations of S^3
#define MAX_LEVELS 20
hypersphere_tessellation_t tessellations[MAX_LEVELS];

#define PROXIMITY_QUEUE_SIZE 100
#define PROXIMITY_QUEUE_MEMORY_LIMIT 1e6

/*
 * Reproject mesh vertices onto the unit hypersphere.
 */
static void reproject_vertices(double **vertices, int nv, int dim)
{
  int i;
  for (i = 0; i < nv; i++) {
    double d = norm(vertices[i], dim);
    mult(vertices[i], vertices[i], 1/d, dim);
  }
}

/*
 * Create the initial (low-res) mesh of S3 (in R4).
 */
static octetramesh_t *init_mesh_S3_octetra()
{
  octetramesh_t *mesh;
  safe_calloc(mesh, 1, octetramesh_t);
  octetramesh_new(mesh, 8, 16, 0, 4);

  double vertices[8][4] = {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1},
			   {-1,0,0,0}, {0,-1,0,0}, {0,0,-1,0}, {0,0,0,-1}};

  int tetrahedra[16][4] = {{0,1,2,3}, {0,1,2,7}, {0,1,6,3}, {0,5,2,3},
			   {0,1,6,7}, {0,5,2,7}, {0,5,6,3}, {0,5,6,7},
			   {4,1,2,3}, {4,1,2,7}, {4,1,6,3}, {4,5,2,3},
			   {4,1,6,7}, {4,5,2,7}, {4,5,6,3}, {4,5,6,7}};

  memcpy(mesh->vertices[0], vertices, 8*4*sizeof(double));
  memcpy(mesh->tetrahedra[0], tetrahedra, 16*4*sizeof(int));

  return mesh;
}

/*
 * Create a tesselation of the 3-sphere (in R4) at a given level (# of subdivisions)
 */
static octetramesh_t *build_octetra(int level)
{
  int i;
  octetramesh_t *mesh = init_mesh_S3_octetra();
  octetramesh_t tmp;

  for (i = 0; i < level; i++) {
    octetramesh_subdivide(&tmp, mesh);
    octetramesh_free(mesh);
    free(mesh);
    mesh = octetramesh_clone(&tmp);
    octetramesh_free(&tmp);
  }

  reproject_vertices(mesh->vertices, mesh->nv, mesh->d);

  return mesh;
}


/*
 * Fill in the fields of a hypersphere_tessellation from an octetramesh.
 */
static void octetramesh_to_tessellation(hypersphere_tessellation_t *T, octetramesh_t *mesh)
{
  // get tetramesh
  T->tetramesh = octetramesh_to_tetramesh(mesh);

  // dimensions
  T->d = 4;
  T->n = T->tetramesh->nt;
  int n = T->n, d = T->d;

  // compute cell centroids and volumes
  T->centroids = new_matrix2(n, d);
  safe_calloc(T->volumes, n, double);
  tetramesh_centroids(T->centroids, T->volumes, T->tetramesh);
  reproject_vertices(T->centroids, n, d);
}

/*
 * Pre-cache some tessellations of hyperspheres.
 */
void hypersphere_init()
{
  int i;
  const int levels = 5; //7;

  memset(tessellations, 0, MAX_LEVELS*sizeof(hypersphere_tessellation_t));

  for (i = 0; i < levels; i++) {
    octetramesh_t *mesh = build_octetra(i);
    octetramesh_to_tessellation(&tessellations[i], mesh);
  }
}

/*
 * Returns a tesselation of the 3-sphere (in R4) with at least n cells.
 */
hypersphere_tessellation_t *tessellate_S3(int n)
{
  int i;
  for (i = 0; i < MAX_LEVELS; i++) {
    if (tessellations[i].tetramesh == NULL) {
      octetramesh_t *mesh = build_octetra(i);
      octetramesh_to_tessellation(&tessellations[i], mesh);
    }
    if (tessellations[i].tetramesh->nt >= n)
      return &tessellations[i];
  }

  return &tessellations[MAX_LEVELS-1];
}

