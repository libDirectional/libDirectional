#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "util.h"
#include "tetramesh.h"


/** Tools for subdividing and smoothing a tetrahedral mesh **/

meshgraph_t *tetramesh_meshgraph(tetramesh_t *mesh)
{
  int i, i0, i1, i2, i3, nv = mesh->nv, nt = mesh->nt;

  meshgraph_t *g = meshgraph_new(nv, 10);

  // add the vertices
  for (i = 0; i < nv; i++)
    g->vertices[i] = i;
  g->nv = nv;

  // add faces (and edges)
  for (i = 0; i < nt; i++) {
    i0 = mesh->tetrahedra[i][0];
    i1 = mesh->tetrahedra[i][1];
    i2 = mesh->tetrahedra[i][2];
    i3 = mesh->tetrahedra[i][3];
    meshgraph_add_face(g, i0, i1, i2);
    meshgraph_add_face(g, i0, i1, i3);
    meshgraph_add_face(g, i0, i2, i3);
    meshgraph_add_face(g, i1, i2, i3);
  }

  return g;
}


/*
 * Convert a tetramesh to a set of centroids and volumes.
 */
void tetramesh_centroids(double **centroids, double *volumes, tetramesh_t *mesh)
{
  int i, d = mesh->d;
  //double x[d], *v0, *v1, *v2, *v3;
  double x[4], *v0, *v1, *v2, *v3;

  for (i = 0; i < mesh->nt; i++) {   // for each tetrahedron, compute the center and volume
    memset(x, 0, d*sizeof(double));
    v0 = mesh->vertices[ mesh->tetrahedra[i][0] ];
    v1 = mesh->vertices[ mesh->tetrahedra[i][1] ];
    v2 = mesh->vertices[ mesh->tetrahedra[i][2] ];
    v3 = mesh->vertices[ mesh->tetrahedra[i][3] ];

    add(x, x, v0, d);
    add(x, x, v1, d);
    add(x, x, v2, d);
    add(x, x, v3, d);
    mult(centroids[i], x, 1/4.0, d);

    volumes[i] = tetrahedron_volume(v0, v1, v2, v3, d);
  }
}


/*
 * Build the graph of a mesh.
 */
graph_t *tetramesh_graph(tetramesh_t *mesh)
{
  int i, i0, i1, i2, i3;
  graph_t *g;
  safe_malloc(g, 1, graph_t);

  g->nv = mesh->nv;
  safe_calloc(g->vertices, g->nv, vertex_t);

  for (i = 0; i < g->nv; i++)
    g->vertices[i].index = i;

  // build neighbor lists
  for (i = 0; i < mesh->nt; i++) {
    i0 = mesh->tetrahedra[i][0];
    i1 = mesh->tetrahedra[i][1];
    i2 = mesh->tetrahedra[i][2];
    i3 = mesh->tetrahedra[i][3];

    vertex_t *v0 = &g->vertices[i0];
    vertex_t *v1 = &g->vertices[i1];
    vertex_t *v2 = &g->vertices[i2];
    vertex_t *v3 = &g->vertices[i3];

    if (!ilist_contains(v0->neighbors, i1))  v0->neighbors = ilist_add(v0->neighbors, i1);
    if (!ilist_contains(v0->neighbors, i2))  v0->neighbors = ilist_add(v0->neighbors, i2);
    if (!ilist_contains(v0->neighbors, i3))  v0->neighbors = ilist_add(v0->neighbors, i3);

    if (!ilist_contains(v1->neighbors, i0))  v1->neighbors = ilist_add(v1->neighbors, i0);
    if (!ilist_contains(v1->neighbors, i2))  v1->neighbors = ilist_add(v1->neighbors, i2);
    if (!ilist_contains(v1->neighbors, i3))  v1->neighbors = ilist_add(v1->neighbors, i3);

    if (!ilist_contains(v2->neighbors, i0))  v2->neighbors = ilist_add(v2->neighbors, i0);
    if (!ilist_contains(v2->neighbors, i1))  v2->neighbors = ilist_add(v2->neighbors, i1);
    if (!ilist_contains(v2->neighbors, i3))  v2->neighbors = ilist_add(v2->neighbors, i3);

    if (!ilist_contains(v3->neighbors, i0))  v3->neighbors = ilist_add(v3->neighbors, i0);
    if (!ilist_contains(v3->neighbors, i1))  v3->neighbors = ilist_add(v3->neighbors, i1);
    if (!ilist_contains(v3->neighbors, i2))  v3->neighbors = ilist_add(v3->neighbors, i2);
  }

  // get the total number of edges:  ne = sum(degree_i) / 2
  g->ne = 0;
  for (i = 0; i < g->nv; i++)
    g->ne += g->vertices[i].neighbors->len;
  g->ne /= 2;

  // allocate space for the vertex edge lists
  for (i = 0; i < g->nv; i++)
    safe_calloc(g->vertices[i].edges, g->vertices[i].neighbors->len, int);

  // get the edges, and build the vertex edge lists
  safe_calloc(g->edges, g->ne, edge_t);
  int cnt = 0;
  for (i = 0; i < g->nv; i++) {

    ilist_t *tmp;
    int ti = 0;  // index of neighbor list
    for (tmp = g->vertices[i].neighbors; tmp; tmp = tmp->next, ti++) {

      if (tmp->x >= i) {
	int j = tmp->x;

	// add the edge to g->edges
	g->edges[cnt].i = i;
	g->edges[cnt].j = j;
	//g->edges[cnt].len = dist(mesh->vertices[i], mesh->vertices[j], mesh->d);

	// add the edge index to the vertex edge lists of vertices i and j
	g->vertices[i].edges[ti] = cnt;
	int tj = ilist_find(g->vertices[j].neighbors, i);
	g->vertices[j].edges[tj] = cnt;

	cnt++;
      }
    }
  }

  return g;
}

/*
 * Compute the midpoints of all the edges in a mesh.
 */
void tetramesh_midpoints(double **dst, tetramesh_t *mesh, graph_t *graph)
{
  int i, j, e;
  for (e = 0; e < graph->ne; e++) {   // for each edge, compute the midpoint
    i = graph->edges[e].i;
    j = graph->edges[e].j;
    avg(dst[e], mesh->vertices[i], mesh->vertices[j], mesh->d);
  }
}

/*
 * Subdivide each tetrahedron in a mesh into 8 smaller tetrahedra.
 */
void tetramesh_subdivide(tetramesh_t *dst, tetramesh_t *src)
{
  int i;
  int p1, p2, p3, p4;
  int q12, q13, q14, q23, q24, q34;

  int nv = src->nv;
  int nt = src->nt;
  int d = src->d;

  // compute the midpoints of each edge in the mesh
  //midpoint_map_t M;
  //tetramesh_midpoints(&M, src);

  graph_t *graph = tetramesh_graph(src);

  int nv2 = nv + graph->ne;  //M.num_edges;  // old vertices plus the midpoints
  int nt2 = 8*nt;

  // allocate space for the new mesh and copy old vertices plus midpoints into dst
  tetramesh_new(dst, nv2, nt2, d);
  memcpy(dst->vertices[0], src->vertices[0], nv*d*sizeof(double));
  //memcpy(dst->vertices[0] + nv*d, M.points[0], (nv2-nv)*d*sizeof(double));
  tetramesh_midpoints(dst->vertices + nv, src, graph);

  for (i = 0; i < nt; i++) {    // for each tetrahedron in the original mesh
    // original point indices
    p1 = src->tetrahedra[i][0];
    p2 = src->tetrahedra[i][1];
    p3 = src->tetrahedra[i][2];
    p4 = src->tetrahedra[i][3];

    // new point indices
    q12 = nv + graph->vertices[p1].edges[ ilist_find(graph->vertices[p1].neighbors, p2) ];  //M.map[p1][p2];
    q13 = nv + graph->vertices[p1].edges[ ilist_find(graph->vertices[p1].neighbors, p3) ];  //M.map[p1][p3];
    q14 = nv + graph->vertices[p1].edges[ ilist_find(graph->vertices[p1].neighbors, p4) ];  //M.map[p1][p4];
    q23 = nv + graph->vertices[p2].edges[ ilist_find(graph->vertices[p2].neighbors, p3) ];  //M.map[p2][p3];
    q24 = nv + graph->vertices[p2].edges[ ilist_find(graph->vertices[p2].neighbors, p4) ];  //M.map[p2][p4];
    q34 = nv + graph->vertices[p3].edges[ ilist_find(graph->vertices[p3].neighbors, p4) ];  //M.map[p3][p4];

    int *t0 = dst->tetrahedra[8*i];
    int *t1 = dst->tetrahedra[8*i+1];
    int *t2 = dst->tetrahedra[8*i+2];
    int *t3 = dst->tetrahedra[8*i+3];
    int *t4 = dst->tetrahedra[8*i+4];
    int *t5 = dst->tetrahedra[8*i+5];
    int *t6 = dst->tetrahedra[8*i+6];
    int *t7 = dst->tetrahedra[8*i+7];

    t0[0] = p1;  t0[1] = q12;  t0[2] = q13;  t0[3] = q14;
    t1[0] = p2;  t1[1] = q12;  t1[2] = q23;  t1[3] = q24;
    t2[0] = p3;  t2[1] = q13;  t2[2] = q23;  t2[3] = q34;
    t3[0] = p4;  t3[1] = q14;  t3[2] = q24;  t3[3] = q34;
    t4[0] = q13;  t4[1] = q12;  t4[2] = q23;  t4[3] = q34;
    t5[0] = q13;  t5[1] = q12;  t5[2] = q14;  t5[3] = q34;
    t6[0] = q24;  t6[1] = q12;  t6[2] = q23;  t6[3] = q34;
    t7[0] = q24;  t7[1] = q12;  t7[2] = q14;  t7[3] = q34;
  }

  graph_free(graph);

  //free_midpoint_map(&M);
}


/*
 * Copy the contents of one mesh into another mesh.
 */
void tetramesh_copy(tetramesh_t *dst, tetramesh_t *src)
{
  memcpy(dst->vertices[0], src->vertices[0], src->nv*src->d*sizeof(double));  // copy raw vertices
  memcpy(dst->tetrahedra[0], src->tetrahedra[0], 4*src->nt*sizeof(int));      // copy raw tetrahedra
}


/*
 * Copy the vertices of one mesh into another mesh.
 */
void tetramesh_copy_vertices(tetramesh_t *dst, tetramesh_t *src)
{
  memcpy(dst->vertices[0], src->vertices[0], src->nv*src->d*sizeof(double));  // copy raw vertices
}


/*
 * Copy the tetrahedra of one mesh into another mesh.
 */
void tetramesh_copy_tetrahedra(tetramesh_t *dst, tetramesh_t *src)
{
  memcpy(dst->tetrahedra[0], src->tetrahedra[0], 4*src->nt*sizeof(int));      // copy raw tetrahedra
}


/*
 * Clone a mesh.
 */
tetramesh_t *tetramesh_clone(tetramesh_t *src)
{
  tetramesh_t *dst;
  safe_calloc(dst, 1, tetramesh_t);
  tetramesh_new(dst, src->nv, src->nt, src->d);
  tetramesh_copy(dst, src);

  return dst;
}


/*
 * Create (allocate) the contents of a mesh.
 */
void tetramesh_new(tetramesh_t *mesh, int nv, int nt, int d)
{
  mesh->nv = nv;
  mesh->nt = nt;
  mesh->d = d;

  mesh->vertices = new_matrix2(nv, d);
  mesh->tetrahedra = new_matrix2i(nt, 4);
}



/*
 * Free the contents of a tetrahedral mesh.
 */
void tetramesh_free(tetramesh_t *mesh)
{
  if (mesh->nv != 0) {
    //free(mesh->vertices[0]);
    //free(mesh->vertices);
    //free(mesh->tetrahedra[0]);  // tetra_raw
    //free(mesh->tetrahedra);

    free_matrix2(mesh->vertices);
    free_matrix2i(mesh->tetrahedra);
  }
}


/*
 *Save a colored tetrahedral mesh to PLY file.
 */
void tetramesh_save_PLY_colors(tetramesh_t *mesh, meshgraph_t *graph, char *filename, int *colors)
{
  FILE *f = fopen(filename, "w");

  int i, j, i0, i1, i2, i3;

  int num_faces = graph->nf;

  fprintf(f, "ply\n");
  fprintf(f, "format ascii 1.0\n");
  fprintf(f, "comment tetramesh model\n");  
  fprintf(f, "element vertex %d\n", mesh->nv);
  fprintf(f, "property float x\n");
  fprintf(f, "property float y\n");
  fprintf(f, "property float z\n");
  fprintf(f, "element face %d\n", num_faces);
  fprintf(f, "property list uchar int vertex_indices\n");
  if (colors) {
    fprintf(f, "property uchar red\n");
    fprintf(f, "property uchar green\n");
    fprintf(f, "property uchar blue\n");
  }
  fprintf(f, "end_header\n");

  for (i = 0; i < mesh->nv; i++) {
    for (j = 0; j < 3 /*mesh->d*/; j++)
      fprintf(f, "%f ", mesh->vertices[i][j]);
    fprintf(f, "\n");
  }

  int *face_colors = NULL;
  if (colors) {
    // determine what color each face should be
    safe_calloc(face_colors, num_faces, int);
    for (i = 0; i < mesh->nt; i++) {
      i0 = mesh->tetrahedra[i][0];
      i1 = mesh->tetrahedra[i][1];
      i2 = mesh->tetrahedra[i][2];
      i3 = mesh->tetrahedra[i][3];

      int face = meshgraph_find_face(graph, i0, i1, i2);
      if (colors[i] > face_colors[face])
	face_colors[face] = colors[i];

      face = meshgraph_find_face(graph, i0, i1, i3);
      if (colors[i] > face_colors[face])
	face_colors[face] = colors[i];

      face = meshgraph_find_face(graph, i0, i2, i3);
      if (colors[i] > face_colors[face])
	face_colors[face] = colors[i];

      face = meshgraph_find_face(graph, i1, i2, i3);
      if (colors[i] > face_colors[face])
	face_colors[face] = colors[i];
    }
  }

  // ~~~ TODO: Use a graph traversal algorithm to avoid double-counting! ~~~
  for (i = 0; i < num_faces; i++) {
    face_t face = graph->faces[i];
    if (colors) {
      color_t color = colormap[ face_colors[i] ];
      fprintf(f, "3 %d %d %d %d %d %d\n", face.i, face.j, face.k, color.r, color.g, color.b);
    }
    else
      fprintf(f, "3 %d %d %d\n", face.i, face.j, face.k);
  }

  /*
  for (i = 0; i < mesh->nt; i++) {
    i0 = mesh->tetrahedra[i][0];
    i1 = mesh->tetrahedra[i][1];
    i2 = mesh->tetrahedra[i][2];
    i3 = mesh->tetrahedra[i][3];
    if (colors) {
      fprintf(f, "3 %d %d %d %d %d %d\n", i0, i1, i2, colors[i].r, colors[i].g, colors[i].b);
      fprintf(f, "3 %d %d %d %d %d %d\n", i0, i1, i3, colors[i].r, colors[i].g, colors[i].b);
      fprintf(f, "3 %d %d %d %d %d %d\n", i0, i2, i3, colors[i].r, colors[i].g, colors[i].b);
      fprintf(f, "3 %d %d %d %d %d %d\n", i1, i2, i3, colors[i].r, colors[i].g, colors[i].b);
    }
    else {
      fprintf(f, "3 %d %d %d\n", i0, i1, i2);
      fprintf(f, "3 %d %d %d\n", i0, i1, i3);
      fprintf(f, "3 %d %d %d\n", i0, i2, i3);
      fprintf(f, "3 %d %d %d\n", i1, i2, i3);
    }
  }
  */

  fclose(f);
}

/*
 *Save a tetrahedral mesh to PLY file.
 */
void tetramesh_save_PLY(tetramesh_t *mesh, meshgraph_t *graph, char *filename)
{
  tetramesh_save_PLY_colors(mesh, graph, filename, 0);
}


/*
 * Compute statistics about tetrahedral areas, edge lengths, etc.
 */
tetramesh_stats_t tetramesh_stats(tetramesh_t *T)
{
  int i, j, e, nt = T->nt, d = T->d;
  tetramesh_stats_t stats;

  graph_t *graph = tetramesh_graph(T);
  stats.num_edges = graph->ne;

  stats.num_vertices = T->nv;
  stats.num_tetrahedra = T->nt;

  stats.min_edge_len = DBL_MAX;
  stats.max_edge_len = 0;
  stats.avg_edge_len = 0;
  stats.std_edge_len = 0;
  stats.min_skewness = DBL_MAX;
  stats.max_skewness = 0;
  stats.avg_skewness = 0;
  stats.std_skewness = 0;
  stats.min_volume = DBL_MAX;
  stats.max_volume = 0;
  stats.avg_volume = 0;
  stats.std_volume = 0;

  for (i = 0; i < nt; i++) {      // for each tetrahedron

    int i0 = T->tetrahedra[i][0];
    int i1 = T->tetrahedra[i][1];
    int i2 = T->tetrahedra[i][2];
    int i3 = T->tetrahedra[i][3];

    double *p0 = T->vertices[i0];
    double *p1 = T->vertices[i1];
    double *p2 = T->vertices[i2];
    double *p3 = T->vertices[i3];

    double d01 = dist(p0, p1, d);
    double d02 = dist(p0, p2, d);
    double d03 = dist(p0, p3, d);
    double d12 = dist(p1, p2, d);
    double d13 = dist(p1, p3, d);
    double d23 = dist(p2, p3, d);

    double d_edge[6] = {d01, d02, d03, d12, d13, d23};
    double dmax = arr_max(d_edge, 6);
    double dmin = arr_min(d_edge, 6);

    double skewness = dmax/dmin;
    if (skewness < stats.min_skewness)
      stats.min_skewness = skewness;
    if (skewness > stats.max_skewness)
      stats.max_skewness = skewness;
    stats.avg_skewness += skewness;

    double volume = tetrahedron_volume(p0, p1, p2, p3, d);
    if (volume < stats.min_volume)
      stats.min_volume = volume;
    if (volume > stats.max_volume)
      stats.max_volume = volume;
    stats.avg_volume += volume;
  }
  stats.avg_skewness /= (double)nt;
  stats.avg_volume /= (double)nt;

  for (i = 0; i < nt; i++) {      // for each tetrahedron

    int i0 = T->tetrahedra[i][0];
    int i1 = T->tetrahedra[i][1];
    int i2 = T->tetrahedra[i][2];
    int i3 = T->tetrahedra[i][3];

    double *p0 = T->vertices[i0];
    double *p1 = T->vertices[i1];
    double *p2 = T->vertices[i2];
    double *p3 = T->vertices[i3];

    double d01 = dist(p0, p1, d);
    double d02 = dist(p0, p2, d);
    double d03 = dist(p0, p3, d);
    double d12 = dist(p1, p2, d);
    double d13 = dist(p1, p3, d);
    double d23 = dist(p2, p3, d);

    double d_edge[6] = {d01, d02, d03, d12, d13, d23};
    double dmax = arr_max(d_edge, 6);
    double dmin = arr_min(d_edge, 6);

    double skewness = dmax/dmin;
    double ds = stats.avg_skewness - skewness;
    stats.std_skewness += ds*ds;

    double volume = tetrahedron_volume(p0, p1, p2, p3, d);
    double dv = stats.avg_volume - volume;
    stats.std_volume += dv*dv;
  }
  stats.std_skewness = sqrt(stats.std_skewness/(double)nt);
  stats.std_volume = sqrt(stats.std_volume/(double)nt);

  for (e = 0; e < graph->ne; e++) {
    i = graph->edges[e].i;
    j = graph->edges[e].j;
    double edge_len = dist(T->vertices[i], T->vertices[j], d);
    if (edge_len < stats.min_edge_len)
      stats.min_edge_len = edge_len;
    if (edge_len > stats.max_edge_len)
      stats.max_edge_len = edge_len;
    stats.avg_edge_len += edge_len;
  }
  stats.avg_edge_len /= (double)stats.num_edges;

  for (e = 0; e < graph->ne; e++) {
    i = graph->edges[e].i;
    j = graph->edges[e].j;
    double edge_len = dist(T->vertices[i], T->vertices[j], d);
    double de = stats.avg_edge_len - edge_len;
    stats.std_edge_len += de*de;
  }
  stats.std_edge_len = sqrt(stats.std_edge_len/(double)stats.num_edges);

  graph_free(graph);

  return stats;
}

/*
 * Print the stats of a mesh.
 */
void tetramesh_print_stats(tetramesh_stats_t stats)
{
  printf("tetramesh stats {\n");
  printf("  nv = %d, ne = %d, nt = %d\n", stats.num_vertices, stats.num_edges, stats.num_tetrahedra);

  printf("  edge_len = [%f, %f], avg: %f, std: %f\n",
	 stats.min_edge_len, stats.max_edge_len, stats.avg_edge_len, stats.std_edge_len);

  printf("  skewness = [%f, %f], avg: %f, std: %f\n",
	 stats.min_skewness, stats.max_skewness, stats.avg_skewness, stats.std_skewness);

  printf("  volume = [%f, %f], avg: %f, std: %f\n",
	 stats.min_volume, stats.max_volume, stats.avg_volume, stats.std_volume);

  printf("}\n");
}
