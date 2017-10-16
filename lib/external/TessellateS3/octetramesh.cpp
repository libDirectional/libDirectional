#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "util.h"
#include "octetramesh.h"
#include <malloc.h>

#ifdef __GNUC__
    #include <alloca.h>
#endif

/*
 * Convert an octa-tetrahedral mesh to a tetrahedral mesh.
 */
tetramesh_t *octetramesh_to_tetramesh(octetramesh_t *mesh)
{
  int i, i0, i1, i2, i3, i4, i5, nv = mesh->nv, nt = mesh->nt, no = mesh->no, d = mesh->d;
  tetramesh_t *tmesh;
  safe_calloc(tmesh, 1, tetramesh_t);
  tetramesh_new(tmesh, nv, nt + 4*no, d);
  memcpy(tmesh->vertices[0], mesh->vertices[0], d*nv*sizeof(double));
  memcpy(tmesh->tetrahedra[0], mesh->tetrahedra[0], 4*nt*sizeof(int));

  for (i = 0; i < mesh->no; i++) {
    i0 = mesh->octahedra[i][0];
    i1 = mesh->octahedra[i][1];
    i2 = mesh->octahedra[i][2];
    i3 = mesh->octahedra[i][3];
    i4 = mesh->octahedra[i][4];
    i5 = mesh->octahedra[i][5];

    // new tetrahedra
    int *t0 = tmesh->tetrahedra[nt+4*i];
    int *t1 = tmesh->tetrahedra[nt+4*i+1];
    int *t2 = tmesh->tetrahedra[nt+4*i+2];
    int *t3 = tmesh->tetrahedra[nt+4*i+3];

    t0[0] = i0;  t0[1] = i1;  t0[2] = i2;  t0[3] = i3;
    t1[0] = i0;  t1[1] = i1;  t1[2] = i3;  t1[3] = i4;
    t2[0] = i5;  t2[1] = i1;  t2[2] = i2;  t2[3] = i4;
    t3[0] = i5;  t3[1] = i2;  t3[2] = i3;  t3[3] = i4;
  }

  return tmesh;
}

/*
 * Build the graph of a mesh.
 */
graph_t *octetramesh_graph(octetramesh_t *mesh)
{
  //printf("octetramesh_graph()\n");

  int i, i0, i1, i2, i3, i4, i5;
  graph_t *g;
  safe_malloc(g, 1, graph_t);

  g->nv = mesh->nv;
  safe_calloc(g->vertices, g->nv, vertex_t);

  for (i = 0; i < g->nv; i++)
    g->vertices[i].index = i;

  // build neighbor lists 
  for (i = 0; i < mesh->nt; i++) {  // tetrahedra
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
  for (i = 0; i < mesh->no; i++) {  // octahedra
    i0 = mesh->octahedra[i][0];
    i1 = mesh->octahedra[i][1];
    i2 = mesh->octahedra[i][2];
    i3 = mesh->octahedra[i][3];
    i4 = mesh->octahedra[i][4];
    i5 = mesh->octahedra[i][5];

    vertex_t *v0 = &g->vertices[i0];
    vertex_t *v1 = &g->vertices[i1];
    vertex_t *v2 = &g->vertices[i2];
    vertex_t *v3 = &g->vertices[i3];
    vertex_t *v4 = &g->vertices[i4];
    vertex_t *v5 = &g->vertices[i5];

    if (!ilist_contains(v0->neighbors, i1))  v0->neighbors = ilist_add(v0->neighbors, i1);
    if (!ilist_contains(v0->neighbors, i2))  v0->neighbors = ilist_add(v0->neighbors, i2);
    if (!ilist_contains(v0->neighbors, i3))  v0->neighbors = ilist_add(v0->neighbors, i3);
    if (!ilist_contains(v0->neighbors, i4))  v0->neighbors = ilist_add(v0->neighbors, i4);

    if (!ilist_contains(v1->neighbors, i0))  v1->neighbors = ilist_add(v1->neighbors, i0);
    if (!ilist_contains(v1->neighbors, i2))  v1->neighbors = ilist_add(v1->neighbors, i2);
    if (!ilist_contains(v1->neighbors, i4))  v1->neighbors = ilist_add(v1->neighbors, i4);
    if (!ilist_contains(v1->neighbors, i5))  v1->neighbors = ilist_add(v1->neighbors, i5);

    if (!ilist_contains(v2->neighbors, i0))  v2->neighbors = ilist_add(v2->neighbors, i0);
    if (!ilist_contains(v2->neighbors, i1))  v2->neighbors = ilist_add(v2->neighbors, i1);
    if (!ilist_contains(v2->neighbors, i3))  v2->neighbors = ilist_add(v2->neighbors, i3);
    if (!ilist_contains(v2->neighbors, i5))  v2->neighbors = ilist_add(v2->neighbors, i5);

    if (!ilist_contains(v3->neighbors, i0))  v3->neighbors = ilist_add(v3->neighbors, i0);
    if (!ilist_contains(v3->neighbors, i2))  v3->neighbors = ilist_add(v3->neighbors, i2);
    if (!ilist_contains(v3->neighbors, i4))  v3->neighbors = ilist_add(v3->neighbors, i4);
    if (!ilist_contains(v3->neighbors, i5))  v3->neighbors = ilist_add(v3->neighbors, i5);

    if (!ilist_contains(v4->neighbors, i0))  v4->neighbors = ilist_add(v4->neighbors, i0);
    if (!ilist_contains(v4->neighbors, i1))  v4->neighbors = ilist_add(v4->neighbors, i1);
    if (!ilist_contains(v4->neighbors, i3))  v4->neighbors = ilist_add(v4->neighbors, i3);
    if (!ilist_contains(v4->neighbors, i5))  v4->neighbors = ilist_add(v4->neighbors, i5);

    if (!ilist_contains(v5->neighbors, i1))  v5->neighbors = ilist_add(v5->neighbors, i1);
    if (!ilist_contains(v5->neighbors, i2))  v5->neighbors = ilist_add(v5->neighbors, i2);
    if (!ilist_contains(v5->neighbors, i3))  v5->neighbors = ilist_add(v5->neighbors, i3);
    if (!ilist_contains(v5->neighbors, i4))  v5->neighbors = ilist_add(v5->neighbors, i4);
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
static void octetramesh_midpoints(double **dst, octetramesh_t *mesh, graph_t *graph)
{
  int i, j, e;
  for (e = 0; e < graph->ne; e++) {   // for each edge, compute the midpoint
    i = graph->edges[e].i;
    j = graph->edges[e].j;
    avg(dst[e], mesh->vertices[i], mesh->vertices[j], mesh->d);
  }
}


/*
 * Compute the midpoints of the selected edges in a mesh.
 */
static void octetramesh_midpoints_select(double **dst, octetramesh_t *mesh, graph_t *graph, int *edgemask)
{
  int i, j, e, cnt = 0;
  for (e = 0; e < graph->ne; e++) {   // for each selected edge, compute the midpoint
    if (edgemask[e]) {
      i = graph->edges[e].i;
      j = graph->edges[e].j;
      avg(dst[cnt], mesh->vertices[i], mesh->vertices[j], mesh->d);
      cnt++;
    }
  }
}


/*
 * Compute the centers of all the octahedra in a mesh.
 */
static void octetramesh_octa_centers(double **dst, octetramesh_t *mesh)
{
  int i, j, v, d = mesh->d;
  //double x[d];
  double x[4];

  for (i = 0; i < mesh->no; i++) {   // for each octahedra, compute the center
    memset(x, 0, d*sizeof(double));
    for (j = 1; j < 5; j++) {
      v = mesh->octahedra[i][j];
      add(x, x, mesh->vertices[v], d);
    }
    mult(dst[i], x, 1/4.0, d);
  }
}


/*
 * Compute the centers of the selected octahedra in a mesh.
 */
static void octetramesh_octa_centers_select(double **dst, octetramesh_t *mesh, int *octmask)
{
  int i, j, v, cnt, d = mesh->d;
  //double x[d];
  double x[4];

  cnt = 0;
  for (i = 0; i < mesh->no; i++) {   // for each octahedra, compute the center
    if (octmask[i]) {
      memset(x, 0, d*sizeof(double));
      for (j = 1; j < 5; j++) {
	v = mesh->octahedra[i][j];
	add(x, x, mesh->vertices[v], d);
      }
      mult(dst[cnt++], x, 1/4.0, d);
    }
  }
}


/*
 * Mark the edges of a subset of tetrahedra and octahedra.
 */
static void select_edges(int *edgemask, octetramesh_t *src, graph_t *graph, int *tetmask, int *octmask)
{
  int i, p0, p1, p2, p3, p4, p5;
  int ne = graph->ne;
  int nt = src->nt;
  int no = src->no;

  memset(edgemask, 0, ne*sizeof(int));

  for (i = 0; i < nt; i++) {
    if (tetmask[i]) {
      p0 = src->tetrahedra[i][0];
      p1 = src->tetrahedra[i][1];
      p2 = src->tetrahedra[i][2];
      p3 = src->tetrahedra[i][3];
      edgemask[ graph_find_edge(graph, p0, p1) ] = 1;
      edgemask[ graph_find_edge(graph, p0, p2) ] = 1;
      edgemask[ graph_find_edge(graph, p0, p3) ] = 1;
      edgemask[ graph_find_edge(graph, p1, p2) ] = 1;
      edgemask[ graph_find_edge(graph, p1, p3) ] = 1;
      edgemask[ graph_find_edge(graph, p2, p3) ] = 1;
    }
  }

  for (i = 0; i < no; i++) {
    if (octmask[i]) {
      p0 = src->octahedra[i][0];
      p1 = src->octahedra[i][1];
      p2 = src->octahedra[i][2];
      p3 = src->octahedra[i][3];
      p4 = src->octahedra[i][4];
      p5 = src->octahedra[i][5];
      edgemask[ graph_find_edge(graph, p0, p1) ] = 1;
      edgemask[ graph_find_edge(graph, p0, p2) ] = 1;
      edgemask[ graph_find_edge(graph, p0, p3) ] = 1;
      edgemask[ graph_find_edge(graph, p0, p4) ] = 1;
      edgemask[ graph_find_edge(graph, p1, p2) ] = 1;
      edgemask[ graph_find_edge(graph, p2, p3) ] = 1;
      edgemask[ graph_find_edge(graph, p3, p4) ] = 1;
      edgemask[ graph_find_edge(graph, p4, p1) ] = 1;
      edgemask[ graph_find_edge(graph, p5, p1) ] = 1;
      edgemask[ graph_find_edge(graph, p5, p2) ] = 1;
      edgemask[ graph_find_edge(graph, p5, p3) ] = 1;
      edgemask[ graph_find_edge(graph, p5, p4) ] = 1;
    }
  }
}


/**
 * Subdivide each specificied octahedron into 6 octahedra and 8 tetrahedra,
 * and subdivide each specified tetrahedron into 4 tetrahedra and an octahedron.
 */
void octetramesh_subdivide_select(octetramesh_t *dst, octetramesh_t *src, int *tetmask, int *octmask)
{
  int i, tcnt, ocnt;
  int p0, p1, p2, p3, p4, p5;
  int q0, q01, q02, q03, q04, q12, q23, q34, q41, q51, q52, q53, q54, q13;

  graph_t *graph = octetramesh_graph(src);

  int ne = graph->ne;
  int nv = src->nv;
  int nt = src->nt;
  int no = src->no;
  int d = src->d;

  // mark edges for subdivision
  //int edgemask[ne];
  int *edgemask = (int *)alloca(sizeof(int) * ne);

  select_edges(edgemask, src, graph, tetmask, octmask);

  // create an edge map from the full edge list to the selected edge list
  //int edgemap[ne];
  int *edgemap = (int *)alloca(sizeof(int) * ne);
  findinv(edgemap, edgemask, ne);

  // create an oct map from the full octahedron list to the selected octahedron list
  //int octmap[no];
  int *octmap = (int *)alloca(sizeof(int) * no);
  findinv(octmap, octmask, no);

  int tdiv = count(tetmask, nt);   // # of tetrahedra to subdivide
  int odiv = count(octmask, no);   // # of octahedra to subdivide
  int ediv = count(edgemask, ne);  // # of edges to subdivide

  int nv2 = nv + ediv + odiv;  // old vertices + midpoints + octahedral centers
  int nt2 = (nt - tdiv) + 4*tdiv + 8*odiv;
  int no2 = (no - odiv) + tdiv + 6*odiv;

  // allocate space for the new mesh and copy old vertices, midpoints, and octahedral centers into dst
  octetramesh_new(dst, nv2, nt2, no2, d);
  memcpy(dst->vertices[0], src->vertices[0], nv*d*sizeof(double));
  octetramesh_midpoints_select(dst->vertices + nv, src, graph, edgemask);
  octetramesh_octa_centers_select(dst->vertices + nv + ediv, src, octmask);

  // copy unchanged tetrahedra into dst
  tcnt = 0;
  for (i = 0; i < nt; i++)
    if (tetmask[i] == 0)
      memcpy(dst->tetrahedra[tcnt++], src->tetrahedra[i], 4*sizeof(int));

  // copy unchanged octahedra into dst
  ocnt = 0;
  for (i = 0; i < no; i++)
    if (octmask[i] == 0)
      memcpy(dst->octahedra[ocnt++], src->octahedra[i], 6*sizeof(int));

  for (i = 0; i < nt; i++) {    // subdivide selected tetrahedra in the original mesh
    if (tetmask[i]) {
      // original point indices
      p0 = src->tetrahedra[i][0];
      p1 = src->tetrahedra[i][1];
      p2 = src->tetrahedra[i][2];
      p3 = src->tetrahedra[i][3];

      // new point indices
      q01 = nv + edgemap[ graph_find_edge(graph, p0, p1) ];
      q02 = nv + edgemap[ graph_find_edge(graph, p0, p2) ];
      q03 = nv + edgemap[ graph_find_edge(graph, p0, p3) ];
      q12 = nv + edgemap[ graph_find_edge(graph, p1, p2) ];
      q13 = nv + edgemap[ graph_find_edge(graph, p1, p3) ];
      q23 = nv + edgemap[ graph_find_edge(graph, p2, p3) ];

      // new tetrahedra
      int *t0 = dst->tetrahedra[tcnt++];
      int *t1 = dst->tetrahedra[tcnt++];
      int *t2 = dst->tetrahedra[tcnt++];
      int *t3 = dst->tetrahedra[tcnt++];

      t0[0] = p0;  t0[1] = q01;  t0[2] = q02;  t0[3] = q03;
      t1[0] = p1;  t1[1] = q01;  t1[2] = q12;  t1[3] = q13;
      t2[0] = p2;  t2[1] = q02;  t2[2] = q12;  t2[3] = q23;
      t3[0] = p3;  t3[1] = q03;  t3[2] = q13;  t3[3] = q23;

      // new octahedron
      int *oct = dst->octahedra[ocnt++];

      oct[0] = q01;  oct[1] = q02;  oct[2] = q03;  oct[3] = q13;  oct[4] = q12;  oct[5] = q23;
    }
  }

  for (i = 0; i < no; i++) {    // subdivide selected octahedra in the original mesh
    if (octmask[i]) {
      // original point indices
      p0 = src->octahedra[i][0];
      p1 = src->octahedra[i][1];
      p2 = src->octahedra[i][2];
      p3 = src->octahedra[i][3];
      p4 = src->octahedra[i][4];
      p5 = src->octahedra[i][5];

      // new point indices
      q0 = nv + ediv + octmap[i];
      q01 = nv + edgemap[ graph_find_edge(graph, p0, p1) ];
      q02 = nv + edgemap[ graph_find_edge(graph, p0, p2) ];
      q03 = nv + edgemap[ graph_find_edge(graph, p0, p3) ];
      q04 = nv + edgemap[ graph_find_edge(graph, p0, p4) ];
      q12 = nv + edgemap[ graph_find_edge(graph, p1, p2) ];
      q23 = nv + edgemap[ graph_find_edge(graph, p2, p3) ];
      q34 = nv + edgemap[ graph_find_edge(graph, p3, p4) ];
      q41 = nv + edgemap[ graph_find_edge(graph, p4, p1) ];
      q51 = nv + edgemap[ graph_find_edge(graph, p5, p1) ];
      q52 = nv + edgemap[ graph_find_edge(graph, p5, p2) ];
      q53 = nv + edgemap[ graph_find_edge(graph, p5, p3) ];
      q54 = nv + edgemap[ graph_find_edge(graph, p5, p4) ];

      // new tetrahedra
      int *t0 = dst->tetrahedra[tcnt++];
      int *t1 = dst->tetrahedra[tcnt++];
      int *t2 = dst->tetrahedra[tcnt++];
      int *t3 = dst->tetrahedra[tcnt++];
      int *t4 = dst->tetrahedra[tcnt++];
      int *t5 = dst->tetrahedra[tcnt++];
      int *t6 = dst->tetrahedra[tcnt++];
      int *t7 = dst->tetrahedra[tcnt++];

      t0[0] = q0;  t0[1] = q01;  t0[2] = q02;  t0[3] = q12;
      t1[0] = q0;  t1[1] = q02;  t1[2] = q03;  t1[3] = q23;
      t2[0] = q0;  t2[1] = q03;  t2[2] = q04;  t2[3] = q34;
      t3[0] = q0;  t3[1] = q04;  t3[2] = q01;  t3[3] = q41;
      t4[0] = q0;  t4[1] = q51;  t4[2] = q52;  t4[3] = q12;
      t5[0] = q0;  t5[1] = q52;  t5[2] = q53;  t5[3] = q23;
      t6[0] = q0;  t6[1] = q53;  t6[2] = q54;  t6[3] = q34;
      t7[0] = q0;  t7[1] = q54;  t7[2] = q51;  t7[3] = q41;

      // new octahedra
      int *oct0 = dst->octahedra[ocnt++];
      int *oct1 = dst->octahedra[ocnt++];
      int *oct2 = dst->octahedra[ocnt++];
      int *oct3 = dst->octahedra[ocnt++];
      int *oct4 = dst->octahedra[ocnt++];
      int *oct5 = dst->octahedra[ocnt++];

      oct0[0] = p0;  oct0[1] = q01;  oct0[2] = q02;  oct0[3] = q03;  oct0[4] = q04;  oct0[5] = q0;
      oct1[0] = q01;  oct1[1] = p1;  oct1[2] = q12;  oct1[3] = q0;  oct1[4] = q41;  oct1[5] = q51;
      oct2[0] = q02;  oct2[1] = p2;  oct2[2] = q23;  oct2[3] = q0;  oct2[4] = q12;  oct2[5] = q52;
      oct3[0] = q03;  oct3[1] = p3;  oct3[2] = q34;  oct3[3] = q0;  oct3[4] = q23;  oct3[5] = q53;
      oct4[0] = q04;  oct4[1] = p4;  oct4[2] = q41;  oct4[3] = q0;  oct4[4] = q34;  oct4[5] = q54;
      oct5[0] = p5;  oct5[1] = q51;  oct5[2] = q52;  oct5[3] = q53;  oct5[4] = q54;  oct5[5] = q0;
    }
  }

  graph_free(graph);
}


/**
 * Subdivide each octahedron into 6 octahedra and 8 tetrahedra,
 * and subdivide each tetrahedron into 4 tetrahedra and an octahedron.
 */
void octetramesh_subdivide(octetramesh_t *dst, octetramesh_t *src)
{
  int i;
  int p0, p1, p2, p3, p4, p5;
  int q0, q01, q02, q03, q04, q12, q23, q34, q41, q51, q52, q53, q54, q13;

  graph_t *graph = octetramesh_graph(src);

  int ne = graph->ne;
  int nv = src->nv;
  int nt = src->nt;
  int no = src->no;
  int d = src->d;

  int nv2 = nv + ne + no;  // old vertices + midpoints + octahedral centers
  int nt2 = 4*nt + 8*no;
  int no2 = nt + 6*no;

  // allocate space for the new mesh and copy old vertices, midpoints, and octahedral centers into dst
  octetramesh_new(dst, nv2, nt2, no2, d);
  memcpy(dst->vertices[0], src->vertices[0], nv*d*sizeof(double));
  octetramesh_midpoints(dst->vertices + nv, src, graph);
  octetramesh_octa_centers(dst->vertices + nv + ne, src);

  for (i = 0; i < nt; i++) {    // subdivide each tetrahedron in the original mesh
    // original point indices
    p0 = src->tetrahedra[i][0];
    p1 = src->tetrahedra[i][1];
    p2 = src->tetrahedra[i][2];
    p3 = src->tetrahedra[i][3];

    // new point indices
    q01 = nv + graph_find_edge(graph, p0, p1);
    q02 = nv + graph_find_edge(graph, p0, p2);
    q03 = nv + graph_find_edge(graph, p0, p3);
    q12 = nv + graph_find_edge(graph, p1, p2);
    q13 = nv + graph_find_edge(graph, p1, p3);
    q23 = nv + graph_find_edge(graph, p2, p3);

    // new tetrahedra
    int *t0 = dst->tetrahedra[4*i];
    int *t1 = dst->tetrahedra[4*i+1];
    int *t2 = dst->tetrahedra[4*i+2];
    int *t3 = dst->tetrahedra[4*i+3];

    t0[0] = p0;  t0[1] = q01;  t0[2] = q02;  t0[3] = q03;
    t1[0] = p1;  t1[1] = q01;  t1[2] = q12;  t1[3] = q13;
    t2[0] = p2;  t2[1] = q02;  t2[2] = q12;  t2[3] = q23;
    t3[0] = p3;  t3[1] = q03;  t3[2] = q13;  t3[3] = q23;

    // new octahedron
    int *oct = dst->octahedra[i];

    oct[0] = q01;  oct[1] = q02;  oct[2] = q03;  oct[3] = q13;  oct[4] = q12;  oct[5] = q23;
  }

  for (i = 0; i < no; i++) {    // subdivide each octahedron in the original mesh
    // original point indices
    p0 = src->octahedra[i][0];
    p1 = src->octahedra[i][1];
    p2 = src->octahedra[i][2];
    p3 = src->octahedra[i][3];
    p4 = src->octahedra[i][4];
    p5 = src->octahedra[i][5];

    // new point indices
    q0 = nv + ne + i;
    q01 = nv + graph_find_edge(graph, p0, p1);
    q02 = nv + graph_find_edge(graph, p0, p2);
    q03 = nv + graph_find_edge(graph, p0, p3);
    q04 = nv + graph_find_edge(graph, p0, p4);
    q12 = nv + graph_find_edge(graph, p1, p2);
    q23 = nv + graph_find_edge(graph, p2, p3);
    q34 = nv + graph_find_edge(graph, p3, p4);
    q41 = nv + graph_find_edge(graph, p4, p1);
    q51 = nv + graph_find_edge(graph, p5, p1);
    q52 = nv + graph_find_edge(graph, p5, p2);
    q53 = nv + graph_find_edge(graph, p5, p3);
    q54 = nv + graph_find_edge(graph, p5, p4);

    // new tetrahedra
    int *t0 = dst->tetrahedra[4*nt+8*i];
    int *t1 = dst->tetrahedra[4*nt+8*i+1];
    int *t2 = dst->tetrahedra[4*nt+8*i+2];
    int *t3 = dst->tetrahedra[4*nt+8*i+3];
    int *t4 = dst->tetrahedra[4*nt+8*i+4];
    int *t5 = dst->tetrahedra[4*nt+8*i+5];
    int *t6 = dst->tetrahedra[4*nt+8*i+6];
    int *t7 = dst->tetrahedra[4*nt+8*i+7];

    t0[0] = q0;  t0[1] = q01;  t0[2] = q02;  t0[3] = q12;
    t1[0] = q0;  t1[1] = q02;  t1[2] = q03;  t1[3] = q23;
    t2[0] = q0;  t2[1] = q03;  t2[2] = q04;  t2[3] = q34;
    t3[0] = q0;  t3[1] = q04;  t3[2] = q01;  t3[3] = q41;
    t4[0] = q0;  t4[1] = q51;  t4[2] = q52;  t4[3] = q12;
    t5[0] = q0;  t5[1] = q52;  t5[2] = q53;  t5[3] = q23;
    t6[0] = q0;  t6[1] = q53;  t6[2] = q54;  t6[3] = q34;
    t7[0] = q0;  t7[1] = q54;  t7[2] = q51;  t7[3] = q41;

    // new octahedra
    int *oct0 = dst->octahedra[nt+6*i];
    int *oct1 = dst->octahedra[nt+6*i+1];
    int *oct2 = dst->octahedra[nt+6*i+2];
    int *oct3 = dst->octahedra[nt+6*i+3];
    int *oct4 = dst->octahedra[nt+6*i+4];
    int *oct5 = dst->octahedra[nt+6*i+5];

    oct0[0] = p0;  oct0[1] = q01;  oct0[2] = q02;  oct0[3] = q03;  oct0[4] = q04;  oct0[5] = q0;
    oct1[0] = q01;  oct1[1] = p1;  oct1[2] = q12;  oct1[3] = q0;  oct1[4] = q41;  oct1[5] = q51;
    oct2[0] = q02;  oct2[1] = p2;  oct2[2] = q23;  oct2[3] = q0;  oct2[4] = q12;  oct2[5] = q52;
    oct3[0] = q03;  oct3[1] = p3;  oct3[2] = q34;  oct3[3] = q0;  oct3[4] = q23;  oct3[5] = q53;
    oct4[0] = q04;  oct4[1] = p4;  oct4[2] = q41;  oct4[3] = q0;  oct4[4] = q34;  oct4[5] = q54;
    oct5[0] = p5;  oct5[1] = q51;  oct5[2] = q52;  oct5[3] = q53;  oct5[4] = q54;  oct5[5] = q0;
  }

  graph_free(graph);
}


/**
 * Multi-resolution subdivision based on approximating
 * a scalar function f to a given resolution.
 */
void octetramesh_subdivide_mres(octetramesh_t *dst, octetramesh_t *src,
				double(*f)(double *, void *), void *fdata, double resolution)
{
  int i;
  int nt = src->nt;
  int no = src->no;

  //int tetmask[nt];
  //int octmask[no];
  int *tetmask = (int *)alloca(sizeof(int) * nt);  
  int *octmask = (int *)alloca(sizeof(int) * no);
  memset(tetmask, 0, nt*sizeof(int));
  memset(octmask, 0, no*sizeof(int));

  for (i = 0; i < nt; i++) {    // determine which tetrahedra to subdivide
    double *v0 = src->vertices[ src->tetrahedra[i][0] ];
    double *v1 = src->vertices[ src->tetrahedra[i][1] ];
    double *v2 = src->vertices[ src->tetrahedra[i][2] ];
    double *v3 = src->vertices[ src->tetrahedra[i][3] ];

    double F[4] = {f(v0, fdata), f(v1, fdata), f(v2, fdata), f(v3, fdata)};

    if (arr_max(F,4) - arr_min(F,4) > resolution)
      tetmask[i] = 1;
  }
  
  for (i = 0; i < no; i++) {    // determine which octahedra to subdivide

    double *v0 = src->vertices[ src->octahedra[i][0] ];
    double *v1 = src->vertices[ src->octahedra[i][1] ];
    double *v2 = src->vertices[ src->octahedra[i][2] ];
    double *v3 = src->vertices[ src->octahedra[i][3] ];
    double *v4 = src->vertices[ src->octahedra[i][4] ];
    double *v5 = src->vertices[ src->octahedra[i][5] ];

    double F[6] = {f(v0, fdata), f(v1, fdata), f(v2, fdata), f(v3, fdata), f(v4, fdata), f(v5, fdata)};

    if (arr_max(F,6) - arr_min(F,6) > resolution)
      octmask[i] = 1;
  }
  
  if (count(tetmask, nt) == 0 && count(octmask, no) == 0) {
    memcpy(dst, src, sizeof(octetramesh_t));
    return;
  }

  // subdivide the selected cells
  octetramesh_subdivide_select(dst, src, tetmask, octmask);
}


/*
 * Create (allocate) the contents of a mesh.
 */
void octetramesh_new(octetramesh_t *mesh, int nv, int nt, int no, int d)
{
  mesh->nv = nv;
  mesh->nt = nt;
  mesh->no = no;
  mesh->d = d;

  mesh->vertices = new_matrix2(nv, d);
  mesh->tetrahedra = new_matrix2i(nt, 4);

  if (no > 0)
    mesh->octahedra = new_matrix2i(no, 6);
}


/*
 * Copy the contents of one mesh into another mesh.
 */
void octetramesh_copy(octetramesh_t *dst, octetramesh_t *src)
{
  memcpy(dst->vertices[0], src->vertices[0], src->nv*src->d*sizeof(double));  // copy raw vertices
  memcpy(dst->tetrahedra[0], src->tetrahedra[0], 4*src->nt*sizeof(int));      // copy raw tetrahedra
  memcpy(dst->octahedra[0], src->octahedra[0], 6*src->no*sizeof(int));        // copy raw octahedra
}

void octetramesh_copy_vertices(octetramesh_t *dst, octetramesh_t *src)
{
  memcpy(dst->vertices[0], src->vertices[0], src->nv*src->d*sizeof(double));  // copy raw vertices
}

void octetramesh_copy_tetrahedra(octetramesh_t *dst, octetramesh_t *src)
{
  memcpy(dst->tetrahedra[0], src->tetrahedra[0], 4*src->nt*sizeof(int));      // copy raw tetrahedra
}

void octetramesh_copy_octahedra(octetramesh_t *dst, octetramesh_t *src)
{
  memcpy(dst->octahedra[0], src->octahedra[0], 6*src->no*sizeof(int));        // copy raw octahedra
}


/*
 * Clone a mesh
 */
octetramesh_t *octetramesh_clone(octetramesh_t *src)
{
  octetramesh_t *dst;
  safe_calloc(dst, 1, octetramesh_t);
  octetramesh_new(dst, src->nv, src->nt, src->no, src->d);
  octetramesh_copy(dst, src);

  return dst;
}


/*
 * Free the contents of a octetrahedral mesh.
 */
void octetramesh_free(octetramesh_t *mesh)
{
  if (mesh->nv > 0) {
    free_matrix2(mesh->vertices);
    free_matrix2i(mesh->tetrahedra);
    if (mesh->no > 0)
      free_matrix2i(mesh->octahedra);
  }
}


/*
 * Save a octetrahedral mesh to PLY file.
 */
void octetramesh_save_PLY(octetramesh_t *mesh, char *filename)
{
  FILE *f = fopen(filename, "w");

  int i, i0, i1, i2, i3, i4, i5;

  // ~~~ TODO: Use a graph traversal algorithm to avoid double-counting! ~~~
  int num_faces = 4*mesh->nt + 8*mesh->no;  //tetramesh_face_count(mesh);

  fprintf(f, "ply\n");
  fprintf(f, "format ascii 1.0\n");
  fprintf(f, "comment tetramesh model\n");  
  fprintf(f, "element vertex %d\n", mesh->nv);
  fprintf(f, "property float x\n");
  fprintf(f, "property float y\n");
  fprintf(f, "property float z\n");
  fprintf(f, "element face %d\n", num_faces);
  fprintf(f, "property list uchar int vertex_indices\n");
  fprintf(f, "end_header\n");

  for (i = 0; i < mesh->nv; i++) {
    double x = mesh->vertices[i][0]; // / (mesh->vertices[i][3]);
    double y = mesh->vertices[i][1]; // / (mesh->vertices[i][3]);
    double z = mesh->vertices[i][2]; // / (mesh->vertices[i][3]);
    fprintf(f, "%f %f %f", x, y, z);
    //for (j = 0; j < 3 /*mesh->d*/; j++)
    //  fprintf(f, "%f ", mesh->vertices[i][j]);
    fprintf(f, "\n");
  }

  // ~~~ TODO: Use a graph traversal algorithm to avoid double-counting! ~~~
  for (i = 0; i < mesh->nt; i++) {
    i0 = mesh->tetrahedra[i][0];
    i1 = mesh->tetrahedra[i][1];
    i2 = mesh->tetrahedra[i][2];
    i3 = mesh->tetrahedra[i][3];
    fprintf(f, "3 %d %d %d\n", i0, i1, i2);
    fprintf(f, "3 %d %d %d\n", i0, i1, i3);
    fprintf(f, "3 %d %d %d\n", i0, i2, i3);
    fprintf(f, "3 %d %d %d\n", i1, i2, i3);
  }
  for (i = 0; i < mesh->no; i++) {
    i0 = mesh->octahedra[i][0];
    i1 = mesh->octahedra[i][1];
    i2 = mesh->octahedra[i][2];
    i3 = mesh->octahedra[i][3];
    i4 = mesh->octahedra[i][4];
    i5 = mesh->octahedra[i][5];
    fprintf(f, "3 %d %d %d\n", i0, i1, i2);
    fprintf(f, "3 %d %d %d\n", i0, i2, i3);
    fprintf(f, "3 %d %d %d\n", i0, i3, i4);
    fprintf(f, "3 %d %d %d\n", i0, i4, i1);
    fprintf(f, "3 %d %d %d\n", i5, i1, i2);
    fprintf(f, "3 %d %d %d\n", i5, i2, i3);
    fprintf(f, "3 %d %d %d\n", i5, i3, i4);
    fprintf(f, "3 %d %d %d\n", i5, i4, i1);
  }

  fclose(f);
}


/*
 * Compute stats of a mesh.
 */
octetramesh_stats_t octetramesh_stats(octetramesh_t *T)
{
  int i, j, e, nt = T->nt, d = T->d;
  octetramesh_stats_t stats;

  graph_t *graph = octetramesh_graph(T);
  stats.num_edges = graph->ne;

  stats.num_vertices = T->nv;
  stats.num_tetrahedra = T->nt;
  stats.num_octahedra = T->no;

  stats.min_edge_len = DBL_MAX;
  stats.max_edge_len = 0;
  stats.avg_edge_len = 0;
  stats.std_edge_len = 0;
  stats.min_tetra_skewness = DBL_MAX;
  stats.max_tetra_skewness = 0;
  stats.avg_tetra_skewness = 0;
  stats.std_tetra_skewness = 0;
  stats.min_tetra_volume = DBL_MAX;
  stats.max_tetra_volume = 0;
  stats.avg_tetra_volume = 0;
  stats.std_tetra_volume = 0;

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
    if (skewness < stats.min_tetra_skewness)
      stats.min_tetra_skewness = skewness;
    if (skewness > stats.max_tetra_skewness)
      stats.max_tetra_skewness = skewness;
    stats.avg_tetra_skewness += skewness;

    double volume = tetrahedron_volume(p0, p1, p2, p3, d);
    if (volume < stats.min_tetra_volume)
      stats.min_tetra_volume = volume;
    if (volume > stats.max_tetra_volume)
      stats.max_tetra_volume = volume;
    stats.avg_tetra_volume += volume;
  }
  stats.avg_tetra_skewness /= (double)nt;
  stats.avg_tetra_volume /= (double)nt;

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
    double ds = stats.avg_tetra_skewness - skewness;
    stats.std_tetra_skewness += ds*ds;

    double volume = tetrahedron_volume(p0, p1, p2, p3, d);
    double dv = stats.avg_tetra_volume - volume;
    stats.std_tetra_volume += dv*dv;
  }
  stats.std_tetra_skewness = sqrt(stats.std_tetra_skewness/(double)nt);
  stats.std_tetra_volume = sqrt(stats.std_tetra_volume/(double)nt);

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
void octetramesh_print_stats(octetramesh_stats_t stats)
{
  printf("octetramesh stats {\n");

  printf("  nv = %d, ne = %d, nt = %d, no = %d\n",
	 stats.num_vertices, stats.num_edges, stats.num_tetrahedra, stats.num_octahedra);

  printf("  edge_len = [%f, %f], avg: %f, std: %f\n",
	 stats.min_edge_len, stats.max_edge_len, stats.avg_edge_len, stats.std_edge_len);

  printf("  tetra_skewness = [%f, %f], avg: %f, std: %f\n",
	 stats.min_tetra_skewness, stats.max_tetra_skewness, stats.avg_tetra_skewness, stats.std_tetra_skewness);

  printf("  tetra_volume = [%f, %f], avg: %f, std: %f\n",\
	 stats.min_tetra_volume, stats.max_tetra_volume, stats.avg_tetra_volume, stats.std_tetra_volume);

  printf("}\n");
}
