 /* gssw.i */
%include "stdint.i"
%include "typemaps.i"

%typemap(in,numinputs=0,noblock=1) size_t *len  {
  size_t templen;
  $1 = &templen;
}

%typemap(out) int* test_wrapper {
  int i;
  $result = PyList_New(templen);
  for (i = 0; i < templen; i++) {
    PyObject *o = PyFloat_FromDouble((double)$1[i]);
    PyList_SetItem($result,i,o);
  }
}


 %module gssw 
 %{
#include "gssw.h"


struct gssw_graph_cigar{
    int32_t length;
    gssw_cigar_element* elements;
};

struct gssw_graph_mapping{
    int32_t position; // position in first node
    int16_t score;
    gssw_graph_cigar cigar;
};


extern int* test_wrapper(gssw_graph* graph,
                                           const char* read,
                                           int32_t readLen,
                                           int8_t* nt_table,
                                           int8_t* score_matrix,
                                           uint8_t gap_open,
                                           uint8_t gap_extension,
                                           int8_t start_full_length_bonus,
                                           int8_t end_full_length_bonus,
                                           size_t *len
                                           );

extern gssw_profile* gssw_init (const int8_t* read, const int32_t readLen, const int8_t* mat, const int32_t n,
                         int8_t start_full_length_bonus, int8_t end_full_length_bonus, const int8_t score_size);

extern int8_t* gssw_create_score_matrix(int32_t match, int32_t mismatch);
extern int8_t* gssw_create_nt_table(void);

extern gssw_node* gssw_node_create(void* data,
                            const uint32_t id,
                            const char* seq,
                            const int8_t* nt_table,
                            const int8_t* score_matrix);
extern void gssw_node_destroy(gssw_node* n);
extern void gssw_node_add_prev(gssw_node* n, gssw_node* m);
extern void gssw_node_add_next(gssw_node* n, gssw_node* m);
extern void gssw_nodes_add_edge(gssw_node* n, gssw_node* m);
extern void gssw_node_del_prev(gssw_node* n, gssw_node* m);
extern void gssw_node_del_next(gssw_node* n, gssw_node* m);
extern void gssw_nodes_del_edge(gssw_node* n, gssw_node* m);
extern void gssw_node_replace_prev(gssw_node* n, gssw_node* m, gssw_node* p);
extern void gssw_node_replace_next(gssw_node* n, gssw_node* m, gssw_node* p);


extern gssw_node* gssw_node_fill (gssw_node* node,
                const gssw_profile* prof,
                const uint8_t weight_gapO,
                const uint8_t weight_gapE,
                const int32_t maskLen,
                bool save_matrixes,
                const gssw_seed* seed);

extern gssw_graph* gssw_graph_fill (gssw_graph* graph,
                 const char* read_seq,
                 const int8_t* nt_table,
                 const int8_t* score_matrix,
                 const uint8_t weight_gapO,
                 const uint8_t weight_gapE,
                 const int8_t start_full_length_bonus,
                 const int8_t end_full_length_bonus,
                 const int32_t maskLen,
                 const int8_t score_size,
                 bool save_matrixes);

extern gssw_graph* gssw_graph_fill_qual_adj(gssw_graph* graph,
                         const char* read_seq,
                         const char* read_qual,
                         const int8_t* nt_table,
                         const int8_t* adj_score_matrix,
                         const uint8_t weight_gapO,
                         const uint8_t weight_gapE,
                         const int8_t start_full_length_bonus,
                         const int8_t end_full_length_bonus,
                         const int32_t maskLen,
                         const int8_t score_size,
                         bool save_matrixes);
    
    
extern gssw_graph* gssw_graph_create(uint32_t size);

extern int32_t gssw_graph_add_node(gssw_graph* graph,
                            gssw_node* node);

extern void gssw_graph_print_score_matrices(gssw_graph* graph,
                                     const char* read,
                                     int32_t readLen,
                                     FILE* out);

extern gssw_graph_mapping* gssw_graph_trace_back (gssw_graph* graph,
                                           const char* read,
                                           int32_t readLen,
                                           int8_t* nt_table,
                                           int8_t* score_matrix,
                                           uint8_t gap_open,
                                           uint8_t gap_extension,
                                           int8_t start_full_length_bonus,
                                           int8_t end_full_length_bonus);                                    

%}

#include "gssw.h"
extern gssw_profile* gssw_init (const int8_t* read, const int32_t readLen, const int8_t* mat, const int32_t n,
                         int8_t start_full_length_bonus, int8_t end_full_length_bonus, const int8_t score_size);

extern int8_t* gssw_create_score_matrix(int32_t match, int32_t mismatch);
extern int8_t* gssw_create_nt_table(void);

extern gssw_node* gssw_node_create(void* data,
                            const uint32_t id,
                            const char* seq,
                            const int8_t* nt_table,
                            const int8_t* score_matrix);
extern void gssw_node_destroy(gssw_node* n);
extern void gssw_node_add_prev(gssw_node* n, gssw_node* m);
extern void gssw_node_add_next(gssw_node* n, gssw_node* m);
extern void gssw_nodes_add_edge(gssw_node* n, gssw_node* m);
extern void gssw_node_del_prev(gssw_node* n, gssw_node* m);
extern void gssw_node_del_next(gssw_node* n, gssw_node* m);
extern void gssw_nodes_del_edge(gssw_node* n, gssw_node* m);
extern void gssw_node_replace_prev(gssw_node* n, gssw_node* m, gssw_node* p);
extern void gssw_node_replace_next(gssw_node* n, gssw_node* m, gssw_node* p);


extern gssw_node* gssw_node_fill (gssw_node* node,
                const gssw_profile* prof,
                const uint8_t weight_gapO,
                const uint8_t weight_gapE,
                const int32_t maskLen,
                bool save_matrixes,
                const gssw_seed* seed);

extern gssw_graph* gssw_graph_fill (gssw_graph* graph,
                 const char* read_seq,
                 const int8_t* nt_table,
                 const int8_t* score_matrix,
                 const uint8_t weight_gapO,
                 const uint8_t weight_gapE,
                 const int8_t start_full_length_bonus,
                 const int8_t end_full_length_bonus,
                 const int32_t maskLen,
                 const int8_t score_size,
                 bool save_matrixes);

extern gssw_graph* gssw_graph_fill_qual_adj(gssw_graph* graph,
                         const char* read_seq,
                         const char* read_qual,
                         const int8_t* nt_table,
                         const int8_t* adj_score_matrix,
                         const uint8_t weight_gapO,
                         const uint8_t weight_gapE,
                         const int8_t start_full_length_bonus,
                         const int8_t end_full_length_bonus,
                         const int32_t maskLen,
                         const int8_t score_size,
                         bool save_matrixes);
    
extern gssw_graph* gssw_graph_create(uint32_t size); 


extern int32_t gssw_graph_add_node(gssw_graph* graph,
                            gssw_node* node);

extern void gssw_graph_print_score_matrices(gssw_graph* graph,
                                     const char* read,
                                     int32_t readLen,
                                     FILE* out);
                                     
extern gssw_graph_mapping* gssw_graph_trace_back (gssw_graph* graph,
                                           const char* read,
                                           int32_t readLen,
                                           int8_t* nt_table,
                                           int8_t* score_matrix,
                                           uint8_t gap_open,
                                           uint8_t gap_extension,
                                           int8_t start_full_length_bonus,
                                           int8_t end_full_length_bonus);



extern int* test_wrapper(gssw_graph* graph,
                                           const char* read,
                                           int32_t readLen,
                                           int8_t* nt_table,
                                           int8_t* score_matrix,
                                           uint8_t gap_open,
                                           uint8_t gap_extension,
                                           int8_t start_full_length_bonus,
                                           int8_t end_full_length_bonus,
                                           size_t *len
                                           );

