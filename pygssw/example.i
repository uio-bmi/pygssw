 /* interface.i */
 %module gssw 
 %{
#include "gssw.h"
extern gssw_profile* gssw_init (const int8_t* read, const int32_t readLen, const int8_t* mat, const int32_t n,
                         int8_t start_full_length_bonus, int8_t end_full_length_bonus, const int8_t score_size);

%}

#include "gssw.h"
extern gssw_profile* gssw_init (const int8_t* read, const int32_t readLen, const int8_t* mat, const int32_t n,
                         int8_t start_full_length_bonus, int8_t end_full_length_bonus, const int8_t score_size);

 
