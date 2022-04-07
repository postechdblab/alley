#ifndef _INTERSECTION_HPP_
#define _INTERSECTION_HPP_

#include <nmmintrin.h>
#include <immintrin.h>
#include <stdint.h>
#include <string.h>
#include <set>
#include <random>
#include "util.h"

#define INITIAL_CAPACITY 10000000
#define CACHE_THRESHOLD 1000

using std::pair;

#ifdef COMPRESS

const int BSR_WIDTH = sizeof(int) * 8;
const int BSR_SHIFT = __builtin_ctzll(BSR_WIDTH); 
const int BSR_MASK  = BSR_WIDTH - 1;

#define BASE(id) (id >> BSR_SHIFT)
#define STATE(id) (1 << (id & BSR_MASK))

void quit();
void align_malloc(void **, size_t, size_t);
    
constexpr int cyclic_shift1 = _MM_SHUFFLE(0,3,2,1); //rotating right
constexpr int cyclic_shift2 = _MM_SHUFFLE(2,1,0,3); //rotating left
constexpr int cyclic_shift3 = _MM_SHUFFLE(1,0,3,2); //between

static const __m128i all_zero_si128 = _mm_setzero_si128();
static const __m128i all_one_si128 = _mm_set_epi32(0xffffffff, 0xffffffff,
        0xffffffff, 0xffffffff);

static const uint8_t shuffle_pi8_array[256] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 255, 255, 255, 255, 255, 255, 255, 255, 
    8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 255, 255, 255, 255, 
    12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 255, 255, 255, 255, 
    8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
};
static const __m128i *shuffle_mask = (__m128i*)(shuffle_pi8_array);

uint32_t* prepare_shuffling_dict_avx();
static const uint32_t *shuffle_mask_avx =prepare_shuffling_dict_avx();

int* prepare_byte_check_mask_dict2();
static const int *byte_check_mask_dict = prepare_byte_check_mask_dict2();

uint8_t * prepare_match_shuffle_dict2();
static const __m128i *match_shuffle_dict = (__m128i*)prepare_match_shuffle_dict2();

static const uint8_t byte_check_group_a_pi8[64] = {
    0, 0, 0, 0, 4, 4, 4, 4, 8, 8, 8, 8, 12, 12, 12, 12,
    1, 1, 1, 1, 5, 5, 5, 5, 9, 9, 9, 9, 13, 13, 13, 13,
    2, 2, 2, 2, 6, 6, 6, 6, 10, 10, 10, 10, 14, 14, 14, 14,
    3, 3, 3, 3, 7, 7, 7, 7, 11, 11, 11, 11, 15, 15, 15, 15,
};
static const uint8_t byte_check_group_b_pi8[64] = {
    0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12,
    1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13,
    2, 6, 10, 14, 2, 6, 10, 14, 2, 6, 10, 14, 2, 6, 10, 14,
    3, 7, 11, 15, 3, 7, 11, 15, 3, 7, 11, 15, 3, 7, 11, 15,
};
static const __m128i *byte_check_group_a_order = (__m128i*)(byte_check_group_a_pi8);
static const __m128i *byte_check_group_b_order = (__m128i*)(byte_check_group_b_pi8);

constexpr int word_check_shuffle_a01 = _MM_SHUFFLE(1,1,0,0); 
constexpr int word_check_shuffle_a23 = _MM_SHUFFLE(3,3,2,2); 
constexpr int word_check_shuffle_b01 = _MM_SHUFFLE(1,0,1,0); 
constexpr int word_check_shuffle_b23 = _MM_SHUFFLE(3,2,3,2); 

uint8_t* prepare_shuffling_dict_u16();
static const __m128i *shuffle_mask16 = (__m128i *)prepare_shuffling_dict_u16();

int intersect_qfilter_bsr_b4_v2(const int* bases_a, const int* states_a, int size_a,
            const int* bases_b, const int* states_b, int size_b,
            int* bases_c, int* states_c);
int intersect_simdgalloping_bsr(const int* bases_a, const int* states_a, int size_a,
        const int* bases_b, const int* states_b, int size_b,
        int* bases_c, int* states_c);

int offline_uint_trans_bsr(int *set_a, int size_a, int *bases_a, int *states_a);

int offline_bsr_trans_uint(int *bases_a, int *states_a, int size_a, int *set_a);

#else

static __m128i shuffle_mask[16];

int getBit(int value, int position);
void prepare_shuffling_dictionary();
inline size_t scalar(const int32_t *a, const size_t lena,
		const int32_t *b, const size_t lenb, int32_t *out);
#endif



///////////////////////////////////////// APIs ////////////////////////////////////////

#ifdef COMPRESS
//assume size_a <= size_b 
int intersect(const int* bases_a, const int* states_a, int size_a,
              const int* bases_b, const int* states_b, int size_b,
              int* bases_c, int* states_c);
//binary
int intersect(const range& RA, const range& RB, int* CB, int* CS);
//assume that R has >= 2 ranges
int intersect(const std::vector<range>& R, int* CB, int* CS, int* TB, int* TS);
#else
int intersect(const range& RA, const range& RB, int* C); 
int intersect(const meta_range& RA, const meta_range& RB, int* C); 
int intersect(int *V, size_t S, const range& RB, int* C);
int intersect(int *V, size_t S, const meta_range& RB, int* C);
#ifdef OPT
int intersect(std::vector<range>& R, int* C, int* temp);
#else
int intersect(std::set<range>& R, int* C, int* temp);
int intersect(std::set<meta_range>& R, int* C, int* temp, AdjElem& failed_edge);
#endif
#endif

int do_union(const range& RA, const range& RB, int* C); 
int do_union(int *V, size_t S, const range& RB, int* C);
int do_union(std::vector<range>& R, int* C);
int do_union(std::set<range>& R, int* C);

#endif
